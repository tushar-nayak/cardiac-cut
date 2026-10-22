#!/usr/bin/env python3
"""Local browser app for reviewing and manually trimming cardiac surface splits."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import re
import threading
import uuid
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np
from plotly.offline import get_plotlyjs


ROOT = Path(__file__).resolve().parents[1]


def _cases(split_dir: Path) -> list[str]:
    return sorted(p.name.removesuffix("_inner.stl")
                  for p in split_dir.glob("*_inner.stl"))


def _path_pair(split_dir: Path, case: str):
    inner = split_dir / f"{case}_inner.stl"
    outer = split_dir / f"{case}_outer_untrimmed.stl"
    if not inner.is_file() or not outer.is_file():
        raise FileNotFoundError(
            f"Expected {inner.name} and {outer.name}; regenerate splits with the updated splitter."
        )
    return inner, outer


def _read_mesh(path: Path):
    suffix = path.suffix.lower()
    if suffix == ".stl":
        from vtkmodules.vtkIOGeometry import vtkSTLReader
        reader = vtkSTLReader()
        reader.SetFileName(str(path))
        reader.Update()
        output = reader.GetOutput()
    else:
        # Convert non-STL input meshes through trimesh so clipping and preview
        # use the same triangle-surface representation for STL, OBJ, and PLY.
        import trimesh
        mesh = trimesh.load(str(path), force="mesh", process=True)
        if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
            raise ValueError(f"No triangle surface found in {path.name}")
        reader = _reader_from_trimesh(mesh)
        output = reader.GetOutput()
    if output.GetNumberOfPoints() == 0 or output.GetNumberOfCells() == 0:
        raise ValueError(f"Empty or unreadable mesh: {path}")
    return reader


def _poly_payload(poly) -> dict:
    if poly.GetNumberOfPoints() == 0 or poly.GetNumberOfCells() == 0:
        return {"x": [], "y": [], "z": [], "i": [], "j": [], "k": [], "n_faces": 0}
    # STL surfaces are triangles already. Convert VTK's contiguous buffers in
    # bulk instead of crossing the Python/VTK boundary once per point and cell.
    from vtkmodules.vtkFiltersCore import vtkTriangleFilter
    from vtkmodules.util.numpy_support import vtk_to_numpy
    tri = vtkTriangleFilter()
    tri.SetInputData(poly)
    tri.PassVertsOff()
    tri.PassLinesOff()
    tri.Update()
    mesh = tri.GetOutput()
    pts = vtk_to_numpy(mesh.GetPoints().GetData()).astype(np.float32, copy=False)
    polys = mesh.GetPolys()
    offsets = vtk_to_numpy(polys.GetOffsetsArray()).astype(np.int64, copy=False)
    connectivity = vtk_to_numpy(polys.GetConnectivityArray()).astype(np.int32, copy=False)
    if not np.all(np.diff(offsets) == 3):
        raise ValueError("Surface contains non-triangular cells")
    faces = connectivity.reshape(-1, 3)
    return {"x": pts[:, 0].tolist(), "y": pts[:, 1].tolist(), "z": pts[:, 2].tolist(),
            "i": faces[:, 0].tolist(), "j": faces[:, 1].tolist(), "k": faces[:, 2].tolist(),
            "n_faces": int(len(faces))}


def _poly_mesh(poly):
    """Convert a VTK triangle surface to a welded mesh for QC and measurements."""
    import trimesh
    from vtkmodules.vtkFiltersCore import vtkTriangleFilter
    from vtkmodules.util.numpy_support import vtk_to_numpy
    triangle = vtkTriangleFilter()
    triangle.SetInputData(poly)
    triangle.PassVertsOff()
    triangle.PassLinesOff()
    triangle.Update()
    mesh = triangle.GetOutput()
    if mesh.GetPoints() is None or mesh.GetNumberOfCells() == 0:
        return trimesh.Trimesh()
    vertices = vtk_to_numpy(mesh.GetPoints().GetData()).astype(np.float64, copy=False)
    polys = mesh.GetPolys()
    faces = vtk_to_numpy(polys.GetConnectivityArray()).astype(np.int64, copy=False).reshape(-1, 3)
    return trimesh.Trimesh(vertices=vertices, faces=faces, process=True)


def _surface_diagnostics(surface):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    components = surface.split(only_watertight=False)
    counts = np.bincount(surface.edges_unique_inverse, minlength=len(surface.edges_unique))
    boundary = surface.edges_unique[counts == 1]
    if len(boundary):
        used = np.unique(boundary)
        remap = np.full(len(surface.vertices), -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        edges = remap[boundary]
        degree = np.bincount(edges.ravel(), minlength=len(used))
        graph = coo_matrix((np.ones(2 * len(edges)),
                            (np.r_[edges[:, 0], edges[:, 1]],
                             np.r_[edges[:, 1], edges[:, 0]])),
                           shape=(len(used), len(used))).tocsr()
        loops = int(connected_components(graph, directed=False, return_labels=False))
        simple_loops = bool(np.all(degree == 2))
    else:
        loops, simple_loops = 0, True
    return {"faces": int(len(surface.faces)), "components": int(len(components)),
            "area_mm2": float(surface.area), "euler": int(surface.euler_number),
            "boundary_edges": int(len(boundary)), "boundary_loops": loops,
            "simple_boundary_loops": simple_loops,
            "watertight": bool(surface.is_watertight)}


def _section_measurement(surface, origin, normal):
    section = surface.section(plane_origin=origin, plane_normal=normal)
    if section is None:
        return {"area_mm2": 0.0, "perimeter_mm": 0.0,
                "major_diameter_mm": None, "minor_diameter_mm": None}
    planar, _ = section.to_2D()
    try:
        polygons = planar.polygons_full
    except (ValueError, TypeError, NotImplementedError):
        return {"area_mm2": None, "perimeter_mm": float(planar.length),
                "major_diameter_mm": None, "minor_diameter_mm": None}
    area = float(sum(p.area for p in polygons))
    perimeter = float(sum(p.length for p in polygons))
    major = minor = None
    if polygons:
        rect = max(polygons, key=lambda p: p.area).minimum_rotated_rectangle
        xy = np.asarray(rect.exterior.coords, dtype=float)
        lengths = np.linalg.norm(np.diff(xy, axis=0), axis=1)
        if len(lengths) >= 4:
            major, minor = float(np.max(lengths)), float(np.min(lengths))
    return {"area_mm2": area, "perimeter_mm": perimeter,
            "major_diameter_mm": major, "minor_diameter_mm": minor}


def _preview_surface(poly):
    """Use the original surface for previews; never decimate its geometry."""
    return poly


def _initial_plane(poly, case: str | None = None, outer_path: Path | None = None,
                   guidance: str = "auto", source_path: Path | None = None) -> dict:
    pts = np.asarray([poly.GetPoint(i) for i in range(poly.GetNumberOfPoints())])
    mean = pts.mean(axis=0)
    source_path = source_path or (ROOT / "dataset" / f"{case}.stl" if case else None)
    method = "manual plane starting position" if guidance == "manual" else "PCA fallback"
    origin = None
    normal = None
    if guidance == "auto" and source_path and source_path.is_file() and outer_path and outer_path.is_file():
        try:
            # Reuse the splitter's own basal-rim plane estimator so the viewer
            # starts with the same normal-based plane used by --top-cut plane.
            import trimesh
            splitter_path = Path(__file__).with_name("split_cardiac_surfaces.py")
            spec = importlib.util.spec_from_file_location("cardiac_surface_splitter", splitter_path)
            splitter = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(splitter)
            main_mesh = trimesh.load_mesh(source_path, process=True)
            outer_mesh = trimesh.load_mesh(outer_path, process=True)
            origin, normal = splitter._basal_plane(main_mesh, outer_mesh, 3.0)
            method = "splitter basal-rim plane (3 mm margin)"
        except Exception:
            # Keep the viewer usable if a case lacks a detectable rim.
            origin = normal = None
    if origin is None or normal is None:
        _, vectors = np.linalg.eigh(np.cov((pts - mean).T))
        normal = vectors[:, -1]
        dominant = int(np.argmax(np.abs(normal)))
        if normal[dominant] < 0:
            normal *= -1
        projection = (pts - mean) @ normal
        origin = mean + normal * (projection.max() - 0.12 * np.ptp(projection))
    else:
        normal = np.asarray(normal, dtype=float)
        normal /= np.linalg.norm(normal) + 1e-12
        origin = np.asarray(origin, dtype=float)
    projection = (pts - mean) @ normal
    return {
        "origin": origin.tolist(), "normal": normal.tolist(),
        "axis_center": mean.tolist(),
        "axis_range": [float(projection.min()), float(projection.max())],
        "method": method,
    }


def _clip(outer_reader, origin, normal, keep_side):
    from vtkmodules.vtkCommonDataModel import vtkPlane
    from vtkmodules.vtkFiltersCore import vtkClipPolyData
    plane = vtkPlane()
    plane.SetOrigin(*origin)
    plane.SetNormal(*normal)
    clip = vtkClipPolyData()
    clip.SetInputConnection(outer_reader.GetOutputPort())
    clip.SetClipFunction(plane)
    clip.SetValue(0.0)
    # VTK's default keeps the positive half-space; InsideOut keeps the negative one.
    clip.SetInsideOut(keep_side == "negative")
    clip.Update()
    return clip


def _reader_from_poly(poly):
    from vtkmodules.vtkCommonExecutionModel import vtkTrivialProducer
    producer = vtkTrivialProducer()
    producer.SetOutput(poly)
    return producer


def _manual_plane_split(source_path: Path, origin, normal, keep_side="negative"):
    """Split a raw closed shell at the clinician's plane, without auto seeds."""
    import trimesh
    splitter_path = Path(__file__).with_name("split_cardiac_surfaces.py")
    spec = importlib.util.spec_from_file_location("cardiac_surface_splitter_manual", splitter_path)
    splitter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(splitter)
    mesh = trimesh.load(str(source_path), force="mesh", process=True)
    main = max(mesh.split(only_watertight=False), key=lambda part: len(part.faces))
    clip_normal = np.asarray(normal) if keep_side == "negative" else -np.asarray(normal)
    clipped = splitter._clip_below(main, np.asarray(origin), clip_normal)
    pieces = sorted(clipped.split(only_watertight=False), key=lambda part: part.area, reverse=True)
    total_area = sum(part.area for part in pieces)
    pieces = [part for part in pieces if part.area >= 0.02 * total_area]
    if len(pieces) != 2:
        raise ValueError(
            f"This plane produced {len(pieces)} substantial surface pieces; adjust its position or angle until it separates the two walls."
        )
    # The plane has disconnected the walls. Use their radial normal direction
    # around the cut axis only to label the two components; no center/normal seed
    # search or geodesic partition is used in this mode.
    def radial_orientation(piece):
        offset = piece.triangles_center - np.asarray(origin)
        radial = offset - np.outer(offset @ normal, normal)
        return float(np.average(np.einsum("ij,ij->i", piece.face_normals, radial),
                                weights=piece.area_faces))
    inner, outer = sorted(pieces, key=radial_orientation)
    return inner, outer


def _reader_from_trimesh(mesh):
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    points = vtkPoints()
    points.SetData(numpy_to_vtk(np.asarray(mesh.vertices, dtype=np.float64), deep=True))
    faces = np.asarray(mesh.faces, dtype=np.int64)
    offsets = np.arange(0, 3 * len(faces) + 1, 3, dtype=np.int64)
    cells = vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(offsets, deep=True),
                  numpy_to_vtkIdTypeArray(faces.reshape(-1), deep=True))
    poly = vtkPolyData()
    poly.SetPoints(points)
    poly.SetPolys(cells)
    return _reader_from_poly(poly)


def _plane_quad(origin, normal, half_extent):
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n) + 1e-12
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(n, ref); u /= np.linalg.norm(u) + 1e-12
    v = np.cross(n, u); v /= np.linalg.norm(v) + 1e-12
    o = np.asarray(origin, dtype=float)
    q = np.asarray([o - half_extent*u - half_extent*v,
                    o + half_extent*u - half_extent*v,
                    o + half_extent*u + half_extent*v,
                    o - half_extent*u + half_extent*v])
    return {"x": q[:, 0].tolist(), "y": q[:, 1].tolist(), "z": q[:, 2].tolist()}


PAGE = r'''<!doctype html>
<html><head><meta charset="utf-8"><title>Cardiac surface review</title>
<script src="/plotly.min.js"></script>
<style>
*{box-sizing:border-box}body{font:14px/1.4 system-ui,-apple-system,Segoe UI,sans-serif;margin:0;color:#18212b;background:#edf2f6}
header{padding:13px 22px;background:#17324d;color:white;display:flex;justify-content:space-between;align-items:center}
h1{font-size:19px;margin:0}header p{margin:2px 0 0;color:#c8d7e5;font-size:12px}
main{padding:12px 16px;display:grid;grid-template-rows:auto auto minmax(410px,1fr) auto auto auto;gap:10px;height:calc(100vh - 62px);min-height:650px}
.card{background:white;border:1px solid #d8e0e8;border-radius:9px;padding:11px 14px;box-shadow:0 1px 2px #132b3d0a}
.topbar{display:flex;align-items:center;gap:12px;flex-wrap:wrap}.topbar .case{flex:1;min-width:250px}.step{font-weight:700;color:#1769aa;margin-right:5px}
select,button,input{font:inherit}select{padding:7px 9px;border:1px solid #b9c7d3;border-radius:6px;background:white;min-width:210px}
button{padding:8px 13px;border:0;border-radius:6px;background:#1769aa;color:white;font-weight:650;cursor:pointer}
button:hover{background:#0f558d}button.secondary{background:#eaf0f5;color:#263746}button.secondary:hover{background:#dce6ef}
button:disabled{opacity:.5;cursor:wait}.adjust{display:grid;grid-template-columns:1.25fr 1fr 1fr 1.3fr auto;gap:16px;align-items:center}
.ctrl{min-width:0}.ctrl label{display:flex;justify-content:space-between;gap:8px;font-size:12px;font-weight:650;color:#34495a;margin-bottom:4px}
.ctrl output{font-variant-numeric:tabular-nums;color:#1769aa}.ctrl input[type=range]{width:100%;accent-color:#1769aa}
.keep{display:flex;align-items:center;gap:8px;white-space:nowrap}.keep select{min-width:190px}
.view-tools{display:flex;align-items:center;gap:8px;flex-wrap:wrap}.view-tools .label{font-size:12px;color:#526272;font-weight:650;margin-right:2px}.view-tools select{min-width:130px;padding:6px 8px}.view-tools button{padding:6px 10px;font-size:12px}.view-tools label{display:flex;align-items:center;gap:5px;font-size:12px;color:#526272;cursor:pointer}.view-tools input[type=checkbox]{accent-color:#1769aa}.view-tools input[type=number]{width:82px;padding:5px 6px;border:1px solid #b9c7d3;border-radius:5px}
.inspection{display:grid;grid-template-columns:1fr 1fr;gap:8px 18px;align-items:center}.inspection h3{font-size:13px;margin:0 0 4px}.inspection .diag{font-size:12px;color:#526272}.inspection .diag.bad{color:#a32626;font-weight:650}.inspection .measure{display:flex;align-items:center;gap:9px;flex-wrap:wrap}.inspection input[type=range]{min-width:160px;flex:1;accent-color:#1769aa}.inspection output{font-variant-numeric:tabular-nums;color:#1769aa}.inspection button{padding:6px 10px;font-size:12px}
.views{display:grid;grid-template-columns:1fr 1fr;gap:10px;min-height:0}.view{background:white;border:1px solid #d8e0e8;border-radius:9px;overflow:hidden;display:flex;flex-direction:column;min-height:0}
.view-head{padding:8px 13px 0;display:flex;align-items:baseline;justify-content:space-between}.view h2{font-size:14px;margin:0}.view-head small{color:#718292}.plot{flex:1;min-height:0}
.footer{display:flex;align-items:center;gap:12px;flex-wrap:wrap}.status{color:#3b5265;flex:1;min-width:220px}.status.error{color:#a32626}.hint{color:#5c6b79;font-size:12px}
details{border-top:1px solid #e3e9ee;margin-top:9px;padding-top:7px;color:#526272;font-size:12px}details summary{cursor:pointer;font-weight:650;color:#405464}
.coords{display:flex;gap:8px;margin-top:8px}.coords label{display:flex;align-items:center;gap:5px}.coords input{width:90px;padding:4px}
.upload-button{display:inline-flex;align-items:center;gap:6px;padding:8px 11px;border:1px solid #b9c7d3;border-radius:6px;background:#eaf0f5;color:#263746;font-weight:650;cursor:pointer;white-space:nowrap}.upload-button:hover{background:#dce6ef}.upload-button input{display:none}
@media(max-width:1000px){main{height:auto;min-height:calc(100vh - 62px);grid-template-rows:auto auto auto auto auto auto}.adjust{grid-template-columns:repeat(2,1fr)}.views{grid-template-columns:1fr}.plot{height:55vh}.inspection{grid-template-columns:1fr}}
</style></head><body><header><div><h1>Cardiac surface review</h1><p>Position the cut plane, check the outer-wall preview, then save the reviewed pair.</p></div></header><main>
<section class="card topbar"><label class="upload-button" title="STL, OBJ, or PLY; coordinates are interpreted as millimeters">Upload mesh<input id="meshUpload" type="file" accept=".stl,.obj,.ply"></label><label class="case"><span class="step">1</span>Choose mesh <select id="case"></select></label><button class="secondary" id="refreshInputs" title="Rescan the local input folder">Refresh inputs</button>
<label class="keep">Split method <select id="planeMode" title="Choose how inner and outer surfaces are separated"><option value="auto" selected>Center + normals</option><option value="manual">Manual plane</option></select></label>
<label class="keep">Keep <select id="side"><option value="negative">Main wall (suggested)</option><option value="positive">Opposite side</option></select></label>
<button class="secondary" id="reset">Reset plane</button><button id="save">Save this mesh</button></section>
<section class="card adjust">
<div class="hint" style="grid-column:1/-1">Center + normals automatically separates inner and outer walls, then suggests a top cut. Manual plane skips that splitter: position the plane to cut the raw mesh into the two walls.</div>
<div class="ctrl"><label><span><span class="step">2</span>Move plane</span><output id="offv">0 mm</output></label><input id="off" type="range" min="-100" max="100" step="0.5"></div>
<div class="ctrl"><label><span>Rotate around Z</span><output id="azv">0°</output></label><input id="az" type="range" min="-180" max="180" step="1"></div>
<div class="ctrl"><label><span>Tilt plane</span><output id="elv">0°</output></label><input id="el" type="range" min="-89" max="89" step="1"></div>
<div class="hint">Drag the mesh to rotate the view. The red sheet marks the cut.</div><button class="secondary" id="preview">Update preview</button>
<details style="grid-column:1/-1"><summary>Precise plane coordinates</summary><div class="coords"><label>X <input id="ox" type="number" step="0.5"></label><label>Y <input id="oy" type="number" step="0.5"></label><label>Z <input id="oz" type="number" step="0.5"></label><span class="hint">Editing a coordinate resets the move slider.</span></div></details>
</section>
<section class="views"><div class="view"><div class="view-head"><h2>Inner surface</h2><small id="innerNote">cut preview · save behavior shown in controls</small></div><div id="inner" class="plot"></div></div>
<div class="view"><div class="view-head"><h2>Outer surface</h2><small>gray = removed · red = cut plane</small></div><div id="outer" class="plot"></div></div></section>
<section class="card view-tools"><span class="label">Mesh</span><select id="displayMode"><option value="cut">Plane preview</option><option id="campOption" value="camp" disabled>CAMP patches</option></select><span class="label">CAMP field</span><select id="campField"><option value="integrated_K">Integrated curvature</option><option value="avg_kg">Mean Gaussian curvature</option><option value="eta">Eta</option><option value="R2">R²</option><option value="patch_id">Patch ID</option><option value="area_mm2">Patch area</option><option value="quality_status">Fit pass/fail</option></select><label><input id="campEdges" type="checkbox" checked> Patch edges</label><label>Remesh <input id="targetFaces" type="number" min="5000" max="200000" step="5000" value="50000"></label><button id="runCamp" disabled>Save cut + Run CAMP</button><span class="label" style="margin-left:8px">Color</span><select id="colorMode" aria-label="Mesh color"><option value="axis">Long-axis</option><option value="anatomical">Anatomical</option></select><label><input id="contours" type="checkbox"> Surface contours</label><label title="For raw-mesh runs, apply the same plane to both saved surfaces"><input id="applyInner" type="checkbox"> Cut inner on save</label><span class="label">View</span><button class="secondary" data-view="fit">Fit</button><button class="secondary" data-view="front">Front</button><button class="secondary" data-view="side">Side</button><button class="secondary" data-view="top">Top</button></section>
<section class="card inspection"><div><h3>Split quality</h3><div id="quality" class="diag">Quality details appear after the preview updates.</div></div><div><h3>Cross-section measurement</h3><div class="measure"><span>From plane</span><input id="measureOffset" type="range" min="0" max="80" step="1" value="10"><output id="measureOffsetValue">10 mm</output><button class="secondary" id="measureButton">Measure</button></div><div id="measurement" class="diag">Move the measurement plane toward the apex and measure inner and outer cross-sections.</div></div><div style="grid-column:1/-1"><details><summary>3D visibility</summary><div class="coords"><label><input id="showInner" type="checkbox" checked> Inner surface</label><label><input id="showOuter" type="checkbox" checked> Outer surface</label><label><input id="showRemoved" type="checkbox" checked> Removed region</label><label><input id="showPlane" type="checkbox" checked> Cut plane</label><label><input id="showMeasurePlane" type="checkbox" checked> Measurement plane</label></div></details></div></section>
<section class="card footer"><span class="step">3</span><span id="status" class="status">Loading meshes…</span><span class="hint">Meshes in input/ appear in the dropdown. Uploads are saved there; reviewed pairs go to results/manual_surface_split.</span></section>
</main><script>
const byId=id=>document.getElementById(id), status=byId('status'); let caseData=null, campData=null, timer=null, baseOrigin=[0,0,0], normal=[1,0,0], bounds=null, latestRequest=0, previewAbort=null, lastPreview=null, measureResult=null, lastCamera=null, lastDragMode='orbit', cameraSyncInstalled=false, cameraSyncBusy=0;
let cameraFollowSource=null,cameraFollowFrame=null,cameraFollowUntil=0,cameraPointerDown=false;
function vecFromAngles(){let a=Number(byId('az').value)*Math.PI/180,e=Number(byId('el').value)*Math.PI/180;return [Math.cos(e)*Math.cos(a),Math.cos(e)*Math.sin(a),Math.sin(e)]}
function currentOrigin(){return [Number(byId('ox').value),Number(byId('oy').value),Number(byId('oz').value)]}
function meshTrace(m,showScale=false,opacity=1,neutral=false){let c=caseData.plane.axis_center,a=caseData.plane.normal;let intensity=m.x.map((x,i)=>(x-c[0])*a[0]+(m.y[i]-c[1])*a[1]+(m.z[i]-c[2])*a[2]);let range=caseData.plane.axis_range,mode=byId('colorMode').value,contourOn=byId('contours').checked;let trace={type:'mesh3d',x:m.x,y:m.y,z:m.z,i:m.i,j:m.j,k:m.k,opacity,flatshading:false,lighting:{ambient:.48,diffuse:.9,specular:.35,roughness:.55,fresnel:.12},lightposition:{x:1000,y:1500,z:2500},showscale:false};if(neutral){trace.color='#718096'}else if(mode==='anatomical'&&!contourOn){trace.color='#c98278'}else{trace.intensity=intensity;trace.intensitymode='vertex';trace.colorscale=mode==='axis'?'Turbo':[[0,'#c98278'],[1,'#c98278']];trace.cmin=range[0];trace.cmax=range[1];trace.showscale=mode==='axis'&&showScale;if(trace.showscale)trace.colorbar={title:{text:'Long-axis<br>position'},thickness:12,len:.55,x:1.02}}if(contourOn&&!neutral){trace.contour={show:true,color:'#354b5d',width:1.5}}return trace}
const CAMP_LABELS={integrated_K:'Integrated Gaussian curvature',avg_kg:'Mean Gaussian curvature',eta:'CAMP eta',R2:'CAMP R²',patch_id:'CAMP patch ID',area_mm2:'Patch area (mm²)',quality_status:'Fit pass/fail'};
function campTrace(m){let field=byId('campField').value,vals=m.values[field],signed=field==='avg_kg'||field==='integrated_K',quality=field==='quality_status',finite=vals.filter(Number.isFinite),lo=finite.reduce((a,b)=>Math.min(a,b),Infinity),hi=finite.reduce((a,b)=>Math.max(a,b),-Infinity);if(!Number.isFinite(lo)||!Number.isFinite(hi)){lo=0;hi=1}if(signed){let ext=Math.max(Math.abs(lo),Math.abs(hi),1e-12);lo=-ext;hi=ext}if(lo===hi){let e=Math.max(Math.abs(lo)*.01,1e-12);lo-=e;hi+=e}let trace={type:'mesh3d',x:m.vertices.map(v=>v[0]),y:m.vertices.map(v=>v[1]),z:m.vertices.map(v=>v[2]),i:m.triangles.map(t=>t[0]),j:m.triangles.map(t=>t[1]),k:m.triangles.map(t=>t[2]),intensity:vals,intensitymode:'cell',flatshading:true,lighting:{ambient:.52,diffuse:.88,specular:.18,roughness:.7},colorscale:quality?[[0,'#2ca25f'],[.499,'#2ca25f'],[.5,'#de2d26'],[1,'#de2d26']]:signed?'RdBu':'Viridis',reversescale:signed,cmin:quality?0:lo,cmax:quality?1:hi,text:m.hover.flatMap(p=>[p,p,p]),hovertemplate:'CAMP patch %{text}<br>'+CAMP_LABELS[field]+': %{intensity:.5g}<extra></extra>',colorbar:{title:{text:CAMP_LABELS[field]},...(quality?{tickvals:[0,1],ticktext:['Pass','Fail']}: {})},showscale:true};let traces=[trace];if(byId('campEdges').checked)traces.push({type:'scatter3d',mode:'lines',x:m.boundary.x,y:m.boundary.y,z:m.boundary.z,line:{color:'#18212b',width:2},hoverinfo:'skip',showlegend:false});return traces}
function planeTrace(origin,n){let norm=Math.hypot(...n)||1;n=n.map(x=>x/norm);let r=Math.abs(n[2])<.9?[0,0,1]:[0,1,0];let cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];let u=cross(n,r);let un=Math.hypot(...u);u=u.map(x=>x/un);let v=cross(n,u);let h=Math.max(...[bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]])*0.7;let q=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([a,b])=>origin.map((x,i)=>x+h*(a*u[i]+b*v[i])));return {type:'mesh3d',x:q.map(p=>p[0]),y:q.map(p=>p[1]),z:q.map(p=>p[2]),i:[0,0],j:[1,2],k:[2,3],color:'#e63946',opacity:.35,showscale:false,hoverinfo:'skip'} }
function measurePlaneTrace(){let n=vecFromAngles(),o=currentOrigin().map((x,i)=>x-n[i]*Number(byId('measureOffset').value)),t=planeTrace(o,n);t.color='#8e24aa';t.opacity=.35;return t}
function copyCamera(camera){return camera?JSON.parse(JSON.stringify(camera)):null}
function cameraFromPlot(id){
    let gd=byId(id),scene=gd?._fullLayout?.scene;
    return copyCamera(scene?._scene?.getCamera()||scene?.camera||gd?.layout?.scene?.camera||null);
}
function rememberCamera(preferred){let source=preferred||cameraFollowSource,camera=(source&&cameraFromPlot(source))||cameraFromPlot('inner')||cameraFromPlot('outer')||lastCamera;if(camera)lastCamera=camera;return copyCamera(lastCamera)}
function rememberDragMode(preferred){let gd=preferred&&byId(preferred), mode=gd?._fullLayout?.scene?.dragmode||gd?.layout?.scene?.dragmode||byId('inner')._fullLayout?.scene?.dragmode||byId('outer')._fullLayout?.scene?.dragmode;if(mode)lastDragMode=mode;return lastDragMode}
async function applyViewState(id,camera,dragMode=lastDragMode){if(!byId(id).data)return;let update={'scene.dragmode':dragMode};if(camera)update['scene.camera']=copyCamera(camera);cameraSyncBusy++;try{await Plotly.relayout(id,update)}finally{cameraSyncBusy--}}
function sameCamera(a,b){
    return a&&b&&['eye','up','center'].every(key=>
        ['x','y','z'].every(axis=>Math.abs(a[key][axis]-b[key][axis])<1e-7))
        &&a.projection?.type===b.projection?.type;
}
function followCamera(now){
    cameraFollowFrame=null;
    if(!cameraFollowSource)return;
    if(!cameraSyncBusy){
        let id=cameraFollowSource,other=id==='inner'?'outer':'inner',camera=cameraFromPlot(id);
        if(camera){
            lastCamera=camera;
            lastDragMode=rememberDragMode(id);
            if(!sameCamera(camera,cameraFromPlot(other))||byId(other)._fullLayout?.scene?.dragmode!==lastDragMode)
                applyViewState(other,camera,lastDragMode);
        }
    }
    // Turntable interpolation can continue after plotly_relayout/mouseup.
    if(cameraPointerDown||now<cameraFollowUntil)cameraFollowFrame=requestAnimationFrame(followCamera);
    else cameraFollowSource=null;
}
function queueCameraFollow(id){
    cameraFollowSource=id;
    cameraFollowUntil=performance.now()+750;
    if(cameraFollowFrame===null)cameraFollowFrame=requestAnimationFrame(followCamera);
}
function onCameraRelayout(id,event){
    if(cameraSyncBusy)return;
    let keys=Object.keys(event||{});
    // Repeated drags can emit {} even though the live camera has changed.
    if(keys.length&&!keys.some(key=>key.startsWith('scene.')))return;
    queueCameraFollow(id);
}
function installCameraSync(){
    if(cameraSyncInstalled||!byId('inner').data||!byId('outer').data)return;
    for(let id of ['inner','outer']){
        byId(id).on('plotly_relayout',event=>onCameraRelayout(id,event));
        byId(id).addEventListener('pointerdown',()=>{cameraPointerDown=true;queueCameraFollow(id)});
    }
    let release=()=>{cameraPointerDown=false;cameraFollowUntil=performance.now()+750};
    window.addEventListener('pointerup',release);
    window.addEventListener('pointercancel',release);
    window.addEventListener('blur',release);
    cameraSyncInstalled=true;
    rememberCamera();
}
// A drag changes the live WebGL camera without necessarily updating gd.layout.
// Carry that camera into every data update, including the immediate plane update.
function layout(camera=rememberCamera()){
    let dragmode=rememberDragMode();
    let spans=[0,1,2].map(i=>Math.max(bounds[i*2+1]-bounds[i*2],.01));
    let scale=Math.cbrt(spans[0]*spans[1]*spans[2]);
    let axis=(title,i)=>{
        let lo=bounds[i*2],hi=bounds[i*2+1],pad=Math.max((hi-lo)*.05,.01);
        return {title,autorange:false,range:[lo-pad,hi+pad]};
    };
    // Use the same fixed coordinate frame in both panes, independent of the plane.
    let scene={uirevision:caseData?caseData.case:'cardiac-review',aspectmode:'manual',
        aspectratio:{x:spans[0]/scale,y:spans[1]/scale,z:spans[2]/scale},dragmode,
        xaxis:axis('X',0),yaxis:axis('Y',1),zaxis:axis('Z',2)};
    if(camera)scene.camera=copyCamera(camera);
    return {uirevision:caseData?caseData.case:'cardiac-review',
        margin:{l:0,r:0,t:0,b:0},showlegend:false,scene};
}
async function renderInner(){let traces=[];if(byId('showInner').checked){if(byId('displayMode').value==='camp'&&campData)traces=campTrace(campData.inner);else if(lastPreview){traces=[meshTrace(lastPreview.inner_kept,false)];if(byId('showRemoved').checked&&lastPreview.inner_removed.n_faces)traces.push(meshTrace(lastPreview.inner_removed,false,.22,true));if(byId('showPlane').checked)traces.push(planeTrace(currentOrigin(),vecFromAngles()))}else if(caseData?.inner)traces=[meshTrace(caseData.inner,false)];if(measureResult&&byId('showMeasurePlane').checked)traces.push(measurePlaneTrace())}await Plotly.react('inner',traces,layout(),{responsive:true,displaylogo:false});installCameraSync()}
async function renderOuter(){let traces=[];if(byId('showOuter').checked){if(byId('displayMode').value==='camp'&&campData)traces=campTrace(campData.outer);else if(lastPreview){traces=[meshTrace(lastPreview.kept,true)];if(byId('showRemoved').checked&&lastPreview.removed.n_faces)traces.push(meshTrace(lastPreview.removed,false,.22,true));if(byId('showPlane').checked)traces.push(planeTrace(currentOrigin(),vecFromAngles()))}else if(caseData?.source){traces=[meshTrace(caseData.source,true)];if(byId('showPlane').checked)traces.push(planeTrace(currentOrigin(),vecFromAngles()))}if(measureResult&&byId('showMeasurePlane').checked)traces.push(measurePlaneTrace())}await Plotly.react('outer',traces,layout(),{responsive:true,displaylogo:false});installCameraSync()}
async function renderBoth(){await renderInner();await renderOuter()}
function showQuality(q){if(!q)return;let describe=s=>`${s.components===1&&s.simple_boundary_loops?'Connected':'Review connectivity'} · ${s.faces.toLocaleString()} faces · ${s.boundary_loops} boundary loops · Euler ${s.euler}`;let bad=q.inner.components!==1||q.outer.components!==1||!q.inner.simple_boundary_loops||!q.outer.simple_boundary_loops;let node=byId('quality');node.textContent=`Inner: ${describe(q.inner)} | Outer: ${describe(q.outer)} | ${q.outer_removed_pct.toFixed(1)}% of outer area removed`;node.classList.toggle('bad',bad)}
function updatePlaneGeometry(id,index,t){
    let camera=cameraFromPlot(id),view={};
    if(camera)view['scene.camera']=camera;
    return Plotly.update(id,{x:[t.x],y:[t.y],z:[t.z]},view,[index]);
}
function moveCutPlane(){if(!byId('showPlane').checked||!bounds)return;let t=planeTrace(currentOrigin(),vecFromAngles());for(let id of ['inner','outer']){let gd=byId(id);if(!gd.data?.length)continue;let i=gd.data.findIndex(trace=>trace.type==='mesh3d'&&trace.color==='#e63946');if(i>=0)updatePlaneGeometry(id,i,t)}}
function moveMeasurePlane(){if(!measureResult||!byId('showMeasurePlane').checked)return;for(let id of ['inner','outer']){let gd=byId(id);if(!gd.data?.length)continue;let i=gd.data.length-1,t=measurePlaneTrace();updatePlaneGeometry(id,i,t)}}
async function measureSection(){if(!caseData)return;let button=byId('measureButton');button.disabled=true;byId('measurement').textContent='Measuring cross-sections…';try{let r=await fetch('/api/measure',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({case:caseData.case,origin:currentOrigin(),normal:vecFromAngles(),keep_side:byId('side').value,split_method:byId('planeMode').value,measure_offset_mm:Number(byId('measureOffset').value)})}),d=await r.json();if(!r.ok)throw new Error(d.error||'Measurement failed');measureResult=true;let f=x=>Number.isFinite(x)?x.toFixed(1):'—',fmt=m=>`${f(m.area_mm2)} mm² · perimeter ${f(m.perimeter_mm)} mm · diameters ${f(m.major_diameter_mm)} × ${f(m.minor_diameter_mm)} mm`;byId('measurement').textContent=`At ${d.offset_mm.toFixed(1)} mm toward the apex — inner: ${fmt(d.inner)} | outer: ${fmt(d.outer)}`;moveMeasurePlane()}catch(e){byId('measurement').textContent='Measurement failed: '+e.message}finally{button.disabled=false}}
async function refresh(){if(!caseData)return;let request=++latestRequest,o=currentOrigin(),n=vecFromAngles();if(previewAbort)previewAbort.abort();previewAbort=new AbortController();status.classList.remove('error');status.textContent='Updating preview…';
try{let res=await fetch('/api/preview',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({case:caseData.case,origin:o,normal:n,keep_side:byId('side').value,split_method:byId('planeMode').value}),signal:previewAbort.signal});let d=await res.json();if(request!==latestRequest)return;if(!res.ok){status.textContent=d.error;status.classList.add('error');return}lastPreview=d;showQuality(d.quality);await renderBoth();status.textContent=byId('planeMode').value==='manual'?`Plane split ready · ${d.inner_kept.n_faces.toLocaleString()} inner and ${d.kept.n_faces.toLocaleString()} outer triangles`:`Preview ready · ${d.kept.n_faces.toLocaleString()} outer triangles kept · start plane: ${caseData.plane.method}`}catch(e){if(e.name!=='AbortError'&&request===latestRequest){status.textContent='Preview failed: '+e.message;status.classList.add('error')}}}
function invalidateCamp(){if(campData){campData=null;byId('campOption').disabled=true;byId('displayMode').value='cut'}}
function schedule(){invalidateCamp();byId('azv').value=byId('az').value+'°';byId('elv').value=byId('el').value+'°';byId('offv').value=byId('off').value+' mm';moveCutPlane();if(measureResult){moveMeasurePlane();byId('measurement').textContent='Cut plane moved; click Measure to update the section readings.'}clearTimeout(timer);timer=setTimeout(refresh,130)}
function updateInnerNote(){byId('innerNote').textContent=byId('applyInner').checked?'cut preview · plane will be saved on inner':'cut preview · inner saved as supplied'}
async function loadCase(name){if(!name)return;status.classList.remove('error');let mode=byId('planeMode').value;status.textContent=mode==='auto'&&sourceMode==='raw_meshes'?'Running center/normal/geodesic split for '+name+'…':'Loading '+name+'…';lastPreview=null;campData=null;measureResult=null;byId('quality').textContent='Quality details appear after the preview updates.';byId('measurement').textContent='Move the measurement plane toward the apex and measure inner and outer cross-sections.';byId('campOption').disabled=true;byId('displayMode').value='cut';let r=await fetch('/api/mesh?case='+encodeURIComponent(name)+'&split_method='+encodeURIComponent(mode));let d=await r.json();if(!r.ok){status.textContent=d.error;status.classList.add('error');return}caseData=d;bounds=d.bounds;byId('applyInner').checked=!!d.apply_cut_to_inner;byId('runCamp').disabled=false;byId('applyInner').disabled=mode==='manual';byId('side').disabled=mode==='manual';updateInnerNote();baseOrigin=d.plane.origin;normal=d.plane.normal;
for(let i=0;i<3;i++){let id=['ox','oy','oz'][i],input=byId(id),lo=bounds[i*2],hi=bounds[i*2+1],span=hi-lo;input.min=lo-span;input.max=hi+span;input.value=baseOrigin[i].toFixed(3)}
let az=Math.atan2(normal[1],normal[0])*180/Math.PI,el=Math.asin(normal[2])*180/Math.PI;byId('az').value=az;byId('el').value=el;let diag=Math.hypot(bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]);byId('off').min=-diag;byId('off').max=diag;byId('off').value=0;byId('side').value='negative';await renderBoth();refresh()}
for(let id of ['ox','oy','oz','az','el','off','side'])byId(id).addEventListener('input',()=>{if(id==='az'||id==='el'){let pivot=currentOrigin(),v=Number(byId('off').value);normal=vecFromAngles();baseOrigin=pivot.map((x,i)=>x-normal[i]*v);/* Keep the plane at the same position while its normal changes. */}else if(id==='off'){let n=vecFromAngles(),v=Number(byId('off').value);for(let i=0;i<3;i++)byId(['ox','oy','oz'][i]).value=baseOrigin[i]+n[i]*v;}else if(id==='ox'||id==='oy'||id==='oz'){baseOrigin=currentOrigin();byId('off').value=0;}schedule()});
byId('case').addEventListener('change',e=>loadCase(e.target.value));byId('planeMode').addEventListener('change',()=>{if(caseData)loadCase(caseData.case)});
function cameraFor(view){if(view==='front')return {eye:{x:2.2,y:0,z:.08},up:{x:0,y:0,z:1},center:{x:0,y:0,z:0}};if(view==='side')return {eye:{x:0,y:2.2,z:.08},up:{x:0,y:0,z:1},center:{x:0,y:0,z:0}};if(view==='top')return {eye:{x:0,y:0,z:2.2},up:{x:0,y:1,z:0},center:{x:0,y:0,z:0}};return {eye:{x:1.6,y:1.6,z:1.25},up:{x:0,y:0,z:1},center:{x:0,y:0,z:0}}}
function setView(view){let camera=cameraFor(view);lastCamera=copyCamera(camera);Promise.all(['inner','outer'].map(id=>applyViewState(id,camera)))}
function updateAppearance(){if(!caseData)return;renderBoth()}
function renderCampStatus(){let rows=campData.summary.rows,parts=rows.map(r=>`${r.stem.endsWith('_inner')?'inner':'outer'}: ${r.n_patches} patches · k fluct ${Number(r.k_fluct).toPrecision(3)} · failed QC ${Number(r.pct_failed_quality).toFixed(1)}%`);status.textContent='CAMP complete · '+parts.join(' | ')}
async function runCamp(){if(!caseData)return;let button=byId('runCamp');button.disabled=true;status.classList.remove('error');status.textContent='Saving reviewed surfaces, then running canonical CAMP on inner and outer…';for(let id of ['case','planeMode','side','off','az','el','ox','oy','oz','applyInner','save'])byId(id).disabled=true;try{let request=await fetch('/api/run_camp',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({case:caseData.case,origin:currentOrigin(),normal:vecFromAngles(),keep_side:byId('side').value,split_method:byId('planeMode').value,apply_cut_to_inner:byId('applyInner').checked,target_faces:Number(byId('targetFaces').value)})}),result=await request.json();if(!request.ok)throw new Error(result.error||'CAMP failed');status.textContent='Loading CAMP patch meshes…';let response=await fetch('/api/camp_mesh?case='+encodeURIComponent(caseData.case)),data=await response.json();if(!response.ok)throw new Error(data.error||'Could not load CAMP meshes');campData=data;byId('campOption').disabled=false;byId('displayMode').value='camp';await renderBoth();renderCampStatus()}catch(e){status.textContent='CAMP failed: '+e.message;status.classList.add('error')}finally{for(let id of ['case','planeMode','side','off','az','el','ox','oy','oz','applyInner','save'])byId(id).disabled=false;byId('side').disabled=byId('planeMode').value==='manual';byId('applyInner').disabled=byId('planeMode').value==='manual';button.disabled=false}}
byId('colorMode').addEventListener('change',updateAppearance);byId('contours').addEventListener('change',updateAppearance);byId('campField').addEventListener('change',updateAppearance);byId('campEdges').addEventListener('change',updateAppearance);byId('displayMode').addEventListener('change',()=>{if(byId('displayMode').value==='camp'&&!campData)byId('displayMode').value='cut';updateAppearance();if(campData&&byId('displayMode').value==='camp')renderCampStatus()});['showInner','showOuter','showRemoved','showPlane','showMeasurePlane'].forEach(id=>byId(id).addEventListener('change',updateAppearance));byId('measureButton').addEventListener('click',measureSection);byId('measureOffset').addEventListener('input',()=>{byId('measureOffsetValue').value=byId('measureOffset').value+' mm';if(measureResult){byId('measurement').textContent='Measurement plane moved; click Measure to update the readings.';moveMeasurePlane()}});byId('runCamp').addEventListener('click',runCamp);byId('applyInner').addEventListener('change',()=>{invalidateCamp();updateInnerNote()});document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>setView(b.dataset.view)));
byId('preview').addEventListener('click',refresh);byId('reset').addEventListener('click',()=>{if(caseData)loadCase(caseData.case)});
byId('save').addEventListener('click',async()=>{if(!caseData)return;byId('save').disabled=true;status.classList.remove('error');status.textContent='Saving reviewed surfaces…';let applyInner=byId('applyInner').checked||byId('planeMode').value==='manual',r=await fetch('/api/save',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({case:caseData.case,origin:currentOrigin(),normal:vecFromAngles(),keep_side:byId('side').value,split_method:byId('planeMode').value,apply_cut_to_inner:applyInner})}),d=await r.json();status.textContent=r.ok?(applyInner?'Saved plane-split inner and outer surfaces · ':'Saved outer cut; inner preserved · ')+d.output:d.error;if(!r.ok)status.classList.add('error');byId('save').disabled=false});
let sourceMode='pre_split';function showCases(cases,selected=''){let node=byId('case');node.replaceChildren(...cases.map(c=>new Option(c,c)));if(selected)node.value=selected}async function refreshInputs(selected=''){let r=await fetch('/api/cases'),d=await r.json();if(!r.ok)throw new Error(d.error||'Could not scan input folder');sourceMode=d.source_mode||'pre_split';byId('planeMode').querySelector('[value="manual"]').disabled=sourceMode!=='raw_meshes';showCases(d.cases,selected);return d}async function start(){try{let d=await refreshInputs();if(d.cases.length)loadCase(byId('case').value||d.cases[0]);else status.textContent='No meshes found in input/. Add an STL, OBJ, or PLY mesh there, or use Upload mesh.'}catch(err){status.textContent=err.message;status.classList.add('error')}}byId('refreshInputs').addEventListener('click',async()=>{try{let current=byId('case').value,d=await refreshInputs(current);if(d.cases.length&&!d.cases.includes(current))loadCase(d.cases[0]);status.textContent=d.cases.length+' mesh(es) found in input/'}catch(err){status.textContent=err.message;status.classList.add('error')}});byId('meshUpload').addEventListener('change',async e=>{let file=e.target.files[0];if(!file)return;let button=byId('meshUpload');button.disabled=true;status.classList.remove('error');status.textContent='Uploading and checking '+file.name+'…';try{let r=await fetch('/api/upload?filename='+encodeURIComponent(file.name),{method:'POST',headers:{'Content-Type':'application/octet-stream'},body:file}),d=await r.json();if(!r.ok)throw new Error(d.error||'Upload failed');await refreshInputs(d.case);byId('planeMode').querySelector('[value="manual"]').disabled=false;await loadCase(d.case);status.textContent='Loaded '+file.name+' · saved in input/'}catch(err){status.textContent='Upload failed: '+err.message;status.classList.add('error')}finally{button.disabled=false;e.target.value=''}});byId('case').addEventListener('change',()=>loadCase(byId('case').value));start();
</script></body></html>'''


def make_handler(split_dir: Path, output_dir: Path, input_dir: Path | None = None,
                 camp_results_dir: Path | None = None,
                 uploads_dir: Path | None = None):
    # Cache each automatically split case for this review session. In raw-mesh
    # mode, the exact project splitter is run once per source STL with no basal
    # cut; the adjustable plane is applied afterward in the viewer.
    mesh_cache = {}
    mesh_cache_lock = threading.Lock()
    input_dir = input_dir or (ROOT / "input")
    input_dir.mkdir(parents=True, exist_ok=True)
    supported_meshes = {".stl", ".obj", ".ply"}
    raw_sources = {}
    def scan_input():
        found = {}
        scan_dirs = list(dict.fromkeys((input_dir, uploads_dir or input_dir)))
        for folder in scan_dirs:
            if not folder.is_dir():
                continue
            for path in sorted(folder.iterdir()):
                if (path.is_file() and path.suffix.lower() in supported_meshes
                        and not path.name.startswith(("._", ".upload-"))):
                    found[path.stem] = path
        raw_sources.clear()
        raw_sources.update(found)
        return found
    scan_input()
    uploads_dir = uploads_dir or input_dir
    uploads_dir.mkdir(parents=True, exist_ok=True)
    upload_names = {}
    for p in uploads_dir.glob("*.stl"):
        if not p.name.startswith("._"):
            raw_sources[p.stem] = p
    camp_results_dir = camp_results_dir or (output_dir.parent / "cardiac_delta_k_gui")

    def run_normals_split(source: Path):
        report_path = split_dir / f"{source.stem}_split_report.json"
        inner_path = split_dir / f"{source.stem}_inner.stl"
        outer_path = split_dir / f"{source.stem}_outer_untrimmed.stl"
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        if report_path.is_file() and inner_path.is_file() and outer_path.is_file():
            try:
                cached = json.loads(report_path.read_text())
                if cached.get("source_sha256") == source_hash:
                    return cached["split_report"]
            except (OSError, ValueError, KeyError):
                pass
        splitter_path = Path(__file__).with_name("split_cardiac_surfaces.py")
        spec = importlib.util.spec_from_file_location("cardiac_surface_splitter", splitter_path)
        splitter = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(splitter)
        report = splitter.split_one(source, split_dir, top_cut="none", fill_holes_mm=0.0)
        report_path.write_text(json.dumps({"source_sha256": source_hash,
                                           "split_report": report}, indent=2) + "\n")
        return report

    def cached_mesh(case: str):
        with mesh_cache_lock:
            if case not in mesh_cache:
                report = None
                apply_cut_to_inner = False
                if case in raw_sources:
                    report = run_normals_split(raw_sources[case])
                    apply_cut_to_inner = True
                inner_path, outer_path = _path_pair(split_dir, case)
                inner_reader = _read_mesh(inner_path)
                outer_reader = _read_mesh(outer_path)
                mesh_cache[case] = (inner_reader, outer_reader,
                                    _preview_surface(outer_reader.GetOutput()),
                                    apply_cut_to_inner, report, inner_path, outer_path)
            return mesh_cache[case]

    class Handler(BaseHTTPRequestHandler):
        def _send(self, status: int, data: bytes, content_type: str):
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)

        def _json(self, status: int, data):
            def clean(value):
                if isinstance(value, float) and not np.isfinite(value):
                    return None
                if isinstance(value, dict):
                    return {key: clean(item) for key, item in value.items()}
                if isinstance(value, (list, tuple)):
                    return [clean(item) for item in value]
                return value
            payload = json.dumps(clean(data), separators=(",", ":"), allow_nan=False).encode()
            if "gzip" in self.headers.get("Accept-Encoding", "") and len(payload) > 4096:
                payload = gzip.compress(payload, compresslevel=1)
                self.send_response(status)
                self.send_header("Content-Type", "application/json; charset=utf-8")
                self.send_header("Content-Encoding", "gzip")
                self.send_header("Content-Length", str(len(payload)))
                self.send_header("Cache-Control", "no-store")
                self.end_headers()
                self.wfile.write(payload)
            else:
                self._send(status, payload, "application/json; charset=utf-8")

        def _body(self):
            length = int(self.headers.get("Content-Length", "0"))
            return json.loads(self.rfile.read(length))

        def _upload(self, url):
            length = int(self.headers.get("Content-Length", "0"))
            max_bytes = 256 * 1024 * 1024
            if length <= 0:
                raise ValueError("The uploaded file is empty")
            if length > max_bytes:
                raise ValueError("Uploads are limited to 256 MB")
            name = parse_qs(url.query).get("filename", [""])[0]
            suffix = Path(name).suffix.lower()
            if suffix not in {".stl", ".obj", ".ply"}:
                raise ValueError("Choose an STL, OBJ, or PLY triangle mesh")
            data = self.rfile.read(length)
            if len(data) != length:
                raise ValueError("The upload ended before the full file arrived")
            temp_path = uploads_dir / f".upload-{uuid.uuid4().hex}{suffix}"
            try:
                import trimesh
                temp_path.write_bytes(data)
                mesh = trimesh.load(str(temp_path), force="mesh", process=True)
                if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
                    raise ValueError("The file does not contain a triangle surface")
                if not np.isfinite(mesh.vertices).all():
                    raise ValueError("The mesh contains non-finite coordinates")
                if len(mesh.faces) > 2_000_000:
                    raise ValueError("Meshes are limited to 2 million triangles for interactive review")
                digest = hashlib.sha256(data).hexdigest()[:10]
                stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(name).stem).strip("._-")
                stem = (stem or "mesh")[:48]
                case = f"upload_{stem}_{digest}"
                target = uploads_dir / f"{case}.stl"
                if suffix == ".stl":
                    temp_path.replace(target)
                else:
                    exported = mesh.export(file_type="stl")
                    target.write_bytes(exported.encode() if isinstance(exported, str) else exported)
                    temp_path.unlink(missing_ok=True)
                with mesh_cache_lock:
                    raw_sources[case] = target
                    upload_names[case] = Path(name).name
                self._json(201, {"case": case, "filename": Path(name).name,
                                 "faces": int(len(mesh.faces)),
                                 "source": "input"})
            finally:
                temp_path.unlink(missing_ok=True)

        def do_GET(self):
            url = urlparse(self.path)
            if url.path == "/":
                self._send(200, PAGE.encode(), "text/html; charset=utf-8")
            elif url.path == "/plotly.min.js":
                self._send(200, get_plotlyjs().encode(), "application/javascript; charset=utf-8")
            elif url.path == "/api/cases":
                sources = scan_input()
                found = sorted(sources) if sources else _cases(split_dir)
                self._json(200, {"cases": found,
                                 "source_mode": "raw_meshes" if sources else "pre_split"})
            elif url.path == "/api/mesh":
                case = parse_qs(url.query).get("case", [""])[0]
                try:
                    split_method = parse_qs(url.query).get("split_method", ["auto"])[0]
                    if split_method not in ("manual", "auto"):
                        raise ValueError("split_method must be manual or auto")
                    if split_method == "manual":
                        if case not in raw_sources:
                            raise ValueError("Manual plane splitting requires a raw STL, OBJ, or PLY mesh in input/")
                        raw_reader = _read_mesh(raw_sources[case])
                        outer = raw_reader.GetOutput()
                        inner_payload = {"x": [], "y": [], "z": [], "i": [], "j": [], "k": [], "n_faces": 0}
                        apply_cut_to_inner, report = True, None
                        outer_path = raw_sources[case]
                    else:
                        inner_reader, outer_reader, _, apply_cut_to_inner, report, inner_path, outer_path = cached_mesh(case)
                        outer = outer_reader.GetOutput()
                        inner_payload = _poly_payload(inner_reader.GetOutput())
                    plane = _initial_plane(outer, case, outer_path,
                                           "manual" if split_method == "manual" else "auto",
                                           source_path=raw_sources.get(case))
                    inner_pts = np.asarray([
                        *[inner_reader.GetOutput().GetPoint(i)
                          for i in range(inner_reader.GetOutput().GetNumberOfPoints())]
                    ]) if split_method == "auto" else np.empty((0, 3))
                    if len(inner_pts):
                        inner_projection = (inner_pts - plane["axis_center"]) @ plane["normal"]
                        plane["axis_range"] = [
                            min(plane["axis_range"][0], float(inner_projection.min())),
                            max(plane["axis_range"][1], float(inner_projection.max())),
                        ]
                    plane["guidance"] = split_method
                    self._json(200, {
                        "case": case, "inner": inner_payload,
                        "source": _poly_payload(outer) if split_method == "manual" else None,
                        "bounds": list(outer.GetBounds()),
                        "plane": plane,
                        "apply_cut_to_inner": apply_cut_to_inner,
                        "split_report": report,
                        "split_method": split_method,
                    })
                except Exception as exc:
                    self._json(400, {"error": str(exc)})
            elif url.path == "/api/camp_mesh":
                case = parse_qs(url.query).get("case", [""])[0]
                try:
                    self._json(200, self._camp_payload(case))
                except Exception as exc:
                    self._json(400, {"error": str(exc)})
            else:
                self._json(404, {"error": "not found"})

        def do_POST(self):
            url = urlparse(self.path)
            if url.path == "/api/upload":
                try:
                    self._upload(url)
                except Exception as exc:
                    self._json(400, {"error": str(exc)})
                return
            try:
                body = self._body()
                case = str(body["case"])
                origin = np.asarray(body["origin"], dtype=float)
                normal = np.asarray(body["normal"], dtype=float)
                if origin.shape != (3,) or normal.shape != (3,) or not np.all(np.isfinite(origin)):
                    raise ValueError("Plane origin and normal must be finite 3D vectors")
                if not np.all(np.isfinite(normal)) or np.linalg.norm(normal) < 1e-8:
                    raise ValueError("Plane normal must be nonzero and finite")
                normal /= np.linalg.norm(normal)
                split_method = body.get("split_method", "auto")
                if split_method not in ("manual", "auto"):
                    raise ValueError("split_method must be manual or auto")
                keep_side = body.get("keep_side", "negative")
                if keep_side not in ("positive", "negative"):
                    raise ValueError("keep_side must be positive or negative")
                if split_method == "manual":
                    if case not in raw_sources:
                        raise ValueError("Manual plane splitting requires a raw STL, OBJ, or PLY mesh in input/")
                    inner_piece, outer_piece = _manual_plane_split(
                        raw_sources[case], origin, normal, keep_side)
                    inner_reader = _reader_from_trimesh(inner_piece)
                    reader = _reader_from_trimesh(outer_piece)
                    preview_poly = reader.GetOutput()
                    inner_path = outer_path = raw_sources[case]
                    default_apply_inner, split_report = True, None
                else:
                    inner_reader, reader, preview_poly, default_apply_inner, split_report, inner_path, outer_path = cached_mesh(case)
                if split_method == "manual":
                    clip = reader
                    inner_clip = inner_reader
                    apply_inner = True
                    preview_reader = reader
                elif url.path == "/api/preview":
                    # Preview the same plane cut on both surfaces.
                    preview_reader = _reader_from_poly(preview_poly)
                    clip = _clip(preview_reader, origin, normal, keep_side)
                    inner_clip = _clip(inner_reader, origin, normal, keep_side)
                else:
                    clip = _clip(reader, origin, normal, keep_side)
                    apply_inner = bool(body.get("apply_cut_to_inner", default_apply_inner))
                    inner_clip = _clip(inner_reader, origin, normal, keep_side) if apply_inner else None
                poly = clip.GetOutput()
                if poly.GetNumberOfCells() == 0:
                    raise ValueError("The selected plane removed the entire outer surface")
                if url.path == "/api/preview":
                    inner_poly = inner_clip.GetOutput()
                    if split_method == "manual":
                        from vtkmodules.vtkCommonDataModel import vtkPolyData
                        empty = vtkPolyData()
                        removed_poly = inner_removed_poly = empty
                    else:
                        removed_clip = _clip(preview_reader, origin, normal,
                                             "positive" if keep_side == "negative" else "negative")
                        inner_removed_clip = _clip(inner_reader, origin, normal,
                                                   "positive" if keep_side == "negative" else "negative")
                        removed_poly = removed_clip.GetOutput()
                        inner_removed_poly = inner_removed_clip.GetOutput()
                    inner_mesh = _poly_mesh(inner_poly)
                    outer_mesh = _poly_mesh(poly)
                    original_area = max(float(_poly_mesh(preview_poly).area), 1e-12)
                    removed_mesh = _poly_mesh(removed_poly)
                    self._json(200, {
                        "kept": _poly_payload(poly),
                        "removed": _poly_payload(removed_poly),
                        "inner_kept": _poly_payload(inner_poly),
                        "inner_removed": _poly_payload(inner_removed_poly),
                        "quality": {
                            "inner": _surface_diagnostics(inner_mesh),
                            "outer": _surface_diagnostics(outer_mesh),
                            "outer_removed_pct": (100.0 * removed_mesh.area / original_area
                                                  if len(removed_mesh.faces) else 0.0),
                        },
                    })
                elif url.path == "/api/measure":
                    measure_offset = float(body.get("measure_offset_mm", 0.0))
                    if not np.isfinite(measure_offset):
                        raise ValueError("Measurement offset must be finite")
                    if split_method == "manual":
                        inner_mesh, outer_mesh = inner_piece, outer_piece
                    else:
                        import trimesh
                        outer_mesh_poly = preview_poly
                        inner_mesh, outer_mesh = _poly_mesh(inner_reader.GetOutput()), \
                            _poly_mesh(outer_mesh_poly)
                        clip_normal = normal if keep_side == "negative" else -normal
                        inner_vertices, inner_faces = trimesh.intersections.slice_faces_plane(
                            inner_mesh.vertices, inner_mesh.faces, -clip_normal, origin)[:2]
                        outer_vertices, outer_faces = trimesh.intersections.slice_faces_plane(
                            outer_mesh.vertices, outer_mesh.faces, -clip_normal, origin)[:2]
                        inner_mesh = trimesh.Trimesh(vertices=inner_vertices,
                                                     faces=inner_faces, process=True)
                        outer_mesh = trimesh.Trimesh(vertices=outer_vertices,
                                                     faces=outer_faces, process=True)
                    measure_origin = origin - normal * measure_offset
                    self._json(200, {
                        "offset_mm": measure_offset,
                        "inner": _section_measurement(inner_mesh, measure_origin, normal),
                        "outer": _section_measurement(outer_mesh, measure_origin, normal),
                    })
                elif url.path in ("/api/save", "/api/run_camp"):
                    if apply_inner and inner_clip.GetOutput().GetNumberOfCells() == 0:
                        raise ValueError("The selected plane removed the entire inner surface")
                    inner_quality_poly = (inner_clip.GetOutput() if apply_inner
                                          else inner_reader.GetOutput())
                    quality = {"inner": _surface_diagnostics(_poly_mesh(inner_quality_poly)),
                               "outer": _surface_diagnostics(_poly_mesh(poly))}
                    folder = self._save(case, inner_path, outer_path, inner_clip, clip,
                                        origin, normal, keep_side, apply_inner, split_report,
                                        split_method, quality)
                    if url.path == "/api/save":
                        self._json(200, {"output": str(folder.resolve()), "quality": quality})
                    else:
                        target_faces = int(body.get("target_faces", 50_000))
                        if target_faces < 5_000 or target_faces > 200_000:
                            raise ValueError("CAMP target faces must be between 5000 and 200000")
                        summary = self._run_camp(case, target_faces)
                        self._json(200, {"output": str(folder.resolve()),
                                         "target_faces": target_faces, "summary": summary})
                else:
                    self._json(404, {"error": "not found"})
            except Exception as exc:
                self._json(400, {"error": str(exc)})

        def _save(self, case, inner_path, outer_path, inner_clip, outer_clip, origin, normal,
                  keep_side, apply_inner, split_report, plane_mode, quality):
            from vtkmodules.vtkIOGeometry import vtkSTLWriter
            folder = output_dir / case
            folder.mkdir(parents=True, exist_ok=True)
            inner_out = folder / f"{case}_inner.stl"
            outer_out = folder / f"{case}_outer.stl"
            if apply_inner:
                writer = vtkSTLWriter()
                writer.SetFileName(str(inner_out))
                writer.SetInputConnection(inner_clip.GetOutputPort())
                if writer.Write() != 1:
                    raise OSError("VTK could not write the clipped inner STL")
            else:
                inner_out.write_bytes(inner_path.read_bytes())
            writer = vtkSTLWriter()
            writer.SetFileName(str(outer_out))
            writer.SetInputConnection(outer_clip.GetOutputPort())
            if writer.Write() != 1:
                raise OSError("VTK could not write the clipped outer STL")
            metadata = {
                "source_inner_stl": str(inner_path.resolve()),
                "source_outer_stl": str(outer_path.resolve()),
                "source_inner_sha256": hashlib.sha256(inner_path.read_bytes()).hexdigest(),
                "source_outer_sha256": hashlib.sha256(outer_path.read_bytes()).hexdigest(),
                "inner_stl": inner_out.name, "outer_stl": outer_out.name,
                "inner_cut_applied": apply_inner,
                "plane_origin": origin.tolist(), "plane_normal": normal.tolist(),
                "kept_side": keep_side,
                "plane_guidance": "manual clinician placement" if plane_mode == "manual" else "automatic basal-rim suggestion",
                "saved_at_utc": datetime.now(timezone.utc).isoformat(),
                "split_method": ("clinician plane cut separated disconnected inner/outer walls"
                                 if plane_mode == "manual" else
                                 "center-of-mass/normal seeds with geodesic partition; clinician-adjusted planar trim"),
                "automatic_split_report": split_report,
                "post_cut_quality": quality,
                "uploaded_source": Path(inner_path).resolve().is_relative_to(uploads_dir.resolve()),
                "original_upload_filename": upload_names.get(case),
                "camp_source": "camp/camp_updated_2026-08-26/campsrc.py",
            }
            (folder / "review.json").write_text(json.dumps(metadata, indent=2) + "\n")
            return folder

        def _run_camp(self, case, target_faces):
            runner_path = Path(__file__).with_name("run_cardiac_delta_k.py")
            spec = importlib.util.spec_from_file_location("cardiac_delta_k_runner", runner_path)
            runner = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(runner)
            run_stl, camp_source = runner.load_camp()
            camp_results_dir.mkdir(parents=True, exist_ok=True)
            rows = []
            for side in ("inner", "outer"):
                stl = output_dir / case / f"{case}_{side}.stl"
                row = run_stl(stl, camp_results_dir, target_faces=target_faces,
                              trim_frac=0.01, save_html=False, save_png=False,
                              verbose=False, campsrc_module=camp_source)
                rows.append(row)
                if row.get("status") != "ok":
                    raise RuntimeError(f"CAMP failed for {case} {side}: {row.get('error', 'unknown error')}")
            summary_path = camp_results_dir / f"{case}_summary.json"
            summary_path.write_text(json.dumps({"case": case, "target_faces": target_faces,
                                                "trim_frac": 0.01, "rows": rows}, indent=2) + "\n")
            return rows

        def _camp_payload(self, case):
            summary_path = camp_results_dir / f"{case}_summary.json"
            if not summary_path.is_file():
                raise FileNotFoundError(f"Run CAMP on {case} first")
            summary = json.loads(summary_path.read_text())
            viewer_path = Path(__file__).with_name("view_cardiac_camp.py")
            spec = importlib.util.spec_from_file_location("cardiac_camp_viewer", viewer_path)
            viewer = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(viewer)
            quality = {row["stem"]: row for row in summary["rows"]}
            surfaces = {}
            for side in ("inner", "outer"):
                stem = f"{case}_{side}"
                npz_path = camp_results_dir / stem / f"{stem}_camp.npz"
                surfaces[side] = viewer.mesh_payload(npz_path, stem, quality)
            return {"case": case, "inner": surfaces["inner"],
                    "outer": surfaces["outer"], "summary": summary}

        def log_message(self, fmt, *args):
            print(f"[review-web] {self.address_string()} - {fmt % args}")

    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "--dataset", dest="input_dir", type=Path,
                        default=ROOT / "input",
                        help="folder scanned for STL, OBJ, and PLY meshes (default: ./input)")
    parser.add_argument("--split-dir", type=Path, default=ROOT / "results" / "gui_normals_split",
                        help="working directory for automatic split outputs")
    parser.add_argument("--pre-split-dir", type=Path,
                        help="review existing *_inner.stl and *_outer_untrimmed.stl pairs instead")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "manual_surface_split")
    parser.add_argument("--camp-results", type=Path,
                        default=ROOT / "results" / "cardiac_delta_k_gui",
                        help="CAMP output directory for per-mesh GUI runs")
    parser.add_argument("--uploads", type=Path, default=None,
                        help="optional upload storage override (default: the input folder)")
    parser.add_argument("--host", default="127.0.0.1", help="bind to localhost by default")
    parser.add_argument("--port", type=int, default=8766)
    parser.add_argument("--no-browser", action="store_true", help="print URL without opening a browser")
    args = parser.parse_args()
    input_dir = args.input_dir
    split_dir = args.pre_split_dir or args.split_dir
    if not any(p.is_file() and p.suffix.lower() in {".stl", ".obj", ".ply"}
               for p in input_dir.iterdir()) if input_dir.is_dir() else True:
        print(f"No STL, OBJ, or PLY meshes found in {input_dir}; browser uploads remain available.")
    if args.pre_split_dir and not args.pre_split_dir.is_dir():
        raise SystemExit(f"Split directory not found: {args.pre_split_dir}")
    split_dir.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer((args.host, args.port),
                                 make_handler(split_dir, args.output, input_dir,
                                              args.camp_results, args.uploads))
    url = f"http://127.0.0.1:{args.port}/"
    print(f"Cardiac split review app running at {url}")
    print("Press Ctrl+C to stop the server.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("Stopping review app.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
