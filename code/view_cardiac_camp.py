#!/usr/bin/env python3
"""Local interactive 3D viewer for the cardiac CAMP results."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import threading
import webbrowser
from collections import defaultdict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np
from plotly.offline import get_plotlyjs
import plotly.graph_objects as go
from plotly.subplots import make_subplots


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "cardiac_delta_k"
LABELS = {
    "quality_status": "CAMP fit QC (green pass / red fail)",
    "patch_id": "CAMP patch ID",
    "avg_kg": "Mean Gaussian curvature",
    "integrated_K": "Integrated Gaussian curvature",
    "eta": "CAMP eta",
    "R2": "CAMP R²",
    "area_mm2": "Patch area (mm²)",
}
SIGNED = {"avg_kg", "integrated_K"}
QUALITY_SCALE = [[0.0, "#2ca25f"], [0.499, "#2ca25f"],
                 [0.5, "#de2d26"], [1.0, "#de2d26"]]
QUALITY_LIMITS = (0.0, 1.0)


def discover_cases(results: Path) -> dict[str, Path]:
    cases = {}
    for path in sorted(results.rglob("*_camp.npz")):
        cases[path.stem.removesuffix("_camp")] = path
    return cases


def load_quality(results: Path) -> dict[str, dict[str, str]]:
    path = results / "summary.csv"
    if not path.exists():
        return {}
    with path.open(newline="") as stream:
        return {row.get("stem", ""): row for row in csv.DictReader(stream)}


def edge_boundaries(vertices: np.ndarray, faces: np.ndarray,
                    labels: np.ndarray) -> tuple[list, list, list]:
    incidence: dict[tuple[int, int], list[int]] = defaultdict(list)
    for tri, patch in zip(faces, labels):
        a, b, c = (int(v) for v in tri)
        incidence[tuple(sorted((a, b)))].append(int(patch))
        incidence[tuple(sorted((b, c)))].append(int(patch))
        incidence[tuple(sorted((c, a)))].append(int(patch))
    xs, ys, zs = [], [], []
    for (a, b), patches in incidence.items():
        if len(patches) == 1 or len(set(patches)) > 1:
            xs.extend((float(vertices[a, 0]), float(vertices[b, 0]), None))
            ys.extend((float(vertices[a, 1]), float(vertices[b, 1]), None))
            zs.extend((float(vertices[a, 2]), float(vertices[b, 2]), None))
    return xs, ys, zs


def mesh_payload(path: Path, case: str, quality: dict[str, dict[str, str]]) -> dict:
    with np.load(path, allow_pickle=False) as data:
        vertices = np.asarray(data["V"], dtype=np.float32)
        faces = np.asarray(data["F"], dtype=np.int32)
        labels = np.asarray(data["face_labels"], dtype=np.int32)
        metrics = {
            "avg_kg": np.asarray(data["avg_kg"], dtype=np.float32),
            "integrated_K": np.asarray(data["K_integrated_patch"], dtype=np.float32),
            "eta": np.asarray(data["eta_patch"], dtype=np.float32),
            "R2": np.asarray(data["R2_patch"], dtype=np.float32),
            "area_mm2": np.asarray(data["patch_area"], dtype=np.float32),
        }
        kept = np.asarray(data["keep"], dtype=bool)

    if len(faces) != len(labels) or (labels.min(initial=0) < 0) or (labels.max(initial=-1) >= len(kept)):
        raise ValueError(f"Invalid face/patch labels in {path}")
    face_values = {"patch_id": labels.astype(float)}
    for name, values in metrics.items():
        face_values[name] = values[labels].astype(float)
    patch_pass = (metrics["R2"] >= 0.90) & (metrics["eta"] <= 0.075)
    face_values["quality_status"] = (~patch_pass[labels]).astype(float)
    bx, by, bz = edge_boundaries(vertices, faces, labels)
    quality_row = quality.get(case, {})
    return {
        "case": case,
        "vertices": vertices[faces].reshape(-1, 3).tolist(),
        "triangles": np.arange(faces.size, dtype=np.int32).reshape(-1, 3).tolist(),
        "patch_ids": labels.tolist(),
        "values": {name: val.tolist() for name, val in face_values.items()},
        "hover": labels.tolist(),
        "boundary": {"x": bx, "y": by, "z": bz},
        "n_faces": int(len(faces)),
        "n_patches": int(len(kept)),
        "n_kept": int(kept.sum()),
        "pct_failed_quality": quality_row.get("pct_failed_quality", ""),
        "pct_at_min_size": quality_row.get("pct_at_min_size", ""),
    }


def make_figure(payload: dict, field: str = "integrated_K", boundaries: bool = True):
    """Create a Plotly figure with one scalar value per CAMP face."""
    vertices = np.asarray(payload["vertices"], dtype=float)
    triangles = np.asarray(payload["triangles"], dtype=int)
    values = np.asarray(payload["values"][field], dtype=float)
    lo, hi = QUALITY_LIMITS if field == "quality_status" else (float(values.min()), float(values.max()))
    if field in SIGNED:
        extent = max(abs(lo), abs(hi), np.finfo(float).eps)
        lo, hi = -extent, extent
    if hi == lo:
        eps = max(abs(lo) * .01, 1e-12)
        lo, hi = lo - eps, hi + eps
    fig = go.Figure(go.Mesh3d(
        x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
        i=triangles[:, 0], j=triangles[:, 1], k=triangles[:, 2],
        intensity=values, intensitymode="cell", flatshading=True,
        colorscale=QUALITY_SCALE if field == "quality_status" else ("RdBu" if field in SIGNED else "Viridis"),
        reversescale=field in SIGNED, cmin=lo, cmax=hi,
        text=np.repeat(payload["hover"], 3),
        hovertemplate=f"CAMP patch %{{text}}<br>{LABELS[field]}: %{{intensity:.6g}}<extra></extra>",
        colorbar={"title": LABELS[field],
                  **({"tickvals": [0, 1], "ticktext": ["Pass", "Fail"]}
                     if field == "quality_status" else {})}, name="CAMP mesh"))
    if boundaries:
        line = payload["boundary"]
        fig.add_trace(go.Scatter3d(
            x=line["x"], y=line["y"], z=line["z"], mode="lines",
            line={"color": "#111820", "width": 2}, hoverinfo="skip",
            showlegend=False))
    buttons = []
    for name, label in LABELS.items():
        vals = np.asarray(payload["values"][name], dtype=float)
        low, high = QUALITY_LIMITS if name == "quality_status" else (float(vals.min()), float(vals.max()))
        if name in SIGNED:
            extent = max(abs(low), abs(high), np.finfo(float).eps)
            low, high = -extent, extent
        if low == high:
            eps = max(abs(low) * .01, 1e-12)
            low, high = low - eps, high + eps
        buttons.append({"label": label, "method": "update", "args": [
            {"intensity": [vals], "colorscale": [QUALITY_SCALE if name == "quality_status" else ("RdBu" if name in SIGNED else "Viridis")],
             "reversescale": [name in SIGNED], "cmin": [low], "cmax": [high],
             "colorbar.title.text": [label],
             "colorbar.tickvals": [[0, 1] if name == "quality_status" else None],
             "colorbar.ticktext": [["Pass", "Fail"] if name == "quality_status" else None],
             "hovertemplate": [f"CAMP patch %{{text}}<br>{label}: %{{intensity:.6g}}<extra></extra>"]},
            {"title.text": f"{payload['case']} · {label}"}, [0]]})
    fig.update_layout(
        title=f"{payload['case']} · {LABELS[field]}", height=800,
        margin={"l": 0, "r": 0, "b": 0, "t": 95},
        scene={"aspectmode": "data", "xaxis": {"visible": False},
               "yaxis": {"visible": False}, "zaxis": {"visible": False}},
        showlegend=False,
        updatemenus=[{"buttons": buttons, "x": 0, "y": 1.13,
                      "xanchor": "left", "direction": "down"}],
        annotations=[{"text": "Rotate: drag · Zoom: scroll · Hover: CAMP patch and values · Color field: dropdown",
                      "x": 0, "y": -0.02, "xref": "paper", "yref": "paper",
                      "showarrow": False}])
    return fig


def make_pair_figure(inner: dict, outer: dict, field: str = "integrated_K",
                     boundaries: bool = True):
    """Show matched inner/outer CAMP surfaces in linked side-by-side panels."""
    fig = make_subplots(
        rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
        subplot_titles=("Inner surface", "Outer surface"), horizontal_spacing=0.02)
    mesh_payloads = (inner, outer)
    mesh_indices = []
    edge_indices = []
    for col, payload in enumerate(mesh_payloads, start=1):
        vertices = np.asarray(payload["vertices"], dtype=float)
        triangles = np.asarray(payload["triangles"], dtype=int)
        vals = np.asarray(payload["values"][field], dtype=float)
        lo, hi = QUALITY_LIMITS if field == "quality_status" else (float(vals.min()), float(vals.max()))
        if field in SIGNED:
            extent = max(abs(lo), abs(hi), np.finfo(float).eps)
            lo, hi = -extent, extent
        if lo == hi:
            eps = max(abs(lo) * .01, 1e-12)
            lo, hi = lo - eps, hi + eps
        mesh_indices.append(len(fig.data))
        fig.add_trace(go.Mesh3d(
            x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
            i=triangles[:, 0], j=triangles[:, 1], k=triangles[:, 2],
            intensity=vals, intensitymode="cell", flatshading=True,
            colorscale=QUALITY_SCALE if field == "quality_status" else ("RdBu" if field in SIGNED else "Viridis"),
            reversescale=field in SIGNED, cmin=lo, cmax=hi,
            text=np.repeat(payload["hover"], 3),
            hovertemplate=f"CAMP patch %{{text}}<br>{LABELS[field]}: %{{intensity:.6g}}<extra></extra>",
            colorbar={"title": LABELS[field], "x": 0.46 if col == 1 else 1.02,
                      **({"tickvals": [0, 1], "ticktext": ["Pass", "Fail"]}
                         if field == "quality_status" else {})},
            name=f"{payload['case']} CAMP mesh"), row=1, col=col)
        if boundaries:
            line = payload["boundary"]
            edge_indices.append(len(fig.data))
            fig.add_trace(go.Scatter3d(
                x=line["x"], y=line["y"], z=line["z"], mode="lines",
                line={"color": "#111820", "width": 2}, hoverinfo="skip",
                showlegend=False, name="CAMP boundaries"), row=1, col=col)

    buttons = []
    for name, label in LABELS.items():
        updates = {"intensity": [], "colorscale": [], "reversescale": [],
                   "cmin": [], "cmax": [], "colorbar.title.text": [],
                   "hovertemplate": []}
        title = f"{inner['case'].removesuffix('_inner')} · {label}"
        for payload in mesh_payloads:
            vals = np.asarray(payload["values"][name], dtype=float)
            lo, hi = QUALITY_LIMITS if name == "quality_status" else (float(vals.min()), float(vals.max()))
            if name in SIGNED:
                extent = max(abs(lo), abs(hi), np.finfo(float).eps)
                lo, hi = -extent, extent
            if lo == hi:
                eps = max(abs(lo) * .01, 1e-12)
                lo, hi = lo - eps, hi + eps
            updates["intensity"].append(vals)
            updates["colorscale"].append(QUALITY_SCALE if name == "quality_status" else ("RdBu" if name in SIGNED else "Viridis"))
            updates["reversescale"].append(name in SIGNED)
            updates["cmin"].append(lo)
            updates["cmax"].append(hi)
            updates["colorbar.title.text"].append(label)
            updates.setdefault("colorbar.tickvals", []).append([0, 1] if name == "quality_status" else None)
            updates.setdefault("colorbar.ticktext", []).append(["Pass", "Fail"] if name == "quality_status" else None)
            updates["hovertemplate"].append(
                f"CAMP patch %{{text}}<br>{label}: %{{intensity:.6g}}<extra></extra>")
        buttons.append({"label": label, "method": "update", "args": [
            updates, {"title.text": title}, mesh_indices]})

    base_name = inner["case"].removesuffix("_inner")
    fig.update_layout(
        title=f"{base_name} · {LABELS[field]}", height=820,
        margin={"l": 0, "r": 0, "b": 0, "t": 110}, showlegend=False,
        scene={"aspectmode": "data", "xaxis": {"visible": False},
               "yaxis": {"visible": False}, "zaxis": {"visible": False}},
        scene2={"aspectmode": "data", "xaxis": {"visible": False},
                "yaxis": {"visible": False}, "zaxis": {"visible": False}},
        updatemenus=[{"buttons": buttons, "x": 0, "y": 1.12,
                      "xanchor": "left", "direction": "down"}],
        annotations=[{"text": "Rotate each surface by dragging · Zoom by scrolling · Hover for patch ID and value",
                      "x": 0.5, "y": -0.02, "xref": "paper", "yref": "paper",
                      "showarrow": False}])
    return fig


PAGE = r"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cardiac CAMP patch viewer</title>
<script>__PLOTLY__</script>
<style>
body{font:15px system-ui,sans-serif;margin:0;background:#f5f7fa;color:#17212b}
header{padding:14px 22px;background:#142638;color:white}h1{font-size:20px;margin:0 0 4px}
.sub{font-size:13px;color:#cbd5df}.controls{display:flex;gap:14px;align-items:center;flex-wrap:wrap;padding:14px 22px;background:white;border-bottom:1px solid #d7dfe7}
label{display:flex;gap:7px;align-items:center;font-size:13px}select{font:inherit;padding:7px;max-width:360px}
#status{padding:8px 22px;color:#465766;font-size:13px}#plot{height:calc(100vh - 168px);min-height:480px}
.hint{margin-left:auto;color:#526272;font-size:12px}
</style></head><body>
<header><h1>Cardiac CAMP patch viewer</h1><div class="sub">Rotate by dragging · zoom by scrolling · hover for CAMP patch and measurements</div></header>
<div class="controls">
 <label>Mesh <select id="case"></select></label>
 <label>Color by <select id="field"></select></label>
 <label><input id="edges" type="checkbox" checked> Show patch boundaries</label>
 <span class="hint">One mesh is loaded at a time.</span>
</div><div id="status">Loading available surfaces…</div><div id="plot"></div>
<script>
const FIELD_LABELS={quality_status:'CAMP fit QC (green pass / red fail)',patch_id:'CAMP patch ID',avg_kg:'Mean Gaussian curvature',integrated_K:'Integrated Gaussian curvature',eta:'CAMP eta',R2:'CAMP R²',area_mm2:'Patch area (mm²)'};
const SIGNED=new Set(['avg_kg','integrated_K']);
const QUALITY_SCALE=[[0,'#2ca25f'],[.499,'#2ca25f'],[.5,'#de2d26'],[1,'#de2d26']];
const caseSel=document.getElementById('case'),fieldSel=document.getElementById('field'),edgeSel=document.getElementById('edges'),statusEl=document.getElementById('status');
let mesh=null, requestId=0;
function limits(values,field){if(field==='quality_status')return [0,1];let finite=values.filter(Number.isFinite),lo=Math.min(...finite),hi=Math.max(...finite);if(SIGNED.has(field)){let a=Math.max(Math.abs(lo),Math.abs(hi),1e-12);lo=-a;hi=a;}if(lo===hi){let d=Math.max(Math.abs(lo)*.01,1e-12);lo-=d;hi+=d;}return [lo,hi];}
function render(){if(!mesh)return;const field=fieldSel.value,vals=mesh.values[field],[lo,hi]=limits(vals,field),patchText=mesh.patch_ids.flatMap(p=>[p,p,p]);
 const data=[{type:'mesh3d',x:mesh.vertices.map(v=>v[0]),y:mesh.vertices.map(v=>v[1]),z:mesh.vertices.map(v=>v[2]),i:mesh.triangles.map(t=>t[0]),j:mesh.triangles.map(t=>t[1]),k:mesh.triangles.map(t=>t[2]),intensity:vals,intensitymode:'cell',flatshading:true,colorscale:field==='quality_status'?QUALITY_SCALE:(SIGNED.has(field)?'RdBu':'Viridis'),reversescale:SIGNED.has(field),cmin:lo,cmax:hi,text:patchText,hovertemplate:'CAMP patch %{text}<br>'+FIELD_LABELS[field]+': %{intensity:.6g}<extra></extra>',colorbar:{title:FIELD_LABELS[field],...(field==='quality_status'?{tickvals:[0,1],ticktext:['Pass','Fail']}: {})},name:'CAMP mesh'}];
 if(edgeSel.checked)data.push({type:'scatter3d',mode:'lines',x:mesh.boundary.x,y:mesh.boundary.y,z:mesh.boundary.z,line:{color:'#111820',width:2},hoverinfo:'skip',showlegend:false});
 Plotly.react('plot',data,{title:{text:mesh.case+' · '+FIELD_LABELS[field],x:.02},margin:{l:0,r:0,b:0,t:42},scene:{aspectmode:'data',xaxis:{visible:false},yaxis:{visible:false},zaxis:{visible:false}},showlegend:false}, {responsive:true,displaylogo:false});
 const q=mesh.pct_failed_quality===''?'':' · failed quality '+Number(mesh.pct_failed_quality).toFixed(1)+'% · at minimum size '+Number(mesh.pct_at_min_size).toFixed(1)+'%';
 statusEl.textContent=`${mesh.n_faces.toLocaleString()} faces · ${mesh.n_patches} patches (${mesh.n_kept} kept)${q}`;
}
async function loadMesh(){const id=++requestId;statusEl.textContent='Loading '+caseSel.value+'…';try{const r=await fetch('/api/mesh?case='+encodeURIComponent(caseSel.value));if(!r.ok)throw new Error(await r.text());const m=await r.json();if(id!==requestId)return;mesh=m;render();}catch(e){statusEl.textContent='Could not load mesh: '+e;}}
async function init(){const r=await fetch('/api/cases');const cases=await r.json();for(const name of cases){const o=document.createElement('option');o.value=name;o.textContent=name;caseSel.append(o);}for(const [key,label]of Object.entries(FIELD_LABELS)){const o=document.createElement('option');o.value=key;o.textContent=label;fieldSel.append(o);}fieldSel.value='quality_status';caseSel.addEventListener('change',loadMesh);fieldSel.addEventListener('change',render);edgeSel.addEventListener('change',render);await loadMesh();}
init().catch(e=>statusEl.textContent='Viewer setup failed: '+e);
</script></body></html>"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--case", help="surface name, e.g. CG10_1_10_MY_inner")
    parser.add_argument("--side-by-side", action="store_true",
                        help="export matched inner and outer surfaces together")
    parser.add_argument("--field", choices=tuple(LABELS), default="integrated_K",
                        help="initial color field for exported HTML (default: integrated_K)")
    parser.add_argument("--export", type=Path,
                        help="write a standalone offline HTML viewer for one surface and exit")
    parser.add_argument("--no-open", action="store_true", help="do not open a browser automatically")
    args = parser.parse_args()
    results = args.results.resolve()
    cases = discover_cases(results)
    if not cases:
        raise SystemExit(f"No CAMP outputs found below {results}; run code/run_cardiac_delta_k.py first.")
    quality = load_quality(results)
    if args.export:
        case = args.case or next(iter(cases))
        if args.side_by_side:
            base = case.removesuffix("_inner").removesuffix("_outer")
            inner_case, outer_case = base + "_inner", base + "_outer"
            if inner_case not in cases or outer_case not in cases:
                raise SystemExit(f"Could not find both surfaces for {base!r} in {results}")
            pair = [mesh_payload(cases[name], name, quality)
                    for name in (inner_case, outer_case)]
            target = args.export.resolve()
            target.parent.mkdir(parents=True, exist_ok=True)
            make_pair_figure(*pair, field=args.field).write_html(target, include_plotlyjs=True, full_html=True)
            print(f"Wrote paired offline viewer: {target} ({target.stat().st_size / 1e6:.1f} MB)")
            return
        if case not in cases:
            raise SystemExit(f"Unknown surface {case!r}. Available surfaces: {', '.join(cases)}")
        payload = mesh_payload(cases[case], case, quality)
        target = args.export.resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        make_figure(payload, field=args.field).write_html(target, include_plotlyjs=True, full_html=True)
        print(f"Wrote offline viewer: {target} ({target.stat().st_size / 1e6:.1f} MB)")
        return
    page = PAGE.replace("__PLOTLY__", get_plotlyjs())

    class Handler(BaseHTTPRequestHandler):
        def send_json(self, payload) -> None:
            raw = json.dumps(payload, separators=(",", ":"), allow_nan=False).encode()
            zipped = gzip.compress(raw)
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Encoding", "gzip")
            self.send_header("Content-Length", str(len(zipped)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(zipped)

        def do_GET(self):
            parsed = urlparse(self.path)
            if parsed.path == "/":
                content = page.encode()
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(content)))
                self.end_headers()
                self.wfile.write(content)
            elif parsed.path == "/api/cases":
                self.send_json(list(cases))
            elif parsed.path == "/api/mesh":
                case = parse_qs(parsed.query).get("case", [""])[0]
                if case not in cases:
                    self.send_error(404, "Unknown CAMP surface")
                    return
                try:
                    self.send_json(mesh_payload(cases[case], case, quality))
                except Exception as exc:  # Return an actionable local viewer error.
                    self.send_error(500, str(exc))
            else:
                self.send_error(404)

        def log_message(self, format, *values):
            print("[viewer] " + format % values)

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    url_host = "localhost" if args.host in ("127.0.0.1", "0.0.0.0") else args.host
    url = f"http://{url_host}:{args.port}/"
    print(f"Serving {len(cases)} CAMP surfaces at {url}", flush=True)
    print("Press Ctrl+C to stop the viewer.", flush=True)
    if not args.no_open:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("Stopping viewer.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
