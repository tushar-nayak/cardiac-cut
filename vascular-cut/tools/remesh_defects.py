"""Every mesh defect the REMESHED surface carries, at each density.

NOT tools/mesh_defects.py, which inventories the SOURCE surfaces in a folder cheaply
and without remeshing. This one is the expensive counterpart: it remeshes each case to
several densities and scans what CAMP actually measures, because a defect that changes
count when the mesh is rebuilt is noise in delta K, while a fixed one is at worst an
offset.

Folds were the one we knew about. This looks for the rest, on the mesh CAMP actually
measures -- after acvd_uniform_remesh, not on the source. A defect that is a fixed
property of the geometry is at worst a constant offset; a defect whose COUNT changes
when the mesh is rebuilt is noise injected straight into delta K, because it changes
which patches carry the number.

    folded edges      dihedral past 90 deg: the surface doubles back, and a quadratic
                      fitted across it returns a curvature that is not one
    slivers           min interior angle under 5 deg: near-degenerate triangles whose
                      normals are numerically unstable
    tiny faces        area under 1% of the median: same problem, different cause
    sharp edges       dihedral over 60 deg but under 90: not folds, but creases that
                      no single quadratic can span
    non-manifold      an edge shared by more than two faces
    winding           inconsistent face orientation, which flips normal signs
    duplicate faces   the same triangle twice
    unreferenced      vertices no face uses

    python tools/remesh_defects.py <shard> <nshards>
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tools"))
import trimesh                                                     # noqa: E402
import shape_metrics as SM                                         # noqa: E402
from cutclean import acvd_uniform_remesh                           # noqa: E402
from camp_cohort import repair_folds                               # noqa: E402

SHARD = int(sys.argv[1]) if len(sys.argv) > 1 else 0
NSH = int(sys.argv[2]) if len(sys.argv) > 2 else 1
DENS = (100_000, 200_000, 400_000)


def scan(m):
    ang = np.degrees(np.asarray(m.face_adjacency_angles))
    tri = np.asarray(m.vertices)[np.asarray(m.faces)]
    e0 = np.linalg.norm(tri[:, 1] - tri[:, 0], axis=1)
    e1 = np.linalg.norm(tri[:, 2] - tri[:, 1], axis=1)
    e2 = np.linalg.norm(tri[:, 0] - tri[:, 2], axis=1)
    E = np.sort(np.stack([e0, e1, e2], axis=1), axis=1)
    # smallest interior angle, by the law of cosines on the shortest side
    num = E[:, 1] ** 2 + E[:, 2] ** 2 - E[:, 0] ** 2
    den = 2 * E[:, 1] * E[:, 2]
    minang = np.degrees(np.arccos(np.clip(num / np.maximum(den, 1e-30), -1, 1)))
    fa = np.asarray(m.area_faces)
    med = float(np.median(fa))
    ee = np.sort(np.asarray(m.edges), axis=1)
    _, cnt = np.unique(ee, axis=0, return_counts=True)
    fs = np.sort(np.asarray(m.faces), axis=1)
    _, fc = np.unique(fs, axis=0, return_counts=True)
    used = np.unique(np.asarray(m.faces))
    return {
        "faces": int(len(m.faces)), "verts": int(len(m.vertices)),
        "folded": int((ang > 90).sum()), "max_dihedral": float(ang.max()) if len(ang) else 0.0,
        "sharp_60_90": int(((ang > 60) & (ang <= 90)).sum()),
        "slivers_5deg": int((minang < 5).sum()), "min_angle": float(minang.min()),
        "tiny_faces": int((fa < 0.01 * med).sum()), "min_area": float(fa.min()),
        "nonmanifold_edges": int((cnt > 2).sum()),
        "winding_consistent": bool(m.is_winding_consistent),
        "duplicate_faces": int((fc > 1).sum()),
        "unreferenced_verts": int(len(m.vertices) - len(used)),
        "euler": int(m.euler_number),
    }


OUT = ROOT / "tools" / f"remesh_defects_{SHARD}.jsonl"
with OUT.open("w", encoding="utf-8") as fh:
    for stl in sorted((ROOT / "output" / "nubs_final").glob("*_nub.stl"))[SHARD::NSH]:
        cid = stl.name.replace("_nub.stl", "")
        w = SM.weld(trimesh.load(str(stl), process=False))
        w, _ = repair_folds(w)
        for tf in DENS:
            try:
                rm = acvd_uniform_remesh(w, target_faces=tf)
                rec = scan(rm)
            except Exception as e:                                 # noqa: BLE001
                rec = {"error": f"{type(e).__name__}: {e}"}
            rec.update(case=cid, target=tf)
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            print(f"  {cid:10s} {tf:>7,}  folded {rec.get('folded','-'):>3} "
                  f"sharp {rec.get('sharp_60_90','-'):>4} slivers {rec.get('slivers_5deg','-'):>4} "
                  f"tiny {rec.get('tiny_faces','-'):>3}", flush=True)
