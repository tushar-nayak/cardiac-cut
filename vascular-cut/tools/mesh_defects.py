"""Inventory the mesh defects in a folder of cut surfaces.

    python tools/mesh_defects.py                    # output/nubs_final
    python tools/mesh_defects.py --indir <dir>

Reports, per mesh, the things that make a local curvature fit meaningless:

  folds        face pairs meeting at more than FOLD_DEG -- the surface doubling
               back on itself. These are the ones that matter: on KP35_1, 46 of
               them within one small region produced a patch at integrated K =
               -41.7 carrying 54% of the whole case's delta K, against a next-worst
               of -5.4.
  slivers      triangles so thin that a normal is unreliable
  nonmanifold  edges shared by more than two faces
  duplicated   coincident vertices left after welding

Folds are counted on the WELDED surface, because an unwelded STL has no face
adjacency at all and every check silently returns zero.

Cheap -- no CAMP, no remesh, about a second per mesh -- so it is worth running
over any cohort before trusting a shape metric computed on it.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import shape_metrics as SM                                  # noqa: E402

# The check itself lives in shape_metrics, beside weld(), so that the web app and the
# desktop tool can warn about a folded input without importing anything from tools/.
FOLD_DEG = SM.FOLD_DEG
SLIVER_Q = SM.SLIVER_Q
defects = SM.surface_defects


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=str(ROOT / "output" / "nubs_final"))
    ap.add_argument("--fold-deg", type=float, default=FOLD_DEG)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    meshes = sorted(Path(args.indir).glob("*.stl"))
    if not meshes:
        raise SystemExit(f"no .stl under {args.indir}")
    rows = []
    print(f"{'case':11s} {'faces':>7} {'folds':>6} {'max dih':>8} {'span':>7} "
          f"{'slivers':>8} {'minQ':>6} {'nonman':>7} {'euler':>6}")
    for p in meshes:
        d = defects(trimesh.load(str(p), process=False), args.fold_deg)
        d["case"] = p.name.replace("_nub.stl", "").replace(".stl", "")
        rows.append(d)
        span = f"{d['fold_span_mm']:6.1f}" if d["fold_span_mm"] is not None else "     -"
        print(f"{d['case']:11s} {d['faces']:7d} {d['folds']:6d} "
              f"{d['max_dihedral']:7.1f}° {span} {d['slivers']:8d} "
              f"{d['min_quality']:6.3f} {d['nonmanifold_edges']:7d} {d['euler']:6d}")

    bad = [r for r in rows if r["folds"]]
    print(f"\n  {len(bad)} of {len(rows)} meshes carry folded geometry")
    if bad:
        worst = sorted(bad, key=lambda r: -r["folds"])[:5]
        print("  worst: " + ", ".join(f"{r['case']} ({r['folds']})" for r in worst))
        print(f"  folds are LOCAL: median cluster radius "
              f"{np.median([r['fold_span_mm'] for r in bad]):.1f} mm")
    nm = [r for r in rows if r["nonmanifold_edges"]]
    if nm:
        print(f"  {len(nm)} mesh(es) have non-manifold edges: "
              + ", ".join(r["case"] for r in nm))
    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=1), encoding="utf-8")
        print(f"  wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
