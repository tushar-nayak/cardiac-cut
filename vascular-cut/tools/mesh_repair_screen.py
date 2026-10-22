"""Which post-remesh repair threshold removes mesh noise without altering geometry?

repair_folds only touches dihedrals past FOLD_DEG = 90. The scan found the population
that correlates with instability is mostly 60-90 degrees -- creases too sharp for one
quadratic to span, invisible to the current threshold. So sweep the threshold.

Two things are measured against each other, because they trade off:

    CV of delta K across 100k/200k/400k   how much mesh noise is left
    moved vertices, mean move, area change  how much geometry was destroyed to get it

A repair that smooths a fifth of the surface is not a repair. The guardrail matters
more than the CV: any threshold that moves area by more than a fraction of a percent
is changing the object being measured.

Winding is fixed first where it is broken -- inconsistent face orientation flips
normal signs, which corrupts both the fan measure and the Monge frame.

    python tools/mesh_repair_screen.py <shard> <nshards>
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tools"))
import trimesh                                                     # noqa: E402
import shape_metrics as SM                                         # noqa: E402
from camp_cohort import repair_folds                               # noqa: E402

SHARD = int(sys.argv[1]) if len(sys.argv) > 1 else 0
NSH = int(sys.argv[2]) if len(sys.argv) > 2 else 1
CASES = ["KK28_11", "KK20_2", "KK32_0", "KR52_1", "KY160_1", "KP35_1", "KY7_1", "KR73_1"]
DENS = (100_000, 200_000, 400_000)
THRESH = (None, 90.0, 75.0, 60.0, 45.0)


def make_repair(deg):
    if deg is None:
        return None, {}
    stats = {}

    def f(m):
        mm = m.copy()
        if not mm.is_winding_consistent:
            mm.fix_normals()
            stats["winding_fixed"] = True
        out, info = repair_folds(mm, fold_deg=deg)
        stats.update(info)
        return out, info
    return f, stats


OUT = ROOT / "tools" / f"mesh_repair_screen_{SHARD}.jsonl"
with OUT.open("w", encoding="utf-8") as fh:
    for cid in CASES[SHARD::NSH]:
        base = SM.weld(trimesh.load(str(ROOT / "output" / "nubs_final" / f"{cid}_nub.stl"),
                                    process=False))
        chi = int(base.euler_number)
        base, _ = repair_folds(base)
        for deg in THRESH:
            for tf in DENS:
                fn, stats = make_repair(deg)
                try:
                    r = SM.delta_k(base, target_faces=tf, trim_frac=0.0,
                                   make_figure=False, expect_euler=chi,
                                   post_remesh_repair=fn)
                    rec = {"case": cid, "deg": deg, "faces": tf,
                           "dk": float(r["k_fluct"]), "n_patches": int(r["n_patches"]),
                           "moved": int(stats.get("moved_vertices", 0)),
                           "mean_move_mm": float(stats.get("mean_move_mm", 0.0)),
                           "area_change_pct": float(stats.get("area_change_pct", 0.0)),
                           "folded_before": int(stats.get("folded_edges", 0)),
                           "folded_after": int(stats.get("folded_after", 0)),
                           "winding_fixed": bool(stats.get("winding_fixed", False))}
                except Exception as e:                             # noqa: BLE001
                    rec = {"case": cid, "deg": deg, "faces": tf,
                           "error": f"{type(e).__name__}: {e}"}
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                print(f"  {cid:9s} deg={str(deg):>5s} {tf:>7,}  "
                      f"dk {rec.get('dk', float('nan')):.5f}  "
                      f"moved {rec.get('moved', 0):5d}  "
                      f"area {rec.get('area_change_pct', 0):+.3f}%", flush=True)
