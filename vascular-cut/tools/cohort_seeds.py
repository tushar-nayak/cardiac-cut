"""The whole nub cohort at 100k, partitioned 8 independent ways.

One partition is one draw: delta K is carried by ~22 patches out of ~850, sited on the
branch ostia, and where k-means puts the boundary there decides the number. Averaging
independent partitions is the matched fix -- CV falls as 1/sqrt(k), so 8 seeds takes a
case from about +/- 16% to +/- 6%.

Density is fixed at 100,000 deliberately. Averaging removes seed noise; it does NOT
remove the systematic drift between densities, so mixing densities here would put back
the one thing this cannot fix.

    python tools/cohort_seeds.py <shard> <nshards>
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tools"))
import trimesh                                                     # noqa: E402
import shape_metrics as SM                                         # noqa: E402
from camp_cohort import repair_folds                               # noqa: E402

SHARD = int(sys.argv[1]) if len(sys.argv) > 1 else 0
NSH = int(sys.argv[2]) if len(sys.argv) > 2 else 1
SRC = sys.argv[3] if len(sys.argv) > 3 else "output/nubs_final"
SUF = sys.argv[4] if len(sys.argv) > 4 else "_nub"
TAG = sys.argv[5] if len(sys.argv) > 5 else ""
SEEDS = list(range(42, 50))
FACES = 100_000

OUT = ROOT / "tools" / f"cohort_seeds{TAG}_{SHARD}.jsonl"
t0 = time.time()
with OUT.open("w", encoding="utf-8") as fh:
    for stl in sorted((ROOT / SRC).glob(f"*{SUF}.stl"))[SHARD::NSH]:
        cid = stl.name.replace(f"{SUF}.stl", "")
        base = SM.weld(trimesh.load(str(stl), process=False))
        chi = int(base.euler_number)
        base, _ = repair_folds(base)
        for sd in SEEDS:
            orig = SM.PARAMS.copy(); SM.PARAMS["seed"] = sd
            try:
                r = SM.delta_k(base, target_faces=FACES, trim_frac=0.0,
                               make_figure=False, expect_euler=chi)
                rec = {"case": cid, "seed": sd, "dk": float(r["k_fluct"]),
                       "area": float(r["A_mm2"]), "n_patches": int(r["n_patches"])}
            except Exception as e:                                 # noqa: BLE001
                rec = {"case": cid, "seed": sd, "error": f"{type(e).__name__}: {e}"}
            finally:
                SM.PARAMS.clear(); SM.PARAMS.update(orig)
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            print(f"  {cid:10s} seed {sd}  dk {rec.get('dk', float('nan')):.5f}  "
                  f"[{(time.time()-t0)/60:.0f} min]", flush=True)
