"""Is delta K's instability partition-sampling noise, and does averaging remove it?

Changing only the k-means seed moves delta K by 40-50% on a byte-identical mesh, which
is more than a 4x change in mesh density does. So the partition is the noisy element:
delta K is carried by ~22 patches out of ~850, and which 22 depends on where k-means
lands.

If those partitions are independent samples of the same underlying quantity, the mean
of n of them has a standard error falling as 1/sqrt(n) -- 4 replicates halve the
spread, 9 third it. If instead the spread plateaus, the partitions are correlated and
averaging cannot rescue it.

Both are measured: the spread ACROSS SEEDS at one density (partition noise on a fixed
mesh) and the spread ACROSS DENSITIES (what the error bars on the cohort plot show).
If averaging seeds shrinks the second, the error bars narrow.

    python tools/seed_average.py <shard> <nshards>
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
CASES = ["KK28_11", "KK20_2", "KK32_0", "KY160_1", "KY7_1", "SA21_1"]
DENS = (100_000, 200_000, 400_000)
SEEDS = [42, 43, 44, 45, 46, 47, 48, 49]

OUT = ROOT / "tools" / f"seed_average_{SHARD}.jsonl"
t0 = time.time()
with OUT.open("w", encoding="utf-8") as fh:
    for cid in CASES[SHARD::NSH]:
        base = SM.weld(trimesh.load(str(ROOT / "output" / "nubs_final" / f"{cid}_nub.stl"),
                                    process=False))
        chi = int(base.euler_number)
        base, _ = repair_folds(base)
        for tf in DENS:
            for sd in SEEDS:
                orig = SM.PARAMS.copy()
                SM.PARAMS["seed"] = sd
                try:
                    r = SM.delta_k(base, target_faces=tf, trim_frac=0.0,
                                   make_figure=False, expect_euler=chi)
                    rec = {"case": cid, "faces": tf, "seed": sd,
                           "dk": float(r["k_fluct"]), "n_patches": int(r["n_patches"])}
                except Exception as e:                             # noqa: BLE001
                    rec = {"case": cid, "faces": tf, "seed": sd,
                           "error": f"{type(e).__name__}: {e}"}
                finally:
                    SM.PARAMS.clear(); SM.PARAMS.update(orig)
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                print(f"  {cid:9s} {tf:>7,} seed {sd}  "
                      f"dk {rec.get('dk', float('nan')):.5f}  "
                      f"[{(time.time()-t0)/60:.0f} min]", flush=True)
