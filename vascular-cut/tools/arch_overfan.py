"""How much of an arch's delta K comes from patches that span an ostium?

An ostium-spanning patch is one whose FACE normals fan past 90 degrees: no height
field describes it, so its Monge fit -- and therefore its K_center, and therefore its
contribution -- is not a curvature measurement of anything. The fan gate cannot catch
these because it judges vertex sets at fit time and these are created afterwards, by
the vertex-to-face label conversion.

Two numbers per case: the share of delta K those patches carry, and what delta K
becomes with them dropped. If the share is small the arch instability lies elsewhere;
if it is large, the arch cohort is largely measuring its own ostia.
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
OUT = ROOT / "tools" / f"arch_overfan_{SHARD}.jsonl"
FACES = 50_000

t0 = time.time()
with OUT.open("w", encoding="utf-8") as fh:
    for stl in sorted((ROOT / "output" / "arch_final").glob("*_arch.stl"))[SHARD::NSH]:
        cid = stl.name.replace("_arch.stl", "")
        w = SM.weld(trimesh.load(str(stl), process=False))
        chi = int(w.euler_number)
        w, _ = repair_folds(w)
        r = SM.delta_k(w, target_faces=FACES, trim_frac=0.0, make_figure=False,
                       expect_euler=chi)
        K = np.asarray(r["integrated_K"], float)
        A = np.asarray(r["patch_area"], float)
        m = r["mesh"]; lab = np.asarray(r["face_labels"], int)
        fn = np.asarray(m.face_normals, float); fa = np.asarray(m.area_faces, float)
        fan = np.zeros(len(K))
        for p in range(len(K)):
            sel = lab == p
            if not sel.any():
                continue
            n = fn[sel]; mm = (n * fa[sel, None]).sum(0)
            mm /= max(np.linalg.norm(mm), 1e-12)
            fan[p] = float(np.degrees(np.arccos(np.clip(n @ mm, -1, 1))).max())
        over = fan > 90.0
        wt = A / A.sum()
        contrib = wt * (K - float((wt * K).sum())) ** 2
        dk = float(contrib.sum())
        share = 100 * float(contrib[over].sum()) / max(dk, 1e-30)
        keep = ~over
        if keep.sum() >= 2:
            w2 = A[keep] / A[keep].sum()
            dk_wo = float(np.sum(w2 * (K[keep] - np.sum(w2 * K[keep])) ** 2))
        else:
            dk_wo = float("nan")
        rec = {"case": cid, "area": float(r["A_mm2"]), "n_patches": int(len(K)),
               "n_over": int(over.sum()), "pct_patches_over": 100 * float(over.mean()),
               "area_frac_over": 100 * float(A[over].sum() / A.sum()),
               "share_dk": share, "dk": dk, "dk_without_over": dk_wo}
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        print(f"  {cid:10s} {int(over.sum()):2d}/{len(K)} patches "
              f"({100*over.mean():4.1f}%)  {A[over].sum()/A.sum()*100:4.1f}% of area  "
              f"-> {share:5.1f}% of dK   [{(time.time()-t0)/60:.0f} min]", flush=True)
