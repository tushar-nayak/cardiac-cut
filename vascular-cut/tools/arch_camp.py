"""CAMP over the locked arch set, at the density the arch sweep picked.

50,000 faces, and NOT the 100,000 the nubs use. The arches were swept separately
because they are a different object: 8,548 mm2 against 45,514, with the same five
openings, so ostia are 5.3x more concentrated. Three criteria were set in advance and
the sweep answered them like this:

    faces    % of vertex ceiling   source triangles per patch   rho vs 2x density
    10,000        48.5%                    40.0                      +0.56
    20,000        36.1%                    25.8                       n/a
    50,000        26.3%                    14.2                      +0.60
   100,000        20.2%                     9.1                       n/a

50k is the only density clearing the source-resolution criterion with margin. The arch
source meshes hold a median 5,019 faces, so at 100k a patch spans 9.1 original
triangles and is largely measuring subdivision; KY186_1 also fails outright at 100k,
the remesher closing an opening (euler -3 -> -2).

READ THE STABILITY CAVEAT BEFORE USING THESE NUMBERS. No arch density reached the
rho >= 0.8 agreement across a doubling that the nubs met at 100k (+0.91). The best
here is +0.60. These values are usable for a within-run comparison; they are not
established as reproducible across a mesh change the way the nub numbers are.

The fan and over-fan counts are kept per case because the arch geometries carry three
branch ostia in a much shorter vessel than the nubs did, so the ostium-spanning patches
that the fan gate does not catch should be commoner here, and that is worth seeing
rather than assuming.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import trimesh                                                     # noqa: E402

import shape_metrics as SM                                         # noqa: E402
from camp_cohort import repair_folds                               # noqa: E402

cu = SM.load_campsrc()
FOLDER = ROOT / "output" / "arch_final"
OUT = ROOT / "tools" / "camp_arch_50k.jsonl"
FACES = 50_000

done = set()
if OUT.exists():
    done = {json.loads(l)["case"] for l in OUT.read_text(encoding="utf-8").splitlines()
            if l.strip()}

t0 = time.time()
with OUT.open("a", encoding="utf-8") as fh:
    for stl in sorted(FOLDER.glob("*_arch.stl")):
        cid = stl.name.replace("_arch.stl", "")
        if cid in done:
            continue
        welded = SM.weld(trimesh.load(str(stl), process=False))
        chi = int(welded.euler_number)
        welded, _ = repair_folds(welded)
        try:
            r = SM.delta_k(welded, target_faces=FACES, trim_frac=0.0,
                           make_figure=False, expect_euler=chi)
        except Exception as e:                                     # noqa: BLE001
            fh.write(json.dumps({"case": cid, "error": f"{type(e).__name__}: {e}"}) + "\n")
            fh.flush()
            print(f"  {cid:10s} {type(e).__name__}: {e}", flush=True)
            continue
        K = np.asarray(r["integrated_K"], float)
        A = np.asarray(r["patch_area"], float)
        m = r["mesh"]
        lab = np.asarray(r["face_labels"], int)
        fn = np.asarray(m.face_normals, float)
        fa = np.asarray(m.area_faces, float)
        nover = 0
        for p in range(len(K)):
            sel = lab == p
            if not sel.any():
                continue
            n = fn[sel]
            mm = (n * fa[sel, None]).sum(0)
            mm /= max(np.linalg.norm(mm), 1e-12)
            if np.degrees(np.arccos(np.clip(n @ mm, -1, 1))).max() > 90:
                nover += 1
        rec = {"case": cid, "area": float(r["A_mm2"]), "n_patches": int(len(K)),
               "over_fan": nover, "faces": int(r["faces"]),
               "pct_at_min_size": float(r["pct_at_min_size"]), "dk": {}}
        for t in (0.0, 0.0025):
            keep = cu._trim_mask_integrated(K, t)
            q, a = K[keep], A[keep]
            w = a / a.sum()
            rec["dk"][str(t)] = float(np.sum(w * (q - np.sum(w * q)) ** 2))
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        print(f"  {cid:10s} area {rec['area']:7.0f}  N {len(K):4d}  over-fan {nover:2d}  "
              f"dk0 {rec['dk']['0.0']:.5f}  [{(time.time() - t0) / 60:.0f} min]", flush=True)
print(f"\nwrote {OUT}")
