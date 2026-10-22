"""Precompute the integrated-K colouring for every geometry, once.

    .venv\Scripts\python tools\camp_patch_cache.py
    .venv\Scripts\python tools\camp_patch_cache.py --dir output/arch_final --suffix _arch

CAMP takes the better part of a minute per case at 100k faces, which is fine for a
cohort run and unusable behind an arrow key. So the partition is computed here and
cached, and the viewer becomes a file read.

What is stored is exactly what campsrc's own figure draws: the decimated mesh, one
integrated-K value PER FACE (its patch's K, as face_values = integrated_K[f_labels]),
and the edges where two patches meet. Nothing is recomputed at draw time, so the
picture in the viewer is the picture CAMP would have made.

Trim is 0 deliberately. A trimmed patch is hidden from campsrc's figure entirely, and
a viewer whose purpose is to check that the number matches the geometry must not start
by hiding the patches carrying the number. The viewer can trim on the way out.
"""
import argparse
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


WANT = ("verts", "faces", "face_K", "face_patch", "patch_K", "patch_area",
        "patch_share", "seam_edges")


def _is_ours(path):
    """True only for a cache file this script wrote, judged by its keys."""
    try:
        with np.load(path) as z:
            return all(k in z.files for k in WANT)
    except Exception:                                              # noqa: BLE001
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="output/nubs_final")
    ap.add_argument("--suffix", default="_nub")
    ap.add_argument("--out", default="output/camp_intk")
    ap.add_argument("--faces", type=int, default=100_000)
    args = ap.parse_args()

    folder = ROOT / args.dir
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    stls = sorted(folder.glob(f"*{args.suffix}.stl"))
    if not stls:
        raise SystemExit(f"no *{args.suffix}.stl in {folder}")

    t0 = time.time()
    for stl in stls:
        cid = stl.name.replace(f"{args.suffix}.stl", "")
        dst = out / f"{cid}.npz"
        # Presence is not enough. output/camp_patches already held npz files of
        # per-patch STATISTICS from other work, and a bare exists() check silently
        # skipped every case and left the viewer reading the wrong schema.
        if dst.exists() and _is_ours(dst):
            print(f"  {cid:10s} cached", flush=True)
            continue
        welded = SM.weld(trimesh.load(str(stl), process=False))
        chi = int(welded.euler_number)
        welded, _ = repair_folds(welded)
        r = SM.delta_k(welded, target_faces=args.faces, trim_frac=0.0,
                       make_figure=False, expect_euler=chi)

        m = r["mesh"]
        lab = np.asarray(r["face_labels"], dtype=np.int32)
        K = np.asarray(r["integrated_K"], dtype=float)
        A = np.asarray(r["patch_area"], dtype=float)
        adj = np.asarray(m.face_adjacency)
        seam = lab[adj[:, 0]] != lab[adj[:, 1]]

        # Each patch's share of delta K, kept per patch so the viewer can rank and
        # highlight without the arrays it would need to recompute them.
        w = A / A.sum()
        contrib = w * (K - float((w * K).sum())) ** 2

        np.savez_compressed(
            dst,
            verts=np.asarray(m.vertices, dtype=np.float32),
            faces=np.asarray(m.faces, dtype=np.int32),
            face_K=K[lab].astype(np.float32),
            face_patch=lab,
            patch_K=K.astype(np.float32),
            patch_area=A.astype(np.float32),
            patch_share=(contrib / max(contrib.sum(), 1e-30)).astype(np.float32),
            seam_edges=np.asarray(m.face_adjacency_edges)[seam].astype(np.int32),
            area=np.float32(r["A_mm2"]),
            dk=np.float32(contrib.sum()),
        )
        print(f"  {cid:10s} {len(K):4d} patches  K [{K.min():+.3f}, {K.max():+.3f}]  "
              f"-> {dst.name} ({dst.stat().st_size/1e6:.1f} MB)  "
              f"[{(time.time()-t0)/60:.0f} min]", flush=True)

    print(f"\ncached {len(list(out.glob('*.npz')))} cases in {out}")


if __name__ == "__main__":
    raise SystemExit(main())
