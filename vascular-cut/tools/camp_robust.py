"""Is the partition robust to a perturbed surface, and to the mesh it is given?

Two questions, and they are not the same one:

    jitter    displace every vertex of the SOURCE surface by Gaussian noise, then run
              the whole pipeline. Sigma is a fraction of that mesh's own mean edge
              length, so the perturbation is comparable across aortas of different
              size. Applied BEFORE decimation deliberately: a shape metric has to
              survive the segmentation being slightly different, and that includes the
              remesh reacting to it, not just CAMP reacting to a fixed remesh.

    density   the same surface at 25k, 50k, 100k and 200k faces. delta K is known to
              depend on face count; what matters for using it is whether the COHORT
              ORDER survives, since every clinical reading is a comparison.

Replicate 0 of the jitter run is unjittered. KMeans is seeded throughout
(random_state=42 globally, depth+1 per split), so it must come back bit-identical to
the reference run -- if it does not, the pipeline is not deterministic and every other
number here is measuring the wrong thing.

    python tools/camp_robust.py --mode jitter  --shard 0 --nshards 6
    python tools/camp_robust.py --mode density --shard 0 --nshards 6
"""
import argparse
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

DENSITIES = (25_000, 50_000, 100_000, 200_000)
DENSITIES_HI = (300_000, 400_000)
JITTER_REPLICATES = 4          # replicate 0 is the unjittered control
JITTER_FRAC = 0.25             # sigma, as a fraction of mean edge length

# Absolute displacements, in mm, for the sweep. Stated in mm rather than as a fraction
# of edge length because the question is clinical, not numerical: these surfaces come
# from CT at 0.5-1.0 mm voxels, so the perturbation that matters is the one a
# segmentation could plausibly be wrong by. 0.05 mm is far below anything real; 0.40 mm
# is most of a voxel.
JITTER_SIGMAS_MM = (0.05, 0.10, 0.20, 0.40)
SWEEP_REPLICATES = 2
TRIMS = (0.0, 0.0025)


def mean_edge(mesh):
    e = np.asarray(mesh.edges_unique)
    v = np.asarray(mesh.vertices)
    return float(np.linalg.norm(v[e[:, 0]] - v[e[:, 1]], axis=1).mean())


def measure(mesh, chi, faces, cu):
    r = SM.delta_k(mesh, target_faces=faces, trim_frac=0.0,
                   make_figure=False, expect_euler=chi)
    K = np.asarray(r["integrated_K"], float)
    A = np.asarray(r["patch_area"], float)
    out = {"n_patches": int(len(K)), "faces": int(r["faces"]),
           "area": float(r["A_mm2"]), "dk": {}}
    for t in TRIMS:
        keep = cu._trim_mask_integrated(K, t)
        q, a = K[keep], A[keep]
        w = a / a.sum()
        out["dk"][str(t)] = float(np.sum(w * (q - np.sum(w * q)) ** 2))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode",
                    choices=("jitter", "density", "sweep", "orphan", "hidensity"),
                    required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--dir", default="output/nubs_final")
    ap.add_argument("--suffix", default="_nub")
    ap.add_argument("--tag", default="",
                    help="suffix for the output file, to keep cohorts apart")
    ap.add_argument("--densities", default="",
                    help="comma-separated face counts, overriding the defaults")
    ap.add_argument("--mps", type=int, default=0,
                    help="override min_patch_size (vertices); 0 keeps the default 18")
    args = ap.parse_args()

    cu = SM.load_campsrc()
    if args.mps:
        SM.PARAMS["min_patch_size"] = args.mps
        print(f"  min_patch_size = {args.mps} vertices "
              f"(the Monge quadratic has 6 coefficients, so this is thin)")
    folder = ROOT / args.dir
    stls = sorted(folder.glob(f"*{args.suffix}.stl"))[args.shard::args.nshards]
    densities = ([int(x) for x in args.densities.split(",")]
                 if args.densities else DENSITIES)
    tag = f"_{args.tag}" if args.tag else ""
    out = ROOT / "tools" / f"camp_robust_{args.mode}{tag}_{args.shard}.jsonl"
    done = set()
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                done.add((r["case"], r.get("rep", r.get("faces"))))

    t0 = time.time()
    with out.open("a", encoding="utf-8") as fh:
        for stl in stls:
            cid = stl.name.replace(f"{args.suffix}.stl", "")
            base = SM.weld(trimesh.load(str(stl), process=False))
            chi = int(base.euler_number)
            sigma = JITTER_FRAC * mean_edge(base)

            if args.mode == "jitter":
                jobs = [("rep", i) for i in range(JITTER_REPLICATES)]
            elif args.mode == "sweep":
                jobs = [("job", f"{sg}:{k}") for sg in JITTER_SIGMAS_MM
                        for k in range(SWEEP_REPLICATES)]
            elif args.mode == "hidensity":
                jobs = [("faces", f) for f in DENSITIES_HI]
            elif args.mode == "orphan":
                # Only the two densities where the floor has stopped binding. Below
                # 100k the partition is decided by the vertex budget, so an orphan rule
                # cannot be read there without confounding the two.
                jobs = [("faces", f) for f in (100_000, 200_000)]
            else:
                jobs = [("faces", f) for f in densities]
            for key, val in jobs:
                if (cid, val) in done:
                    continue
                m = base.copy()
                use_sigma = 0.0
                if args.mode == "jitter" and val > 0:
                    use_sigma = sigma
                elif args.mode == "sweep":
                    use_sigma = float(str(val).split(":")[0])
                if use_sigma > 0:
                    # Seeded per (case, job) so a re-run reproduces exactly.
                    rng = np.random.default_rng(abs(hash((cid, val))) % (2 ** 32))
                    m.vertices = np.asarray(m.vertices) + rng.normal(
                        0.0, use_sigma, size=np.asarray(m.vertices).shape)
                m, _ = repair_folds(m)
                faces = 100_000 if args.mode in ("jitter", "sweep") else val
                if args.mode == "orphan":
                    SM.PARAMS["orphan_adjacency"] = True
                try:
                    rec = measure(m, chi, faces, cu)
                except Exception as e:                             # noqa: BLE001
                    rec = {"error": f"{type(e).__name__}: {e}"}
                rec.update({"case": cid, key: val, "mode": args.mode})
                if args.mode == "jitter":
                    rec["sigma_mm"] = sigma if val > 0 else 0.0
                elif args.mode == "sweep":
                    rec["sigma_mm"] = use_sigma
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                tag = rec.get("error", f"dk0 {rec.get('dk', {}).get('0.0', float('nan')):.5f}")
                print(f"  [{args.shard}] {cid:10s} {key}={val:<7} {tag}  "
                      f"[{(time.time() - t0) / 60:.0f} min]", flush=True)
    print(f"shard {args.shard} done -> {out}")


if __name__ == "__main__":
    raise SystemExit(main())
