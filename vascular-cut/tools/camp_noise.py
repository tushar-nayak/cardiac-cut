"""Does repairing folded geometry change delta K by more than the noise?

    python tools/camp_noise.py run          # 30 cases x 3 densities x 2 variants
    python tools/camp_noise.py report

WHY THIS EXISTS

delta K is measured on a partition, and the partition is rebuilt from scratch
whenever the mesh changes at all. Repairing a fold moves ~0.05% of the surface --
enough to re-cluster everything -- so a before/after difference is NOT evidence
that the repair did anything. Twice already a change was attributed to a repair
that was really the partition landing differently.

So each case is measured at three face counts (98k / 100k / 102k). The spread
across those, holding the geometry fixed, IS that case's noise band: it is the same
surface, so anything it produces is partition sensitivity and nothing else. A
repair is then only believed where its effect is bigger than that band.

Runtime is ~180 CAMP runs, about 40 minutes at 6 workers.
"""
import argparse
import json
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

HERE = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "output" / "nubs_final"
DEFAULT_OUT = HERE / "camp_noise.jsonl"

FACES = (98_000, 100_000, 102_000)
TRIMS = (0.0, 0.001, 0.0025, 0.005, 0.01)


def _one(path_str, faces, repaired):
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = "1"
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(HERE))
    import time

    import numpy as np
    import trimesh

    import shape_metrics as SM
    from camp_cohort import repair_folds

    t0 = time.time()
    path = Path(path_str)
    mesh = trimesh.load(str(path), process=False)
    welded = SM.weld(mesh)
    chi = int(welded.euler_number)
    fix = {}
    if repaired:
        welded, fix = repair_folds(welded)
        if int(welded.euler_number) != chi:
            raise RuntimeError("repair changed topology")
    res = SM.delta_k(welded, target_faces=faces, trim_frac=0.0,
                     make_figure=False, expect_euler=chi)
    iK = np.asarray(res["integrated_K"], float)
    area = np.asarray(res["patch_area"], float)
    from camp import trim_mask_integrated
    dks = {}
    for t in TRIMS:
        keep = trim_mask_integrated(iK, t)
        q, a = iK[keep], area[keep]
        w = a / a.sum()
        dks[f"{t:g}"] = float(np.sum(w * (q - np.sum(w * q)) ** 2))
    return {
        "case": path.name.replace("_nub.stl", "").replace(".stl", ""),
        "file": path.name, "faces_target": faces, "repaired": bool(repaired),
        "dk": dks, "A_mm2": res["A_mm2"], "n_patches": res["n_patches"],
        "min_K": float(iK.min()), "max_K": float(iK.max()),
        "faces": res["faces"], "euler": res["euler_number"],
        "pct_at_min_size": res["pct_at_min_size"],
        **{f"repair_{k}": v for k, v in fix.items()},
        "seconds": round(time.time() - t0, 1),
    }


def cmd_run(args):
    meshes = sorted(Path(args.indir).glob("*.stl"))
    out = Path(args.out)
    done = set()
    if out.exists() and not args.fresh:
        for ln in out.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done.add((r["file"], r["faces_target"], r["repaired"]))
    jobs = [(m, f, rep) for m in meshes for f in FACES for rep in (False, True)
            if (m.name, f, rep) not in done]
    print(f"{len(meshes)} mesh(es) x {len(FACES)} densities x 2 variants = "
          f"{len(meshes) * len(FACES) * 2} runs, {len(jobs)} to do, "
          f"{args.workers} workers")
    if not jobs:
        return cmd_report(args)
    fails = []
    with out.open("a", encoding="utf-8") as fh, \
            ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_one, str(m), f, rep): (m, f, rep)
                for m, f, rep in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            m, f, rep = futs[fut]
            try:
                rec = fut.result()
            except Exception as e:                              # noqa: BLE001
                fails.append((m.name, f, rep, f"{type(e).__name__}: {e}"))
                print(f"  [{i}/{len(jobs)}] {m.stem} {f//1000}k "
                      f"{'repaired' if rep else 'raw':8s}: FAILED {e}", flush=True)
                if args.verbose:
                    traceback.print_exc()
                continue
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print(f"  [{i}/{len(jobs)}] {rec['case']:10s} {f//1000}k "
                  f"{'repaired' if rep else 'raw':8s} dK0={rec['dk']['0']:.5f} "
                  f"minK={rec['min_K']:7.2f} {rec['seconds']:.0f}s", flush=True)
    for name, f, rep, why in fails:
        print(f"  FAILED {name} {f} repaired={rep}: {why[:120]}")
    return cmd_report(args)


def cmd_report(args):
    import numpy as np

    rows = [json.loads(ln) for ln in
            Path(args.out).read_text(encoding="utf-8").splitlines() if ln.strip()]
    by = {}
    for r in rows:
        by[(r["case"], r["faces_target"], r["repaired"])] = r
    cases = sorted({r["case"] for r in rows})
    trim = args.trim

    print(f"\ndelta K at trim {trim:g}. 'noise' is the spread over 98k/100k/102k "
          f"on the SAME geometry;\n'effect' is repaired/raw at 100k. A repair is "
          f"believed only where effect is outside the noise.\n")
    print(f"{'case':10s} {'folds':>5} {'raw noise':>10} {'rep noise':>10} "
          f"{'raw dK':>9} {'rep dK':>9} {'effect':>8}  verdict")
    verdicts = {}
    for c in cases:
        raw = [by[(c, f, False)]["dk"][f"{trim:g}"] for f in FACES
               if (c, f, False) in by]
        rep = [by[(c, f, True)]["dk"][f"{trim:g}"] for f in FACES
               if (c, f, True) in by]
        if len(raw) < len(FACES) or len(rep) < len(FACES):
            continue
        folds = by[(c, 100_000, True)].get("repair_folded_edges", 0)
        n_raw = max(raw) / min(raw)
        n_rep = max(rep) / min(rep)
        r0 = by[(c, 100_000, False)]["dk"][f"{trim:g}"]
        p0 = by[(c, 100_000, True)]["dk"][f"{trim:g}"]
        eff = p0 / r0
        band = max(n_raw, n_rep)
        if not folds:
            v = "no folds"
        elif eff < 1.0 / band:
            v = "REAL improvement"
        elif eff > band:
            v = "REAL worsening"
        else:
            v = "within noise"
        verdicts[c] = v
        print(f"{c:10s} {folds:5d} {n_raw:9.2f}x {n_rep:9.2f}x {r0:9.5f} "
              f"{p0:9.5f} {eff:7.2f}x  {v}")

    print()
    for v in ("REAL improvement", "REAL worsening", "within noise", "no folds"):
        n = [c for c, x in verdicts.items() if x == v]
        print(f"  {v:18s} {len(n):2d}   {', '.join(n) if len(n) <= 8 else ''}")

    # How big is the noise itself, and does it track folding?
    print(f"\n  noise band (raw geometry) against fold count:")
    for c in sorted(cases, key=lambda c: -by[(c, 100_000, True)].get(
            "repair_folded_edges", 0))[:10]:
        raw = [by[(c, f, False)]["dk"][f"{trim:g}"] for f in FACES
               if (c, f, False) in by]
        if len(raw) == len(FACES):
            print(f"    {c:10s} {by[(c, 100_000, True)].get('repair_folded_edges', 0):3d} "
                  f"fold(s)   noise {max(raw)/min(raw):5.2f}x")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "report"):
        p = sub.add_parser(name)
        p.add_argument("--indir", default=str(DEFAULT_IN))
        p.add_argument("--out", default=str(DEFAULT_OUT))
        p.add_argument("--workers", type=int, default=6)
        p.add_argument("--trim", type=float, default=0.0)
        p.add_argument("--fresh", action="store_true")
        p.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    return {"run": cmd_run, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
