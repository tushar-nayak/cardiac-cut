"""What the two CAMP partition fixes did to the UNBRANCHED cohort.

The connectivity split and the fan gate were both measured on the BRANCHED nub
cohort, where they mattered a great deal: the normals band went from 0.14-3.47
(CV 1.04) to 0.66-1.37 (CV 0.23) at trim 0. Every unbranched run in this repository
was made AFTER both fixes, so nothing on disk answers the obvious question --

    did the fixes also move the unbranched cohort?

It matters for the write-up. Both defects are things that happen where two surfaces
come close through the air but far apart across the surface, and where a patch can
wrap around something: on a branched aorta that means the branch ostia. An unbranched
tube has none. So the EXPECTATION is little or no change -- but that is a mechanism
argument, and a mechanism argument is not a measurement.

This runs all 30 unbranched meshes twice at a fixed 100k faces and a fixed seed,
changing only the two flags:

    OFF   split_disconnected=False, fan_thresh=None     the pre-fix partition
    ON    the current defaults                          what we ship

Same mesh, same remesh target, same seed, so the ONLY difference is the partition
rules. Reports per case and, more importantly, whether the cohort's shape survives:
the normals band, the normal-vs-diseased separation, and the rank correlation.

    .venv\\Scripts\\python tools\\unbranched_fix_effect.py [--faces 100000]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import shape_metrics as SM                                          # noqa: E402

OUT = ROOT / "tools" / "unbranched_fix_effect.jsonl"
SRC = ROOT / "output" / "unbranched_src"
# KY and SA are the healthy controls, as everywhere else in this cohort.
NORMAL = ("KY", "SA")


def is_normal(case):
    return case.startswith(NORMAL)


def run_one(stl, faces, flags):
    """delta K for one mesh under one set of partition flags."""
    saved = dict(SM.PARAMS)
    try:
        SM.PARAMS.update(flags)
        r = SM.delta_k(str(stl), target_faces=faces, trim_frac=0.0,
                       make_figure=False)
    finally:
        SM.PARAMS.clear()
        SM.PARAMS.update(saved)
    return {
        "dk": float(r["k_fluct"]),
        "area": float(r["A_mm2"]),
        "n_patches": int(r["n_patches"]),
        "faces": int(r["faces"]),
        "disconnected": int(r.get("n_disconnected", -1)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faces", type=int, default=100_000)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    out = Path(args.out)

    # Both fixes off is the PRE-FIX partition; the defaults are what we ship.
    CONDITIONS = {
        "off": dict(split_disconnected=False, fan_thresh=None),
        "on": dict(split_disconnected=True, fan_thresh=90.0),
    }

    stls = sorted(SRC.glob("*_camp_input.stl"))
    if not stls:
        raise SystemExit(f"no unbranched meshes under {SRC}")
    done = set()
    if out.exists():
        for line in out.open(encoding="utf-8"):
            if line.strip():
                r = json.loads(line)
                done.add((r["case"], r["cond"]))
        print(f"resuming: {len(done)} already measured")

    print(f"{len(stls)} meshes x 2 conditions at {args.faces:,} faces, seed 42")
    t0 = time.time()
    with out.open("a", encoding="utf-8") as fh:
        for i, stl in enumerate(stls, 1):
            case = stl.name.replace("_camp_input.stl", "")
            for cond, flags in CONDITIONS.items():
                if (case, cond) in done:
                    continue
                try:
                    rec = run_one(stl, args.faces, flags)
                    rec.update(case=case, cond=cond)
                except Exception as e:                              # noqa: BLE001
                    rec = {"case": case, "cond": cond,
                           "error": f"{type(e).__name__}: {e}"}
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                el = time.time() - t0
                print(f"  [{i:2d}/{len(stls)}] {case:<10} {cond:<3} "
                      f"{rec.get('dk', float('nan')):.6f}  ({el / 60:.1f} min)")
    report(out)


def report(out):
    by = {}
    for line in Path(out).open(encoding="utf-8"):
        if not line.strip():
            continue
        r = json.loads(line)
        if "error" in r:
            continue
        by.setdefault(r["case"], {})[r["cond"]] = r
    both = {c: v for c, v in by.items() if "off" in v and "on" in v}
    if not both:
        print("nothing to report yet")
        return

    cases = sorted(both, key=lambda c: both[c]["on"]["dk"])
    off = np.array([both[c]["off"]["dk"] for c in cases])
    on = np.array([both[c]["on"]["dk"] for c in cases])
    ratio = on / np.where(off == 0, np.nan, off)

    print()
    print("=" * 78)
    print(f"UNBRANCHED COHORT, {len(cases)} cases, trim 0, both fixes OFF vs ON")
    print("=" * 78)
    print(f"{'case':<12}{'group':<10}{'OFF':>11}{'ON':>11}{'ON/OFF':>9}"
          f"{'disc':>7}")
    for c in cases:
        d = both[c]
        print(f"{c:<12}{'normal' if is_normal(c) else 'diseased':<10}"
              f"{d['off']['dk']:>11.6f}{d['on']['dk']:>11.6f}"
              f"{d['on']['dk'] / d['off']['dk']:>9.3f}"
              f"{d['off'].get('disconnected', -1):>7}")

    print()
    print(f"  median ON/OFF ratio : {np.nanmedian(ratio):.4f}")
    print(f"  worst move          : {np.nanmax(np.abs(ratio - 1)) * 100:.1f}% "
          f"({cases[int(np.nanargmax(np.abs(ratio - 1)))]})")
    print(f"  cases moving >10%   : {int(np.nansum(np.abs(ratio - 1) > 0.10))} "
          f"of {len(cases)}")
    print(f"  cases moving >2%    : {int(np.nansum(np.abs(ratio - 1) > 0.02))} "
          f"of {len(cases)}")

    # The thing that actually matters: does the cohort keep its shape and order?
    from scipy.stats import spearmanr                               # noqa: PLC0415
    rho, p = spearmanr(off, on)
    print(f"  rank correlation    : rho {rho:+.4f}  (p {p:.2g})")

    for label, arr in (("OFF", off), ("ON", on)):
        nrm = np.array([v for c, v in zip(cases, arr) if is_normal(c)])
        dis = np.array([v for c, v in zip(cases, arr) if not is_normal(c)])
        if len(nrm) and len(dis):
            band = nrm / nrm.mean()
            print(f"  {label:<4} normals {nrm.mean():.6f} (band "
                  f"{band.min():.2f}-{band.max():.2f}, CV {nrm.std(ddof=1) / nrm.mean():.3f})"
                  f"  diseased {np.median(dis) / np.median(nrm):.2f}x normals")


if __name__ == "__main__":
    if "--report-only" in sys.argv:
        report(OUT)
    else:
        main()
