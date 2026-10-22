"""Do the zones tile the aortic trunk exactly once? Checked on the centerline.

    .venv\\Scripts\\python tools\\zone_tiling.py                  # thoracic
    .venv\\Scripts\\python tools\\zone_tiling.py Y:/aorta/abdominal/MK/JK

The surface is the wrong place to ask this. A face's zone depends on the station field,
which is a projection and has its own failure modes -- the arch hairpin above all -- so
a fault there could be the zones or could be the projection, and the two are not
separable by looking.

The CENTERLINE is where the question is clean. Every node of the aortic trunk has one
arc-length station and each zone is an interval, so "how many zones claim this node"
is arithmetic with nothing in between. Every node from the first zone's start to the
last zone's end must be claimed exactly once:

    claimed twice   the intervals overlap -- two zones would both contain the same
                    stretch of aorta, and any measurement of either includes the other
    claimed never   the intervals leave a hole, so a stretch of aorta belongs to no
                    zone and would be missed by a sweep that trusted the set to cover

A hole is expected in exactly one place and is not a fault: where a zone is not offered
at all, because its landmarks are missing or coincide. Those are reported separately
from holes inside the offered set.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import branch_cut as BC                                      # noqa: E402
from anat_bench import find_cases                            # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
TILING = ("Ishimaru zone 0", "Ishimaru zone 1", "Ishimaru zone 2",
          "Ishimaru zone 3", "Ishimaru zone 4")
STEP_MM = 0.25


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    cases = find_cases(args.cohort)
    if args.limit:
        cases = cases[:args.limit]

    n_ok = n_over = n_hole = 0
    print(f"{len(cases)} case(s), sampling the trunk every {STEP_MM} mm\n")
    for cid, cl, mp in cases:
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            route = BC.aorta_route(s.tree)
            zs = sorted((p for p in s.presets() if p["name"] in TILING),
                        key=lambda p: p["s0"])
            if not zs:
                print(f"  {cid:12s} no tiling zones offered")
                continue
            L = float(route.length)
            x = np.arange(0.0, L, STEP_MM)
            claims = np.zeros(len(x), int)
            for p in zs:
                lo = float(p["s0"])
                hi = L if p["s1"] is None else float(p["s1"])
                claims += ((x >= lo) & (x < hi)).astype(int)

            lo0 = min(float(p["s0"]) for p in zs)
            hi0 = max(L if p["s1"] is None else float(p["s1"]) for p in zs)
            span = (x >= lo0) & (x < hi0)
            over = int((claims[span] > 1).sum())
            hole = int((claims[span] == 0).sum())
            n_over += bool(over)
            n_hole += bool(hole)
            n_ok += (not over) and (not hole)

            tag = ""
            if over:
                tag += f"   OVERLAP {over * STEP_MM:.1f} mm"
            if hole:
                # Where -- a hole between two offered zones is a different thing from
                # the stretch a missing zone would have covered.
                gaps = []
                d = np.diff(np.concatenate(([0], (claims == 0) & span, [0])).astype(int))
                for a, b in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]):
                    gaps.append(f"{x[a]:.0f}-{x[min(b, len(x) - 1)]:.0f}")
                tag += f"   HOLE {hole * STEP_MM:.1f} mm at {', '.join(gaps[:3])}"
            missing = [z for z in TILING if z not in [p["name"] for p in zs]]
            print(f"  {cid:12s} {len(zs)} zone(s) covering "
                  f"{lo0:6.1f}-{hi0:6.1f} of {L:6.1f} mm{tag}"
                  + (f"   [not offered: "
                     f"{', '.join(m.replace('Ishimaru ', '') for m in missing)}]"
                     if missing else ""))
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:12s} ERROR {type(e).__name__}: {e}")

    print(f"\n  {n_ok} of {len(cases)} tile their covered span exactly once")
    print(f"  {n_over} with an overlap, {n_hole} with a hole inside the covered span")


if __name__ == "__main__":
    main()
