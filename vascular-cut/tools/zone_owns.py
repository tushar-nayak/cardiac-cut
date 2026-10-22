"""What branches does each zone own, across a cohort?

    .venv\\Scripts\\python tools\\zone_owns.py                    # thoracic
    .venv\\Scripts\\python tools\\zone_owns.py Y:/aorta/abdominal/MK/JK

A zone owns the vessel it is NAMED for, and that is the whole of what it should own:

    zone 0   the brachiocephalic, and the daughters that share its ostium
    zone 1   the left common carotid
    zone 2   the left subclavian
    zone 3   NOTHING -- its distal bound is a distance, 20 mm past the subclavian,
             not a vessel, so no ostium should fall in it
    zone 4   nothing, unless the scan traced an intercostal

So zone 3 and zone 4 are the sharp test: anything they own is either a vessel that
genuinely arises there and was traced, or an ostium in the wrong place. Both are worth
knowing about and they look identical in a count, so the names are printed.

Membership comes from region_station.branch_inside -- the same function the cut uses.
A checker with its own copy of the rule can only ever agree with itself.
"""
import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import branch_cut as BC                                      # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
TILING = ("Ishimaru zone 0", "Ishimaru zone 1", "Ishimaru zone 2",
          "Ishimaru zone 3", "Ishimaru zone 4")
EXPECT_EMPTY = ("Ishimaru zone 3", "Ishimaru zone 4")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    cases = find_cases(args.cohort)
    if args.limit:
        cases = cases[:args.limit]
    owns = defaultdict(list)

    for cid, cl, mp in cases:
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            route = BC.aorta_route(s.tree)
            jg = BC._junction_graph(s.tree)
            L = float(route.length)
            zs = sorted((p for p in s.presets() if p["name"] in TILING),
                        key=lambda p: p["s0"])
            _p, _t, sta, owner, _r = RS._segment_station_samples(s.tree, route, jg)
            seg_st = {}
            for nm, v in zip(owner, sta):
                if nm not in route.chain:
                    seg_st.setdefault(nm, float(v))
            land = BC.aorta_landmark_stations(s.tree, route, jg)
            for p in zs:
                segs = [nm for nm, v in seg_st.items()
                        if RS.branch_inside(np.array([v]), p["s0"], p["s1"], L)[0]]
                names = sorted({k for k, v in land.items()
                                if RS.branch_inside(np.array([float(v)]), p["s0"],
                                                    p["s1"], L)[0]})
                owns[p["name"]].append((cid, len(segs), names))
            print(f"  {cid:12s} " + "  ".join(
                f"{p['name'][-6:]}:{len([1 for nm, v in seg_st.items() if RS.branch_inside(np.array([v]), p['s0'], p['s1'], L)[0]])}"
                for p in zs))
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:12s} ERROR {type(e).__name__}: {e}")

    print("\n  zone                offered  owns nothing  usual contents")
    for z in TILING:
        rows = owns.get(z, [])
        if not rows:
            continue
        empty = sum(1 for _, n, _ in rows if n == 0)
        common = Counter(nm for _, _, names in rows for nm in names)
        print(f"    {z:18s} {len(rows):5d} {empty:11d}  "
              f"{', '.join(n for n, _ in common.most_common(3)) or '-'}")

    print("\n  cases where a zone owns something it should not")
    bad = 0
    for z in EXPECT_EMPTY:
        for cid, n, names in owns.get(z, []):
            if n:
                bad += 1
                print(f"    {cid:12s} {z:18s} {n} segment(s): "
                      f"{', '.join(names) or 'unnamed'}")
    if not bad:
        print("    none - zones 3 and 4 own no branch on any case")


if __name__ == "__main__":
    main()
