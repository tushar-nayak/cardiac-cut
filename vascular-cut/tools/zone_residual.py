"""What streak is LEFT after each stub is given to its ostium, and is it really stub?

    .venv\\Scripts\\python tools\\zone_residual.py

The crease region stops at the mouth, so it does not cover a collar the older frontier
flood would have taken. Scored on the frontier's wider notion of a stub, 27 of 89 stubs
still look split after the fix -- but that set runs PAST the mouth onto the aorta, and
wall out there is supposed to carry the next zone. Counting it as a streak would be
counting the fix against the very over-reach it exists to avoid.

So the residue is split three ways, using the aorta's own radius profile (read off its
cross sections, so branches cannot inflate it) to say which side of the mouth a vertex
is on:

  ON THE AORTA     ratio <= 1.0, inside the aortic envelope. Wall. Carrying the next
                   zone there is correct, not a streak.
  OFF THE AORTA    ratio > 1.0 and disagreeing with its stub's zone. A real streak the
                   crease region failed to reach, and the only number here worth acting
                   on.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
from zone_stub_bench import label_case, own_by_stub          # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"


def ratio_field(mesh, route):
    ss, rr = RS.aortic_radius_profile(mesh, route)
    P = np.array([route.point(float(t)) for t in ss])
    d, i = cKDTree(P).query(np.asarray(mesh.vertices, float))
    return d / np.maximum(rr[i], 1e-9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    cases = find_cases(args.cohort)
    if args.only:
        cases = [c for c in cases if c[0] in set(args.only.split(","))]

    tot_wall = tot_off = tot_stub = 0
    worst = []
    for cid, cl, mp in cases:
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            mesh, route, _z, lab, isb, adj, _f = label_case(s)
            reg = RS.stub_regions(mesh, route, adj=adj)
            new = own_by_stub(mesh, lab, reg, isb, adj)
            fr = RS.branch_regions(mesh, route)
            ratio = ratio_field(mesh, route)
            n_wall = n_off = n_stub = 0
            for j in range(int(fr.max()) + 1):
                sel = fr == j
                v = new[sel]
                v = v[v >= 0]
                if len(v) < 10:
                    continue
                n_stub += 1
                z = int(np.bincount(v).argmax())
                res = sel & (new != z) & (new >= 0)
                off = int((res & (ratio > 1.0)).sum())
                n_wall += int((res & (ratio <= 1.0)).sum())
                n_off += off
                if off:
                    worst.append((off, cid, j))
            tot_wall += n_wall
            tot_off += n_off
            tot_stub += n_stub
            print(f"  {cid:10s} residue off the aorta {n_off:4d}   "
                  f"on the aorta {n_wall:4d}   ({n_stub} stub(s))", flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:10s} ERROR {type(e).__name__}: {e}", flush=True)

    print(f"\n  {tot_stub} stubs")
    print(f"  residue OFF the aorta (a real streak)  {tot_off:5d} vertices")
    print(f"  residue ON the aorta (correct)         {tot_wall:5d} vertices")
    for n, cid, j in sorted(worst, reverse=True)[:10]:
        print(f"    {cid:10s} stub {j}: {n} vertices off the aorta disagree")


if __name__ == "__main__":
    main()
