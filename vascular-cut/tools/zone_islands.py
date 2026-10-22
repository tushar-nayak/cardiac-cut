"""No zone may sit as an island inside another zone.

    .venv\\Scripts\\python tools\\zone_islands.py
    .venv\\Scripts\\python tools\\zone_islands.py --only KR47_1

A zone is a stretch of aorta. A patch of zone 1 marooned in the middle of zone 0 is not
a stretch of anything -- it is a labelling artefact, and any area or delta K taken on
the zone averages across it. So each zone's connected pieces are counted, and for every
piece that is not the largest, the zone that SURROUNDS it is named: that is the pair
that has to be explained, and it is what makes the difference between "a branch whose
mouth straddles a boundary", which is anatomy, and "eight vertices left behind", which
is not.

Reported per case so a defect can be looked at, and totalled so a change can be judged.
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
from zone_stub_bench import final_labels                     # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"


def islands(mesh, lab):
    """Every non-largest piece of every zone, with the zone that surrounds it.

    Returns a list of (zone, size, surrounding_zone_or_None).
    """
    e = np.asarray(mesh.edges_unique, int)
    nv = len(lab)
    same = lab[e[:, 0]] == lab[e[:, 1]]
    ee = e[same]
    g = sp.coo_matrix((np.ones(len(ee)), (ee[:, 0], ee[:, 1])), shape=(nv, nv))
    k, comp = connected_components(g, directed=False)
    sizes = np.bincount(comp, minlength=k)
    best = {}
    for c in range(k):
        v = int(np.flatnonzero(comp == c)[0])
        L = int(lab[v])
        if L < 0:
            continue
        if L not in best or sizes[c] > sizes[best[L]]:
            best[L] = c
    out = []
    for c in range(k):
        idx = np.flatnonzero(comp == c)
        L = int(lab[idx[0]])
        if L < 0 or best.get(L) == c:
            continue
        on = np.zeros(nv, bool)
        on[idx] = True
        touch = Counter()
        for u, w in ((e[:, 0], e[:, 1]), (e[:, 1], e[:, 0])):
            take = on[u] & ~on[w] & (lab[w] >= 0)
            touch.update(lab[w[take]].tolist())
        around = touch.most_common(1)[0][0] if touch else None
        out.append((L, int(len(idx)), None if around is None else int(around)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--only", default="")
    ap.add_argument("--no-stubs", action="store_true",
                    help="label without stub ownership, to attribute a defect")
    args = ap.parse_args()
    cases = find_cases(args.cohort)
    if args.only:
        cases = [c for c in cases if c[0] in set(args.only.split(","))]

    n_bad = n_isl = 0
    for cid, cl, mp in cases:
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            # THROUGH final_labels, which is the whole labelling the tool ships.
            # Rebuilding the last few steps here is how this checker once reported an
            # island on KY160_1 that does not exist in the map.
            mesh, route, zones, lab, isb, adj, _f, reg, oz = final_labels(s)
            isl = islands(mesh, lab)
            if isl:
                n_bad += 1
                n_isl += len(isl)
                for z, n, around in isl:
                    print(f"  {cid:10s} zone {z} island of {n:5d} vertices "
                          f"inside zone {around}", flush=True)
            else:
                print(f"  {cid:10s} clean, {len(zones)} zone(s)", flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:10s} ERROR {type(e).__name__}: {e}", flush=True)
    print(f"\n  {n_isl} island(s) on {n_bad} of {len(cases)} case(s)")


if __name__ == "__main__":
    main()
