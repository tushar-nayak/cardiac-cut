"""Does giving each stub to its own ostium fix the streak without moving the trunk?

    .venv\\Scripts\\python tools\\zone_stub_bench.py
    .venv\\Scripts\\python tools\\zone_stub_bench.py --alpha 1.2,1.5,2.0

THE DEFECT. A great vessel comes out one colour for almost its whole length and then
carries a streak of the next zone down one flank. It happens where the stub curls back
and hovers over the next zone: the wall rule labels by the harmonic coordinate, and on
the part of the stub that the branch mask missed, that coordinate has already moved on.

WHY THIS IS MEASURED HERE AND NOT IN zone_map. Every previous attempt was reported with
the branch mask that was itself under test, and one of those reports -- "OPT 3 changes 0
wall vertices" -- was false: an independent check found 168 of 193 changed vertices on
KR6_1 to be aortic. So the two numbers here come from different places on purpose:

  TRUNK DAMAGE   counted with geometry the stub finder never sees. A vertex nearer the
                 aortic centerline than 0.95 x the local aortic radius is aortic wall,
                 whatever else it is. That is one-sided -- a stub vertex near its own
                 mouth can sit inside that envelope -- so it undercounts, and every
                 vertex it reports is real.
  STREAK         counted per stub as the share of it not carrying the stub's own
                 majority zone. Reported over BOTH stub definitions, the concavity
                 regions and the frontier flood, because a fix scored on its own notion
                 of where the branch is can only agree with itself.

Zone shape is guarded too: pieces per zone, and boundary flatness, before and after.
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import branch_cut as BC                                      # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
# The labelling itself now lives in region_station, so that the web app and the
# desktop tool can import it without putting this bench directory on their path.
# Re-exported here because a dozen bench scripts and two test suites import it
# from this module by name. Setting a flag on THIS module no longer reaches the
# code that reads it -- set it on region_station instead.
from region_station import (                                # noqa: E402,F401
    TILING,
    WALL_K,
    CROSS_SECTION_PLACEMENT,
    CONTOUR_TOL_MM,
    RIM_PLACEMENT,
    LANDMARK_LAG_FIX,
    LANDMARK_LAG_MIN_MM,
    PLANE_BOUNDARY,
    ZONE3_LENGTH_MM,
    wall_mask,
    label_case,
    own_by_stub,
    ostium_zones,
    final_labels,
    plane_boundaries,
)


def streak(lab, regions):
    """Per stub: how much of it does not carry the stub's own majority zone."""
    n_split, n_stub, frac = 0, 0, []
    for j in range(int(regions.max()) + 1):
        sel = regions == j
        v = lab[sel]
        v = v[v >= 0]
        if len(v) < 10:
            continue
        n_stub += 1
        share = 1.0 - np.bincount(v).max() / len(v)
        frac.append(float(share))
        if share > 0.02:
            n_split += 1
    return n_stub, n_split, frac


def pieces(mesh, lab):
    from scipy.sparse.csgraph import connected_components
    adj = RS.surface_graph(mesh)
    tot = 0
    for zi in range(int(lab.max()) + 1):
        idx = np.flatnonzero(lab == zi)
        if len(idx) < 2:
            continue
        n, _ = connected_components(adj[idx][:, idx], directed=False)
        tot += int(n)
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--alpha", default="1.5")
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    alphas = [float(x) for x in args.alpha.split(",")]

    cases = find_cases(args.cohort)
    if args.only:
        cases = [c for c in cases if c[0] in set(args.only.split(","))]

    agg = {a: dict(chg=0, wall=0, stubs=0, split0=0, split1=0,
                   fsplit0=0, fsplit1=0, fstubs=0, pc0=0, pc1=0) for a in alphas}
    for cid, cl, mp in cases:
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            mesh, route, zones, lab, isb, adj, flat = label_case(s)
            wall = wall_mask(mesh, route)
            fr = RS.branch_regions(mesh, route)
            p0 = pieces(mesh, lab)
            line = f"  {cid:10s}"
            for a in alphas:
                reg = RS.stub_regions(mesh, route, alpha=a, adj=adj)
                # The ostium decides, as it does in tools/zone_map.py. Measuring the
                # proxy instead would be measuring a program that does not ship.
                new = own_by_stub(mesh, lab, reg, isb, adj,
                                  region_zone=ostium_zones(s, mesh, route, reg))
                chg = new != lab
                ns, s0, _ = streak(lab, reg)
                _n, s1, _ = streak(new, reg)
                fn, f0, _ = streak(lab, fr)
                _f, f1, _ = streak(new, fr)
                p1 = pieces(mesh, new)
                d = agg[a]
                d["chg"] += int(chg.sum())
                d["wall"] += int((chg & wall).sum())
                d["stubs"] += ns
                d["split0"] += s0
                d["split1"] += s1
                d["fstubs"] += fn
                d["fsplit0"] += f0
                d["fsplit1"] += f1
                d["pc0"] += p0
                d["pc1"] += p1
                line += (f"  a{a}: chg {int(chg.sum()):4d} wall {int((chg & wall).sum()):3d}"
                         f" split {s0}->{s1} pieces {p0}->{p1}")
            print(line, flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:10s} ERROR {type(e).__name__}: {e}", flush=True)

    print()
    for a in alphas:
        d = agg[a]
        print(f"  alpha {a}")
        print(f"    labels changed        {d['chg']:6d}  of which aortic wall {d['wall']:5d}")
        print(f"    stubs split, own set  {d['split0']:4d} -> {d['split1']:4d}"
              f"  of {d['stubs']}")
        print(f"    stubs split, frontier {d['fsplit0']:4d} -> {d['fsplit1']:4d}"
              f"  of {d['fstubs']}")
        print(f"    zone pieces           {d['pc0']:5d} -> {d['pc1']:5d}")


if __name__ == "__main__":
    main()
