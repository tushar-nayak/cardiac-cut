"""Does zone 3 begin at the left subclavian, or is there aorta in between?

    .venv\\Scripts\\python tools\\zone3_start.py
    .venv\\Scripts\\python tools\\zone3_start.py --only KK32_0

Zone 3 runs from the left subclavian to 20 mm beyond it, so its proximal boundary IS
that vessel's ostium. The preset says so exactly -- no clearance offset is applied -- so
if zone 3 looks late, the boundary is landing distal of the mouth rather than the number
being wrong.

MEASURED ALONG THE WALL, IN MILLIMETRES. The first version of this compared the
boundary's harmonic coordinate with the ostium rim's and converted the difference
through the station curve. It reported a median gap of 0.5 mm and flagged six cases;
asked to mark the bad ones by eye, the reviewer marked twenty, and barely any were the
six. The metric was wrong in a way worth recording: a coordinate difference run through
a calibration says little about a distance on a bent, dilated arch, which is exactly
where the defect lives.

So measure the thing being looked at -- the WALL between the subclavian's mouth and the
start of zone 3. Geodesic distance across the surface from the ostium rim to the zone 2
/ zone 3 boundary, over wall vertices only and forbidden to cut through the stub, so it
runs along the aorta rather than up and over the vessel. Zone 3 is meant to begin
immediately after the vessel, so anything much beyond a millimetre or two is zone 2
holding aorta that belongs to zone 3.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
from zone_stub_bench import label_case, ostium_zones         # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"


def wall_gap(mesh, adj, reg, lab, lsa_regions):
    """(median, nearest) mm of wall between the subclavian's mouth and zone 3."""
    e = np.asarray(mesh.edges_unique, int)
    sel = np.isin(reg, lsa_regions)
    if not sel.any():
        return None
    cut = sel[e[:, 0]] != sel[e[:, 1]]
    rim = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
    rim_wall = rim[~sel[rim]]
    on = (((lab[e[:, 0]] == 2) & (lab[e[:, 1]] == 3))
          | ((lab[e[:, 0]] == 3) & (lab[e[:, 1]] == 2)))
    ring = np.unique(e[on])
    ring = ring[reg[ring] < 0]
    if len(rim_wall) < 6 or len(ring) < 8:
        return None
    keep = np.flatnonzero(reg < 0)
    pos = -np.ones(len(reg), int)
    pos[keep] = np.arange(len(keep))
    seed = pos[rim_wall]
    seed = seed[seed >= 0]
    tgt = pos[ring]
    tgt = tgt[tgt >= 0]
    if not len(seed) or not len(tgt):
        return None
    d = dijkstra(adj[keep][:, keep], directed=False, indices=seed, min_only=True)
    dd = d[tgt]
    dd = dd[np.isfinite(dd)]
    if not len(dd):
        return None
    # THE CLOSEST APPROACH IS THE SLACK. A boundary is a ring around the whole vessel
    # and the mouth is a patch on one side of it, so most of the ring is half a
    # circumference away no matter where the boundary sits -- which is why the median
    # of this distance says almost nothing. What says something is how much wall lies
    # between the mouth and the NEAREST part of the ring: that is the gap zone 3 has to
    # cross before it starts, and it is what a nudge has to close.
    return float(np.percentile(dd, 10)), float(dd.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    cases = find_cases(args.cohort)
    if args.only:
        cases = [c for c in cases if c[0] in set(args.only.split(","))]

    rows = []
    for cid, cl, mp in cases:
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            mesh, route, zones, lab, isb, adj, _f = label_case(s)
            reg = RS.stub_regions(mesh, route, adj=adj)
            if reg.max() < 0:
                continue
            oz = ostium_zones(s, mesh, route, reg)
            lab = RS.own_stubs(mesh, lab, reg, isb, adj=adj, region_zone=oz)
            znum = [int(p["name"].split()[-1]) for p in zones]
            if 2 not in znum or 3 not in znum:
                print(f"  {cid:10s} no zone 2/3 pair")
                continue
            lsa = [j for j, z in oz.items() if z == 2]
            if not lsa:
                print(f"  {cid:10s} no stub owned by zone 2")
                continue
            got = wall_gap(mesh, adj, reg, lab, lsa)
            if got is None:
                print(f"  {cid:10s} could not measure")
                continue
            med, near = got
            rows.append((med, near, cid))
            print(f"  {cid:10s} slack between the subclavian and zone 3: "
                  f"p10 {med:5.1f} mm   nearest {near:5.1f} mm", flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:10s} ERROR {type(e).__name__}: {e}", flush=True)

    if not rows:
        return
    g = np.array([r[0] for r in rows])
    print(f"\n  {len(rows)} cases: median {np.median(g):.1f} mm of wall between the "
          f"subclavian and zone 3  (p90 {np.percentile(g, 90):.1f}, max {g.max():.1f})")
    print("  worst")
    for med, near, cid in sorted(rows, reverse=True)[:14]:
        print(f"    {cid:10s} slack {med:5.1f} mm  (nearest {near:5.1f} mm)")


if __name__ == "__main__":
    main()
