"""Cut a plane normal to the trunk at every node. It must touch exactly one zone.

    .venv\\Scripts\\python tools\\zone_rings.py                   # thoracic
    .venv\\Scripts\\python tools\\zone_rings.py --limit 5

tools/zone_tiling.py checks the same property in one dimension -- do the intervals
overlap -- and that is arithmetic, so it can only ever catch a mistake in the interval
arithmetic. This checks the SURFACE. At each node of the aortic trunk, take the ring of
wall the normal plane there passes through and count the distinct zones it belongs to.
One is correct. Two means the wall at one station is split between zones, which the
interval test cannot see, because it happens when the STATION FIELD assigns part of a
circumference differently from the rest -- the arch hairpin above all, where the far
limb of the fold lies a few millimetres away and can capture wall that is not its own.

THE FAR LIMB IS EXCLUDED GEOMETRICALLY, not by station. A plane through the arch really
does pass through the descending aorta as well, and that surface really does belong to
a different zone, so counting it would report a fault that is not one. The ring is
taken as the connected piece of the plane's intersection nearest the centerline point,
found by proximity between the candidate vertices themselves -- which never consults
the field being tested.
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

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
STEP_MM = 1.0                 # how often along the trunk to cut a plane
# The slab has to be thicker than the mesh is coarse, or it falls between the
# vertices and the node is skipped. Fixed at 0.8 mm it found a ring on 4 nodes of
# 380 on KK15_1 and on 366 of 366 on DJ18_2 -- the same check, silently measuring
# almost nothing on half the cohort. Scaled to the mesh's own edge length instead.
SLAB_EDGES = 1.1
LINK_EDGES = 3.0
# How close to a zone boundary a split stops being a fault. Two zones meet along a
# level set of the station field, which is near the normal plane but not identical to
# it, so the rings within a whisker of a boundary legitimately see both sides.
BOUNDARY_GRACE_MM = 3.0


def ring_labels(V, label, origin, normal, tree, slab, link):
    """Zone labels on the near ring of the plane at *origin*, or None if too thin."""
    near = tree.query_ball_point(origin, r=60.0)
    if len(near) < 12:
        return None
    idx = np.asarray(near)
    rel = V[idx] - origin
    along = rel @ normal
    band = idx[np.abs(along) < slab]
    if len(band) < 8:
        return None
    pts = V[band]
    pairs = cKDTree(pts).query_pairs(link, output_type="ndarray")
    if len(pairs) == 0:
        return None
    g = sp.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                      shape=(len(pts), len(pts)))
    k, lab = connected_components(g, directed=False)
    # The piece whose closest point to the centerline node is nearest: the far limb of
    # a fold is a different piece and sits further away.
    best, bestd = 0, np.inf
    for c in range(k):
        d = np.linalg.norm(pts[lab == c] - origin, axis=1).min()
        if d < bestd:
            best, bestd = c, d
    sel = band[lab == best]
    if len(sel) < 8:
        return None
    return label[sel]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    cases = find_cases(args.cohort)
    if args.limit:
        cases = cases[:args.limit]

    clean = 0
    worst = []
    for cid, cl, mp in cases:
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            mesh = s.mesh
            route = BC.aorta_route(s.tree)
            L = float(route.length)
            station, _own, on_route = RS.station_field(mesh, s.tree, route)
            zs = sorted((p for p in s.presets() if p["name"] in TILING),
                        key=lambda p: p["s0"])
            V = np.asarray(mesh.vertices, float)
            label = -np.ones(len(V), int)
            for zi, p in enumerate(zs):
                f = RS.zone_field(station, on_route, p["s0"], p["s1"], route_end=L)
                lo = p["s0"]
                hi = p["s1"] if p["s1"] is not None else np.inf
                label[(f >= lo) & (f < hi)] = zi
            label, _moved = RS.enforce_connected_labels(mesh, label)

            e = mesh.edges_unique
            edge_len = float(np.median(np.linalg.norm(
                V[e[:, 0]] - V[e[:, 1]], axis=1)))
            slab = max(0.8, SLAB_EDGES * edge_len)
            link = max(2.5, LINK_EDGES * edge_len)
            bounds = sorted({float(p["s0"]) for p in zs}
                            | {float(p["s1"]) for p in zs if p["s1"] is not None})
            tree = cKDTree(V)
            checked = split = 0
            offenders = defaultdict(float)
            for st in np.arange(1.0, L - 1.0, STEP_MM):
                name, local, _fwd = route.locate(float(st))
                segment = s.tree.segments[name]
                origin, tangent = segment.at_arclength(float(local))
                lb = ring_labels(V, label, np.asarray(origin, float),
                                 np.asarray(tangent, float), tree, slab, link)
                if lb is None:
                    continue
                checked += 1
                u = set(int(x) for x in np.unique(lb) if x >= 0)
                near_bound = any(abs(st - b) <= BOUNDARY_GRACE_MM for b in bounds)
                if len(u) > 1 and not near_bound:
                    split += 1
                    offenders[tuple(sorted(u))] += STEP_MM
            frac = 100.0 * split / max(checked, 1)
            clean += (split == 0)
            tag = ""
            if split:
                worst.append((frac, cid, dict(offenders)))
                pieces = ", ".join(
                    f"{'+'.join(zs[i]['name'][-6:] for i in k)} {v:.0f} mm"
                    for k, v in sorted(offenders.items(), key=lambda t: -t[1])[:3])
                tag = f"   SPLIT at {split} of {checked} nodes ({frac:.1f}%): {pieces}"
            print(f"  {cid:12s} {checked:4d} of {int((L - 2) / STEP_MM)}"
                  f" nodes checked (edge {edge_len:.2f} mm){tag}")
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:12s} ERROR {type(e).__name__}: {e}")

    print(f"\n  {clean} of {len(cases)} have every node's ring in exactly one zone")
    for frac, cid, off in sorted(worst, reverse=True)[:8]:
        print(f"    {cid:12s} {frac:5.1f}% of nodes split")


if __name__ == "__main__":
    main()
