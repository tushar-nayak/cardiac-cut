"""How ragged is each zone boundary -- not how tilted, how WANDERING.

    .venv\\Scripts\\python tools\\zone_ragged.py
    .venv\\Scripts\\python tools\\zone_ragged.py --only KR47_1

Flatness -- distance off the best-fit plane -- says nothing about this. A boundary can
sit dead flat and still zigzag across that plane, and a boundary that zigzags is what
"zone 0 looks glitchy" means when you see it: lobes of one colour interlocking with the
next. On KR47_1 the 0/1 boundary is 2.64 mm off plane, better than KK15_1's 2.92, and
it is the one that looks wrong.

So measure the boundary against the shortest ring that could replace it. The cut length
is the total length of the mesh edges the boundary crosses; the ring is the circle of
the same enclosed girth, 2 pi R with R the median in-plane distance from the boundary's
own centroid. A smooth boundary crossing a triangulated tube gives a constant a little
over 1 -- the crossings zigzag at the scale of one triangle, and nothing can be done
about that. Anything well above it is wandering at a scale bigger than the mesh.

Reported with the boundary's flatness beside it, because the two are independent and a
change that improves one usually costs the other: fitting the full ring took flatness
from 2.35 to 3.20 median, and majority-label smoothing took raggedness from 1.63 to
1.59 while flatness went from 1.89 to 11.12.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
from zone_stub_bench import label_case                       # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"


def ragged(mesh, lab, zi, zj, frozen=None):
    """(raggedness, flatness_mm, n_vertices) for the boundary between two zones.

    ON WALL EDGES ONLY, matching flatten_boundary. A stub that has taken its ostium's
    zone puts its mouth on the boundary ring, and those vertices sit off in the
    branch's direction: including them drags the fitted plane and reports a boundary
    as bent when what is bent is the vessel leaving it. Measured on KR47_1, that alone
    is the difference between 13.3 mm and 2.6 mm.
    """
    V = np.asarray(mesh.vertices, float)
    e = np.asarray(mesh.edges_unique, int)
    if frozen is not None:
        e = e[~frozen[e[:, 0]] & ~frozen[e[:, 1]]]
    a, b = lab[e[:, 0]], lab[e[:, 1]]
    cross = ((a == zi) & (b == zj)) | ((a == zj) & (b == zi))
    if cross.sum() < 6:
        return None
    cut_len = float(np.linalg.norm(V[e[cross, 0]] - V[e[cross, 1]], axis=1).sum())
    ring = np.unique(np.concatenate([e[cross, 0], e[cross, 1]]))
    P = V[ring]
    c = P.mean(axis=0)
    # Plane of the ring: its smallest principal direction.
    _u, _s, vt = np.linalg.svd(P - c, full_matrices=False)
    n = vt[2]
    off = (P - c) @ n
    flat = float(np.abs(off).max())
    inplane = (P - c) - np.outer(off, n)
    R = float(np.median(np.linalg.norm(inplane, axis=1)))
    if R <= 1e-6:
        return None
    return cut_len / (2.0 * np.pi * R), flat, len(ring)


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
            if reg.max() >= 0:
                lab = RS.own_stubs(mesh, lab, reg, isb, adj=adj)
            out = []
            for zi in range(len(zones) - 1):
                r = ragged(mesh, lab, zi, zi + 1, frozen=(reg >= 0))
                if r is None:
                    out.append("     -")
                    continue
                rows.append((r[0], cid, zi, r[1]))
                out.append(f"{r[0]:5.2f}")
            print(f"  {cid:10s} " + "  ".join(
                f"{zi}/{zi + 1} {v}" for zi, v in enumerate(out)), flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:10s} ERROR {type(e).__name__}: {e}", flush=True)

    if not rows:
        return
    vals = np.array([r[0] for r in rows])
    print(f"\n  {len(rows)} boundaries: raggedness median {np.median(vals):.2f}  "
          f"p90 {np.percentile(vals, 90):.2f}  max {vals.max():.2f}")
    print("  worst")
    for v, cid, zi, fl in sorted(rows, reverse=True)[:10]:
        print(f"    {cid:10s} boundary {zi}/{zi + 1}: raggedness {v:5.2f}  "
              f"flatness {fl:5.2f} mm")


if __name__ == "__main__":
    main()
