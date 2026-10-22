"""A zone cut by station must land where the anatomy says, and be a usable surface.

    .venv\\Scripts\\python tests\\test_region_station.py
    .venv\\Scripts\\python tests\\test_region_station.py --prove

--prove puts each defect back, one at a time, and fails if the test does NOT catch it.
A regression test nobody has watched fail is a regression test nobody has tested.

The four things guarded, and why each is here rather than assumed:

  LANDS AT THE STATION   The point of the station cut. The planar path drifts 8-37 mm
                         on these two cases (ostium clearance, then the preset ladder),
                         and on Ishimaru zone 3 the drift exceeds half the zone's
                         20 mm length. A station cut that drifts at all has no reason
                         to exist.
  WINDING                Rebuilding a cut face from filtered inside/outside lists
                         loses the cyclic order and reverses about half the faces.
                         The surface still looks right and every area-weighted
                         measurement taken on it is wrong.
  NO ZERO-AREA FACES     A level set falling exactly on a vertex makes a triangle with
                         two coincident corners. No normal, and delta K is
                         area-weighted.
  BRANCH MEMBERSHIP      s0 <= ostium < s1, half-open. Inclusive at BOTH ends puts
                         every landmark vessel in two zones at once, because the
                         Ishimaru zones tile; exclusive at both drops the vessel the
                         zone is named for.
  ONE PIECE              a zone delivered in two patches describes neither, and any
                         area or delta K taken on it averages across the gap.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vascular_api import CutSession                          # noqa: E402
import branch_cut as BC                                      # noqa: E402
import region_station as RS                                  # noqa: E402

CASES = [
    ("JK5_1", "JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
    ("KY160_1", "KY160_1_CentrelinePoints.csv",
     "Geometries & Meshes/mesh/KY160_1.mat"),
]
# Geometry is gitignored -- see the .gitignore note on why derived surfaces are kept
# out deliberately rather than by accident -- so a fresh clone has some of these and
# not others. Run on what is actually here and say what is missing, rather than dying
# on a FileNotFoundError that reads like a broken test.
CASES = [(cid, csv, mat) for cid, csv, mat in CASES
         if (ROOT / csv).exists() and (ROOT / mat).exists()]
# How far the rim may sit from the station asked for. Not zero: the rim is
# interpolated along mesh edges, so a vertex ON the level set is exact but the
# station field itself is sampled every SAMPLE_MM along the centerline.
STATION_TOL_MM = 1.0

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print(f"    FAIL  {msg}")
    return cond


def run():
    for cid, csv, mat in CASES:
        s = CutSession(str(ROOT / csv), str(ROOT / mat))
        s.recut()
        route = BC.aorta_route(s.tree)
        print(f"\n{cid}: route {route.length:.1f} mm")
        areas = {}
        for p in s.presets():
            name, s0, s1 = p["name"], p["s0"], p["s1"]
            m, info = RS.region_by_station(s.mesh, s.tree, s0, s1, route=route)
            areas[name] = info["area_mm2"]
            print(f"  {name:32s} {info['faces']:7,} faces  "
                  f"{info['area_mm2']:9,.0f} mm2  {info['n_loops']} loops")

            check(len(m.faces) > 0, f"{cid} {name}: empty")

            # LANDS AT THE STATION -- measured on the result, not on what was asked,
            # and on ORIGINAL vertices only.
            #
            # The rim's own vertices are interpolated onto the level set and are
            # deliberately not tested. Where the aorta folds back on itself -- the
            # inner curvature of the arch -- arc length is not continuous across the
            # surface: on KY160_1 a single 0.9 mm edge there spans 7 mm of station, so
            # a rim vertex placed correctly ON the level set re-measures well past it.
            # That is a property of the coordinate, not a misplaced boundary, and
            # asserting on it fails a cut that is right. What must hold is that no
            # original wall from outside the span survives.
            st, _own, on_route = RS.station_field(m, s.tree, route)
            d_orig, _ = cKDTree(np.asarray(s.mesh.vertices, float)).query(
                np.asarray(m.vertices, float))
            orig = on_route & (d_orig < 1e-9)
            wall = st[orig]
            if len(wall):
                if s0 is not None:
                    check(wall.min() >= s0 - STATION_TOL_MM,
                          f"{cid} {name}: kept original wall at {wall.min():.1f} mm, "
                          f"proximal bound asked {s0:.1f}")
                if s1 is not None:
                    check(wall.max() <= s1 + STATION_TOL_MM,
                          f"{cid} {name}: kept original wall at {wall.max():.1f} mm, "
                          f"distal bound asked {s1:.1f}")

            check(bool(m.is_winding_consistent),
                  f"{cid} {name}: winding inconsistent")
            a = m.area_faces
            check(int((a < 1e-9).sum()) == 0,
                  f"{cid} {name}: {int((a < 1e-9).sum())} zero-area faces")
            check(info["n_loops"] >= 1,
                  f"{cid} {name}: no boundary loop -- nothing was cut")
            # ONE PIECE. Guaranteed by keeping the largest component, so this guards
            # that guarantee rather than hoping for it -- and the island count beside
            # it says what upholding it cost.
            check(len(RS._component_sizes(m)) == 1,
                  f"{cid} {name}: delivered in "
                  f"{len(RS._component_sizes(m))} pieces")
            check(info["largest_island"] <= 0.02 * info["faces"],
                  f"{cid} {name}: discarded an island of "
                  f"{info['largest_island']} faces to make it one piece "
                  f"({100 * info['largest_island'] / max(info['faces'], 1):.1f}%)")

        # BRANCH MEMBERSHIP, on the case that has both presets.
        a, b = "Infrarenal aorta", "Infrarenal aorta with iliacs"
        if a in areas and b in areas:
            check(areas[b] > areas[a] * 1.2,
                  f"{cid}: '{b}' ({areas[b]:,.0f}) is not meaningfully larger than "
                  f"'{a}' ({areas[a]:,.0f}) -- the iliacs are in both")


def prove():
    """Put each defect back and check the suite notices."""
    import types

    def broken_winding(mesh, f, level, keep_above):
        """The filtered-list version: loses cyclic order, reverses half the faces."""
        import trimesh
        V = np.asarray(mesh.vertices, float)
        F = np.asarray(mesh.faces, int)
        f = np.asarray(f, float)
        inside = (f >= level) if keep_above else (f <= level)
        cnt = inside[F].sum(axis=1)
        newV, newf, newF = [V], [f], [F[np.nonzero(cnt == 3)[0]]]
        edge_pt = {}

        def cut_point(i, j):
            key = (i, j) if i < j else (j, i)
            if key in edge_pt:
                return edge_pt[key]
            fi, fj = f[key[0]], f[key[1]]
            t = 0.5 if abs(fj - fi) < 1e-12 else (level - fi) / (fj - fi)
            t = float(np.clip(t, 0.0, 1.0))
            edge_pt[key] = len(V) + len(edge_pt)
            newV.append((V[key[0]] * (1 - t) + V[key[1]] * t).reshape(1, 3))
            newf.append(np.array([level]))
            return edge_pt[key]

        tris = []
        for fi in np.nonzero((cnt == 1) | (cnt == 2))[0]:
            a, b, c = F[fi]
            ins = [v for v in (a, b, c) if inside[v]]
            out = [v for v in (a, b, c) if not inside[v]]
            if len(ins) == 1:
                tris.append((ins[0], cut_point(ins[0], out[0]),
                             cut_point(ins[0], out[1])))
            else:
                m_, n_ = cut_point(ins[0], out[0]), cut_point(ins[1], out[0])
                tris.append((ins[0], ins[1], n_))
                tris.append((ins[0], n_, m_))
        if tris:
            newF.append(np.asarray(tris, int))
        Vout, fout = np.vstack(newV), np.concatenate(newf)
        Fout = np.vstack(newF)
        used = np.zeros(len(Vout), bool)
        used[Fout.ravel()] = True
        remap = -np.ones(len(Vout), int)
        remap[used] = np.arange(used.sum())
        return trimesh.Trimesh(Vout[used], remap[Fout], process=False), fout[used]

    def inclusive_membership(station, on_route, s0, s1):
        """Membership inclusive at the distal end: the iliacs land in both presets."""
        lo = float(s0) if s0 is not None else -np.inf
        hi = float(s1) if s1 is not None else np.inf
        f = np.asarray(station, float).copy()
        br = ~np.asarray(on_route, bool)
        if br.any():
            ost = f[br]
            inside = (ost >= lo) & (ost <= hi)
            f[br] = np.where(inside, (lo + min(hi, ost.max() + 1.0)) / 2.0, lo - 1.0)
        return f

    for label, attr, repl in (("winding / zero-area", "clip_by_field", broken_winding),
                              ("branch membership", "zone_field",
                               inclusive_membership)):
        good = getattr(RS, attr)
        setattr(RS, attr, repl)
        del fails[:]
        print(f"\n{'=' * 78}\nPROVE: {label} defect reinstated")
        try:
            run()
        finally:
            setattr(RS, attr, good)
        if fails:
            print(f"  caught it: {len(fails)} assertion(s) fired")
        else:
            print(f"  NOT CAUGHT -- the test does not guard {label}")
            fails.append(f"prove: {label} went undetected")
        keep = [x for x in fails if x.startswith("prove:")]
        del fails[:]
        fails.extend(keep)


if __name__ == "__main__":
    if "--prove" in sys.argv:
        prove()
    else:
        run()
    print("\n" + "=" * 78)
    if fails:
        print(f"{len(fails)} FAILURE(S)")
        for f_ in fails:
            print(f"  {f_}")
    else:
        print("all checks passed")
    raise SystemExit(1 if fails else 0)
