"""A branch stub must be found from its own mouth, and must stay off the aorta.

    .venv\\Scripts\\python tests\\test_stub_regions.py
    .venv\\Scripts\\python tests\\test_stub_regions.py --prove

--prove puts each defect back, one at a time, and fails if the test does NOT catch it.
A regression test nobody has watched fail is a regression test nobody has tested.

WHAT IS BEING GUARDED. A great vessel that curls back over the arch used to come out
one zone for most of its length with a streak of the next zone down the flank that
hovers over the following zone. Two families of fix were measured and both failed:
distance from the aortic centerline (the populations overlap, d/r 1.05-1.45 for a
curled stub against 1.00 with p90 1.06 for wall) and frontier width (a pure trade-off,
418 wall vertices dragged with 56 stubs still split at one end of the sweep, 2934 with
4 still split at the other). What works is the crease at the mouth, and what it must
never do is what those two did:

  STAYS OFF THE AORTA   judged with geometry the finder never sees -- a vertex nearer
                        the aortic centerline than 0.95 x the local aortic radius is
                        aortic wall, the radius read off the vessel's own cross
                        sections so that branches cannot inflate it. This is the
                        assertion a previous attempt was reported as passing when it
                        did not: OPT 3 was announced as changing no wall, and an
                        independent check found 168 of 193 changed vertices on KR6_1
                        to be aortic.

                        The MOUTH ITSELF IS EXEMPT, and measurably has to be. An
                        ostium sits in a recess: on JK5_1 every vertex the raw test
                        flagged lay a median 0.0 mm and at most 1.5 mm from the
                        region's own boundary, at 0.78-0.95 of the median radius.
                        That is the fillet, which belongs to the branch. So the test
                        looks only DEEPER than one mouth radius inside the region.

                        NO LONGER FULLY INDEPENDENT, and it should be read that way.
                        stub_regions now stops a growth step whose new vertices have a
                        median radius ratio below 1.10, so finder and audit draw on the
                        same aortic_radius_profile. The assertion still catches a region
                        that runs somewhere the step rule did not veto, but it is no
                        longer the arms-length check it was. The guards that ARE
                        independent of the finder are the topological ones below --
                        bounded, one piece, no pockets, no zone islanded inside another
                        -- and, in the end, a person looking at 30 cases.
  BOUNDED               a stub is not a fifth of the aorta. The leak this catches is
                        specific: raise the level far enough and the region is the
                        whole mesh, whose boundary is nothing at all, so a rule that
                        takes the LAST level satisfying the neck bound is satisfied
                        again on the far side of the leak.
  WHOLE                 the region holds its entire cut loop and is one piece, with no
                        enclosed pockets. A stray concave vertex inside a stub makes
                        the flood go round it, and the region comes out freckled -- and
                        a freckle is a vertex that keeps the wall's zone, which is the
                        streak this exists to remove.
  ONE ZONE PER STUB     the requirement itself, checked on the labelling that ships.
  GREW                  a region that is its own cut loop and nothing else passes every
                        guard above -- it is bounded, whole, unfreckled and spills no
                        wall, because it is nothing -- while its rim is the cut opening
                        partway up the vessel instead of the ostium, so every rule that
                        reads the mouth reads it in the wrong place. KNOWN RED on
                        JK5_1's sixth stub, whose opening is crowded against a
                        neighbour's so the retry cannot grow it. That region was equally
                        wrong before this guard existed; the guard is what makes it
                        visible.
  NO ZONE INSIDE A ZONE a patch of zone 1 marooned in the middle of zone 0 is not a
                        stretch of aorta, it is a labelling artefact, and any area or
                        delta K taken on the zone averages across it. Taking a stub
                        back can strand a few vertices of the zone it was borrowed
                        from -- on KR52_1 an island of 8 -- so own_stubs repairs with
                        the stubs frozen, and this is the guard on that.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import branch_cut as BC                                      # noqa: E402
import region_station as RS                                  # noqa: E402
from zone_stub_bench import label_case, own_by_stub, wall_mask, streak  # noqa: E402
from zone_islands import islands                              # noqa: E402

LOCAL = [
    ("JK5_1", "JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
    ("KY160_1", "KY160_1_CentrelinePoints.csv",
     "Geometries & Meshes/mesh/KY160_1.mat"),
]
# The stranded-island defect needs an arch whose brachiocephalic mouth straddles the
# 0/1 boundary, and neither case in the repo has one -- so the guard on it could be
# written but not PROVEN. KR52_1 does, and is used when the cohort is mounted; when it
# is not, the prove run says so rather than passing quietly.
COHORT = r"Z:/Martino/branches_CAMP/GT"
COHORT_CASE = "KR52_1"
# KK16_1 is here for the GREW guard: its left subclavian flares just below where the
# nub cut crosses it, so the first rung of the ladder already carries a rim 1.68 times
# the cut loop and the length bound breaks before the flood has reported anything. The
# region came out as the seed loop plus one vertex, which is not a small stub but no
# stub at all -- and nothing else in the suite noticed, because a 60-vertex region is
# bounded, whole, unfreckled and spills nothing.
COHORT_CASES = (COHORT_CASE, "KK16_1")


def cases():
    # Geometry is gitignored on purpose, so a fresh clone carries some of these and
    # not others. Run on what is here rather than dying on a missing file.
    out = [(cid, ROOT / csv, ROOT / mat) for cid, csv, mat in LOCAL
           if (ROOT / csv).exists() and (ROOT / mat).exists()]
    try:
        from anat_bench import find_cases
        for cid, cl, mp in find_cases(COHORT):
            if cid in COHORT_CASES:
                out.append((cid, Path(cl), Path(mp)))
    except Exception:                                        # noqa: BLE001
        pass
    return out
MAX_STUB_FRAC = 0.20      # a stub is never a fifth of the aorta
WALL_TOL = 0.02           # a mouth vertex may sit inside the aortic envelope
# Vertices a region must hold beyond its own cut loop before the flood counts as having
# run. Held here rather than read from RS, so that --prove can disable the retry in the
# finder without also disabling the assertion that catches it.
MIN_GROWTH = 5

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print(f"    FAIL  {msg}")
    return cond


def mouth(e, sel):
    """The region's own boundary: its vertices with a neighbour outside it."""
    cut = sel[e[:, 0]] != sel[e[:, 1]]
    m = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
    return m[sel[m]]


def deep_spill(adj, e, ew, sel, wall):
    """Aortic wall inside the region, past the fillet at its mouth."""
    m = mouth(e, sel)
    if not len(m):
        return int((sel & wall).sum())
    lmouth = float(ew[(sel[e[:, 0]] != sel[e[:, 1]])].sum())
    g = dijkstra(adj, directed=False, indices=m, min_only=True)
    return int((sel & wall & (g > lmouth / (2.0 * np.pi))).sum())


def pockets(adj, sel, keep_out):
    """Islands of unselected surface wholly enclosed by the region."""
    out = ~sel
    if not out.any():
        return 0
    idx = np.flatnonzero(out)
    ncmp, lab = connected_components(adj[idx][:, idx], directed=False)
    if ncmp <= 1:
        return 0
    pos = np.searchsorted(idx, np.intersect1d(idx, np.asarray(keep_out, int)))
    return int(ncmp - len(set(lab[pos].tolist())))


def run():
    for cid, csv, mat in cases():
        s = CutSession(str(csv), str(mat))
        if not s.is_thoracic:
            # Since ad00bc1 the nub-in cut -- and the Ishimaru zones these stubs serve --
            # are thoracic-only, and nub_preset refuses an abdominal case by design. That
            # refusal used to crash this test on JK5_1, the repo's abdominal sample and
            # its first case, so the thoracic cases after it were never checked. JK5_1
            # was this test's one permanently red case (its crowded sixth stub fails the
            # GREW guard); it is out of scope now rather than a failure.
            print(f"{cid}: reads as abdominal; the nub-in cut is thoracic-only, skipped")
            continue
        s.nub_preset()
        s.recut()
        mesh = s.mesh
        route = BC.aorta_route(s.tree)
        adj = RS.surface_graph(mesh)
        reg = RS.stub_regions(mesh, route, adj=adj)
        wall = wall_mask(mesh, route)
        nv = len(mesh.vertices)
        loops = BC.boundary_loops(mesh)
        V = np.asarray(mesh.vertices, float)
        a, b = route.point(0.0), route.point(route.length)
        cen = [V[np.asarray(L_)].mean(axis=0) for L_ in loops]
        ia = int(np.argmin([np.linalg.norm(c - a) for c in cen]))
        ib = int(np.argmin([np.linalg.norm(c - b) for c in cen]))
        stubs = [np.asarray(L_, int) for k, L_ in enumerate(loops)
                 if k not in (ia, ib)]
        far = np.concatenate([np.asarray(loops[ia], int),
                              np.asarray(loops[ib], int)])
        e = np.asarray(mesh.edges_unique, int)
        ew = np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1)
        n = int(reg.max()) + 1
        print(f"\n{cid}: {nv:,} vertices, {len(stubs)} cut stub(s), "
              f"{n} region(s)")
        check(n > 0, f"{cid}: no stub region found at all")

        for j in range(n):
            sel = reg == j
            size = int(sel.sum())
            spill = deep_spill(adj, e, ew, sel, wall)
            pk = pockets(adj, sel, far)
            ncmp, _ = connected_components(
                adj[np.flatnonzero(sel)][:, np.flatnonzero(sel)], directed=False)
            print(f"  region {j}: {size:5d} verts  {100 * size / nv:4.1f}% of mesh  "
                  f"spill {spill:4d}  pieces {ncmp}  pockets {pk}")

            check(size >= 10, f"{cid} region {j}: {size} vertices -- nothing found")
            check(size <= MAX_STUB_FRAC * nv,
                  f"{cid} region {j}: {100 * size / nv:.0f}% of the mesh -- "
                  f"that is the aorta, not a stub")
            check(spill <= WALL_TOL * size,
                  f"{cid} region {j}: {spill} of {size} vertices are aortic wall "
                  f"deeper than the mouth, by the radius test")
            check(ncmp == 1, f"{cid} region {j}: delivered in {ncmp} pieces")
            check(pk == 0, f"{cid} region {j}: {pk} enclosed pocket(s) -- freckled")
            # WHOLE CUT LOOP. The stub is defined as the surface behind its own
            # opening, so an opening split between region and trunk is not a stub.
            hit = [int(sel[L_].sum()) / len(L_) for L_ in stubs]
            check(max(hit) > 0.99 or max(hit) == 0.0,
                  f"{cid} region {j}: best cut loop only "
                  f"{100 * max(hit):.0f}% inside the region")
            # GREW PAST ITS OWN CUT LOOP. Every guard above is satisfied by a region
            # that is simply the loop the flood started from -- it is bounded, whole,
            # unfreckled and spills nothing, because it is nothing. But its rim is then
            # the cut opening partway up the vessel rather than the ostium, so the
            # mouth every later rule reads is in the wrong place: on KK16_1 that put
            # the subclavian's lip 11 mm downstream and the zone 2/3 boundary with it.
            own = int(np.argmax(hit))
            if hit[own] > 0.99:
                grew = size - len(stubs[own])
                check(grew >= MIN_GROWTH,
                      f"{cid} region {j}: only {grew} vertices beyond its own "
                      f"{len(stubs[own])}-vertex cut loop -- the flood never ran")

        # ONE ZONE PER STUB, on the labelling that ships.
        mesh2, route2, zones, lab, isb, adj2, _flat = label_case(s)
        reg2 = RS.stub_regions(mesh2, route2, adj=adj2)
        new = own_by_stub(mesh2, lab, reg2, isb, adj2)
        chg = new != lab
        wall2 = wall_mask(mesh2, route2)
        ns, before, _ = streak(lab, reg2)
        _n, after, _ = streak(new, reg2)
        print(f"  zones: {chg.sum()} label(s) changed, {int((chg & wall2).sum())} "
              f"of them aortic wall; stubs split {before} -> {after} of {ns}")
        check(after == 0, f"{cid}: {after} stub(s) still carry two zones")
        check(int((chg & wall2).sum()) == 0,
              f"{cid}: giving stubs to their ostium moved "
              f"{int((chg & wall2).sum())} aortic wall vertices")
        isl = islands(mesh2, new)
        for z, n, around in isl:
            print(f"    island: zone {z}, {n} vertices, inside zone {around}")
        check(not isl, f"{cid}: {len(isl)} zone island(s) inside another zone")


def prove():
    """Put each defect back and check the suite notices."""
    good_regions = RS.stub_regions
    good_fill = RS._fill_pockets
    good_own = RS.own_stubs

    def last_success(mesh, route, alpha=RS.NECK_ALPHA, smooth=2, adj=None):
        """Take the LAST level that satisfies the neck bound, not the first failure.

        Past the leak the region is the whole mesh and its boundary is nothing, so the
        bound comes good again and the rule sails past the mouth.
        """
        return _regions_variant(mesh, route, alpha, smooth, adj, stop_early=False)

    def no_fill(adj, sel, keep_out):
        return sel

    def no_repair(mesh, label, regions, is_branch, adj=None, region_zone=None):
        """Give the stub its ostium's zone and stop -- no connectivity repair.

        This is what strands the vertices the stub was taken back from.
        """
        lab = np.asarray(label, int).copy()
        reg, isb = np.asarray(regions, int), np.asarray(is_branch, bool)
        for j in range(int(reg.max()) + 1):
            sel = reg == j
            if not sel.any():
                continue
            cand = lab[sel & isb]
            cand = cand[cand >= 0]
            if not len(cand):
                cand = lab[sel][lab[sel] >= 0]
            if len(cand):
                lab[sel] = int(np.bincount(cand).argmax())
        return lab

    good_growth = RS.MIN_STUB_GROWTH
    # WHAT THE SUITE ALREADY SAYS WITH NOTHING BROKEN. A defect used to count as
    # caught if any assertion fired at all, and that stopped being a proof the moment
    # one case went permanently red: JK5_1's crowded stub fails the GREW guard on every
    # run, so "did anything fail" is true whatever the defect does. A defect is caught
    # only if it fires an assertion that the healthy code does NOT.
    del fails[:]
    print(f"\n{'=' * 78}\nPROVE: baseline, nothing broken")
    run()
    baseline = set(fails)
    del fails[:]
    for f in sorted(baseline):
        print(f"  known red: {f}")

    for label, attr, repl in (("stop at last success", "stub_regions", last_success),
                              ("no pocket fill", "_fill_pockets", no_fill),
                              ("no connectivity repair", "own_stubs", no_repair),
                              # A negative threshold makes the size test unreachable,
                              # which is the code before the retry existed.
                              ("no retry for a failed flood", "MIN_STUB_GROWTH", -1)):
        setattr(RS, attr, repl)
        del fails[:]
        print(f"\n{'=' * 78}\nPROVE: {label} defect reinstated")
        try:
            run()
        finally:
            setattr(RS, attr, {"stub_regions": good_regions,
                               "_fill_pockets": good_fill,
                               "own_stubs": good_own,
                               "MIN_STUB_GROWTH": good_growth}[attr])
        fresh = [x for x in fails if x not in baseline]
        if fresh:
            print(f"  caught it: {len(fresh)} assertion(s) fired that do not fire on "
                  f"healthy code")
            for f in fresh[:4]:
                print(f"    {f}")
        elif label == "no connectivity repair" and len(cases()) < 3:
            print(f"  NOT EXERCISED -- {COHORT_CASE} is not mounted, so no case here "
                  f"strands a vertex. The guard is untested, not passed.")
        else:
            print(f"  NOT CAUGHT -- the test does not guard {label}")
            fails.append(f"prove: {label} went undetected")
        keep = [x for x in fails if x.startswith("prove:")]
        del fails[:]
        fails.extend(keep)


def _regions_variant(mesh, route, alpha, smooth, adj, stop_early):
    """stub_regions with the stopping rule made switchable, for --prove only."""
    from branch_cut import boundary_loops
    V = np.asarray(mesh.vertices, float)
    out = -np.ones(len(V), int)
    loops = boundary_loops(mesh)
    if len(loops) < 3:
        return out
    a, b = route.point(0.0), route.point(route.length)
    cen = [V[np.asarray(L_)].mean(axis=0) for L_ in loops]
    ia = int(np.argmin([np.linalg.norm(c - a) for c in cen]))
    ib = int(np.argmin([np.linalg.norm(c - b) for c in cen]))
    conc = RS.surface_concavity(mesh, smooth)
    e = np.asarray(mesh.edges_unique, int)
    ew = np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1)
    adjlist = [[] for _ in range(len(V))]
    for i, j in e:
        adjlist[i].append(j)
        adjlist[j].append(i)
    adj = adj if adj is not None else RS.surface_graph(mesh)
    far = np.concatenate([np.asarray(loops[ia], int), np.asarray(loops[ib], int)])
    for j, k in enumerate([x for x in range(len(loops)) if x not in (ia, ib)]):
        seed = np.asarray(loops[k], int)
        key = RS._bottleneck(adjlist, conc, seed)
        fin = np.isfinite(key)
        if fin.sum() < len(seed) * 2:
            continue
        bs = np.sort(key[fin])
        onloop = np.zeros(len(V), bool)
        onloop[seed] = True
        lopen = float(np.linalg.norm(V[seed] - V[np.roll(seed, 1)], axis=1).sum())
        best = None
        for t in np.linspace(bs[0], bs[int(RS.NECK_WINDOW * len(bs))],
                             RS.NECK_LEVELS):
            sel = onloop | (fin & (key <= t))
            cut = sel[e[:, 0]] != sel[e[:, 1]]
            over = float(ew[cut].sum()) > alpha * max(lopen, 1e-9)
            if over and stop_early:
                break
            if not over:
                best = sel
        if best is None:
            continue
        out[RS._fill_pockets(adj, best, far)] = j
    return out


if __name__ == "__main__":
    if "--prove" in sys.argv:
        prove()
    else:
        run()
    print("\n" + "=" * 78)
    if fails:
        print(f"{len(fails)} FAILURE(S)")
        for f in fails:
            print(f"  {f}")
        sys.exit(1)
    print("all checks passed")
