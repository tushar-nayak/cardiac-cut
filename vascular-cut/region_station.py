"""Isolate an anatomical zone by CENTERLINE STATION, not with a pair of planes.

A branch cut can insist on a flat rim: it is placed where the branch is a clean tube,
and if it will not come out flat the answer is to move it a few millimetres and try
again. A ZONE cannot. Its boundary is fixed by anatomy -- the lowest renal, the left
subclavian, 20 mm past it -- and on a curved or aneurysmal aorta there is often no
planar cross-section at that station at all. The tool's answer until now was to walk
the boundary along a ladder of alternative stations until one cut flat, which buys a
planar rim by moving the zone off the anatomy that defines it. That is the wrong trade
for a measurement: a zone whose end is 10 mm from where it should be is not the zone.

So a zone boundary here is not a plane. Every vertex is given a station along the
aorta and the zone is the set of vertices whose station lies in [s0, s1]. The
boundary is the level set of that field: it follows the surface, it exists on every
geometry, and it sits exactly where the anatomy says.

WHAT THIS COSTS. The rim is a closed curve on the surface, not a disc, so the result
is not a CFD inlet and must not be described as one. What it is fit for is measuring
-- area, curvature, delta K -- where the rim's shape does not enter the answer.

WHAT A BRANCH DOES. Its station is CONSTANT along its whole length: the station of
its ostium on the aorta. So a branch is kept or dropped whole, by where it leaves the
aorta, and no isocontour ever runs up a renal artery. That is the same membership rule
region_plan already applies (see its "Membership is decided by each branch's OSTIUM
STATION" note), expressed as a field instead of as a pair of planes.

A ZONE KEEPS ITS BRANCHES, AND THE ZONES DO NOT OVERLAP. Those two requirements
together fix the rule exactly, and neither alone does.

A zone exists to say what the geometry does at one place on the aorta, and at the
places these zones are drawn the interesting geometry IS the junction: Ishimaru zones 1
and 2 are defined by the great vessels bounding them, so a zone that drops its branches
is a bare tube with the branch flare -- the thing that makes that stretch different from
any other -- cut away. So branches are kept.

But the Ishimaru zones TILE the aorta: zone 1 ends exactly where zone 2 begins. Keep a
branch at both of its zone's bounds and every landmark vessel lands in two zones at
once -- the left carotid is zone 1's distal bound and zone 2's proximal bound, so both
claim it, and a zone that should show one vessel shows several.

So membership is the HALF-OPEN interval s0 < ostium <= s1 -- out at the proximal
bound, IN at the distal one. That is not an arbitrary choice between two half-open
conventions: it is what the zone names mean. A zone contains the vessel at its far
end, and covering that zone means covering that vessel's origin:

  zone 0   ascending, to the brachiocephalic    THE BRACHIOCEPHALIC, and with it the
                                                daughters that share its ostium
  zone 1   brachiocephalic to left carotid      the left carotid
  zone 2   left carotid to left subclavian      the left subclavian
  zone 3   left subclavian, +20 mm              none
  zone 4   beyond that                          none

Anything else a zone contains is a vessel that genuinely arises inside it, which is
worth seeing rather than hiding.

The mirror convention -- in at the proximal bound, out at the distal -- also tiles and
was built first. It puts the brachiocephalic in zone 1 rather than zone 0, which is
not what zone 0 means.
"""
import re

import numpy as np
import trimesh
from scipy.spatial import cKDTree

import branch_cut as BC
from branch_cut import (
    _junction_graph, _station_on_route, aorta_route, keep_largest_component,
)

# Spacing the centerlines are resampled to before the nearest-point lookup. The
# station field is only as sharp as this: at 2 mm a vertex can take the station of a
# point 1 mm away along the curve, which puts the same jitter on the rim. 0.4 mm is
# well under the mesh's own edge length on every case here and costs nothing that
# matters -- an aorta's centerlines come to ~1500 samples.
SAMPLE_MM = 0.4

# How close to a bound an ostium has to be to count as sitting ON it. A landmark's
# station and the zone bound derived from it should agree exactly, but they are
# computed by different code -- aorta_landmark_stations for one, the junction table
# for the other -- so the comparison is given a millimetre and a half of room rather
# than being decided by the last bit of a float.
BOUND_TOL_MM = 1.5


def _segment_station_samples(tree, route, jg, mesh=None):
    """Every centerline point, with the route station it stands for.

    On-route segments carry their own arc length. Everything else carries the station
    of the junction where its subtree meets the route -- one value for the whole
    branch, which is what makes a branch indivisible by a zone boundary.
    """
    jstation = {}
    for i, name in enumerate(route.chain):
        ns, ne = jg[name]
        seg = tree.segments[name]
        a, b = (ns, ne) if route.forward[i] else (ne, ns)
        if a is not None:
            jstation.setdefault(a, route.starts[i])
        if b is not None:
            jstation.setdefault(b, route.starts[i] + seg.length)

    pts, tans, sta, owner = [], [], [], []
    for name, seg in tree.segments.items():
        n = max(int(seg.length / SAMPLE_MM) + 1, 2)
        local = np.linspace(0.0, seg.length, n)
        p = np.column_stack([np.interp(local, seg.s, seg.pts[:, k]) for k in range(3)])
        tg = np.column_stack([np.interp(local, seg.s, seg.tangents[:, k])
                              for k in range(3)])
        tg /= np.maximum(np.linalg.norm(tg, axis=1, keepdims=True), 1e-12)
        if name in route.chain:
            s = np.array([route.station_of(name, float(x)) for x in local])
        else:
            ns, ne = jg.get(name, (None, None))
            st = None
            for j in (ns, ne):
                st = _station_on_route(jg, jstation, j)
                if st is not None:
                    break
            if st is None:
                # Not reachable from the route -- a stray component. Give it a station
                # outside every zone rather than a wrong one inside one.
                st = -1e9
            s = np.full(n, float(st))
        pts.append(p)
        tans.append(tg)
        sta.append(s)
        owner += [name] * n
    P, T, S, O = (np.vstack(pts), np.vstack(tans), np.concatenate(sta),
                  np.array(owner))
    # LOCAL RADIUS PER SAMPLE, measured off the surface rather than assumed. The
    # distance from a centerline point to the nearest wall is that vessel's radius
    # there, near enough, and it is what turns a raw distance into a comparable one.
    if mesh is None:
        R = np.ones(len(P))
    else:
        d, _ = cKDTree(np.asarray(mesh.vertices, float)).query(P)
        R = np.maximum(d, 0.5)
    return P, T, S, O, R


# How far from perpendicular the line from a centerline point to a vertex may be
# before that point is disqualified from owning it. A vertex on a vessel wall sits
# broadly radial to its own centerline; one that a folded-back neighbour happens to
# be nearer to sits nearly ALONG that neighbour's tangent.
PERP_COS_MAX = 0.55                       # |cos| 0.55 is about 33 degrees off radial
NEIGHBOURS = 48


def station_field(mesh, tree, route=None, jg=None):
    """Route station for every vertex of *mesh*.

    NEAREST CENTERLINE POINT IS NOT ENOUGH, and an aortic arch is where that shows.
    At the top of the arch the ascending and descending limbs run within a few
    millimetres of each other, so a vertex on one limb is often physically closest to
    the OTHER limb's centerline and takes its station -- 30 to 70 mm wrong. Zones then
    reach past their own bounds: measured on KY160_1, Ishimaru zone 3 came out
    spanning 129.7-158.9 mm for a span asked as 131.3-151.3.

    Two things a candidate has to be, not one:

    RADIAL     the line from the centerline point out to the vertex is roughly
               perpendicular to that centerline. A vertex looking across a fold at the
               far limb sees it nearly end-on.
    INSIDE     the centerline point lies on the INWARD side of the wall at that vertex,
               (p - v) . n < 0 for the outward vertex normal n. This is the one that
               actually settles the hairpin. A centerline belongs to the lumen its wall
               encloses, so the point that owns a vertex is always behind that vertex's
               surface; the other limb of a fold, however close, is in a different
               lumen and is therefore in FRONT of this wall. Radial alone does not
               catch it -- where the two limbs run parallel, the far limb is seen
               broadside and looks perfectly radial.

    The proper test is whether the segment from p to v stays inside the lumen, which is
    a ray cast; with no embree in this environment that is far too slow for 30k
    vertices times 48 candidates. The normal test is one dot product, and on a closed
    wound surface it answers the same question everywhere except within a hair of the
    surface itself.

    Among the candidates that survive, the nearest wins; if none do, the nearest
    overall is used rather than leaving the vertex unassigned.

    Returns (station, owner_segment_name_per_vertex, on_route_mask).
    """
    jg = jg or _junction_graph(tree)
    route = route or aorta_route(tree, jg)
    if route is None:
        raise ValueError("no aorta route: the family call has to succeed first")
    pts, tans, sta, owner, rad = _segment_station_samples(tree, route, jg, mesh)
    on_route = np.isin(owner, list(route.chain))

    V = np.asarray(mesh.vertices, float)
    k = int(min(NEIGHBOURS, len(pts)))
    d, idx = cKDTree(pts).query(V, k=k)
    if k == 1:
        d, idx = d[:, None], idx[:, None]
    # NORMALISED BY THE VESSEL'S OWN RADIUS. Raw distance says a wall vertex belongs
    # to whatever centerline runs nearest, which on the arch hands great swathes of
    # aortic wall to a great vessel passing a few millimetres away. Divided by the
    # local radius the question becomes "whose lumen is this the wall of", and a
    # vertex sits at about ONE radius from its own centerline and several from any
    # other. Candidates are re-ordered on that ratio.
    ratio = d / rad[idx]
    ordr = np.argsort(ratio, axis=1)
    rows0 = np.arange(len(V))[:, None]
    idx = idx[rows0, ordr]
    d = d[rows0, ordr]

    rel = V[:, None, :] - pts[idx]                      # vertex minus candidate
    n = np.maximum(np.linalg.norm(rel, axis=2), 1e-12)
    cos = np.abs(np.einsum("ijk,ijk->ij", rel / n[..., None], tans[idx]))
    # (p - v) . n < 0  <=>  (v - p) . n > 0, with n the OUTWARD vertex normal.
    nrm = np.asarray(mesh.vertex_normals, float)
    inward = np.einsum("ijk,ik->ij", rel, nrm) > 0
    ok = (cos <= PERP_COS_MAX) & inward
    # d comes back sorted, so the first acceptable column IS the nearest acceptable.
    first = np.argmax(ok, axis=1)
    none_ok = ~ok.any(axis=1)
    first[none_ok] = 0
    row = np.arange(len(V))
    pick = idx[row, first]

    # PROJECT ONTO THE TANGENT. Taking the chosen sample's arc length wholesale
    # quantises the field to whichever sample won, and neighbouring vertices often
    # pick samples a few millimetres apart -- so the level set, which is the zone
    # boundary, inherits that as a ragged edge. Zone 3 sits in the descending aorta
    # where a cross-section IS very nearly a plane, and its boundaries came out
    # jagged for no geometric reason at all.
    #
    # s = s(p) + (v - p) . t(p) is the first-order correction: it measures how far
    # along the centerline the vertex actually sits rather than which sample was
    # nearest. For a true perpendicular foot the correction is zero, so this changes
    # nothing where the pick was already exact, and removes the staircase everywhere
    # else. The level set becomes the normal plane at that station, which is what a
    # cross-section means.
    rel_pick = V - pts[pick]
    along = np.einsum("ij,ij->i", rel_pick, tans[pick])
    station = sta[pick] + along
    # A BRANCH KEEPS ITS OSTIUM STATION, uncorrected. Its whole length stands for one
    # station by construction; projecting would spread a renal artery across 30 mm of
    # aorta and let a zone boundary cut through the middle of it.
    on_r = on_route[pick]
    station = np.where(on_r, station, sta[pick])
    return station, owner[pick], on_r


def branch_inside(ostium, s0, s1, route_end=None):
    """Does a branch with this ostium station belong to the zone [s0, s1]?

    HALF-OPEN, OUT AT THE PROXIMAL BOUND AND IN AT THE DISTAL ONE. A zone contains the
    vessel at its far end, which is what the zone names mean: zone 0 includes the
    brachiocephalic, zone 1 the left carotid, zone 2 the left subclavian, and covering
    a zone means covering that vessel's origin. Half-open either way makes the zones
    tile; only this direction puts each vessel in the zone it is named for.

    ...except at the route's own end, where what leaves is the aorta's continuation
    rather than a branch off the last zone. Without that "Infrarenal aorta" keeps the
    iliacs at the bifurcation and becomes "Infrarenal aorta with iliacs", which is the
    one distinction those two presets exist to draw.

    Kept as its own function because the cohort checks ask the same question, and a
    checker with its own copy of the rule can only ever agree with itself.
    """
    ost = np.asarray(ostium, float)
    lo = float(s0) if s0 is not None else -np.inf
    hi = float(s1) if s1 is not None else np.inf
    inside = (ost > lo + BOUND_TOL_MM) & (ost <= hi + BOUND_TOL_MM)
    if route_end is not None and s1 is not None and s1 >= route_end - BOUND_TOL_MM:
        inside &= ost < route_end - BOUND_TOL_MM
    return inside


def zone_field(station, on_route, s0, s1, route_end=None):
    """The station field rewritten so ONE level-set cut does the whole zone.

    The aorta keeps its true arc length, because the boundary has to land at the
    station asked for and nowhere else. A branch is not cut at all -- it is kept or
    dropped whole -- so every vertex of an INCLUDED branch is moved to a value safely
    inside the span and every vertex of an excluded one to a value just outside it.
    No isocontour can then pass through a branch, and the clip needs no special case.

    Just outside, not far outside: the cut point on an edge is interpolated, so a
    sentinel of -1e9 against a neighbour at 100 lands the new vertex on top of that
    neighbour and makes a zero-area triangle. One millimetre past the bound puts it
    where it belongs, at the branch mouth.
    """
    lo = float(s0) if s0 is not None else -np.inf
    hi = float(s1) if s1 is not None else np.inf
    f = np.asarray(station, float).copy()
    br = ~np.asarray(on_route, bool)
    if not br.any():
        return f
    ost = f[br]
    inside = branch_inside(ost, s0, s1, route_end)
    mid = (max(lo, ost.min() - 1.0) + min(hi, ost.max() + 1.0)) / 2.0
    if not np.isfinite(mid):
        mid = (lo + 1.0) if np.isfinite(lo) else (hi - 1.0)
    out = np.where(ost <= lo, lo - 1.0, hi + 1.0)
    f[br] = np.where(inside, mid, out)
    return f


def clip_by_field(mesh, f, level, keep_above):
    """Cut a mesh along a level set of a per-vertex scalar, exactly.

    Faces wholly on one side are kept or dropped whole; a face the level crosses is
    split on its edges by linear interpolation, so the rim lands ON the level set
    rather than following the nearest run of triangle edges. That is the difference
    between a boundary at the station asked for and a boundary up to one edge length
    away from it, which on a 100k-face aorta is about 1 mm of anatomy.

    Returns (mesh, field) -- the field is carried through so a second cut can be made
    against the same coordinate without recomputing it.
    """
    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    f = np.asarray(f, float)
    inside = (f >= level) if keep_above else (f <= level)

    cnt = inside[F].sum(axis=1)
    keep_all = np.nonzero(cnt == 3)[0]
    cross = np.nonzero((cnt == 1) | (cnt == 2))[0]

    newV = [V]
    newf = [f]
    newF = [F[keep_all]]
    # One interpolated vertex per crossed EDGE, shared by both faces on it, or the
    # two sides of the rim come out as separate points and the surface is torn along
    # the boundary it was supposed to be cut at.
    edge_pt = {}

    def cut_point(i, j):
        key = (i, j) if i < j else (j, i)
        hit = edge_pt.get(key)
        if hit is not None:
            return hit
        fi, fj = f[key[0]], f[key[1]]
        t = 0.5 if abs(fj - fi) < 1e-12 else (level - fi) / (fj - fi)
        t = float(np.clip(t, 0.0, 1.0))
        p = V[key[0]] * (1.0 - t) + V[key[1]] * t
        n = len(V) + len(edge_pt)
        edge_pt[key] = n
        newV.append(p.reshape(1, 3))
        newf.append(np.array([level]))
        return n

    # CYCLIC ORDER IS THE WHOLE THING. Filtering a face into "inside" and "outside"
    # lists throws away which vertex followed which, and rebuilding from those lists
    # reverses roughly half the faces -- the mesh comes back with every zone
    # winding-inconsistent, normals pointing both ways, and no way to tell inside from
    # out. So rotate the face until the split sits in a known place and emit against
    # the original cycle.
    tris = []
    for fi in cross:
        tri = F[fi]
        ins = inside[tri]
        if ins.sum() == 1:
            k = int(np.argmax(ins))                 # rotate the lone inside one to 0
            v0, v1, v2 = tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]
            tris.append((v0, cut_point(v0, v1), cut_point(v2, v0)))
        else:
            k = int(np.argmin(ins))                 # rotate the lone OUTSIDE one to 2
            v0, v1, v2 = tri[(k + 1) % 3], tri[(k + 2) % 3], tri[k]
            m, n = cut_point(v1, v2), cut_point(v2, v0)
            tris.append((v0, v1, m))
            tris.append((v0, m, n))
    if tris:
        newF.append(np.asarray(tris, int))

    Vout = np.vstack(newV)
    fout = np.concatenate(newf)
    Fout = np.vstack(newF) if newF else np.zeros((0, 3), int)

    # A level that falls exactly on a vertex interpolates to that vertex, which makes
    # a triangle with two coincident corners. It is harmless to look at and poison to
    # measure with -- a zero-area face has no normal, and area-weighted statistics
    # (delta K among them) divide by what it contributes.
    tv = Vout[Fout]
    area2 = np.linalg.norm(np.cross(tv[:, 1] - tv[:, 0], tv[:, 2] - tv[:, 0]), axis=1)
    Fout = Fout[area2 > 1e-9]

    used = np.zeros(len(Vout), bool)
    used[Fout.ravel()] = True
    remap = -np.ones(len(Vout), int)
    remap[used] = np.arange(used.sum())
    out = trimesh.Trimesh(Vout[used], remap[Fout], process=False)
    return out, fout[used]


def _component_sizes(mesh):
    """Face counts of every surface-connected piece, largest first."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    fa = np.asarray(mesh.face_adjacency)
    n = len(mesh.faces)
    if n == 0:
        return []
    if len(fa) == 0:
        return [1] * n
    g = sp.coo_matrix((np.ones(len(fa)), (fa[:, 0], fa[:, 1])), shape=(n, n))
    k, lab = connected_components(g, directed=False)
    return sorted(np.bincount(lab, minlength=k).tolist(), reverse=True)


def enforce_connected_labels(mesh, label, max_iter=12, frozen=None):
    """Make every zone one connected piece by REASSIGNING strays, not deleting them.

    keep_largest_component answers the same requirement by throwing the small pieces
    away, which is right when one zone is being isolated on its own -- there is nowhere
    else for them to go -- and wrong when the whole aorta is being labelled, because
    every one of those faces belongs to SOME zone and dropping it leaves a hole in a
    partition that is supposed to cover the vessel.

    So each label's components are found on the surface graph, the largest is kept, and
    every other one is given to whichever label its own boundary neighbours most. That
    converges quickly -- an island is small and is surrounded by one neighbour almost
    by definition -- and it is capped anyway, because a pathological case that keeps
    handing a patch back and forth should stop rather than spin.

    Returns (label, moved) with moved the number of vertices reassigned.
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    lab = np.asarray(label, int).copy()
    e = np.asarray(mesh.edges_unique, int)
    nv = len(lab)
    moved = 0
    # A branch's zone comes from its ostium and is not a matter of connectivity: a
    # great vessel whose mouth sits on a zone boundary is its own island by
    # construction, and reassigning it to whichever neighbour it touches is exactly
    # the "zone 1 grabbed zone 0's branch" fault.
    frozen = (np.zeros(nv, bool) if frozen is None else np.asarray(frozen, bool))
    for _ in range(max_iter):
        same = lab[e[:, 0]] == lab[e[:, 1]]
        ee = e[same]
        g = sp.coo_matrix((np.ones(len(ee)), (ee[:, 0], ee[:, 1])), shape=(nv, nv))
        k, comp = connected_components(g, directed=False)
        # Largest component per label; everything else is an island.
        order = np.lexsort((np.arange(nv), comp))
        sizes = np.bincount(comp, minlength=k)
        best = {}
        for c in range(k):
            v = int(np.argmax(comp == c)) if sizes[c] else -1
            if v < 0:
                continue
            L = int(lab[v])
            if L < 0:
                continue
            if L not in best or sizes[c] > sizes[best[L]]:
                best[L] = c
        island = np.zeros(nv, bool)
        for c in range(k):
            v = int(np.argmax(comp == c)) if sizes[c] else -1
            if v < 0:
                continue
            L = int(lab[v])
            if L >= 0 and best.get(L) != c:
                island |= (comp == c)
        if not island.any():
            break
        # Neighbour vote, counting only neighbours OUTSIDE the island being replaced.
        a, b = e[:, 0], e[:, 1]
        cand = np.full(nv, -2, int)
        for u, w in ((a, b), (b, a)):
            take = island[u] & ~island[w] & (lab[w] >= 0)
            cand[u[take]] = lab[w[take]]
        fix = island & (cand != -2) & ~frozen
        if not fix.any():
            break
        lab[fix] = cand[fix]
        moved += int(fix.sum())

    # UNLABELLED WALL GOES TO THE NEIGHBOUR IT TOUCHES. A vertex can end up with no
    # zone -- the harmonic coordinate at the very inlet and outlet rings sits a hair
    # outside the first and last zone's range -- and leaving it grey puts a hole in a
    # partition that is meant to cover the vessel. It is wall, it is connected to
    # exactly one zone, and that is where it belongs. Flood outward until nothing
    # unlabelled has a labelled neighbour left.
    a, b = e[:, 0], e[:, 1]
    for _ in range(200):
        gap = lab < 0
        if not gap.any():
            break
        cand = np.full(nv, -2, int)
        for u, w in ((a, b), (b, a)):
            take = gap[u] & (lab[w] >= 0)
            cand[u[take]] = lab[w[take]]
        fix = gap & (cand != -2) & ~frozen
        if not fix.any():
            break
        lab[fix] = cand[fix]
        moved += int(fix.sum())
    return lab, moved


def region_by_station(mesh, tree, s0, s1, route=None, jg=None, largest=True):
    """The zone between two stations, cut along the surface rather than by planes.

    s1=None keeps everything distal to s0, matching region_plan's open distal end.

    Returns (mesh, info) where info carries what a caller has to be able to report:
    the stations actually used, the rim loops produced, and the fact that the
    boundaries are NOT planar -- which the sidecar has to record, or a reader will
    assume they are.
    """
    jg = jg or _junction_graph(tree)
    route = route or aorta_route(tree, jg)
    station, _owner, on_route = station_field(mesh, tree, route, jg)
    f = zone_field(station, on_route, s0, s1,
                   route_end=float(route.length))

    out = mesh
    if s0 is not None:
        out, f = clip_by_field(out, f, float(s0), keep_above=True)
    if s1 is not None:
        out, f = clip_by_field(out, f, float(s1), keep_above=False)
    # A ZONE IS ONE PIECE. That is a standard, not a preference: a zone reported as
    # a single area that is really two patches on opposite walls describes neither of
    # them, and every per-zone number downstream -- area, delta K -- silently averages
    # across the gap. So the components are COUNTED before the largest is kept, and
    # what was discarded is reported, because a zone that arrives in pieces is
    # evidence about the station field and not a tidying job.
    parts = _component_sizes(out)
    if largest:
        before = len(out.faces)
        out = keep_largest_component(out)
        dropped = before - len(out.faces)
    else:
        dropped = 0

    from branch_cut import boundary_loops
    loops = boundary_loops(out)
    info = {
        "kind": "station",
        "n_components": len(parts),
        "island_faces": int(sum(sorted(parts, reverse=True)[1:])),
        "largest_island": int(sorted(parts, reverse=True)[1]) if len(parts) > 1 else 0,
        "s0": None if s0 is None else float(s0),
        "s1": None if s1 is None else float(s1),
        "planar": False,
        "n_loops": len(loops),
        "faces": int(len(out.faces)),
        "area_mm2": float(out.area),
        "euler": int(out.euler_number),
        "faces_dropped_as_islands": int(dropped),
    }
    return out, info


# ---------------------------------------------------------------------------
# A coordinate that does not project.
#
# Everything above asks, for each wall vertex, which centerline point owns it. That
# question has no answer where the vessel bends tighter than it is wide: nearest-point
# projection onto a curve is single-valued only inside the REACH, 1/curvature, and on
# the arch k*d runs to 1.4 median and 13 peak around Ishimaru zone 3. Past k*d = 1 the
# normal planes have crossed before they reach the wall, so a wall point has several
# perpendicular feet and its station is genuinely multivalued. That is why the zone 2
# to zone 3 boundary is ragged, and why filtering candidates, projecting onto
# tangents, repairing components and snapping to planes all failed to fix it -- none
# of them changed the question being asked.
#
# So ask a different one. Solve Laplace's equation on the wall with the inlet ring
# held at 0 and the outlet ring at 1. The solution is smooth by construction, monotone
# from one end to the other, and its level sets are closed curves around the vessel --
# cross-sections in every sense that matters here, with no reach condition anywhere.
# Branch mouths are left as natural (insulated) boundaries, which is why a branch stub
# comes out at very nearly the value of the mouth it leaves from, exactly as a branch
# should sit at one station.

def _cotan_laplacian(V, F):
    """Cotangent-weighted Laplacian. Uniform weights blur where triangles are thin."""
    import scipy.sparse as sp
    i0, i1, i2 = F[:, 0], F[:, 1], F[:, 2]
    e0, e1, e2 = V[i2] - V[i1], V[i0] - V[i2], V[i1] - V[i0]

    def cot(a, b):
        cross = np.linalg.norm(np.cross(a, b), axis=1)
        return np.einsum("ij,ij->i", a, b) / np.maximum(cross, 1e-12)

    c0, c1, c2 = cot(-e1, e2), cot(-e2, e0), cot(-e0, e1)
    I = np.concatenate([i1, i2, i2, i0, i0, i1])
    J = np.concatenate([i2, i1, i0, i2, i1, i0])
    W = np.concatenate([c0, c0, c1, c1, c2, c2]) * 0.5
    n = len(V)
    A = sp.coo_matrix((W, (I, J)), shape=(n, n)).tocsr()
    A = (A + A.T) * 0.5                        # symmetric, and kills sign noise
    d = np.asarray(A.sum(axis=1)).ravel()
    return sp.diags(d) - A


def harmonic_coordinate(mesh, route):
    """Smooth 0..1 coordinate along the vessel, from the inlet ring to the outlet.

    Returns None when the two aortic ends cannot be identified, so a caller can fall
    back rather than get a coordinate built on the wrong holes.
    """
    import scipy.sparse.linalg as spl
    from branch_cut import boundary_loops

    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    loops = boundary_loops(mesh)
    if len(loops) < 2:
        return None
    # The aorta's own two ends are the loops nearest the route's endpoints. Branch
    # mouths are nearer to their own ostia and are left alone.
    a, b = route.point(0.0), route.point(route.length)
    cen = [V[np.asarray(L_)].mean(axis=0) for L_ in loops]
    ia = int(np.argmin([np.linalg.norm(c - a) for c in cen]))
    ib = int(np.argmin([np.linalg.norm(c - b) for c in cen]))
    if ia == ib:
        return None

    L = _cotan_laplacian(V, F)
    fixed = np.zeros(len(V), bool)
    val = np.zeros(len(V))
    fixed[np.asarray(loops[ia])] = True
    fixed[np.asarray(loops[ib])] = True
    val[np.asarray(loops[ib])] = 1.0
    free = ~fixed
    if not free.any():
        return None
    try:
        rhs = -L[free][:, fixed] @ val[fixed]
        u = np.zeros(len(V))
        u[fixed] = val[fixed]
        u[free] = spl.spsolve(L[free][:, free].tocsc(), rhs)
    except Exception:                                        # noqa: BLE001
        return None
    if not np.all(np.isfinite(u)):
        return None
    return np.clip(u, 0.0, 1.0)


def calibrate_u_to_station(mesh, route, u, step_mm=1.0):
    """Monotone map from the harmonic coordinate to arc length along the trunk.

    Sampled with a BALL around each centerline point, never a plane: a ball is
    well-defined wherever the centerline goes, and the whole reason this coordinate
    exists is that planes are not.
    """
    from scipy.spatial import cKDTree
    V = np.asarray(mesh.vertices, float)
    tree = cKDTree(V)
    ss, uu = [], []
    for s in np.arange(0.0, float(route.length) + 1e-9, step_mm):
        p = route.point(float(s))
        d, _i = tree.query(p)
        near = tree.query_ball_point(p, r=max(1.6 * float(d), 3.0))
        if len(near) < 6:
            continue
        ss.append(s)
        uu.append(float(np.median(u[np.asarray(near)])))
    if len(ss) < 4:
        return None, None
    ss, uu = np.asarray(ss), np.asarray(uu)
    # Force monotone, so inverting it is unambiguous.
    uu = np.maximum.accumulate(uu)
    keep = np.concatenate(([True], np.diff(uu) > 1e-9))
    return ss[keep], uu[keep]


def harmonic_gradient(mesh, u):
    """Per-vertex |grad u|, area-weighted from the faces.

    A branch stub carries no flux -- its far end is an insulated hole -- so u is very
    nearly constant along it and |grad u| collapses. The trunk carries all of it. That
    is a way to tell wall from branch that owes nothing to nearest-point projection,
    which is the point: the projection's trunk/branch call fails in the same place its
    station does, and using it to decide which vertices get the branch rule put whole
    patches of arch wall into a great vessel's zone.
    """
    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    u0, u1, u2 = u[F[:, 0]], u[F[:, 1]], u[F[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    a2 = np.maximum(np.linalg.norm(n, axis=1), 1e-12)
    nh = n / a2[:, None]
    g = (np.cross(nh, (p2 - p1)) * u0[:, None]
         + np.cross(nh, (p0 - p2)) * u1[:, None]
         + np.cross(nh, (p1 - p0)) * u2[:, None]) / a2[:, None]
    gm = np.linalg.norm(g, axis=1)
    out = np.zeros(len(V))
    wt = np.zeros(len(V))
    for k in range(3):
        np.add.at(out, F[:, k], gm * a2)
        np.add.at(wt, F[:, k], a2)
    return out / np.maximum(wt, 1e-12)


def _ordered_boundary_loops(F, nv):
    """Boundary loops as ORDERED vertex rings, walking the one-sided edges.

    An unordered loop is useless for capping: a triangle fan built from a shuffled
    ring folds over itself, and the cotangent weights that come out of those folds
    are meaningless.
    """
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    uq, cnt = np.unique(e, axis=0, return_counts=True)
    be = uq[cnt == 1]
    if len(be) == 0:
        return []
    nb = {}
    for a, b in be:
        nb.setdefault(int(a), []).append(int(b))
        nb.setdefault(int(b), []).append(int(a))
    seen, loops = set(), []
    for start in nb:
        if start in seen:
            continue
        loop, cur, prev = [start], start, None
        seen.add(start)
        while True:
            nxt = [x for x in nb.get(cur, []) if x != prev]
            nxt = [x for x in nxt if x not in seen or x == start]
            if not nxt:
                break
            step = nxt[0]
            if step == start:
                break
            loop.append(step)
            seen.add(step)
            prev, cur = cur, step
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def harmonic_coordinate_capped(mesh, route, is_branch):
    """The harmonic coordinate, solved on the TRUNK with its branch mouths capped.

    A branch mouth is a hole, and a hole carries a no-flux condition, so the level
    sets of the solution have to bend around it. Where several great-vessel mouths sit
    close together on a tight arch that squeezes the level sets into wedges -- which is
    exactly the zone shape that looks wrong on those cases, and it is the field doing
    it, not the zone definition.

    So strip the branches, cap each mouth flat, and solve on a plain tube from inlet to
    outlet. With nothing to bend around, the level sets come back as rings. The caps
    are discarded afterwards; they exist only so the solve sees a closed pipe.

    Returns u on the ORIGINAL vertex numbering, or None if it cannot be built.
    """
    import scipy.sparse.linalg as spl

    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    trunk_f = ~is_branch[F].any(axis=1)
    if trunk_f.sum() < 100:
        return None
    Ft = F[trunk_f]
    keep = np.zeros(len(V), bool)
    keep[Ft.ravel()] = True
    remap = -np.ones(len(V), int)
    remap[keep] = np.arange(int(keep.sum()))
    Vt = V[keep]
    Ftr = remap[Ft]

    loops = _ordered_boundary_loops(Ftr, len(Vt))
    if len(loops) < 2:
        return None
    a, b = route.point(0.0), route.point(route.length)
    cen = [Vt[np.asarray(L_)].mean(axis=0) for L_ in loops]
    ia = int(np.argmin([np.linalg.norm(c - a) for c in cen]))
    ib = int(np.argmin([np.linalg.norm(c - b) for c in cen]))
    if ia == ib:
        return None

    # Cap every hole that is not one of the aorta's own two ends.
    Vc = [Vt]
    Fc = [Ftr]
    n = len(Vt)
    for li, loop in enumerate(loops):
        if li in (ia, ib):
            continue
        ring = np.asarray(loop, int)
        Vc.append(Vt[ring].mean(axis=0).reshape(1, 3))
        c = n
        n += 1
        fan = np.column_stack([np.full(len(ring), c), ring, np.roll(ring, -1)])
        Fc.append(fan)
    Vc = np.vstack(Vc)
    Fc = np.vstack(Fc)

    L = _cotan_laplacian(Vc, Fc)
    fixed = np.zeros(len(Vc), bool)
    val = np.zeros(len(Vc))
    fixed[np.asarray(loops[ia], int)] = True
    fixed[np.asarray(loops[ib], int)] = True
    val[np.asarray(loops[ib], int)] = 1.0
    free = ~fixed
    try:
        rhs = -L[free][:, fixed] @ val[fixed]
        uu = np.zeros(len(Vc))
        uu[fixed] = val[fixed]
        uu[free] = spl.spsolve(L[free][:, free].tocsc(), rhs)
    except Exception:                                        # noqa: BLE001
        return None
    if not np.all(np.isfinite(uu)):
        return None
    uu = np.clip(uu[:len(Vt)], 0.0, 1.0)

    # Back onto the original numbering. Branch vertices take the value of the nearest
    # trunk vertex, which is their mouth -- they are labelled by the ostium rule
    # anyway, so this only has to be sane, not exact.
    u = np.zeros(len(V))
    u[keep] = uu
    if (~keep).any():
        d, j = cKDTree(Vt).query(V[~keep])
        u[~keep] = uu[j]
    return u


def plane_contour_collar(mesh, origin, normal, radius_mm, adj=None):
    """Vertices within *radius_mm* ALONG THE SURFACE of the plane's cross-section.

    Returns (collar_mask, signed_distance) or (None, None) if the plane misses.

    The plane is cut with the wall and only the intersection component nearest
    *origin* is kept -- that is the cross-section at this landmark. A plane runs
    through all of space and on a folded arch meets the wall in more than one ring;
    the others belong to a different stretch of the vessel and must not be allowed to
    decide anything here.

    The collar is measured GEODESICALLY, out along the wall from that contour, never
    as a ball in space. A Euclidean neighbourhood on the arch reaches straight across
    the fold to wall that is 200 mm away along the vessel, and re-deciding those by
    which side of the plane they fall is what shattered the boundary into interleaved
    patches the first time this was tried.
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components, dijkstra

    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    sd = (V - np.asarray(origin, float)) @ np.asarray(normal, float)
    pos = sd > 0
    cross = pos[F].any(axis=1) & (~pos[F]).any(axis=1)
    if not cross.any():
        return None, None

    fa = np.asarray(mesh.face_adjacency)
    ok = cross[fa[:, 0]] & cross[fa[:, 1]]
    ci = np.nonzero(cross)[0]
    pos_of = -np.ones(len(F), int)
    pos_of[ci] = np.arange(len(ci))
    pairs = fa[ok]
    g = sp.coo_matrix((np.ones(len(pairs)),
                       (pos_of[pairs[:, 0]], pos_of[pairs[:, 1]])),
                      shape=(len(ci), len(ci)))
    k, lab = connected_components(g, directed=False)
    cen = np.asarray(mesh.triangles_center)[ci]
    best, bestd = 0, np.inf
    for c in range(k):
        d = np.linalg.norm(cen[lab == c] - origin, axis=1).min()
        if d < bestd:
            best, bestd = c, d
    seed = np.unique(F[ci[lab == best]])

    if adj is None:
        adj = surface_graph(mesh)
    dist = dijkstra(adj, directed=False, indices=seed, limit=float(radius_mm),
                    min_only=True)
    return np.isfinite(dist), sd


def surface_graph(mesh):
    """Vertex adjacency weighted by edge length, for geodesic-ish distances."""
    import scipy.sparse as sp
    V = np.asarray(mesh.vertices, float)
    e = np.asarray(mesh.edges_unique, int)
    w = np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1)
    n = len(V)
    return sp.coo_matrix((np.concatenate([w, w]),
                          (np.concatenate([e[:, 0], e[:, 1]]),
                           np.concatenate([e[:, 1], e[:, 0]]))),
                         shape=(n, n)).tocsr()


# ---------------------------------------------------------------------------
# THE ARCH IS AN ARC, SO USE THE ANGLE.
#
# The clinical zone map draws every boundary as a line converging on one point in the
# concavity of the arch and fanning out to the brachiocephalic, the left carotid, the
# left subclavian. Those are not cross-sections of a tube -- they are RAYS FROM THE
# CENTRE OF CURVATURE, and the zones between them are angular sectors.
#
# That is also the way out of the problem measured everywhere else in this module.
# Nearest-point projection is multivalued once the vessel bends tighter than it is
# wide (k*d > 1 on every arch here), and a plane at a landmark cuts the wall in two or
# three places rather than one. An angle around a fixed point has neither failure: it
# is defined for every point in space, exactly once, however tight the bend.

def fit_arch_circle(route, s_lo=None, s_hi=None, n=200):
    """Best-fit circle through the arch centerline: (centre, plane normal, radius).

    The centre is the point the clinical figure's boundaries converge on. Fitted by
    least squares in the plane of the arch, which is itself taken as the smallest
    principal direction of the centerline points -- an arch is close to planar and the
    third eigenvector is the direction it does not use.
    """
    lo = 0.0 if s_lo is None else float(s_lo)
    hi = float(route.length) if s_hi is None else float(s_hi)
    ss = np.linspace(lo, hi, n)
    P = np.array([route.point(float(s)) for s in ss])
    c0 = P.mean(axis=0)
    _u, _s, vt = np.linalg.svd(P - c0)
    N = vt[2] / np.linalg.norm(vt[2])
    e1, e2 = vt[0], vt[1]
    x = (P - c0) @ e1
    y = (P - c0) @ e2
    # algebraic circle fit: x^2+y^2 + D x + E y + F = 0
    A = np.column_stack([x, y, np.ones(len(x))])
    b = -(x ** 2 + y ** 2)
    try:
        D, E, F = np.linalg.lstsq(A, b, rcond=None)[0]
    except Exception:                                        # noqa: BLE001
        return None, None, None
    cx, cy = -D / 2.0, -E / 2.0
    r2 = cx ** 2 + cy ** 2 - F
    if not np.isfinite(r2) or r2 <= 0:
        return None, None, None
    centre = c0 + cx * e1 + cy * e2
    return centre, N, float(np.sqrt(r2))


def sweep_angle(points, centre, normal, ref):
    """Angle of each point around *centre* in the plane with *normal*, from *ref*.

    Continuous, single-valued, and indifferent to how tightly the vessel bends -- the
    property the projection coordinate does not have.
    """
    r = np.asarray(points, float) - np.asarray(centre, float)
    n = np.asarray(normal, float)
    r = r - np.outer(r @ n, n)
    e1 = np.asarray(ref, float)
    e1 = e1 - (e1 @ n) * n
    e1 /= max(np.linalg.norm(e1), 1e-12)
    e2 = np.cross(n, e1)
    return np.arctan2(r @ e2, r @ e1)


def vessel_end_collar(mesh, route, adj=None, radii=1.0):
    """The trumpet at each end of the aorta, which no zone boundary may re-decide.

    The inlet flares: the aortic root is wider than the vessel above it and its wall
    runs back toward the arch, so it is geodesically CLOSE to a boundary that is 30 mm
    away up the vessel. On KK40_3 that put the root inside the collar of the 0/1
    boundary, its flare fell on the distal side of the fitted plane, and it was handed
    to zone 1 and then to zone 2 -- a band of the wrong zone around the aortic inlet.

    The vessel's own ends belong to the first and last zones by the coordinate, and no
    plane fitted a zone away has anything to say about them.
    """
    from scipy.sparse.csgraph import dijkstra
    from branch_cut import boundary_loops

    V = np.asarray(mesh.vertices, float)
    out = np.zeros(len(V), bool)
    loops = boundary_loops(mesh)
    if len(loops) < 2:
        return out
    a, b = route.point(0.0), route.point(route.length)
    cen = [V[np.asarray(L_)].mean(axis=0) for L_ in loops]
    adj = adj if adj is not None else surface_graph(mesh)
    for end in (a, b):
        k = int(np.argmin([np.linalg.norm(c - end) for c in cen]))
        seed = np.asarray(loops[k], int)
        r = float(np.median(np.linalg.norm(V[seed] - V[seed].mean(axis=0), axis=1)))
        d = dijkstra(adj, directed=False, indices=seed, limit=radii * r,
                     min_only=True)
        out |= np.isfinite(d)
    return out


def flatten_boundary(mesh, label, zi, zj, adj=None, iters=3,
                     collar_radii=1.2, frozen=None, target=None,
                     coord=None, coord_half_width=None, coord_calib=None):
    """Pull the boundary between two zones onto its own best-fit plane.

    The boundary already comes out as ONE closed ring -- measured, on every case -- it
    is simply not flat: 5 to 11 mm off its best-fit plane on the arch, against 2 mm on
    the descending aorta where the vessel is straight. That tilt is the harmonic level
    set meeting the wall obliquely on a bend.

    Fitting a plane to the RING and re-deciding a collar around it is safe in a way
    that the centerline-tangent plane was not. That plane was built from the landmark
    and cut the wall in two or three places; this one is built from the boundary
    itself, so it sits on the wall where the boundary already is, and the collar is
    geodesic so it cannot reach the far limb of a fold.

    Iterated, because moving the ring changes its best fit. Three passes is enough --
    the deviation stops falling after two on every case tried.

    THE COLLAR IS A LOCAL ADJUSTMENT AND HAS TO STAY ONE. Its width is 1.2 times the
    ring's own radius, which is fine on a long zone and is most of a short one: on
    KK40_3 the ascending aorta is 30.8 mm to the brachiocephalic and flares wide at the
    root, so the collar around the 0/1 boundary reached the root, the root's flare fell
    on the distal side of the plane, and 182 vertices of it were handed to zone 1 --
    then to zone 2 by the next boundary's collar, which is why that case had a band of
    zone 2 around the aortic inlet.

    Passing *coord* (the harmonic coordinate) confines the collar to vertices within
    reach of the ring ALONG THE VESSEL rather than across its surface, which is the
    distinction that fails at a flared root: the wall there bulges toward the boundary,
    so it is geodesically close while being 30 mm away up the vessel.

    The width is the same collar_radii ring radii, converted through *coord_calib*, the
    (station, coord) curve from calibrate_u_to_station -- so this is not a second
    tightening knob, it is the existing one measured in the right direction. Giving
    *coord_half_width* directly overrides it; confining the collar instead to half the
    narrower neighbouring zone was tried and starves it, taking KK40_3's 1/2 boundary
    from 2.53 to 15.86 mm off plane.
    """
    from scipy.sparse.csgraph import dijkstra

    V = np.asarray(mesh.vertices, float)
    e = np.asarray(mesh.edges_unique, int)
    if adj is None:
        adj = surface_graph(mesh)
    lab = np.asarray(label, int).copy()
    # BRANCHES ARE NOT UP FOR RE-DECISION. A branch belongs to one zone by its ostium,
    # and that is the rule -- but the collar around a boundary ring reaches across a
    # nearby mouth, and re-deciding those vertices by which side of the plane they
    # fall hands a great vessel to the neighbouring zone. Measured on KR52_1: 425
    # branch vertices moved, zone 0 losing 191 of them to zone 1, and zones 3 and 4
    # acquiring 89 between them when they own no branch at all.
    frozen = (np.zeros(len(V), bool) if frozen is None
              else np.asarray(frozen, bool))

    for _ in range(iters):
        # The ring is fitted on WALL edges. Including the branch vertices on it was
        # tried and measured worse -- flatness median 2.35 -> 3.20 mm, p95 3.31 ->
        # 10.07, max 4.02 -> 23.78 -- because a mouth's vertices sit off in the
        # branch's direction and drag the plane with them.
        wall_e = e[~frozen[e[:, 0]] & ~frozen[e[:, 1]]]
        on = (((lab[wall_e[:, 0]] == zi) & (lab[wall_e[:, 1]] == zj))
              | ((lab[wall_e[:, 0]] == zj) & (lab[wall_e[:, 1]] == zi)))
        ring = np.unique(wall_e[on])
        if len(ring) < 8:
            return lab, None
        P = V[ring]
        if target is not None:
            # THE PLANE OF THE LANDMARK'S OWN CROSS-SECTION, not of wherever the
            # boundary currently sits. Fitting to the current ring makes it flat but
            # locks in whatever displacement it already had: measured against the true
            # cross-sections, flattening moved boundaries from a median +1.1 mm distal
            # to +2.7 mm. The cross-section's plane is flat AND in the right place, so
            # the boundary gets both instead of trading one for the other.
            c = np.asarray(target[0], float)
            nrm = np.asarray(target[1], float)
            nrm = nrm / max(np.linalg.norm(nrm), 1e-12)
        else:
            c = P.mean(axis=0)
            _u, _s, vt = np.linalg.svd(P - c)
            nrm = vt[2]
        # Orient so +n means the zj side, or the two labels would swap.
        side = (V - c) @ nrm
        if np.median(side[lab == zj]) < np.median(side[lab == zi]):
            nrm, side = -nrm, -side
        r = float(np.median(np.linalg.norm(P - c, axis=1)))
        dist = dijkstra(adj, directed=False, indices=ring,
                        limit=collar_radii * r, min_only=True)
        collar = (np.isfinite(dist) & ((lab == zi) | (lab == zj)) & ~frozen)
        if coord is not None:
            cv = np.asarray(coord, float)
            c0 = float(np.median(cv[ring]))
            w = coord_half_width
            if w is None and coord_calib is not None:
                cs, cu = coord_calib
                sb = float(np.interp(c0, cu, cs))
                w = 0.5 * abs(float(np.interp(sb + collar_radii * r, cs, cu))
                              - float(np.interp(sb - collar_radii * r, cs, cu)))
            if w is not None:
                collar &= np.abs(cv - c0) <= w
        lab[collar & (side <= 0)] = zi
        lab[collar & (side > 0)] = zj

    wall_e = e[~frozen[e[:, 0]] & ~frozen[e[:, 1]]]
    on = (((lab[wall_e[:, 0]] == zi) & (lab[wall_e[:, 1]] == zj))
          | ((lab[wall_e[:, 0]] == zj) & (lab[wall_e[:, 1]] == zi)))
    ring = np.unique(wall_e[on])
    if len(ring) < 8:
        return lab, None
    P = V[ring]
    c = P.mean(axis=0)
    _u, _s, vt = np.linalg.svd(P - c)
    return lab, float(np.abs((P - c) @ vt[2]).max())


def smooth_labels(mesh, label, frozen=None, iters=6):
    """Majority vote over surface neighbours, to take the zigzag off a boundary.

    Flattening decides each vertex by which side of a plane it falls, and where that
    plane runs oblique to the wall the sign flips back and forth between neighbours --
    the boundary lands in the plane but arrives there as a staircase. Flat and smooth
    are different properties and the flattening only buys the first.

    A majority filter fixes the second without touching the first: a vertex whose
    neighbours mostly say otherwise is outvoted, which removes single-vertex spurs and
    leaves a boundary that is still on the plane. Frozen vertices never vote and are
    never moved.
    """
    lab = np.asarray(label, int).copy()
    e = np.asarray(mesh.edges_unique, int)
    nv = len(lab)
    frozen = (np.zeros(nv, bool) if frozen is None else np.asarray(frozen, bool))
    live = ~frozen
    for _ in range(iters):
        best = lab.copy()
        counts = {}
        for z in np.unique(lab[lab >= 0]):
            c = np.zeros(nv, int)
            for a, b in ((e[:, 0], e[:, 1]), (e[:, 1], e[:, 0])):
                np.add.at(c, a, (lab[b] == z).astype(int))
            counts[int(z)] = c
        deg = np.zeros(nv, int)
        for a in (e[:, 0], e[:, 1]):
            np.add.at(deg, a, 1)
        top = np.full(nv, -1)
        topn = np.zeros(nv, int)
        for z, c in counts.items():
            take = c > topn
            top[take] = z
            topn[take] = c[take]
        move = live & (top >= 0) & (top != lab) & (topn * 2 > deg)
        if not move.any():
            break
        lab[move] = top[move]
    return lab


def branch_regions(mesh, route, max_mm=40.0, jump=2.2, pullback=1.0):
    """Which vertices belong to which branch stub, found topologically.

    Every distance-based classification tried here fails in the same place: near the
    arch, a great vessel's centerline runs within millimetres of aortic wall, so
    "whose centerline is nearest" -- normalised by radius or not -- hands patches of
    aorta to a branch and patches of branch to the aorta. Verifying anything with that
    mask is circular, and it produced a confident false negative: OPT 3 was reported as
    changing 0 wall vertices when a geometric check found 168 of 193 changed vertices
    on KR6_1 to be aortic.

    The cut already knows where the branches are. Each stub was cut across and so ends
    in its OWN boundary loop, and the stub is the surface between that loop and the
    aorta. Flooding outward from the loop, the frontier stays the width of the tube --
    one circumference per step -- until it reaches the mouth, where it spills onto the
    aortic wall and the frontier suddenly grows. That jump is the mouth, and finding it
    needs no centerline at all.

    Returns an int array: branch index per vertex, -1 for the trunk.
    """
    from scipy.sparse.csgraph import dijkstra
    from branch_cut import boundary_loops

    V = np.asarray(mesh.vertices, float)
    loops = boundary_loops(mesh)
    if len(loops) < 3:
        return -np.ones(len(V), int)

    # The aorta's own two ends are not branches.
    a, b = route.point(0.0), route.point(route.length)
    cen = [V[np.asarray(L_)].mean(axis=0) for L_ in loops]
    ia = int(np.argmin([np.linalg.norm(c - a) for c in cen]))
    ib = int(np.argmin([np.linalg.norm(c - b) for c in cen]))

    adj = surface_graph(mesh)
    out = -np.ones(len(V), int)
    for k, loop in enumerate(loops):
        if k in (ia, ib):
            continue
        seed = np.asarray(loop, int)
        d = dijkstra(adj, directed=False, indices=seed, limit=max_mm, min_only=True)
        fin = np.isfinite(d)
        if fin.sum() < len(seed) * 2:
            continue
        # Frontier size per half-millimetre of geodesic distance. In the tube it is
        # flat; at the mouth it steps up.
        step = 0.5
        edges = np.arange(0.0, max_mm + step, step)
        hist, _ = np.histogram(d[fin], bins=edges)
        base = np.median(hist[1:6]) if len(hist) > 6 else hist[0]
        cut = max_mm
        for i in range(4, len(hist)):
            if base > 0 and hist[i] > jump * base:
                cut = edges[i]
                break
        # STOP SHORT OF THE MOUTH, not past it. The frontier only looks wider once
        # the flood is already out on the aortic wall, so the jump is detected a step
        # late and the region carries a collar of aorta -- which is what moved zones
        # 2 and 3 when the region was forced to one zone. A branch is the part that
        # leaves the aorta, so the region is pulled back by roughly the mouth's own
        # radius. Cutting instead at "one aortic radius from the aortic centerline"
        # was tried and is far worse: the radius estimate is itself inflated by branch
        # vertices, and detection fell from 81 of 90 stubs to 3.
        r_open = float(np.median(np.linalg.norm(V[seed] - V[seed].mean(axis=0),
                                                axis=1)))
        cut = max(cut - pullback * r_open, 2.0 * step)
        out[fin & (d < cut)] = k
    return out


# ---------------------------------------------------------------------------
# THE OSTIUM IS A SADDLE.
#
# Two families of discriminator have been measured to death here. Distance from the
# aortic centerline cannot separate a stub from the wall: a stub that curls back over
# the arch sits at d/r 1.05-1.45 where the wall sits at 1.00 with p90 1.06, so the
# populations overlap and no threshold splits them. Frontier width cannot either --
# swept, it trades one defect for the other along a single curve, from 418 wall
# vertices dragged with 56 stubs still split, to 2934 dragged with 4 split.
#
# Both fail for the same reason: once a stub curls back, nothing about WHERE it is or
# HOW WIDE it is says it is still a branch. What does is the shape of the surface at
# its base. A vessel wall is convex; the join where a branch leaves the aorta is the
# union of two convex tubes, and such a union turns inward all the way round. So the
# mouth is a closed ring where the surface is at its LEAST convex, and it is there
# whether the branch stands up off the arch or lies back along it.
#
# Least convex, not concave: on a smoothed segmentation the fillet rarely goes past
# zero. KR47_1's wall sits at -0.050 with a p90 of -0.013, and the ring its
# brachiocephalic has to be crossed at is about -0.045. Only the RELATIVE maximum
# matters, which is what the bottleneck flood below reads, so nothing here depends on
# the crease being concave in absolute terms.
#
# That is a purely local property of the surface. It needs no centerline, no radius,
# and no distance -- which is what both dead families needed.

def surface_concavity(mesh, smooth=2):
    """Per-vertex concavity, 1/mm: positive where the surface is creased inward.

    Discrete mean curvature normal from the cotangent Laplacian, divided by the
    barycentric vertex area and read against the outward normal, then negated so
    concave is positive. Checked against a sphere of radius 10, where it returns
    -0.1995 against a true -2H of -0.2.

    Smoothed a little because segmentation noise puts a stray concave vertex in the
    middle of a perfectly convex tube, and one such vertex is enough to stop a flood.
    """
    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    L = _cotan_laplacian(V, F)
    fa = np.asarray(mesh.area_faces, float)
    area = np.zeros(len(V))
    for c in range(3):
        np.add.at(area, F[:, c], fa / 3.0)
    hn = (L @ V) / np.maximum(area, 1e-9)[:, None]
    h = np.einsum("ij,ij->i", hn, np.asarray(mesh.vertex_normals, float))
    e = np.asarray(mesh.edges_unique, int)
    deg = np.maximum(np.bincount(e.ravel(), minlength=len(V)).astype(float), 1.0)
    for _ in range(int(smooth)):
        acc = np.zeros(len(V))
        np.add.at(acc, e[:, 0], h[e[:, 1]])
        np.add.at(acc, e[:, 1], h[e[:, 0]])
        h = 0.5 * h + 0.5 * acc / deg
    return -h


def _bottleneck(adjlist, cost, seeds):
    """Least worst crossing: min over paths from *seeds* of the max cost on the path.

    Ordinary shortest path adds cost up and so grows with length, which is distance
    again under another name. This does not: a vertex is cheap if it can be reached
    WITHOUT crossing anything concave, however far away it is, and expensive the moment
    every route to it climbs the ostium ring. Every aortic vertex has to cross that one
    ring, so they all inherit its height and the sorted values step.
    """
    import heapq
    n = len(cost)
    key = np.full(n, np.inf)
    key[seeds] = -np.inf
    pq = [(-np.inf, int(v)) for v in seeds]
    heapq.heapify(pq)
    done = np.zeros(n, bool)
    while pq:
        k, v = heapq.heappop(pq)
        if done[v]:
            continue
        done[v] = True
        for w in adjlist[v]:
            if done[w]:
                continue
            nk = k if cost[w] < k else cost[w]
            if nk < key[w]:
                key[w] = nk
                heapq.heappush(pq, (nk, int(w)))
    return key


NECK_ALPHA = 1.5          # a mouth may be half again the stub's own cut loop
NECK_LEVELS = 80
NECK_WINDOW = 0.25        # never look past a quarter of the mesh for one stub
# How far off the aortic centerline, in local aortic radii, a growth step's vertices
# must sit on average for the step to count as more vessel rather than a spill.
STEP_RATIO_MIN = 1.10
# Below this many vertices added beyond its own cut loop, a stub's flood did not run at
# all and is retried without the length bound. The cohort separates cleanly: the four
# failures add exactly 1 vertex, and the smallest genuine stub adds 17.
MIN_STUB_GROWTH = 5


def _radius_ratio(mesh, route):
    """Each vertex's distance from the aortic centerline, in local aortic radii."""
    ss, rr = aortic_radius_profile(mesh, route)
    P = np.array([route.point(float(s)) for s in ss])
    d, i = cKDTree(P).query(np.asarray(mesh.vertices, float))
    return d / np.maximum(rr[i], 1e-9)


def _fill_pockets(adj, sel, keep_out):
    """Absorb islands of unselected surface enclosed by the region.

    Segmentation noise leaves single concave vertices inside a stub and the flood goes
    round them, so the region comes out freckled. A pocket of unselected surface that
    cannot reach the rest of the mesh is interior by definition.
    """
    from scipy.sparse.csgraph import connected_components
    out = ~sel
    if not out.any():
        return sel
    idx = np.flatnonzero(out)
    ncmp, lab = connected_components(adj[idx][:, idx], directed=False)
    if ncmp <= 1:
        return sel
    pos = np.searchsorted(idx, np.intersect1d(idx, np.asarray(keep_out, int)))
    outside = set(lab[pos].tolist())
    if not outside:
        return sel
    sel = sel.copy()
    sel[idx[~np.isin(lab, list(outside))]] = True
    return sel


def stub_regions(mesh, route, alpha=NECK_ALPHA, smooth=2, adj=None, ratio=None):
    """Which vertices belong to which branch stub, decided by the ostium crease.

    Each stub was cut across, so it ends in its own boundary loop. Flooding out from
    that loop over the concavity field orders the surface by the worst crease that has
    to be climbed to reach it: the stub is everything reachable without climbing its
    own mouth.

    Where to stop is not a threshold on concavity -- that number varies with calibre
    and with mesh resolution -- but a statement about the region's own shape. While the
    region is the stub, its boundary is a ring the size of the stub; the instant it
    spills, the boundary is an aortic ring several times longer. So the level is raised
    until the boundary first exceeds *alpha* times the stub's own cut loop, and the
    last level before that is kept. Both lengths belong to the stub, so the rule
    carries no length scale of its own.

    FIRST FAILURE, NOT LAST SUCCESS. Raise the level far enough and the region is the
    whole mesh, whose boundary is nothing at all, so the bound comes good again on the
    far side of the leak. The mouth is where it is first broken.

    AND THE STEP MUST ADD VESSEL. The length bound asks whether the rim is short; it
    cannot ask whether the rim is in the right place, and on a vessel with a wide
    opening -- KR47_1's brachiocephalic cut loop is 80 mm round -- half again is
    generous enough to let a splash of aortic wall through underneath it. So each step
    is also judged by WHAT IT ADDS: the median distance of its new vertices from the
    aortic centerline, in local aortic radii.

    That is the radius test, which is dead as a per-vertex classifier and was measured
    dead -- a curled stub sits at 1.05-1.45 and the wall at 1.00 with p90 1.06, so the
    populations overlap and no threshold splits them. A MEDIAN OVER A STEP is not the
    same statistic. Each step adds dozens of vertices at once, and on KR47_1 every
    legitimate step adds at a median of 1.18 to 1.51 while the step that lays down the
    splash adds at 1.02. It is a question about a set, not about a vertex.

    Returns an int array: stub index per vertex, -1 for the trunk.
    """
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

    conc = surface_concavity(mesh, smooth)
    e = np.asarray(mesh.edges_unique, int)
    ew = np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1)
    adjlist = [[] for _ in range(len(V))]
    for i, j in e:
        adjlist[i].append(j)
        adjlist[j].append(i)
    adj = adj if adj is not None else surface_graph(mesh)
    ratio = _radius_ratio(mesh, route) if ratio is None else np.asarray(ratio, float)
    far = np.concatenate([np.asarray(loops[ia], int), np.asarray(loops[ib], int)])

    for j, k in enumerate([x for x in range(len(loops)) if x not in (ia, ib)]):
        seed = np.asarray(loops[k], int)
        key = _bottleneck(adjlist, conc, seed)
        fin = np.isfinite(key)
        if fin.sum() < len(seed) * 2:
            continue
        bs = np.sort(key[fin])
        onloop = np.zeros(len(V), bool)
        onloop[seed] = True
        lopen = float(np.linalg.norm(V[seed] - V[np.roll(seed, 1)], axis=1).sum())
        # LEVELS BY RANK, NOT BY VALUE. Concavity is nearly constant through the tube
        # and then climbs steeply at the mouth, so a ladder spaced evenly in
        # concavity spends itself on the flat part: measured on KR47_1's
        # brachiocephalic, 25 of 40 rungs left the region at 33 vertices and the last
        # two rungs took it from 154 to 2704. The mouth was between them and could not
        # be seen. Spacing the rungs evenly in RANK puts a roughly equal number of new
        # vertices on each, so the resolution is where the region is actually growing.
        hi = max(2, int(NECK_WINDOW * len(bs)))
        rungs = np.unique(np.linspace(0, hi - 1, NECK_LEVELS).astype(int))

        def climb(bound, forbid=None):
            """Raise the level until a test stops it, and keep the last level that held.

            *forbid* are openings this region may not reach: the aorta's own two ends
            and every OTHER stub's cut loop. Reaching one is not a matter of degree --
            a stub that has arrived at another vessel's opening has left its own vessel,
            whatever its rim happens to measure.
            """
            prev = onloop.copy()
            kept = None
            for i in rungs:
                sel = onloop | (fin & (key <= bs[i]))
                add = sel & ~prev
                if add.sum() >= 5 and float(np.median(ratio[add])) < STEP_RATIO_MIN:
                    break
                if forbid is not None and sel[forbid].any():
                    break
                cut = sel[e[:, 0]] != sel[e[:, 1]]
                if float(ew[cut].sum()) > bound * max(lopen, 1e-9):
                    break
                kept = sel
                prev = sel
            return kept

        best = climb(alpha)
        # A REGION THE SIZE OF ITS OWN CUT LOOP IS NOT A SMALL STUB, IT IS A FAILED
        # FLOOD. The length bound weighs the region's rim against the loop the flood
        # started from, and on a vessel that flares just below where it was cut the
        # very first rung already carries a rim half again as long -- so the ladder
        # breaks before it has reported anything and the region is the seed loop
        # itself. Its rim is then the cut opening partway up the vessel instead of the
        # ostium, and everything downstream reads the wrong mouth: on KK16_1 the
        # subclavian's distal lip came out 11 mm past where it is, which is exactly how
        # far the zone 2/3 boundary then sat beyond the branch.
        #
        # That failure announces itself -- the region did not grow -- so it can be
        # retried rather than met by loosening the bound for every stub. Over the
        # cohort 31 of 90 stubs stop on the length test and NOT ONE of them stops on a
        # step whose new vertices look like aortic wall; the ratio test catches all 59
        # leaks unaided. So the retry drops the length bound and leaves the ratio test
        # to stop the flood. It fires on four stubs in thirty cases -- KK16_1, KK32_0,
        # KP35_1, KR6_1 -- and every other region is untouched.
        #
        # The retry needs a bound of its own, because the ratio test cannot supply one
        # here: it asks whether a step's new vertices stand off the aortic centerline
        # like a vessel, and a NEIGHBOURING VESSEL answers yes. Left unbounded on
        # JK5_1, one stub ran 54 vertices to 846, straight over another stub and
        # erasing its region. So the retry may not reach another opening.
        if best is None or int(best.sum()) - len(seed) < MIN_STUB_GROWTH:
            others = np.concatenate(
                [far] + [np.asarray(loops[x], int)
                         for x in range(len(loops)) if x not in (ia, ib, k)])
            grown = climb(float("inf"), forbid=others)
            if grown is not None and int(grown.sum()) - len(seed) >= MIN_STUB_GROWTH:
                best = grown
        if best is None:
            continue
        out[_fill_pockets(adj, best, far)] = j
    return out


def stub_ostium_zone(tree, route, mesh, regions, zones, jg=None):
    """The zone each stub region's OWN OSTIUM puts it in.

    own_stubs used to read a region's zone off whichever of its vertices station_field
    had called a branch, and fall back to the surrounding wall when there were none.
    Both are proxies for the rule, and the fallback contradicts it outright: it gives a
    vessel the zone of the wall around its mouth, which is the zone the streak came
    from. On a bovine arch that is not a corner case but the normal one -- the shared
    trunk's daughters lie along the arch, station_field calls them trunk, and both
    KK5_1 and KY244_1 handed a daughter to zone 2 when its ostium is in zone 0.

    So ask the centerline directly. Every off-route segment already carries the station
    of the junction where its subtree meets the route -- one value for the whole branch,
    which is what makes a branch indivisible by a zone boundary -- and the nearest such
    segment to a region is the vessel that region is part of. The median over the
    region's vertices, so a few vertices near a neighbouring mouth cannot decide it.

    Returns {region index: zone number}, missing where no off-route centerline is near.
    """
    jg = jg or _junction_graph(tree)
    P, _T, S, O, _R = _segment_station_samples(tree, route, jg, mesh)
    off = ~np.isin(O, list(route.chain))
    if not off.any():
        return {}
    kd = cKDTree(P[off])
    s_off = S[off]
    V = np.asarray(mesh.vertices, float)
    L = float(route.length)
    reg = np.asarray(regions, int)
    out = {}
    for j in range(int(reg.max()) + 1):
        sel = reg == j
        if not sel.any():
            continue
        _d, i = kd.query(V[sel])
        st = float(np.median(s_off[i]))
        if st < -1e8:
            continue
        for p in zones:
            if branch_inside(np.array([st]), p["s0"], p["s1"], L)[0]:
                m = re.search(r"(\d+)\s*$", p["name"] or "")
                out[j] = int(m.group(1)) if m else 0
                break
    return out


def landmark_lag(mesh, coord, calib, regions, region_zones, zone=2):
    """How far the landmark station for a zone's vessel is PAST that vessel's mouth.

    A zone boundary is asked for at the station where the branch's centerline meets the
    aorta's. On most arches that is the ostium: measured over this cohort the left
    subclavian's landmark sits a median 0.1 mm from its own mouth. On a few it is not,
    and it is late by up to 14 mm -- the branch centerline runs alongside the aorta
    before it diverges, so the junction is placed well past where the vessel actually
    leaves the wall. Zone 2 then runs on past the subclavian and zone 3 starts late,
    which is exactly what that looks like on the screen.

    The mouth is not an estimate: the crease regions find it on the wall. So return the
    lag, landmark minus mouth, for the caller to take out of the bounds.

    Positive means the landmark is distal of the mouth, which is the failure. Returns
    None when there is no stub for the zone to measure against.
    """
    reg = np.asarray(regions, int)
    js = [j for j, z in region_zones.items() if z == zone]
    if not js or reg.max() < 0:
        return None
    e = np.asarray(mesh.edges_unique, int)
    sel = np.isin(reg, js)
    if not sel.any():
        return None
    cut = sel[e[:, 0]] != sel[e[:, 1]]
    rim = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
    rim = rim[~sel[rim]]                    # the wall side: the ring it leaves from
    if len(rim) < 6:
        return None
    ss, uu = calib
    return float(np.median(np.interp(np.asarray(coord, float)[rim], uu, ss)))


# ANCHOR ON THE MOUTH'S PROXIMAL LIP, NOT ITS DISTAL ONE. The distal lip is where the
# vessel finally lets go of the aorta, and on a shallow takeoff it trails a long way
# downstream: measured across the marked cases, the mouth spans 7.0, 7.9 and 8.8 mm of
# centerline where the takeoff is square, and 15.0 and 19.7 mm where it is oblique. A
# cut at the distal lip is therefore exact on a compact mouth and late in proportion to
# the elongation -- overshooting the marks by 0.5, 1.2 and 1.0 mm on the first three,
# and by 5.5 and 10.6 mm on the last two. That is the whole of the "zone 3 starts too
# late" complaint, and it is not a constant, so no constant nudge ever fixed it.
#
# The PROXIMAL lip does not move with obliquity -- it is where the vessel first parts
# from the aorta. Measured from there the marks sit at +6.5, +6.8, +7.8, +9.2 and
# +9.4 mm: one number, sd 1.3, against sd 4.2 for the distal lip. So the cut is placed
# a mouth's width past where the vessel leaves, which is what "immediately after the
# LSA" means on an arch where the ostium is drawn out along the wall.
#
# Scored against 167 marks: proximal, mean 3.0 mm / trimmed 1.7 / median 1.5, against
# 5.1 / 3.5 / 2.8 for the distal lip. No marked case is worse; the two the plane
# already fitted are unchanged to the decimal, and the two worst roughly halve.
# SUPERSEDED, AND LEFT HERE ONLY SO THE SWEEPS IN tools/zone_score.py STILL RUN. The
# shipped path is mouth_clearing_plane, which takes no percentile and no offset; these
# two constants and mouth_plane() below no longer affect any map. Both anchors this pair
# can express were measured against eight cases judged by eye and each failed in both
# directions at once, because a lip plus a constant cannot track a mouth whose own span
# runs from 1 to 20 mm of centerline.
MOUTH_PCT = 10.0
MOUTH_OFFSET_MM = 8.0
PLANE_REACH_MM = 70.0     # how far the plane is allowed to decide anything


def mouth_plane(mesh, route, tree, coord, regions, region_ids,
                       pct=MOUTH_PCT, offset_mm=MOUTH_OFFSET_MM):
    """The CENTERLINE-NORMAL PLANE just past a vessel's mouth.

    This is the boundary the clinical figure draws: a cut across the vessel, with the
    branch sitting at the distal end of its zone and the next zone beginning immediately
    past it. Everything tried before moved a LEVEL SET of the harmonic coordinate
    instead, and a level set follows the wall -- on a bend it runs oblique, so shifting
    it along the vessel changes where it sits without ever making it a transverse cut.
    That is why moving it never produced what the figure shows, whichever way it moved.

    So take a lip of the mouth, project it onto the centerline, and cut there. A plane
    needs no coordinate, no calibration and no level set, and it is the same plane the
    branch cuts already use.

    WHICH lip is the whole question, and it is the PROXIMAL one -- see MOUTH_PCT. The
    distal lip trails downstream as the takeoff gets shallower, so a cut there is late
    by however oblique the ostium happens to be; the proximal lip stays put. *pct* is
    which rim vertex stands for it, by coordinate rather than by the single extreme, so
    one stray vertex cannot place the boundary. *offset_mm* then carries the cut a
    mouth's width past it, which is where the marks put the line.

    Returns (origin, normal, station) with the normal pointing DISTAL, or None.
    """
    reg = np.asarray(regions, int)
    sel = np.isin(reg, list(region_ids))
    if not sel.any():
        return None
    e = np.asarray(mesh.edges_unique, int)
    cut = sel[e[:, 0]] != sel[e[:, 1]]
    rim = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
    rim = rim[~sel[rim]]                       # the ring of WALL the vessel leaves from
    if len(rim) < 6:
        return None
    V = np.asarray(mesh.vertices, float)
    cv = np.asarray(coord, float)[rim]
    # The harmonic coordinate is monotone along the vessel, so it ranks the rim from
    # proximal to distal without any of the projection's failure modes. Only the ONE
    # vertex it picks is then projected, and a single point projects cleanly even where
    # a whole ring would not.
    lip = rim[int(np.argmin(np.abs(cv - np.percentile(cv, pct))))]
    ss = np.arange(0.0, float(route.length) + 1e-9, 0.25)
    P = np.array([route.point(float(x)) for x in ss])
    st = float(ss[int(np.argmin(np.linalg.norm(P - V[lip], axis=1)))])
    st = min(st + float(offset_mm), float(route.length) - 1e-3)
    name, local, _f = route.locate(st)
    c, t = tree.segments[name].at_arclength(float(local))
    t = np.asarray(t, float)
    n = t / max(np.linalg.norm(t), 1e-12)
    # route.point runs from inlet to outlet, so +tangent is distal. Confirm it against
    # the route itself rather than trusting the segment's own direction, which follows
    # however the centerline was traced.
    ahead = route.point(min(st + 1.0, float(route.length)))
    if float(np.dot(np.asarray(ahead, float) - np.asarray(c, float), n)) < 0:
        n = -n
    return np.asarray(c, float), n, st


MOUTH_CLEARANCE_MM = 0.0   # how far past the mouth the cut sits, once it has cleared it


def mouth_clearing_plane(mesh, route, regions, region_ids,
                         clearance_mm=MOUTH_CLEARANCE_MM, step=0.25):
    """The most proximal centerline-normal plane that leaves the WHOLE mouth behind it.

    Two anchors were tried before this and both were a lip of the mouth plus a constant:
    the distal lip (exact on a square takeoff, late in proportion to obliquity) and the
    proximal lip plus 8 mm (right only where the mouth happens to span 8 mm). The second
    was fitted to five hand-marked cases and then failed on eight judged by eye, in both
    directions at once -- overshooting the branch where the mouth is short and cutting
    through it where the mouth is long. A constant cannot track a mouth whose own span
    runs from 1 to 20 mm of centerline.

    The quantity that does track it was found by measuring the eight judgements against
    every candidate: the plane's SIGNED CLEARANCE over the ostium. Called correct at
    +0.05 mm; "encroaching on the LSA" at -17.6; branch faces showing through at -5.2
    and -7.3; "way too far" at +9.2; "should be closer to that last branch" at +13.7 and
    +18.6. One number, monotone, no overlap between the verdicts.

    So ask for clearance zero directly instead of computing a lip and hoping a constant
    lands there. Sweep the centerline and keep the first station whose normal plane puts
    every rim vertex on its proximal side. This needs no percentile, no offset, and no
    harmonic coordinate -- which also means it still works on the cases where the capped
    solve returns nothing and the old rule silently did not run at all.

    THE PLANE THEREFORE NEVER CUTS THE OSTIUM, which is the same statement as the branch
    keeping its own zone. The zone-3 faces that showed through on a branch were faces
    with one corner in the mouth and one corner out past a plane that had sliced between
    them; a plane that clears the mouth has no such face to make.

    Returns (origin, normal, station) with the normal pointing DISTAL, or None.
    """
    reg = np.asarray(regions, int)
    sel = np.isin(reg, list(region_ids))
    if not sel.any():
        return None
    e = np.asarray(mesh.edges_unique, int)
    cut = sel[e[:, 0]] != sel[e[:, 1]]
    if not cut.any():
        return None
    rim = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
    rim = rim[~sel[rim]]                       # the ring of WALL the vessel leaves from
    if len(rim) < 6:
        return None
    V = np.asarray(mesh.vertices, float)
    Q = V[rim]
    ss = np.arange(0.0, float(route.length) + 1e-9, float(step))
    P = np.array([route.point(float(x)) for x in ss])
    # Tangents by central difference on the sampled route. The station itself is never
    # used as a measurement here -- only as the parameter being swept -- so the
    # projection that is unreliable on a tight arch is not in the loop at all.
    T = np.empty_like(P)
    T[1:-1] = P[2:] - P[:-2]
    T[0] = P[1] - P[0]
    T[-1] = P[-1] - P[-2]
    T /= np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-12)
    # signed distance of every rim vertex from the plane at every candidate station
    sd = np.einsum("ij,ikj->ik", T, Q[None, :, :] - P[:, None, :])
    ok = (sd <= -float(clearance_mm)).all(axis=1)
    if not ok.any() or bool(ok[0]):
        # ok[0] true would mean the whole mouth is behind the inlet plane, which is not
        # a mouth on this route; better to decline than to place a confident wrong cut.
        return None
    i = int(np.argmax(ok))
    return P[i].copy(), T[i].copy(), float(ss[i])


def split_by_plane(mesh, label, za, zb, origin, normal, adj=None, frozen=None,
                   near=None, reach_mm=PLANE_REACH_MM):
    """Re-decide the za/zb boundary as a plane cut: proximal is za, distal is zb.

    A plane is global and an aorta is folded, so on a tight arch one plane meets the
    wall in more than one place. Two things keep that from mattering. Only vertices
    ALREADY in za or zb are re-decided, and those are a band about the arch rather than
    the whole vessel; and *near* plus *reach_mm* confine it to the geodesic
    neighbourhood of the mouth, so a second intersection somewhere down the descending
    aorta has nothing there to change.

    Branches are frozen: a stub belongs to the zone its ostium is in, and a plane
    sweeping past its mouth must not take it.
    """
    from scipy.sparse.csgraph import dijkstra

    lab = np.asarray(label, int).copy()
    V = np.asarray(mesh.vertices, float)
    e = np.asarray(mesh.edges_unique, int)
    frozen = (np.zeros(len(lab), bool) if frozen is None
              else np.asarray(frozen, bool))
    live = ((lab == za) | (lab == zb)) & ~frozen
    # SEEDED FROM THE BOUNDARY ITSELF, not only from the mouth. A boundary is a ring
    # around the whole vessel and the mouth is a patch on one side of it, so the far
    # side of the ring is most of a circumference away -- past any sane reach. Seeding
    # only from the mouth left that far side on its old level set, and the "plane" came
    # out 12 mm off plane instead of one edge length.
    on = (((lab[e[:, 0]] == za) & (lab[e[:, 1]] == zb))
          | ((lab[e[:, 0]] == zb) & (lab[e[:, 1]] == za)))
    seeds = [np.unique(e[on])] if on.any() else []
    if near is not None and len(np.atleast_1d(near)):
        seeds.append(np.asarray(near, int).ravel())
    if seeds:
        adj = adj if adj is not None else surface_graph(mesh)
        d = dijkstra(adj, directed=False, indices=np.concatenate(seeds),
                     limit=float(reach_mm), min_only=True)
        live &= np.isfinite(d)
    if not live.any():
        return lab, 0
    side = (V - np.asarray(origin, float)) @ np.asarray(normal, float)
    before = lab.copy()
    lab[live & (side <= 0)] = za
    lab[live & (side > 0)] = zb
    return lab, int((lab != before).sum())


def shift_zone_bounds(zones, mm, from_zone=3):
    """Slide a zone and everything after it *mm* PROXIMALLY along the centerline.

    The complaint this answers is that zone 3 starts too late -- there is aorta between
    the left subclavian and the start of the zone that is meant to begin at it. What is
    wanted is the band moved a little way back, keeping its shape and its 20 mm length.

    TWO EARLIER SHAPES OF THIS WERE WRONG, both because they changed the band instead of
    moving it. Re-placing only the zone 2 / zone 3 boundary GREW zone 3 rather than
    sliding it, and since a boundary is a level set of the harmonic coordinate, a fixed
    change of coordinate is a fixed change of nothing in particular: it moved a hair on
    the outer curvature and swept the whole underside of KK32_0's arch. Eroding zone 2
    by a geodesic collar was uniform in millimetres of wall but not in what it removed:
    a zone on a bend is a wedge, wide on the outer curvature and narrow on the inner, so
    a constant collar ate the whole inner side and left zone 2 with no bottom.

    Moving the BOUNDS moves the band. Zone 2 ends earlier, zone 3 starts and ends
    earlier by the same amount so it keeps its length, and zone 4 starts earlier -- all
    in arc length along the vessel, which is the direction the zones are defined in.
    Every boundary is still placed and flattened exactly as before; they are simply
    asked for at a different station.
    """
    out = []
    for p in zones:
        q = dict(p)
        m = re.search(r"(\d+)\s*$", q.get("name") or "")
        z = int(m.group(1)) if m else 0
        if z >= from_zone and q.get("s0") is not None:
            q["s0"] = float(q["s0"]) - mm
        if z >= from_zone - 1 and q.get("s1") is not None:
            q["s1"] = float(q["s1"]) - mm
        out.append(q)
    return out


def nudge_boundary(mesh, label, za, zb, mm, adj=None, frozen=None):
    """Move the za/zb boundary by *mm* of WALL, zb growing into za.

    A zone boundary is a level set of the harmonic coordinate, so shifting it by a fixed
    amount OF COORDINATE moves it by wildly different amounts of wall: barely anything
    on the outer curvature where the level sets are far apart, a great deal on the inner
    curvature where they crowd. On KK32_0 a shift worth 8.7 mm of arc length swept zone 3
    across the whole underside of the arch.

    A geodesic collar does not do that. Every vertex of za within *mm* of the boundary
    goes to zb, so the boundary moves the same distance everywhere, which is what
    "move the band over" means.

    Branches are left alone: a stub belongs to the zone its ostium is in, and sliding a
    boundary past its mouth must not take it.
    """
    from scipy.sparse.csgraph import dijkstra

    lab = np.asarray(label, int).copy()
    if mm <= 0:
        return lab, 0
    e = np.asarray(mesh.edges_unique, int)
    adj = adj if adj is not None else surface_graph(mesh)
    frozen = (np.zeros(len(lab), bool) if frozen is None
              else np.asarray(frozen, bool))
    on = (((lab[e[:, 0]] == za) & (lab[e[:, 1]] == zb))
          | ((lab[e[:, 0]] == zb) & (lab[e[:, 1]] == za)))
    ring = np.unique(e[on])
    ring = ring[(lab[ring] == za) | (lab[ring] == zb)]
    if len(ring) < 6:
        return lab, 0
    d = dijkstra(adj, directed=False, indices=ring, limit=float(mm), min_only=True)
    move = np.isfinite(d) & (lab == za) & ~frozen
    lab[move] = zb
    return lab, int(move.sum())


def stub_rim_coord(mesh, coord, regions, pct=95.0, wall_side=True):
    """The coordinate value at each stub region's own rim -- i.e. at its ostium.

    A zone boundary is not an arc length, it is a VESSEL: zone 0 ends at the
    brachiocephalic, zone 1 at the left carotid, zone 2 at the left subclavian. Placing
    it by interpolating the station->coordinate curve asks the centerline where that
    vessel is and then asks the curve where that station is, and both steps are biased
    on a tight arch: measured on KK40_3, the 1/2 boundary landed 15 mm distal of the
    left carotid with a 48 mm spread, so zone 1 never closed around the trunk.

    The stub region's rim IS the ostium, found from the crease and nothing else. Its
    median coordinate is where the boundary belongs, with no curve in between.

    *pct* picks where on the rim: 95 is the DISTAL LIP, which is what a zone boundary
    wants. A zone owns the vessel at its distal end and the next zone begins immediately
    after it, so the boundary belongs just past the far edge of the mouth. Taking the
    rim's median instead puts it through the middle of the mouth, and every boundary
    placed that way came out distal of where it should be.

    *wall_side* takes the rim vertices belonging to the AORTA rather than to the stub.
    They are the ring of wall the vessel leaves from, and they carry the coordinate the
    surrounding wall is labelled by.

    Returns {region index: coordinate on its rim}.
    """
    V = np.asarray(mesh.vertices, float)
    e = np.asarray(mesh.edges_unique, int)
    reg = np.asarray(regions, int)
    cv = np.asarray(coord, float)
    out = {}
    for j in range(int(reg.max()) + 1):
        sel = reg == j
        if not sel.any():
            continue
        cut = sel[e[:, 0]] != sel[e[:, 1]]
        rim = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
        rim = rim[~sel[rim]] if wall_side else rim[sel[rim]]
        if len(rim) < 6:
            continue
        out[j] = float(np.percentile(cv[rim], pct))
    return out


def own_stubs(mesh, label, regions, is_branch, adj=None, region_zone=None):
    """Give every stub the zone its own ostium already chose, and tidy up after.

    THE STUB ALREADY HAS A ZONE. Its ostium station decided it, through branch_inside,
    and all this does is spread that answer over the whole stub, so a vessel cannot end
    up partly in a zone its mouth is not in. Nothing new is decided here.

    *region_zone* is that answer taken straight from the centerline, one zone per
    region, and is what should be passed: see stub_ostium_zone. Without it the zone is
    read off whichever of the region's vertices station_field called a branch, which is
    a proxy, and where there are none it falls back to the wall around the mouth --
    which is the zone the streak came from, and is how both bovine cases handed a
    daughter of the shared trunk to zone 2 when its ostium is in zone 0.

    Called AFTER the boundaries are flattened, so it can only ever repaint a stub: the
    wall between two zones has already been settled and no vertex of it is in a region.

    THEN REPAIR, WITH THE STUBS FROZEN. Taking a stub back can strand a few vertices of
    the zone it was borrowed from -- on KR52_1 an island of 8 -- and an island of one
    zone sitting inside another is a defect whoever made it. The repair is the ordinary
    connectivity pass, with the stub regions frozen so that it cannot undo the very
    ownership it is cleaning up after; without that guard it hands the stub straight
    back to whichever zone surrounds its mouth.
    """
    lab = np.asarray(label, int).copy()
    reg = np.asarray(regions, int)
    isb = np.asarray(is_branch, bool)
    for j in range(int(reg.max()) + 1):
        sel = reg == j
        if not sel.any():
            continue
        if region_zone is not None and j in region_zone:
            lab[sel] = int(region_zone[j])
            continue
        cand = lab[sel & isb]
        cand = cand[cand >= 0]
        if not len(cand):
            cand = lab[sel][lab[sel] >= 0]
        if len(cand):
            lab[sel] = int(np.bincount(cand).argmax())
    lab, _moved = enforce_connected_labels(mesh, lab, frozen=(reg >= 0))
    return lab


def stub_zones(mesh, route, regions, openings, land, zones):
    """Zone per branch region: its own opening names it, the landmark table places it.

    Deriving the station from the stub's own geometry was tried and is biased distal
    by 4-10 mm, which is a whole zone on an arch. The zone BOUNDS come from
    aorta_landmark_stations, so membership has to be read from that same table or the
    two disagree by construction.

    Regions whose opening matches no landmark -- a trim group, an unnamed vessel --
    return None and are left to whatever labelled them, rather than being forced
    somewhere on a guess.
    """
    V = np.asarray(mesh.vertices, float)
    if not openings:
        return {}
    cen = np.array([o.get("centroid") or o.get("plane_origin")
                    for o in openings], float)
    L = float(route.length)
    out = {}
    for k in sorted(set(np.asarray(regions).tolist()) - {-1}):
        sel = np.asarray(regions) == k
        if not sel.any():
            continue
        j = int(np.argmin(np.linalg.norm(cen - V[sel].mean(axis=0), axis=1)))
        nm = openings[j].get("name") or ""
        st = land.get(nm)
        if st is None:
            hits = [v for kk, v in land.items()
                    if kk.lower() in nm.lower() or nm.lower() in kk.lower()]
            st = hits[0] if hits else None
        if st is None:
            out[k] = None
            continue
        z = [zi for zi, p in enumerate(zones)
             if branch_inside(np.array([float(st)]), p["s0"], p["s1"], L)[0]]
        out[k] = z[0] if z else None
    return out


def boundary_u_at(mesh, route, u, station, seg=None, prior=None):
    """The harmonic value ON the landmark's own cross-section.

    The calibration curve answers "what is u at station s" by sampling balls along the
    centerline and interpolating. On the arch that is biased: measured against the true
    cross-sections, boundaries placed from it sit a median +3.0 mm DISTAL of their
    landmark and up to 12.6 mm, so every zone grows into the next -- while the one
    boundary in the straight descending aorta lands at +0.1 mm.

    The cross-section is a direct answer to the same question. Take the plane at the
    landmark, keep the intersection component nearest that centerline point -- the
    others belong to a different stretch of a folded vessel -- and read the median u on
    it. The boundary then lies on the landmark's own cross-section by construction
    rather than by interpolation.

    Returns None if the plane misses, so the caller can fall back to the curve.
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    name, local, _fwd = route.locate(float(station))
    segment = (route.tree.segments[name] if seg is None else seg)
    c, t = segment.at_arclength(float(local))
    c = np.asarray(c, float)
    t = np.asarray(t, float)

    sd = (V - c) @ t
    pos = sd > 0
    cross = pos[F].any(axis=1) & (~pos[F]).any(axis=1)
    if not cross.any():
        return None
    fa = np.asarray(mesh.face_adjacency)
    ok = cross[fa[:, 0]] & cross[fa[:, 1]]
    ci = np.nonzero(cross)[0]
    po = -np.ones(len(F), int)
    po[ci] = np.arange(len(ci))
    pr = fa[ok]
    g = sp.coo_matrix((np.ones(len(pr)), (po[pr[:, 0]], po[pr[:, 1]])),
                      shape=(len(ci), len(ci)))
    k, lab = connected_components(g, directed=False)
    cen = np.asarray(mesh.triangles_center)[ci]
    # WHICH CONTOUR IS THE CROSS-SECTION. Nearest in space is wrong on a fold: at
    # station 111.6 on KK20_2 it picked a ring on the descending limb and returned
    # u = 0.376 where the next landmark 20 mm further on returned 0.278 -- u running
    # backwards, so that zone's span inverted and collapsed. u must increase along the
    # vessel, so when a prior is available the component whose value is CLOSEST TO IT
    # is the cross-section, and the prior only has to be good to a few millimetres.
    best, bd = 0, np.inf
    for q in range(k):
        ring_q = np.unique(F[ci[lab == q]])
        if len(ring_q) < 8:
            continue
        if prior is None:
            score = float(np.linalg.norm(cen[lab == q] - c, axis=1).min())
        else:
            score = abs(float(np.median(u[ring_q])) - float(prior))
        if score < bd:
            best, bd = q, score
    ring = np.unique(F[ci[lab == best]])
    if len(ring) < 8:
        return None
    return float(np.median(u[ring]))


def cross_section_plane(mesh, route, tree, station):
    """(centroid, normal) of the landmark's own cross-section contour, or None."""
    V = np.asarray(mesh.vertices, float)
    name, local, _f = route.locate(float(station))
    c, t = tree.segments[name].at_arclength(float(local))
    u_dummy = np.zeros(len(V))
    ring = _cross_section_ring(mesh, np.asarray(c, float), np.asarray(t, float))
    if ring is None or len(ring) < 8:
        return None
    P = V[ring]
    cq = P.mean(axis=0)
    _a, _b, vt = np.linalg.svd(P - cq)
    return cq, vt[2]


def _cross_section_ring(mesh, c, t):
    """Vertices of the plane's intersection component nearest c."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    sd = (V - c) @ t
    pos = sd > 0
    cr = pos[F].any(axis=1) & (~pos[F]).any(axis=1)
    if not cr.any():
        return None
    fa = np.asarray(mesh.face_adjacency)
    ok = cr[fa[:, 0]] & cr[fa[:, 1]]
    ci = np.nonzero(cr)[0]
    po = -np.ones(len(F), int)
    po[ci] = np.arange(len(ci))
    pr = fa[ok]
    g = sp.coo_matrix((np.ones(len(pr)), (po[pr[:, 0]], po[pr[:, 1]])),
                      shape=(len(ci), len(ci)))
    k, lab = connected_components(g, directed=False)
    cen = np.asarray(mesh.triangles_center)[ci]
    best, bd = 0, np.inf
    for q in range(k):
        d = np.linalg.norm(cen[lab == q] - c, axis=1).min()
        if d < bd:
            best, bd = q, d
    return np.unique(F[ci[lab == best]])


def aortic_radius_profile(mesh, route, step=1.0):
    """Aortic radius along the route, read off the vessel's own cross-sections.

    Estimating it as "distance from the centerline to the nearest wall" is what broke
    the earlier attempt to tell stub from aorta: branch vertices sit near the
    centerline too, so the estimate is contaminated by exactly the surface it is meant
    to exclude, and detection fell from 81 of 90 stubs to 3.

    A cross-section cannot be contaminated that way. The plane's intersection component
    nearest the centerline point is the aorta's own ring, and the branches are other
    components. Where no ring is found -- the tightest part of an arch, where a plane
    stops cutting cross-sections at all -- the profile is interpolated from its
    neighbours rather than guessed.
    """
    V = np.asarray(mesh.vertices, float)
    ss = np.arange(0.0, float(route.length) + 1e-9, step)
    r = np.full(len(ss), np.nan)
    # A SECTION THROUGH AN OPEN END IS NOT A SECTION. At station 0 the plane passes
    # through the inlet itself and the component it finds wraps the root's flare:
    # measured, 103 mm on KK4_2 and 132 on KR47_1 where the vessel is near 20. Any
    # ring touching a boundary loop -- the aorta's own two ends, or a branch's cut --
    # is dropped and interpolated over like a station with no ring at all.
    from branch_cut import boundary_loops
    on_edge = np.zeros(len(V), bool)
    for L_ in boundary_loops(mesh):
        on_edge[np.asarray(L_, int)] = True
    for i, x in enumerate(ss):
        name, local, _f = route.locate(float(x))
        c, t = route.tree.segments[name].at_arclength(float(local))
        ring = _cross_section_ring(mesh, np.asarray(c, float), np.asarray(t, float))
        if ring is None or len(ring) < 8:
            continue
        d = np.linalg.norm(V[ring] - np.asarray(c, float), axis=1)
        # A CROSS SECTION IS ROUND. Where the plane grazes the wall instead of cutting
        # across it, the component it finds is a strip lying IN the plane but running
        # along the vessel, and its "radius" runs from nothing to several calibres. A
        # true section, even an oval one, holds its radius to within about a fifth.
        # Rejecting on that is a shape test, so it does not need to know how wide this
        # aorta is -- which the size tests it replaces did, and got wrong on KK4_2
        # where the profile read 103 mm at the root and 71-83 across the arch.
        med = float(np.median(d))
        if med <= 1e-6 or float(np.percentile(d, 90)) / med > 1.6:
            continue
        if on_edge[ring].any():
            continue
        r[i] = med
    ok = np.isfinite(r)
    if ok.sum() < 2:
        return ss, np.full(len(ss), np.nan)
    # AND REJECT THE RINGS THAT ARE NOT CROSS SECTIONS. Where a normal plane grazes
    # the wall instead of cutting across it -- the aortic root's trumpet, the tightest
    # part of an arch -- the component it finds runs ALONG the vessel and its median
    # radius is several times the true one. Measured on KK4_2: 103 mm at the root and
    # 71-83 mm across the arch against a true radius near 20. Those are not radii, and
    # they are worse than a gap, because a gap is interpolated over and these are
    # believed. An aorta does not change calibre by more than a factor of two over its
    # length, so anything outside that band of the case's own median is dropped and
    # interpolated like any other station with no ring.
    # BACKSTOP AT THE ENDS. The shape test catches a grazing plane, but at station 0
    # the plane can cut a perfectly round contour around the root's flare and be
    # believed: 103 mm on KK4_2, 132 on KR47_1. Nothing about an aorta is twice its own
    # median calibre, so that is dropped too, on a second pass once the median is known
    # from the sections that survived the shape test.
    med = float(np.median(r[ok]))
    ok &= r < 2.0 * med
    if ok.sum() < 2:
        return ss, np.full(len(ss), med)
    r = np.interp(np.arange(len(ss)), np.nonzero(ok)[0], r[ok])
    return ss, r


def off_the_aorta(mesh, route, k=1.35, step=1.0):
    """True where a vertex is NOT part of the aorta's own tube.

    The user's definition of a branch: the part that goes off and is not part of the
    greater curve of the aorta. A vertex about one aortic radius from the aortic
    centerline is the aorta; further out, it is on something leaving it.
    """
    from scipy.spatial import cKDTree
    V = np.asarray(mesh.vertices, float)
    ss, r = aortic_radius_profile(mesh, route, step)
    if not np.all(np.isfinite(r)):
        return np.zeros(len(V), bool)
    P = np.array([route.point(float(x)) for x in ss])
    d, i = cKDTree(P).query(V)
    return d > k * r[i]


# ---------------------------------------------------------------------------
# THE ZONE LABELLING, MOVED HERE FROM tools/zone_stub_bench.py.
#
# It lived in tools/ because that is where it was developed, and tools/ is a bench
# directory: it imports anat_bench, it hard-codes a cohort path, and nothing that
# ships is allowed to depend on it. That made the zone map unreachable from the two
# front ends -- webapp/cut_app.py puts only the repo root and webapp/ on sys.path, so
# the browser could not import the pipeline at all, and branch_gui.py never tried.
#
# So the labelling lives in the library layer, beside the ostium finder it is built
# from, and tools/zone_stub_bench.py imports it back. Every bench script and test that
# already says `from zone_stub_bench import label_case` keeps working unchanged.
# ---------------------------------------------------------------------------

TILING = ("Ishimaru zone 0", "Ishimaru zone 1", "Ishimaru zone 2",
          "Ishimaru zone 3", "Ishimaru zone 4")
WALL_K = 0.95
# Place each boundary on its landmark's own CROSS SECTION rather than by interpolating
# the station->u curve. MUST MATCH tools/zone_map.py's flag of the same name: this file
# exists to measure what ships, and a reference implementation that ships something else
# measures nothing.
CROSS_SECTION_PLACEMENT = False
CONTOUR_TOL_MM = 8.0
# Place each boundary at the distal lip of the mouth of the vessel it names.
# OFF until reviewed: it is right by definition and it moves cases in BOTH directions,
# because the landmark station is sometimes proximal of the real mouth and sometimes
# distal of it -- KK32_0's subclavian lip is 8.7 mm proximal of its station, KK20_2's is
# 4.3 mm distal. tools/zone_map.py does not implement it at all, so what ships is
# unchanged while this is decided.
RIM_PLACEMENT = False
# Correct the subclavian's landmark to its own mouth where the two disagree.
LANDMARK_LAG_FIX = False
LANDMARK_LAG_MIN_MM = 2.0


def wall_mask(mesh, route, k=WALL_K):
    """Vertices certainly on the aortic wall, from the aorta's own radius profile."""
    ss, rr = aortic_radius_profile(mesh, route)
    P = np.array([route.point(float(t)) for t in ss])
    d, i = cKDTree(P).query(np.asarray(mesh.vertices, float))
    return d < k * rr[i]


def label_case(s, freeze_stubs=False, cap_crease=True, cross_section=None,
               zone3_nudge=0.0):
    """The labelling tools/zone_map.py produces, stopping short of stub ownership.

    freeze_stubs passes the crease regions to flatten_boundary as its *frozen* set --
    the argument that function has always had and that the shipped path could never
    supply, because the only branch mask available was the one that fails at the arch.
    Without it the collar around a boundary ring reaches across a nearby mouth and
    re-decides those vertices by which side of the fitted plane they fall on.

    OFF, because it was measured and it is a wash: raggedness median 2.326 -> 2.294,
    flatness median 2.74 -> 2.82 mm, islands 0 either way. It trades the two against
    each other in the third decimal place, so the labelling stays as reviewed and this
    stays here as a settled question rather than an untried idea.
    """
    mesh = s.mesh
    route = BC.aorta_route(s.tree)
    L = float(route.length)
    station, _own, onr = station_field(mesh, s.tree, route)
    is_branch = ~onr
    adj = surface_graph(mesh)
    # CAP ON THE CREASE REGIONS, not on the projection's branch call. The capped solve
    # strips the branches and closes their mouths so the holes cannot bend the level
    # sets; feeding it a mask that is wrong at the arch means cutting wall out instead.
    # Measured on KK40_3 it was removing 2377 vertices where the stubs are 511, and the
    # field paid for it: |grad u| p99 0.01824 against 0.00516, hot spots where the
    # level sets bunch and a boundary runs lengthwise down the vessel. On KY68_1 that
    # was a tongue of zone 2 down the whole ascending aorta.
    reg = stub_regions(mesh, route, adj=adj)
    cap = (reg >= 0) if (cap_crease and reg.max() >= 0) else is_branch
    u = harmonic_coordinate_capped(mesh, route, cap)
    if u is None:
        u = harmonic_coordinate(mesh, route)
    ss, uu = calibrate_u_to_station(mesh, route, u)
    zones = sorted((p for p in s.presets() if p["name"] in TILING),
                   key=lambda p: p["s0"])
    # THE SUBCLAVIAN'S LANDMARK, CHECKED AGAINST ITS OWN MOUTH. Median lag over this
    # cohort is 0.1 mm, so on most arches nothing moves; where the branch centerline
    # runs alongside the aorta before diverging the junction lands well past the
    # ostium -- 14 mm on KK32_0 -- and zone 2 then runs past the vessel it ends at.
    lag = 0.0
    if LANDMARK_LAG_FIX and reg.max() >= 0:
        oz0 = stub_ostium_zone(s.tree, route, mesh, reg, zones)
        mouth = landmark_lag(mesh, u, (ss, uu), reg, oz0, zone=2)
        lsa = [p["s1"] for p in zones
               if int(re.search(r"(\d+)\s*$", p["name"]).group(1)) == 2]
        if mouth is not None and lsa and lsa[0] is not None:
            lag = float(lsa[0]) - mouth
            lag = lag if lag > LANDMARK_LAG_MIN_MM else 0.0
    if zone3_nudge or lag:
        zones = shift_zone_bounds(zones, float(zone3_nudge) + lag)
    # BY NUMBER, NOT BY POSITION. A bovine arch offers zones 0, 2, 3, 4 -- the left
    # carotid shares the brachiocephalic's origin so zone 1 has zero length -- and
    # indexing by position calls zone 2 "zone 1".
    znum = [int(re.search(r"(\d+)\s*$", p["name"]).group(1)) for p in zones]
    # A BOUNDARY IS A VESSEL, SO PUT IT ON THAT VESSEL'S MOUTH. By the half-open rule a
    # zone owns the vessel at its distal end and the next begins immediately after it,
    # so the boundary closing zone k belongs at the far lip of the mouth of the stub
    # zone k owns. That lip is read off the wall from the crease, with no station->u
    # curve in between -- and the curve is what put it wrong. Measured as millimetres of
    # wall between the subclavian's mouth and the start of zone 3, the cohort was
    # carrying a median 22 mm of zone 2 past the vessel, which is most of zone 3's 20 mm.
    #
    # An earlier attempt used the rim's MEDIAN, on the stub side, which sits in the
    # middle of the mouth and moved every boundary further out instead.
    xsec = CROSS_SECTION_PLACEMENT if cross_section is None else cross_section
    rim_ = (stub_rim_coord(mesh, u, reg)
            if (RIM_PLACEMENT and reg.max() >= 0) else {})
    oz_ = stub_ostium_zone(s.tree, route, mesh, reg, zones) if rim_ else {}
    by_zone = {}
    for _j, _z in oz_.items():
        if _j in rim_:
            by_zone.setdefault(_z, []).append(rim_[_j])

    def _ub(st, zone=None):
        if zone is not None and zone in by_zone:
            return float(max(by_zone[zone]))
        base = float(np.interp(st, ss, uu))
        if not xsec:
            return base
        seg = s.tree.segments[route.locate(float(st))[0]]
        v = boundary_u_at(mesh, route, u, st, seg, prior=base)
        if v is None:
            return base
        back = float(np.interp(v, uu, ss))
        return v if abs(back - float(st)) <= CONTOUR_TOL_MM else base

    lab = -np.ones(len(u), int)
    for zi, p in enumerate(zones):
        lo = -1.0 if zi == 0 else _ub(p["s0"], znum[zi - 1])
        hi = (2.0 if (p["s1"] is None or zi == len(zones) - 1)
              else _ub(p["s1"], znum[zi]))
        m = (~is_branch) & (u >= lo) & (u < hi)
        inb = branch_inside(np.where(is_branch, station, -1e9),
                               p["s0"], p["s1"], L)
        lab[m | (is_branch & inb)] = znum[zi]
    lab, _ = enforce_connected_labels(mesh, lab)
    frozen = (reg >= 0) if (freeze_stubs and reg.max() >= 0) else np.zeros(len(u), bool)
    # The aorta's own two ends are not up for re-decision by a boundary a zone away.
    frozen = frozen | vessel_end_collar(mesh, route, adj)
    # HALF OF THE NARROWER NEIGHBOUR. A boundary may tidy its own collar; it may not
    # reach across a zone. Zone 0 on KK40_3 is 30.8 mm long and the collar was 1.2 ring
    # radii, so it took in the aortic root.
    ub = [(-1.0 if zi == 0 else float(np.interp(p["s0"], ss, uu)),
           2.0 if (p["s1"] is None or zi == len(zones) - 1)
           else float(np.interp(p["s1"], ss, uu)))
          for zi, p in enumerate(zones)]
    flat = []
    for zi in range(len(zones) - 1):
        za, zb = znum[zi], znum[zi + 1]
        # ...and flattened ONTO that cross section, not onto the plane of wherever the
        # ring currently sits. Fitting the current ring makes it flat while locking in
        # its displacement. zone_map has always done this; this file had not, so it was
        # measuring a boundary the tool does not produce.
        tgt = None
        bb = zones[zi]["s1"]
        if bb is not None and xsec:
            base = float(np.interp(bb, ss, uu))
            seg = s.tree.segments[route.locate(float(bb))[0]]
            v = boundary_u_at(mesh, route, u, bb, seg, prior=base)
            if (v is not None
                    and abs(float(np.interp(v, uu, ss)) - float(bb)) <= CONTOUR_TOL_MM):
                tgt = cross_section_plane(mesh, route, s.tree, bb)
        # Protecting the two ends is enough. Confining the collar along the vessel as
        # well was measured: it fixes the same root and costs KR27_2's 1/2 boundary
        # 2.74 -> 11.29 mm off plane, KK32_0's 0/1 3.17 -> 9.58.
        lab, dev = flatten_boundary(mesh, lab, za, zb, adj, frozen=frozen,
                                       target=tgt)
        if dev is not None:
            flat.append(float(dev))
    lab, _ = enforce_connected_labels(mesh, lab)
    return mesh, route, zones, lab, is_branch, adj, flat


def own_by_stub(mesh, lab, regions, is_branch, adj=None, region_zone=None):
    """What ships, called by name. A checker with its own copy of the rule can only
    ever agree with itself, and this file exists to disagree when the rule is wrong."""
    return own_stubs(mesh, lab, regions, is_branch, adj=adj,
                        region_zone=region_zone)


def ostium_zones(s, mesh, route, regions):
    """Each stub's zone straight from its own ostium on the centerline."""
    zones = sorted((p for p in s.presets() if p["name"] in TILING),
                   key=lambda p: p["s0"])
    return stub_ostium_zone(s.tree, route, mesh, regions, zones)


# ZONE 3 IS THE 20 MM BAND JUST PAST THE SUBCLAVIAN, CUT BY TWO PLANES. Its proximal
# edge is a centerline-normal plane at the distal lip of that vessel's mouth; its distal
# edge is the same plane 20 mm further on. Every earlier attempt moved a LEVEL SET of the
# harmonic coordinate instead, and a level set follows the wall -- on a bend it is
# oblique, so no amount of moving it makes it the transverse cut the clinical figure
# draws. That is why the boundary looked wrong however it was placed.
PLANE_BOUNDARY = True
ZONE3_LENGTH_MM = 20.0


def final_labels(s, **kw):
    """The labelling that ships, start to finish, in one call.

    Three separate times a checker in this directory was found to be scoring a program
    the tool does not run, because each rebuilt the last few steps for itself. There is
    now one path: everything -- the map, the gate, the renders -- comes through here.

    Returns (mesh, route, zones, lab, is_branch, adj, flatness, regions, ostium_zones).
    """
    mesh, route, zones, lab, isb, adj, flat = label_case(s, **kw)
    reg = stub_regions(mesh, route, adj=adj)
    oz = ostium_zones(s, mesh, route, reg) if reg.max() >= 0 else {}
    if reg.max() >= 0:
        lab = own_stubs(mesh, lab, reg, isb, adj=adj, region_zone=oz)
    if PLANE_BOUNDARY and reg.max() >= 0:
        lab = plane_boundaries(s, mesh, route, lab, reg, oz, adj)
    return mesh, route, zones, lab, isb, adj, flat, reg, oz


def plane_boundaries(s, mesh, route, lab, reg, oz, adj):
    """Cut zone 3 out with two centerline-normal planes, off the subclavian's mouth."""
    lsa = [j for j, z in oz.items() if z == 2]
    zs = sorted(set(lab[lab >= 0].tolist()))
    if not lsa or 2 not in zs or 3 not in zs:
        return lab
    # NO COORDINATE. The anchor is the first plane that clears the mouth, which is a
    # question about the rim's geometry alone -- so the capped harmonic solve is no
    # longer on this path. That matters beyond tidiness: it used to return None on
    # KK5_1 and KY183_1 (a singular Laplacian from stray disconnected trunk pieces) and
    # plane_boundaries then returned the labelling untouched, so those two cases kept
    # the old level-set boundary and no change to the rule could reach them. One of
    # them was among the cases the user complained about.
    got = mouth_clearing_plane(mesh, route, reg, lsa,
                                  clearance_mm=MOUTH_CLEARANCE_MM)
    if got is None:
        return lab
    origin, normal, st = got
    e = np.asarray(mesh.edges_unique, int)
    sel = np.isin(reg, lsa)
    cut = sel[e[:, 0]] != sel[e[:, 1]]
    rim = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
    frozen = reg >= 0
    # The plane is frozen off the stubs AND their first ring of neighbouring wall. A
    # face with one corner inside a mouth and one corner beyond the plane renders as
    # zone 3 sitting on the branch even though no stub vertex was ever repainted --
    # 677 such faces over the cohort, on 16 of 27 cases. Freezing the 1-ring removes
    # 99% of them. The REPAIR pass below keeps the narrow mask: widening that one
    # strands vertices it can no longer absorb and leaves new islands.
    ring1 = np.zeros(len(frozen), bool)
    ring1[e[frozen[e[:, 0]] | frozen[e[:, 1]]].ravel()] = True
    lab, _n = split_by_plane(mesh, lab, 2, 3, origin, normal, adj=adj,
                               frozen=frozen | ring1, near=rim)
    # ...and the far edge, the same plane 20 mm on, so zone 3 is the band the
    # definition says it is rather than whatever the coordinate left behind.
    if 4 in zs:
        far = min(st + ZONE3_LENGTH_MM, float(route.length) - 1e-3)
        name, local, _f = route.locate(far)
        c2, t2 = s.tree.segments[name].at_arclength(float(local))
        n2 = np.asarray(t2, float)
        n2 = n2 / max(np.linalg.norm(n2), 1e-12)
        ahead = route.point(min(far + 1.0, float(route.length)))
        if float(np.dot(np.asarray(ahead, float) - np.asarray(c2, float), n2)) < 0:
            n2 = -n2
        lab, _m = split_by_plane(mesh, lab, 3, 4, np.asarray(c2, float), n2,
                                    adj=adj, frozen=frozen | ring1, near=rim)
    lab, _mv = enforce_connected_labels(mesh, lab, frozen=frozen)
    return lab


# ---------------------------------------------------------------------------
# One palette, shared. tools/zone_map.py drew the cohort renders, zone_viewer.py the
# desktop window and webapp/cut_app.py the browser, and three copies of a colour list
# drift -- a zone that is gold in the paper and orange in the tool is a bug report
# waiting to happen. The numbering matters as much as the colours: on a bovine arch
# the left carotid arises from the brachiocephalic, zone 1 has zero length and is
# never offered, and the case runs 0, 2, 3, 4. Colour by the zone's NUMBER, never by
# its position in the list.
# ---------------------------------------------------------------------------
ZONE_COLOURS = ["#1d4ed8", "#0e9488", "#b8860b", "#c2410c", "#7d3ac1"]
ZONE_UNCLAIMED = "#d9d4cf"
ZONE_NAMES = {
    0: "ascending, to the brachiocephalic",
    1: "brachiocephalic to left carotid",
    2: "left carotid to left subclavian",
    3: "subclavian, and 20 mm on",
    4: "beyond zone 3",
}


def hex_to_rgb(h):
    h = str(h).lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float)


def zone_face_colours(mesh, label, unclaimed=ZONE_UNCLAIMED):
    """Per-FACE uint8 RGB for a per-vertex zone labelling.

    Per face rather than per vertex because that is what the viewers take, and because
    a face whose three corners disagree should read as the earlier zone rather than be
    smeared between two: the boundary is then a clean edge on the mesh instead of a
    gradient a millimetre wide, which is what it is supposed to be.
    """
    lab = np.asarray(label, int)
    F = np.asarray(mesh.faces, int)
    zf = lab[F]
    zf = np.where(zf < 0, np.iinfo(np.int32).max, zf).min(axis=1)
    out = np.tile(hex_to_rgb(unclaimed), (len(F), 1))
    for z in range(len(ZONE_COLOURS)):
        out[zf == z] = hex_to_rgb(ZONE_COLOURS[z])
    return out.astype(np.uint8)
