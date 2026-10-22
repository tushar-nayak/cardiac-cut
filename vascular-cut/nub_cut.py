"""Nub-in cuts: every aortic branch left as a short stub proud of the aortic wall.

The preset for the branches paper. Instead of cutting each branch somewhere along
its useful length, cut it just after it clears the aorta, so what survives is a
"nub" -- a stub of vessel standing proud of the aortic surface, the same length
rule applied to every branch of every case.

THE RULE

    n = s_wall + k * (D_aorta / 2)

*n* is measured ALONG THE BRANCH CENTERLINE, starting at the junction node, which
sits on the aortic centerline. So the walk begins on the aortic axis and goes
outward up the branch, and *n* is where the cut plane lands.

    s_wall     how far up the branch you must go before you are clear of the aorta
    D_aorta    the aorta's local diameter where that branch leaves it
    k          how much stub you want past the wall, in aortic radii

WHY s_wall IS MEASURED AND NOT ASSUMED

The rule was first written as n = D/2 + k*(D/2): travel one aortic radius to reach
the wall, then k of another for the nub. That is right whenever a branch leaves
perpendicular, and wrong by the obliquity whenever it does not -- and the arch
vessels do not. Measured on KY160_1, walking each branch out from its junction
node and asking at every station whether it is more than one local aortic radius
off the aortic axis:

    vessel      D_aorta/2 predicts     actually clear at
    Branch 1         13.3 mm                18.0 mm
    Branch 3         13.3 mm                10.0 mm
    Branch 2         18.3 mm                21.0 mm
    Line1            22.0 mm                21.0 mm

The estimate is out by up to 4.7 mm. At k = 0.1 the intended nub is 1.3-2.2 mm --
shorter than the error -- so on two of those four the plane landed INSIDE the
aortic lumen and there was no nub to have. Measuring s_wall makes the proud length
mean what it says at every k, and where a take-off really is perpendicular the two
forms agree exactly.

Two other signals were tried and rejected, both on KY160_1:

  the branch's own faces      Fires 3-5 mm early. A branch OWNS its flare -- the
                              funnel where its mouth meets the aorta -- and that
                              flare lies inside the aortic envelope. Right answer
                              to "where does this vessel's wall begin", wrong one
                              for "where am I clear of the aorta".
  perpendicular section size  Confounded by obliquity. Near an ostium the plane
                              slices the aorta at a shallow angle, so the section
                              collapses as the branch turns away rather than as it
                              leaves. It put Branch 1's wall at 9 mm, inside a
                              13.3 mm-radius aorta.
"""
import numpy as np
from scipy.spatial import cKDTree

from branch_cut import (
    FOREIGN_TOLERANCE_MM2,
    MAX_OBLIQUITY_DEG,
    MIN_OFFSET_MM,
    RIM_FLATNESS_MM,
    _loops_by_terminal,
    _junction_graph,
    aorta_route,
    boundary_loops,
    cut_terminal,
    band_width,
    keep_largest_component,
    obliquity_deg,
    preflight,
    radius_at,
    radius_profiles_from_faces,
    segment_radius,
    snap_tolerance,
    suggest_pullback,
    vessel_attribution,
    validate_cut,
    vessel_coords,
)
from branch_tree import Terminal

# How much stub past the aortic wall, in aortic radii.
#
# The figure this started from was 0.1, which on a thoracic aorta is 1.3-2.2 mm of
# nub -- shorter than the error in the term it was added to, and on KY160 it put two
# of four planes inside the aortic lumen.
#
# 0.75 is the cohort's own answer rather than a preference. A branch leaves its aorta
# as a funnel and only becomes a tube some way out, and a plane through a funnel does
# not come back flat, so there is a shortest nub each vessel can have. Over 82 vessels
# in 30 thoracic cases that shortest length is k = 0.50 at the median and k = 0.76 at
# the 90th percentile. Ask for less and a quarter of the cohort is silently lengthened
# to what it can take, which is the one thing this preset exists to prevent: at 0.5,
# 20 of 82 vessels were grown; at 0.75, 8.
NUB_K = 0.75

# The wall walk. 0.25 mm is well under the finest thing being resolved and the walk is
# at most a few hundred steps per branch, so there is no reason to be coarser.
WALL_STEP_MM = 0.25

# A crossing only counts if the branch STAYS clear for this far. Without it a single
# station where the aortic radius profile happens to dip reads as the ostium.
WALL_HOLD_MM = 2.0

# How far up a branch to bother looking. Beyond a few aortic radii the branch has
# either left or it is running alongside the aorta, and neither wants a longer walk.
WALL_SEARCH_RADII = 3.0
WALL_SEARCH_PAD_MM = 15.0

# Only aortic centerline within this much arc length of the ostium is a candidate for
# "the aorta near this branch". Stops a branch that curves back toward a different part
# of the aorta from being measured against the wrong station.
AORTA_WINDOW_PAD_MM = 40.0

# Keep the plane this far short of whatever ends the vessel -- its tip, or the split
# at the far end of a trunk. A cut sitting exactly on a bifurcation has no vessel to
# leave behind on one side.
NUB_MARGIN_MM = 1.5

# Below this a "nub" is not one, whatever the arithmetic says. Only ever reached by a
# clamp against a short trunk, and it flags rather than rejects.
MIN_USEFUL_PROUD_MM = 2.0

# A first-order vessel feeding this many great vessels is a SHARED ORIGIN, not a
# brachiocephalic, and its daughters carry the stubs instead of it.
#
# A normal arch has three vessels leaving the aorta and the brachiocephalic feeds two,
# so every stub stands for at most two great vessels. A bovine arch has two leaving the
# aorta and the shared stem feeds three -- so a stub on the stem stands for three, and
# the left common carotid disappears inside it. Eight of the thirty thoracic cases are
# bovine by this test, and before this rule they were being handled two different ways:
# six had a stem too short to cut and were handed down anyway, while KY160_1 and
# KY215_1 had a stem long enough and kept it, coming out with two branch stubs where
# every other case has three. Same anatomy, different geometry, decided by whether the
# stem happened to clear the aorta before it split.
#
# Deciding it on what the vessel FEEDS rather than on whether it happens to be
# cuttable makes every stub in the cohort stand for at most two great vessels.
SHARED_ORIGIN_VESSELS = 3

# How much further out to push a nub whose rim came out torn once the whole model was
# assembled. Two millimetres clears it on both cohort cases that needed it.
RETRY_STEP_MM = 2.0

# Searching for the shortest nub that cuts cleanly. Coarse first because most vessels
# need no growth at all, then finer, because a millimetre of vessel is worth keeping.
# The limit is generous: the most any vessel in the thoracic cohort needed was 18.8 mm.
GROW_STEP_MM = 2.0
GROW_FINE_MM = 0.5
GROW_LIMIT_MM = 40.0

# How far the cut reaches around its plane, as multiples of the width the tool picks
# automatically. That width is scaled from the segment's measured radius, and where
# that radius is wrong the band is too narrow to close a ring -- the cut then leaves
# part of the original surface in the rim and the opening reads as torn at EVERY
# offset. KK40_3's blocking vessel measures 7.43 mm in radius against a section 33 to
# 42 mm across, so its 18.6 mm band could never span it; at 30 mm the same vessel cuts
# at 0.000 mm off plane. Widening the plane is what a person does with the Plane size
# slider, so the preset tries it before giving up, and records what it used.
BAND_TRIES = (1.0, 1.6, 2.4, 3.2)

# Tilting the plane off centerline-perpendicular, tried only when NOTHING flat works.
# A vessel with almost no tube past the aortic wall is still in its flare where the
# plane lands, and a perpendicular slice through a flare cannot come back flat.
# KY244_1's bovine brachiocephalic has 2.5 mm of runway: flat, its best rim over 56
# cuts spanning every offset, both modes and bands from 9.6 to 63.7 mm is 1.798 mm.
# Tilted, it cuts at 0.000 mm, and not on a knife edge -- 165 combinations work,
# clustered around tilt B -20 to -30 with tilt A about -5.
TILT_TRIES_DEG = (0, -10, 10, -20, 20, -30, 30)

# The aorta's own two ends: how far in from each end to cut, and on what plane. Not
# nubs -- these are the model's inlet and outlet, squared off horizontally so the
# ascending and descending ends look the same across the cohort.
TRUNK_OFFSET_MM = 8.0
TRUNK_MODE = "z"

# Offsets are placed on the desktop tool's own slider grid. Its sliders move in tenths
# of a millimetre, so an offset off that grid cannot be represented there: opening a
# case in the window to adjust it would shift every untouched cut by up to 0.05 mm and
# quietly hand back a different mesh -- 66 faces different on KY160, from a plane
# moved 0.04 mm. A tenth of a millimetre is far below anything these meshes resolve,
# so nothing is lost by rounding, and the batch result and the window's result become
# the same geometry rather than two nearly identical ones.
OFFSET_GRID_MM = 0.1


def _grid(mm):
    """An offset the desktop tool's sliders can hold exactly."""
    return float(round(float(mm) / OFFSET_GRID_MM) * OFFSET_GRID_MM)


# A segment shorter than this many of its neighbour's radii cannot measure its own
# radius: segment_radius samples the middle of the segment, and within a couple of
# radii of the mesh's end the nearest surface is the END CAP rather than the wall.
SPLIT_STUB_RADII = 2.0


# A connected piece smaller than this fraction of the model is debris, not anatomy.
# KK40_3 carries a 4-face tetrahedron 0.2 mm across, floating clear of a 45,658-face
# aorta. It is dropped at the end anyway by keep_largest_component, but until then it
# is part of the surface every probe measures against.
DEBRIS_FRACTION = 0.005


def drop_debris(mesh, verbose=False):
    """The model without any stray specks. Returns *mesh* itself when there are none."""
    try:
        parts = mesh.split(only_watertight=False)
    except Exception:                                              # noqa: BLE001
        return mesh
    if len(parts) < 2:
        return mesh
    keep = [p for p in parts if len(p.faces) >= DEBRIS_FRACTION * len(mesh.faces)]
    if len(keep) == len(parts) or not keep:
        return mesh
    dropped = len(parts) - len(keep)
    out = keep[0] if len(keep) == 1 else __import__("trimesh").util.concatenate(keep)
    if verbose:
        print(f"  dropped {dropped} stray piece(s) of mesh debris")
    return out


class NubError(Exception):
    """A nub plan that cannot be built at all, as opposed to one with bad vessels."""


def continuation_segments(tree, jg=None):
    """seg -> the other pieces of the SAME vessel, joined at non-junctions.

    A node where two segments meet is a break in the exported line, not a branch
    point, so the pieces either side are one vessel. Anything asking "did this cut
    take surface off a vessel it was not cutting" has to know that, or it answers no
    when the cut was perfect: DJ18_2's descending aorta cuts with a rim 0.000 mm off
    plane and is rejected for taking 10.43 mm2 of another vessel -- 8 faces, all of
    them the far half of the aorta being cut.
    """
    jg = jg or _junction_graph(tree)
    members = {}
    for name in tree.segments:
        for node in jg[name]:
            if node is not None:
                members.setdefault(node, []).append(name)

    group = {n: {n} for n in tree.segments}          # union-find, small enough to be naive
    for segs in members.values():
        if len(segs) != 2:
            continue
        a, b = segs
        merged = group[a] | group[b]
        for n in merged:
            group[n] = merged
    return {n: (g - {n}) for n, g in group.items()}


def heal_split_segments(tree, mesh=None, verbose=False):
    """Repair centerlines that arrive broken into pieces at a NON-junction.

    A node where exactly TWO segments meet is not a bifurcation -- it is a break in
    the exported line, and the two pieces are one vessel. Seven of the thirty
    thoracic cases have one, and where the piece left over is a short stub at the end
    of the model it poisons everything downstream of its radius:

    DJ18_2's descending aorta is split into a 22.2 mm stub and a 224.9 mm remainder.
    segment_radius measures the stub by sampling its middle, but the middle of a
    22 mm stub sits within a few millimetres of the mesh's end cap, so the nearest
    surface it finds is that cap and not the wall. It came out 7.84 mm on a vessel
    the cut plane measures 47-53 mm across. Band width is scaled from that radius, so
    the cut was told to reach 19.6 mm around a 50 mm vessel: it cannot close a ring,
    and the opening came out 15 mm out of plane. Every offset failed, in both plane
    modes -- which reads as "this vessel cannot be cut" when the truth is that it was
    measured wrong.

    Giving the stub the radius of the vessel it continues into fixes it outright:
    the same end then cuts with a rim 0.000 mm off plane.

    Radius only. The segments are not merged and no geometry is touched, because
    merging them would change what every other part of the tool sees.
    """
    jg = _junction_graph(tree)
    members = {}
    for name in tree.segments:
        for node in jg[name]:
            if node is not None:
                members.setdefault(node, []).append(name)

    healed = []
    for node, segs in members.items():
        if len(segs) != 2:
            continue                          # a real junction, or a free end
        a, b = (tree.segments[s] for s in segs)
        if mesh is not None:
            for s in (a, b):
                if getattr(s, "_radius_mm", None) is None:
                    s._radius_mm = segment_radius(mesh, s)
        ra = float(getattr(a, "_radius_mm", None) or 0.0)
        rb = float(getattr(b, "_radius_mm", None) or 0.0)
        if ra <= 0 or rb <= 0:
            continue
        # The longer piece is the one that could measure itself; the shorter inherits
        # only if it was too short to have measured its own radius honestly.
        short, long_ = (a, b) if a.length <= b.length else (b, a)
        r_long = float(long_._radius_mm)
        if short.length >= SPLIT_STUB_RADII * r_long:
            continue
        if abs(float(short._radius_mm) - r_long) < 0.5:
            continue
        healed.append((short.name, float(short._radius_mm), r_long, long_.name))
        short._radius_mm = r_long

    if healed:
        # Ownership is decided by distance NORMALISED BY EACH SEGMENT'S RADIUS, and
        # its cache is keyed by segment names and lengths -- not radii. So changing a
        # radius leaves a cache that silently contradicts it, and the difference only
        # shows up wherever the cache happens to have been filled first.
        #
        # That is exactly how this went wrong. The batch computes ownership after the
        # repair and gets DJ18_2's aortic root rim at 0.478 mm, inside the 0.5 mm
        # standard. The window cuts once as it loads, BEFORE the repair, so it had
        # already cached ownership from the 6.28 mm radius; the repair then changed
        # the radius and not the cache, and the same case came up in the window at
        # 0.7 mm and refused to export -- the batch and the window disagreeing about
        # a cut they had both been given identically.
        for attr in ("_attr_full",):
            if hasattr(tree, attr):
                delattr(tree, attr)

    if verbose:
        for nm, was, now, src in healed:
            print(f"  split-line repair: {nm} radius {was:.2f} -> {now:.2f} mm "
                  f"(inherited from {src}, which it continues into)")
    return healed


def aorta_reference(mesh, tree, prof=None, route=None, jg=None):
    """The aortic centerline as points, each carrying the local lumen radius.

    One array rather than a per-segment lookup, because the walk up a branch has to
    ask "how far am I from the aorta, and how wide is the aorta there" without
    knowing in advance which route segment the answer comes from -- a branch off the
    arch is nearest the arch at its ostium and can be nearest the ascending aorta a
    centimetre later.
    """
    jg = jg or _junction_graph(tree)
    route = route or aorta_route(tree, jg)
    if route is None:
        raise NubError("no aortic route: the tree has no trunk to hang branches off")

    if prof is None:
        owner, s_val, radial = vessel_coords(mesh, tree, return_radial=True)
        prof = radius_profiles_from_faces(tree, owner, s_val, radial)

    pts, radii, station = [], [], []
    for i, name in enumerate(route.chain):
        seg = tree.segments[name]
        r = radius_at(prof, name, seg.s)
        s_route = route.starts[i] + (seg.s if route.forward[i]
                                     else (seg.length - seg.s))
        pts.append(seg.pts)
        radii.append(np.asarray(r, float))
        station.append(np.asarray(s_route, float))

    pts = np.vstack(pts)
    return {"kd": cKDTree(pts), "pts": pts,
            "radii": np.concatenate(radii),
            "station": np.concatenate(station),
            "route": route, "prof": prof, "jg": jg}


def _windowed(ref, near_point):
    """A KD tree over just the stretch of aorta around *near_point*."""
    _, i = ref["kd"].query(np.asarray(near_point, float))
    s0 = float(ref["station"][i])
    span = WALL_SEARCH_RADII * float(ref["radii"][i]) + AORTA_WINDOW_PAD_MM
    keep = np.flatnonzero(np.abs(ref["station"] - s0) <= span)
    return cKDTree(ref["pts"][keep]), ref["radii"][keep]


def wall_crossing(seg, s_origin, direction, ref):
    """Where the branch centerline gets clear of the aorta, walking out from the node.

    Returns a dict with *s_wall_mm* (arc length from the junction node), the local
    aortic radius there, and whether a crossing was found at all. Distances are to
    the aortic centerline, compared against the aorta's own measured radius at the
    nearest station, so an oblique take-off is handled by construction.
    """
    p0, _ = seg.at_arclength(s_origin)
    kd, radii = _windowed(ref, p0)

    d0, i0 = kd.query(p0)
    r_here = float(radii[i0])
    reach = min(WALL_SEARCH_RADII * r_here + WALL_SEARCH_PAD_MM, seg.length)

    ss = np.arange(0.0, reach + 1e-9, WALL_STEP_MM)
    pts = np.array([seg.at_arclength(s_origin + direction * s)[0] for s in ss])
    d, idx = kd.query(pts)
    r = radii[idx]
    clear = d > r

    # The first station clear that STAYS clear for WALL_HOLD_MM. A run test rather
    # than a first crossing: the aortic radius is a percentile over binned wall faces,
    # so it wobbles, and a lone station over the line is noise rather than an ostium.
    hold = max(int(round(WALL_HOLD_MM / WALL_STEP_MM)), 1)
    run = 0
    hit = -1
    for j in range(len(clear) - 1, -1, -1):
        run = run + 1 if clear[j] else 0
        if run >= hold:
            hit = j
    if hit < 0:
        return {"found": False, "s_wall_mm": float("nan"),
                "r_aorta_mm": r_here, "d_aorta_mm": float(d0),
                "reach_mm": float(reach)}

    # Sub-step placement: the crossing is where d - r changes sign, not at whichever
    # sample happened to land past it.
    s_wall = float(ss[hit])
    if hit > 0:
        g0, g1 = float(d[hit - 1] - r[hit - 1]), float(d[hit] - r[hit])
        if g1 != g0:
            s_wall = float(ss[hit - 1] + WALL_STEP_MM * (0.0 - g0) / (g1 - g0))
        s_wall = float(np.clip(s_wall, ss[hit - 1], ss[hit]))

    return {"found": True, "s_wall_mm": s_wall,
            "r_aorta_mm": float(r[hit]), "d_aorta_mm": float(d[hit]),
            "reach_mm": float(reach)}


def first_order_vessels(tree, jg=None, route=None, mesh=None):
    """Every vessel leaving the aorta directly, whether or not it ends in an opening.

    A brachiocephalic trunk carries no terminal of its own -- it splits -- so a list
    built from terminals misses exactly the vessel the preset was asked to include.
    Interior segments are collected here on the same footing, and the terminals they
    swallow are recorded so the caller can drop them from the cut list.
    """
    jg = jg or _junction_graph(tree)
    route = route or aorta_route(tree, jg)
    if route is None:
        return []
    chain = set(route.chain)

    # every junction node that lies on the route
    on_route = set()
    for name in route.chain:
        on_route.update(n for n in jg[name] if n is not None)

    out = []
    for name, seg in tree.segments.items():
        if name in chain:
            continue
        ns, ne = jg[name]
        if ns in on_route:
            node, s_origin, direction, at_end = ns, 0.0, +1.0, True
        elif ne in on_route:
            node, s_origin, direction, at_end = ne, seg.length, -1.0, False
        else:
            continue                      # second order or further; its parent is cut

        if mesh is not None and getattr(seg, "_radius_mm", None) is None:
            seg._radius_mm = segment_radius(mesh, seg)

        own = [t for t in tree.terminals if t.segment.name == name]
        interior = (ns is not None and ne is not None)
        sub, swallowed = set(), []
        if interior:
            from branch_cut import _subtree_segments
            far = ne if at_end else ns
            sub = _subtree_segments(jg, far, name)
            swallowed = [t for t in tree.terminals if t.segment.name in sub]

        out.append({"segment": name, "node": node, "s_origin": s_origin,
                    "direction": direction, "at_end": at_end,
                    "interior": interior, "subtree": sub,
                    "terminal": own[0] if own and not interior else None,
                    "swallowed": swallowed, "depth": 0, "parent_vessel": None})
    out.sort(key=lambda c: c["segment"])
    return out


def daughters_of(tree, cand, jg=None, mesh=None):
    """The vessels a trunk splits into, as nub candidates in their own right.

    Wanted when a trunk cannot hold a nub at all. Four of the thoracic cohort's
    bovine arches and two other short trunks SPLIT BEFORE THEY LEAVE THE AORTA --
    7 to 25 mm of vessel inside an aorta 14 to 28 mm in radius, so there is no
    station on them that is outside it and no nub to cut. Their daughters are 34 to
    185 mm long and clear it easily, so the great vessels are still represented;
    they are just represented one generation down, which is recorded.
    """
    if not cand["interior"]:
        return []
    jg = jg or _junction_graph(tree)
    far = None
    ns, ne = jg[cand["segment"]]
    far = ne if cand["at_end"] else ns
    out = []
    for name in cand["subtree"]:
        seg = tree.segments[name]
        ks, ke = jg[name]
        if ks == far:
            s_origin, direction, at_end = 0.0, +1.0, True
        elif ke == far:
            s_origin, direction, at_end = seg.length, -1.0, False
        else:
            continue                       # a granddaughter, reached via its own parent
        if mesh is not None and getattr(seg, "_radius_mm", None) is None:
            seg._radius_mm = segment_radius(mesh, seg)
        interior = ks is not None and ke is not None
        own = [t for t in tree.terminals if t.segment.name == name]
        sub, swallowed = set(), []
        if interior:
            from branch_cut import _subtree_segments
            sub = _subtree_segments(jg, ke if at_end else ks, name)
            swallowed = [t for t in tree.terminals if t.segment.name in sub]
        out.append({"segment": name, "node": far, "s_origin": s_origin,
                    "direction": direction, "at_end": at_end,
                    "interior": interior, "subtree": sub,
                    "terminal": own[0] if own and not interior else None,
                    "swallowed": swallowed, "depth": cand["depth"] + 1,
                    "parent_vessel": cand["segment"]})
    out.sort(key=lambda c: c["segment"])
    return out


def _evaluate(tree, v, k, ref, names, clamp):
    """One vessel's row: where its nub goes, or why it cannot have one."""
    seg = tree.segments[v["segment"]]
    wall = wall_crossing(seg, v["s_origin"], v["direction"], ref)

    label = (v["terminal"].name if v["terminal"] is not None
             else names.get(v["segment"]) or f"Trunk {v['segment']}")
    if v["depth"]:
        label = f"{label} (daughter of {v['parent_vessel']})"

    row = {"segment": v["segment"], "name": label, "depth": int(v["depth"]),
           "interior": v["interior"], "length_mm": float(seg.length),
           "radius_mm": float(getattr(seg, "_radius_mm", 0.0) or 0.0),
           "k": float(k), **wall}

    # The runway: from the wall out to wherever the vessel stops being cuttable. For
    # a trunk that is its split, and it is usually the binding constraint -- KY160's
    # brachiocephalic clears the aorta at 21.2 mm and splits at 29.3, leaving 6.6 mm
    # of nub against the 10.1 mm that k = 0.5 asks for.
    n_max = float(seg.length - NUB_MARGIN_MM)
    row["n_max_mm"] = n_max

    if not wall["found"]:
        row.update(ok=False, n_mm=float("nan"), offset_mm=float("nan"),
                   k_max=float("nan"), can_descend=True,
                   why=f"splits or ends {wall['reach_mm']:.0f} mm along, still "
                       "inside the aorta -- no station on it is outside")
        return row

    s_wall = wall["s_wall_mm"]
    r_ao = wall["r_aorta_mm"]
    row["k_max"] = float((n_max - s_wall) / r_ao) if r_ao > 0 else float("nan")

    if s_wall >= n_max:
        row.update(ok=False, n_mm=float("nan"), offset_mm=float("nan"),
                   proud_mm=float("nan"), can_descend=True,
                   why=f"clears the aorta at {s_wall:.1f} mm but the vessel ends at "
                       f"{seg.length:.1f} mm -- no room for any nub")
        return row

    n = s_wall + k * r_ao
    if n > n_max:
        if not clamp:
            row.update(ok=False, n_mm=float(n), offset_mm=float("nan"),
                       proud_mm=float(k * r_ao), can_descend=True,
                       why=f"k={k:g} wants {n:.1f} mm, vessel allows {n_max:.1f} mm "
                           f"(k_max={row['k_max']:.2f})")
            return row
        # A vessel with SOME runway keeps its nub, shortened, rather than being
        # handed down to its daughters: a 6.6 mm stub of brachiocephalic is still a
        # brachiocephalic, and descending would drop the vessel from the set
        # entirely and replace it with the two it feeds.
        n = n_max
        row["clamped"] = True
        row["why"] = (f"shortened: k={k:g} wants {s_wall + k * r_ao:.1f} mm, "
                      f"vessel allows {n_max:.1f} mm (k_max={row['k_max']:.2f})")

    row["n_mm"] = float(n)
    row["proud_mm"] = float(n - s_wall)
    row["k_used"] = float((n - s_wall) / r_ao) if r_ao > 0 else float("nan")
    # A clamp can leave a stub too short to be one. KY244's brachiocephalic clears
    # the aorta 1.0 mm before it splits, so it is placeable and still not worth
    # keeping without someone looking at it. Placed, cut, and flagged -- deciding for
    # the reviewer which marginal vessels to drop is not this function's job.
    row["marginal"] = bool(row["proud_mm"] < MIN_USEFUL_PROUD_MM)
    off = _grid(seg.length - n)
    row.update(ok=True, offset_mm=max(off, MIN_OFFSET_MM), can_descend=False)
    row.setdefault("why", "")
    return row


def nub_plan(mesh, tree, k=NUB_K, ref=None, vessels=None, names=None,
             clamp=True, descend=True, max_depth=2):
    """Where every branch's nub cut goes, and whether it fits.

    One row per vessel that ends up carrying a nub. *offset_mm* is what the cutter
    wants -- distance back from the far end -- and equals length - n in every case,
    terminal-bearing or interior, because n is measured from the aortic end and the
    cutter measures from the other one.

    A vessel that cannot hold a nub AT ALL hands the job to its daughters, and the
    row for the vessel itself is kept alongside them so the reason is visible rather
    than inferred from an absence.
    """
    jg = _junction_graph(tree)
    route = aorta_route(tree, jg)
    ref = ref or aorta_reference(mesh, tree, route=route, jg=jg)
    queue = list(vessels if vessels is not None else first_order_vessels(
        tree, jg=jg, route=route, mesh=mesh))
    names = names or {}

    rows, used = [], []
    seen = set()
    while queue:
        v = queue.pop(0)
        if v["segment"] in seen:
            continue
        seen.add(v["segment"])
        row = _evaluate(tree, v, k, ref, names, clamp)

        # A shared origin hands down on principle, not because it failed to cut: a stub
        # on it would stand for three great vessels where every other stub stands for
        # at most two.
        if (row["ok"] and descend and v["depth"] < max_depth
                and len(v.get("swallowed") or []) >= SHARED_ORIGIN_VESSELS):
            kids = daughters_of(tree, v, jg=jg, mesh=mesh)
            if kids:
                row.update(ok=False, can_descend=True, shared_origin=True,
                           handed_to=[c["segment"] for c in kids],
                           why=f"a shared origin feeding {len(v['swallowed'])} great "
                               f"vessels, so its {len(kids)} daughters carry the stubs")
                queue.extend(kids)
                rows.append(row)
                continue

        # A trunk that can only give a stub too short to be one hands down instead.
        # KY244_1's bovine brachiocephalic clears the aorta 2.5 mm before it splits, so
        # its own nub is 1.0 mm proud against a cohort median of 13.0 -- placeable, and
        # not a nub. Its two daughters are already clear of the aorta where they leave
        # it and take 11 to 13 mm stubs without complaint. Two real stubs beat one that
        # is barely there.
        if row["ok"] and row.get("marginal") and descend and v["depth"] < max_depth:
            kids = daughters_of(tree, v, jg=jg, mesh=mesh)
            if kids:
                row.update(ok=False, can_descend=True,
                           handed_to=[c["segment"] for c in kids],
                           why=f"only {row['proud_mm']:.1f} mm proud, too short to be a "
                               f"nub; its {len(kids)} daughters take one instead")
                queue.extend(kids)
                rows.append(row)
                continue

        if row["ok"]:
            rows.append(row)
            used.append(v)
            continue

        kids = (daughters_of(tree, v, jg=jg, mesh=mesh)
                if descend and row.get("can_descend") and v["depth"] < max_depth
                else [])
        if kids:
            row["handed_to"] = [c["segment"] for c in kids]
            row["why"] += f"; nubbing its {len(kids)} daughters instead"
            queue.extend(kids)
        rows.append(row)

    rows.sort(key=lambda r: (not r["ok"], r["depth"], r["name"]))
    return rows, used


def nub_terminals(tree, plan, vessels, anat=None, trunk_offset_mm=TRUNK_OFFSET_MM,
                  trunk_mode=TRUNK_MODE):
    """(terminals to cut, offsets) for a plan: trunk ends squared off, branches nubbed.

    The aorta's own two ends are NOT nubbed -- they are the model's inlet and outlet.
    They are cut *trunk_offset_mm* in from each end on a HORIZONTAL plane, which is
    what the great vessels' own exports do and what the paper wants the ascending and
    descending ends to look like. Openings downstream of a nubbed trunk are dropped
    from the list, since that one cut removes all of them.
    """
    by_seg = {v["segment"]: v for v in vessels}
    offsets, terms = {}, []
    swallowed = set()

    for row in plan:
        if not row["ok"]:
            continue
        v = by_seg[row["segment"]]
        seg = tree.segments[row["segment"]]
        if v["interior"]:
            t = Terminal(seg, v["at_end"], label=f"{seg.name}.nub", mode="normal",
                         name=row["name"])
            t.is_interior = True
            t.subtree_segments = set(v["subtree"])
            swallowed.update(x.label for x in v["swallowed"])
        else:
            t = v["terminal"]
            if t is None:
                continue
        terms.append(t)
        offsets[t.label] = float(row["offset_mm"])

    # The aorta's own ends, and any branch the plan could not place, keep the cut they
    # already had. An interior segment never carries a terminal of its own -- both its
    # ends are junctions -- so a nub on one can only collide with the openings it
    # swallows, which are already excluded.
    for t in tree.terminals:
        if t.label in swallowed or t.label in offsets:
            continue
        if getattr(t, "is_trunk", False):
            # The mode is set rather than assumed. classify_terminals already picks
            # "z" for the two aorta ends, but only when it is asked to, and a case
            # where it did not would otherwise be cut perpendicular to a centerline
            # that is 17-25 degrees off horizontal there -- a different-looking end
            # on one case out of thirty, which is exactly the kind of thing nobody
            # spots in a cohort until the figures are drawn.
            t.mode = trunk_mode
            offsets[t.label] = float(trunk_offset_mm)
        terms.append(t)
    return terms, offsets


def cut_nubs(mesh, tree, k=NUB_K, default_offset_mm=6.0, snap_mm=None,
             rescue=True, passes=2, anat=None, trunk_offset_mm=TRUNK_OFFSET_MM,
             trunk_mode=TRUNK_MODE, verbose=False):
    """Plan and perform a whole nub-in cut.

    Returns (mesh, plan, notes, info), where info carries the offsets and terminal
    list actually cut -- neither is reproducible from the plan once a rescue has
    rewritten an offset or a refusal has dropped a terminal.

    *rescue* lets a nub that will not cut cleanly move OUTWARD to the deepest offset
    that does. Outward is the safe direction here: a smaller offset is a plane
    further from the take-off, so a rescued nub is longer than asked for, never
    shorter, and never ends up back inside the aorta. Every rescue is recorded,
    because a nub that had to move is a case for the reviewer to look at.

    *passes* is how many times the assembled model may be re-cut to fix a rim that
    only comes out torn once its neighbours exist.
    """
    # The desktop tool snaps near-coincident rim vertices unless its "clean" box is
    # cleared, and the batch did not, so the same case cut both ways came out 66 faces
    # apart on KY160 with every offset identical. Whichever is preferable, they cannot
    # differ: the workbench hands these cases to that window, and a mesh that changes
    # the moment it is opened is a mesh nobody can vouch for. Follow the window.
    if snap_mm is None:
        snap_mm = snap_tolerance(mesh)

    jg = _junction_graph(tree)
    route = aorta_route(tree, jg)
    ref = aorta_reference(mesh, tree, route=route, jg=jg)
    vessels = first_order_vessels(tree, jg=jg, route=route, mesh=mesh)
    plan, used = nub_plan(mesh, tree, k=k, ref=ref, vessels=vessels, names=anat)
    terms, offsets = nub_terminals(tree, plan, used, anat=anat,
                                   trunk_offset_mm=trunk_offset_mm,
                                   trunk_mode=trunk_mode)

    by_label = {r["segment"]: r for r in plan}
    notes = []
    dropped = []             # terminals that refused; they have no opening to check

    # A horizontal plane is right for an aortic end until the vessel is steeply
    # oblique there, when a z cut stops going across the vessel and starts running
    # along it. Three cases in thirty: 52, 56 and 58 degrees against the 50 the
    # cutter allows, and because preflight refuses before anything is cut, each one
    # took its whole case down. Those ends fall back to a centerline-perpendicular
    # cut at the same 8 mm.
    cont = continuation_segments(tree, jg)
    attr0 = vessel_attribution(mesh, tree)
    for t in terms:
        if not getattr(t, "is_trunk", False) or t.label not in offsets:
            continue
        if t.mode != trunk_mode:
            continue
        want = float(offsets[t.label])
        if float(obliquity_deg(t, want)) <= MAX_OBLIQUITY_DEG:
            continue                       # horizontal is fine where it was asked for

        # Obliquity is not a property of the vessel, it is a property of the STATION:
        # KK40_3's ascending meets a horizontal plane at 49 degrees 2 mm in and at 17
        # degrees 40 mm in. So look for somewhere along the vessel where a horizontal
        # cut is both within the limit and actually clean, and only concede a
        # perpendicular plane when there is nowhere. The ends are wanted horizontal.
        L = float(t.segment.length)
        cands = np.arange(MIN_OFFSET_MM + OFFSET_GRID_MM, 0.95 * L, GROW_FINE_MM)
        found = None
        for cand in sorted(cands, key=lambda c: abs(c - want)):
            if float(obliquity_deg(t, float(cand))) > MAX_OBLIQUITY_DEG:
                continue
            if nub_works(mesh, tree, t, float(cand), attr0,
                         own_extra=cont.get(t.segment.name)):
                found = float(_grid(cand))
                break
        if found is not None:
            offsets[t.label] = found
            notes.append(f"{t.name}: kept horizontal, moved from {want:.1f} to "
                         f"{found:.1f} mm where the plane meets the vessel at "
                         f"{obliquity_deg(t, found):.0f} deg")
        else:
            ob = float(obliquity_deg(t, want))
            t.mode = "normal"
            notes.append(f"{t.name}: no horizontal cut works anywhere on this vessel "
                         f"({ob:.0f} deg at {want:.0f} mm), cut perpendicular instead")

    # Recorded, not raised. Preflight catches configurations worth knowing about, but
    # every individual cut is defended below and every opening is measured afterwards,
    # so a complaint here is a thing to show the reviewer -- not a reason to return no
    # geometry at all for a case whose branches are all perfectly cuttable.
    pre = preflight(tree, offsets, default_offset_mm, terminals=terms)
    if pre is not None:
        notes.append(f"preflight: {pre}")

    out = _assemble(mesh, tree, terms, offsets, default_offset_mm, snap_mm,
                    by_label, notes, dropped, rescue, cont)

    # Verify the ASSEMBLED model, and grow whatever came out torn. The probe above
    # judges each nub against the uncut mesh, one cut at a time, and that is not
    # quite the model that gets built: a neighbouring cut can already have taken
    # faces this one's band needed, so a nub that probes clean can still assemble
    # ragged. Measured over the thoracic cohort it is 2 nubs in 80, both on trunks
    # whose cut also removes a subtree, and both come good by moving 2 mm further
    # out. Only the torn ones pay for the retry.
    for _ in range(max(int(passes) - 1, 0)):
        torn = _torn(out, tree, offsets, terms)
        if not torn:
            break
        moved = False
        for t in torn:
            row = by_label.get(t.segment.name)
            if row is None or t.label not in offsets:
                continue
            step = max(RETRY_STEP_MM, 0.0)
            new_off = offsets[t.label] - step
            if new_off < MIN_OFFSET_MM:
                continue
            offsets[t.label] = new_off
            row["offset_mm"] = new_off
            row["n_mm"] = float(row["length_mm"] - new_off)
            row["proud_mm"] = float(row["n_mm"] - row["s_wall_mm"])
            row["retried_mm"] = float(row.get("retried_mm", 0.0) + step)
            notes.append(f"{row['name']}: rim came out torn in the assembled model, "
                         f"grown a further {step:.1f} mm")
            moved = True
        if not moved:
            break
        out = _assemble(mesh, tree, terms, offsets, default_offset_mm, snap_mm,
                        by_label, notes, dropped, rescue=False, cont=cont)

    # A terminal that refused made no opening, so leaving it in the list tells
    # validate_cut to expect one and it reports a hole count that is short by
    # exactly the cuts it was already told about.
    terms = [t for t in terms if t not in dropped]
    problem = validate_cut(out, tree, offsets, default_offset_mm, terminals=terms)
    if problem:
        notes.append(f"validate: {problem}")
    if verbose:
        print(f"  {len(boundary_loops(out))} openings, euler={int(out.euler_number)}"
              f"{'' if problem is None else '  ** ' + problem + ' **'}")
    # The offsets and terminal list come back because they are not reproducible from
    # the plan alone: a rescue rewrites an offset, a refusal drops a terminal, and
    # anything scoring the result has to measure the planes that were actually cut.
    return out, plan, notes, {"offsets": offsets, "terminals": terms,
                              "validate": problem}


def _base_loops(mesh):
    """How many openings this mesh has before anything is cut. Cached per mesh."""
    hit = getattr(mesh, "_nub_base_loops", None)
    if hit is None:
        hit = len(boundary_loops(mesh))
        try:
            mesh._nub_base_loops = hit
        except Exception:                                          # noqa: BLE001
            pass
    return hit


def nub_works(mesh, tree, terminal, off, attr=None, band_scale=2.5, own_extra=None,
              band_mm=None, tilt=None):
    """Would a cut at this offset produce a nub worth keeping?

    The same three questions the cohort benchmark asks, in the cheap order: does the
    cut take any surface off a vessel it is not cutting, does it leave exactly one
    opening, and does that opening's rim lie on the plane that made it.
    """
    if off < MIN_OFFSET_MM:
        return False
    if band_mm is not None or tilt is not None:
        was_band = terminal.band_mm
        was_tilt = (terminal.tilt_a_deg, terminal.tilt_b_deg)
        if band_mm is not None:
            terminal.band_mm = float(band_mm)
        if tilt is not None:
            terminal.tilt_a_deg, terminal.tilt_b_deg = float(tilt[0]), float(tilt[1])
        try:
            return nub_works(mesh, tree, terminal, off, attr, band_scale, own_extra)
        finally:
            terminal.band_mm = was_band
            terminal.tilt_a_deg, terminal.tilt_b_deg = was_tilt
    attr = attr if attr is not None else vessel_attribution(mesh, tree)
    owner = attr[0]
    own = owner == terminal.segment.name
    sub = getattr(terminal, "subtree_segments", None)
    if sub:
        own = own | np.isin(owner, list(sub))
    # The rest of the same vessel, where the exported line was broken at a
    # non-junction. Not passed through subtree_segments, which would mean "remove all
    # of this" rather than "this is still me".
    if own_extra:
        own = own | np.isin(owner, list(own_extra))
    try:
        band, delete, _ = cut_terminal(mesh, tree, terminal, off,
                                       band_scale=band_scale, _probe=True, _attr=attr)
    except ValueError:
        return False
    if not band.any():
        return False
    if float(mesh.area_faces[(delete | band) & ~own].sum()) > FOREIGN_TOLERANCE_MM2:
        return False
    try:
        out = cut_terminal(mesh, tree, terminal, off, band_scale=band_scale,
                           snap_mm=snap_tolerance(mesh), _attr=attr, _suggest=False)
    except ValueError:
        return False
    # ONE MORE OPENING THAN THE MESH ALREADY HAD -- not "euler == 1", which is only
    # that statement for a plain closed surface of genus 0. Two cohort meshes are not:
    # KY238_1 is watertight with euler 0, a tunnel through it, and KK40_3 is watertight
    # with euler 4, two separate closed components. Cutting either gives an euler that
    # is right for the shape and wrong for the test, so every candidate was rejected --
    # KY238_1's two branches were being refused while producing a rim 0.000 mm off
    # plane. Counting openings says what was meant and is true of any topology.
    loops = boundary_loops(out)
    if len(loops) != _base_loops(mesh) + 1:
        return False
    pt, nrm = terminal.cut_plane(off)
    verts = np.asarray(out.vertices)
    rim = min(loops, key=lambda L_:
              float(np.linalg.norm(verts[np.asarray(L_)].mean(axis=0) - pt)))
    return float(np.abs((verts[np.asarray(rim)] - pt) @ nrm).max()) <= RIM_FLATNESS_MM


def place_end(mesh, tree, terminal, want_off, attr=None, own_extra=None, limit=None):
    """The offset nearest *want_off* that cuts this vessel END cleanly, or None.

    Separate from grow_to_fit because the safe direction is the opposite one. A nub
    is made safer by moving AWAY from its ostium, which is a smaller offset; a
    trunk end is made safer by moving deeper INTO the vessel, which is a larger one.
    suggest_pullback only ever reduces, so asked to rescue DJ18_2's aortic root it
    walked to 0.8 mm -- onto the vessel tip, where validation then refused it for
    having nothing left to remove. So this searches outward in BOTH directions and
    takes whichever working offset is closest to the one asked for.
    """
    attr = attr if attr is not None else vessel_attribution(mesh, tree)
    L = float(terminal.segment.length)
    hi = float(limit if limit is not None else 0.95 * L)
    # Down to the smallest offset validation will accept, in the finest step the
    # sliders can hold. Borrowing the nub margin here (1.5 mm, in half-millimetre
    # steps) missed DJ18_2's aortic root entirely: the only offsets that cut it
    # cleanly are 1.2 and 1.4 mm, and the search reported that no offset worked at
    # all. An aortic end is not a nub -- nothing has to survive past it -- so the
    # only floor that applies is the one that stops a cut landing on the tip itself.
    lo = MIN_OFFSET_MM + OFFSET_GRID_MM
    auto = float(band_width(terminal))
    cands = sorted(np.arange(lo, hi + 1e-9, OFFSET_GRID_MM),
                   key=lambda c: abs(c - want_off))
    for f in BAND_TRIES:
        bm = None if f == 1.0 else _grid(auto * f)
        if nub_works(mesh, tree, terminal, want_off, attr, own_extra=own_extra,
                     band_mm=bm):
            return {"offset": want_off, "grow": 0.0, "band_mm": bm, "tilt": None}
        for cand in cands:
            if nub_works(mesh, tree, terminal, float(cand), attr,
                         own_extra=own_extra, band_mm=bm):
                return {"offset": float(_grid(cand)),
                        "grow": float(cand - want_off), "band_mm": bm, "tilt": None}
    return None


def grow_to_fit(mesh, tree, terminal, off, attr=None, own_extra=None):
    """The SHORTEST nub at or past *off* that cuts cleanly, or None.

    Searched over nub length in fixed millimetre steps rather than by
    suggest_pullback's ladder of offset fractions. The ladder is the right tool for
    "pull this back until it works" from a hand-driven slider, and the wrong one
    here, because where it lands depends on where it started and the acceptance test
    is not monotone in the offset: KY160's left subclavian grew 6.9 mm when asked
    for k = 0.5 and 43.0 mm when asked for k = 0.75, a LONGER nub from a longer
    request. A preset that promises a uniform stub cannot have that in it.

    Coarse then fine, because most vessels need no growth at all and the ones that
    do are usually within a few millimetres.
    """
    attr = attr if attr is not None else vessel_attribution(mesh, tree)
    floor = max(NUB_MARGIN_MM, MIN_OFFSET_MM)
    auto = float(band_width(terminal))
    # Narrowest plane first: a wider one is only ever needed where the automatic width
    # is wrong, and widening it needlessly reaches further around the vessel than the
    # cut has any business doing.
    for f in BAND_TRIES:
        bm = None if f == 1.0 else _grid(auto * f)
        if nub_works(mesh, tree, terminal, off, attr, own_extra=own_extra, band_mm=bm):
            return {"offset": off, "grow": 0.0, "band_mm": bm, "tilt": None}

    coarse, band = None, None
    for f in BAND_TRIES:
        bm = None if f == 1.0 else _grid(auto * f)
        for grow in np.arange(GROW_STEP_MM, GROW_LIMIT_MM + 1e-9, GROW_STEP_MM):
            cand = off - float(grow)
            if cand < floor:
                break
            if nub_works(mesh, tree, terminal, cand, attr, own_extra=own_extra,
                         band_mm=bm):
                coarse, band = cand, bm
                break
        if coarse is not None:
            break

    if coarse is None:
        # Nothing flat, at any offset or width. Tilt the plane to meet the vessel's
        # own cross-section, nearest-to-upright first so a needless tilt is never
        # taken. Only reached by a vessel that would otherwise be left uncut.
        for ta, tb in sorted(((a, b) for a in TILT_TRIES_DEG for b in TILT_TRIES_DEG),
                             key=lambda p: abs(p[0]) + abs(p[1])):
            if ta == 0 and tb == 0:
                continue
            for f in BAND_TRIES:
                bm = None if f == 1.0 else _grid(auto * f)
                for cand in np.arange(floor, off + 1e-9, GROW_FINE_MM):
                    if nub_works(mesh, tree, terminal, float(cand), attr,
                                 own_extra=own_extra, band_mm=bm, tilt=(ta, tb)):
                        return {"offset": float(_grid(cand)),
                                "grow": float(off - cand), "band_mm": bm,
                                "tilt": (float(ta), float(tb))}
        return None

    # Back toward the ostium in finer steps: the coarse hit is at most one step past
    # the shortest nub that works, and a millimetre of vessel is worth keeping.
    best = coarse
    for back in np.arange(GROW_FINE_MM, GROW_STEP_MM, GROW_FINE_MM):
        cand = coarse + float(back)
        if cand >= off:
            break
        if nub_works(mesh, tree, terminal, cand, attr, own_extra=own_extra,
                     band_mm=band):
            best = cand
    return {"offset": best, "grow": off - best, "band_mm": band, "tilt": None}


def _torn(out, tree, offsets, terms):
    """Which cut terminals ended up with a rim that does not lie on its plane."""
    loops = _loops_by_terminal(out, tree, offsets, terminals=terms)
    verts = np.asarray(out.vertices)
    bad = []
    for t in terms:
        if t.label not in offsets:
            continue
        loop = loops.get(t.label)
        if loop is None or len(loop) == 0:
            continue
        pt, nrm = t.cut_plane(float(offsets[t.label]))
        if float(np.abs((verts[np.asarray(loop)] - pt) @ nrm).max()) > RIM_FLATNESS_MM:
            bad.append(t)
    return bad


def _assemble(mesh, tree, terms, offsets, default_offset_mm, snap_mm,
              by_label, notes, dropped, rescue, cont=None):
    """Cut every terminal, from the original mesh, at whatever offsets now say."""
    cont = cont if cont is not None else continuation_segments(tree)
    out = mesh.copy()
    for t in terms:
        if t in dropped:
            continue
        off = float(offsets.get(t.label, default_offset_mm))
        row = by_label.get(t.segment.name)
        is_nub = row is not None and t.label in offsets

        # Every nub is PROBED before it is cut, not only the ones that raise. A cut
        # only refuses outright past TORN_REFUSE_MM; a 2 mm tear cuts happily and
        # scores as a valid mesh, and at nub depth -- 85-95% along the branch, the
        # worst regime this cutter has -- that is the common failure, not the rare
        # one. suggest_pullback holds a candidate to the rim being ON its plane and
        # to taking no more than a square millimetre of any other vessel, which is
        # the standard the paper needs, and it returns the offset asked for
        # unchanged when that offset is already clean.
        if is_nub and rescue:
            try:
                # Probed against the ORIGINAL mesh, never the partly cut one. The
                # acceptance test includes euler == 1, which is a statement about a
                # surface carrying exactly one opening -- true of a single cut on the
                # closed mesh, false the moment a second cut exists. Run on the
                # running mesh it rejected every nub after the first, four of seven
                # on the first three cases, for having neighbours. Cut order is
                # already known not to matter here, so the base mesh is the honest
                # thing to ask anyway, and its attribution is cached.
                found = grow_to_fit(mesh, tree, t, off,
                                    own_extra=cont.get(t.segment.name))
            except Exception:                                      # noqa: BLE001
                found = None
            if found is None:
                row.update(ok=False, why="no offset on this vessel cuts cleanly")
                notes.append(f"{row['name']}: no clean offset, left uncut")
                offsets.pop(t.label, None)
                dropped.append(t)
                continue
            # A wider plane where the automatic one could not span the vessel, and a
            # tilt where nothing flat would cut at all. Both are set on the terminal so
            # the cut uses them and recorded so the window can be given the same ones:
            # they are its Plane size and Tilt A/B controls.
            if found.get("band_mm") is not None:
                t.band_mm = float(found["band_mm"])
                if row.get("band_mm") != t.band_mm:
                    row["band_mm"] = float(t.band_mm)
                    notes.append(f"{row['name']}: plane widened to {t.band_mm:.1f} mm "
                                 f"(automatic width {band_width(t):.1f} mm could not "
                                 f"span this vessel)")
            if found.get("tilt") is not None:
                t.tilt_a_deg, t.tilt_b_deg = found["tilt"]
                row["tilt_a_deg"], row["tilt_b_deg"] = found["tilt"]
                notes.append(f"{row['name']}: plane tilted {t.tilt_a_deg:+.0f}/"
                             f"{t.tilt_b_deg:+.0f} deg; no flat plane cuts this vessel "
                             f"anywhere along it")
            new_off = float(found["offset"])
            if new_off < off - 1e-6:
                grew = off - new_off
                # The offsets dict is not bookkeeping -- validate_cut and the opening
                # matcher rebuild every plane from it, so a rescue recorded only in
                # the plan gets validated against the plane already rejected.
                offsets[t.label] = new_off
                row["offset_mm"] = new_off
                row["n_mm"] = float(row["length_mm"] - new_off)
                row["rescued_mm"] = float(grew)
                row["proud_mm"] = float(row["n_mm"] - row["s_wall_mm"])
                notes.append(f"{row['name']}: nub grown {grew:.1f} mm to cut cleanly "
                             f"({row['proud_mm']:.1f} mm proud)")
                off = new_off

        try:
            out = cut_terminal(out, tree, t, off, snap_mm=snap_mm)
            continue
        except Exception as exc:                                   # noqa: BLE001
            first = exc

        if row is not None:
            row.update(ok=False, why=f"refused: {first}")
            notes.append(f"{row['name']}: refused, left uncut ({first})")
            offsets.pop(t.label, None)
            dropped.append(t)
            continue

        # Not a nub -- an aorta end, or a branch the plan left at its default. These
        # used to re-raise, which cost the WHOLE case: DJ18_2's distal aorta comes out
        # 15 mm off plane at every offset from 2 to 10 mm, a property of that case's
        # mesh and nothing to do with the preset, and it took all three of its nubs
        # down with it. Try to place it somewhere workable, and if there is nowhere,
        # leave that one end closed and say so. Three nubs and an unopened aorta is a
        # case for the reviewer; no case at all is a case nobody sees.
        try:
            found = place_end(mesh, tree, t, off,
                              own_extra=cont.get(t.segment.name))
        except Exception:                                          # noqa: BLE001
            found = None
        if found is not None:
            if found.get("band_mm") is not None:
                t.band_mm = float(found["band_mm"])
                notes.append(f"{t.name}: plane widened to {t.band_mm:.1f} mm to span "
                             f"this vessel")
            if found.get("tilt") is not None:
                t.tilt_a_deg, t.tilt_b_deg = found["tilt"]
                notes.append(f"{t.name}: plane tilted {t.tilt_a_deg:+.0f}/"
                             f"{t.tilt_b_deg:+.0f} deg to cut cleanly")
            if abs(float(found["offset"]) - off) > 1e-6:
                offsets[t.label] = float(found["offset"])
                notes.append(f"{t.name}: moved from {off:.1f} to "
                             f"{found['offset']:.1f} mm to cut cleanly")
            try:
                out = cut_terminal(out, tree, t, float(offsets[t.label]),
                                   snap_mm=snap_mm)
                continue
            except Exception:                                      # noqa: BLE001
                pass
        notes.append(f"{t.name}: NO offset on this vessel cuts cleanly, left closed "
                     f"({first})")
        offsets.pop(t.label, None)
        dropped.append(t)

    return keep_largest_component(out)
