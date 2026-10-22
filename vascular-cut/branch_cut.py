# Cutting a branched aorta mesh open at its centerline terminals

# Each terminal gets a plane taken from the centerline (point + tangent) and the
# cap beyond it is clipped away, leaving a clean planar opening.

# Why the cut has to be LOCAL: an infinite plane placed at the arch height slices
# every vessel it happens to cross -- a global slice_plane there opens five holes
# at once and disconnects branches from the trunk. So for each terminal we first
# isolate the cap region (the faces reaching past the plane that are surface-
# connected to that terminal's tip) and clip only those, splicing the result back
# onto the untouched remainder.

import numpy as np
import trimesh
from pathlib import Path
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from trimesh.intersections import slice_faces_plane

from branch_tree import Terminal


def _surface_from_cells(cells, faces_per_cell):
    """Outer surface of a solid mesh: the faces used by exactly ONE cell.

    An interior face is shared by two cells and cancels; what is left is the wall.
    """
    faces = np.vstack([cells[:, list(f)] for f in faces_per_cell])
    _, idx, counts = np.unique(np.sort(faces, axis=1), axis=0,
                               return_index=True, return_counts=True)
    return faces[idx[counts == 1]]


# A tetrahedron's four faces, and a hexahedron's six (as two triangles each) in the
# standard Simpleware/Abaqus node order.
_TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
_HEX_FACES = [(0, 1, 2), (0, 2, 3), (4, 5, 6), (4, 6, 7),
              (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
              (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)]


def load_mat_mesh(path, key=None):
    """Load a Simpleware-exported MATLAB mesh as a triangular SURFACE.

    Handles both kinds of export seen in the wild:
      * a surface mesh   (*_surface_vertices / *_surface_triangles)
      * a solid FE mesh  (*_FE_vertices / *_FE_tetrahedra or *_FE_hexahedra),
        from which the wall is extracted -- the cut tool only ever needs the
        surface, and refusing an FE file outright just blocks the user for a
        reason they cannot act on.

    The variable prefix is discovered rather than assumed ("Model1" in the sample
    file, but there is no guarantee another export uses the same name).
    """
    import scipy.io as sio

    m = sio.loadmat(str(path))
    if key is None:
        keys = [k[: -len("_surface_vertices")] for k in m
                if k.endswith("_surface_vertices")]
        if keys:
            key = keys[0]
        else:
            for suffix, cell_faces in (("_FE_tetrahedra", _TET_FACES),
                                       ("_FE_hexahedra", _HEX_FACES)):
                soli = [k[: -len(suffix)] for k in m if k.endswith(suffix)]
                if not soli:
                    continue
                k0 = soli[0]
                vkey = f"{k0}_FE_vertices"
                if vkey not in m:
                    continue
                V = np.asarray(m[vkey], dtype=float)
                C = np.asarray(m[f"{k0}{suffix}"]).astype(np.int64)
                if C.size and C.min() == 1:          # MATLAB is 1-based
                    C = C - 1
                mesh = trimesh.Trimesh(
                    vertices=V, faces=_surface_from_cells(C, cell_faces),
                    process=False)
                mesh.remove_unreferenced_vertices()   # drop the solid's interior
                mesh.fix_normals()                    # cells carry no consistent winding
                return mesh
            raise ValueError(
                f"{Path(path).name}: no surface or FE mesh found. Expected "
                f"'*_surface_vertices' + '*_surface_triangles', or '*_FE_vertices' "
                f"+ '*_FE_tetrahedra'/'*_FE_hexahedra'. This file holds: "
                f"{[k for k in m if not k.startswith('__')]}"
            )

    V = np.asarray(m[f"{key}_surface_vertices"], dtype=float)
    F = np.asarray(m[f"{key}_surface_triangles"]).astype(np.int64)
    if F.min() == 1:                      # MATLAB is 1-based
        F = F - 1
    return trimesh.Trimesh(vertices=V, faces=F, process=False)


# How much of the centerline must lie within the mesh's own bounds before the two are
# accepted as the same case. Measured, and the separation is total: over 79 correctly
# paired cases across both lab cohorts the coverage is 100.0% on every single one --
# minimum, not median. The four known-bad pairings in the thoracic cohort, where a
# branched centerline was matched with the trunk-only mesh beside it, come out at
# 61.3, 72.0, 72.7 and 81.1%.
#
# This is the check that works everywhere. The radius sign only catches a mismatch on
# a WATERTIGHT surface: on an open one segment_radius falls back to unsigned distance
# and returns cheerfully positive nonsense -- a cropped KY160 measured 79 mm on a
# 22 mm aorta and was flagged by nothing.
#
# It is necessary, not sufficient: a wall or solid FE model has the centerline inside
# the bounding box and outside the material, which is what _no_radius_reason is for.
MIN_CENTERLINE_COVERAGE = 0.98


def centerline_coverage(tree, mesh, pad_frac=0.02):
    """Fraction of centerline points inside the mesh's bounds, with a little slack.

    The slack is a wall thickness, not a free pass: a centerline reaching the very end
    of the vessel it describes can sit a fraction outside the surface's bounding box
    without anything being wrong.
    """
    try:
        pts = np.vstack([s.pts for s in tree.segments.values()])
        lo, hi = np.asarray(mesh.bounds, float)
        pad = float(pad_frac) * (hi - lo)
        return float(((pts >= lo - pad) & (pts <= hi + pad)).all(axis=1).mean())
    except Exception:                                     # noqa: BLE001
        return 1.0            # unmeasurable is not evidence of a mismatch


def _no_radius_reason(mesh, tree):
    """Say WHY no terminal could be measured, distinguishing the two real causes.

    Both present as "zero radius everywhere", but the fixes are opposite:
      * a wall / solid FE model -- the centerline is right, it runs down the lumen,
        but the lumen is a HOLE in the supplied geometry. The user needs the inner
        surface, not the solid.
      * a genuine mismatch -- different case, different frame, or wrong units.
    """
    generic = ("Check the CSV and the mesh are the same case, in the same "
               "coordinate frame and the same units.")
    try:
        pts = np.vstack([s.pts for s in tree.segments.values()])
        lo, hi = mesh.bounds
        mid = pts[len(pts) // 2]
        # Sitting comfortably inside the mesh's own bounding box, yet outside its
        # material, means the geometry wraps around the centerline: a wall model.
        span = np.maximum(hi - lo, 1e-9)
        frac = (mid - lo) / span
        enclosed = bool(np.all(frac > 0.05) and np.all(frac < 0.95))
        if enclosed and not bool(mesh.contains([mid])[0]):
            return (
                "every opening measures zero radius, but the centerline runs "
                "INSIDE the mesh's bounding box, this looks like a wall or solid "
                "FE model, where the lumen is a hole rather than the surface. "
                "This tool cuts a lumen surface: export the inner (blood-facing) "
                "surface, or the segmentation the centerline was built from."
            )
    except Exception:
        pass
    return ("every terminal measures zero radius, the centerline appears to lie "
            "outside the mesh. " + generic)


def boundary_loops(mesh):
    """Vertex-index loops of the open boundary (edges used by exactly one face)."""
    lone = trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)
    if len(lone) == 0:
        return []
    edges = mesh.edges_sorted[lone]

    n = len(mesh.vertices)
    g = coo_matrix(
        (np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n, n)
    )
    _ncomp, labels = connected_components(g, directed=False)

    # Group the boundary vertices by component in one sort. The graph's nodes are
    # EVERY vertex, so every interior vertex is a component of its own, and the loop
    # this replaced -- one pass over the boundary per component -- ran once per
    # vertex: 6.3 of the 17 s in the first cut of a 159k-face aorta. Same loops, same
    # order (by label, members ascending); tests/test_boundary_loops.py holds both.
    used = np.unique(edges)
    lab = labels[used]
    order = np.lexsort((used, lab))
    lab, used = lab[order], used[order]
    groups = np.split(used, np.flatnonzero(np.diff(lab)) + 1)
    return [m.copy() for m in groups if len(m) >= 3]


def opening_stats(mesh):
    """Per-opening centroid / perimeter / equivalent diameter."""
    out = []
    for loop in boundary_loops(mesh):
        pts = mesh.vertices[loop]
        c = pts.mean(axis=0)
        r = np.linalg.norm(pts - c, axis=1)
        out.append(
            dict(
                n=len(loop),
                centroid=c,
                radius_mean=float(r.mean()),
                diameter=float(2 * r.mean()),
            )
        )
    return out


# How many rings of straddling faces the rim may be completed by, outwards from the
# region the cut removes. It exists because ownership by nearest centerline goes wrong
# for a face or two either side of an ostium, and the rim has to be clipped there
# anyway. Small on purpose: this is a tolerance for that error, not a reach. Raising
# it lets one cut walk along the plane onto a neighbouring vessel, which is the
# failure it was introduced to stop.
RIM_GROW_LAYERS = 2


def band_width(terminal, band_scale=2.5, radius_mm=None):
    """How far around its plane a cut reaches, the 'plane size' the UI exposes.

    Shrinking this is a real remedy, not just a knob: it rescues cuts that the
    automatic width gets wrong (a plane that reaches too far on a tortuous branch
    grabs a stretch that doubles back and clips it too, punching a second hole).

    *radius_mm* is the vessel's radius WHERE THE CUT IS, when the caller knows it. A
    branch flares as it leaves its parent -- KY160's left subclavian is 6 mm at its
    origin and 3.2 mm in its middle -- so the width derived from a segment-wide median
    is too small exactly where the cross-section is biggest, and the rim runs out of
    band and tears at precisely `band_w` from the cut point. That was the 12 mm wall.
    """
    if terminal.band_mm:
        return float(terminal.band_mm)
    r = radius_mm if radius_mm else getattr(terminal.segment, "_radius_mm", None)
    r = max(float(r or 5.0), 1.0)
    if getattr(terminal, "is_interior", False):
        # No 12 mm floor for a pre-bifurcation trim. That floor exists to guarantee
        # ring coverage where a plane meets a FAT trunk obliquely; on a thin, hairpin
        # vessel it instead reaches a stretch that doubles back near the cut, the
        # plane clips that too, and one trim punches two holes.
        return max(band_scale * r, 3.0)
    # The same floor misfires the same way on a SHORT branch: it is a fixed distance
    # that takes no account of how much vessel there is to reach across. Measured on
    # 40 real aortas, every "one cut punched two holes" failure was a thin branch
    # (r 1.6-3.5 mm) on a stub of 15-29 mm, where the 12 mm floor spans most of the
    # vessel and necessarily reaches the parent and its neighbours. Past half a
    # segment the band is at the far junction anyway, so bound the floor there --
    # a fat trunk is long enough that this never binds and it keeps its 12 mm.
    seg_len = float(getattr(terminal.segment, "length", 0.0) or 0.0)
    floor = 12.0 if seg_len <= 0.0 else min(12.0, max(0.5 * seg_len, 3.0))
    return max(band_scale * r, floor)


# How far a region boundary's reach may be scaled, nearest-first and alternating
# outward. 1.0 comes first so a region that already cuts keeps exactly the plane it
# would have had and pays one probe for the check -- 167 of 167 working cuts in the
# abdominal cohort were left untouched by this.
#
# It has to go BOTH ways, because the two failures pull in opposite directions:
#   torn rim   "the opening comes out N mm out of plane" -- the plane ran out of band
#              before it finished spanning the ring. Wants a WIDER reach.
#   severed    an extra opening -- the plane reached sideways through a bystander, on
#              an abdominal aorta a renal running parallel rather than leaving
#              laterally. Wants a NARROWER one.
# Measured over 182 region cuts: 10 of the 15 failures are rescued, JK4_2 at 0.60x
# and JK17_2 and JK3_1 at 1.30x.
REGION_REACH_SCALES = (1.0, 1.3, 0.75, 1.6, 0.6, 2.0, 0.5, 2.5, 0.4)
# The ladder stops at 0.4 and going lower was tried and measured. Extending it to
# 0.2 lets JK3_1's aortic neck stop severing a neighbour and moves the BOUNDS proxy
# by one -- and moves the number that matters, presets cutting end to end, by
# nothing at all: 146 of 182 either way, thoracic unchanged on both counts. A 0.2x
# band is about 7.5 mm on a 15 mm aorta, narrower than the vessel's own radius, so
# it stops severing by failing to reach rather than by cutting well, and validate_cut
# still refuses what it produces. A proxy point bought with an implausible plane.


def fit_region_reach(mesh, tree, bounds, snap_mm=None,
                     scales=REGION_REACH_SCALES):
    """Give a region's boundary planes a reach that actually cuts. Returns the scale.

    Sets ``band_mm`` on each bound and returns the multiplier used, or None if no
    reach works and the bands are left as they were. A bound whose ``band_mm`` was
    already set by hand is never touched -- a person who typed a plane size means it.

    Probed on a COPY with only the two bounds cut, because that is the only state in
    which the test is unambiguous: two planes through a closed surface leave exactly
    two openings, and anything more is a plane that severed a vessel it was not aimed
    at. Once the branches inside the region are also being cut, the count no longer
    says which plane did what.
    """
    if len(bounds) != 2 or mesh is None:
        return None
    if any(getattr(t, "band_mm", None) for t in bounds):
        return None                      # a hand-set plane size is the user's answer
    auto = [band_width(t) for t in bounds]
    if not all(a and a > 0 for a in auto):
        return None
    snap = snap_tolerance(mesh) if snap_mm is None else snap_mm
    for k in scales:
        for t, a in zip(bounds, auto):
            t.band_mm = float(a) * float(k)
        m = mesh
        try:
            for t in bounds:
                m = cut_terminal(m, tree, t, getattr(t, "_region_offset", 6.0),
                                 snap_mm=snap)
        except Exception:                                 # noqa: BLE001
            continue                     # refused at this reach; try the next
        if len(boundary_loops(m)) == 2:
            return float(k)
    for t in bounds:
        t.band_mm = None
    return None


# How far a region boundary may be MOVED along the aorta to find a cuttable
# cross-section, when no plane reach works where it was asked for.
#
# 8 mm is not a round number picked for comfort: region_plan already moves a boundary
# by OSTIUM_CLEAR_R renal radii to clear a branch mouth, which is about 8 mm on a real
# renal, so this stays within a distance the planner already moves boundaries by for
# its own reasons.
#
# The cap is the point. Measured over the five region cuts no reach rescues, three are
# rescued by moving: JK10_3's visceral segment by 1 mm, JK14_1's by 8, and JK6_1's
# AORTIC NECK by 12. The first two are conventions -- "around the celiac/SMA/renal
# group" is already a +/-15 mm convention -- and moving them changes nothing anyone
# reports. The third is not: an aortic neck BEGINS at the lowest renal by definition,
# and a neck starting 12 mm below the renal is a different measurement wearing the
# same name. So the cap deliberately falls between them, and JK6_1 stays refused
# rather than being silently answered with the wrong region.
REGION_NUDGE_MM = 8.0

# ...and it is OFF, because it was measured and it does not work.
#
# Judged on opening count alone, moving the boundary rescued three of the five cuts no
# reach saves: JK10_3 by 1 mm, JK14_1 by 8, JK6_1's aortic neck by 12. Judged on what
# the tool ACCEPTS -- validate_cut, which also checks rim flatness and where the
# opening landed -- it rescues none of them. All three were artifacts of a probe
# measuring something easier than the goal, and shipping it would have moved JK14_1's
# boundary 8 mm along the aorta to produce a cut the very next recut still refused.
#
# The mechanism is kept and the parameter is live, so a different cohort can be tried
# without rebuilding it. It stays at "do not move" until something measures otherwise.
REGION_NUDGES = (0.0,)
REGION_NUDGES_TRIED = (0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 5.0, -5.0, 8.0, -8.0)


def fit_region(mesh, tree, route, s0, s1, snap_mm=None, jg=None,
               nudges=REGION_NUDGES, scales=REGION_REACH_SCALES):
    """Stations and a plane reach that actually cut. Returns a dict or None.

    ``{"s0", "s1", "scale", "moved_mm"}``. The region is tried where it was asked for
    first, so anything that already cuts is returned unchanged with moved_mm 0.0.

    Both ends move TOGETHER. A neck moved 3 mm along the aorta is still a neck; a neck
    shortened by 3 mm is a different measurement, and silently changing a region's
    LENGTH to make it cuttable would corrupt the very thing being measured.
    """
    jg = jg or _junction_graph(tree)
    snap = snap_tolerance(mesh) if snap_mm is None else snap_mm
    open_distal = s1 is None
    for d in nudges:
        a = s0 + d
        b = None if open_distal else s1 + d
        if a < 0 or (b is not None and b > route.length):
            continue
        try:
            bounds, _in, _rm, _nt = region_plan(tree, a, b, route, jg, mesh=mesh,
                                                clear_ostia=True)
        except Exception:                                 # noqa: BLE001
            continue
        if not bounds:
            continue                # nothing to place; the vessel's ends bound it
        scale = (fit_region_reach(mesh, tree, bounds, snap_mm=snap, scales=scales)
                 if len(bounds) == 2 else None)
        m = mesh
        offs = {}
        try:
            for t in bounds:
                off = float(getattr(t, "_region_offset", 6.0))
                offs[t.label] = off
                m = cut_terminal(m, tree, t, off, snap_mm=snap)
        except Exception:                                 # noqa: BLE001
            continue
        if len(boundary_loops(m)) != len(bounds):
            continue                                      # severed a bystander
        # And it has to pass what the TOOL accepts, not merely what cut_terminal
        # allows. validate_cut checks rim flatness and where the opening ended up;
        # cut_terminal lets through cuts it then rejects, so a fitter judged on
        # opening COUNT alone reports a station that the very next recut refuses --
        # which is exactly what happened on JK14_1 and JK10_3, moved 8 mm and 1 mm to
        # no purpose because the probe was measuring something easier than the goal.
        if validate_cut(m, tree, offs, terminals=bounds) is None:
            return {"s0": a, "s1": b, "scale": scale or 1.0, "moved_mm": float(d)}
    return None


def _centerline_kd(tree):
    """KD-tree over every centerline point, cached on the tree (built once)."""
    cache = getattr(tree, "_kd_cache", None)
    if cache is None:
        names, pts, arcs = [], [], []
        for name, seg in tree.segments.items():
            names.extend([name] * len(seg.pts))
            pts.append(seg.pts)
            arcs.append(seg.s)
        cache = (cKDTree(np.vstack(pts)), np.array(names), np.concatenate(arcs))
        tree._kd_cache = cache
    return cache


def _segment_kd(seg):
    """KD-tree over ONE segment's centerline points, cached on the segment."""
    kd = getattr(seg, "_kd_own", None)
    if kd is None:
        kd = cKDTree(np.asarray(seg.pts, dtype=float))
        seg._kd_own = kd
    return kd


def _point_spacing(seg):
    """Median gap between this segment's centerline points, cached on the segment.

    Used as the tolerance on "is this point past the end", because it is the
    resolution at which the segment's extent is known in the first place.
    """
    sp = getattr(seg, "_pt_spacing", None)
    if sp is None:
        P = np.asarray(seg.pts, float)
        sp = (float(np.median(np.linalg.norm(np.diff(P, axis=0), axis=1)))
              if len(P) > 1 else 0.0)
        seg._pt_spacing = sp
    return sp


def _wall_distance(seg, pts, dd, ii, stops=(True, True)):
    """How far each point is from where this segment's WALL runs.

    *stops* says whether the start and the end are places the vessel actually stops
    -- a junction -- rather than a free end where the centerline merely runs out.

    A KD-tree over a segment's centerline points CLAMPS: ask it about a point past
    the end of the line and it answers with the distance to the last point, which is
    mostly a distance ALONG the vessel rather than across it. Comparing that to the
    segment's radius -- which is what the wall contest in vessel_coords does -- asks
    whether a face sits on a wall the segment does not have there, because the
    segment stopped.

    That is not hypothetical. On JK1_2 the aorta is split at z=517.6 into Line09
    below and Line06 above. The 86 faces on Line06's inlet stub sit 16 mm ABOVE
    Line09's last centerline point, and every one of them clamps to it: 18.74 mm
    away, against a 15.67 mm radius, for a wall miss of 3.06 -- better than the
    5.11 mm miss against Line06, whose stub they are actually on. So Line09 won a
    quarter of the surface the inlet cut has to remove, the cutter would not remove
    another vessel's wall, and the leftover tongue held the rim 25-35% off plane at
    every offset, plane size and mode. A plain plane cut of that stub is flat to
    0.000 mm: there is nothing above the plane but the aorta's own end, in one
    connected piece.

    So take the component perpendicular to the segment's axis, which is the only
    part that means anything as a radius. Between the endpoints nothing changes --
    the nearest point is interior and the distance is already perpendicular. Off the
    end, 18.74 mm of which 16.4 is axial becomes 9.06 mm across, the miss goes to
    6.6, and the stub stays with the vessel it belongs to.

    Correcting the distance is not enough on its own, and the correction alone makes
    things worse. Collapse the axial part of a big overshoot and what is left is
    small, which FITS a fat vessel: with only the correction in place, JK1_2's stub
    went to Line10, whose centerline stops 36 mm further down still, because 36 mm
    of overshoot leaves 8 mm across and its wall is 15.6 mm. WALL_REACH does not
    catch that either -- it scales with the vessel, and on a 15 mm aorta it reaches
    47 mm.

    So the distance is measured properly AND a claim past the end is refused
    outright. "More past the end than beside it" was tried as the refusal and is too
    weak on a fat vessel: at the default 6 mm the plane sits 11 mm above Line09's
    last point, the aorta is 13 mm across there, and 11 < 13 let a quarter of the
    stub go anyway. The rule that holds is the plain one -- past its endpoint the
    segment is not there, so it may not claim anything.

    The tolerance is the segment's own point spacing, which is the resolution at
    which "where this segment runs" is defined at all, rather than a number chosen
    to make a case pass. It leaves the one claim past an end that is legitimate --
    the collar of a shallow ostium, a fraction of a millimetre behind its own s=0 --
    and removes the ones that are 11 and 16 mm out.

    ONLY AT A JUNCTION. An endpoint means two unrelated things. At a junction the
    vessel really does stop and another one carries on, so a face past it belongs to
    the segment that continues. At a FREE end the vessel does not stop at all -- the
    centerline stops short of the mesh, and the cap past it is this vessel's own
    surface, which it must stay free to claim.

    Measured, the distinction costs and buys nothing: gating both kinds of end and
    gating only junctions give byte-identical cut_bench runs over 4,145 cuts on 79
    cases. It is kept because it is the true statement of when a segment has ended,
    not because a number moved -- past a free end there is only the cap, and the
    nearest-centerline assignment has already given that to the right vessel, so the
    contest never arises there. Do not read the equality as licence to drop it: it
    says this cohort never exercises the case, not that the case cannot happen.

    Both values are returned: the caller needs the honest distance for the fit and
    the overshoot for the refusal.
    """
    n = len(seg.pts)
    over = np.zeros(len(dd))
    if n < 2:
        return dd, over
    out = dd.copy()
    P = np.asarray(seg.pts, float)
    for idx, nb, stop in ((0, 1, stops[0]), (n - 1, n - 2, stops[1])):
        if not stop:                        # a free end: the vessel runs on past it
            continue
        at = np.flatnonzero(ii == idx)      # only the CLAMPED points can be off the end
        if not len(at):
            continue
        tang = P[idx] - P[nb]
        norm = float(np.linalg.norm(tang))
        if norm < 1e-9:
            continue
        t = np.maximum((pts[at] - P[idx]) @ (tang / norm), 0.0)
        over[at] = t
        out[at] = np.sqrt(np.maximum(dd[at] ** 2 - t ** 2, 0.0))
    return out, over


def _segment_radii(mesh, tree):
    """Each segment's median lumen radius, filled in and cached on the segments."""
    out = {}
    for name, seg in tree.segments.items():
        r = getattr(seg, "_radius_mm", None)
        if r is None:
            r = segment_radius(mesh, seg)
            seg._radius_mm = r
        out[name] = max(float(r), 1.0)
    return out


# Rays cast around the centerline to measure the lumen at one station, and the
# percentile of their lengths taken as the radius. A low percentile on purpose: at an
# ostium or a confluence some rays leave through the opening and run the length of the
# parent -- on KY160's arch they return 90 to 100 mm against a wall 20 mm away -- so
# the mean and even the median are dragged upward by whichever direction happens to be
# open. The quartile ignores those and still averages over most of the circumference.
RADIUS_PCTL = 25.0
RADIUS_STEP_MM = 4.0
# A station whose reading is this many times the vessel's own typical one is not
# measuring the vessel. It happens at both ends: the last stretch of a branch's
# centerline lies inside its parent, so every ray leaves through the ostium and
# returns the parent's far wall -- KY160's left subclavian reads 59 mm at its origin
# against 3 to 6 mm along its length. Those stations are dropped and interpolated
# across rather than believed.
RADIUS_OUTLIER = 2.5


def radius_profiles_from_faces(tree, owner, s_val, radial, step_mm=RADIUS_STEP_MM):
    """Radius along every centerline, read off the surface that is already there.

    The faces themselves are the measurement. Each one sits on some vessel's wall and
    the KD query has already said which centerline it is nearest and how far, so
    binning by arc length and taking a percentile per bin gives the lumen radius at
    every station for the cost of a sort. Rays were the obvious way to do this and the
    wrong one: casting them took four seconds a model, which is four seconds of a
    frozen window, and they answer badly at exactly the places that matter, escaping
    through an ostium and returning whatever wall lies beyond it.

    A percentile rather than a mean because a bin near a junction holds a minority of
    faces from the vessel next door, sitting further out. The wall is what most of the
    faces agree on.
    """
    prof = {}
    for name, seg in tree.segments.items():
        L = float(seg.length)
        n = max(int(np.ceil(L / max(step_mm, 0.5))) + 1, 3)
        ss = np.linspace(0.0, L, n)
        m = owner == name
        r = np.full(n, np.nan)
        if m.any():
            idx = np.clip(np.searchsorted(ss, s_val[m]) - 1, 0, n - 1)
            d = radial[m]
            order = np.argsort(idx, kind="stable")
            idx, d = idx[order], d[order]
            starts = np.searchsorted(idx, np.arange(n), side="left")
            ends = np.searchsorted(idx, np.arange(n), side="right")
            for i in range(n):
                if ends[i] - starts[i] >= 3:
                    r[i] = np.percentile(d[starts[i]:ends[i]], RADIUS_PCTL)
        ok = ~np.isnan(r)
        if not ok.any():
            r[:] = 1.0
        else:
            mid = r[ok]
            if len(mid) >= 5:
                core = np.median(mid)
                if np.isfinite(core) and core > 0:
                    ok &= ~(r > RADIUS_OUTLIER * core)
            r = (np.interp(ss, ss[ok], r[ok]) if ok.any()
                 else np.full(n, float(np.nanmedian(r))))
        if n >= 5:
            r = np.convolve(np.pad(r, 1, mode="edge"), np.ones(3) / 3.0, mode="valid")
        prof[name] = (ss, np.maximum(r, 0.3))
    return prof


def radius_at(prof, name, s):
    """The vessel's radius at arc length *s*, for one segment or an array of faces."""
    ss, rr = prof[name]
    return np.interp(s, ss, rr)


# When to doubt a face's owner, and what it takes to overturn it.
#
# A genuine wall face sits about ONE of its own vessel's radii from that vessel's
# centerline: measured over KY160_1's seven segments, the median ratio is 1.0 to 1.2
# and a fat trunk never exceeds 1.5. Past 1.5 the attribution is worth re-checking
# against every segment; below it, checking would only cost time.
#
# Overturning takes a clear margin, because the re-check is itself an estimate: a
# vessel's radius is a median over its middle, so it understates a flared ostium and
# overstates a stenosed mid-section. Requiring the alternative to fit a third better
# keeps the ordinary ostium flare with the vessel it belongs to.
RADIAL_SUSPECT = 1.5
RADIAL_MARGIN = 0.75

# A face is worth re-checking once it sits this far from where its owner's wall runs,
# and a rival must beat that miss by this much to take it. The margin is what keeps
# attribution stable: without it, faces on a shared ostium flip between two vessels on
# a tenth of a millimetre and the boundary between them comes out ragged.
WALL_MISS_MM = 1.5
WALL_WIN_MM = 0.5
# ...and this fraction of the vessel's own radius, whichever is larger. A wall is not
# smooth, and how far it wanders scales with the vessel: 1.5 mm of wobble on a 20 mm
# aorta is a normal wall, on a 3 mm branch it is a different vessel. Judging both by
# the absolute figure alone flagged 26% of every face as worth re-checking, and each
# one is re-queried against every segment -- most of the cost of attribution, spent on
# faces that already sit where they should.
WALL_MISS_FRAC = 0.35
# However well a vessel's wall happens to fit, a face this many radii from its axis is
# not on it. Stops a thin branch three vessels away from winning on a coincidence.
WALL_REACH = 3.0


def vessel_coords(mesh, tree, return_radial=False):
    """For every face: which centerline segment owns it, and at what arc length.

    This is the surface's "vessel coordinate" -- it says *which vessel* a face lies on
    and *how far along* it sits. A cut plane can tell us neither, which is why the
    cutter decides what to remove from this and not from the plane.

    Ownership is by NORMALISED distance: distance to a segment's centerline divided by
    that segment's own radius, not raw nearest point. Raw nearest point fails wherever
    a thin vessel runs alongside a fat one, because a fat trunk's axis lies one trunk
    radius deep in the lumen -- so a patch of trunk WALL can be nearer to a thin
    branch's centerline than to its own.

    That is not a corner case, it is the bug behind the torn openings. Measured on
    KY160_1: 304 mm2 of descending-aorta wall was claimed by a 3.7 mm branch, sitting
    11 to 18 mm out -- three to five of that branch's radii, where a real wall face
    sits at one. The cutter took those faces as part of the branch, deleted them, and
    punched a hole in the aorta; the opening came out 17 mm from its plane on a cut
    32 mm clear of any junction, and every explanation in terms of the cut plane was
    looking in the wrong place.

    Normalising makes the comparison fair. Those faces are 1.0 aorta-radii from the
    aorta's centerline and 3 to 5 branch-radii from the branch's, so the aorta wins on
    the only measure that means anything.

    With *return_radial*, also returns each face's distance to the centerline that
    ends up owning it.
    """
    kd, names, arcs = _centerline_kd(tree)
    centers = np.asarray(mesh.triangles_center)
    # Every core, and there are tens of thousands of faces. This one query is half of
    # what a re-cut costs, and it is embarrassingly parallel.
    d, idx = kd.query(centers, workers=-1)
    owner, s_val, radial = names[idx].copy(), arcs[idx].copy(), d.copy()

    radii = _segment_radii(mesh, tree)
    uniq = np.array(sorted(radii))
    r_owner = np.array([radii[n] for n in uniq])[np.searchsorted(uniq, owner)]

    # Only the doubtful faces are re-checked against every segment. On KY160_1 that is
    # a few per cent of them, which keeps this a single extra KD query on a small
    # subset rather than one per segment over the whole mesh -- the cutter runs this
    # once per opening on every slider release, so the difference is felt.
    # WHOSE WALL IS THIS. The nearest centerline is not the right question and never
    # was: a trunk's axis lies a full radius deep in its own lumen, so a branch leaving
    # it has wall that is genuinely closer to the trunk's axis than to its own. What
    # distinguishes them is not distance but WHERE EACH VESSEL'S WALL IS. A face on the
    # aorta sits about one aortic radius from the aortic axis. A face on a branch by
    # its ostium sits FURTHER out than that -- it is outside the aortic tube, which is
    # what "it is clearly the branch" means when you look at it.
    #
    # So the face goes to whichever vessel it is nearest the WALL of, |d - r(s)|, with
    # r read off that vessel's radius profile at the station in question rather than a
    # single median. On KY160's left subclavian at its flare: the face is 6.5 mm out
    # where the subclavian's wall is 6.4 mm (miss: 0.1) and 23.6 mm out where the
    # aorta's is 20 (miss: 3.6). The subclavian keeps its own wall.
    # Cached against the segments themselves, not the mesh. A cut returns a new mesh
    # but the same vessels, so keying by mesh rebuilt these on every opening; keying by
    # the segment names and lengths still notices region planning splitting a segment
    # underneath us, which is the case that made a mesh-keyed cache unsafe.
    shape = tuple((n, round(float(s.length), 3)) for n, s in tree.segments.items())
    hit = getattr(tree, "_rprof_key", None)
    if hit == shape:
        prof = tree._rprof
    else:
        prof = radius_profiles_from_faces(tree, owner, s_val, radial)
        tree._rprof_key, tree._rprof = shape, prof
    r_local = np.empty(len(centers))
    for name in tree.segments:
        m = owner == name
        if m.any():
            r_local[m] = radius_at(prof, name, s_val[m])
    miss = np.abs(radial - r_local)
    tree._rprof = prof              # the cut reads its band width off this

    nodes = _segment_nodes(tree)    # which ends are junctions, for _wall_distance
    suspect = np.flatnonzero(
        (radial > RADIAL_SUSPECT * r_owner)
        | (miss > np.maximum(WALL_MISS_MM, WALL_MISS_FRAC * r_local)))
    if len(suspect):
        pts = centers[suspect]
        best = radial[suspect] / r_owner[suspect]
        best_miss = miss[suspect]
        # The contest, decided on the wall. A candidate takes the face when it misses
        # by less -- when the face sits closer to where that vessel's surface actually
        # runs. The old ratio test stays as a floor so a vessel cannot win a face that
        # is absurdly far from its axis by having a conveniently sized wall there.
        for name, seg in tree.segments.items():
            # Single-threaded on purpose. These run once per segment over the few per
            # cent of faces worth re-checking, and spawning a thread pool for each was
            # costing more in setup than the query itself: 96 thread starts a re-cut,
            # a quarter of a second of pure overhead.
            dd, ii = _segment_kd(seg).query(pts)
            s_here = np.asarray(seg.s)[ii]
            # Across the vessel, not along it -- a challenger reaching past its own
            # last centerline point is judged on the part of the distance that is
            # actually a radius. See _wall_distance. The reach test below keeps the
            # TRUE distance, so collapsing the perpendicular part cannot buy a face
            # that is simply far away.
            nd = nodes.get(name) or {}
            across, over = _wall_distance(
                seg, pts, dd, ii,
                (nd.get("start") is not None, nd.get("end") is not None))
            miss_here = np.abs(across - radius_at(prof, name, s_here))
            take = ((miss_here < best_miss - WALL_WIN_MM)
                    & (dd < WALL_REACH * radii[name])
                    & (over <= _point_spacing(seg)))   # where this segment RUNS
            if not take.any():
                continue
            sel = suspect[take]
            owner[sel] = name
            s_val[sel] = s_here[take]
            radial[sel] = dd[take]
            best_miss[take] = miss_here[take]

    if return_radial:
        return owner, s_val, radial
    return owner, s_val


def vessel_attribution(mesh, tree):
    """vessel_coords, plus the radius of whichever vessel ended up owning each face.

    The radius is what turns a raw radial distance into a judgement. On its own, "this
    face is 9 mm from the centerline that claims it" says nothing; against a 3 mm
    branch it says the attribution is wrong, and against an 18 mm aorta it says the
    face is exactly where a wall should be.
    """
    # CARRIED FORWARD FROM THE MODEL AS LOADED, not recomputed per mesh.
    #
    # Which vessel a face belongs to is a fact about the face, and cutting does not
    # change it -- a face that survives a cut is the same face on the same wall. But
    # this ran in full on every intermediate mesh, and it is the most expensive thing
    # in a cut: on an 11-opening abdominal case, 108 ms of which 86 is re-deciding
    # 8,683 doubtful faces against all 19 segments, eleven times over, for answers
    # that were settled the first time.
    #
    # So the full result is kept, and a later mesh is attributed by matching each of
    # its faces to the face it came from. Survivors match exactly; a face the plane
    # clipped has moved by a fraction of its own edge and matches its parent. If
    # anything lands further away than that the mesh is not a descendant of this one
    # and the whole thing is recomputed, which is also what happens when region
    # planning splits a segment underneath us -- the stale-cache failure that made
    # two suites accept cuts they exist to refuse.
    shape = tuple((n, round(float(s.length), 3)) for n, s in tree.segments.items())
    hit = getattr(tree, "_attr_full", None)
    if hit is not None and hit["shape"] == shape:
        c = np.asarray(mesh.triangles_center)
        d, i = hit["kd"].query(c, workers=-1)
        if len(d) and float(d.max()) <= hit["tol"]:
            return (hit["owner"][i], hit["s"][i], hit["radial"][i], hit["r_owner"][i])

    owner, s_val, radial = vessel_coords(mesh, tree, return_radial=True)

    # The radius WHERE THE FACE IS, not the vessel's average. Everything that consumes
    # this is asking "is this face about where its vessel's wall should be", and a
    # segment-wide median answers that wrongly at both ends of every branch: KY160's
    # left subclavian is 6 mm at its origin against a 3.69 mm median, so its own wall
    # there reads as nearly two radii out and any cut is entitled to treat it as
    # nobody's. Profiles are already built by the attribution above.
    # Tried reading this off the profile too, so the "nobody owns this convincingly"
    # test would scale with the vessel where the face is. It cost a fifth of a second
    # per re-cut and fixed nothing measurable, so it is not here. The profile is used
    # where it earns its keep: deciding whose wall a face is on.
    radii = _segment_radii(mesh, tree)
    uniq = np.array(sorted(radii))
    r_owner = np.array([radii[n] for n in uniq])[np.searchsorted(uniq, owner)]

    # Keep it, with the tolerance set by the mesh's own triangles: a clipped face's
    # centroid moves by a fraction of an edge, so an edge and a half is generous
    # enough to match every descendant and tight enough to reject a different model.
    centroids = np.asarray(mesh.triangles_center)
    edge = float(np.median(mesh.edges_unique_length)) if len(mesh.faces) else 1.0
    tree._attr_full = {"shape": shape, "kd": cKDTree(centroids), "owner": owner,
                       "s": s_val, "radial": radial, "r_owner": r_owner,
                       "tol": max(1.5 * edge, 1e-6)}
    return owner, s_val, radial, r_owner


def _segment_nodes(tree):
    """For every segment: the junction position at its start / end, or None.

    Cached on the tree. A junction endpoint is a *branch node*, where this
    segment's centerline joins the rest of the tree.
    """
    cache = getattr(tree, "_segnodes_cache", None)
    if cache is None:
        cache = {name: {"start": None, "end": None} for name in tree.segments}
        for jpos, incident in tree.junctions:
            for name, at_end in incident:
                if name in cache:
                    cache[name]["end" if at_end else "start"] = np.asarray(jpos, float)
        tree._segnodes_cache = cache
    return cache


def _terminal_name_for_segment(tree, seg_name, at_end=None):
    """Human name of the opening on *seg_name*.

    When *at_end* is given, prefer the terminal on that specific end, this is
    what an unbranched vessel needs, where BOTH terminals sit on the one segment
    and a plain "first match" always returns the proximal name even when the
    cursor is at the distal end. A real branch has a single terminal, so the end
    hint is harmless there.
    """
    cands = [t for t in tree.terminals if t.segment.name == seg_name]
    if at_end is not None:
        for t in cands:
            if t.at_end == at_end:
                return t.name
    return cands[0].name if cands else None


def measure_along_branch(tree, world_point):
    """Distance along the centerline from *world_point* to its branch node.

    The point (a spot on the mesh surface, e.g. under the cursor) is projected to
    the nearest centerline point; that fixes which vessel it is on and how far
    along. We then return the arc-length distance, in mm, back to the node the
    vessel comes off of, the junction at the near end of its segment.

    For a branch (a segment with one junction end and one free tip) this is
    unambiguous: the only junction on the segment IS its take-off, and walking the
    other way just reaches the tip. For an interior segment (two junction ends) the
    nearer junction is used. For a lone unbranched vessel (no junctions at all) it
    falls back to the nearer vessel end, flagged as such.

    Returns a dict, or None if the tree carries no centerline.
    """
    if not tree.segments:
        return None
    kd, names, arcs = _centerline_kd(tree)

    p = np.asarray(world_point, float)
    radial, idx = kd.query(p)
    seg_name = str(names[idx])
    s = float(arcs[idx])
    seg = tree.segments[seg_name]
    L = seg.length

    nodes = _segment_nodes(tree)[seg_name]
    cands = []                                   # (distance_mm, node_pos, at_end, kind)
    if nodes["start"] is not None:
        cands.append((s, nodes["start"], False, "junction"))
    if nodes["end"] is not None:
        cands.append((L - s, nodes["end"], True, "junction"))
    if not cands:                                # no junction on this segment
        cands.append((s, seg.pts[0], False, "end"))
        cands.append((L - s, seg.pts[-1], True, "end"))

    dist_mm, node_pos, at_end, kind = min(cands, key=lambda c: c[0])

    proj = seg.at_arclength(s)[0]                # projection back onto the centerline
    keep = (seg.s >= s) if at_end else (seg.s <= s)
    mid = seg.pts[keep]
    # Always ordered proj -> node, so a caller can rely on path[0]/path[-1].
    path = np.vstack([proj, mid if at_end else mid[::-1], node_pos])

    # Name the vessel: for a branch, its own opening; for the no-junction fallback,
    # the end actually being measured to (so an unbranched vessel reads "distal"
    # near its distal end, not always "proximal").
    name_hint = at_end if kind == "end" else None
    branch_name = _terminal_name_for_segment(tree, seg_name, name_hint) or seg_name

    return {
        "segment": seg_name,
        "branch_name": branch_name,
        "s_mm": s,
        "seg_length_mm": L,
        "distance_mm": float(dist_mm),
        "proj": np.asarray(proj, float),
        "node_pos": np.asarray(node_pos, float),
        "node_kind": kind,                       # "junction" or "end"
        "radial_mm": float(radial),              # how far the point sat off the centerline
        "path": np.asarray(path, float),
    }


def section_diameter(mesh, point, normal, near_point=None):
    """Cross-section size where a cut plane meets the vessel.

    A vessel cross-section is rarely a circle, so there is no single "diameter".
    The headline number is the EQUIVALENT-AREA DIAMETER, D = 2*sqrt(A/pi): the
    diameter of the circle with the same cross-sectional area. It is the standard
    vascular measure because area is what governs flow, and it is stable regardless
    of how oval the section is. Alongside it we report the min and max caliper
    widths (the shortest and longest the section measures across), which give the
    ovality directly, and the circularity 4*pi*A/P^2 (1.0 = a perfect circle).

    The plane can graze more than one vessel, so the loop nearest *near_point*
    (default: the plane point) is the one measured. Returns a dict, or None if the
    plane misses the mesh.
    """
    point = np.asarray(point, float)
    n = _unit_np(np.asarray(normal, float))
    near = np.asarray(point if near_point is None else near_point, float)

    try:
        sec = mesh.section(plane_origin=point, plane_normal=n)
    except Exception:
        sec = None
    if sec is None:
        return None

    best = None                                  # (dist_to_near, loop_pts)
    for ent in sec.entities:
        loop = np.asarray(sec.vertices[ent.points], float)
        if len(loop) < 4:
            continue
        d = np.linalg.norm(loop.mean(axis=0) - near)
        if best is None or d < best[0]:
            best = (d, loop)
    if best is None:
        return None
    loop = best[1]

    a = _unit_np(np.cross(n, [0.0, 0.0, 1.0]))
    if np.linalg.norm(np.cross(n, [0.0, 0.0, 1.0])) < 1e-6:
        a = _unit_np(np.cross(n, [1.0, 0.0, 0.0]))
    b = _unit_np(np.cross(n, a))
    c = loop.mean(axis=0)
    uv = np.column_stack([(loop - c) @ a, (loop - c) @ b])

    x, y = uv[:, 0], uv[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    perim = float(np.sum(np.linalg.norm(
        np.diff(np.vstack([uv, uv[0]]), axis=0), axis=1)))
    if area < 1e-6:
        return None

    ang = np.linspace(0.0, np.pi, 90, endpoint=False)
    widths = np.array([np.ptp(uv @ np.array([np.cos(t), np.sin(t)])) for t in ang])

    return {
        "area_mm2": float(area),
        "ead_mm": float(2.0 * np.sqrt(area / np.pi)),
        "dmin_mm": float(widths.min()),
        "dmax_mm": float(widths.max()),
        "perimeter_mm": perim,
        "circularity": float(4.0 * np.pi * area / max(perim ** 2, 1e-9)),
        "centroid": np.asarray(c, float),
        "loop": loop,
        "n_loops": len(sec.entities),
    }


def _endpoint_graph(tree):
    """Junction/terminal graph for along-centerline distances. Cached on the tree.

    Returns (node_of, adj) where node_of maps (seg_name, at_end) -> node id and
    adj maps node -> list of (other_node, length_mm, seg_name).
    """
    cache = getattr(tree, "_endgraph_cache", None)
    if cache is not None:
        return cache

    node_of = {}
    for i, (_, incident) in enumerate(tree.junctions):
        for name, at_end in incident:
            node_of[(name, at_end)] = ("J", i)
    for i, t in enumerate(tree.terminals):
        node_of.setdefault((t.segment.name, t.at_end), ("T", i))
    # any endpoint not yet placed (e.g. a lone segment) becomes its own node
    for name in tree.segments:
        for at_end in (False, True):
            node_of.setdefault((name, at_end), ("E", f"{name}:{int(at_end)}"))

    adj = {}
    for name, seg in tree.segments.items():
        a = node_of[(name, False)]
        b = node_of[(name, True)]
        adj.setdefault(a, []).append((b, seg.length, name))
        adj.setdefault(b, []).append((a, seg.length, name))

    cache = (node_of, adj)
    tree._endgraph_cache = cache
    return cache


def _graph_dist(adj, start):
    """Dijkstra distances from *start* over the small endpoint graph."""
    import heapq
    dist = {start: 0.0}
    pq = [(0.0, start)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, np.inf):
            continue
        for v, w, _ in adj.get(u, ()):
            nd = d + w
            if nd < dist.get(v, np.inf):
                dist[v] = nd
                heapq.heappush(pq, (nd, v))
    return dist


def centerline_geodesic(tree, world_a, world_b):
    """Distance ALONG the centerline between two points (each projected to it).

    Same-vessel is the common case (measuring a length along one vessel) and is
    exact: |s_a - s_b|, with the traversed centerline returned for drawing. Across
    a bifurcation it walks the junction graph (shortest path), and the straight
    endpoints are returned for drawing since the full multi-segment path is not
    reconstructed. Returns a dict, or None.
    """
    if not tree.segments:
        return None
    kd, names, arcs = _centerline_kd(tree)
    ia = kd.query(np.asarray(world_a, float))[1]
    ib = kd.query(np.asarray(world_b, float))[1]
    na, sa = str(names[ia]), float(arcs[ia])
    nb, sb = str(names[ib]), float(arcs[ib])
    seg_a, seg_b = tree.segments[na], tree.segments[nb]
    pa, pb = seg_a.at_arclength(sa)[0], seg_b.at_arclength(sb)[0]

    if na == nb:
        lo, hi = sorted((sa, sb))
        keep = (seg_a.s >= lo) & (seg_a.s <= hi)
        path = np.vstack([seg_a.at_arclength(lo)[0], seg_a.pts[keep],
                          seg_a.at_arclength(hi)[0]])
        return {"distance_mm": float(abs(sa - sb)), "same_segment": True,
                "path": path, "proj_a": pa, "proj_b": pb}

    node_of, adj = _endpoint_graph(tree)
    La, Lb = seg_a.length, seg_b.length
    ends_a = {node_of[(na, False)]: sa, node_of[(na, True)]: La - sa}
    ends_b = {node_of[(nb, False)]: sb, node_of[(nb, True)]: Lb - sb}

    best = np.inf
    for ea, ra in ends_a.items():
        dist = _graph_dist(adj, ea)
        for eb, rb in ends_b.items():
            total = ra + dist.get(eb, np.inf) + rb
            if total < best:
                best = total
    if not np.isfinite(best):
        return None
    return {"distance_mm": float(best), "same_segment": False,
            "path": np.vstack([pa, pb]), "proj_a": pa, "proj_b": pb}


def _unit_np(v, eps=1e-12):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / max(n, eps)


def _component(mesh, mask, seeds):
    """Faces of *mask* surface-connected to any face in *seeds* (walking only
    through *mask*)."""
    seeds = np.atleast_1d(seeds)
    if len(seeds) == 0 or not mask.any():
        return np.zeros(len(mesh.faces), dtype=bool)

    fa = np.asarray(mesh.face_adjacency)
    link = fa[mask[fa[:, 0]] & mask[fa[:, 1]]] if len(fa) else np.empty((0, 2), int)

    nf = len(mesh.faces)
    g = coo_matrix((np.ones(len(link)), (link[:, 0], link[:, 1])), shape=(nf, nf))
    _, labels = connected_components(g, directed=False)
    return mask & np.isin(labels, np.unique(labels[seeds]))


def _tears(mesh, delete, untouched, where=None):
    """Whether the removed region has a boundary that is not the cut plane.

    A deleted face touching a face that stays leaves the raw edge between them in the
    result; a clipped face leaves its edge on the plane. So this is exactly the
    question "would this opening come out torn", answered without doing the slice --
    which is what makes it affordable to ask of a dozen candidate offsets.
    """
    adj = np.asarray(mesh.face_adjacency)
    if not len(adj):
        return False
    if where is not None:
        delete = delete & where
    a, b = adj[:, 0], adj[:, 1]
    return bool(((delete[a] & untouched[b]) | (delete[b] & untouched[a])).any())


# Offsets to try when a cut has nowhere to land, as a fraction of the one asked for.
# Only pulling BACK is offered: the vessel gets its own cross-section back by moving
# away from the junction, and a suggestion that moved the cut closer to it would be
# advice to make the same mistake harder.
PULLBACK_TRIES = (0.85, 0.7, 0.55, 0.4, 0.3, 0.2, 0.1)

# ...and how far to push the cut IN, when pulling it back does not help. Absolute
# millimetres rather than fractions: a cut that fails at 6 mm usually needs a few
# more, not a proportion of a vessel whose length has nothing to do with it. Nearest
# step first, so the answer stays as close as possible to what was asked for.
DEEPEN_TRIES = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 11.0, 15.0, 20.0, 26.0)

# Surface a cut may take from a vessel it is not cutting before it stops counting as
# a clean opening. An ostium is fuzzy and a triangle straddling one is nobody's fault;
# a square millimetre is the width of that excuse.
FOREIGN_TOLERANCE_MM2 = 1.0

# Past this, an opening is not an imperfect cut, it is a hole in the wrong shape, and
# handing one back is worse than saying no. The validator already refuses to export it
# and the cut has usually taken a piece of a neighbour to make it, so the surface it
# removed is gone from a vessel nobody asked it to touch.
#
# Measured over 2,927 cuts on 65 aortas: refusing above 3 mm drops 501 of them and
# 55,096 mm2 of foreign surface, and not one of those was ever a clean cut. The bar is
# deliberately well above the 0.5 mm a cut has to meet to count as clean -- between the
# two the opening is flawed but recognisable, the rim is close to its plane, and it is
# the user's call whether to keep it or pull back.
TORN_REFUSE_MM = 3.0

# ...and the bar a rim has to meet to be CLEAN, which is the other end of that band.
# Between the two the cut is handed back, deliberately -- but until now it was handed
# back SILENTLY, and that is the part that was wrong. A cut at 2.15 mm has taken a
# piece of a neighbouring vessel (31 mm2 of it on JK24_1's Branch 4) and looks like
# any other opening in the list, while a clean cut sits 0.8 mm further back and
# suggest_pullback can find it. Refusing is policy and stays the user's call; not
# MENTIONING it is a defect. Same value cut_bench scores against, so a cut this
# module calls flawed is one the benchmark would not count as clean.
RIM_FLAW_MM = 0.5


# Fractions of the asked-for offset to try when a MODEL END will not come out flat.
# Descending, so the answer stays as close as possible to what was asked for and as
# much of the vessel is kept as the geometry allows.
TRUNK_BACKOFF = (0.75, 0.5, 0.35, 0.25, 0.15)


def is_model_end(tree, terminal):
    """Is this a free end of the TRUNK -- the end of the model itself?

    True only where nothing lies past the plane but this vessel's own stub: a trunk
    terminal, on an end that is not a junction, and not one of the three things that
    look like a model end and are not -- a trim (which is meant to swallow the
    openings below it), a region boundary (a plane inside the vessel), or an interior
    cut.
    """
    if not getattr(terminal, "is_trunk", False):
        return False                      # a branch TIP is a free end too, not this
    if getattr(terminal, "is_interior", False):
        return False
    if getattr(terminal, "subtree_segments", None):
        return False
    if getattr(terminal, "_region_offset", None) is not None:
        return False
    nd = _segment_nodes(tree).get(terminal.segment.name) or {}
    return nd.get("end" if terminal.at_end else "start") is None


def flat_trunk_offset(mesh, tree, terminal, want, snap_mm=0.0, attr=None):
    """The largest offset at or below *want* giving ONE flat opening, or None.

    The offset at a vessel's FREE END is arbitrary in a way an interior plane's is
    not: it says how much of the end stub to trim off, and a shallower cut simply
    keeps more vessel. A region boundary is the opposite -- it is a station chosen on
    the aorta, and moving it changes which stretch was measured -- so bounds are never
    sent here; they have clear_bound_offset, which moves them ALONG the vessel and
    reports the distance.

    Used for the aorta's own end and, when a count failure points at one, for a branch
    tip: JK30_4 refuses with "10 openings, expected 9" because Branch 4's plane runs
    through a neighbour, and 4.5 mm instead of 6.0 fixes it. The name says trunk
    because that is what it was written for.

    Four abdominal inlets refused at the default 6 mm and cut EXACTLY flat a few
    millimetres shallower, the deviation growing monotonically with depth:

        offset      1     2     3     4     5   6 (default)  7
        JK14_3   0.00  0.00  0.00  0.98  1.15    1.40       2.38
        JK23_3   0.00  0.00  0.00  0.70  1.26    2.25        --
        JK30_3   0.00  0.00  0.00  0.91  1.73    2.54        --
        JK6_2    0.00  0.00  0.00  0.00  0.00    0.00       1.02

    The cause is geometric and no better choice of faces fixes it. These inlet
    segments are 15 to 22 mm long, and the deeper the plane travels the more likely
    it crosses a NEIGHBOURING branch that runs upward. That branch's collar borders
    the stub, so removing it cuts into the branch and leaving it puts its above-plane
    corners into the rim -- one to six faces, 0.3 to 1.6 mm proud, which is exactly
    the deviation reported. Four separate face-selection rules were built against
    that and every one traded one case for another; the worst raised foreign surface
    from 19,375 to 840,623 mm2. See tools/ABDOMINAL.md.

    suggest_pullback searches the same way and is deliberately NOT used here: it also
    requires the cut to take under a square millimetre from any other vessel, and at
    these shallow offsets attribution still hands part of the stub to the neighbour,
    so it rejects offsets whose rim measures 0.000 mm. The bar that decides whether a
    user can export is the one validate_cut applies, so that is the bar used.
    """
    attr = attr or vessel_attribution(mesh, tree)
    for frac in TRUNK_BACKOFF:
        off = max(float(want) * frac, MIN_OFFSET_MM)
        try:
            out = cut_terminal(mesh, tree, terminal, off, snap_mm=snap_mm,
                               _attr=attr, _suggest=False)
        except ValueError:
            if off <= MIN_OFFSET_MM:
                break
            continue
        loops = boundary_loops(out)
        # ONE opening, and flat. Flatness alone was the first test and it misses the
        # other way this cut fails: on JK14_2 the inlet plane runs straight through
        # the CELIAC TRUNK -- 109 of its faces sit a median 1.01 mm off the plane --
        # so the cut opens the celiac as well as the aorta and comes back with two
        # openings where one was asked for. That is the "12 openings, expected 11"
        # that case refuses with on the plain path. Both are the same fault, a plane
        # crossing a neighbouring vessel, and both answer to the same remedy.
        if len(loops) == 1:
            pt, nrm = terminal.cut_plane(off)
            v = np.asarray(out.vertices)
            rim = min(loops, key=lambda L: float(
                np.linalg.norm(v[np.asarray(L)].mean(axis=0) - pt)))
            if float(np.abs((v[np.asarray(rim)] - pt) @ nrm).max()) <= RIM_FLATNESS_MM:
                return off
        if off <= MIN_OFFSET_MM:
            break
    return None


# How far a region boundary may be nudged to clear the ostium it is undercutting,
# and in what order. Nearest first, so the region stays as close to what was asked
# for as the geometry allows, and BOTH WAYS -- outward only was the first version and
# it left rescuable bounds refusing, because which side is clear depends on which
# side the branch is on. Measured over the eight that would not cut at all:
#
#     JK25_3 aortic neck    -6 mm   rim 0.00      (+1 to +2 refuse outright)
#     JK25_3 visceral       +1 mm   rim 0.10
#     JK6_1  visceral       -8 mm   rim 0.00      (every positive step refuses)
#
# Eight millimetres is the limit and it is not a small move on a 40 mm region, which
# is why it is reported rather than absorbed. The alternative on these is refusing
# the region entirely.
# 0.0 comes FIRST so a boundary that already cuts cleanly is left exactly where it
# was planned. Without it the search moves every bound it is handed, including ones
# with nothing wrong -- JK30_5's region end was being shifted a millimetre for no
# reason at all.
BOUND_NUDGE_MM = (0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0,
                  5.0, -5.0, 6.0, -6.0, 7.0, -7.0, 8.0, -8.0)


def clear_bound_offset(mesh, tree, terminal, want, snap_mm=0.0, attr=None):
    """Nudge a region boundary until its cut leaves ONE opening, or None.

    A boundary sitting a millimetre or two from a branch's ostium undercuts it: the
    removed side takes aortic wall from around the base of a vessel the region KEEPS,
    and the hole that leaves is an opening nobody asked for. It is made entirely of
    faces the cut legitimately owns, which is why no ownership, band or subtree rule
    touches it -- measured, the cut removes zero foreign surface and still comes back
    with two openings where one was expected.

    What it looks like on JK17_2's infrarenal bound: the aortic cross-section is
    163 points, 15.46 mm across and 0.00 mm off its plane -- flawless -- beside a
    19-point hole 1.79 mm across, 35 mm away. That second loop is a branch's ostium.

    Openings left by the bound's own cut, against the nudge applied:

        shift      0    +1   +2   +3   +5
        JK17_2     2     1    1    1    1
        JK30_5     2     2    1    1    1
        JK4_2      2     2    2    1    1

    One to three millimetres is enough in every case. Probed with a single cut per
    step rather than a full re-cut, because the count from the bound alone is the
    whole signal.

    Ruled out first, and recorded so they are not re-tried: the plane's reach (swept
    55 mm down to 6 mm, no effect on any of six presets), ownership of the faces
    taken (zero foreign), the subtree being removed in part (taking it whole changes
    nothing), and the subtree wrongly containing a branch that leaves inside the
    region (none does).
    """
    attr = attr or vessel_attribution(mesh, tree)
    L = float(getattr(terminal.segment, "length", 0.0) or 0.0)
    for d in BOUND_NUDGE_MM:
        off = float(want) + d
        if off > 0.95 * L or off < MIN_OFFSET_MM:
            continue
        try:
            out = cut_terminal(mesh, tree, terminal, off, snap_mm=snap_mm,
                               _attr=attr, _suggest=False)
        except ValueError:
            continue
        loops = boundary_loops(out)
        if len(loops) != 1:
            continue
        # Flat as well as single. The opening count alone was the first test and it
        # is not the bar the tool applies: validate_cut refuses a rim over
        # RIM_FLATNESS_MM whatever its count, so a nudge judged on count can hand
        # back an offset the very next re-cut rejects. Same mistake as scoring a
        # boundary fitter on opening count, recorded in tools/LOOP.md rule 7.
        pt, nrm = terminal.cut_plane(off)
        v = np.asarray(out.vertices)
        rim = min(loops, key=lambda Lp: float(
            np.linalg.norm(v[np.asarray(Lp)].mean(axis=0) - pt)))
        if float(np.abs((v[np.asarray(rim)] - pt) @ nrm).max()) <= RIM_FLATNESS_MM:
            return off
    return None


def suggest_pullback(mesh, tree, terminal, offset_mm, band_scale=2.5):
    """The largest offset below the one asked for that still has a clean cross-section.

    Returned as (offset_mm, fraction_of_vessel) or None. "Pull it back" is the right
    advice and a useless amount of it; this is the same advice with a number on it.
    The candidates are probed with the real rule and judged by the tear test, so a
    suggestion is one the cutter can actually honour, not an estimate.
    """
    L = float(getattr(terminal.segment, "length", 0.0) or 0.0)
    attr = vessel_attribution(mesh, tree)   # computed once, handed to every probe
    owner = attr[0]
    own = owner == terminal.segment.name
    sub = getattr(terminal, "subtree_segments", None)
    if sub:
        own = own | np.isin(owner, list(sub))

    def works(off):
        if off < 0.5:
            return False
        try:
            band, delete, untouched = cut_terminal(
                mesh, tree, terminal, off, band_scale=band_scale,
                _probe=True, _attr=attr)
        except ValueError:
            return False
        if not band.any():
            return False
        # Flat is not enough to recommend: an offset that closes its rim by taking a
        # piece of the parent is the failure this whole exercise is about. Held to the
        # same square millimetre a cut is, and checked first because it is cheap and
        # rejects without doing the slice.
        if float(mesh.area_faces[(delete | band) & ~own].sum()) > FOREIGN_TOLERANCE_MM2:
            return False

        # Then measure the opening this offset would actually make. The tear test was
        # used here and is too blunt for the job: it asks whether ANY deleted face
        # touches one that stays, which is true wherever a vessel touches its
        # neighbour down its length, whatever this cut does. On KY160's left
        # subclavian it rejected offsets whose rim comes out at 0.00 mm and walked the
        # answer back 40 mm past cuts that were already perfect. The rim is the
        # question, so ask the rim.
        try:
            out = cut_terminal(mesh, tree, terminal, off, band_scale=band_scale,
                               snap_mm=snap_tolerance(mesh), _attr=attr,
                               _suggest=False)   # never search from inside a search
        except ValueError:
            return False
        loops = boundary_loops(out)
        if not loops or int(out.euler_number) != 1:
            return False
        pt, nrm = terminal.cut_plane(off)
        verts = np.asarray(out.vertices)
        rim = min(loops, key=lambda L_:
                  float(np.linalg.norm(verts[np.asarray(L_)].mean(axis=0) - pt)))
        return float(np.abs((verts[np.asarray(rim)] - pt) @ nrm).max()) <= RIM_FLATNESS_MM

    if works(offset_mm):
        # Asked from a cut that is already fine. Say so rather than march it back to
        # a "safer" offset it never needed -- the whole point is to keep as much
        # vessel as the geometry allows.
        return offset_mm, (offset_mm / L if L > 0 else 0.0)

    # DEEPER FIRST, then shallower. The search only ever went back toward the tip,
    # and the message it produced when that failed -- "no offset along this vessel
    # has a clean cross-section" -- was simply untrue whenever a deeper one worked.
    # JK1_2's inlet refuses at the default 6 mm and cuts perfectly at 10; the tool
    # told the user nothing would work while sitting 4 mm from an answer.
    #
    # It is tried in step order so the nearest working offset wins, which keeps the
    # cut as close as possible to where it was asked for. Bounded at 90% of the
    # segment: past that the plane is at the far junction and has no vessel left to
    # cross. Measured over the 16 abdominal cases the plain path refuses, searching
    # both ways rescues 4 that searching one way does not.
    for step in DEEPEN_TRIES:
        off = offset_mm + step
        if off > 0.9 * L:
            break
        if works(off):
            return float(off), (off / L if L > 0 else 0.0)

    lo, hi = None, offset_mm      # hi is the nearest offset known NOT to work
    for frac in PULLBACK_TRIES:
        off = offset_mm * frac
        if off < 0.5:
            break
        if works(off):
            lo = off
            break
        hi = off
    if lo is None:
        return None

    # Close in on the largest offset that still works. The ladder lands wherever its
    # steps happen to fall, and the difference is not cosmetic: on KY160's left
    # subclavian the ladder's first hit is 146 mm when 160 mm cuts perfectly, so the
    # advice would be to give up 14 mm of vessel for no reason. Whoever is cutting
    # wants to keep as much of it as they can.
    for _ in range(6):
        if hi - lo < 1.0:
            break
        mid = 0.5 * (lo + hi)
        if works(mid):
            lo = mid
        else:
            hi = mid
    return lo, (lo / L if L > 0 else 0.0)


def _no_offset_reason(mesh, tree, terminal, band_scale=2.5):
    """Why NO offset works here -- the rim, or the neighbour, or a real mismatch.

    Three different faults present as "nothing cuts", and only one of them is
    anything the person at the screen can act on:

      * every plane spills into another vessel. The vessel is too short for the
        reach its own calibre needs, so a plane wide enough to span the ring also
        crosses the junction. Nothing on this row fixes it.
      * every plane tears. The surface here has no planar cross-section at all.
      * the centerline is outside the mesh. A different fault entirely, and
        preflight already says so, but it lands here too.
    """
    seg = terminal.segment
    L = float(getattr(seg, "length", 0.0) or 0.0)
    r = float(getattr(seg, "_radius_mm", None) or 0.0)
    reach = band_width(terminal, band_scale)
    if r > 0 and L > 0 and reach > L:
        return (f"No offset works: this vessel is {L:.0f} mm long and needs a "
                f"{reach:.0f} mm plane to span its own {2 * r:.0f} mm width, so every "
                f"plane wide enough to cut it also reaches past its far end into the "
                f"vessel beyond. That is the segmentation, not a setting -- trim the "
                f"neighbour first, or cut this end on a longer stretch of vessel.")
    return ("No offset along this vessel has a clean cross-section: shallower planes "
            "tear, and deeper ones cut into the vessel beyond it. This usually means "
            "the mesh and the centerline disagree here.")


def _with_suggestion(msg, mesh, tree, terminal, offset_mm, band_scale, probing):
    """Append a workable offset to a refusal, when one can be found."""
    if probing:
        return msg                      # never search from inside a search
    found = suggest_pullback(mesh, tree, terminal, offset_mm, band_scale)
    if found is None:
        # Say WHICH way it is impossible. "No offset has a clean cross-section" sent
        # people hunting for a slider position that does not exist, when the usual
        # cause is the opposite: an offset deep enough to give a flat rim reaches past
        # this vessel and into the next one. JK1_2's inlet is 17.8 mm of segment on an
        # 8.6 mm radius, so the plane needs 21.5 mm of reach to span the ring -- more
        # than the whole segment -- and every plane that spans it lands on a
        # neighbouring branch. No knob fixes that, and the message should not imply
        # one might.
        return msg + " " + _no_offset_reason(mesh, tree, terminal, band_scale)
    off, frac = found
    return msg + (f" Cutting at {off:.0f} mm instead of {offset_mm:.0f} mm "
                  f"({frac * 100:.0f}% along this vessel) does have one.")


def cut_terminal(mesh, tree, terminal, offset_mm, band_scale=2.5, snap_mm=0.0,
                 _probe=False, _attr=None, _suggest=True):
    """Open one terminal: remove this vessel past the cut, leaving a planar face.

    Three roles, deliberately separated -- conflating them is what broke earlier
    versions:

      ARC LENGTH decides what is REMOVED. A cut plane is an infinite half-space
        and knows nothing about the vessel; a tortuous branch re-enters its own
        half-space (Line5's centerline crosses its own plane three times at large
        offsets), so "beyond the plane" stops meaning "past the cut". Arc length
        is monotone along the vessel by construction and never lies.

      A LOCAL BAND around the cut point decides what gets CLIPPED FLAT. Only
        faces near the cut need the plane; restricting to them means the infinite
        plane cannot touch anything far away that it happens to slice.

      OWNERSHIP keeps the removed core on this vessel, so a cut cannot chew
        through an ostium into the parent. Ownership is NOT applied to the band:
        near a junction the cut ring is partly owned by the parent's centerline,
        and excluding it left the ring unclipped (the aortic root came out 9%
        undersized).

    Everything distal by arc length but outside the band is deleted outright.
    That is what kills the floating scraps a wandering branch leaves behind when
    it dips back across its own cut plane.
    """
    if len(mesh.faces) == 0:
        raise ValueError("mesh has no faces left to cut")

    s_cut = terminal.s_cut(offset_mm)
    # The band centre follows the in-plane shift; the slice itself is unchanged
    # (the shift lies in the plane), so the opening still lands at the true point.
    point, outward = terminal.cut_anchor(offset_mm)
    seg = terminal.segment

    if getattr(seg, "_radius_mm", None) is None:    # cached: rebuilding a proximity
        seg._radius_mm = segment_radius(mesh, seg)  # query per cut is far too slow

    # How far from the cut point the rim may reach. This is the "plane size" the UI
    # exposes, and with the band defined below it is a SAFETY BOUND rather than the
    # selection itself: the automatic value sits well clear of any legitimate ring
    # (2.5 radii against a ring spanning r/cos(obliquity), so it does not bind below
    # 66 degrees), and shrinking it by hand is the escape hatch for a cut sitting so
    # close to a junction that the plane's intersection curve has merged with the
    # parent's.
    owner, s_face, radial, r_owner = _attr or vessel_attribution(mesh, tree)

    # Sizing the band on the LOCAL radius was tried here and made deep cuts worse, not
    # better. Near an ostium the vessel is at its widest and its parent is at its
    # closest, so the wider band reaches the parent before it finishes spanning the
    # ring: KY160's left subclavian at 178 mm went from a rim torn 11.9 mm out of plane
    # to a cut that left the vessel altogether. No single radius works there, because
    # the ring and the parent are the same distance away. Left on the median.
    band_w = band_width(terminal, band_scale)
    centers = np.asarray(mesh.triangles_center)
    verts = np.asarray(mesh.vertices, float)
    tip = np.asarray(terminal.tip, float)

    # --- core: everything the cut removes -----------------------------------
    distal = (owner == seg.name) & terminal.distal_mask(s_face, s_cut, 0.0)
    # A pre-bifurcation cut removes the whole downstream subtree as well.
    if getattr(terminal, "subtree_segments", None):
        distal = distal | np.isin(owner, list(terminal.subtree_segments))

    # TWO ATTEMPTS TO FIX THE AORTIC INLET LIVED HERE AND BOTH FAILED cut_bench.
    #
    # The fault is real and measured: on twelve of the 49 abdominal scans the
    # segment holding the aortic inlet has a median radius 1.88x smaller than its
    # local calibre, attribution hands about a quarter of the surface near the tip
    # to the next AORTIC segment, that surface stays above the plane, and the rim
    # comes out 25-35% ragged at every offset, plane size and mode. Median rim
    # deviation is 0.000: most of the rim is perfect and a few arcs are stuck.
    #
    # Taking those faces anyway does fix it, and costs more than it is worth:
    #   every terminal, beyond the plane and inside the band
    #       foreign 14,275 -> 121,942 mm2, CLEAN 2162 -> 2092, 71 cuts broken
    #   trunk ends only, aortic segments only, same bounds
    #       foreign 14,275 -> 54,122 mm2, CLEAN 2162 -> 2099, 64 cuts broken
    # for one extra case on the plain path either way.
    #
    # A half-space is the wrong guard on a curved vessel: an aorta bends back over
    # its own cut plane, so "beyond the plane and within the band" reaches around
    # the arch and takes wall that other cuts need. The next thing to try is
    # CONNECTIVITY -- only faces contiguous with the piece actually being removed --
    # which cannot jump a gap the way a half-space does. See tools/ABDOMINAL.md.
    if not distal.any():
        # Nothing beyond the plane to remove. Usually the centerline simply runs
        # past the end of the mesh -- a very common export mismatch -- so the cut
        # lands in empty space BEFORE the geometry begins. Saying how far in the
        # surface actually starts turns a dead end into a number the user can dial.
        own = owner == seg.name
        hint = ""
        if own.any():
            near, far = float(s_face[own].min()), float(s_face[own].max())
            gap = (seg.length - far) if terminal.at_end else near
            if gap > 0.5:
                edge = "tip" if terminal.at_end else "start"
                hint = (f", the centerline runs about {gap:.0f} mm past the end of "
                        f"the mesh at this {edge}, so a {offset_mm:.0f} mm cut falls "
                        f"outside the surface entirely. Try at least "
                        f"{gap + 2.0:.0f} mm.")
        raise ValueError(f"no distal region for {terminal.name}{hint}")
    dcand = np.flatnonzero(distal)
    tip_face = int(dcand[np.argmin(np.linalg.norm(centers[dcand] - tip, axis=1))])
    core = _component(mesh, distal, [tip_face])
    if not core.any():
        raise ValueError(f"empty core for {terminal.name}")

    # --- what the cut actually removes, and where the plane gets to decide ----
    #
    # This used to be a BALL of radius band_w around the cut point: everything inside
    # it was handed to the plane, everything outside was judged by arc length alone.
    # The ball was wrong in both directions at once, and the two errors hid each
    # other.
    #
    #   Too generous. Faces owned by ANOTHER segment had to be admitted on distance
    #     alone -- near a junction part of the ring is claimed by the parent's
    #     centerline, and excluding it left the ring unclipped (the aortic root came
    #     out 9% undersized). But that exemption also put the parent's far wall and
    #     the siblings around the junction inside the ball, and the infinite plane
    #     sliced those too. Measured on KY160_1: a cut 95% down the descending aorta
    #     deleted 5,693 mm2 of two OTHER vessels, and the hole that left merged with
    #     the opening into a rim 30 mm out of plane.
    #
    #   Too mean. An arc-length window inside the ball kept it from swallowing a
    #     stretch that doubles back, but that window also cut through the ring itself
    #     near a junction, so part of the boundary was never clipped and stayed as
    #     raw triangle edges. That is what "torn, not sliced" actually was: not the
    #     plane tearing anything, but the plane never being applied there.
    #
    # Neither is a tuning problem, so this drops the ball and states the thing
    # directly. The removed region is ONE CONNECTED PIECE, reached from the tip, of
    # everything that is
    #
    #     not entirely on the keep side of the PLANE          (so the rim is flat)
    #     and either past the cut by ARC LENGTH               (so "past" survives
    #                                                          a tortuous vessel)
    #     or near the cut and on a face this cut may touch    (so the rim can close
    #                                                          across the ostium)
    #
    # and the faces the plane straddles on its boundary are the rim. All three parts
    # are load bearing. Arc length is what makes "past the cut" mean anything on a
    # vessel that re-enters its own half-space. The plane is what makes the boundary
    # flat -- no set defined by ownership has a flat boundary, and a boundary that
    # follows triangle edges is exactly what a torn opening is. And the walk being
    # confined to faces this cut may touch is what stops it crossing to the vessel
    # next door.
    #
    # A pre-bifurcation cut needs no branch of its own here. It used to take
    # `(core | own) & ball` with no component walk, because letting the band swallow
    # the whole core left `delete` empty and every subtree face proximal to the plane
    # survived -- still attached, still capped, so it added no boundary loop and
    # nothing could see it, and a 40 mm blind-ended vessel shipped while the status
    # bar said "Ready". Here `delete` is the removed region minus its own rim, which
    # is only empty when there is genuinely nothing to remove.
    d_vert = (verts - point) @ outward
    fd = d_vert[mesh.faces]
    keep_side = fd.max(axis=1) <= 0.0          # entirely on the surviving side
    straddle = (fd.min(axis=1) < 0.0) & (fd.max(axis=1) > 0.0)
    near_cut = np.linalg.norm(centers - point, axis=1) <= band_w

    # Two things stop the walk running away, and it needs both.
    #
    # WHOSE VESSEL. Left free, the walk leaves through the ostium, crosses the strip
    # of aortic wall between two ostia, and climbs the vessel NEXT DOOR. That is what
    # was tearing openings nowhere near a junction: on JK5_1 the rim of branch 9's cut
    # came back with 23 of its 98 points sitting on branch 10, up to 8.9 mm off the
    # plane, and branch 10's cut returned the favour with 22 points on branch 9. So a
    # face may be walked through only if this vessel's centerline claims it, or if
    # nobody's does convincingly -- the ostium fuzz, where a wall face genuinely sits
    # several radii from every centerline and attribution is a coin toss. A face
    # sitting one radius from ANOTHER vessel's centerline is that vessel's wall, and
    # this cut has no business there.
    #
    # HOW FAR. Beyond band_w only the core goes, so ambiguity near the cut cannot
    # become a licence to travel. A cut that is going to work never reaches this
    # bound; it binds where the cut was already doomed, and turns a runaway into a
    # torn rim the validator can name.
    ambiguous = radial > RADIAL_SUSPECT * r_owner
    mine = (owner == seg.name) | ambiguous
    if getattr(terminal, "subtree_segments", None):
        mine = mine | np.isin(owner, list(terminal.subtree_segments))
    # SEGMENTS THAT ARE THE SAME VESSEL, named by whoever built this cut.
    #
    # RADIAL_SUSPECT cuts through the middle of a DILATED vessel's wall. A normal
    # aorta's wall sits about one radius from its own centerline and the segment owns
    # it outright; an aneurysmal one drifts out to around 1.5, which is the threshold
    # itself -- so part of the ring reads as nobody's and is cut, and the rest is
    # confidently assigned to the next centerline along and refused. One ring, split
    # by a constant running through it, which is why the tear is a contiguous sector
    # and why no plane size changes it.
    #
    # Measured on KP20_1, a 53 mm ascending aorta: the sliced two-thirds of the rim
    # sit at radial/r_owner 1.67 -- over the line, so "nobody's" -- and the torn third
    # at 1.48, under it, owned by Line1. Line1 is the NEXT SEGMENT OF THE SAME AORTA.
    # The bound was declining to cut the aorta because a different piece of the
    # aorta's own centerline claimed it.
    #
    # So the fix is not a bigger threshold, which would loosen every cut in the tool
    # to rescue two. It is to say what the caller already knows: for a region
    # boundary, every segment of the aortic route is one vessel, and the split into
    # Line4 / Line5_B is centerline bookkeeping rather than anatomy. Kept separate
    # from subtree_segments because that means "the group this trim removes", which a
    # boundary does not have.
    if getattr(terminal, "_own_segments", None):
        mine = mine | np.isin(owner, list(terminal._own_segments))

    # Refuse rather than approximate. Both messages below are the cut declining to
    # produce something it cannot make flat, which is the whole point: the alternative
    # is a cut that reaches onto a neighbour to close its rim and reports success.
    seeds = np.flatnonzero(core & ~keep_side)
    if len(seeds) == 0:
        raise ValueError(_with_suggestion(
            f"the cut plane for {terminal.name} does not cross this vessel: every "
            f"part of it past the cut sits on the keep side of the plane. The plane "
            f"is lying ALONG the vessel rather than across it. Switch this opening "
            f"to Centerline-normal, or reduce its tilt.",
            mesh, tree, terminal, offset_mm, band_scale, _probe))
    # EVERY core face is a starting point, not just the one by the tip.
    #
    # Growing the removed region from a single seed assumes the vessel past the cut is
    # one piece once the keep side is taken away, and near an ostium it is not: the
    # vessel curves, a stretch of it dips back across its own plane, and the surface
    # beyond that dip is cut off from the tip. Those faces are this vessel's own, sit
    # inside the band, on the removed side, and were still left behind -- their
    # neighbours were deleted around them, so the boundary between the two became raw
    # triangle edge instead of the plane. On KY160's left subclavian that is the whole
    # difference between a cut at 36 mm from the junction, which comes out exactly
    # flat, and one at 35 mm, which strays 12.4 mm: 25 stranded faces of its own wall.
    #
    # Seeding from the whole core costs nothing -- it is the same flood fill from more
    # starting points -- and cannot reach further than the candidate set already
    # allowed, which is still gated by ownership.
    # Five cuts in 3,371 go from clean to not-clean here, all thin branches (2 to
    # 4 mm) taking 8 to 46 mm2 of a neighbour or shedding an extra piece. Restricting
    # what the extra seeds may add to this vessel's own wall was tried and changes
    # none of them: the surface they take is admitted by OWNERSHIP, upstream of the
    # walk, not by the seeding. Left alone rather than carrying a guard that costs a
    # component pass and fixes nothing.
    region = _component(mesh, ~keep_side & (core | (mine & near_cut)),
                        seeds.tolist())

    # ARC LENGTH's verdict, for what the plane cannot see. A plane is an infinite
    # half-space and a tortuous vessel re-enters its own: Line5's centerline crosses
    # its own cut plane three times at large offsets, so past the second crossing
    # "beyond the plane" stops meaning "past the cut". Away from the cut, arc length
    # is the only thing that still knows, and it takes the vessel whichever side of
    # the plane it happens to be lying on.
    region = region | (core & ~near_cut)

    band = region & straddle       # the rim: these get clipped flat
    delete = region & ~straddle    # wholly past the plane: these just go
    if not band.any():
        # The region has no rim of its own. Every face the plane crosses here belongs
        # to a vessel this cut is not cutting, which means the cut has run into the
        # junction: there is no cross-section of THIS vessel left to cut.
        gap = s_cut if terminal.at_end else (seg.length - s_cut)
        raise ValueError(_with_suggestion(
            f"{terminal.name} has no cross-section left to cut here. The cut sits "
            f"{gap:.0f} mm from the junction this vessel springs from, close enough "
            f"that the only surface the plane meets belongs to its parent. Pull it "
            f"back.", mesh, tree, terminal, offset_mm, band_scale, _probe))

    untouched = ~(delete | band)
    if _probe:
        return band, delete, untouched

    if snap_mm > 0:
        # Prevent slivers rather than repair them. A sliver appears when the plane
        # grazes a triangle a hair from one of its corners, spawning a rim vertex
        # almost on top of an existing one: the leftover triangle is a NEEDLE --
        # 1-3 mm long, razor thin, and with no short edge, so welding coincident
        # points cannot remove it. Nudging the offending vertex ONTO the plane
        # first means the clip lands on it cleanly and no needle is ever created.
        # Only band vertices are eligible, so a far-off vessel that the infinite
        # plane happens to pass through is never touched.
        verts = verts.copy()
        bv = np.unique(mesh.faces[band])
        d = (verts[bv] - point) @ outward
        sel = np.abs(d) < float(snap_mm)
        verts[bv[sel]] -= np.outer(d[sel], outward)

    kept_v, kept_f, _ = slice_faces_plane(
        verts,
        np.asarray(mesh.faces[band], np.int64),
        plane_normal=-np.asarray(outward, float),   # keeps the proximal side
        plane_origin=np.asarray(point, float),
    )

    rest = trimesh.Trimesh(verts, mesh.faces[untouched], process=False)
    rest.remove_unreferenced_vertices()
    clipped = trimesh.Trimesh(kept_v, kept_f, process=False)

    out = trimesh.util.concatenate([rest, clipped])
    out.merge_vertices()                 # stitches the clipped strip back on
    out.remove_unreferenced_vertices()
    out = keep_largest_component(out)

    # LOOK AT WHAT WAS MADE. Every rule above chooses faces, and none of them can see
    # the opening that results -- which is the only thing the caller asked for. So
    # measure it, and refuse to hand back a hole in the wrong shape rather than leave
    # that for someone downstream to notice.
    #
    # Gated on the cheap test first, because measuring costs 30 ms a cut and most cuts
    # are fine. If no deleted face touches one that stays, the opening's whole boundary
    # came from clipping and every point of it lies on the plane by construction -- an
    # exact answer, not an approximation, so nothing is missed by not looking. The test
    # errs the other way sometimes, where a vessel touches its neighbour somewhere down
    # its length, and that costs only the measurement it triggers.
    if not _tears(mesh, delete, untouched):
        # Exact, not assumed: the whole boundary came from clipping, so every point
        # of it is on the plane. Recorded rather than left unset so a reader can tell
        # "measured and flat" from "never measured".
        terminal._rim_dev_mm = 0.0
        return out

    pt, nrm = terminal.cut_plane(offset_mm)
    loops = boundary_loops(out)
    if loops:
        vo = np.asarray(out.vertices)
        rim = min(loops, key=lambda L_:
                  float(np.linalg.norm(vo[np.asarray(L_)].mean(axis=0) - pt)))
        dev = float(np.abs((vo[np.asarray(rim)] - pt) @ nrm).max())
        # Keep the measurement whatever it says. It was computed and thrown away
        # unless it crossed the refusal bar, which left the middle of the band --
        # cut handed back, rim visibly off its plane -- with nothing for a caller to
        # report. See RIM_FLAW_MM.
        terminal._rim_dev_mm = dev
        if dev > TORN_REFUSE_MM:
            raise ValueError(_with_suggestion(
                # No name here: every caller prefixes one, and saying it twice reads
                # like a stutter in the status bar.
                f"cannot be cut here, the opening comes out {dev:.0f} mm out of "
                f"plane, torn rather than sliced, because this vessel has no clean "
                f"cross-section left at this depth.",
                mesh, tree, terminal, offset_mm, band_scale,
                _probe or not _suggest))
    return out


def flawed_openings(tree, terminals=None):
    """Openings whose rim came out off its plane, worst first.

    [(terminal, dev_mm)] for every cut rim above RIM_FLAW_MM. Empty when every
    opening is clean, and empty before anything has been cut.

    This is the middle of the band described at TORN_REFUSE_MM: above 3 mm the cut is
    refused and the message carries a workable offset, below 0.5 mm it is clean. In
    between it was handed back with no mention at all, so a flawed opening looked
    exactly like a good one in every list the tool shows -- and on the two cases this
    was traced from, a CLEAN cut sat 0.6 and 0.8 mm further back and suggest_pullback
    finds it in both. The refusal path has told people that for a while; the accepted
    path never did.
    """
    out = []
    for t in (terminals if terminals is not None else tree.terminals):
        dev = getattr(t, "_rim_dev_mm", None)
        if dev is not None and float(dev) > RIM_FLAW_MM:
            out.append((t, float(dev)))
    return sorted(out, key=lambda p: -p[1])


def _ray_radius(mesh, segment, ss, tangents, k=12):
    """Lumen radius measured ACROSS the vessel, by rays that cannot see a cap.

    The last resort for a segment so short that its whole middle is beside the
    model's end cap, where there is no wall sample left to take a median of. Rays
    are cast perpendicular to the tangent, k of them evenly around the axis at each
    station, and the median first-hit distance is the radius. A cap is broadside to
    these rays rather than in front of them, so it cannot be mistaken for the wall.

    Costs 20-80 ms for a segment and runs for about one segment per case, so it is
    the fallback and never the main path.
    """
    origins, dirs = [], []
    for s, t in zip(ss, tangents):
        p = segment.at_arclength(float(s))[0]
        a = np.array([1.0, 0.0, 0.0])
        if abs(float(t @ a)) > 0.9:
            a = np.array([0.0, 1.0, 0.0])
        u = np.cross(t, a)
        u = u / max(float(np.linalg.norm(u)), 1e-9)
        v = np.cross(t, u)
        for th in np.linspace(0.0, 2.0 * np.pi, k, endpoint=False):
            origins.append(p)
            dirs.append(np.cos(th) * u + np.sin(th) * v)
    if not origins:
        return float("nan")
    origins = np.asarray(origins, float)
    try:
        loc, idx, _tri = mesh.ray.intersects_location(
            origins, np.asarray(dirs, float), multiple_hits=False)
    except Exception:                                   # noqa: BLE001
        return float("nan")
    if not len(loc):
        return float("nan")
    return float(np.median(np.linalg.norm(loc - origins[idx], axis=1)))


def segment_radius(mesh, segment, lo=0.2, hi=0.8, n=25):
    """Median lumen radius over the MIDDLE of a segment.

    Sampled away from the ends on purpose: signed_distance near a terminal
    measures the distance to the nearby *cap*, not to the vessel wall, which
    badly understates the radius (the aortic root reads ~7 mm at the tip but
    ~21 mm mid-segment). The middle is what actually separates trunk from branch.

    SAMPLING THE MIDDLE IS NOT ENOUGH, and for years this quietly returned the
    wrong number on exactly the segments that matter most. On a SHORT segment at a
    model end the middle is still beside the cap, so the dodge above fails and what
    comes back is the distance to the end of the model. The cleanest demonstration
    is the same centerline and the same 25 stations measured twice, differing only
    in whether the cap is there: JK1_2's Line06 reads 8.586 mm on the input mesh
    and 14.673 mm once the inlet has been cut away. 8.586 is what everything
    downstream used. Across the abdominal cohort it is 39 of 891 segments, always
    about one per case, always the short segment carrying the inlet.

    So each sample is asked WHICH SURFACE it found: split the vector to the closest
    point into the part along the centerline and the part across it, and a sample
    whose nearest surface is more along than across measured a cap, not a wall.
    Take the median over the wall samples. When too few survive -- the short
    inlet segments, where every station is beside the cap -- fall back to rays cast
    across the vessel, which a cap cannot answer.

    The fallback only ever raises a radius that was reading a cap. A NEGATIVE
    median is left exactly alone, because that is the signal that the centerline
    lies outside the mesh, and mismatched() and _no_radius_reason() are built on it;
    rays fired from outside the surface would come back with a confident positive
    number and hide a mesh that does not match its centerline at all.

    signed_distance needs a CLOSED surface to tell inside from outside. On an open
    one the sign is meaningless and the radius comes back NEGATIVE, which sorts the
    vessel to the bottom of every ranking that follows. Two ordinary things produce
    an open mesh: a segmentation that arrived with holes (3 of 44 abdominal cases
    here), and -- unavoidably -- a mesh this tool cut and the user loaded back in to
    carry on. Reloading JK1_1's own export turns its inlet from +7.98 mm into
    -14.61 mm, drops the aorta out of the trunk ranking, and takes every anatomical
    name with it.

    The DISTANCE is still right when the sign is not, so measure that instead. A
    closed mesh takes the original path and reads exactly as it did before.
    """
    L = segment.length
    ss = np.linspace(lo * L, hi * L, n)
    at = [segment.at_arclength(s) for s in ss]
    pts = np.array([a[0] for a in at])
    tans = np.array([a[1] for a in at])
    # WHICH surface each station found. closest_point is asked for the point as well
    # as the distance, because the direction to it is the whole test: a wall lies
    # across the vessel, a cap lies along it.
    cp, d_un, tid = trimesh.proximity.closest_point(mesh, pts)
    axial = np.abs(np.einsum("ij,ij->i", cp - pts, tans))
    across = np.sqrt(np.maximum(d_un ** 2 - axial ** 2, 0.0))
    on_wall = axial <= across
    if mesh.is_watertight:
        # INSIDE OR OUT, from the face we already found rather than from a ray cast.
        #
        # signed_distance re-derives this same distance and adds a sign by firing rays
        # through the whole mesh to count crossings. Without embree that is pure
        # Python and it dominates everything: 0.24 s per segment, 5.5 s of the 9.7 s
        # a 13-opening case takes to open, more than every cut in it put together. The
        # sign was the only thing being bought.
        #
        # closest_point already returns the triangle it landed on, and a point is
        # inside when it sits behind that triangle's outward normal. Decided ONCE per
        # segment by majority rather than per station: the two tests disagree on 27 of
        # 22,750 stations across the cohort (0.119%), always where the nearest point
        # falls on an edge between faces whose normals diverge, and a majority vote
        # cannot be moved by that. Checked against signed_distance over 910 segments:
        # not one median changes sign, and cut_bench is byte-identical over 503 cuts.
        #
        # Getting this wrong is not a rounding error. A valid centerline that reads as
        # OUTSIDE comes back with a negative radius, drops out of the trunk ranking,
        # and takes every anatomical name with it -- see the note below on reloading
        # JK1_1's own export.
        behind = np.einsum("ij,ij->i", pts - cp, np.asarray(mesh.face_normals)[tid])
        d = d_un if float((behind < 0.0).mean()) >= 0.5 else -d_un
    else:
        d = d_un                       # already measured above, do not query twice
    if on_wall.sum() >= max(5, int(0.4 * n)):
        return float(np.median(d[on_wall]))
    med = float(np.median(d))
    if med <= 0.0:
        return med                     # outside the mesh: that signal must survive
    r = _ray_radius(mesh, segment, ss, tans)
    return r if r == r else med        # NaN means the rays found nothing either


def detect_cut_mode(mesh, terminal, bias_deg=5.0):
    """Pick Z-plane vs centerline-normal by LOOKING at the cap the mesh already has.

    The segmentation was cut somehow, and the flat face it left behind records how.
    On KY160 the aorta's caps sit 7-10 deg off horizontal but 17-25 deg off
    centerline-perpendicular, so a Z-plane restores the original cut and a normal
    one does not. That is a fact about that mesh, not a law: another case may have
    been cut perpendicular to the vessel, and forcing a Z-plane on it puts a
    horizontal plane through an arch, where the vessel runs horizontal and the cut
    degenerates.

    So measure. Faces near the tip whose normals agree with the outward direction
    ARE the cap; compare its orientation against both candidates and take the
    closer, with a small bias toward the (safer) centerline-normal cut.

    Rounded tips have no flat cap to read, and fall back to centerline-normal.
    """
    tip = np.asarray(terminal.tip, float)
    _, outward = terminal.cut_plane(0.0)

    r = max(float(getattr(terminal.segment, "_radius_mm", None) or 5.0), 1.0)
    centers = np.asarray(mesh.triangles_center)
    near = np.linalg.norm(centers - tip, axis=1) < max(2.5 * r, 8.0)
    if near.sum() < 5:
        return "normal", None

    normals = np.asarray(mesh.face_normals)
    cap = near.copy()
    cap[np.flatnonzero(near)] = normals[near] @ outward > 0.75
    if cap.sum() < 5:
        return "normal", None                 # rounded tip: nothing flat to read

    areas = np.asarray(mesh.area_faces)
    n = np.average(normals[cap], axis=0, weights=areas[cap])
    n /= max(np.linalg.norm(n), 1e-12)

    tilt_z = np.degrees(np.arccos(min(abs(float(n[2])), 1.0)))
    tilt_cl = np.degrees(np.arccos(min(abs(float(n @ outward)), 1.0)))
    mode = "z" if tilt_z + bias_deg < tilt_cl else "normal"
    return mode, (tilt_z, tilt_cl)


# A second terminal is a trunk END only if it is comparable in calibre to the
# fattest one. Kept as the FALLBACK for a vessel with no junction to read -- an
# unbranched tube, or a tree whose structure could not be walked. It is a poor
# primary test and is no longer used as one: see PAIR_AREA_FRACTION.
TRUNK_RATIO = 0.6

# An aorta either runs END TO END through the model -- a thoracic one, entering at
# the root and leaving by the distal descending end -- or it enters once and
# TERMINATES BY DIVIDING, which is what an abdominal one does at the iliac
# bifurcation. Everything downstream turns on which: the number of trunk ends
# decides whether the openings are named ascending/descending or inlet/iliac, and
# region_presets offers Ishimaru zones to one and the aortic neck to the other.
#
# Deciding it on the second-fattest terminal's CALIBRE reads the wrong answer on a
# fifth of real cases, in both directions. A common iliac is legitimately
# comparable to the aorta feeding it, so an abdominal case trips the "two ends"
# test (JK1_1: 0.657 against a 0.60 cut-off); and an aneurysm makes a thoracic
# aorta's two ends wildly unequal, so it fails the same test the other way
# (KK40_3: 0.435). Over 44 abdominal and 30 thoracic cases the calibre rule got
# 58 of 74 right -- 14 of the 44 abdominals were treated as thoracic.
#
# The structural question is better posed and easier to measure. Walking down the
# trunk, is there a junction where the outflow DIVIDES, rather than a branch merely
# leaving? Flow divides with cross-sectional AREA, not with radius, so the two
# sides are compared by summed area over the terminals each leads to. A true
# bifurcation sends a quarter or more of the area down its lesser side; a branch
# takeoff never does.
#
# Measured on those 74 cases: thoracic reaches 0.214 at worst (with a single
# outlier at 0.309, KP35_1, whose first-order vessel feeds three great vessels),
# abdominal 0.269 at least. A quarter sits inside that gap and states the anatomy
# plainly. Together with the open-surface fix in segment_radius -- which is what
# rescues JK6_1 -- it gets 73 of 74. The one residual, KP35_1, is reported as a
# near call rather than left silent. See tools/ANATOMY.md.
#
# THAT GAP DOES NOT SURVIVE THE SMOOTHING THE TOOL ITSELF APPLIES. Re-measured over
# 79 cases on the centerlines actually cut -- smoothed when jagged, which is 52 of
# them -- the two families OVERLAP: thoracic reaches 0.249 and abdominal falls to
# 0.225. No cut-off can separate those, so a threshold is not the thing to tune.
# 0.25 still sits in the best place available (between the 0.249 thoracic and the
# 0.270 abdominal) and scores 77 of 79; the two it misses, KP35_1 at 0.309 and
# JK25_1 at 0.225, are on the wrong side of EVERY possible cut-off relative to each
# other. This is the rule's ceiling rather than a tuning failure, which is why both
# front ends report a near call and let the family be set by hand.
PAIR_AREA_FRACTION = 0.25


def _through(jg, junction, seg_name):
    """Segments reachable from *junction* by leaving through *seg_name*."""
    a, b = jg[seg_name]
    far = b if a == junction else a
    if far is None:
        return {seg_name}
    return {seg_name} | _subtree_segments(jg, far, seg_name)


def trunk_walk(tree, rad, jg=None, start=None):
    """Follow the trunk from *start*, always continuing into the greater outflow.

    Returns (chain, splits). *chain* is the ordered segment names from *start*'s own
    segment to wherever the trunk runs out; *splits* has one entry per junction
    passed, (fraction, junction_index, (near_seg, far_seg)), where the fraction is
    the lesser side's share of the greater side's cross-sectional area. The junction
    in splits[i] sits between chain[i] and chain[i + 1].

    Always taking the greater area is what makes this the TRUNK rather than an
    arbitrary path through the tree: at every junction the aorta carries more area
    than the branch leaving it, right down to the iliac bifurcation, where the two
    sides are comparable and the walk simply continues into the larger iliac.

    Starting from every terminal in turn was tried and is much worse: begin at a
    renal and the first junction offers "up the aorta" against "down the aorta", two
    large and comparable sides, which reads as a bifurcation on a thoracic case --
    59 of 74 against 72 in the same sweep.
    """
    jg = jg or _junction_graph(tree)
    term_of_seg = {}
    for t in tree.terminals:
        term_of_seg.setdefault(t.segment.name, []).append(t)

    def area(names):
        # Negative radii mean an open surface (see segment_radius); they carry no
        # area rather than subtracting it.
        return float(sum(max(float(rad.get(t.label, 0.0)), 0.0) ** 2
                         for n in names for t in term_of_seg.get(n, [])))

    if start is None:
        start = max(tree.terminals, key=lambda t: float(rad.get(t.label, 0.0)))
    seg = start.segment.name
    ns, ne = jg[seg]
    j = ne if ns is None else ns
    came, chain, splits = seg, [seg], []
    for _ in range(len(tree.segments) + 1):          # cannot revisit: a tree
        if j is None:
            break
        opts = [n for n, _e in tree.junctions[j][1] if n != came]
        if not opts:
            break
        scored = sorted(((area(_through(jg, j, n)), n) for n in opts),
                        key=lambda s: -s[0])
        top = scored[0]
        rival = scored[1] if len(scored) > 1 else (0.0, None)
        frac = (rival[0] / top[0]) if top[0] > 0 else 0.0
        splits.append((frac, j, (top[1], rival[1])))
        came = top[1]
        chain.append(top[1])
        a, b = jg[top[1]]
        j = b if a == j else a
    return chain, splits


def outflow_split(tree, rad, jg=None, start=None):
    """The most even division of outflow found walking down the trunk.

    Returns (fraction, junction_index, (near_seg, far_seg)). The fraction is the
    lesser side's share of the greater side's cross-sectional area, so 1.0 is a
    perfectly even bifurcation and ~0.1 an ordinary branch leaving the aorta.

    A junction that divides nothing -- fraction exactly zero, which is what an
    unmeasured surface gives everywhere -- is not an answer, so it reports no
    junction at all rather than the first one it walked past.
    """
    _chain, splits = trunk_walk(tree, rad, jg, start)
    best = (0.0, None, (None, None))
    for s in splits:
        if s[0] > best[0]:
            best = s
    return best


def classify_terminals(tree, mesh, z_cut_ends=True, n_trunk_ends=None, verbose=True):
    """Label the trunk ends vs the branch tips and set each terminal's cut mode.

    Nothing here is specific to one dataset. The trunk ends are the fattest
    terminals by a wide margin, so ranking on mid-segment radius and taking the
    top ones is unambiguous. A THORACIC aorta has two (the root and the distal
    descending end, the proximal being the superior of the pair); an ABDOMINAL
    aorta has one (the inlet) and leaves through its iliac and visceral branches.
    Which it is, is decided by calibre rather than assumed -- see TRUNK_RATIO --
    and *n_trunk_ends* pins it by hand if a case ever needs it.

    An UNBRANCHED aorta simply has two terminals and no branches, and falls out of
    the same logic -- the two cut planes then just isolate a length of trunk.

    The separation ratio is reported so a marginal call is visible rather than
    silent; the GUI lets any terminal be re-labelled or re-moded by hand.
    """
    rad = {}
    for t in tree.terminals:
        t.segment._radius_mm = segment_radius(mesh, t.segment)   # cached for cut_terminal
        rad[t.label] = t.segment._radius_mm
    order = sorted(tree.terminals, key=lambda t: -rad[t.label])

    if len(order) < 2:
        raise ValueError(
            f"need at least 2 terminals, found {len(order)}, "
            f"the centerline must run end-to-end through the vessel"
        )

    r0, r1 = rad[order[0].label], rad[order[1].label]
    if r0 <= 0.0:
        raise ValueError(_no_radius_reason(mesh, tree))
    ratio = r1 / r0

    # Which vessels' centerlines lie OUTSIDE the surface. preflight already refuses
    # to CUT on this and says why, but it runs long after the anatomy has been
    # decided and put on screen: the family call, every vessel name and every region
    # preset come off these radii first.
    #
    # Measured on the thoracic cohort, pairing a branched centerline with the
    # trunk-only mesh sitting beside it: 4 to 5 of 6 terminals came out negative, one
    # at -68 mm, and only 58-75% of the centerline lay inside the mesh's bounding box
    # at all. outflow_split then finds no dividing junction -- every subtree measures
    # zero area -- so the family quietly falls back to the 2nd/1st radius rule, which
    # is right on 58 of 74. KK40_3, KR47_1 and KR52_1 were all called abdominal that
    # way, and all three read correctly the moment the matching mesh is used.
    #
    # So it is recorded here, and the call it produced is marked untrustworthy rather
    # than presented as a reading of the anatomy.
    tree.radius_outside = [t.label for t in tree.terminals if rad[t.label] <= 0.0]
    tree.coverage = centerline_coverage(tree, mesh)
    tree.mesh_matches = (tree.coverage >= MIN_CENTERLINE_COVERAGE
                         and not tree.radius_outside)

    # Does the trunk end by DIVIDING (an iliac bifurcation) or by running out at a
    # tip (a descending aorta)? That, not calibre, is what one-vs-two trunk ends
    # actually asks -- see PAIR_AREA_FRACTION.
    jg = _junction_graph(tree)
    split, split_j, split_segs = 0.0, None, (None, None)
    if len(order) > 2 and tree.junctions:
        split, split_j, split_segs = outflow_split(tree, rad, jg)
    tree.outflow_split = split if split_j is not None else None
    # A branched tree whose junctions divide NOTHING has no measurable outflow, which
    # on real anatomy only happens when the radii are wrong. Together with a vessel
    # outside the surface that is the signature of a mismatched mesh, and it is
    # exactly when the fallback rule is least worth trusting.
    tree.family_trusted = bool(tree.mesh_matches) and not (
        len(order) > 2 and tree.junctions and split_j is None)

    if n_trunk_ends is None:
        # No branches at all means an unbranched vessel: both ends ARE the trunk.
        if len(order) == 2 or split_j is None:
            n_trunk_ends = 2 if (len(order) == 2 or ratio >= TRUNK_RATIO) else 1
        else:
            n_trunk_ends = 1 if split >= PAIR_AREA_FRACTION else 2

    n_trunk = min(max(int(n_trunk_ends), 1), len(order))
    trunk, branches = order[:n_trunk], order[n_trunk:]

    # Inlet first -- read off the anatomy, not off the z axis (see _proximal_end).
    prox_end, how = _proximal_end(tree, trunk, jg, rad)
    trunk = [prox_end] + [t for t in trunk if t is not prox_end]
    tree.proximal_from = how
    if len(trunk) >= 2:
        trunk[0].name, trunk[1].name = "Aorta root (proximal)", "Aorta distal"
    else:
        trunk[0].name = "Aorta inlet"
    for t in tree.terminals:
        t.is_trunk = t in trunk               # the rest are branches
        t.is_inlet = (t is trunk[0])

    # Read each trunk cap off the mesh rather than assuming. Forcing a Z-plane on
    # a mesh that was cut perpendicular to its vessel drops a horizontal plane
    # into an arch, where the vessel runs horizontal and the cut degenerates.
    caps = {}
    for t in trunk:
        if z_cut_ends:
            t.mode, caps[t.label] = detect_cut_mode(mesh, t)
            # The aorta's OWN ends are wanted horizontal. That is how the great-vessel
            # exports were cut and what the paper's figures show, and detect_cut_mode
            # reads it off the cap correctly on most cases -- but a rounded or
            # remeshed tip has no flat cap to read and it falls back to
            # centerline-normal, so the same anatomy came out two different ways
            # depending on how the segmentation happened to finish. Asked for from the
            # tool: the ascending and descending ends should start on Z.
            #
            # Except where a horizontal plane would not cut ACROSS the vessel. Past
            # MAX_OBLIQUITY_DEG it starts running along it and the opening tears --
            # which is the same fallback nub_cut applies to these ends, measured
            # there on three cases in thirty at 52, 56 and 58 degrees.
            # Measured WITH the horizontal plane in place: obliquity_deg reads the
            # angle of whatever mode the terminal currently has, and a
            # centerline-normal plane is perpendicular by construction, so asking
            # before switching always answers 0 and always switches. That cost two
            # thoracic cases before it was caught.
            was = t.mode
            t.mode = "z"
            if float(obliquity_deg(t, 6.0)) > MAX_OBLIQUITY_DEG - OBLIQUITY_MARGIN_DEG:
                t.mode = was
        else:
            t.mode, caps[t.label] = "normal", None

    # Margin of the call: the leanest trunk end against the FATTEST branch, i.e.
    # the closest competitor. (Ranking the branches by height first and comparing
    # against the topmost one flatters the result, it read 5.4x when the real
    # margin was 3.0x.)
    sep = (rad[trunk[-1].label] / max(max(rad[b.label] for b in branches), 1e-6)
           if branches else float("inf"))
    tree.confidence = sep

    branches = sorted(branches, key=lambda t: -t.tip[2])
    for i, t in enumerate(branches, start=1):
        t.name = f"Branch {i}"
        t.mode = "normal"

    if verbose:
        kind = "unbranched" if not branches else f"{len(branches)} branch(es)"
        ends = "1 trunk end (inlet)" if len(trunk) == 1 else f"{len(trunk)} trunk ends"
        print(f"terminal classification ({kind}, {ends}, by mid-segment radius):")
        for t in order:
            extra = ""
            c = caps.get(t.label)
            if c:
                extra = (f"   cap is {c[0]:.0f}° off horizontal, "
                         f"{c[1]:.0f}° off perpendicular")
            print(f"  {t.name:22s} {t.label:12s} r={rad[t.label]:6.2f} mm  "
                  f"tip_z={t.tip[2]:8.1f}  cut={t.mode}{extra}")
        if branches:
            print(f"  trunk/branch radius separation: {sep:.1f}x"
                  + ("   <-- MARGINAL, check the labels" if sep < 1.5 else ""))
        else:
            print("  no branches: the two planes isolate a length of trunk")
        # Show the evidence the 1-vs-2 call was made on, and say when it was close.
        # An aorta read as the wrong family is not a subtle error -- it renames
        # every opening and changes which regions can be isolated -- so the number
        # behind it belongs on screen rather than in the code.
        if split_j is not None:
            near, far = split_segs
            print(f"  outflow divides most evenly at junction {split_j} "
                  f"({near} vs {far}): the lesser side carries {split:.2f} of the "
                  f"greater side's area, cut-off {PAIR_AREA_FRACTION:.2f} "
                  f"-> {n_trunk} trunk end(s)")
            if abs(split - PAIR_AREA_FRACTION) < 0.08:
                other = 1 if n_trunk == 2 else 2
                print(f"  NOTE: that is close to the cut-off, so {other} trunk "
                      f"end(s) is arguable. Set the family by hand if this is wrong.")
        elif len(order) > 2:
            print(f"  no junction could be walked; fell back to 2nd/1st terminal "
                  f"radius = {ratio:.2f} against {TRUNK_RATIO:.2f} "
                  f"-> {n_trunk} trunk end(s)")
        if tree.coverage < MIN_CENTERLINE_COVERAGE:
            print(f"  WARNING: only {100 * tree.coverage:.0f}% of the centerline lies "
                  f"within this mesh's bounds. Every correctly paired case in both "
                  f"lab cohorts measures 100%. This mesh and this centerline are "
                  f"most likely different exports of different extent -- the anatomy "
                  f"below is measured off them and should not be believed.")
        if tree.radius_outside:
            print(f"  WARNING: {len(tree.radius_outside)} of {len(tree.terminals)} "
                  f"centerlines run OUTSIDE this mesh "
                  f"({', '.join(tree.radius_outside[:4])}"
                  f"{' ...' if len(tree.radius_outside) > 4 else ''}). The usual "
                  f"cause is a trunk-only mesh paired with a branched centerline. "
                  f"Every name and region below is measured off those radii, so "
                  f"treat the anatomy as unread until the pair is fixed.")
        # A radius can only come out non-positive on an open surface, where the
        # inside/outside sign is undefined (see segment_radius). The measurement
        # falls back to unsigned distance, but a vessel whose centerline genuinely
        # leaves the mesh is a data fault worth naming, not absorbing.
        if not mesh.is_watertight:
            # Deliberately no opening COUNT here: a surface loaded without merging
            # its vertices reports one boundary per edge (34,783 on a 34,783-face
            # export), and a number that wrong is worse than no number.
            print("  NOTE: this surface is NOT closed, so radii are measured as "
                  "unsigned distance rather than signed. If this is a mesh the tool "
                  "already cut, its centerlines now run past the new rims, so the "
                  "radii read long there and the vessel names may not be reliable.")

    tree.terminals = list(trunk) + branches
    return tree


def _junction_graph(tree):
    """seg name -> (start_junction | None, end_junction | None)."""
    jof = {}
    for i, (pos, inc) in enumerate(tree.junctions):
        for (n, e) in inc:
            jof[(n, e)] = i
    return {name: (jof.get((name, False)), jof.get((name, True)))
            for name in tree.segments}


def _subtree_segments(seg_nodes, start_junction, via_seg):
    """Segment names reachable from *start_junction* without crossing *via_seg*."""
    adj = {}
    for name, (ns, ne) in seg_nodes.items():
        if ns is not None:
            adj.setdefault(ns, []).append((name, ne))
        if ne is not None:
            adj.setdefault(ne, []).append((name, ns))
    seen_seg, seen_node, stack = set(), {start_junction}, [start_junction]
    while stack:
        j = stack.pop()
        for (name, other) in adj.get(j, []):
            if name == via_seg:
                continue
            seen_seg.add(name)
            if other is not None and other not in seen_node:
                seen_node.add(other)
                stack.append(other)
    return seen_seg


def tree_root_terminal(tree):
    """The vessel's inlet: the most superior trunk end, else the fattest terminal.

    Rooting the tree here is what gives every interior segment a single, unambiguous
    DOWNSTREAM direction -- the way blood leaves. Without it, "trim this vessel
    before it splits" has two readings, and the wrong one proposes trimming an iliac
    by deleting the entire aorta above it.
    """
    if not tree.terminals:
        return None
    flagged = [t for t in tree.terminals if getattr(t, "is_inlet", False)]
    if flagged:
        return flagged[0]            # decided by classify_terminals, z-free
    trunks = [t for t in tree.terminals if getattr(t, "is_trunk", False)]
    if trunks:
        return max(trunks, key=lambda t: float(t.tip[2]))
    # Fallback for a tree that has not been classified yet. Radii are cached by
    # classify_terminals, so they are unset in exactly this case -- fall back to the
    # most superior terminal, which needs no mesh, rather than to an arbitrary one.
    rad = [float(getattr(t.segment, "_radius_mm", None) or 0.0) for t in tree.terminals]
    if max(rad) > 0.0:
        return tree.terminals[int(np.argmax(rad))]
    return max(tree.terminals, key=lambda t: float(t.tip[2]))


def pre_bifurcation_candidates(tree, mesh=None):
    """Interior segments that can be cut BEFORE they split -- "trim a branch group".

    A candidate is an interior segment (a junction at BOTH ends) whose DOWNSTREAM
    side -- the side facing away from the inlet -- carries two or more openings.
    Cutting there removes that whole group and leaves ONE opening on the feeding
    vessel: a common iliac trimmed before it splits into internal and external, or
    a celiac trunk before its branches. One opening downstream is not a candidate,
    since that is just an ordinary cut on the branch itself.

    Sorted smallest-first, so the tight anatomical trims come before the drastic
    ones (an aorta cut just above its bifurcation legitimately removes every iliac).

    Returns dicts: {segment, at_end (the split-facing end), subtree (segment names
    removed), terminals (the openings collapsed into one), junction (its position)}.
    """
    root = tree_root_terminal(tree)
    if root is None:
        return []
    seg_nodes = _junction_graph(tree)
    out = []
    for name in tree.segments:
        ns, ne = seg_nodes[name]
        if ns is None or ne is None:
            continue                          # a free end -> not an interior segment
        # Downstream is the side whose subtree does NOT lead back to the inlet. Cutting
        # this segment must SEPARATE the two sides: if the inlet is reachable from both
        # (a cycle in the centerline graph) or from neither (a disconnected component),
        # there is no downstream and the "trim" would delete the model.
        sub_end = _subtree_segments(seg_nodes, ne, name)
        sub_start = _subtree_segments(seg_nodes, ns, name)
        in_end, in_start = (root.segment.name in sub_end), (root.segment.name in sub_start)
        if in_end == in_start:
            continue
        at_end = not in_end
        sub = sub_end if at_end else sub_start
        terms = [t for t in tree.terminals if t.segment.name in sub]
        if len(terms) < 2 or root in terms:
            continue                          # one opening past it = an ordinary cut
        # A trim needs somewhere to put the plane, and that used to be tested as
        # "at least four radii long". The radius it tested was wrong in exactly the
        # case the rule mattered. An interior segment runs from ITS PARENT'S AXIS to
        # the split, so its first stretch lies inside the parent -- and a segment
        # leaving a 20 mm arch has 20 mm of centerline inside the aorta before it
        # reaches its own ostium. Every face near that stretch is aortic wall, so the
        # measured radius is the aorta's: KY160's brachiocephalic trunk came out at
        # 12.6 mm, needed 50.6 mm of length, had 29.3, and was never offered. Cut 2 mm
        # before its split it comes out with a rim 0.00 mm off plane over 91 points.
        #
        # So the check is only that there is a vessel here at all. Whether a given
        # offset works is a question about that offset, which the cut answers when it
        # is asked, and the pull-back button answers before it is.
        seg = tree.segments[name]
        if getattr(seg, "_radius_mm", None) is None and mesh is not None:
            seg._radius_mm = segment_radius(mesh, seg)   # interior segments are uncached
        if seg.length < 2.0 * MIN_OFFSET_MM:
            continue
        jidx = ne if at_end else ns
        out.append({"segment": name, "at_end": at_end, "subtree": sub,
                    "terminals": terms,
                    "junction": np.asarray(tree.junctions[jidx][0], float)})
    # Fewest openings swallowed first; anything that would take a trunk end with it
    # (cutting the aorta itself, say) is legitimate but drastic, so it sorts last
    # rather than heading the list. Ties go to the roomiest vessel -- the one most
    # likely to cut clean -- then by name, so the menu order is reproducible rather
    # than following CSV row order.
    out.sort(key=lambda c: (any(getattr(t, "is_trunk", False) for t in c["terminals"]),
                            len(c["terminals"]),
                            -tree.segments[c["segment"]].length,
                            c["segment"]))
    return out


# ---------------------------------------------------------------------------
# Anatomical naming
#
# The frame these names assume. KY160's .mhd declares AnatomicalOrientation = RAI,
# and both sample datasets agree with it independently: the aortic root sits ABOVE
# the distal descending aorta (+z superior); the vessels off the brachiocephalic
# trunk -- right-sided by definition -- point to -x (+x left); and the celiac and
# SMA sit at far lower y than the aorta they leave (-y anterior).
#
# A mirrored export would flip left from right without changing anything else, and
# nothing in the geometry could tell. So the assumption is RECORDED and reported
# rather than buried, and only sides that follow from it are ever stated.
ANATOMY_FRAME = "+x = patient left, +y = posterior, +z = superior (RAI)"

# Every "I could not confirm handedness" note ends with this. Both families emit one
# -- the arch when no brachiocephalic trunk is found, the abdomen when no
# splenic/hepatic pair is -- and when one is present the Left/Right names rest on the
# frame ALONE rather than on a witness in the anatomy. A front end has to be able to
# find either without matching the whole sentence, because that difference is the
# difference between a reading and a guess, and the user is entitled to know which
# one they are looking at.
SIDES_ASSUMED_NOTE = f"sides assume {ANATOMY_FRAME}"

# How far a region boundary must stay from a branch ostium, in that branch's radii.
# Measured on JK5: a renal (r~2.5-2.9 mm) on a 25 mm aorta needs ~8 mm of clearance
# before its mouth stops being clipped, i.e. about three of its own radii.
OSTIUM_CLEAR_R = 3.0

# Side/lateral thresholds are expressed in AORTIC RADII, not millimetres, and are
# measured from the aorta's own axis rather than from x = 0. Absolute world
# coordinates are meaningless here: the same arch exported with a different origin
# scrambled every label (a carotid read as a subclavian) until this was relative.
SIDE_TOL_R = 0.35        # off-axis by this many aortic radii before a side is called
CAROTID_MAX_R = 1.6      # ascends near the midline
SUBCLAVIAN_MIN_R = 2.5   # reaches far laterally


# The artery names that can be assigned BY HAND, per family.
#
# They live here rather than in the interface because region_presets matches on the
# words in them -- "renal", "celiac", "sma", "mesenteric" -- so a name typed in
# different words looks assigned and silently builds no region. Offering a list is
# what makes a hand-named renal actually restore the aortic neck.
#
# This is not a fallback for a weak rule. Over 47 abdominal scans name_anatomy names
# no renal on 18, and on 14 of those because SEVERAL lateral pairs fit the renal
# pattern equally well and naming either would be a guess presented as a fact. No
# geometric rule should settle that; a person looking at the model can. Accessory and
# polar renals alone are ~30% of people. See tools/ABDOMINAL.md.
ANATOMY_NAMES_ABDOMINAL = (
    "Celiac trunk",
    "Common hepatic artery",
    "Splenic artery",
    "Left gastric artery",
    "Superior mesenteric artery (SMA)",
    "Inferior mesenteric artery (IMA)",
    "Left renal artery",
    "Right renal artery",
    # Accessory and polar renals carry "renal" deliberately: the infrarenal neck
    # begins below the LOWEST renal artery whichever one that is, so they have to
    # count toward that boundary rather than be excluded from it.
    "Left accessory renal artery",
    "Right accessory renal artery",
    "Left common iliac",
    "Right common iliac",
    "Left internal iliac",
    "Right internal iliac",
    "Left external iliac",
    "Right external iliac",
    "Inferior phrenic artery",
    "Lumbar artery",
    "Median sacral artery",
)

ANATOMY_NAMES_THORACIC = (
    "Ascending aorta",
    "Descending thoracic aorta",
    "Brachiocephalic (innominate) trunk",
    "Right common carotid artery",
    "Left common carotid artery",
    "Right subclavian artery",
    "Left subclavian artery",
    "Right vertebral artery",
    "Left vertebral artery",
    "Intercostal artery",
    "Bronchial artery",
)


def anatomy_choices(tree):
    """The artery names offered for hand-assignment on this anatomy."""
    n_trunk = sum(1 for t in tree.terminals if getattr(t, "is_trunk", False))
    return (ANATOMY_NAMES_ABDOMINAL if n_trunk == 1 else ANATOMY_NAMES_THORACIC)


def set_manual_name(tree, label, name):
    """Assign (or clear, with a falsy *name*) one opening's artery by hand."""
    manual = dict(getattr(tree, "manual_names", None) or {})
    if name:
        manual[str(label)] = str(name)
    else:
        manual.pop(str(label), None)
    tree.manual_names = manual
    return manual


def _apply_manual(tree, names, vessels, notes):
    """Hand-assigned names win, and displace any automatic name they duplicate.

    The displacement is not tidiness. region_presets places the aortic neck at the
    LOWEST landmark matching "renal", so an automatically named renal left sitting
    beside a hand-named one would put that boundary on whichever is lower -- quite
    possibly the vessel the person was correcting.
    """
    manual = {str(k): str(v) for k, v in
              (getattr(tree, "manual_names", None) or {}).items() if v}
    if not manual:
        return names, vessels, notes
    known = {t.label for t in tree.terminals}
    used = 0
    for label, name in manual.items():
        if label not in known:
            continue                     # left over from a different model
        for lbl, nm in list(names.items()):
            if nm == name and lbl != label:
                names.pop(lbl, None)
        for seg, nm in list(vessels.items()):
            if nm == name:
                vessels.pop(seg, None)
        names[label] = name
        used += 1
    if used:
        notes.append(f"{used} opening(s) named by hand -- those names override "
                     f"anything read from the geometry")
    return names, vessels, notes


def _aorta_chain(tree, jg, seg_a, seg_b):
    """Ordered segment names along the path from *seg_a* to *seg_b* -- the trunk."""
    if seg_a == seg_b:
        return [seg_a]
    adj = {}
    for name, (ns, ne) in jg.items():
        for j in (ns, ne):
            if j is not None:
                adj.setdefault(j, set()).add(name)
    prev, queue = {seg_a: None}, [seg_a]
    while queue:
        cur = queue.pop(0)
        if cur == seg_b:
            chain = []
            while cur is not None:
                chain.append(cur)
                cur = prev[cur]
            return chain[::-1]
        for j in jg[cur]:
            if j is None:
                continue
            for nxt in adj.get(j, ()):
                if nxt not in prev:
                    prev[nxt] = cur
                    queue.append(nxt)
    return [seg_a, seg_b]


def _aorta_path(tree, jg, seg_a, seg_b):
    """Segment names on the path between two segments -- i.e. the trunk itself."""
    return set(_aorta_chain(tree, jg, seg_a, seg_b))


def _aorta_stations(tree, jg, seg_a, seg_b):
    """Junction index -> how far along the aorta it sits, from *seg_a*'s free end.

    This is what "the first branch off the arch" actually means. Ranking takeoffs by
    HEIGHT instead gets a type III arch backwards, because there the brachiocephalic
    origin lies below the vessels distal to it.
    """
    chain = _aorta_chain(tree, jg, seg_a, seg_b)
    stations, d = {}, 0.0
    for i, name in enumerate(chain[:-1]):
        d += tree.segments[name].length
        shared = ({j for j in jg[name] if j is not None}
                  & {j for j in jg[chain[i + 1]] if j is not None})
        for j in shared:
            stations.setdefault(j, d)
    return stations


def _branch_takeoffs(tree, jg, trunk, chain, stations):
    """Where each branch leaves the aorta, in mm along it from the first trunk end.

    A branch counts at the junction where its WHOLE VESSEL TREE joins the aorta, not
    only when its own segment touches it. Most arches put the right subclavian and the
    right carotid on the brachiocephalic, so neither touches the aorta itself; and a
    centerline that draws the left carotid and left subclavian off a short common stem
    leaves no branch touching it at all. Counting only direct attachments then found no
    takeoffs, _proximal_end fell through to the taper rule, and on ZF_Example's KK16_6
    (_2 centerline, 2 Oct 2026) -- whose descending end is the fatter, r 21.6 against
    14.6 mm -- the descending aorta was named the root and zone 0 ran 316 mm down it.

    ONE VOTE PER TAKEOFF, NOT PER TIP. What says "this is the ascending end" is where the
    branches leave the aorta, not how many tips each group ends in. Counting tips let a
    bushy group pull the mean: a celiac trunk with three tips low on a thoracoabdominal
    descending aorta outvoted the arch and flipped a root the direct-only rule had right.
    So the result is the set of distinct takeoff stations.

    A branch whose vessels reach the aorta at MORE THAN ONE junction -- a loop, such as a
    dissection's two lumens drawn as parallel segments -- has no single takeoff and is
    skipped, as name_anatomy skips cycles; otherwise segment order would decide.
    tests/test_proximal_end.py guards all three.
    """
    on_chain = set(chain)
    link = {}                         # junction -> junctions one off-aorta segment away
    for name, (ns, ne) in jg.items():
        if name in on_chain or ns is None or ne is None:
            continue
        link.setdefault(ns, []).append(ne)
        link.setdefault(ne, []).append(ns)
    out = set()
    for t in tree.terminals:
        if t in trunk:
            continue
        ns, ne = jg[t.segment.name]
        start = ns if ns is not None else ne
        if start is None:
            continue
        # Every aortic junction this branch's vessels reach without running along the
        # aorta. An aortic junction is an end point of the search, not a way through it.
        seen, stack, reached = {start}, [start], set()
        while stack:
            j = stack.pop()
            if j in stations:
                reached.add(j)
                continue
            for o in link.get(j, ()):
                if o not in seen:
                    seen.add(o)
                    stack.append(o)
        if len(reached) == 1:
            out.add(stations[reached.pop()])
    return sorted(out)


def _proximal_end(tree, trunk, jg, rad=None):
    """Which trunk end the blood enters by -- decided WITHOUT trusting the z axis.

    "Higher z is proximal" holds only while the export runs head-first with +z
    superior. A feet-first acquisition, or a pipeline that writes z the other way,
    silently swaps the aortic root with the distal descending end -- and every name
    that follows is then wrong. Two z-free signals settle it instead:

      * the great vessels all leave from the ASCENDING end, so the trunk end nearer
        the branch takeoffs (measured ALONG the aorta, each branch counted where its
        vessel tree joins it -- see _branch_takeoffs) is the proximal one;
      * failing that, the aorta TAPERS from root to periphery, so the fatter end is
        the proximal one.

    Only when neither is decisive does it fall back to z, and it says so.
    Returns (terminal, how) with how in {"single", "branches", "taper", "z"}.
    """
    if len(trunk) < 2:
        return (trunk[0] if trunk else None), "single"
    a, b = trunk[0], trunk[1]
    chain = _aorta_chain(tree, jg, a.segment.name, b.segment.name)
    stations = _aorta_stations(tree, jg, a.segment.name, b.segment.name)
    total = float(sum(tree.segments[n].length for n in chain))

    takeoffs = _branch_takeoffs(tree, jg, trunk, chain, stations)
    if takeoffs and total > 0:
        mean = float(np.mean(takeoffs))
        if abs(mean - 0.5 * total) > 0.05 * total:      # decisively one side
            return (a if mean < 0.5 * total else b), "branches"

    if rad:
        ra, rb = float(rad.get(a.label, 0.0)), float(rad.get(b.label, 0.0))
        if max(ra, rb) > 0.0 and abs(ra - rb) > 0.05 * max(ra, rb):
            return (a if ra > rb else b), "taper"

    return max(trunk, key=lambda t: float(t.tip[2])), "z"


class AortaRoute:
    """The aorta as ONE continuous ruler, from the inlet to its far end.

    Every clinically interesting region -- an aortic neck, an Ishimaru landing zone,
    the arch -- is a stretch of this ruler between two stations, so having a single
    arc-length coordinate that runs across segment boundaries is what lets those be
    expressed at all. Cuts up to now could only be placed relative to a vessel TIP,
    which cannot say "from the lowest renal to 15 mm below it".
    """

    def __init__(self, tree, chain, forward, starts, total):
        self.tree, self.chain = tree, chain
        self.forward, self.starts, self.length = forward, starts, total

    def locate(self, s):
        """Station along the route -> (segment name, arc length within it, forward)."""
        s = float(np.clip(s, 0.0, self.length))
        i = max(0, min(int(np.searchsorted(self.starts, s, "right")) - 1,
                       len(self.chain) - 1))
        seg = self.tree.segments[self.chain[i]]
        local = s - self.starts[i]
        return self.chain[i], (local if self.forward[i] else seg.length - local), \
            self.forward[i]

    def station_of(self, seg_name, s_local):
        """The inverse: in-segment arc length -> station along the route."""
        i = self.chain.index(seg_name)
        seg = self.tree.segments[seg_name]
        return self.starts[i] + (s_local if self.forward[i] else seg.length - s_local)

    def point(self, s):
        name, local, _ = self.locate(s)
        return self.tree.segments[name].at_arclength(local)[0]


def _abdominal_chain(tree, jg, root):
    """The aorta's own segments, inlet first, ending at the iliac bifurcation.

    The aorta ends at its bifurcation, and there is no aorta below that -- so of the
    divisions found on the trunk, the right one is the MOST DISTAL. Two independent
    ways of finding a division are tried and the further answer wins:

      structural  an interior segment with exactly two openings past it, which is a
                  common iliac splitting into internal and external. This is the
                  original rule and it is right on most scans.
      outflow     walking down the trunk, a junction where the outflow DIVIDES rather
                  than a branch merely leaving -- the same test that tells an
                  abdominal aorta from a thoracic one, see PAIR_AREA_FRACTION.

    Neither alone is enough. The structural rule takes the FIRST entry of a list
    ordered by fewest-openings then longest segment, so where the iliacs never
    bifurcate -- a stent, a femoro-femoral bypass, poor contrast, or the CT stopping
    short -- it lands on whatever else splits in two and stops there: JK20_1 got a
    47 mm ruler on a 206 mm aorta. The outflow rule misses in the other direction,
    because a fat visceral branch can carry a quarter of what is left: on JK5 the SMA
    reaches 0.324 against the whole infrarenal aorta, so taking the first crossing
    ends the ruler at 72 mm. Measured over 47 abdominal scans, the further of the two
    is right wherever either is, and neither over-reaches past the bifurcation.

    Where a scan genuinely has no division at all, the route runs to the far end of
    the trunk rather than refusing. Losing every region tool because a patient has a
    stent is not an acceptable answer.
    """
    def _structural():
        cands = [c for c in pre_bifurcation_candidates(tree)
                 if len(c["terminals"]) == 2]
        if not cands:
            return []
        ch = _aorta_chain(tree, jg, root.segment.name, cands[0]["segment"])
        return ch[:-1] if len(ch) > 1 else ch

    rad = {t.label: float(getattr(t.segment, "_radius_mm", None) or 0.0)
           for t in tree.terminals}
    if max(rad.values(), default=0.0) <= 0.0:
        # No calibre measured yet, so there is no outflow to compare: every junction
        # would score zero and the walk would follow CSV row order. The structural
        # rule needs no radii, so it answers alone rather than badly.
        return _structural()

    chain, splits = trunk_walk(tree, rad, jg, root)
    over = [i for i, (frac, _j, _s) in enumerate(splits)
            if frac >= PAIR_AREA_FRACTION]
    # The walk does not stop at the aortic bifurcation, it carries on into the larger
    # common iliac -- whose own internal/external split is frequently the MORE even of
    # the two. So the last division that still has a branching vessel beyond it is the
    # aortic one; an iliac dividing into two tips has a single opening each way.
    rich = [i for i in over
            if splits[i][2][0] is not None
            and sum(1 for t in tree.terminals
                    if t.segment.name in _through(jg, splits[i][1],
                                                  splits[i][2][0])) >= 2]
    end = (rich or over or [None])[-1]
    # splits[i] sits between chain[i] and chain[i + 1], so keeping through chain[end]
    # ends the ruler AT the division rather than one vessel past it.
    outflow = chain[:end + 1] if end is not None else chain

    def _length(ch):
        return sum(tree.segments[n].length for n in ch) if ch else -1.0

    structural = _structural()
    return max((structural, outflow), key=_length)


def aorta_route(tree, jg=None):
    """Build the aorta's arc-length ruler, inlet first. None if it cannot be traced."""
    jg = jg or _junction_graph(tree)
    root = tree_root_terminal(tree)
    if root is None:
        return None
    trunk = [t for t in tree.terminals if getattr(t, "is_trunk", False)]
    other = [t for t in trunk if t is not root]
    if other:                                   # thoracic: run to the far trunk end
        chain = _aorta_chain(tree, jg, root.segment.name, other[0].segment.name)
    else:                                       # abdominal: run down to the bifurcation
        chain = _abdominal_chain(tree, jg, root)
    if not chain:
        return None

    # Which way each segment is walked: entered by its start junction => forward.
    forward, starts, total = [], [], 0.0
    prev_j = None
    for i, name in enumerate(chain):
        ns, ne = jg[name]
        if i == 0:
            fwd = not root.at_end               # leave the inlet tip heading inward
        else:
            fwd = (ns == prev_j)
        forward.append(bool(fwd))
        starts.append(total)
        total += tree.segments[name].length
        prev_j = ne if fwd else ns
    return AortaRoute(tree, chain, forward, np.asarray(starts, float), total)


def aorta_landmarks(tree, route=None, names=None):
    """Named branch ostia as stations along the route: {name: station_mm}."""
    route = route or aorta_route(tree)
    if route is None:
        return {}
    jg = _junction_graph(tree)
    # station of each junction on the route
    jstation = {}
    for i, name in enumerate(route.chain):
        ns, ne = jg[name]
        seg = tree.segments[name]
        a, b = (ns, ne) if route.forward[i] else (ne, ns)
        if a is not None:
            jstation.setdefault(a, route.starts[i])
        if b is not None:
            jstation.setdefault(b, route.starts[i] + seg.length)
    out = {}
    for t in tree.terminals:
        if getattr(t, "is_trunk", False):
            continue
        ns, ne = jg[t.segment.name]
        j = ns if ns is not None else ne
        if j in jstation:
            out[(names or {}).get(t.label, t.name)] = float(jstation[j])
    return out


def region_terminals(tree, s_prox, s_dist, route=None, jg=None):
    """Two cuts bounding the stretch of aorta between two stations.

    The proximal one removes everything back toward the inlet, the distal one
    everything beyond -- so what survives is the region itself, with whatever
    branches leave from inside it still attached. This is the trim generalised: a
    trim removes ONE side of an interior point, and a region is simply two of them
    facing each other.
    """
    jg = jg or _junction_graph(tree)
    route = route or aorta_route(tree, jg)
    if route is None or s_dist - s_prox < 1e-6:
        return []
    out = []
    for s, keep_upstream in ((s_prox, False), (s_dist, True)):
        name, local, fwd = route.locate(s)
        seg = tree.segments[name]
        # at_end picks which side of the plane is discarded; forward traversal means
        # "toward the inlet" is the segment's start.
        at_end = (not fwd) if keep_upstream is False else fwd
        ns, ne = jg[name]
        far_j = ne if at_end else ns
        sub = _subtree_segments(jg, far_j, name) if far_j is not None else set()
        t = Terminal(seg, at_end, label=f"{name}.region{'B' if keep_upstream else 'A'}",
                     mode="normal",
                     name="Region start" if not keep_upstream else "Region end")
        t.is_interior = True
        # Say what this is ON THE OBJECT. A region boundary is not a vessel end: its
        # depth comes from _region_offset just below, and any caller that sets an
        # offset on it is overruled. vascular_api.offset_for_mark already guarded
        # against picking one with `getattr(x, "is_region_bound", False)` -- but
        # nothing ever set the attribute, so the guard was always False and inert.
        # "Send cut to mark" could therefore land on a boundary, set an offset that
        # _cut_list overwrites, and report success having moved nothing.
        #
        # It is set here rather than derived by the callers because this is the one
        # place a boundary is made; deriving it elsewhere is what produced two
        # different answers to the same question (see openings(), which has to
        # exclude trims separately because "not in tree.terminals" is a DIFFERENT
        # question and answering it that way disabled the slider on every trim).
        t.is_region_bound = True
        t.subtree_segments = sub
        t._region_offset = (seg.length - local) if at_end else local
        out.append(t)
    return out


def aorta_landmark_stations(tree, route=None, jg=None, names=None, vessels=None):
    """Every named ostium as a station along the aorta -- branches AND named trunks.

    aorta_landmarks only reports terminals, but the landmark that bounds an Ishimaru
    zone is often a vessel that was itself followed past a split (the brachiocephalic),
    so those have to be here too or the zones cannot be placed.
    """
    jg = jg or _junction_graph(tree)
    route = route or aorta_route(tree, jg)
    if route is None:
        return {}
    jstation = {}
    for i, name in enumerate(route.chain):
        ns, ne = jg[name]
        seg = tree.segments[name]
        a, b = (ns, ne) if route.forward[i] else (ne, ns)
        if a is not None:
            jstation.setdefault(a, route.starts[i])
        if b is not None:
            jstation.setdefault(b, route.starts[i] + seg.length)

    # WALKED OUT TO THE ROUTE, not required to sit on it.
    #
    # A vessel that leaves the aorta directly attaches at a junction on the route. One
    # that leaves THROUGH another vessel does not: a carotid springing from the
    # brachiocephalic has its junction on the brachiocephalic, which is nowhere on the
    # route. Tested with a bare `j in jstation` those simply vanish, and a landmark
    # that vanishes takes its region with it.
    #
    # _station_on_route already exists for exactly this and is already used when
    # deciding what a region CONTAINS -- it was written after "after isolating the arch
    # I can no longer trim the brachiocephalic trunk". Landmark discovery never got the
    # same treatment, so the two disagreed about where a vessel meets the aorta.
    #
    # Measured on KK28_11: only the brachiocephalic was found, at 84.4 mm of a 402 mm
    # route. The arch preset is built from the first and last landmark, so it came out
    # 30 mm long, collapsed once its ostia were cleared, and was dropped -- the case
    # offered no aortic arch at all. Every other thoracic case has at least three
    # landmarks and never noticed.
    out = {}
    for t in tree.terminals:
        if t.segment.name in route.chain:
            continue
        ns, ne = jg[t.segment.name]
        st = _station_on_route(jg, jstation, ns if ns is not None else ne)
        if st is not None:
            out[(names or {}).get(t.label, t.name)] = float(st)
    for seg_name, label in (vessels or {}).items():
        ns, ne = jg.get(seg_name, (None, None))
        for j in (ns, ne):
            st = _station_on_route(jg, jstation, j)
            if st is not None:
                out.setdefault(label, float(st))
                break
    return out


# Where an aortic neck ends: the first station whose calibre exceeds the reference
# at the lowest renal by this factor. 1.10 is the reporting convention -- the neck
# ends at the first slice more than 10% wider than baseline at the lowest renal.
# This was 1.20, which is more conservative than the literature and reports a neck
# longer than a reader following the standard would measure.
NECK_DILATION_FACTOR = 1.10


def _dilation_station(tree, mesh, route, s_from, factor=NECK_DILATION_FACTOR,
                      step=2.0):
    """Where the aorta first exceeds its calibre at *s_from* by *factor*.

    This is what ends an aortic neck: the neck runs from the lowest renal to where
    the vessel starts to balloon, so the distal boundary is a measurement, not a
    guess. Returns None if it never dilates.

    The reference is taken just distal to the lowest renal, matching the convention
    of baselining at that artery's distal edge.
    """
    nm, loc, _ = route.locate(min(s_from + 1.0, route.length))
    ref = section_diameter(mesh, *tree.segments[nm].at_arclength(loc))
    if not ref:
        return None
    for s in np.arange(s_from + step, route.length, step):
        nm, loc, _ = route.locate(s)
        sd = section_diameter(mesh, *tree.segments[nm].at_arclength(loc))
        if sd and sd["ead_mm"] > factor * ref["ead_mm"]:
            return float(s)
    return None


def centerline_arc_angle(tree, route, s0, s1, step_mm=1.0):
    """How much the centerline turns over a stretch, from the circle that fits it.

    theta = L / R, the angle a circle of radius R subtends over arc length L. It is
    the natural formulation for neck angulation: the answer does not depend on where
    someone chooses to draw two tangents, which is what makes hand-measured
    angulation disagree between observers.

    Returns None if the stretch is too short to fit. Otherwise a dict with:
      angle_deg   theta, in degrees
      radius_mm   the osculating circle's radius
      length_mm   centerline length used
      rms_mm      how far the points sit from that circle, on average
      planarity   out-of-plane spread over in-plane spread, 0 for a flat curve

    Read rms_mm and planarity before trusting angle_deg. A stretch that is not
    close to a circular arc has an R, and therefore an angle, that means little.
    """
    if route is None or s1 - s0 < 2.0:
        return None
    step = max(float(step_mm), 0.25)
    stations = list(np.arange(float(s0), float(s1), step))
    # arange stops short of the end, and the shortfall is measured length: at 1 mm
    # steps over a 31.4 mm arc it dropped 0.4 mm and under-reported the angle by
    # 1.3%. Always finish at s1.
    if not stations or stations[-1] < float(s1) - 1e-9:
        stations.append(float(s1))
    pts = []
    for s in stations:
        try:
            nm, loc, _ = route.locate(float(s))
            pos, _tan = tree.segments[nm].at_arclength(loc)
        except Exception:
            continue
        pts.append(np.asarray(pos, dtype=float))
    if len(pts) < 5:
        return None
    P = np.vstack(pts)

    # Arc length along the sampled points, not s1 - s0: the route's own stations
    # already measure that, and using them keeps the two consistent.
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    L = float(seg.sum())
    if L < 2.0:
        return None

    # Best-fit plane through the points; the circle is fitted inside it.
    c = P.mean(axis=0)
    U, S, Vt = np.linalg.svd(P - c, full_matrices=False)
    e1, e2, n = Vt[0], Vt[1], Vt[2]
    planarity = float(S[2] / S[0]) if S[0] > 0 else 0.0

    # A straight stretch has no circle. The fit is degenerate on collinear points
    # and returns whatever radius least squares lands on: on a dead-straight line it
    # produced R = 17.6 mm and an angle of 195 degrees, which is worse than useless
    # because it looks like an answer. Curvature shows up as spread across the
    # second axis, so test that first.
    straightness = float(S[1] / S[0]) if S[0] > 0 else 0.0
    if straightness < 1e-3:
        return {
            "angle_deg": 0.0, "radius_mm": None, "length_mm": round(L, 1),
            "rms_mm": 0.0, "planarity": round(planarity, 3),
            "n_points": int(len(P)), "straight": True,
            "centre": None, "normal": n.tolist(),
        }

    x = (P - c) @ e1
    y = (P - c) @ e2

    # Algebraic circle fit: x^2 + y^2 = 2ax + 2by + d, linear in (a, b, d).
    A = np.column_stack([2.0 * x, 2.0 * y, np.ones_like(x)])
    b = x ** 2 + y ** 2
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    ax, by, d = sol
    r2 = d + ax ** 2 + by ** 2
    if r2 <= 0:
        return None
    R = float(np.sqrt(r2))
    if not np.isfinite(R) or R <= 1e-6:
        return None

    rms = float(np.sqrt(np.mean(
        (np.sqrt((x - ax) ** 2 + (y - by) ** 2) - R) ** 2)))
    theta = L / R                                   # radians
    return {
        "angle_deg": round(float(np.degrees(theta)), 1),
        "radius_mm": round(R, 1),
        "length_mm": round(L, 1),
        "rms_mm": round(rms, 2),
        "planarity": round(planarity, 3),
        "n_points": int(len(P)),
        "straight": False,
        "centre": (c + ax * e1 + by * e2).tolist(),
        "normal": n.tolist(),
    }


# The shortest zone worth offering. Below this the two landmarks are effectively
# the same station -- a bovine arch puts the left carotid ON the brachiocephalic --
# and what would be isolated is a ring of triangles rather than a segment of aorta.
MIN_ZONE_SPAN_MM = 2.0


def region_presets(tree, mesh=None, names=None, vessels=None, planar=False):
    """Standard regions this particular anatomy can actually support.

    Each is {name, s0, s1, why}. Only offered when its landmarks exist -- an arch
    with no left common carotid cannot bound Ishimaru zones 1 and 2, and saying so
    by omission is better than placing a plane at a vessel that was never segmented.
    """
    route = aorta_route(tree)
    if route is None:
        return []
    land = aorta_landmark_stations(tree, route, names=names, vessels=vessels)
    L = route.length
    out = []

    def find(*words):
        hits = [(k, v) for k, v in land.items()
                if any(w in k.lower() for w in words)]
        return hits

    trunk = [t for t in tree.terminals if getattr(t, "is_trunk", False)]
    if len(trunk) == 1:                                   # ---- abdominal ----
        renals = find("renal")
        if renals:
            lo = max(v for _, v in renals)
            end = (_dilation_station(tree, mesh, route, lo) if mesh is not None
                   else None)
            out.append({"name": "Aortic neck", "s0": lo,
                        "s1": end if end is not None else min(lo + 40.0, L),
                        "why": (f"lowest renal to where the aorta first exceeds "
                                f"that calibre by "
                                f"{(NECK_DILATION_FACTOR - 1) * 100:.0f}%"
                                if end is not None else
                                "lowest renal, 40 mm distally (load a mesh to find "
                                "the dilation)")})
            out.append({"name": "Infrarenal aorta", "s0": lo, "s1": L,
                        "why": "lowest renal down to the bifurcation"})
            # s1=None: one plane, everything below it kept. The two-plane version
            # stops a millimetre short of the route end, and the iliac ostia sit AT
            # the bifurcation, so they fall outside it and cannot be kept or
            # trimmed. This keeps them as ordinary openings.
            out.append({"name": "Infrarenal aorta with iliacs", "s0": lo,
                        "s1": None,
                        "why": ("everything above the lowest renal removed, the "
                                "iliacs kept so they can be cut or trimmed as a "
                                "group")})
        if find("celiac", "sma", "mesenteric"):
            top = min(v for _, v in find("celiac", "sma", "mesenteric"))
            out.append({"name": "Visceral segment", "s0": max(top - 15.0, 0.0),
                        "s1": (max(v for _, v in renals) + 15.0) if renals else L,
                        "why": "around the celiac / SMA / renal group"})
        out.append({"name": "Whole abdominal aorta", "s0": 0.0, "s1": L,
                    "why": "inlet to the bifurcation"})
    else:                                                 # ---- thoracic ----
        innom = find("brachiocephalic", "innominate")
        lcc = [(k, v) for k, v in find("common carotid") if "left" in k.lower()]
        lsa = [(k, v) for k, v in find("subclavian") if "left" in k.lower()]
        if innom:
            i = min(v for _, v in innom)
            out.append({"name": "Ishimaru zone 0", "s0": 0.0, "s1": i,
                        "why": "ascending aorta, proximal to the brachiocephalic"})
            if lcc:
                out.append({"name": "Ishimaru zone 1", "s0": i,
                            "s1": min(v for _, v in lcc),
                            "why": "brachiocephalic to left common carotid"})
        if lcc and lsa:
            out.append({"name": "Ishimaru zone 2", "s0": min(v for _, v in lcc),
                        "s1": min(v for _, v in lsa),
                        "why": "left common carotid to left subclavian"})
        if lsa:
            s = min(v for _, v in lsa)
            out.append({"name": "Ishimaru zone 3", "s0": s, "s1": min(s + 20.0, L),
                        "why": "the first 20 mm distal to the left subclavian"})
            out.append({"name": "Ishimaru zone 4", "s0": min(s + 20.0, L), "s1": L,
                        "why": "beyond 20 mm distal to the left subclavian"})
        if land:
            # The distal end is a CONVENTION (last great vessel + 20 mm), not a
            # clinical definition -- and on half of 20 real arches that convention
            # lands inside the hairpin, where there is no planar cross-section to
            # cut: the arch preset failed 17 times, 15 of them "rim not flat".
            # Pushing the bound further along the descending aorta rescues them
            # (measured: KY160_1 strays 8.6 mm at +20, 2.4 at +25, clean from +30;
            # 8 of the 18 offered cases recover this way). So carry alternatives and
            # let the caller take the first that actually cuts, saying how far it
            # had to go -- rather than silently redefining the arch for everyone.
            # ANCHORED ON THE LEFT SUBCLAVIAN, not on whichever landmark happens to
            # sit furthest along. The distal end of the arch is defined off the LSA
            # -- it is where Ishimaru zone 3 begins -- and "the last great vessel" is
            # only the same thing when the LSA is the last one.
            #
            # It often is not. On four of thirty thoracic cases the left common
            # carotid is mapped DISTAL to the subclavian, so "last + 20 mm" put the
            # arch's end 36 to 61 mm past the subclavian instead of 20: KK40_3 +41.3,
            # DJ8_1 +39.0, KY68_1 +21.9, KY7_1 +16.3 of overshoot, every one of them
            # traceable to the anchor rather than to the cut refusing.
            #
            # Whether a carotid really does arise distal to the subclavian on those
            # four, or whether they are named wrongly, is a question about the scans
            # and not about this rule. Either way the arch ends a fixed distance past
            # the SUBCLAVIAN, so that is what it is measured from, and the old
            # behaviour survives only where no subclavian was named at all.
            lsa_st = [v for k, v in land.items() if "left subclavian" in k.lower()]
            anchor = min(lsa_st) if lsa_st else max(land.values())
            # THE ASCENDING END, anchored on the BRACHIOCEPHALIC and not on whichever
            # landmark sorts first. The two agree on all 30 arches measured -- every
            # vessel arising from the trunk is walked out to where the TRUNK meets the
            # aorta, so a right carotid or a bovine left carotid reports the trunk's
            # station rather than one of its own -- but they agree by construction of
            # something else, and the distal end has already been broken once by an
            # anchor that was "whichever is last". Name the vessel the rule means.
            #
            # 10 mm proximal to it, because that is where a person put it: over 23
            # arches placed by hand the ascending bound sat a median 11.0 mm before
            # the trunk, IQR 1.05 mm wide, and every one cut with a flat rim. A fixed
            # distance also fits far better than a clearance scaled by the trunk's
            # calibre (CV 0.27 against 0.59) -- KR27_2 has a 20.9 mm-radius bovine
            # trunk, where three radii of clearance predicts 63 mm and the hand-placed
            # answer was 13.
            bct_st = [v for k, v in land.items()
                      if "brachioceph" in k.lower() or "innominate" in k.lower()]
            head = min(bct_st) if bct_st else min(land.values())
            # Target first, then outward but INSIDE +/-5 mm, and only then past it.
            # Measured: DJ8_1, KK20_2 and KY68_1 take 10, KK32_0 takes 12, KK40_3
            # takes 15, and KK28_11 and KP20_1 refuse every station in the band --
            # their ascending aorta has no planar cross-section there at all, which
            # is why the rungs continue past it rather than giving up.
            out.append({"name": "Aortic arch",
                        "s0": max(head - 10.0, 0.0),
                        "s0_options": [max(head - m, 0.0)
                                       for m in (10, 11, 12, 9, 13, 8, 14, 7, 15,
                                                 6, 5, 16, 18, 20, 22, 25, 28, 30)],
                        "s1": min(anchor + 20.0, L),
                        # 20 first, then the rest of the target band, and only then
                        # out into the hairpin. A case that can cut at 25 or 15 should
                        # not be pushed to 30 because 20 happened to refuse.
                        "s1_options": [min(anchor + m, L)
                                       for m in (20, 25, 15, 22, 18, 30, 35, 40,
                                                 50, 60)
                                       if anchor + m <= L + 1e-9 or m == 20],
                        "why": "across all the great vessels, from 10 mm before the "
                               "brachiocephalic to 20 mm past the left subclavian"})
    # WHICH ZONES ARE REAL, and separately, which ones a PLANE could have cut.
    #
    # This filter used to plan each zone with planar boundaries and drop the ones whose
    # ends did not survive being nudged clear of every ostium. That is a sound test for
    # a planar cut and a wrong one for a zone: Ishimaru zone 1 is the gap between the
    # brachiocephalic and the left carotid, often 10-15 mm, and OSTIUM_CLEAR_R radii of
    # clearance at each end eats it. Measured over the 30-case thoracic cohort, the
    # filter was the binding constraint on almost every missing zone --
    #
    #     zone 1   offered  1/30   ->  18/30 without it
    #     zone 2   offered  8/30   ->  29/30
    #     zone 3   offered 26/30   ->  30/30
    #
    # -- and a station cut needs no clearance at all, because a branch is kept or
    # dropped whole and no boundary ever slices a mouth. So the default test is now
    # whether the ZONE exists: its landmarks were found and it spans something.
    #
    # The 12 cases that still cannot offer zone 1 are anatomy and not tooling: 8 are
    # the structurally detected bovine arches, where the left carotid arises from the
    # brachiocephalic so the zone has no extent, and on the other 4 it arises within
    # MIN_ZONE_SPAN_MM of it.
    #
    # planar=True restores the old behaviour for a caller that really is going to cut
    # with planes and would rather not be offered a zone it cannot make.
    jg = _junction_graph(tree)
    keep = []
    for p in out:
        open_distal = p["s1"] is None
        if not open_distal and p["s1"] - p["s0"] <= MIN_ZONE_SPAN_MM:
            continue
        if not planar:
            keep.append(p)
            continue
        try:
            bounds, _inside, _removed, _notes = region_plan(
                tree, p["s0"], p["s1"], route, jg, mesh)
        except Exception:
            keep.append(p)      # could not plan it -- not evidence it is impossible
            continue
        # Any bound at all is a region that can be cut. It used to demand exactly
        # two (or one for an open distal end), which silently discarded every region
        # whose boundary IS the vessel's own end -- the whole abdominal aorta,
        # Ishimaru zone 4 -- once those stopped placing a plane there. That dropped
        # 33 presets from the abdominal cohort while looking like an improved pass
        # rate, because the ones removed were the ones that had been failing.
        #
        # Zero bounds is still refused: that is a span with nothing to cut at either
        # end, which is the whole vessel and not a region.
        if bounds:
            keep.append(p)
    return keep


# NOT capped. Capping the ostium clearance at 12 mm was tried and measured worse:
# on 18 thoracic arches it cut the cases whose 'Aortic arch' region could be made
# to cut from 16 to 10, because the proximal bound stopped being pushed clear of
# the brachiocephalic mouth and tore there instead -- and unlike a distal bound,
# a torn START cannot be rescued by extending the region. The wide clearance a fat
# trunk claims is doing real work on the PRESET path. What a user drags is handled
# separately: region_plan(clear_ostia=False) leaves their boundary alone entirely.


def _ostium_radius(seg, mesh):
    """Radius used to size an ostium's clearance, measured rather than guessed.

    classify_terminals only caches _radius_mm on segments that carry a TERMINAL, so
    interior vessels (a brachiocephalic) have none and this used to fall back to a
    flat 2.0 mm. That made the clearance -- and therefore the region's boundaries --
    depend on whether pre_bifurcation_candidates happened to run first and fill the
    cache as a side effect: measured on DJ8_1 and KK16_1, Ishimaru zone 1 planned
    two clean bounds from a fresh load and collapsed entirely once that function had
    run. Same case, same region, different answer by call order. Measure it here
    when a mesh is available so the plan does not depend on history.
    """
    if getattr(seg, "_radius_mm", None) is None and mesh is not None:
        try:
            seg._radius_mm = segment_radius(mesh, seg)
        except Exception:
            pass
    return max(float(getattr(seg, "_radius_mm", None) or 2.0), 1.0)


# A direction change past this is the centerline doubling back, not a bend. An aortic
# arch turns through 180 degrees over its length and never approaches this between two
# consecutive points.
CENTERLINE_REVERSE_DEG = 90.0
# Two points this close together but this far apart ALONG the line means the path
# revisits itself -- it is lying on top of itself rather than following a vessel.
CENTERLINE_TOUCH_MM = 1.0
CENTERLINE_TOUCH_GAP = 10


def centerline_defects(tree):
    """Segments whose centerline folds back on itself: {name: (reversals, touches)}.

    Reported, never repaired. A path that reverses direction or lies on top of itself
    is not describing a vessel, and everything measured from it -- length, radius,
    every cut plane's normal -- is measured from the wrong thing. On DJ18_2 the
    descending aorta reverses four times and touches itself at fourteen places, which
    is what makes it draw as a thick tangle; the same case is the one thoracic aortic
    end that has to give up a horizontal cut plane as too oblique.

    Rare: 2 of 79 cases across the two lab cohorts, both on pieces the loader had to
    split, and no abdominal case at all. Rare enough to be worth naming rather than
    working around, since a user handed one has no other way to know.
    """
    from scipy.spatial import cKDTree

    out = {}
    for name, seg in tree.segments.items():
        pts = np.asarray(seg.pts, float)
        if len(pts) < 12:
            continue
        d = np.diff(pts, axis=0)
        n = np.linalg.norm(d, axis=1, keepdims=True)
        d = d / np.maximum(n, 1e-9)
        turn = np.degrees(np.arccos(np.clip((d[:-1] * d[1:]).sum(axis=1), -1.0, 1.0)))
        rev = int((turn > CENTERLINE_REVERSE_DEG).sum())
        touch = sum(1 for i, j in cKDTree(pts).query_pairs(CENTERLINE_TOUCH_MM)
                    if abs(i - j) > CENTERLINE_TOUCH_GAP)
        if rev or touch:
            out[name] = (rev, touch)
    return out


def _station_on_route(jg, jstation, junction):
    """Where a vessel hanging off *junction* meets the AORTA, in route mm.

    A branch that leaves the aorta directly attaches at a junction on the route and
    this returns it immediately. One that leaves through another vessel does not: the
    right common carotid springs from the brachiocephalic, so its own junction sits on
    the brachiocephalic and is nowhere on the route.

    Those were silently dropped, and the consequence is not small. Every great vessel
    off the brachiocephalic fell out of an isolated region: KY160's "Aortic arch" kept
    only the left subclavian, and the left carotid, right carotid and right subclavian
    -- the vessels the region exists to contain -- were neither inside it nor cut.
    The desktop tool then hid the brachiocephalic trim, because a trim is only offered
    when its whole group is in view. Reported as "after isolating the arch I can no
    longer trim the brachiocephalic trunk".

    So walk outward until the route is reached. It is a tree, so the first junction
    found on the route is the only one, and it is where this vessel's GROUP leaves the
    aorta -- which is the station a region boundary has to respect anyway.
    """
    if junction is None:
        return None
    if junction in jstation:
        return float(jstation[junction])
    adj = {}
    for name, (a, b) in jg.items():
        if a is not None:
            adj.setdefault(a, []).append((name, b))
        if b is not None:
            adj.setdefault(b, []).append((name, a))
    seen, queue = {junction}, [junction]
    while queue:
        j = queue.pop(0)
        for _name, other in adj.get(j, ()):
            if other is None or other in seen:
                continue
            if other in jstation:
                return float(jstation[other])
            seen.add(other)
            queue.append(other)
    return None


def region_plan(tree, s_prox, s_dist, route=None, jg=None, mesh=None,
                clear_ostia=True, clear_start=True):
    """Everything needed to isolate the aorta between two stations.

    Returns (bounds, inside, removed, notes):
      bounds  -- the two cuts facing each other across the region
      inside  -- branches that leave from WITHIN it, which stay and keep their own
                 cuts (an isolated arch should still have its great vessels)
      removed -- the segment names the bounds discard
      notes   -- any boundary that had to move to clear a branch ostium

    Membership is decided by each branch's OSTIUM STATION along the route, not by
    which segments survive: the aorta's own end terminals sit on the same segments
    the bounds cut, so a set test would re-cut a vessel the boundary has already
    taken, and the second cut then finds nothing left to remove.

    clear_start=False keeps the clearance on the DISTAL bound and takes it off the
    proximal one. That combination exists because the two bounds answer to different
    evidence. Clearance is worth its cost distally -- capping it was measured to take
    the arches that cut from 16 to 10. Proximally it is over-generous: on 30 arches it
    pushed the ascending bound to a median 17.5 mm before the brachiocephalic, and as
    far as 42 mm, where a person placing the same bound by hand chose 11 mm on 23 of
    them and every one of those cut with a dead-flat rim. So a caller that intends to
    walk its own ladder of proximal stations turns it off, tries them, and falls back
    to the cleared position only if none of them cut.
    """
    jg = jg or _junction_graph(tree)
    route = route or aorta_route(tree, jg)
    if route is None:
        return [], [], set(), []

    # station of every junction on the route, so a branch can be placed along it
    jstation = {}
    for i, name in enumerate(route.chain):
        ns, ne = jg[name]
        seg = tree.segments[name]
        a, b = (ns, ne) if route.forward[i] else (ne, ns)
        if a is not None:
            jstation.setdefault(a, route.starts[i])
        if b is not None:
            jstation.setdefault(b, route.starts[i] + seg.length)

    # Every vessel that leaves the aorta, whether it ends in an opening or goes on to
    # split. A brachiocephalic trunk carries no terminal of its own, so a terminals-only
    # list misses its ostium entirely -- and an Ishimaru zone 0 boundary is placed
    # exactly there, which tore the opening until this included it.
    ostia = []                       # (owner, station, clearance)
    seen_seg = set()
    for t in tree.terminals:
        if t.segment.name in route.chain:
            continue        # an END of the aorta itself, not a branch off it: the
                            # region's own bounds already cut that vessel, and cutting
                            # it twice leaves the second cut nothing to remove
        ns, ne = jg[t.segment.name]
        j = ns if ns is not None else ne
        st = _station_on_route(jg, jstation, j)
        if st is not None:
            ostia.append((t, st, OSTIUM_CLEAR_R * _ostium_radius(t.segment, mesh)))
            seen_seg.add(t.segment.name)
    for seg_name, (ns, ne) in jg.items():
        if seg_name in route.chain or seg_name in seen_seg:
            continue
        for j in (ns, ne):
            if j in jstation:
                seg = tree.segments[seg_name]
                ostia.append((seg, float(jstation[j]),
                              OSTIUM_CLEAR_R * _ostium_radius(seg, mesh)))
                break

    # A boundary must not land IN an ostium. The mouth is an ellipse several mm long
    # on the aorta wall -- a renal on a 25 mm aorta measures ~8 mm -- so a plane there
    # slices the mouth and leaves a ragged hole that is neither a clean opening nor a
    # sealed wall. Nudge the boundary just clear, INTO the region, and say so: the
    # region is a few mm shorter than asked, which is visible and correctable, whereas
    # a torn ostium is neither.
    # A boundary is nudged to whichever clear side it was ALREADY on, so the region
    # keeps the branches it was meant to contain: an arch boundary placed just before
    # the subclavian moves further before it, not past it.
    # A boundary AT the very tip has nothing beyond it to remove, so keep both ends
    # a little inside the vessel -- the same minimum a normal cut respects.
    # s_dist=None means "keep everything distal": one plane, and the vessels' own
    # ends downstream of it stay as ordinary openings. A two-plane region cannot
    # express that, because its distal plane is clamped a millimetre short of the
    # route end -- so on an abdominal aorta the iliacs, whose ostia sit AT the
    # bifurcation, fell outside every region and could never be kept or trimmed.
    open_distal = s_dist is None
    lo, hi = MIN_OFFSET_MM, route.length - MIN_OFFSET_MM
    notes = []
    s0 = float(np.clip(s_prox, lo, hi))
    s1 = float("inf") if open_distal else float(np.clip(s_dist, lo, hi))

    eps = 1e-6                               # a nudged boundary sits ON the clear edge
    if not clear_ostia:
        # The user placed this boundary. Moving it silently is what makes the slider
        # feel broken: on KY160_1 a bovine trunk (r 12.7 mm, so 3 radii = 38 mm of
        # clearance) pinned 62 mm of slider travel onto two positions. Leave the
        # boundary where it was put and say what it is near, so the choice is theirs.
        for t, st, clear in ostia:
            for tag, s in ((("start", s0),) if open_distal
                           else (("start", s0), ("end", s1))):
                if abs(s - st) < clear - eps:
                    notes.append(
                        f"{tag} sits {abs(s - st):.1f} mm from the {t.name} ostium, "
                        f"a cut through a branch mouth can tear; move it if the rim "
                        f"comes out ragged")
    else:
        for _ in range(4):                   # a nudge can bring another ostium into range
            moved = False
            for t, st, clear in ostia:
                if clear_start and abs(s0 - st) < clear - eps:
                    s0 = (st + clear) if s0 >= st else (st - clear)
                    notes.append(f"start moved to {s0:.1f} mm to clear the {t.name} ostium")
                    moved = True
                if not open_distal and abs(s1 - st) < clear - eps:
                    s1 = (st + clear) if s1 > st else (st - clear)
                    notes.append(f"end moved to {s1:.1f} mm to clear the {t.name} ostium")
                    moved = True
            if not moved:
                break
    notes = list(dict.fromkeys(notes))
    s0 = float(np.clip(s0, lo, hi))
    if not open_distal:
        s1 = float(np.clip(s1, lo, hi))
        if s1 - s0 < 1.0:
            return [], [], set(), ["the region collapses once the ostia are cleared"]

    if open_distal:
        # One plane. region_terminals always builds a facing pair, so ask for the
        # pair and keep only the proximal one; index 0 is the cut that removes
        # everything back toward the inlet.
        pair = region_terminals(tree, s0, hi, route, jg)
        bounds = pair[:1]
    else:
        bounds = region_terminals(tree, s0, s1, route, jg)
    if len(bounds) != (1 if open_distal else 2):
        return [], [], set(), notes

    # A boundary that lands on the vessel's OWN END gets no boundary plane. The
    # region does not stop there because something was cut -- it stops because the
    # vessel does. A plane a millimetre in is trying to slice a wafer off the end of
    # the model, and preflight refuses it: "a cut of 1.0 mm sits on the vessel tip".
    #
    # This is not rare. "Whole abdominal aorta" runs inlet to bifurcation, so its
    # proximal bound sits exactly on the inlet; Ishimaru zone 4 runs to the distal
    # end. Across both cohorts 101 presets place a bound within 2 mm of a vessel end,
    # and dropping that plane fixes 20 of them while changing nothing on the other 79.
    #
    # The END STILL GETS CUT, by its own terminal at its own offset. The input mesh is
    # closed, so an uncut vessel end is a CAP, not an opening -- dropping the plane
    # and stopping there left "Whole abdominal aorta" with 4 openings where it had 5,
    # sealing the inlet. That is the same trade the tool makes everywhere else: the
    # aorta's ends are cut as ends, not as region boundaries.
    # Decided on the STATION along the route, never on the in-segment offset.
    #
    # _region_offset is a distance within one SEGMENT, and it goes to zero whenever a
    # boundary lands near a junction -- which happens constantly in the middle of an
    # aorta, since every branch takeoff is one. Reading that as "the vessel's own end"
    # dropped the plane there and left nothing to remove the far side. On JK1_1 the
    # Infrarenal aorta start sits exactly on the renal junction, so its plane vanished
    # and the region kept the whole aorta ABOVE it -- the half it exists to discard --
    # while a 2 mm nudge moved the boundary off the junction and flipped it back.
    # Reported from the tool as "decrease the start by 1 mm and it grabs the wrong
    # segment of the anatomy", which is exactly what it was doing.
    on_end = []
    if bounds:
        if s0 <= lo + 1e-6:
            on_end.append(bounds[0])
        if not open_distal and len(bounds) > 1 and s1 >= hi - 1e-6:
            on_end.append(bounds[-1])
    if on_end:
        for t in on_end:
            side = "start" if t is bounds[0] else "end"
            notes.append(f"region {side} is the vessel's own end, so no plane is cut "
                         f"there -- it is already an opening")
        bounds = [t for t in bounds if t not in on_end]
    # "Automatic" has to mean measured. A region boundary sits on an INTERIOR aortic
    # segment, and classify_terminals only caches a radius for the segments that carry
    # terminals -- so band_width found none here and fell through to its own r = 5.0
    # default, handing every region boundary in the cohort a 12.50 mm reach whatever
    # the aorta actually measured. True radii here run 14.1 to 25.6 mm, so the plane
    # was being asked to span a 51 mm ring with a 25 mm disc, which is the "opening
    # comes out N mm out of plane" refusal.
    #
    # Worse than wrong, it was ORDER DEPENDENT: anything that happened to call
    # pre_bifurcation_candidates(tree, mesh) first filled the cache, and the same
    # region then cut with a 35 mm band instead. Same input, two different cuts.
    if mesh is not None:
        for t in bounds:
            if getattr(t.segment, "_radius_mm", None) is None:
                t.segment._radius_mm = segment_radius(mesh, t.segment)
    # NB: the band is otherwise left alone. Clamping it to the distance to the
    # nearest ostium looks right but is not: a boundary nudged clear of a branch sits
    # exactly that clearance away, and the clearance is sized by the BRANCH's radius,
    # so the clamp shrinks the band to a branch-sized disc on a fat aorta and tears
    # every opening. Measured: it broke all four working JK5 regions.
    # The aorta is ONE vessel, however many segments the centerline broke it into, and
    # a boundary cutting across it may clip any of them. Ownership otherwise stops at
    # the segment the plane happens to sit on, which tears the far side of a wide
    # vessel where the neighbouring piece of the same route claims it.
    for t in bounds:
        t._own_segments = set(route.chain)
    removed = set().union(*(t.subtree_segments for t in bounds)) if bounds else set()
    # Only TERMINALS become cuts of their own; the ostia list also carries interior
    # vessels (a brachiocephalic), which are there to keep a boundary clear of them,
    # not to be cut.
    #
    # A dropped boundary widens membership on that side: with no plane at the
    # proximal end, everything from the inlet down is inside the region, because
    # nothing removes it.
    lo_s, hi_s = s0, s1
    for t in on_end:
        if getattr(t, "name", "") == "Region start":
            lo_s = -np.inf
        else:
            hi_s = np.inf
    inside = [t for t, st, _ in ostia
              if isinstance(t, Terminal) and lo_s < st < hi_s]
    # With a boundary dropped, membership stops being a station test ENTIRELY --
    # exactly the reasoning the open-distal case already uses below. Everything the
    # remaining cut does not discard survives and needs its own opening, including
    # the vessel's own end terminal, which is otherwise left as a cap on a closed
    # mesh.
    #
    # A station test was tried here and is WRONG. It needs a station for every
    # terminal, and a terminal hanging off a junction the route ends at rather than
    # passes through has none -- so it took a default, and a permissive default let
    # every off-route terminal in. On JK1_1 "Infrarenal aorta" then kept the inlet,
    # the renals and the branch above them: the stretch ABOVE the boundary, which is
    # the half the region exists to remove. Nudging the start 2 mm flipped it back to
    # the iliacs, so the region appeared to grab a different piece of anatomy for a
    # millimetre of slider travel. Reported from the tool, not from a benchmark.
    if on_end:
        bsegs = {t.segment.name for t in bounds}
        inside = [t for t in tree.terminals
                  if t.segment.name not in removed and t.segment.name not in bsegs]
    if open_distal and bounds:
        # With one plane, membership is not a station test at all: every terminal
        # the cut does not discard survives and needs its own opening. A station
        # test misses them, because a route that runs down one iliac puts that
        # iliac's terminal ON the chain, and the others hang off a junction the
        # route ends at rather than passes through.
        bseg = bounds[0].segment.name
        inside = [t for t in tree.terminals
                  if t.segment.name not in removed and t.segment.name != bseg]
    return bounds, inside, removed, notes


def _leaves_along(tree, seg, j_pos, span_mm=12.0):
    """Unit direction a segment heads as it leaves the junction at *j_pos*."""
    a, b = seg.pts[0], seg.pts[-1]
    from_start = np.linalg.norm(a - j_pos) <= np.linalg.norm(b - j_pos)
    s = min(span_mm, 0.5 * seg.length)
    far = seg.at_arclength(s if from_start else seg.length - s)[0]
    return _unit_np(np.asarray(far, float) - np.asarray(j_pos, float))


def name_anatomy(tree, mesh=None):
    """Name each opening (and trimmable vessel) after the artery it actually is.

    Only where the geometry says so unambiguously. A branch whose signature is not
    decisive keeps its Branch N label: a guess that reads like a fact is worse than
    no name, and anatomical variants are common (accessory renals ~30%, replaced
    hepatics ~15-20%). Everything here is measured -- takeoff level, the direction
    the vessel leaves in, calibre, length, and whether it splits -- against the
    frame in ANATOMY_FRAME.

    Returns (terminal_names {label: name}, vessel_names {segment: name}, notes).
    """
    names, vessels, notes = {}, {}, []
    # Labels named by WHERE THEY BRANCH rather than by what they look like -- the
    # brachiocephalic trunk's own daughters. Only the arch fills this; an abdominal
    # aorta leaves it empty, which is why it has to exist before that branch of the
    # code decides whether to run at all.
    claimed = set()
    root = tree_root_terminal(tree)
    if root is None or len(tree.terminals) < 2:
        return names, vessels, notes

    def _unique(names, vessels, notes, protected=()):
        """Drop any artery name that matched more than one vessel.

        Every rule here is a pattern match, and a body carries plenty of vessels that
        resemble a great one: a vertebral looks like a carotid, an intercostal like a
        subclavian, a lumbar pair like the renals. When a name lands twice, neither
        match was distinctive, so BOTH are withdrawn -- a confident wrong label is the
        one outcome worse than no label.

        EXCEPT WHERE ONE OF THEM IS NOT A PATTERN MATCH. A vessel named because it
        arises FROM THE BRACHIOCEPHALIC TRUNK is named by the tree's own branching, not
        by climbing more than it spreads, and that is a fact about the anatomy rather
        than a resemblance. Withdrawing it and keeping the guess had a specific cost:
        on DJ8_1, KK20_2 and KY7_1 the trunk's own carotid daughter lost "Right common
        carotid" to a collision, the four-vessel elimination pass then handed it "Left
        common carotid", and the left carotid came out sharing the right subclavian's
        station -- so zone 1 had zero length and those three cases lost it entirely.
        """
        prot = set(protected)
        seen = {}
        for src, key in ((names, "opening"), (vessels, "vessel")):
            for lbl, nm in src.items():
                seen.setdefault(nm, []).append((src, lbl))
        for nm, hits in seen.items():
            if len(hits) <= 1:
                continue
            keep = [(src, lbl) for src, lbl in hits if lbl in prot]
            if len(keep) == 1:
                for src, lbl in hits:
                    if (src, lbl) != keep[0]:
                        src.pop(lbl, None)
                notes.append(f"“{nm}” fitted {len(hits)} vessels; it was kept on the "
                             "one arising from the brachiocephalic trunk, which is "
                             "named by where it branches rather than by its shape")
                continue
            for src, lbl in hits:
                src.pop(lbl, None)
            notes.append(f"“{nm}” fitted {len(hits)} different vessels, none was "
                         "named, since the match was not distinctive")
        return names, vessels, notes
    jg = _junction_graph(tree)

    # Root the tree at the inlet: what lies downstream of each interior segment.
    down = {}
    for name, (ns, ne) in jg.items():
        if ns is None or ne is None:
            continue
        se, ss = _subtree_segments(jg, ne, name), _subtree_segments(jg, ns, name)
        in_e, in_s = root.segment.name in se, root.segment.name in ss
        if in_e == in_s:
            continue                                  # cycle / disconnected
        at_end = not in_e
        jidx = ne if at_end else ns
        upidx = ns if at_end else ne                  # the end facing the inlet
        sub = se if at_end else ss
        down[name] = {"jpos": np.asarray(tree.junctions[jidx][0], float),
                      "upos": np.asarray(tree.junctions[upidx][0], float),
                      "upidx": upidx,
                      "sub": sub,
                      "terms": [t for t in tree.terminals if t.segment.name in sub]}

    def attach_idx(t):
        """Index of the junction where terminal *t*'s vessel joins the tree."""
        ns, ne = jg[t.segment.name]
        return ns if ns is not None else ne

    def attach(t):
        """Where terminal *t*'s vessel joins the tree."""
        j = attach_idx(t)
        return None if j is None else np.asarray(tree.junctions[j][0], float)

    trunk = [t for t in tree.terminals if getattr(t, "is_trunk", False)]
    branches = [t for t in tree.terminals if not getattr(t, "is_trunk", False)]

    # ---- thoracic: two trunk ends, an arch, great vessels leaving superiorly ----
    if len(trunk) == 2 and branches:
        # Which end is the root comes from the anatomy (_proximal_end), not from z.
        prox = next((t for t in trunk if getattr(t, "is_inlet", False)),
                    max(trunk, key=lambda t: t.tip[2]))
        dist = next(t for t in trunk if t is not prox)
        names[prox.label] = "Ascending aorta"
        names[dist.label] = "Descending thoracic aorta"

        # Midline and scale come from the AORTA, so a volume with a different origin
        # names its vessels the same way.
        aorta = _aorta_path(tree, jg, prox.segment.name, dist.segment.name)
        axis = np.vstack([tree.segments[n].pts for n in aorta])
        mid_x = float(np.median(axis[:, 0]))
        R = max(float(getattr(prox.segment, "_radius_mm", None) or 0.0), 1.0)
        # "Up" the body = back along the aorta toward the root. Used instead of +z so
        # a carotid still reads as ASCENDING in a z-flipped export.
        far_pt = axis[int(np.argmax(np.linalg.norm(
            axis - np.asarray(prox.tip, float), axis=1)))]
        up = _unit_np(np.asarray(prox.tip, float) - far_pt)

        # The brachiocephalic is the FIRST vessel off the arch, and it splits: two
        # openings (right subclavian + right common carotid), or three in a BOVINE
        # arch (~15%) where the left common carotid shares its origin.
        #
        # "First" has to mean first ALONG THE AORTA, not highest. On a type III arch
        # (~25-30% of people) the innominate origin sits BELOW the more distal arch
        # vessels -- KY160's own arch is already inverted that way. Electing by height
        # there crowns the left subclavian instead, and since the brachiocephalic is
        # treated as right-sided BY DEFINITION, that silently inverts left and right
        # for the whole arch.
        stations = _aorta_stations(tree, jg, prox.segment.name, dist.segment.name)
        cand = [n for n, d in down.items()
                if 2 <= len(d["terms"]) <= 3
                and not any(x.is_trunk for x in d["terms"])
                and d.get("upidx") in stations]
        innom = []
        if cand:
            cand.sort(key=lambda n: stations[down[n]["upidx"]])
            first = cand[0]
            # It must also be the first branch off the arch FULL STOP -- if some other
            # vessel leaves earlier, this splitting one is a subdivided distal vessel
            # (a subclavian followed past its vertebral, a carotid past its
            # bifurcation), not the brachiocephalic.
            here = stations[down[first]["upidx"]]
            others = [stations[j] for j in
                      (attach_idx(t) for t in branches) if j in stations]
            earliest_other = min(others, default=float("inf"))
            next_cand = min([stations[down[n]["upidx"]] for n in cand[1:]],
                            default=float("inf"))
            if here > earliest_other + 1e-6:
                notes.append("the first vessel off the arch does not split, so no "
                             "brachiocephalic trunk could be confirmed, sides were "
                             "not read from it")
            elif next_cand - here < 0.25 * R:
                notes.append("two splitting vessels leave the arch at the same point, "
                             "the brachiocephalic could not be identified with "
                             "confidence, so sides were not assumed")
            else:
                innom = [first]

        claimed, left_sign, bovine = set(), None, False
        if len(innom) == 1:
            kids = down[innom[0]]["terms"]
            # This vessel is RIGHT-sided by definition, so whichever way its branches
            # point IS the patient's right. That makes the arch its own handedness
            # witness -- the sides no longer rest on trusting the file's orientation.
            subclav = max(kids, key=lambda t: abs(float(t.tip[0]) - mid_x))
            rights = [subclav] + [t for t in kids
                                  if t is not subclav and float(t.tip[2]) > float(prox.tip[2])]
            bovine = len(kids) == 3
            vessels[innom[0]] = ("Brachiocephalic (innominate) trunk"
                                 + (", bovine arch" if bovine else ""))
            dx_sub = float(subclav.tip[0]) - mid_x
            if abs(dx_sub) > SIDE_TOL_R * R:
                left_sign = -1.0 if dx_sub > 0 else 1.0     # subclavian side is RIGHT
            names[subclav.label] = "Right subclavian"
            claimed.add(subclav.label)
            carotids = sorted((t for t in kids if t is not subclav),
                              key=lambda t: (float(t.tip[0]) - mid_x))
            if bovine and len(carotids) == 2 and left_sign is not None:
                # Two carotids off one trunk: the one on the patient's LEFT is the LCC.
                order = carotids if left_sign > 0 else carotids[::-1]
                names[order[0].label] = "Right common carotid"
                names[order[1].label] = "Left common carotid"
                claimed |= {c.label for c in carotids}
                notes.append("bovine arch: the left common carotid shares the "
                             "brachiocephalic origin (a normal variant, ~15%)")
            elif len(carotids) == 1:
                names[carotids[0].label] = "Right common carotid"
                claimed.add(carotids[0].label)

        if left_sign is None:
            left_sign = 1.0            # fall back to the frame in ANATOMY_FRAME
            notes.append("no brachiocephalic trunk found to confirm left from right, "
                         + SIDES_ASSUMED_NOTE)

        # The remaining arch vessels. Two things separate them, and neither depends on
        # where the volume's origin sits:
        #   a carotid ASCENDS to the neck  -- it climbs much more than it spreads;
        #   a subclavian REACHES OUT to the arm -- it spreads much more than it climbs.
        # Both are differences from the vessel's own takeoff, so they survive any
        # translation. Side comes from the coordinates when they are decisive, and
        # otherwise from anatomy: everything on the arch past the brachiocephalic is
        # left-sided (an aberrant right subclavian would be the exception, and is
        # rare enough to leave to a hand rename).
        def arch_side(far):
            dx = (float(far[0]) - mid_x) * left_sign
            return ("Left" if dx > SIDE_TOL_R * R else
                    ("Right" if dx < -SIDE_TOL_R * R else "Left"))

        def arch_shape(j, far):
            """(rise, lateral) for a vessel leaving the aorta at *j* toward *far*."""
            return (float(np.dot(np.asarray(far, float) - np.asarray(j, float), up)),
                    abs(float(far[0]) - float(j[0])))

        def arch_name(j, far):
            """Great-vessel name for something leaving the aorta at *j* toward *far*."""
            rise, lateral = arch_shape(j, far)
            side = arch_side(far)
            if lateral > rise and lateral > SUBCLAVIAN_MIN_R * R:
                return f"{side} subclavian"
            if rise > lateral and rise > CAROTID_MAX_R * R:
                return f"{side} common carotid"
            return None

        # Only vessels arising DIRECTLY from the aorta are great vessels. Without this
        # a carotid followed past its bifurcation names its internal AND external
        # branches "common carotid", and a subclavian followed past the vertebral does
        # the same -- duplicate names for vessels that are not what they claim.
        arch_pool = []
        for t in branches:
            if t.label in claimed or attach_idx(t) not in stations:
                continue
            nm = arch_name(attach(t), np.asarray(t.tip, float))
            if nm:
                names[t.label] = nm
                arch_pool.append((t, nm))

        # When a great-vessel name lands twice, decide between them rather than
        # withdrawing both. _unique's blanket withdrawal is right where a match is
        # merely undistinctive, and wrong here, because the two candidates ARE the
        # two vessels and only their labels are swapped: on DJ8_1 three vessels score
        # "carotid", two land on "Right common carotid", and both were dropped.
        #
        # A carotid ascends to the head and a subclavian reaches out to the arm and
        # finishes LOWER -- so on the same side the higher of the pair is the carotid
        # and the lower is the subclavian. That is the lab's own reading of these
        # cases, and it is the one signal that separates them when the absolute
        # thresholds cannot: those are scaled to the AORTIC radius, which on a 21 mm
        # aorta asks 53 mm of lateral reach before a subclavian is called at all.
        by_name = {}
        for t, nm in arch_pool:
            by_name.setdefault(nm, []).append(t)
        for nm, group in by_name.items():
            if len(group) < 2:
                continue
            side = nm.split()[0]
            ranked = sorted(
                group, key=lambda t: -arch_shape(attach(t), np.asarray(t.tip, float))[0])
            names[ranked[0].label] = f"{side} common carotid"
            names[ranked[-1].label] = f"{side} subclavian"
            for t in ranked[1:-1]:
                names.pop(t.label, None)       # a third vessel here is neither
            notes.append(
                f"two vessels both fitted “{nm}”; the one climbing further was taken "
                f"as the carotid and the one finishing lower as the subclavian, since "
                f"a subclavian reaches for the arm rather than the head")

        # A great vessel that was itself followed past its own split is an interior
        # segment: name the VESSEL (so the trim menu reads "Left common carotid"),
        # not its children.
        for n, d in down.items():
            if n in innom or d["upidx"] not in stations or n in vessels:
                continue
            nm = arch_name(d["upos"], d["jpos"])
            if nm:
                vessels[n] = nm

        # Withdraw the undistinctive matches BEFORE filling, deliberately. _unique
        # drops a name that landed on two vessels, and until it has run those two
        # still look named -- so a fill that asks "which of the four names are
        # missing" sees none missing, does nothing, and the withdrawal then leaves
        # both openings as Branch N. DJ8_1 is exactly that: two vessels tied on
        # "Right common carotid", both dropped, and the fill never got a chance.
        names, vessels, notes = _unique(names, vessels, notes, claimed)

        # A FOUR-VESSEL ARCH IS FULLY DETERMINED, so fill in whatever is left.
        #
        # 29 of the 30 thoracic cases carry exactly four branch openings, and four
        # openings off an arch are the right and left common carotids and the right
        # and left subclavians -- there is no fifth thing they can be. Two a side,
        # and within a side the one reaching further out is the subclavian, because
        # it is going to the arm while the carotid goes to the head.
        #
        # Ranking beats a threshold here, and the measurements say so. Over the
        # vessels this namer is already confident about, lateral reach separates the
        # two cleanly -- 0.59 aortic radii for a carotid against 4.56 for a
        # subclavian -- while the absolute floors are scaled to the AORTA, so a 21 mm
        # aorta demands 53 mm of reach before a subclavian is named at all. Ranked by
        # lateral reach this agrees with 47 of 49 confident names; ranked by climb it
        # agrees with 55 of 94, because a side holding one vessel ranks it a carotid
        # whatever it is doing.
        #
        # Only UNNAMED openings are filled. The two disagreements are both on KP20_1
        # and look like the existing pair being swapped, but "looks like" is not
        # ground truth, and overwriting a confident name on that basis is the kind of
        # unvalidated tuning this file keeps warning about.
        if len(branches) == 4 and any(t.label not in names for t in branches):
            def _lateral(t):
                j = attach(t)
                return -1.0 if j is None else abs(float(t.tip[0]) - float(j[0]))

            # Work from the NAMES that are missing, not from a split of the vessels.
            # The four names are fixed, so whatever the confident rules already
            # placed, the rest of the openings take the rest of the names -- and a
            # side test is only needed for the ones still in doubt. Splitting the
            # vessels two and two instead needs the sides to come out even, and they
            # do not: DJ18_2 puts three tips right of the aorta's median x and one
            # left, and a case with one named vessel in each half fills neither.
            canon = ["Right subclavian", "Right common carotid",
                     "Left subclavian", "Left common carotid"]
            missing = [n for n in canon if n not in names.values()]
            spare = sorted((t for t in branches if t.label not in names),
                           key=_lateral, reverse=True)
            # Subclavians reach further out than carotids, so the widest-reaching
            # openings take the subclavian names that are still going.
            for want in sorted(missing, key=lambda n: "subclavian" not in n):
                if not spare:
                    break
                t = spare.pop(0) if "subclavian" in want else spare.pop()
                names[t.label] = want
                notes.append(
                    f"“{t.name}” named {want} by elimination: a four-vessel arch "
                    f"carries both carotids and both subclavians and nothing else, "
                    f"and of what was left this one reaches "
                    f"{'furthest out' if 'subclavian' in want else 'least far out'}")

        # A CAROTID CANNOT SIT DISTAL TO THE SUBCLAVIAN ON ITS OWN SIDE.
        #
        # The rule above separates them by shape -- a carotid ascends to the neck, a
        # subclavian reaches out to the arm -- measured from takeoff to TIP. That is
        # sound where both vessels were traced to a comparable extent and wrong where
        # they were not: a subclavian whose centerline stops early rises more than it
        # spreads and reads as a carotid, and a carotid followed far laterally reads
        # as a subclavian. Nothing in the shape test can see that it has produced an
        # arch that does not exist.
        #
        # Order along the aorta can. The great vessels leave in a fixed sequence, so
        # on either side the carotid is PROXIMAL to the subclavian; the reverse is not
        # a variant, it is a mistake. Measured over the thoracic cohort it fires on
        # four of thirty, and on every one the vessel it renames subclavian is both
        # the more distal AND the wider of the two -- KY68_1 r=4.58 against 2.40,
        # KK40_3 4.50 against 3.08 -- which is the independent check, since a
        # subclavian is the fatter vessel.
        #
        # It swaps rather than withdraws: the two candidates ARE the two vessels and
        # only their labels are the wrong way round.
        # LEFT ONLY. On the left the order is fixed -- the carotid leaves before the
        # subclavian and no common variant reverses it. On the RIGHT it is not: an
        # aberrant right subclavian arises from the distal arch, past the carotid and
        # often past the left subclavian too, which this rule would "correct" into
        # nonsense. The comment on arch_side already sets that expectation, and DJ8_1
        # shows why -- applied to the right side the swap named the 5.43 mm vessel a
        # carotid and the 3.23 mm one a subclavian, against their calibre.
        for side in ("Left",):
            car = next((t for t in branches
                        if names.get(t.label) == f"{side} common carotid"), None)
            sub = next((t for t in branches
                        if names.get(t.label) == f"{side} subclavian"), None)
            if car is None or sub is None:
                continue
            # STATIONS, not junction indices. attach_idx returns an index into
            # tree.junctions, whose order says nothing about position on the aorta --
            # comparing those swapped the cases that were already right (KY160_1,
            # SA21_1) and left the wrong ones alone. `stations` is the map that
            # actually answers "how far along the aorta does this junction sit".
            # WALKED OUT TO THE AORTA, because a great vessel may leave through
            # another one. On three of the four the vessel named subclavian hangs off
            # the BRACHIOCEPHALIC, so its junction is not an aortic junction at all
            # and a plain lookup returns nothing -- the guard then skipped exactly the
            # cases it exists for. _station_on_route answers where the vessel's group
            # meets the aorta, which is the station that orders them.
            sc = _station_on_route(jg, stations, attach_idx(car))
            ss = _station_on_route(jg, stations, attach_idx(sub))
            if sc is None or ss is None or sc <= ss:
                continue
            names[car.label], names[sub.label] = (f"{side} subclavian",
                                                  f"{side} common carotid")
            notes.append(
                f"{side.lower()} carotid and subclavian swapped: the one named "
                f"carotid left the aorta DISTAL to the one named subclavian, which "
                f"no arch does. The distal vessel is the subclavian")

        if len([t for t in branches if t.label in names]) < len(branches):
            notes.append("some arch vessels were left as Branch N, their signature "
                         "was not decisive")
        return _apply_manual(tree, names, vessels, notes)

    # ---- abdominal: one inlet, visceral branches, a terminal iliac bifurcation ----
    if len(trunk) == 1:
        inlet = trunk[0]
        names[inlet.label] = "Supraceliac aorta (inlet)"

        # Midline and scale from the aorta itself: inlet down to the iliac bifurcation.
        pairs = [(n, d) for n, d in down.items() if len(d["terms"]) == 2]
        iliacs = []
        for i in range(len(pairs)):
            for k in range(i + 1, len(pairs)):
                (n1, d1), (n2, d2) = pairs[i], pairs[k]
                if d1["sub"] & d2["sub"]:
                    continue                      # nested, not siblings
                iliacs.append((n1, d1, n2, d2))
        # The aortic bifurcation: the two iliacs must SHARE the junction they leave
        # from. Without that a renal that happens to divide in two reads as a "common
        # iliac", and the real contralateral iliac is dropped.
        iliacs = [p for p in iliacs
                  if np.linalg.norm(p[1]["upos"] - p[3]["upos"]) < 1e-6]
        # The aortic bifurcation is the sibling pair FURTHEST DOWN THE AORTA from the
        # inlet -- a distance along the vessel, not "lowest z", so a z-flipped export
        # still finds the iliacs rather than the visceral branches.
        def _far(pair):
            st = _aorta_stations(tree, jg, inlet.segment.name, pair[0])
            return st.get(pair[1]["upidx"], -1.0)
        bifurc = max(iliacs, key=_far) if iliacs else None
        far_seg = bifurc[0] if bifurc else inlet.segment.name
        aorta = _aorta_path(tree, jg, inlet.segment.name, far_seg)
        axis = np.vstack([tree.segments[n].pts for n in aorta])
        mid_x = float(np.median(axis[:, 0]))
        R = max(float(getattr(inlet.segment, "_radius_mm", None) or 0.0), 1.0)
        stations = _aorta_stations(tree, jg, inlet.segment.name, far_seg)
        # "Up" the body = back along the aorta toward the inlet.
        far_pt = axis[int(np.argmax(np.linalg.norm(
            axis - np.asarray(inlet.tip, float), axis=1)))]
        up = _unit_np(np.asarray(inlet.tip, float) - far_pt)

        # Celiac trunk: an ANTERIOR vessel that fans out into three openings. Named
        # FIRST, because its splenic/hepatic pair is the handedness witness -- the
        # spleen is on the left and the liver on the right, whatever the file says.
        left_sign = None
        # Measured from where it LEAVES THE AORTA, heading downstream -- read the other
        # way round a visceral trunk looks like it points backwards. Named only when
        # exactly one anterior trifurcation exists.
        celiacs = [n for n, d in down.items()
                   if len(d["terms"]) == 3
                   and _leaves_along(tree, tree.segments[n], d["upos"])[1] <= -0.30]
        for n in celiacs[:1] if len(celiacs) == 1 else []:
            d = down[n]
            vessels[n] = "Celiac trunk"
            kids = d["terms"]
            # The gastric climbs back toward the stomach: furthest along "up the body",
            # measured from the celiac's own origin rather than in raw z.
            gastric = max(kids, key=lambda t: float(
                np.dot(np.asarray(t.tip, float) - d["upos"], up)))
            rest = [t for t in kids if t is not gastric]
            if len(rest) == 2:
                names[gastric.label] = "Left gastric artery"
                spleen = max(rest, key=lambda t: t.segment.length)  # long and tortuous
                liver = [t for t in rest if t is not spleen][0]
                names[spleen.label] = "Splenic artery"
                names[liver.label] = "Common hepatic artery"
                dx = float(spleen.tip[0]) - float(liver.tip[0])
                if abs(dx) > SIDE_TOL_R * R:
                    left_sign = 1.0 if dx > 0 else -1.0   # spleen side IS the left
            break

        if left_sign is None:
            left_sign = 1.0
            notes.append("no splenic/hepatic pair found to confirm left from right, "
                         + SIDES_ASSUMED_NOTE)

        def side_of(x):
            dx = (float(x) - mid_x) * left_sign
            return "" if abs(dx) < SIDE_TOL_R * R else ("Left" if dx > 0 else "Right")

        if bifurc:
            n1, d1, n2, d2 = bifurc
            # Judge the two limbs AGAINST EACH OTHER, not against mid_x.
            #
            # side_of compares a limb to the median x of the aortic path -- and that
            # path is built out to `far_seg`, which IS one of these two limbs, so the
            # midline is dragged toward whichever sorted first (measured on JK5_1:
            # 18.05 with the limb in, 20.54 without, a 2.5 mm bias against a 4.15 mm
            # deadband). It is also a GLOBAL statistic being used as a LOCAL midline
            # at the far end of the aorta.
            #
            # The two limbs of a bifurcation are on opposite sides of each other by
            # construction, so their separation is the measurement that was wanted and
            # it needs no midline at all. This is the same relative comparison the
            # celiac already uses to fix handedness from the splenic/hepatic pair.
            #
            # It closes both ways the old test failed: one limb inside the deadband
            # was skipped entirely and left as "Branch N" with no note, and a midline
            # outside BOTH limbs named them the same side, after which _unique
            # withdrew all six iliac names (seen on JK36_2).
            dx = (float(d1["jpos"][0]) - float(d2["jpos"][0])) * left_sign
            if abs(dx) < SIDE_TOL_R * R:
                notes.append("the two limbs of the aortic bifurcation are not far "
                             "enough apart to tell left from right, so they were "
                             "left unnamed -- name them by hand")
                pairs = ()
            else:
                pairs = tuple(zip(((n1, d1), (n2, d2)),
                                  ("Left", "Right") if dx > 0 else ("Right", "Left")))
            for (n, d), side in pairs:
                vessels[n] = f"{side} common iliac"
                # internal (hypogastric) runs posteriorly into the pelvis; external
                # carries on inferolaterally to become the femoral.
                kids = d["terms"]
                if len(kids) == 2:
                    post = max(kids, key=lambda t: float(t.tip[1]))
                    ext = min(kids, key=lambda t: float(t.tip[1]))
                    if post is not ext:
                        names[post.label] = f"{side} internal iliac"
                        names[ext.label] = f"{side} external iliac"

        # SMA: a long anterior branch of its own, below the celiac. It is named only
        # if it is the ONLY thing that fits -- an abdominal aorta also carries an IMA,
        # and lumbar and median sacral branches, any of which can look similar enough.
        # Picking the "best" match among several is how an IMA ends up called the SMA.
        # "Below the celiac" = further along the aorta than it, not lower in z.
        celiac_st = min((stations.get(d["upidx"], float("inf"))
                         for n, d in down.items() if vessels.get(n) == "Celiac trunk"),
                        default=None)
        sma = []
        for t in branches:
            if t.label in names:
                continue
            j = attach(t)
            st = stations.get(attach_idx(t))
            if j is None or (celiac_st is not None
                             and (st is None or st <= celiac_st)):
                continue
            if t.segment.length < 2.0 * R:        # scale-relative, not a fixed 25 mm
                continue
            if _leaves_along(tree, t.segment, j)[1] < -0.45:
                sma.append(t)
        if len(sma) == 1:
            names[sma[0].label] = "Superior mesenteric artery (SMA)"
        elif len(sma) > 1:
            notes.append(f"{len(sma)} anterior branches below the celiac could each be "
                         "the SMA (an IMA or a lumbar looks similar), none named")

        # Renals: a lateral pair leaving at essentially the same level, on opposite
        # sides. Several pairs can fit -- an accessory renal is ~30% of people, and a
        # lumbar pair looks the same -- so when they do, the pair that is MOST LEVEL
        # is taken and the choice is reported.
        #
        # That tie-break is measured, not assumed. Over the lab's 47 abdominal scans,
        # 27 have exactly one pair and everyone accepts it as the answer. Widening the
        # level window on those manufactures distractors while the answer stays known,
        # which is the only held-out test available here. Across windows of 3, 5 and 8
        # aortic radii the most-level pair finds the accepted answer 8/8, 9/9 and
        # 11/11; the fattest pair manages 5, the most proximal 3, "the first pair
        # below the visceral group" 3. Scoring against the UNWIDENED candidates proves
        # nothing at all -- with one candidate every rule picks it, and six mutually
        # contradictory rules all scored 27/27 that way.
        #
        # The renals arise as a true pair at one aortic level, which is why levelness
        # separates them: an accessory paired with a real one, or a lumbar pair, is
        # not as level. See tools/ABDOMINAL.md.
        rest = [t for t in branches if t.label not in names]
        pairs_r = []
        for i in range(len(rest)):
            for k in range(i + 1, len(rest)):
                a, b = rest[i], rest[k]
                ja, jb = attach(a), attach(b)
                if ja is None or jb is None or abs(ja[2] - jb[2]) > 1.5 * R:
                    continue
                sa, sb = side_of(a.tip[0]), side_of(b.tip[0])
                if not sa or not sb or sa == sb:
                    continue
                pairs_r.append((abs(float(ja[2] - jb[2])), a, b, sa, sb))
        if pairs_r:
            pairs_r.sort(key=lambda p: p[0])
            _dz, a, b, sa, sb = pairs_r[0]
            names[a.label], names[b.label] = f"{sa} renal artery", f"{sb} renal artery"
            if len(pairs_r) > 1:
                # Named, but as a choice rather than a reading. The interface shows
                # notes, and any opening can be re-declared by hand.
                notes.append(
                    f"{len(pairs_r)} lateral branch pairs could be the renals (an "
                    f"accessory renal or a lumbar pair looks the same); the most "
                    f"level pair was taken, {_dz:.1f} mm apart against "
                    f"{pairs_r[1][0]:.1f} mm for the next -- worth confirming, and "
                    f"any opening can be named by hand")

        if len([t for t in branches if t.label in names]) < len(branches):
            notes.append("some branches were left as Branch N, their signature was "
                         "not decisive (accessory or unsegmented vessels are common)")
    return _apply_manual(tree, *_unique(names, vessels, notes, claimed))


# A correct planar clip lands its rim ON the plane (KY160 measures 0.000 mm), so
# anything past a fraction of a millimetre means the opening was torn, not sliced.
RIM_FLATNESS_MM = 0.5


def _loops_by_terminal(mesh, tree, offsets=None, default_offset_mm=6.0, terminals=None):
    """Terminal label -> the vertex loop of the opening it made.

    *terminals* overrides which cuts to match (the GUI passes its active set, which
    can include a pre-bifurcation cut in place of the branches it removed).
    """
    offsets = dict(offsets or {})
    terms = tree.terminals if terminals is None else terminals
    loops = boundary_loops(mesh)
    if not loops:
        return {}
    cents = np.array([mesh.vertices[l].mean(axis=0) for l in loops])
    out = {}
    for t in terms:
        point, _ = t.cut_plane(float(offsets.get(t.label, default_offset_mm)))
        out[t.label] = loops[int(np.argmin(np.linalg.norm(cents - point, axis=1)))]
    return out


def match_openings_to_terminals(mesh, tree, offsets=None, default_offset_mm=6.0,
                                terminals=None):
    """Map each terminal label -> the opening it produced.

    boundary_loops() returns loops in whatever order the connectivity search
    finds them, which need NOT line up with tree.terminals -- so match each
    terminal to the opening whose centroid is nearest its cut point rather than
    zipping the two lists and hoping. *terminals* overrides the cut list.
    """
    offsets = dict(offsets or {})
    terms = tree.terminals if terminals is None else terminals
    stats = opening_stats(mesh)
    out = {}
    if not stats:
        return out

    cents = np.array([s["centroid"] for s in stats])
    for t in terms:
        point, _ = t.cut_plane(float(offsets.get(t.label, default_offset_mm)))
        j = int(np.argmin(np.linalg.norm(cents - point, axis=1)))
        out[t.label] = stats[j]
    return out


# How far a cut plane may lean away from perpendicular-to-the-vessel before the
# opening stops being usable. A centerline-normal cut is 0 deg by construction; a
# Z-plane's obliquity depends on which way the vessel happens to run there.
#
# Measured, not guessed. Sweeping both Z-cut terminals of KY160 across their whole
# range separates perfectly:
#     FLAT openings : obliquity  3 - 50 deg
#     TORN openings : obliquity 53 - 71 deg
# Past ~50 deg the plane is grazing ALONG the vessel instead of across it: the rim
# stops being planar (up to 8.9 mm off), the opening stretches (roundness 0.34) and
# its apparent diameter inflates (59 mm on a 33 mm vessel).
MAX_OBLIQUITY_DEG = 50.0

# Head-room kept when CHOOSING a horizontal plane for an aortic end, as opposed to
# refusing one. The angle is measured at the default offset, and a user who then moves
# the cut a few millimetres changes it -- so a mode picked at exactly the limit is one
# slider nudge from refusing. Preflight still allows anything up to the limit.
OBLIQUITY_MARGIN_DEG = 8.0


def obliquity_deg(terminal, offset_mm):
    """Angle between the cut plane's normal and the vessel it is cutting.

    0 deg = perpendicular (a clean ring). 90 deg = the plane lies along the vessel.
    """
    s = terminal.s_cut(offset_mm)
    _, tan = terminal.segment.at_arclength(s)
    _, n = terminal.cut_plane(offset_mm)
    return float(np.degrees(np.arccos(min(abs(float(np.dot(tan, n))), 1.0))))

# A cut AT the tip removes nothing: the terminal sits on the surface, so there is
# no vessel past it and no opening to make. 1 mm clears it (and is well under the
# ~2 mm mesh edge, so the plane still meets real faces).
MIN_OFFSET_MM = 1.0

# Sliver snapping is a MESH-DENSITY question, not an absolute-distance one. A
# sliver appears when the cut plane grazes a triangle close to one of its corners,
# and "close" only means anything relative to that triangle's size. A fixed 0.1 mm
# would be half an edge on a 0.2 mm mesh (far too aggressive) and invisible on a
# 5 mm one (far too timid). So express it as a fraction of the local edge length:
# 0.05 x the median edge reproduces the 0.1 mm that was measured to work on the
# sample mesh (median edge 2.08 mm), and carries to any other density.
SNAP_FRAC_OF_EDGE = 0.05


def snap_tolerance(mesh, frac=SNAP_FRAC_OF_EDGE):
    """Sliver-snap tolerance scaled to this mesh's own tessellation."""
    return float(frac) * float(np.median(mesh.edges_unique_length))


def z_cut_hint(tree, offsets=None, default_offset_mm=6.0):
    """Kept for callers; the obliquity check in preflight() now does the work."""
    return ""


def preflight(tree, offsets=None, default_offset_mm=6.0, terminals=None):
    """Check the requested cuts make sense BEFORE any mesh work. Returns a problem
    string, or None.

    Two terminals can share one segment -- an unbranched aorta is exactly that,
    and the two planes then trim the SAME tube from opposite ends. If the
    start-side cut lands beyond the end-side cut they have crossed and NOTHING is
    left between them; the cutter then dies with an opaque numpy error.

    The comparison has to respect which end each terminal belongs to. Merely
    sorting the two cuts by position and measuring the distance between them
    looks fine in both cases -- a crossed pair reports the same healthy "27 mm
    gap" as a valid one -- because sorting throws away the very thing that
    distinguishes them.

    *terminals* overrides the cut list (the GUI's active set, which may hold a
    pre-bifurcation cut in place of the branches it removes).
    """
    offsets = dict(offsets or {})
    terms = tree.terminals if terminals is None else terminals

    def off(t):
        return float(offsets.get(t.label, default_offset_mm))

    # Does this centerline belong to this mesh? segment_radius is the median signed
    # distance from a segment's centerline into the surface, so a value at or below
    # zero says that stretch of centerline runs OUTSIDE the vessel it is supposed to
    # be inside. Everything downstream is then meaningless: the cut band is sized from
    # that radius, ownership is normalised by it, and preflight's own "two radii apart"
    # test uses it.
    #
    # classify_terminals only refuses when the LARGEST terminal is outside, which
    # catches a wholesale mismatch and nothing else. Measured across the two lab
    # cohorts, the common case is partial: the trunk pairs correctly and the branch
    # centerlines do not, because the mesh supplied was the trunk-only export rather
    # than the branched one. Seven of nine segments outside the surface, and the tool
    # cut it without a word.
    outside = [t for t in terms
               if float(getattr(t.segment, "_radius_mm", None) or 0.0) <= 0.0]
    if outside:
        names = ", ".join(f"“{t.name}”" for t in outside[:4])
        more = f" and {len(outside) - 4} more" if len(outside) > 4 else ""
        return (f"{names}{more}: the centerline for {'these' if len(outside) > 1 else 'this'} "
                f"opening{'s' if len(outside) > 1 else ''} runs OUTSIDE the mesh, so "
                f"there is no vessel here to cut. The usual cause is a mesh and a "
                f"centerline from different exports, most often a trunk-only surface "
                f"paired with a branched centerline. Check they are the same case, in "
                f"the same frame and the same units.")

    for t in terms:
        if off(t) < MIN_OFFSET_MM:
            return (f"“{t.name}”: a cut of {off(t):.1f} mm sits on the vessel tip "
                    f"itself, so there is nothing past it to remove. "
                    f"Use at least {MIN_OFFSET_MM:.0f} mm.")

    # The plane must cut ACROSS the vessel, not along it. This is what tears an
    # opening, and it has nothing to do with how deep the cut is: the aortic root
    # stays clean to 85% of its length, while the descending aorta tears at 25%
    # simply because it happens to run horizontal there.
    for t in terms:
        ob = obliquity_deg(t, off(t))
        if ob > MAX_OBLIQUITY_DEG:
            fix = ("Switch it to Centerline-normal, which is perpendicular by "
                   "construction." if t.mode == "z" else
                   "Reduce its tilt.")
            return (f"“{t.name}”: the cut plane meets the vessel at {ob:.0f}°, it "
                    f"is slicing ALONG the vessel rather than across it, which "
                    f"tears the opening instead of cutting it. {fix}")

    by_seg = {}
    for t in terms:
        by_seg.setdefault(t.segment.name, []).append(t)

    for name, ts in by_seg.items():
        starts = [t for t in ts if not t.at_end]
        ends = [t for t in ts if t.at_end]
        if not starts or not ends:
            continue
        for t_s in starts:
            for t_e in ends:
                # Kept vessel is the stretch between them: s_cut(start) .. s_cut(end).
                lo, hi = t_s.s_cut(off(t_s)), t_e.s_cut(off(t_e))
                r = max(float(getattr(t_s.segment, "_radius_mm", 5.0)), 1.0)
                # Two cuts facing each other across an ISOLATED REGION are meant to be
                # close -- an Ishimaru zone 3 is 20 mm on a 36 mm aorta by definition,
                # and demanding 2 radii between them would reject it as "almost met".
                # They still must not cross.
                gap = 0.0 if (getattr(t_s, "is_interior", False)
                              and getattr(t_e, "is_interior", False)) else 2.0 * r
                if hi - lo < gap:
                    if hi <= lo:
                        return (f"“{t_s.name}” and “{t_e.name}” have cut past each "
                                f"other, no vessel would be left between them. "
                                f"Pull one of them back.")
                    return (f"“{t_s.name}” and “{t_e.name}” have almost met, only "
                            f"{hi - lo:.0f} mm of vessel would be left between "
                            f"them. Pull one of them back.")
    return None


def _deepest_cut_hint(tree, offsets=None, default_offset_mm=6.0, terminals=None):
    """Name the cut that has travelled furthest into its vessel.

    "4 openings, expected 5" tells you something broke but not what to do about
    it. With the sliders opened up to nearly the whole vessel, a cut that has run
    too close to its junction is far and away the likeliest culprit, and it is the
    one to pull back.
    """
    offsets = dict(offsets or {})
    worst, frac = None, 0.0
    for t in (tree.terminals if terminals is None else terminals):
        f = float(offsets.get(t.label, default_offset_mm)) / max(t.segment.length, 1e-9)
        if f > frac:
            worst, frac = t, f
    if worst is None or frac < 0.60:
        return ""
    return (f" ,  “{worst.name}” is trimmed {100*frac:.0f}% of the way down its "
            f"vessel, close to the junction it springs from. Pull it back.")


def validate_cut(mesh, tree, offsets=None, default_offset_mm=6.0, terminals=None):
    """Check an open mesh is actually correct. Returns a problem string, or None.

    Topology alone is NOT enough, and that is the whole lesson of this module: a
    cut that leaks into the parent vessel still removes one connected disk, so it
    still leaves one boundary loop and the Euler number stays exactly right. The
    leaked openings were only caught by geometry -- their centroids sat 39-55 mm
    from the plane that was supposed to have made them.

    So we check both: the topology AND that every opening sits on its own cut
    plane, at a radius consistent with its vessel. *terminals* overrides the cut
    list (a pre-bifurcation cut counts as ONE opening in place of its branches).
    """
    offsets = dict(offsets or {})
    terms = tree.terminals if terminals is None else terminals
    want = len(terms)

    pre = preflight(tree, offsets, default_offset_mm, terminals=terms)
    if pre is not None:
        return pre

    hint = z_cut_hint(tree, offsets, default_offset_mm) or _deepest_cut_hint(
        tree, offsets, default_offset_mm, terminals=terms)

    loops = boundary_loops(mesh)
    if len(loops) != want:
        return f"{len(loops)} openings, expected {want}{hint}"

    euler = int(mesh.euler_number)
    if euler != 2 - want:
        return f"euler {euler}, expected {2 - want}{hint}"

    matched = match_openings_to_terminals(mesh, tree, offsets, default_offset_mm,
                                          terminals=terms)
    loops_by_label = _loops_by_terminal(mesh, tree, offsets, default_offset_mm,
                                        terminals=terms)

    for t in terms:
        s = matched.get(t.label)
        if s is None:
            return f"{t.name}: no opening{hint}"
        off = float(offsets.get(t.label, default_offset_mm))
        point, normal = t.cut_plane(off)
        r = max(float(getattr(t.segment, "_radius_mm", 5.0)), 1.0)

        drift = float(np.linalg.norm(np.asarray(s["centroid"]) - point))
        if drift > max(2.0 * r, 12.0):
            return (f"{t.name}: opening is {drift:.0f} mm off its cut plane "
                    f"(cut escaped this vessel){hint}")

        # The opening must actually be FLAT. A clip puts its rim exactly on the
        # plane (KY160 measures 0.000 mm), so any real deviation means the cut
        # tore rather than sliced -- and a ragged, non-planar inlet is useless to
        # a solver. Topology cannot see this: the hole still has one loop and the
        # euler number is still right. It went unnoticed until a tortuous vessel
        # produced an 11 mm-out-of-plane rim that sailed through every other check.
        loop = loops_by_label.get(t.label)
        if loop is not None and len(loop) >= 3:
            dev = float(np.abs((mesh.vertices[loop] - point) @ normal).max())
            if dev > RIM_FLATNESS_MM:
                # Say what actually helps. Measured over 400 cuts on the two sample
                # aortas, EVERY torn opening sat 69-90% of the way down its vessel --
                # within a few mm of the junction where it merges into its parent,
                # where the vessel has no clean cross-section left to cut. Neither
                # widening nor shrinking the plane rescues those: the tear scales
                # with the plane's reach, and shrinking it far enough to matter
                # leaves nothing to clip. Pulling the cut back is what works.
                # The junction is at the END the terminal is NOT on, so the distance
                # to it is measured from the cut station to that far end.
                s_c = t.s_cut(off)
                gap = s_c if t.at_end else (t.segment.length - s_c)
                deep = off > 0.6 * t.segment.length
                cure = (f"It sits {gap:.0f} mm from the junction this vessel springs "
                        f"from, where there is no longer a clean cross-section to "
                        f"cut. Pull it back." if deep else
                        "Move this cut, or shrink its plane size.")
                return (f"{t.name}: the opening is not flat, its rim strays "
                        f"{dev:.1f} mm off the cut plane (torn, not sliced). "
                        f"{cure}{hint}")
    return None


def keep_largest_component(mesh):
    """Largest surface-connected piece.

    Done on the face-adjacency graph rather than via mesh.split(), which rebuilds
    a full Trimesh per component -- far too slow to run inside every cut (it alone
    was most of a ~1.1 s re-cut).
    """
    fa = np.asarray(mesh.face_adjacency)
    nf = len(mesh.faces)
    if nf == 0 or len(fa) == 0:
        return mesh

    g = coo_matrix((np.ones(len(fa)), (fa[:, 0], fa[:, 1])), shape=(nf, nf))
    n, labels = connected_components(g, directed=False)
    if n <= 1:
        return mesh

    keep = labels == int(np.argmax(np.bincount(labels, minlength=n)))
    out = trimesh.Trimesh(mesh.vertices, mesh.faces[keep], process=False)
    out.remove_unreferenced_vertices()
    return out


def open_terminals(mesh, tree, offsets=None, default_offset_mm=6.0, snap_mm=0.0,
                   verbose=True):
    """Open every terminal of *tree*. *offsets* maps terminal label -> mm to cut
    back from the tip; anything unlisted uses *default_offset_mm*."""
    offsets = dict(offsets or {})

    pre = preflight(tree, offsets, default_offset_mm)
    if pre is not None:
        raise ValueError(pre)

    out = mesh.copy()

    for t in tree.terminals:
        off = float(offsets.get(t.label, default_offset_mm))
        before = len(out.faces)
        out = cut_terminal(out, tree, t, off, snap_mm=snap_mm)
        if verbose:
            print(
                f"  cut {t.name:22s} [{t.mode:6s}] offset={off:6.1f} mm  "
                f"faces {before} -> {len(out.faces)}"
            )

    out = keep_largest_component(out)
    if verbose:
        loops = boundary_loops(out)
        euler = int(out.euler_number)
        problem = validate_cut(out, tree, offsets, default_offset_mm)
        print(
            f"  result: {len(out.faces)} faces, {len(loops)} openings, "
            f"euler={euler} "
            f"{'OK' if problem is None else '** ' + problem + ' **'}"
        )
        matched = match_openings_to_terminals(out, tree, offsets, default_offset_mm)
        for t in tree.terminals:
            s = matched.get(t.label)
            if s:
                print(
                    f"    {t.name:22s} diameter~{s['diameter']:5.1f} mm  "
                    f"@ {np.round(s['centroid'], 1)}"
                )
    return out


def cut_open(csv_path, mat_path, offsets=None, default_offset_mm=6.0,
             target_faces=None, z_cut_ends=True, verbose=True):
    """Headless end-to-end: centerline CSV + closed mesh -> open mesh."""
    from branch_tree import load_tree

    tree = load_tree(csv_path)
    mesh = load_mat_mesh(mat_path)
    if verbose:
        print(tree.summary())
        print(f"\nmesh: {len(mesh.faces)} faces, watertight={mesh.is_watertight}\n")

    classify_terminals(tree, mesh, z_cut_ends=z_cut_ends, verbose=verbose)
    if verbose:
        print()

    out = open_terminals(mesh, tree, offsets=offsets,
                         default_offset_mm=default_offset_mm, verbose=verbose)
    if target_faces:
        out = remesh_uniform(out, target_faces=target_faces)
        if verbose:
            print(f"  remeshed -> {len(out.faces)} faces, "
                  f"euler={int(out.euler_number)}")
    return out, tree


def rim_sliver_count(mesh, q_thresh=0.01):
    """How many near-degenerate triangles the clip left on the rim.

    Quality = 4*sqrt(3)*area / sum(edge^2): 1 for equilateral, ~0 for a sliver.
    """
    tri = mesh.vertices[mesh.faces]
    e = [np.linalg.norm(tri[:, (i + 1) % 3] - tri[:, i], axis=1) for i in range(3)]
    denom = np.maximum(e[0] ** 2 + e[1] ** 2 + e[2] ** 2, 1e-30)
    q = 4.0 * np.sqrt(3.0) * mesh.area_faces / denom
    return int((q < q_thresh).sum())


def min_triangle_quality(mesh):
    """Worst triangle quality in the mesh (1 = equilateral, ~0 = degenerate)."""
    tri = mesh.vertices[mesh.faces]
    e = [np.linalg.norm(tri[:, (i + 1) % 3] - tri[:, i], axis=1) for i in range(3)]
    denom = np.maximum(e[0] ** 2 + e[1] ** 2 + e[2] ** 2, 1e-30)
    return float((4.0 * np.sqrt(3.0) * mesh.area_faces / denom).min())


def remesh_uniform(mesh, target_faces=50_000):
    """Uniform ACVD remesh.

    NOT part of the default pipeline and never called by the GUI -- the cut mesh
    is the original geometry with its caps removed, and that is normally what you
    want. Only reachable by passing target_faces to cut_open().

    This used to carry its own copy of the clustering, which meant it also carried
    its own copy of the bug where a small mesh silently returns fewer faces than
    asked for. One implementation now, imported late so that pyacvd stays optional
    for everyone who never remeshes.
    """
    from cutclean import acvd_uniform_remesh

    return acvd_uniform_remesh(mesh, target_faces=target_faces)


if __name__ == "__main__":
    import sys
    from pathlib import Path

    root = Path(__file__).parent
    csv = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "KY160_1_CentrelinePoints.csv"
    mat = Path(sys.argv[2]) if len(sys.argv) > 2 else \
        root / "Geometries & Meshes" / "mesh" / "KY160_1.mat"

    out_mesh, _ = cut_open(csv, mat)
    out_dir = root / "output"
    out_dir.mkdir(exist_ok=True)
    dest = out_dir / f"{mat.stem}_open.stl"
    out_mesh.export(dest)
    print(f"\nwrote {dest}")
