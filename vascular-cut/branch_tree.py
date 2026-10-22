# Branched-centerline tree

# Parses a Simpleware per-data-point spline export (one "Spline name" per
# centerline segment) and assembles the segments into a tree by matching
# coincident endpoints. Exposes, for each terminal, a cut plane parameterised by
# how far back along the centerline the cut should sit.

import io
import json
import re

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import cKDTree

POS_COLS = ["Position x (global)", "Position y", "Position z"]

# Two centerline dialects are in circulation, and they differ only in bookkeeping:
#
#   Simpleware spline export : "Spline name" + "First derivative x/y/z", one row
#                              block per branch, a provenance comment above the
#                              header.  (KY160_1_CentrelinePoints.csv)
#   cl.py's save_centerline_csv : "Tangent x/y/z", header on line 0, no spline
#                              names, one unbranched path per file. (KK2_5.csv)
#
# Both carry the same physics (a position and a direction per point), so the
# loader sniffs which one it is rather than demanding a particular flavour.
TANGENT_COLS = [
    ["First derivative x", "First derivative y", "First derivative z"],
    ["Tangent x", "Tangent y", "Tangent z"],
]
D1_COLS = TANGENT_COLS[0]        # kept for callers that reference it by name

# A file may hold several polylines back to back with no marker between them.
# A genuine break shows up as a step far larger than the sampling interval.
BREAK_FACTOR = 8.0

# Endpoints exported from the same spline network are bit-identical, but keep a
# tolerance so hand-edited or re-exported files still stitch.
JOIN_TOL_MM = 0.5
# Shortest stub a T-junction repair may leave behind. Splitting a host nearer than
# this to one of its own ends makes a sliver whose length is under a normal cut
# depth, so the cut would sit past the end of its own vessel.
MIN_SPLIT_MM = 8.0


def _unit(v, eps=1e-12):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, eps)


def _rotate(v, axis, angle_rad):
    """Rodrigues rotation of *v* about a unit *axis*."""
    v = np.asarray(v, float)
    k = _unit(axis)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return v * c + np.cross(k, v) * s + k * float(np.dot(k, v)) * (1.0 - c)


# Above this much tangent turning, a centerline is a staircase rather than a curve
# and its own tangent column cannot be trusted for a cut normal.
JAGGED_DEG_PER_MM = 15.0

# Chord length used to re-derive tangents on a jagged centerline.
CHORD_SPAN_MM = 8.0

# How far from each end the endpoint-pin correction is blended back to zero.
PIN_WINDOW_MM = 10.0


def _chord_tangents(pts, s, span_mm):
    """Tangents from a chord spanning +/- span/2 mm, instead of point-to-point.

    A voxel-geodesic path (cl.py's exporter) runs dead straight along voxel edges
    and then turns a corner: its per-point tangent is either 0 or enormous. Taking
    the direction between two points a real distance apart averages the staircase
    out and recovers the direction the vessel is actually going.
    """
    pts = np.asarray(pts, float)
    half = 0.5 * float(span_mm)
    lo = np.searchsorted(s, s - half, side="left")
    hi = np.searchsorted(s, s + half, side="right") - 1
    lo = np.clip(lo, 0, len(pts) - 1)
    hi = np.clip(hi, 0, len(pts) - 1)
    same = hi <= lo                                   # degenerate at the very ends
    hi[same] = np.minimum(lo[same] + 1, len(pts) - 1)
    lo[same] = np.maximum(hi[same] - 1, 0)
    return _unit(pts[hi] - pts[lo])


def tangent_turning(pts, tangents, s):
    """Degrees the tangent turns per mm, how jagged a centerline is."""
    if len(pts) < 3:
        return 0.0
    dot = np.clip(np.einsum("ij,ij->i", tangents[:-1], tangents[1:]), -1.0, 1.0)
    turn = np.degrees(np.arccos(dot))
    step = np.maximum(np.diff(s), 1e-9)
    return float(np.percentile(turn / step, 95))


def smooth_polyline(pts, rms_mm, n_out=None):
    """Fit a smoothing spline through *pts*, allowing ~rms_mm of deviation.

    This is what de-staircases a voxel-geodesic path. Keep rms_mm small: on the
    sample, 1 mm removed the quantisation AND left the line deeper inside the
    lumen (min wall distance rose 3.67 -> 4.74 mm), while 4 mm started cutting
    corners and hugging walls (min wall distance collapsed to 1.54 mm).
    """
    from scipy.interpolate import splprep, splev

    pts = np.asarray(pts, float)
    n = len(pts)
    if n < 4 or rms_mm <= 0:
        return pts

    tck, _ = splprep([pts[:, 0], pts[:, 1], pts[:, 2]], s=n * float(rms_mm) ** 2, k=3)
    u = np.linspace(0.0, 1.0, n_out or n)
    q = np.stack(splev(u, tck), axis=1)

    # resample to uniform arc length so downstream arc-length maths stays honest
    qs = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))))
    if qs[-1] <= 0:
        return pts
    uu = np.linspace(0.0, qs[-1], n_out or n)
    q = np.stack([np.interp(uu, qs, q[:, k]) for k in range(3)], axis=1)

    # PIN THE ENDS. A branched centerline is a set of segments that meet at shared
    # endpoints, and build_tree reconstitutes the tree by matching those within
    # JOIN_TOL_MM. Smoothing moves an endpoint by 0.6-2.4 mm on the samples -- past
    # that tolerance -- so every junction dissolves and the tree falls apart into
    # loose segments (JK5 went from 11 openings to 28).
    #
    # The correction decays over a short window at each end rather than across the
    # whole segment: spread over the full length it drags the MIDDLE of a thin branch
    # sideways too, eating into its clearance from the vessel wall (measured -1.2 mm
    # on JK5). Confined to a few mm it fixes the endpoints and leaves the interior --
    # the part smoothing exists to straighten -- where the fit put it.
    s = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))))
    win = max(min(PIN_WINDOW_MM, 0.5 * s[-1]), 1e-6)
    w0 = np.clip(1.0 - s / win, 0.0, 1.0)[:, None]
    w1 = np.clip(1.0 - (s[-1] - s) / win, 0.0, 1.0)[:, None]
    return q + w0 * (pts[0] - q[0]) + w1 * (pts[-1] - q[-1])


class Segment:
    """One centerline curve: positions, unit tangents, cumulative arc length.

    *smooth_mm* > 0 fits a smoothing spline through the points, which is how a
    voxel-geodesic staircase (cl.py's exporter) is turned back into a curve. It is
    OFF by default and must be asked for: a Simpleware spline export is already a
    fitted curve (0.6-2.8 deg/mm of turning) and smoothing it would only degrade
    it. Nothing is ever smoothed behind the user's back.
    """

    def __init__(self, name, pts, tangents, smooth_mm=0.0):
        self.name = str(name)
        raw_pts = np.asarray(pts, dtype=float)
        raw_tans = _unit(tangents)

        seg = np.linalg.norm(np.diff(raw_pts, axis=0), axis=1)
        raw_s = np.concatenate(([0.0], np.cumsum(seg)))

        # Measured on the file AS SUPPLIED, this is what decides whether the tool
        # suggests smoothing, so it must never see a smoothed curve.
        self.raw_turning = tangent_turning(raw_pts, raw_tans, raw_s)
        self.jagged = self.raw_turning > JAGGED_DEG_PER_MM
        self.smooth_mm = float(smooth_mm)

        if self.smooth_mm > 0:
            self.pts = smooth_polyline(raw_pts, self.smooth_mm)
            d = np.linalg.norm(np.diff(self.pts, axis=0), axis=1)
            self.s = np.concatenate(([0.0], np.cumsum(d)))
            self.tangents = _chord_tangents(self.pts, self.s, CHORD_SPAN_MM)
        else:
            self.pts = raw_pts
            self.s = raw_s
            self.tangents = raw_tans

        self.length = float(self.s[-1])
        self.turning = tangent_turning(self.pts, self.tangents, self.s)

        # Span over which cut_plane averages the tangent to orient its plane. A
        # rough centerline needs a longer one, smoothed or not.
        self.tangent_span_mm = 10.0 if self.turning > 6.0 else 4.0

    def __len__(self):
        return len(self.pts)

    def at_arclength(self, s):
        """Position and unit tangent at arc length *s* from the segment start.

        The tangent always points in the direction of increasing s.
        """
        s = float(np.clip(s, 0.0, self.length))
        pos = np.array([np.interp(s, self.s, self.pts[:, k]) for k in range(3)])
        tan = np.array([np.interp(s, self.s, self.tangents[:, k]) for k in range(3)])
        return pos, _unit(tan)

    def endpoint(self, at_end):
        return self.pts[-1] if at_end else self.pts[0]


class Terminal:
    """A free end of the tree, one of the openings we want to cut.

    *at_end* says which end of the parent segment is free, which fixes the sign
    of the outward direction: tangents run start -> end, so a terminal at the
    end points outward along +tangent and one at the start along -tangent.

    *mode* picks the cut plane's orientation:
      "normal", perpendicular to the centerline (right for branch tips, whose
                 caps sit within a few degrees of centerline-perpendicular).
      "z"     , horizontal, normal = +/-z (right for the two aorta ends, which
                 Simpleware cut on a roughly horizontal plane; measured 7-10 deg
                 from z but 17-25 deg from centerline-perpendicular).
    The +/- sign in "z" mode is taken from the centerline's own outward tangent,
    so it never has to be specified by hand.
    """

    def __init__(self, segment, at_end, label=None, mode="normal", name=None):
        self.segment = segment
        self.at_end = bool(at_end)
        self.label = label or f"{segment.name}.{'end' if at_end else 'start'}"
        self.mode = mode
        self.name = name or self.label

        # Manual overrides, all "unset" by default so the automatic behaviour is
        # unchanged until someone deliberately reaches for them.
        self.tilt_a_deg = 0.0     # rotate the plane about its first in-plane axis
        self.tilt_b_deg = 0.0     # ...and about the second
        self.shift_a_mm = 0.0     # slide the cut anchor along that first axis...
        self.shift_b_mm = 0.0     # ...and the second, within the plane
        self.band_mm = None       # how far the cut reaches around the plane; None = auto

        # Set by classify_terminals: the two fattest ends are the trunk, the rest
        # are branches. Lets the GUI act on "all branches" without name-matching.
        # is_inlet marks the end blood ENTERS by, worked out from the anatomy rather
        # than from which end sits higher in z.
        self.is_trunk = False
        self.is_inlet = False

        # A pre-bifurcation cut is an ordinary Terminal placed on an INTERIOR
        # segment, at the end that faces the split, carrying the set of downstream
        # segment names whose whole surface it removes. None on a real opening.
        self.is_interior = False
        self.subtree_segments = None

    def plane_frame(self, normal):
        """A stable pair of in-plane axes for *normal*.

        Deterministic: the same normal always yields the same two axes, so a tilt
        angle means the same thing from one redraw to the next.
        """
        n = _unit(np.asarray(normal, float))
        seed = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(n, seed))) > 0.9:
            seed = np.array([1.0, 0.0, 0.0])
        a = _unit(np.cross(n, seed))
        b = _unit(np.cross(n, a))
        return a, b

    def reset_plane(self):
        self.tilt_a_deg = 0.0
        self.tilt_b_deg = 0.0
        self.shift_a_mm = 0.0
        self.shift_b_mm = 0.0
        self.band_mm = None

    @property
    def tip(self):
        return self.segment.endpoint(self.at_end)

    def s_cut(self, offset_mm):
        """Arc length along the parent segment at which the cut sits."""
        L = self.segment.length
        off = float(max(offset_mm, 0.0))
        return float(np.clip((L - off) if self.at_end else off, 0.0, L))

    def distal_mask(self, s_values, s_cut, margin_mm=0.0):
        """Which arc lengths lie on the cap side of *s_cut* (within a margin).

        This -- not the plane's half-space -- is what defines "past the cut". A
        tortuous vessel re-enters its own cut half-space (Line5's centerline
        crosses its own plane 3 times at large offsets), so a side-of-plane test
        cannot say what is distal. Arc length always can.
        """
        if self.at_end:
            return s_values >= (s_cut - margin_mm)
        return s_values <= (s_cut + margin_mm)

    def cut_plane(self, offset_mm, tangent_smooth_mm=None):
        """Plane *offset_mm* back from the tip, measured along the centerline.

        Returns (point, outward_normal). The normal points from the mesh toward
        the tip, i.e. toward the cap being removed.

        The tangent is averaged over a short span rather than sampled at a point:
        spline derivatives get twitchy right at an endpoint (the export shows
        second derivatives blowing up there), and a wobbly normal tilts the cut.
        """
        offset_mm = float(max(offset_mm, 0.0))
        seg = self.segment
        L = seg.length

        if tangent_smooth_mm is None:
            tangent_smooth_mm = getattr(seg, "tangent_span_mm", 4.0)

        s_cut = (L - offset_mm) if self.at_end else offset_mm
        s_cut = float(np.clip(s_cut, 0.0, L))

        half = 0.5 * float(max(tangent_smooth_mm, 1e-3))
        s_lo, s_hi = max(s_cut - half, 0.0), min(s_cut + half, L)
        n_samp = 7
        tans = np.array(
            [seg.at_arclength(s)[1] for s in np.linspace(s_lo, s_hi, n_samp)]
        )
        tan = _unit(tans.mean(axis=0))

        pos, _ = seg.at_arclength(s_cut)
        outward = tan if self.at_end else -tan

        if self.mode == "z":
            # Keep the point on the centerline; force the plane horizontal. The
            # sign comes from which way the vessel actually exits.
            outward = np.array([0.0, 0.0, 1.0 if outward[2] >= 0 else -1.0])

        outward = _unit(outward)

        # Manual tilt, applied last so it works on top of either mode.
        if self.tilt_a_deg or self.tilt_b_deg:
            a, b = self.plane_frame(outward)
            outward = _rotate(outward, a, np.radians(self.tilt_a_deg))
            outward = _rotate(outward, b, np.radians(self.tilt_b_deg))

        # NB: cut_plane returns the TRUE on-centerline point, this is where the
        # opening actually lands, since an in-plane shift never moves the plane. The
        # shifted anchor is a separate thing (cut_anchor) used only to place the
        # cut's band/footprint, so opening-location consumers (matching, rim
        # flatness, drift, live diameter) must stay on this true point.
        return pos, _unit(outward)

    def cut_anchor(self, offset_mm, tangent_smooth_mm=None):
        """The cut's band centre: the true point slid sideways WITHIN the plane.

        Same plane and slice as cut_plane (the two shift axes lie in the plane), but
        the returned point is moved by the in-plane shift. Only the band selection
        and the on-screen disc use this, everything that needs where the OPENING
        ends up must use cut_plane, because the opening lands at the true point.
        """
        pos, outward = self.cut_plane(offset_mm, tangent_smooth_mm)
        if self.shift_a_mm or self.shift_b_mm:
            a, b = self.plane_frame(outward)
            pos = pos + float(self.shift_a_mm) * a + float(self.shift_b_mm) * b
        return pos, outward


class Tree:
    def __init__(self, segments, junctions, terminals):
        self.segments = segments          # name -> Segment
        self.junctions = junctions        # list of (position, [(name, at_end), ...])
        self.terminals = terminals        # list of Terminal

    def terminal(self, label):
        for t in self.terminals:
            if t.label == label:
                return t
        raise KeyError(f"no terminal {label!r}; have {[t.label for t in self.terminals]}")

    def summary(self):
        lines = [f"Tree: {len(self.segments)} segments, "
                 f"{len(self.junctions)} junctions, {len(self.terminals)} terminals"]
        for name, s in self.segments.items():
            lines.append(f"  segment {name:8s} n={len(s):5d}  length={s.length:7.2f} mm")
        for i, (pos, inc) in enumerate(self.junctions):
            who = ", ".join(f"{n}.{'end' if e else 'start'}" for n, e in inc)
            lines.append(f"  junction J{i+1} @ {np.round(pos, 2)}  degree={len(inc)}  [{who}]")
        for t in self.terminals:
            lines.append(f"  terminal {t.label:12s} @ {np.round(t.tip, 2)}")
        return "\n".join(lines)


def _split_on_breaks(pts, tans, factor=BREAK_FACTOR):
    """Split a point run wherever the step jumps, that is a new polyline.

    Files without a 'Spline name' column can still hold several curves, laid end
    to end with nothing to mark the seam. The sampling interval is regular, so a
    step many times the median is a break, not a sample.
    """
    if len(pts) < 3:
        return [(pts, tans)]
    step = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    med = float(np.median(step))
    if med <= 0:
        return [(pts, tans)]
    cuts = np.flatnonzero(step > factor * med) + 1
    if len(cuts) == 0:
        return [(pts, tans)]
    out = []
    for a, b in zip(np.r_[0, cuts], np.r_[cuts, len(pts)]):
        if b - a >= 2:
            out.append((pts[a:b], tans[a:b]))
    return out


def _norm_col(c):
    """Column name reduced to letters and digits, lowercased, for fuzzy matching."""
    return re.sub(r"[^a-z0-9]", "", str(c).strip().lower())


# What a column has to look like to be read as a coordinate, a direction, or a
# grouping key. Written as suffix/exact sets rather than exact spellings so a file
# from any of the usual tools is understood without being told which tool wrote it.
_POS_KEYS = {
    "x": ("x", "px", "xmm", "positionx", "positionxglobal", "coordx", "pointx",
          "worldx", "xcoord", "xpos", "posx", "centerlinex"),
    "y": ("y", "py", "ymm", "positiony", "positionyglobal", "coordy", "pointy",
          "worldy", "ycoord", "ypos", "posy", "centerliney"),
    "z": ("z", "pz", "zmm", "positionz", "positionzglobal", "coordz", "pointz",
          "worldz", "zcoord", "zpos", "posz", "centerlinez"),
}
_TAN_KEYS = {
    "x": ("firstderivativex", "tangentx", "tx", "dirx", "directionx", "dx",
          "derivativex", "tangentxglobal"),
    "y": ("firstderivativey", "tangenty", "ty", "diry", "directiony", "dy",
          "derivativey", "tangentyglobal"),
    "z": ("firstderivativez", "tangentz", "tz", "dirz", "directionz", "dz",
          "derivativez", "tangentzglobal"),
}
_GROUP_KEYS = ("splinename", "branchid", "cellid", "groupid", "group", "segmentid",
               "segment", "branch", "pathid", "path", "lineid", "line", "vesselid",
               "vessel", "centerlineid", "trackid", "name", "id")


def _match_cols(columns, keys):
    """Pick the column whose normalised name matches one of *keys*, best first."""
    norm = {c: _norm_col(c) for c in columns}
    for want in keys:
        for c, n in norm.items():
            if n == want:
                return c
    return None


def derive_tangents(pts):
    """Tangents from the points themselves.

    Only POSITIONS are essential to a centerline -- arc length, tangents and cut
    planes are all computed from them. Requiring a file to carry its own derivative
    columns rejects perfectly good centerlines (a bare x,y,z export is the commonest
    thing anyone will hand us), so when they are absent we measure them over a chord
    instead, which is the same thing done for a staircased path.
    """
    pts = np.asarray(pts, float)
    if len(pts) < 2:
        return np.tile([1.0, 0.0, 0.0], (len(pts), 1))
    s = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))))
    return _chord_tangents(pts, s, CHORD_SPAN_MM)


def _segments_from_runs(runs, stem, smooth_mm):
    """Build {name: Segment} from [(name|None, pts, tans|None), ...]."""
    segments, auto = {}, 0
    for name, pts, tans in runs:
        pts = np.asarray(pts, float)
        if len(pts) >= 2:
            # Drop repeated points. Exports produce them routinely -- a stitched edge
            # soup repeats the shared vertex, and some tools emit a duplicate at a
            # join -- and a zero-length step carries no direction, so it only muddies
            # the tangents and the arc length.
            keep = np.concatenate(
                ([True], np.linalg.norm(np.diff(pts, axis=0), axis=1) > 1e-9))
            pts = pts[keep]
            tans = np.asarray(tans, float)[keep] if tans is not None else None
        if len(pts) < 2:
            continue
        if tans is None or len(tans) != len(pts):
            tans = derive_tangents(pts)
        if name is None:
            auto += 1
            name = stem if len(runs) == 1 else f"{stem}_{auto}"
        key = str(name)
        n = 2
        while key in segments:                     # never silently drop a duplicate
            key, n = f"{name}_{n}", n + 1
        segments[key] = Segment(key, pts, tans, smooth_mm=smooth_mm)
    return segments


def _read_table(path):
    """Read a delimited centerline table, sniffing the header row and delimiter."""
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"{path.name}: the file is empty")

    def delim_of(ln):
        for d in (",", ";", "\t", "|"):
            if ln.count(d) >= 2:
                return d
        return r"\s+"

    # The header is the first line that names a coordinate; a provenance comment may
    # sit above it (Simpleware writes one), and some exports have no header at all.
    header_row, delim = None, None
    for i, ln in enumerate(lines[:40]):
        d = delim_of(ln)
        cells = re.split(d, ln.strip()) if d == r"\s+" else ln.split(d)
        cells = [c.strip().strip('"') for c in cells]
        if len(cells) >= 3 and _match_cols(cells, _POS_KEYS["x"]) is not None:
            header_row, delim = i, d
            break
    if header_row is None:                          # headerless: pure numbers
        delim = delim_of(lines[0])
        df = pd.read_csv(io.StringIO("\n".join(lines)), sep=delim, header=None,
                         engine="python", comment="#")
        num = df.select_dtypes("number")
        if num.shape[1] < 3:
            raise ValueError(
                f"{path.name}: could not find x, y, z. No header naming them, and "
                f"only {num.shape[1]} numeric column(s) were found. Expected either "
                f"a header row naming the coordinates, or at least 3 numeric columns.")
        num = num.iloc[:, :3]
        num.columns = ["x", "y", "z"]
        return num

    df = pd.read_csv(io.StringIO("\n".join(lines[header_row:])), sep=delim,
                     engine="python", comment=None)
    df.columns = [str(c).strip().strip('"') for c in df.columns]
    return df


def _runs_from_table(df, path):
    """(name, pts, tans) runs from a table, grouped by whatever column groups it."""
    cols = list(df.columns)
    pos = [_match_cols(cols, _POS_KEYS[k]) for k in ("x", "y", "z")]
    if any(c is None for c in pos):
        raise ValueError(
            f"{path.name}: could not identify the x, y and z columns. Found "
            f"{cols[:12]}. Name them x/y/z (or Position x (global), Point.x, "
            f"coordX ...), or supply a file with at least 3 numeric columns.")
    tan = [_match_cols(cols, _TAN_KEYS[k]) for k in ("x", "y", "z")]
    tan = tan if all(c is not None for c in tan) else None

    grp = _match_cols(cols, _GROUP_KEYS)
    if grp is not None and df[grp].nunique(dropna=True) < 2:
        grp = None                                  # a constant column groups nothing

    def block(g):
        p = g[pos].to_numpy(float)
        t = g[tan].to_numpy(float) if tan else None
        return p, t

    runs = []
    if grp is not None:
        d = df.dropna(subset=[grp])
        for name in dict.fromkeys(d[grp]):          # keep file order
            p, t = block(d[d[grp] == name])
            runs.append((str(name), p, t))
    else:
        p, t = block(df)
        # No grouping column: several curves may still sit end to end, marked only by
        # a jump far larger than the sampling interval.
        for pp, tt in _split_on_breaks(p, t if t is not None
                                       else np.zeros_like(p)):
            runs.append((None, pp, tt if t is not None else None))
    return runs


def _runs_from_polydata(path):
    """(name, pts, None) runs from a .vtp/.vtk/.vtu polyline file."""
    import pyvista as pv

    pd_ = pv.read(str(path))
    if getattr(pd_, "n_lines", 0) == 0 and hasattr(pd_, "extract_surface"):
        pd_ = pd_.extract_surface()
    lines = np.asarray(getattr(pd_, "lines", []), dtype=np.int64)
    if lines.size == 0:
        raise ValueError(
            f"{path.name}: no polyline cells in this file, a centerline is expected "
            f"to be stored as lines, not points or triangles.")

    cells, i = [], 0
    while i < len(lines):
        n = int(lines[i])
        cells.append(lines[i + 1:i + 1 + n])
        i += 1 + n
    # Some tools write every edge as its own 2-point cell; stitch those into curves.
    if cells and np.mean([len(c) for c in cells]) < 2.5:
        try:
            pd_ = pd_.strip(join=True)
            lines = np.asarray(pd_.lines, dtype=np.int64)
            cells, i = [], 0
            while i < len(lines):
                n = int(lines[i])
                cells.append(lines[i + 1:i + 1 + n])
                i += 1 + n
        except Exception:
            pass

    pts = np.asarray(pd_.points, float)
    return [(None, pts[c], None) for c in cells if len(c) >= 2]


def _runs_from_json(path):
    """(name, pts, None) runs from 3D Slicer markups, or a plain list of points."""
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    runs = []

    def as_pts(seq):
        arr = np.asarray(seq, float)
        return arr if arr.ndim == 2 and arr.shape[1] >= 3 else None

    if isinstance(data, dict) and "markups" in data:          # Slicer .mrk.json
        for i, m in enumerate(data["markups"], 1):
            cp = m.get("controlPoints") or []
            p = as_pts([c["position"] for c in cp if "position" in c])
            if p is not None:
                runs.append((m.get("name") or f"markup_{i}", p[:, :3], None))
    elif isinstance(data, dict):
        for k, v in data.items():
            p = as_pts(v)
            if p is not None:
                runs.append((str(k), p[:, :3], None))
    elif isinstance(data, list):
        p = as_pts(data)
        if p is not None:
            runs.append((None, p[:, :3], None))
        else:                                                  # list of curves
            for i, v in enumerate(data, 1):
                q = as_pts(v)
                if q is not None:
                    runs.append((f"curve_{i}", q[:, :3], None))
    if not runs:
        raise ValueError(f"{path.name}: no point list found in this JSON.")
    return runs


TABLE_SUFFIXES = (".csv", ".tsv", ".txt", ".dat", ".asc", ".pts", ".xyz", "")
POLYDATA_SUFFIXES = (".vtp", ".vtk", ".vtu", ".obj", ".ply")


def load_spline_csv(path, smooth_mm=0.0):
    """Read a centerline into {name: Segment}, from whatever wrote it.

    The essential content of a centerline is a sequence of POSITIONS; tangents, arc
    length and cut planes are all derived. So the loader asks only for coordinates,
    and works out the rest:

      * delimited tables (.csv/.tsv/.txt/.dat) with a header naming the coordinates
        in any of the usual spellings, with a provenance comment above it or not,
        or with no header at all (the first three numeric columns);
      * polyline files (.vtp/.vtk/.vtu), including exports that store every edge as
        its own two-point cell;
      * 3D Slicer markup JSON, and plain JSON point lists.

    Several curves in one file are separated by whatever the file uses to mark them
    -- a name/branch/cell column, separate polyline cells, separate markups -- and
    failing that by a jump in the sampling interval.
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"{path.name}: file not found")
    suf = path.suffix.lower()

    if suf in POLYDATA_SUFFIXES:
        runs = _runs_from_polydata(path)
    elif suf == ".json":
        runs = _runs_from_json(path)
    else:
        runs = _runs_from_table(_read_table(path), path)

    segments = _segments_from_runs(runs, path.stem, smooth_mm)
    if not segments:
        raise ValueError(
            f"{path.name}: no usable centerline found, every curve had fewer than "
            f"two points.")
    return segments


load_centerline = load_spline_csv          # the name the format-agnostic loader earns


def build_tree(segments, tol_mm=JOIN_TOL_MM):
    """Group coincident segment endpoints into nodes; degree>=2 is a junction,
    degree 1 is a terminal (an opening to cut)."""
    ends = []                                   # (position, seg_name, at_end)
    for name, seg in segments.items():
        ends.append((seg.pts[0], name, False))
        ends.append((seg.pts[-1], name, True))

    unassigned = list(range(len(ends)))
    clusters = []
    while unassigned:
        i = unassigned.pop(0)
        group = [i]
        pi = ends[i][0]
        for j in list(unassigned):
            if np.linalg.norm(ends[j][0] - pi) <= tol_mm:
                group.append(j)
                unassigned.remove(j)
        clusters.append(group)

    junctions, terminals = [], []
    for group in clusters:
        pos = np.mean([ends[i][0] for i in group], axis=0)
        incident = [(ends[i][1], ends[i][2]) for i in group]
        if len(group) == 1:
            name, at_end = incident[0]
            terminals.append(Terminal(segments[name], at_end))
        else:
            junctions.append((pos, incident))

    return Tree(segments, junctions, terminals)


def _components(tree):
    """Group segment names into connected pieces, via shared junctions."""
    adj = {n: set() for n in tree.segments}
    for _pos, incident in tree.junctions:
        names = [n for n, _e in incident]
        for a in names:
            for b in names:
                if a != b:
                    adj[a].add(b)
    seen, comps = set(), []
    for start in adj:
        if start in seen:
            continue
        stack, comp = [start], []
        seen.add(start)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        comps.append(comp)
    return comps


def repair_disconnected(segments, tol_mm=JOIN_TOL_MM, max_weld_mm=6.0):
    """Weld an orphaned piece onto the vessel it actually branches from.

    Segments are joined by clustering COINCIDENT ENDPOINTS. A branch that leaves
    the aorta partway along it is a T-junction: its start sits near the MIDDLE of
    the parent's polyline, not near either parent endpoint, so no cluster ever
    forms and the piece is orphaned -- while looking perfectly continuous in
    ScanIP. Measured on real exports: JK3_1 has a branch whose end is 0.11 mm from
    the aorta and 12.55 mm from the nearest endpoint. No join tolerance can fix
    that; the geometry is a T, not a gap.

    So: split the parent where the orphan meets it, and snap the orphan's endpoint
    onto the split, which ordinary clustering then joins.

    This is deliberately NOT the same as raising tol_mm. Raising the tolerance acts
    everywhere, and at the 1.7 mm one real case needs it merges a genuine distal
    aortic opening into a junction (11 openings become 9). This only ever touches
    a tree that is ALREADY in more than one piece, and only the orphaned side, so
    it cannot weld two openings of a healthy tree together.

    Returns (segments, notes); notes is empty when nothing was done.
    """
    notes = []
    for _round in range(8):                       # each pass welds one orphan
        tree = build_tree(segments, tol_mm=tol_mm)
        comps = _components(tree)
        if len(comps) <= 1:
            break
        comps.sort(key=len, reverse=True)
        main = set(comps[0])
        # every point of the main piece, so an orphan can meet it mid-vessel
        host_pts, host_of, host_idx = [], [], []
        for name in main:
            seg = segments[name]
            for i, p in enumerate(seg.pts):
                host_pts.append(p)
                host_of.append(name)
                host_idx.append(i)
        if not host_pts:
            break
        kd = cKDTree(np.asarray(host_pts))

        best = None                               # (dist, orphan, at_end, host, i)
        for comp in comps[1:]:
            for name in comp:
                seg = segments[name]
                for at_end, p in ((False, seg.pts[0]), (True, seg.pts[-1])):
                    d, k = kd.query(p)
                    if best is None or d < best[0]:
                        best = (float(d), name, at_end, host_of[k], host_idx[k])
        if best is None or best[0] > max_weld_mm:
            break

        d, oname, at_end, hname, hi = best
        host = segments[hname]
        # Splitting close to a host END leaves a sliver, and a sliver is worse than
        # no split: the cut depth then exceeds its own vessel ("trimmed 142% of the
        # way down"). Near an end this is just a gap the tolerance missed, so snap
        # to that end and leave the host whole.
        near_start = float(host.s[hi])
        near_end = float(host.s[-1] - host.s[hi])
        if 0 < hi < len(host.pts) - 1 and min(near_start, near_end) >= MIN_SPLIT_MM:
            a_pts, a_tan = host.pts[:hi + 1], host.tangents[:hi + 1]
            b_pts, b_tan = host.pts[hi:], host.tangents[hi:]
            if len(a_pts) >= 2 and len(b_pts) >= 2:
                del segments[hname]
                segments[f"{hname}_A"] = Segment(f"{hname}_A", a_pts, a_tan)
                segments[f"{hname}_B"] = Segment(f"{hname}_B", b_pts, b_tan)
                notes.append(f"{hname} split at {d:.2f} mm to receive {oname}")
            target = np.asarray(host.pts[hi], float)
        else:
            hi = 0 if near_start <= near_end else len(host.pts) - 1
            target = np.asarray(host.pts[hi], float)
            d = float(np.linalg.norm(
                target - (segments[oname].pts[-1 if at_end else 0])))
            if d > max_weld_mm:
                break                       # too far to be the same junction

        orphan = segments[oname]
        pts = np.array(orphan.pts, float)
        pts[-1 if at_end else 0] = target          # snap the loose end onto the host
        segments[oname] = Segment(oname, pts, orphan.tangents)
        notes.append(f"{oname}.{'end' if at_end else 'start'} welded to {hname} "
                     f"across {d:.2f} mm")
    return segments, notes


def load_tree(path, tol_mm=JOIN_TOL_MM, smooth_mm=0.0):
    """Load a centerline. *smooth_mm* > 0 de-staircases it; 0 leaves it untouched."""
    segs = load_spline_csv(path, smooth_mm=smooth_mm)
    segs, repairs = repair_disconnected(segs, tol_mm=tol_mm)
    tree = build_tree(segs, tol_mm=tol_mm)
    tree.repairs = repairs
    tree.smooth_mm = float(smooth_mm)
    tree.jagged = any(s.jagged for s in tree.segments.values())
    tree.worst_turning = max(s.raw_turning for s in tree.segments.values())
    return tree


if __name__ == "__main__":
    import sys

    csv = sys.argv[1] if len(sys.argv) > 1 else "KY160_1_CentrelinePoints.csv"
    print(load_tree(csv).summary())
