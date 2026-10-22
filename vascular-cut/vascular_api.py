# Headless orchestration for the cut tool.
#
# branch_cut.py holds the geometry; branch_gui.py holds the Qt app that drives it.
# Everything in between -- which terminals are being cut, at what offset, whether a
# region is isolated, what the sidecar says -- lived only inside the widget, so the
# only way to cut a mesh was to open a window.
#
# This module is that middle layer with no Qt in it, so a web request, a script or a
# batch run can do exactly what the window does. The window remains the reference
# implementation: where the two could drift, this file follows branch_gui.CutWindow.

import io
import json
import re
from pathlib import Path

import numpy as np
import trimesh

from branch_cut import (
    PAIR_AREA_FRACTION,
    ANATOMY_FRAME,
    RIM_FLAW_MM,
    aorta_route,
    band_width,
    boundary_loops,
    _component,
    vessel_attribution,
    centerline_arc_angle,
    centerline_defects,
    classify_terminals,
    cut_terminal,
    keep_largest_component,
    load_mat_mesh,
    match_openings_to_terminals,
    measure_along_branch,
    min_triangle_quality,
    MIN_OFFSET_MM,
    anatomy_choices,
    name_anatomy,
    set_manual_name,
    pre_bifurcation_candidates,
    preflight,
    region_plan,
    fit_region_reach,
    fit_region,
    region_presets,
    rim_sliver_count,
    flat_trunk_offset,
    clear_bound_offset,
    is_model_end,
    section_diameter,
    segment_radius,
    snap_tolerance,
    suggest_pullback,
    validate_cut,
)
from branch_tree import load_tree, Terminal

DEFAULT_OFFSET_MM = 6.0

# What the window puts in its smoothing box, and it ticks the box itself when the
# centerline reads as jagged. Matching both defaults matters: a centerline smoothed
# differently gives slightly different cut planes, so a session that quietly skipped
# smoothing would not reproduce the same mesh the window produced.
DEFAULT_SMOOTH_MM = 1.0

# Kept in step with branch_gui: delta K is extensive, so a curvature comparison is
# only meaningful when cut extent, face count, trim and topology all match the
# reference cohort. The sidecar records ours rather than leaving them to be guessed.
CAMP_TARGET_FACES = 100_000
CAMP_REF_AREA_MM2 = (26_252.0, 98_459.0)

FAMILY = {"auto": None, "abdominal": 1, "thoracic": 2}

# How near the cut-off counts as a CLOSE call worth reporting. Matches the desktop
# window, which has reported these since the family rule went in.
FAMILY_CLOSE_BAND = 0.08

# Matches branch_gui: past 60 degrees a cut plane is so oblique to its vessel that
# the rim tears rather than closing, and the validator refuses it anyway.
MAX_TILT_DEG = 60.0
MIN_BAND_MM = 3.0
# A cut may not travel past its own vessel. The honest limit is per-vessel, so this
# is a fraction rather than a distance -- same value the window uses.
MAX_FRAC_OF_SEGMENT = 0.95

CUT_MODES = {"z": "Z-plane (horizontal)", "normal": "Centerline-normal"}


class CutError(RuntimeError):
    """Raised when a session cannot be built at all (bad file, no centerline)."""


class CutSession:
    """One mesh + centerline, and the cutting decisions made about them.

    Construct it, read `openings`, change offsets or isolate a region, call
    `recut()`, then read `problem` / `mesh` / `sidecar()`. Nothing renders and
    nothing is written until you ask for it.
    """

    def __init__(self, centerline, mesh, family="auto", smooth_mm=None,
                 default_offset_mm=DEFAULT_OFFSET_MM, clean_rims=False):
        self.csv_path = Path(centerline)
        self.mesh_path = Path(mesh)
        self.family = family
        self.default_offset = float(default_offset_mm)
        # Off by default because that is what the window's checkbox does, and this
        # module is meant to be its headless twin rather than a second opinion.
        # Worth turning on for curvature work: on KY160_1 it takes the rim from 3
        # slivers to 0 and the worst triangle quality from 0.004 to 0.06, and a
        # near-degenerate triangle is exactly where a curvature fit misbehaves.
        self.clean_rims = bool(clean_rims)

        self._offsets = {}
        self._trims = []          # pre-bifurcation cuts, each standing for a group
        # Tilt / shift / plane size set on a REGION BOUNDARY, keyed by the boundary's
        # role. Those Terminals are rebuilt by region_plan on every cut, so this is
        # the only place such a setting can survive. See _region_terminals.
        self._bound_planes = {}
        self.region = None
        self.region_name = ""
        self._region_manual = False
        self._region_clear_start = True
        self._preset_start_mm = None
        self.mesh = None
        self.problem = None
        self.matched = {}

        self._load(smooth_mm)
        self.recut()

    # ── loading ──────────────────────────────────────────────────────────────
    def _load(self, smooth_mm):
        suffix = self.mesh_path.suffix.lower()
        try:
            if suffix == ".mat":
                self.base = load_mat_mesh(self.mesh_path)
            else:
                # From the BYTES, never the path. Handed a path, trimesh builds a
                # resolver that opens files NAMED INSIDE the mesh -- an OBJ's mtllib, a
                # PLY's TextureFile -- joined onto the mesh's folder with no traversal
                # check. So an uploaded file saying `mtllib ../../x` makes the server
                # open x. From a stream there is no folder to resolve against and
                # nothing is opened; checked with an audit hook on every open().
                # Passing the type explicitly also keeps the server's own path out of
                # trimesh's error text, which a hosted page would show to a stranger.
                self.base = trimesh.load(io.BytesIO(self.mesh_path.read_bytes()),
                                         file_type=suffix.lstrip("."), process=False)
                if suffix == ".stl" and isinstance(self.base, trimesh.Trimesh):
                    # STL stores every triangle's three corners separately, so read
                    # without processing it is a triangle soup, not a surface: KY160
                    # comes out as 119,406 vertices and 39,802 "openings" instead of
                    # 19,903 and five. Merged, it cuts exactly as the .mat does. No
                    # suite had ever cut an STL, and STL is the format a stranger has.
                    self.base.merge_vertices()
        except Exception as e:
            raise CutError(f"Could not read the mesh: {e}") from e
        # A file trimesh cannot make sense of comes back as an empty Scene rather than
        # an error, and the failure then surfaces much later as an AttributeError that
        # means nothing to whoever uploaded it. Say it here, in words.
        if not isinstance(self.base, trimesh.Trimesh) or len(self.base.faces) == 0:
            raise CutError("Could not read the mesh: the file holds no triangle surface.")

        try:
            # Load raw first and measure it as supplied, exactly as the window does,
            # so the smoothing decision is made against the real centerline.
            probe = load_tree(self.csv_path, smooth_mm=0.0)
            self.jagged = bool(probe.jagged)
            # smooth_mm=None means "do what the window would do"; pass a number to
            # override, including 0.0 to force the centerline through untouched.
            smooth = (DEFAULT_SMOOTH_MM if self.jagged else 0.0) \
                if smooth_mm is None else float(smooth_mm)
            self.tree = probe if smooth <= 0 else load_tree(self.csv_path,
                                                            smooth_mm=smooth)
            self.smooth_mm = float(smooth)
        except Exception as e:
            raise CutError(f"Could not read the centerline: {e}") from e

        if self.family not in FAMILY:
            raise CutError(f"family must be one of {sorted(FAMILY)}")
        classify_terminals(self.tree, self.base, n_trunk_ends=FAMILY[self.family],
                           verbose=False)
        self.anat_names, self.anat_vessels, self.anat_notes = name_anatomy(
            self.tree, self.base)
        for t in self.tree.terminals:
            if t.label in self.anat_names:
                t.name = self.anat_names[t.label]
        self.route = aorta_route(self.tree)

    # ── what is being cut ────────────────────────────────────────────────────
    def _region_terminals(self):
        """(bounds, inside, notes) for the active region, or ([], [], []) if none."""
        if self.region is None or self.route is None:
            return [], [], []
        _a, _b = self.region
        s0, s1 = ((float(_a), None) if _b is None
                  else (min(_a, _b), max(_a, _b)))
        # A boundary the caller chose is theirs: clearing ostia would move it behind
        # their back. A preset keeps the automatic clearance, since no station was
        # picked by hand.
        bounds, inside, _removed, notes = region_plan(
            self.tree, s0, s1, self.route, mesh=self.base,
            clear_ostia=not self._region_manual,
            clear_start=getattr(self, "_region_clear_start", True))
        # A reach the fitter settled on for THIS region, reapplied to the fresh
        # Terminal objects region_plan builds on every call.
        scale = getattr(self, "_region_reach", None)
        if scale:
            for t in bounds:
                if not getattr(t, "band_mm", None):
                    t.band_mm = float(band_width(t)) * float(scale)
            notes = list(notes) + [
                f"boundary plane reach set to {scale:.2f}x the automatic size, "
                f"which is what cuts here"]
        # A nudge settled on for THIS region, reapplied to the fresh Terminals
        # region_plan builds on every call. Mutating the objects instead does nothing
        # at all: they are rebuilt from scratch here, which is the same trap that ate
        # hand-set boundary tilts (tests/test_region_planes.py).
        nudge = getattr(self, "_bound_nudge", None)
        if nudge:
            for t in bounds:
                if t.name in nudge:
                    t._region_offset = float(nudge[t.name])
        # Tilt, shift and plane size the CALLER set on a boundary, reapplied for the
        # same reason as the two above: region_plan rebuilt these objects, so a value
        # written straight onto one is gone by the next re-cut. The window has carried
        # these since tests/test_region_planes.py -- which is a CutWindow test, so it
        # never covered this module, and set_plane on a boundary silently did nothing
        # here. The web draws those controls, so they were live and inert.
        #
        # Keyed by the bound's ROLE, not its label: the label carries whichever
        # segment the plane happens to land on, and that changes under the caller as
        # the boundary moves along the aorta.
        saved = getattr(self, "_bound_planes", None) or {}
        for t in bounds:
            for attr, val in saved.get(t.name, {}).items():
                setattr(t, attr, val)
        return bounds, inside, notes

    def _region_bound_role(self, label):
        """The role name of the region boundary with this label, or None.

        A boundary Terminal is not in tree.terminals -- it is manufactured by
        region_plan for this cut only -- which is exactly what makes it the thing
        that has to be remembered on the session instead.
        """
        if self.region is None:
            return None
        if any(t.label == label for t in self.tree.terminals):
            return None
        for t in self._region_terminals()[0]:
            if t.label == label:
                return t.name
        return None

    def _visible_terminals(self):
        """The tree's openings, with any trimmed group replaced by its one cut.

        A trim is one cut on the vessel FEEDING a group -- the brachiocephalic trunk
        before it splits, a common iliac before internal and external -- and it
        stands in for every opening it swallows. Cutting the tree's terminals as
        well would re-cut branches that no longer exist.
        """
        swallowed = set()
        for tr in self._trims:
            swallowed |= set(tr["subtree"])
        out = [t for t in self.tree.terminals if t.segment.name not in swallowed]
        return out + [tr["terminal"] for tr in self._trims]

    def _trims_within(self, terminals):
        """The trims whose whole group lies inside *terminals*.

        Mirrors branch_gui.CutWindow.trims_within, but keyed on the SEGMENTS a trim
        swallows rather than on the terminal labels the window carries: a trim here is
        {"segment", "subtree", "terminal"} and has no candidate record to read the
        group's labels back out of.

        A trim only means something where every opening it swallows is present; one
        offered while looking at another stretch of the vessel would remove openings
        that are not in view.
        """
        segs = {t.segment.name for t in terminals}
        return [tr for tr in self._trims
                if tr["subtree"] and set(tr["subtree"]) <= segs]

    def _apply_trims(self, terminals, trims):
        """Replace each trim's group of openings with its single interior cut."""
        gone = set()
        for tr in trims:
            gone |= set(tr["subtree"])
        return ([t for t in terminals if t.segment.name not in gone]
                + [tr["terminal"] for tr in trims])

    def _cut_list(self):
        """(terminals to cut, offsets) -- the single source of truth for a re-cut."""
        if self.region is None:
            active = self._visible_terminals()
            offs = {t.label: float(self._offsets.get(t.label, self.default_offset))
                    for t in active}
            return active, offs

        bounds, inside, _ = self._region_terminals()
        # TRIMS STILL APPLY INSIDE A REGION, and this file had drifted from the window
        # on exactly that. branch_gui.region_rows() replaces each trimmed group with
        # its single interior cut before the region's rows are built -- it says so, and
        # records that a region used to drop them entirely, so trimming the iliacs
        # while looking at the infrarenal aorta recorded the trim and changed nothing.
        # Here the same call was simply missing, so the API cut a DIFFERENT MESH from
        # the window for the same session: every trimmed group came back as separate
        # openings. This module follows the window where the two could drift, which is
        # the whole contract at the top of this file.
        #
        # A trim only counts where every opening it swallows is inside the region: one
        # offered while looking elsewhere would remove openings that are not on screen.
        inside = self._apply_trims(inside, self._trims_within(inside))
        active = list(inside) + list(bounds)
        offs = {t.label: float(self._offsets.get(t.label, self.default_offset))
                for t in inside}
        # A region bound is not a vessel end, so the default offset means nothing to
        # it -- region_plan already worked out where each bound has to sit to clear
        # the nearby ostia. Overriding that with 6 mm puts the cut somewhere else and
        # leaves a stray opening behind.
        for t in bounds:
            offs[t.label] = float(getattr(t, "_region_offset", self.default_offset))
        return active, offs

    # ── controls ─────────────────────────────────────────────────────────────
    def set_offset(self, label, mm):
        self._offsets[str(label)] = float(mm)

    # ── trimming a branch group ──────────────────────────────────────────────
    #
    # One cut on the vessel feeding a group, placed before it splits, leaving ONE
    # opening where there were several: the brachiocephalic trunk ahead of the right
    # carotid and subclavian, a common iliac ahead of internal and external.

    def trim_candidates(self):
        """Vessels that can be cut before they split, and what each would swallow."""
        if self.tree is None or self.base is None:
            return []
        taken = set()
        for tr in self._trims:
            taken |= set(tr["subtree"]) | {tr["segment"]}
        # Inside a region, only offer trims whose WHOLE group is in view -- the same
        # rule branch_gui._refresh_trim_menu applies, and for the same reason it gives:
        # the menu used to list every candidate in the tree, so the visceral segment
        # offered to trim the iliacs, which are not part of what is being cut.
        #
        # It is not cosmetic here. _cut_list only applies trims via _trims_within, so
        # a trim taken on a group outside the region was RECORDED and then silently
        # dropped: the caller got a success, a trim in the session, and no change to
        # the mesh. That is the same fault the cohort reported against the region
        # sliders -- a control that reports success and does nothing.
        visible = None
        if self.region is not None:
            _bounds, inside, _notes = self._region_terminals()
            visible = {t.label for t in inside}
        out = []
        for c in pre_bifurcation_candidates(self.tree, self.base):
            if c["segment"] in taken:
                continue
            if visible is not None and not {t.label for t in c["terminals"]} <= visible:
                continue
            names = [t.name for t in c["terminals"]]
            group = (" + ".join(names) if len(names) <= 3
                     else f"{names[0]} + {len(names) - 1} more")
            vessel = self.anat_vessels.get(c["segment"], c["segment"])
            out.append({"segment": c["segment"], "vessel": vessel, "group": group,
                        "length_mm": round(float(self.tree.segments[c["segment"]].length), 1),
                        "n_openings": len(names)})
        return out

    def add_trim(self, segment):
        """Replace a group of openings with one cut on the vessel feeding them."""
        segment = str(segment)
        cand = next((c for c in pre_bifurcation_candidates(self.tree, self.base)
                     if c["segment"] == segment), None)
        if cand is None:
            raise CutError(f"{segment} cannot be cut before it splits.")
        seg = self.tree.segments[segment]
        if getattr(seg, "_radius_mm", None) is None:
            seg._radius_mm = segment_radius(self.base, seg)   # band width needs it
        t = Terminal(seg, cand["at_end"], label=f"{seg.name}.trim", mode="normal",
                     name=self.anat_vessels.get(segment) or f"Trim before {segment}")
        t.is_interior = True
        t.subtree_segments = set(cand["subtree"])

        # Open it where it actually cuts. A trim sits on a short interior segment and
        # the usual default is often already past the point where its parent takes
        # over, so the first thing the user would see is a refusal on a cut they were
        # just told they could make.
        start = min(self.default_offset, max(seg.length * 0.5, MIN_OFFSET_MM))
        try:
            found = suggest_pullback(self.base, self.tree, t, start)
        except Exception:                                       # noqa: BLE001
            found = None
        self._offsets[t.label] = float(found[0]) if found else start

        # A trim already inside another trim's removed subtree is redundant.
        self._trims = [tr for tr in self._trims
                       if tr["segment"] not in cand["subtree"]]
        self._trims.append({"segment": segment, "subtree": set(cand["subtree"]),
                            "terminal": t})
        self._cut_steps = None          # the cut list changed shape
        return t.label

    def remove_trim(self, segment):
        """Undo a trim: the openings it swallowed come back as individual cuts."""
        segment = str(segment)
        before = len(self._trims)
        self._trims = [tr for tr in self._trims if tr["segment"] != segment]
        if len(self._trims) == before:
            raise CutError(f"{segment} is not trimmed.")
        self._cut_steps = None
        return True

    def pull_back(self, label):
        """Move one opening to the deepest offset that still cuts cleanly.

        Returns (offset_mm, note). "Pull it back" is the right advice with no number
        attached; the cutter can work the number out by trying, and the search it uses
        is the real rule -- a candidate is accepted only when the rim it produces lies
        on its plane and the cut takes no more than a square millimetre of any vessel
        it is not cutting. So this cannot offer something the cut will then refuse.

        Never moves an opening CLOSER to the junction, and leaves one alone if it is
        already at the deepest clean cut.
        """
        t = self.terminal(label)
        cur = float(self._offsets.get(str(label), DEFAULT_OFFSET_MM))
        found = suggest_pullback(self.base, self.tree, t, cur)
        if found is None:
            raise CutError(
                f"{t.name}: no offset along this vessel has a clean cross-section of "
                f"its own, so there is nowhere to pull back to.")
        off, frac = found
        if off >= cur - 0.05:
            return cur, f"{t.name} is already at the deepest clean cut ({cur:.0f} mm)."
        self.set_offset(label, off)
        return off, (f"{t.name} pulled back to {off:.0f} mm ({frac * 100:.0f}% along "
                     f"the vessel), the deepest cut here that leaves a flat opening.")

    # ── the cut plane itself ─────────────────────────────────────────────────
    #
    # Depth alone does not rescue every cut. A plane that meets a tortuous branch
    # at the wrong angle, or reaches further around than the vessel is wide, punches
    # a second hole no matter where along the vessel you put it. Tilt, shift and
    # plane size are the remedies for exactly that, and without them the web tool
    # can only cut the cases that happened to work first time.

    def terminal(self, label):
        for t in self.tree.terminals:
            if t.label == label:
                return t
        for t in self._cut_list()[0]:
            if t.label == label:
                return t
        raise CutError(f"no opening labelled {label!r}")

    def plane_limits(self, label):
        """The ranges the controls may offer, derived per opening like the window."""
        t = self.terminal(label)
        auto_band = float(band_width(t))
        r = max(float(getattr(t.segment, "_radius_mm", None) or 5.0), 1.0)
        # Past roughly (reach - radius) the band stops covering the cut ring at all,
        # so offering more would just be offering failure.
        shift_range = max(auto_band - r, 6.0)
        return {
            "tilt": MAX_TILT_DEG,
            "shift": round(shift_range, 2),
            "band_min": MIN_BAND_MM,
            "band_max": round(max(4.0 * auto_band, 30.0), 2),
            "band_auto": round(auto_band, 2),
        }

    def plane_controls(self, label):
        t = self.terminal(label)
        lim = self.plane_limits(label)
        return {
            "tilt_a_deg": round(float(t.tilt_a_deg), 2),
            "tilt_b_deg": round(float(t.tilt_b_deg), 2),
            "shift_a_mm": round(float(t.shift_a_mm), 2),
            "shift_b_mm": round(float(t.shift_b_mm), 2),
            "band_mm": round(float(band_width(t)), 2),
            "band_auto": t.band_mm is None,
            "mode": getattr(t, "mode", "z"),
            **{f"limit_{k}": v for k, v in lim.items()},
        }

    def set_plane(self, label, tilt_a_deg=None, tilt_b_deg=None,
                  shift_a_mm=None, shift_b_mm=None, band_mm=None):
        """Set any subset of the plane controls. Values are clamped to the same
        ranges the window offers rather than trusted."""
        t = self.terminal(label)
        lim = self.plane_limits(label)
        changed = {}
        if tilt_a_deg is not None:
            changed["tilt_a_deg"] = float(
                np.clip(tilt_a_deg, -lim["tilt"], lim["tilt"]))
        if tilt_b_deg is not None:
            changed["tilt_b_deg"] = float(
                np.clip(tilt_b_deg, -lim["tilt"], lim["tilt"]))
        if shift_a_mm is not None:
            changed["shift_a_mm"] = float(
                np.clip(shift_a_mm, -lim["shift"], lim["shift"]))
        if shift_b_mm is not None:
            changed["shift_b_mm"] = float(
                np.clip(shift_b_mm, -lim["shift"], lim["shift"]))
        if band_mm is not None:
            changed["band_mm"] = float(
                np.clip(band_mm, lim["band_min"], lim["band_max"]))
        for attr, val in changed.items():
            setattr(t, attr, val)
        # A region BOUNDARY is rebuilt from scratch on every re-cut, so setting it on
        # the object alone is a no-op the caller cannot see. Record it against the
        # boundary's role; _region_terminals puts it back. See that method.
        role = self._region_bound_role(label)
        if role and changed:
            self._bound_planes.setdefault(role, {}).update(changed)

    def reset_plane(self, label):
        """Back to the automatic plane: no tilt, no shift, automatic size."""
        t = self.terminal(label)
        t.reset_plane()
        t.band_mm = None
        role = self._region_bound_role(label)
        if role:
            # Drop the record too, or _region_terminals would reapply what was
            # just reset and the control would appear dead in the other direction.
            self._bound_planes.pop(role, None)

    def set_mode(self, label, mode):
        if mode not in CUT_MODES:
            raise CutError(f"mode must be one of {sorted(CUT_MODES)}")
        self.terminal(label).mode = mode

    def rename(self, label, name):
        """The automatic anatomical name is a guess; let it be corrected."""
        self.terminal(label).name = str(name).strip() or label

    # ── naming a vessel by hand ──────────────────────────────────────────────
    #
    # rename() above is a LABEL: it changes what an opening is called and nothing
    # else. Declaring an artery is different -- region_presets matches on the words
    # in these names, so a declared renal is what brings back the aortic neck and the
    # infrarenal aorta. Over 47 abdominal scans name_anatomy names no renal on 18,
    # and on 14 of those because several lateral pairs fit equally well and choosing
    # between them would be a guess presented as a fact. See tools/ABDOMINAL.md.

    def anatomy_choices(self):
        """The artery names that can be declared on this anatomy."""
        return list(anatomy_choices(self.tree))

    def declare_artery(self, label, name):
        """Say which artery an opening IS. A falsy *name* hands it back to the rules.

        Re-derives every name, because a declaration displaces the automatic name of
        the same artery -- without that, region_presets would see two "renal"
        landmarks and place the aortic neck at whichever sat lower.
        """
        label = str(label)
        self.terminal(label)                          # raises if it is not ours
        choices = set(self.anatomy_choices())
        if name and str(name) not in choices:
            raise CutError(f"{name!r} is not an artery this anatomy can carry; "
                           f"use rename() for a free label")
        set_manual_name(self.tree, label, name or None)
        self.anat_names, self.anat_vessels, self.anat_notes = name_anatomy(
            self.tree, self.base)
        for t in self.tree.terminals:
            if t.label in self.anat_names:
                t.name = self.anat_names[t.label]
        self.route = aorta_route(self.tree)
        return dict(self.anat_names)

    def declared_arteries(self):
        """{label: artery} for the openings named by hand rather than by the rules."""
        return dict(getattr(self.tree, "manual_names", None) or {})

    def reset_all(self):
        """Back to the model as loaded: no region, no trims, no hand names, no planes.

        The same scope branch_gui.CutWindow.reset uses, and its reasoning applies
        here verbatim: "Reset is where a stuck user goes... someone who presses Reset
        wants the model they loaded, and leaving three branches collapsed into one
        makes the button look broken in precisely the situation it exists for."

        Lives here rather than in either front end because the two had already
        drifted -- the web's reset cleared offsets, planes and the region but left
        trims and hand-declared arteries standing, so the same button did two
        different things depending on which one you were using.
        """
        self.clear_region()
        self._trims = []
        self._offsets.clear()
        for t in self.tree.terminals:
            self.reset_plane(t.label)
        if getattr(self.tree, "manual_names", None):
            self.tree.manual_names = {}
            try:
                self.anat_names, self.anat_vessels, self.anat_notes = name_anatomy(
                    self.tree, self.base)
                for t in self.tree.terminals:
                    if t.label in self.anat_names:
                        t.name = self.anat_names[t.label]
                self.route = aorta_route(self.tree)
            except Exception:                                   # noqa: BLE001
                # A reset that cannot re-derive names is still a useful reset; the
                # window swallows this the same way rather than refusing to reset.
                pass
        self.recut()
        return True

    # ── the nub-in preset ────────────────────────────────────────────────────

    def nub_preset(self, k=None):
        """Cut every branch back to a short stub, and square off the aorta's ends.

        Each stub ends 0.75 of the AORTA's local radius past the aortic wall, where
        the branch leaves it (nub_cut.py: n = s_wall + k * D_aorta/2) -- measured over
        82 vessels in 30 thoracic cases, not chosen. An earlier wording here said "each
        vessel's own radius"; the rule has always been the aorta's. Returns the notes, which
        matter: a nub that would not cut cleanly is moved OUTWARD until it does, and
        a trunk end too oblique for a horizontal plane falls back to perpendicular.
        Both leave a case worth looking at, and neither is visible in the result.

        Refused while a region is isolated: the preset squares off the aorta's own
        ends, and a region has already replaced those with its own bounds.

        Thoracic aortas only. That is where it was measured and what it is for: the
        arch with its great vessels left as nubs. An abdominal case is refused rather
        than cut by a rule nobody has checked there.
        """
        import nub_cut as NC

        if self.tree is None or self.base is None:
            raise CutError("load a mesh and a centerline first")
        if not self.is_thoracic:
            raise CutError("the nub-in preset is for thoracic aortas, and this case "
                           "reads as abdominal. If that call is wrong, set the anatomy "
                           "to thoracic and load it again.")
        if self.region is not None:
            raise CutError("clear the region first; the nub preset cuts the aorta's "
                           "own ends, which a region has already replaced")
        notes = [f"split-line repair: {nm} radius {was:.2f} -> {now:.2f} mm"
                 for nm, was, now, _s in NC.heal_split_segments(self.tree, self.base)]
        try:
            _out, _plan, cut_notes, info = NC.cut_nubs(
                self.base, self.tree, k=NC.NUB_K if k is None else float(k),
                anat=dict(self.anat_vessels or {}))
        except Exception as e:                                   # noqa: BLE001
            raise CutError(f"the nub preset failed: {e}") from e
        notes += list(cut_notes)

        # Interior nubs are trims: one cut on a feeding vessel standing for the group
        # past it. They have to go in first, since adding one changes what is offered.
        for t in info["terminals"]:
            if not getattr(t, "is_interior", False):
                continue
            seg = t.segment.name
            if any(tr["segment"] == seg for tr in self._trims):
                continue
            try:
                self.add_trim(seg)
            except CutError:
                notes.append(f"{t.name}: {seg} is not offered as a trim here, its "
                             f"nub is left as separate openings")

        placed = 0
        for t in info["terminals"]:
            off = info["offsets"].get(t.label)
            if off is None:
                continue
            # An interior nub is called "<seg>.nub" here and "<seg>.trim" once the
            # trim exists, so match on the segment rather than on the label.
            mine = next((x for x in self._cut_list()[0]
                         if x.label == t.label
                         or (getattr(x, "is_interior", False)
                             == bool(getattr(t, "is_interior", False))
                             and x.segment.name == t.segment.name)), None)
            if mine is None:
                notes.append(f"{t.name}: no opening here to place it on")
                continue
            mine.mode = t.mode
            if getattr(t, "band_mm", None) is not None:
                mine.band_mm = float(t.band_mm)
            for attr in ("tilt_a_deg", "tilt_b_deg"):
                val = float(getattr(t, attr, 0.0) or 0.0)
                if abs(val) > 1e-9:
                    setattr(mine, attr, val)
            self.set_offset(mine.label, float(off))
            placed += 1

        self.recut()
        return [f"nub-in cut applied: {placed} opening(s) placed"
                f"{f', {len(self._trims)} trimmed group(s)' if self._trims else ''}"
                ] + notes

    def set_offsets(self, mapping):
        for k, v in dict(mapping).items():
            self.set_offset(k, v)

    def presets(self, planar=False):
        """Standard regions this anatomy can actually support.

        TWO CALLERS WANT TWO DIFFERENT LISTS, and conflating them breaks one of them.

        The zone work asks what zones the ANATOMY has: it labels by station, so a zone
        exists wherever its two landmarks do, whether or not a flat cross-section
        happens to bound it. That is the default, and region_station's labelling reads
        this list to decide which zones to draw -- filter it and zones 1 and 2 simply
        vanish from the map on most arches.

        A CUT wants what it can actually cut. Pass planar=True for that: measured over
        six cohort arches, the unfiltered list offers Ishimaru zones 1 and 2 nine more
        times and apply_preset refuses every one of them, and a preset a user can pick
        and cannot cut is worse than one never offered. Both front ends and
        apply_preset pass it.
        """
        return region_presets(self.tree, self.base, self.anat_names,
                              self.anat_vessels, planar=bool(planar))

    def set_region(self, s0, s1, name="", manual=True, clear_start=True):
        """Isolate the vessel between two stations, in mm along the aortic route.

        `manual=True` means the stations came from a person, so region_plan is told
        to leave them alone rather than nudge them clear of a branch ostium.
        Presets pass manual=False.

        `clear_start=False` asks for the PROXIMAL station exactly as given while the
        distal one keeps its clearance. apply_preset uses it to walk a ladder of
        ascending stations; nothing else should need it.
        """
        # s1 None is an open distal end: one plane, everything below it kept. That
        # is the only way to keep the iliacs, whose ostia sit at the bifurcation and
        # therefore fall outside any two-plane region.
        self._region_reach = None      # a new region needs its own answer
        self._region_moved_mm = 0.0
        self.region = (float(s0), None if s1 is None else float(s1))
        self.region_name = name or f"{s0:.0f}-{s1:.0f} mm"
        self._region_manual = bool(manual)
        self._region_clear_start = bool(clear_start)

    def apply_preset(self, name):
        """Take a preset, and if it carries alternatives, the first that CUTS.

        Mirrors branch_gui.CutWindow._apply_preset, which this had drifted from. The
        aortic arch's distal end is a convention -- the last great vessel plus 20 mm --
        and region_presets says in as many words that on half of twenty real arches
        that convention lands inside the hairpin, where there is no planar
        cross-section to cut. So it carries s1_options, "and let the caller take the
        first that actually cuts, saying how far it had to go -- rather than silently
        redefining the arch for everyone."

        This caller never looked. The API took s1 and refused, so KK20_2 and KK4_2
        came back "Region end: cannot be cut here" from a preset that was already
        offering the answer.

        How far it went is recorded rather than absorbed. Where the arch ENDS is a
        judgement about anatomy, not about cutting, so a span the tool had to extend
        must be visible to whoever reads the result.
        """
        for p in self.presets(planar=True):
            if p["name"] != name:
                continue
            options = list(p.get("s1_options") or [p["s1"]])
            if p["s1"] not in options:
                options.insert(0, p["s1"])
            for i, s1 in enumerate(options):
                self.set_region(p["s0"], s1, p["name"], manual=False)
                st = self.recut()
                if st["ok"]:
                    self._preset_extended_mm = (
                        None if i == 0 else round(float(s1) - float(p["s1"]), 2))
                    self._settle_start(p, s1)
                    return p
            # Nothing cut. Left on the conventional bound with its own complaint,
            # which is the honest answer rather than the least-bad arbitrary span.
            self.set_region(p["s0"], p["s1"], p["name"], manual=False)
            self._preset_extended_mm = None
            return p
        raise CutError(f"no preset named {name!r}")

    def _settle_start(self, p, s1):
        """With the distal end settled, bring the ASCENDING bound in to the target.

        Two phases rather than one nested search: the distal ladder is walked first
        with the default start, then the start is walked with that distal station
        held. Nesting them is 10 x 18 recuts for the same answer, and the two bounds
        do not interact -- an ascending station that cuts, cuts whichever end the
        arch was given.

        The ostium clearance is off for the proximal bound here, which is the whole
        point: leaving it on, every station asked for is pushed straight back out to
        where it already was. It is not abandoned, though. If no rung cuts, the last
        act is to restore the cleared placement, so a case whose ascending aorta has
        no planar cross-section near the trunk keeps the behaviour that at least gets
        it cut -- KK28_11 and KP20_1 both refuse every station from 5 to 30 mm.
        """
        options = list(p.get("s0_options") or [])
        if not options or self.region is None:
            return
        if p["s0"] not in options:
            options.insert(0, p["s0"])
        for i, s0 in enumerate(options):
            self.set_region(s0, s1, p["name"], manual=False, clear_start=False)
            try:
                st = self.recut()
            except CutError:
                continue
            if st["ok"]:
                self._preset_start_mm = (
                    None if i == 0 else round(float(s0) - float(p["s0"]), 2))
                return
        self.set_region(p["s0"], s1, p["name"], manual=False)
        self.recut()
        self._preset_start_mm = None

    def clear_region(self):
        self._region_reach = None
        # Leaving a region entirely drops its boundary planes, as the window does:
        # carrying a tilt into the NEXT region would apply it to a plane somewhere
        # else. Moving a boundary WITHIN a region keeps them, which is the point.
        self._bound_planes = {}
        self._region_moved_mm = 0.0
        self._preset_extended_mm = None
        self._preset_start_mm = None
        self.region = None
        self.region_name = ""
        self._region_manual = False
        self._region_clear_start = True

    # ── the cut ──────────────────────────────────────────────────────────────
    def recut(self):
        """Re-apply every cut. Sets .mesh and .problem; returns the status dict.

        A REGION whose boundary planes will not cut gets one retry at a different
        plane reach. Only on failure, and only once: probing costs two extra cuts per
        attempt, and doing it up front made a region slider unusable -- the fitter ran
        on every nudge and the whole test suite stopped finishing. This is the same
        shape as suggest_pullback, which also only runs once a cut has refused.

        Measured over 182 region cuts: 10 of the 15 failures are rescued and the 6
        that severed a neighbouring vessel drop to none, with the 167 that already
        worked paying nothing at all.
        """
        st = self._recut_once()
        if not st["ok"]:
            st = self._retry_trunk_offset(st)
        if not st["ok"]:
            st = self._retry_bound_clearance(st)
        if st["ok"] or self.region is None or getattr(self, "_reach_tried", False):
            return st
        bounds, _inside, _notes = self._region_terminals()
        if any(getattr(t, "band_mm", None) for t in bounds):
            return st                         # a hand-set plane size is the answer
        self._reach_tried = True
        try:
            a, b = self.region
            got = fit_region(self.base, self.tree, self.route, a, b)
            if got is None:
                return st
            self._region_reach = got["scale"]
            if abs(got["moved_mm"]) > 1e-9:
                # The region MOVED. Recorded on the session, because a boundary that
                # is no longer where the caller asked for it changes what is being
                # measured, and a sidecar that does not say so is not reproducible.
                self._region_moved_mm = got["moved_mm"]
                self.region = (got["s0"], got["s1"])
            elif abs(got["scale"] - 1.0) < 1e-9:
                return st                     # nothing to change
            self._cut_steps = None            # the planes changed; the cache has not
            return self._recut_once()
        finally:
            self._reach_tried = False

    def _double_opening_cut(self, active, offs):
        """The one cut that leaves two openings on its own, or None.

        Costs one cut per opening and runs only after a re-cut has already failed on
        an opening count, so an eleven-opening model pays about two seconds for an
        answer that would otherwise be a guess. Region bounds are excluded: they have
        their own remedy in _retry_bound_clearance, which moves them along the aorta
        rather than trimming a vessel end.
        """
        attr = vessel_attribution(self.base, self.tree)
        snap = snap_tolerance(self.base) if self.clean_rims else 0.0
        for t in active:
            if getattr(t, "_region_offset", None) is not None:
                continue
            try:
                out = cut_terminal(self.base, self.tree, t, float(offs[t.label]),
                                   snap_mm=snap, _attr=attr, _suggest=False)
            except Exception:                                   # noqa: BLE001
                continue
            if len(boundary_loops(out)) > 1:
                return t
        return None

    def _retry_bound_clearance(self, st):
        """A region boundary undercutting an ostium gets nudged clear of it.

        See clear_bound_offset for the evidence. Only fires when the cut came back
        with MORE openings than were asked for, which is what an undercut ostium
        looks like, and only on a boundary -- every other opening sits where the
        anatomy puts it.

        The move is recorded and reported. A boundary that is not where the caller
        put it changes which stretch of aorta was measured, the same reason
        _region_moved_mm is carried.

        A BOUNDARY THE USER DRAGGED IS THEIRS. region_plan already refuses to clear
        ostia once _region_manual is set -- "a dragged boundary is the user's, do not
        nudge" -- and this retry has to obey the same rule or it undoes it: moving a
        hand-placed bound up to 8 mm to make it cut is exactly the slider-snapping
        this repository fixed once already. Reported from the tool as the start
        jumping when it is moved a millimetre.

        So it fires only on a boundary a PRESET placed. If a station someone chose by
        hand will not cut, the honest answer is the refusal and the reason.
        """
        if (self.tree is None or self.base is None or self.region is None
                or not self.problem
                or self._region_manual
                or getattr(self, "_bound_clear_tried", False)):
            return st
        active, offs = self._cut_list()
        bounds = [t for t in active if getattr(t, "_region_offset", None) is not None]
        # TWO failures answer to the same remedy, and both are the boundary sitting
        # somewhere it cannot make a clean opening.
        #
        #   one opening too many -- it is undercutting a branch's ostium, and the
        #     hole left at that branch's base is the extra loop;
        #   the boundary itself refusing, or coming back not flat -- there is no
        #     clean cross-section where it was put.
        #
        # Too FEW openings is not one of them: that is a region collapsing, and
        # moving a plane will not conjure the missing cut.
        count = re.search(r"(\d+) openings, expected (\d+)", str(self.problem))
        if count:
            if int(count.group(1)) <= int(count.group(2)):
                return st
        elif not any(str(self.problem).startswith(str(t.name)) for t in bounds):
            return st
        if not bounds:
            return st
        attr = vessel_attribution(self.base, self.tree)
        snap = snap_tolerance(self.base) if self.clean_rims else 0.0
        moved = []
        for t in bounds:
            was = float(offs[t.label])
            found = clear_bound_offset(self.base, self.tree, t, was,
                                       snap_mm=snap, attr=attr)
            if found is not None and abs(found - was) > 1e-9:
                moved.append((str(t.name), was, found))
        if not moved:
            return st

        self._bound_clear_tried = True
        keep = dict(getattr(self, "_bound_nudge", None) or {})
        try:
            self._bound_nudge = dict(keep)
            self._bound_nudge.update({n: off for n, _w, off in moved})
            self._cut_steps = None        # the planes moved; the cache has not
            out = self._recut_once()
            if out["ok"]:
                self._bound_moved = [{"name": n, "from_mm": round(w, 2),
                                      "to_mm": round(f, 2)}
                                     for n, w, f in moved]
                # status() was built before the move was recorded, so ask again --
                # otherwise the mesh is right and the line describing it is silent.
                return self.status()
            self._bound_nudge = keep      # no better: put them back as planned
            self._cut_steps = None
            return self._recut_once()
        finally:
            self._bound_clear_tried = False

    def _retry_trunk_offset(self, st):
        """A vessel end that will not cut cleanly gets one retry at an offset that does.

        See flat_trunk_offset for why the offset at a free end may be moved when an
        interior plane's may not: it only says how much of the end stub to trim, and a
        shallower cut keeps MORE vessel.

        Fires on either shape of failure -- the terminal named in the message, or an
        opening count nobody is named in, where the offending cut is FOUND rather than
        guessed at (see _double_opening_cut).

        Recorded on the session, and reported, because a cut that is not where the
        caller asked for it changes what is being measured -- the same reason
        _region_moved_mm exists.
        """
        if self.tree is None or self.base is None or not self.problem:
            return st
        active, offs = self._cut_list()
        ends = [x for x in active if is_model_end(self.tree, x)]
        # NAMED in the message, or blamed by an opening count nobody is named in.
        #
        # A model end whose plane runs through a neighbouring vessel opens that vessel
        # too, and the failure then reads "N openings, expected N-1" -- which names no
        # terminal at all, so a rule keyed on the name never fired. JK14_2's inlet
        # plane passes through the celiac trunk, whose faces sit a median 1.01 mm off
        # it; the cut makes two openings and the case refuses on the plain path.
        # Same fault as a rim that will not come out flat, same remedy.
        t = next((x for x in ends
                  if str(self.problem).startswith(str(x.name))), None)
        count = re.search(r"(\d+) openings, expected (\d+)", str(self.problem))
        if t is None and count and int(count.group(1)) > int(count.group(2)):
            # NOBODY IS NAMED, so find the cut that is at fault instead of guessing.
            #
            # A plane running through a neighbouring vessel opens that vessel too, and
            # the failure reads "N openings, expected N-1", naming no terminal. Which
            # cut did it is not a guess: cut each one alone and see which produces two
            # openings by itself. On JK14_2 that is the inlet, whose plane crosses the
            # celiac trunk; on JK30_4 it is Branch 4, so keying this to model ends
            # only -- as the first version did -- missed it.
            t = self._double_opening_cut(active, offs)
        if t is None or getattr(self, "_trunk_backoff_tried", False):
            return st

        want = float(offs[t.label])
        found = flat_trunk_offset(
            self.base, self.tree, t, want,
            snap_mm=snap_tolerance(self.base) if self.clean_rims else 0.0)
        if found is None:
            return st

        self._trunk_backoff_tried = True
        try:
            self.set_offset(t.label, found)
            self._cut_steps = None            # the plane moved; the cache has not
            got = self._recut_once()
            if got["ok"]:
                self._trunk_moved = {"label": t.label, "name": str(t.name),
                                     "from_mm": round(want, 2),
                                     "to_mm": round(found, 2)}
                return self.status()      # rebuilt, so the note is in the message
            self.set_offset(t.label, want)    # no better: put it back as asked
            self._cut_steps = None
            return self._recut_once()
        finally:
            self._trunk_backoff_tried = False

    def _recut_once(self):
        active, offs = self._cut_list()

        # A region can collapse to nothing once its branch ostia are cleared. With
        # no cuts to make every check below passes on the untouched watertight mesh,
        # and an export would write the CLOSED original under an _open name.
        if not active:
            self.mesh = self.base.copy()
            self.matched = {}
            self.problem = (
                "Nothing is being cut. This region collapses to nothing once its "
                "branch ostia are cleared. Widen it, or go back to the whole vessel."
                if self.region is not None else
                "Nothing is being cut. No openings are selected.")
            return self.status()

        pre = preflight(self.tree, offs, terminals=active)
        if pre is not None:
            self.mesh = self.base.copy()
            self.problem = pre
            return self.status()

        snap = snap_tolerance(self.base) if self.clean_rims else 0.0

        # RESUME rather than redo, and cut the opening that moved LAST.
        #
        # Re-applying every cut from the base mesh means a five-opening model pays for
        # five cuts to change one, which over a network is the difference between a
        # slider that answers and one that does not. Each cut depends only on the cuts
        # before it, so the mesh after each is kept and a re-cut restarts at the first
        # opening whose offset actually moved.
        #
        # The order is free to change: cutting these in any sequence gives a
        # byte-identical mesh -- checked over eight permutations each on a 5-opening
        # and an 11-opening model -- because each cut removes a different vessel. So a
        # changed opening moves to the end and the next change to it resumes from
        # everything before. Dropped whole if any cut threw, since replaying from a
        # step that failed would silently skip it.
        # EVERY setting that changes the plane, not just the offset. The resume cache
        # existed to skip re-cutting openings whose OFFSET had not moved, and it keyed
        # on the offset alone -- so switching a cut from centerline-normal to Z-plane,
        # tilting it, or resizing its band left the key unchanged and the cached mesh
        # was reused. The plane visibly moved and the mesh did not, and the only way
        # out was to nudge the slider and put it back, which is what people were
        # doing. Reported from the tool.
        def _plane_key(t):
            return (round(float(offs[t.label]), 6), round(snap, 6), str(t.mode),
                    round(float(getattr(t, "tilt_a_deg", 0.0) or 0.0), 4),
                    round(float(getattr(t, "tilt_b_deg", 0.0) or 0.0), 4),
                    round(float(getattr(t, "shift_a_mm", 0.0) or 0.0), 4),
                    round(float(getattr(t, "shift_b_mm", 0.0) or 0.0), 4),
                    round(float(getattr(t, "band_mm", 0.0) or 0.0), 4))

        want = {t.label: _plane_key(t) for t in active}
        by_label = {t.label: t for t in active}
        cache = getattr(self, "_cut_steps", None)
        reuse = 0
        if (cache is not None and cache["base_id"] == id(self.base)
                and set(cache["order"]) == set(want)):
            for lbl in cache["order"]:
                if want.get(lbl) == cache["steps"].get(lbl):
                    reuse += 1
                else:
                    break
            rest = cache["order"][reuse:]
            order = (cache["order"][:reuse]
                     + [l for l in rest if want[l] == cache["steps"].get(l)]
                     + [l for l in rest if want[l] != cache["steps"].get(l)])
        else:
            order = [t.label for t in active]
        active = [by_label[l] for l in order]

        m = cache["meshes"][reuse].copy() if reuse else self.base.copy()
        self.problem = None
        meshes = (cache["meshes"][:reuse] if reuse else []) + [m]   # [k] = after k cuts
        for t in active[reuse:]:
            try:
                m = cut_terminal(m, self.tree, t, offs[t.label], snap_mm=snap)
            except Exception as e:
                self.problem = f"{t.name}: {e}"
            meshes.append(m)
        self._cut_steps = ({"base_id": id(self.base), "order": order,
                            "steps": want, "meshes": meshes}
                           if self.problem is None else None)
        self.mesh = keep_largest_component(m)
        self.matched = match_openings_to_terminals(self.mesh, self.tree, offs,
                                                   terminals=active)
        if self.problem is None:
            # Topology alone cannot see a cut that escapes into the parent vessel:
            # it still removes one disk, so the euler number stays right. Only the
            # opening's position catches it.
            self.problem = validate_cut(self.mesh, self.tree, offs, terminals=active)
        return self.status()

    # ── reading the result ───────────────────────────────────────────────────
    @property
    def ok(self):
        return self.problem is None

    def status(self):
        if self.mesh is None:
            return {"ok": False, "message": "nothing cut yet"}
        loops = len(boundary_loops(self.mesh))
        if self.problem is not None:
            return {"ok": False, "message": self.problem, "n_openings": loops}
        slivers = rim_sliver_count(self.mesh)
        # A cut the tool MOVED is said out loud. The offset at a model end is
        # arbitrary and backing it off is how four inlets come out flat at all, but
        # the result is not the mesh that was asked for and a status line that does
        # not say so is not reproducible.
        moved = getattr(self, "_trunk_moved", None)
        note = ""
        if moved:
            note = (f". {moved['name']} trimmed at {moved['to_mm']:g} mm instead of "
                    f"{moved['from_mm']:g} mm, which is where it comes out flat")
        bmoved = getattr(self, "_bound_moved", None)
        if bmoved:
            note += ". " + "; ".join(
                f"{b['name']} moved {abs(b['to_mm'] - b['from_mm']):g} mm to clear a "
                f"branch it was undercutting" for b in bmoved)
        return {
            "ok": True,
            # Wording matches branch_gui exactly. The two status lines say the same
            # thing about the same mesh, so they should not read differently.
            "message": (f"Ready, {loops} openings, "
                        f"euler {int(self.mesh.euler_number)}"
                        f"{f', {slivers} rim sliver(s)' if slivers else ', no slivers'}"
                        f"{note}"),
            "trunk_offset_moved": moved,
            "bounds_moved": bmoved,
            "n_openings": loops,
            "euler_number": int(self.mesh.euler_number),
            "n_faces": int(len(self.mesh.faces)),
            "area_mm2": round(float(self.mesh.area), 2),
            "rim_slivers": slivers,
        }

    def openings(self):
        """One row per opening being cut: name, offset, and measured diameter."""
        active, offs = self._cut_list()
        real = set(id(t) for t in self.tree.terminals)
        rows = []
        for t in active:
            m = self.matched.get(t.label) or {}
            row = {
                "label": t.label,
                "name": t.name,
                "anatomical_name": self.anat_names.get(t.label),
                "is_trunk": bool(getattr(t, "is_trunk", False)),
                "cut_mode": getattr(t, "mode", None),
                "offset_mm": offs[t.label],
                "limit_mm": round(self.offset_limit(t), 1),
                "diameter_mm": (round(float(m["diameter"]), 2)
                                if m.get("diameter") is not None else None),
                # A trim is not in tree.terminals either, so "not one of the tree's
                # own" is not the same question as "is a region bound" -- answering
                # it that way disabled the slider on every trimmed opening.
                "is_region_bound": (id(t) not in real
                                    and not str(t.label).endswith(".trim")),
                "is_trim": str(t.label).endswith(".trim"),
                "trim_segment": (t.segment.name
                                 if str(t.label).endswith(".trim") else None),
            }
            row.update(self.plane_controls(t.label))
            rows.append(row)
        return rows

    def offset_limit(self, terminal):
        """How far back this cut's slider may travel from its tip. This is the window's
        _row_limit, and it must match it exactly.

        An interior (trim) cut is hemmed in at BOTH ends: run it too far back and
        the plane arrives at the junction feeding the vessel, where it tears. So
        stop a couple of radii short of that rather than at the segment's end.

        There is deliberately no absolute ceiling. An earlier version here capped
        the travel at 120 mm, which silently shortened every vessel longer than
        ~126 mm, the left subclavian among them, so the web tool could not reach
        depths the desktop tool offers.
        """
        L = float(getattr(terminal.segment, "length", 0.0) or 0.0)
        if L <= 0:
            return 60.0
        if getattr(terminal, "is_interior", False):
            r = max(float(getattr(terminal.segment, "_radius_mm", None) or 5.0), 1.0)
            return float(np.clip(L - 2.0 * r, MIN_OFFSET_MM + 1.0,
                                 MAX_FRAC_OF_SEGMENT * L))
        return float(MAX_FRAC_OF_SEGMENT * L)

    # ── measurement ──────────────────────────────────────────────────────────
    def measure(self, step_mm=2.0):
        """Equivalent-area diameter along the aorta, station by station.

        EAD is 2*sqrt(A/pi): the diameter of the circle with the same
        cross-sectional area. A vessel section is rarely a circle, so there is no
        single honest "diameter" -- EAD is the standard vascular measure because
        area is what governs flow, and it does not depend on which way you happen
        to lay a caliper across an oval.

        Restricted to the isolated region when there is one, because that is what
        the person is looking at.
        """
        if self.route is None or self.base is None:
            return None
        if self.region is None:
            s0, s1 = 0.0, float(self.route.length)
        else:
            _a, _b = self.region
            # An open region runs to the end of the aortic route, so measure to
            # there rather than refusing.
            s0 = float(_a) if _b is None else min(_a, _b)
            s1 = float(self.route.length) if _b is None else max(_a, _b)
        step = max(float(step_mm), 0.25)

        prof = []
        for s in np.arange(s0, s1 + 1e-9, step):
            try:
                nm, loc, _ = self.route.locate(float(s))
                pos, tan = self.tree.segments[nm].at_arclength(loc)
                sd = section_diameter(self.base, pos, tan, near_point=pos)
            except Exception:
                continue
            if sd and sd.get("ead_mm"):
                prof.append((round(float(s), 2), round(float(sd["ead_mm"]), 2)))
        if not prof:
            return None

        ds = [d for _, d in prof]
        i_max = int(np.argmax(ds))
        i_min = int(np.argmin(ds))
        return {
            "n_stations": len(prof),
            "step_mm": step,
            "length_mm": round(s1 - s0, 1),
            "max_mm": ds[i_max],
            "max_at_mm": prof[i_max][0],
            "min_mm": ds[i_min],
            "min_at_mm": prof[i_min][0],
            "mean_mm": round(float(np.mean(ds)), 2),
            "start_mm": prof[0][1],
            "end_mm": prof[-1][1],
            "region": self.region_name if self.region is not None else None,
            "profile": prof,
        }

    # ── colouring the surface by which branch owns it ────────────────────────
    #
    # Lifted out of branch_gui so the two front-ends cannot drift: "which branch is
    # this face" ought to have one answer, and having the window and the web tool
    # each carry their own copy of the rule is how they stop agreeing.

    MESH_RGB = (201, 66, 74)          # the aorta, and anything unclaimed
    PALETTE = ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#d55e00",
               "#56b4e9", "#8a6d3b", "#7b2fbf", "#2ca02c", "#c0392b"]
    # A face is tinted as a branch only if it lies within this many of that branch's
    # own radii of its centerline. A real branch-wall face sits about one radius
    # out; a patch of fat aorta wall that the nearest-point rule mis-captured sits
    # many branch radii out, and stays red.
    TINT_RADII = 2.0

    def opening_colours(self):
        """label -> hex colour, branches only, in row order. Trunks stay red."""
        out, i = {}, 0
        for t in self._cut_list()[0]:
            if not getattr(t, "is_trunk", False):
                out[t.label] = self.PALETTE[i % len(self.PALETTE)]
                i += 1
        return out

    def face_colours(self, mesh=None):
        """Per-face RGB for the cut surface, or None if there is nothing to tint.

        Two gates, and both are load-bearing. RADIAL rejects the fat-aorta-wall
        patches the nearest-centerline rule mis-captures at an ostium. CONTINUITY
        keeps only the one surface-connected piece grown from the branch's own tip,
        so a scrap that passes the radial test while sitting over on another branch
        is not painted this branch's colour.
        """
        shown = self.mesh if mesh is None else mesh
        if shown is None or self.tree is None:
            return None
        colours = self.opening_colours()
        if not colours:
            return None
        owner, _s, radial, _r = vessel_attribution(shown, self.tree)
        centers = np.asarray(shown.triangles_center)

        by_seg = {}
        for t in self._cut_list()[0]:
            if t.label in colours:
                by_seg.setdefault(t.segment.name, []).append(t)

        out = np.tile(np.array(self.MESH_RGB, np.uint8), (len(owner), 1))
        for seg, terms in by_seg.items():
            if len(terms) != 1:
                continue                  # shared segment: ambiguous, leave it red
            t = terms[0]
            R = max(float(getattr(t.segment, "_radius_mm", None) or 5.0), 1.0)
            cand = (owner == seg) & (radial <= self.TINT_RADII * R)
            if not cand.any():
                continue
            ci = np.flatnonzero(cand)
            tip_face = int(ci[np.argmin(
                np.linalg.norm(centers[ci] - np.asarray(t.tip, float), axis=1))])
            keep = _component(shown, cand, [tip_face])
            h = colours[t.label].lstrip("#")
            out[keep] = [int(h[i:i + 2], 16) for i in (0, 2, 4)]
        return out

    def measure_at(self, world_point):
        """How far a point on the surface sits from the node its vessel comes off.

        The desktop shows this under the cursor; here it is a click, because every
        hover would be a round trip. Same underlying answer either way.
        """
        if self.tree is None:
            return None
        res = measure_along_branch(self.tree, np.asarray(world_point, float))
        if not res:
            return None
        # Name it the way the opening list does. anat_vessels only covers the named
        # trunk segments, so a branch came back as "Line6" -- a label nobody reading
        # the screen has any use for.
        seg = res["segment"]
        named = next((t.name for t in self._cut_list()[0] if t.segment.name == seg),
                     None)
        return {
            "mm_to_node": round(float(res["distance_mm"]), 1),
            "vessel": named or self.anat_vessels.get(seg, seg),
            "segment": res["segment"],
            "s_mm": round(float(res["s_mm"]), 1),
            "at_end": bool(res.get("at_end")),
            "from_vessel_end": bool(res.get("from_vessel_end")),
            "point": [float(v) for v in res["proj"]],
        }

    def offset_for_mark(self, mark):
        """The cut depth that would put a cut AT a pinned mark, and on which opening.

        A mark is a station on a vessel; a cut is a depth measured from that vessel's
        free end. Turning one into the other is the whole point of pinning it -- you
        measure where the disease stops and then cut there, rather than reading a
        number off the screen and typing it into a slider.
        """
        seg = str(mark.get("segment") or "")
        t = next((x for x in self._cut_list()[0]
                  if x.segment.name == seg and not getattr(x, "is_region_bound", False)),
                 None)
        if t is None:
            return None
        L = float(t.segment.length)
        s = float(mark.get("s_mm", 0.0))
        off = (L - s) if t.at_end else s
        off = float(np.clip(off, MIN_OFFSET_MM, self.offset_limit(t)))
        return {"label": t.label, "name": t.name, "offset_mm": round(off, 1)}

    def measure_station(self, s_mm):
        """Where one station along the measured stretch is, and how wide.

        The summary says "widest 41.2 mm at 86 mm" and gives no way to see whether
        86 mm is where you meant. This is that station as a plane -- origin, normal
        and a radius to draw it at -- so the placement can be looked at and moved,
        which is the one thing the numbers cannot show.

        Returns None outside the measurable stretch, or where the mesh has no
        cross-section to take.
        """
        if self.route is None or self.base is None:
            return None
        lo, hi = self.measured_range()
        s = float(np.clip(float(s_mm), lo, hi))
        try:
            nm, loc, _ = self.route.locate(s)
            pos, tan = self.tree.segments[nm].at_arclength(loc)
        except Exception:                                       # noqa: BLE001
            return None
        sd = section_diameter(self.base, pos, tan, near_point=pos)
        dia = float(sd["ead_mm"]) if sd and sd.get("ead_mm") else 0.0
        return {
            "s_mm": round(s, 1),
            "origin": [float(v) for v in pos],
            "normal": [float(v) for v in tan],
            # A little wider than the lumen, so the disc's edge reads against the
            # wall rather than hiding inside it.
            "radius": round(max(0.8 * dia, 2.0), 2),
            "ead_mm": round(dia, 2) if dia else None,
            "min_mm": round(float(sd["min_mm"]), 2) if sd and sd.get("min_mm") else None,
            "max_mm": round(float(sd["max_mm"]), 2) if sd and sd.get("max_mm") else None,
        }

    def measured_range(self):
        """(start, end) in mm of what `measure` covers: the region, or the vessel."""
        if self.route is None:
            return (0.0, 0.0)
        if self.region is None:
            return (0.0, float(self.route.length))
        a, b = self.region
        if b is None:
            return (float(a), float(self.route.length))
        return (float(min(a, b)), float(max(a, b)))

    def angulation(self, step_mm=1.0):
        """How much the centerline turns over the isolated region.

        Fits the circle that best matches the centerline there and reports the arc
        it subtends: theta = L / R. For a seal zone that is neck angulation, and it
        is reproducible in a way a hand-drawn pair of tangents is not, because it
        does not depend on where the tangents are drawn.

        Whole vessel when nothing is isolated, which is rarely meaningful: an
        aorta is not one arc. Isolate the zone first.
        """
        if self.route is None:
            return None
        if self.region is None:
            s0, s1 = 0.0, float(self.route.length)
            whole = True
        else:
            _a, _b = self.region
            s0 = float(_a) if _b is None else min(_a, _b)
            s1 = float(self.route.length) if _b is None else max(_a, _b)
            whole = False
        res = centerline_arc_angle(self.tree, self.route, s0, s1, step_mm=step_mm)
        if res is None:
            return None
        res["region"] = None if whole else (self.region_name or "custom")
        res["whole_vessel"] = whole
        # A jagged centerline inflates arc length, and theta = L/R inflates with it.
        res["smoothed_mm"] = self.smooth_mm
        res["centerline_was_jagged"] = bool(self.jagged)
        return res

    @property
    def is_thoracic(self):
        """The anatomy call for this case -- one trunk end is an abdominal aorta, more
        is a thoracic one. Cheap, so a front end can ask it on every redraw. False
        until something is loaded."""
        if self.tree is None:
            return False
        return sum(1 for t in self.tree.terminals if getattr(t, "is_trunk", False)) != 1

    def anatomy_summary(self):
        kind = "thoracic" if self.is_thoracic else "abdominal"
        how = "you set this" if FAMILY[self.family] else "read from the geometry"
        # Whether any of this can be believed. A centerline running outside the mesh
        # makes the radii meaningless, and the family call, the names and the region
        # presets are all measured off them -- so the caller is told before it shows
        # any of them, not after the cut refuses.
        outside = list(getattr(self.tree, "radius_outside", None) or [])
        # A centerline that doubles back on itself describes no vessel, and every
        # length, radius and plane normal taken from it is measured off the wrong
        # thing. Reported here because a user handed one has no other way to know.
        try:
            folded = centerline_defects(self.tree)
        except Exception:                                       # noqa: BLE001
            folded = {}
        # How close the abdominal/thoracic call was. The desktop window has always
        # said so and this did not, so the web app asserted a marginal reading as a
        # fact. It matters more than it looks: measured on the centerlines the tool
        # actually cuts -- smoothed when jagged -- the two families OVERLAP. An
        # abdominal aorta reaches 0.225 and a thoracic one 0.249, so no cut-off can
        # separate them and 77 of 79 is this rule's ceiling rather than a tuning
        # failure. A marginal call is exactly when the user has to be told.
        split = getattr(self.tree, "outflow_split", None)
        close = (split is not None and not FAMILY[self.family]
                 and abs(float(split) - PAIR_AREA_FRACTION) < FAMILY_CLOSE_BAND)
        return {"kind": kind, "declared": bool(FAMILY[self.family]), "how": how,
                "close_call": bool(close),
                "outflow_split": (None if split is None else round(float(split), 3)),
                "family_cutoff": PAIR_AREA_FRACTION,
                "named_openings": len(self.anat_names), "notes": list(self.anat_notes),
                "trusted": bool(getattr(self.tree, "family_trusted", True)),
                "mesh_matches": bool(getattr(self.tree, "mesh_matches", True)),
                "coverage": round(float(getattr(self.tree, "coverage", 1.0)), 3),
                "outside_mesh": outside,
                "folded_centerlines": {k: list(v) for k, v in folded.items()},
                "n_terminals": len(self.tree.terminals)}

    # ── output ───────────────────────────────────────────────────────────────
    def _region_geometry(self, s0, s1, active, offs):
        """Where the region sits on the aorta and how wide it is along the way.

        Written so a measurement can be reproduced, or compared across cases, without
        reloading the mesh: the two stations in aortic mm, the route they are a span
        of, the calibre at each end and every 2 mm between, which branches fell inside
        it, and anything region_plan wanted to say about the boundaries.

        THE STATIONS ARE WHERE THE PLANES ENDED UP, not what was asked for.

        These used to be read straight off self.region, which is the REQUEST. On the
        preset path region_plan then nudges a boundary clear of any ostium it
        undercuts, and the two are not close: measured over 30 arches the start ends
        up 1.9 to 32.1 mm proximal of the recorded station, median 7.5. DJ8_1 recorded
        a start of 73.65 while its own boundary note in the same file said the start
        had moved to 66.2, and a reader comparing cases on that field was comparing
        plans rather than cuts. A hand-placed region sets clear_ostia=False and is
        never moved, which is why this went unnoticed: 23 of the 30 agreed exactly.

        So take the stations from the BOUNDS. The request is kept alongside, under its
        own name, whenever the two differ -- that a boundary was not honoured is worth
        recording, not worth hiding.

        Best effort. A diameter is a section fitted to the surface and it can fail on
        a boundary that sits in a hairpin; a missing number is left None rather than
        allowed to stop an export that is otherwise sound.
        """
        req0, req1 = float(s0), (None if s1 is None else float(s1))
        got = []
        for t in active:
            if getattr(t, "_region_offset", None) is None:
                continue
            try:
                sc = float(t.s_cut(float(offs[t.label])))
                got.append(float(self.route.station_of(t.segment.name, sc)))
            except Exception:                                       # noqa: BLE001
                continue
        got.sort()
        # Two bounds is a closed region and one is an open distal end, where the
        # single plane is the START. Anything else -- a bound that could not be
        # located -- falls back to the request rather than inventing a station.
        if len(got) == 2 and s1 is not None:
            s0, s1 = got
        elif len(got) == 1 and s1 is None:
            s0 = got[0]

        end = float(self.route.length) if s1 is None else float(s1)
        out = {"start_mm_along_aorta": round(float(s0), 2),
               "end_mm_along_aorta": None if s1 is None else round(float(s1), 2),
               "length_mm": round(end - float(s0), 2),
               "aorta_route_length_mm": round(float(self.route.length), 2)}
        if abs(float(s0) - req0) > 0.01:
            out["start_requested_mm_along_aorta"] = round(req0, 2)
        if req1 is not None and s1 is not None and abs(float(s1) - req1) > 0.01:
            out["end_requested_mm_along_aorta"] = round(req1, 2)

        def calibre(s):
            try:
                nm, loc, _ = self.route.locate(float(s))
                pos, tan = self.tree.segments[nm].at_arclength(loc)
                sd = section_diameter(self.base, pos, tan, near_point=pos)
            except Exception:                                       # noqa: BLE001
                return None
            return round(float(sd["ead_mm"]), 2) if sd and sd.get("ead_mm") else None

        prof = []
        for s in np.arange(float(s0), end, 2.0):
            d = calibre(s)
            if d is not None:
                prof.append((round(float(s - float(s0)), 2), d))
        out["start_diameter_mm"] = calibre(s0)
        out["end_diameter_mm"] = None if s1 is None else calibre(s1)
        out["max_diameter_mm"] = max((d for _, d in prof), default=None)
        try:
            _bounds, _inside, notes = self._region_terminals()
        except Exception:                                           # noqa: BLE001
            notes = []
        # The branches as they were CUT, not as region_plan found them. A trim
        # collapses a group into one opening -- a brachiocephalic standing for both
        # right-sided vessels -- and listing the group's members here would name
        # vessels that have no opening in the mesh this file describes.
        out["branches_inside"] = [str(t.name) for t in active
                                  if getattr(t, "_region_offset", None) is None]
        out["boundary_notes"] = [str(n) for n in notes]
        out["diameter_profile_mm"] = prof
        return out

    def sidecar(self):
        """The JSON that travels with the STL. Mirrors branch_gui.CutWindow.export."""
        if self.mesh is None:
            raise CutError("nothing has been cut yet")
        n_loops = len(boundary_loops(self.mesh))
        area = float(self.mesh.area)
        lo, hi = CAMP_REF_AREA_MM2
        active, offs = self._cut_list()
        meta = {
            # File NAMES, not paths -- the name says which case this is, which is all
            # a reader needs; the path would additionally publish wherever it lived.
            "mesh": self.mesh_path.name,
            "centerline": self.csv_path.name,
            "euler_number": int(self.mesh.euler_number),
            "n_openings": n_loops,
            "n_faces": int(len(self.mesh.faces)),
            "area_mm2": round(area, 2),
            "snap_mm": round(snap_tolerance(self.base), 4) if self.clean_rims else 0.0,
            "rim_slivers": rim_sliver_count(self.mesh),
            "min_triangle_quality": round(min_triangle_quality(self.mesh), 5),
            "anatomy": {
                "frame": ANATOMY_FRAME,
                "left_right_verified": False,
                "named_openings": len(self.anat_names),
                # THE SAME KEYS THE WINDOW WRITES -- this file has the same name and
                # downstream readers cannot tell which front end produced it. Both
                # were missing here, so a web or headless export recorded only a
                # COUNT of named openings: not which opening got which name, and not
                # which of them a person declared.
                #
                # That second one is not bookkeeping. Every side name inferred from
                # the geometry rests on the frame recorded just above, and a mirrored
                # export would swap all of them -- but a HAND-DECLARED name does not
                # rest on it and is not subject to that doubt. A reader checking the
                # left/right convention, or placing a measurement against a
                # hand-named renal, has to be able to tell the two apart.
                "named_by_hand": dict(getattr(self.tree, "manual_names", None) or {}),
                "names": dict(self.anat_names),
                "notes": list(self.anat_notes),
            },
            "camp": {
                "target_faces": CAMP_TARGET_FACES,
                "pre_remesh_euler": int(self.mesh.euler_number),
                "pre_remesh_loops": n_loops,
                "reference_area_mm2": list(CAMP_REF_AREA_MM2),
                "area_in_reference_range": bool(lo <= area <= hi),
            },
            "openings": [],
        }
        if self.region is not None:
            _a, _b = self.region
            s0, s1 = ((float(_a), None) if _b is None
                      else (min(_a, _b), max(_a, _b)))
            # s0_mm / s1_mm are the stations to ASK FOR to reproduce this cut, which
            # is what nub_session restores a region from. They are deliberately the
            # request: feeding back a nudged station would nudge it again. Where the
            # planes actually landed is start_mm_along_aorta / end_mm_along_aorta,
            # added below -- that pair is the one to measure from.
            meta["region"] = {
                "name": self.region_name, "s0_mm": round(s0, 2),
                # None records an open distal end: one plane, everything below kept.
                "s1_mm": None if s1 is None else round(s1, 2),
                "open_distal": s1 is None,
                "manual": self._region_manual}
            # THE SAME KEYS THE WINDOW WRITES, because this file has the same name.
            #
            # Both writers produce <case>_cut_planes.json and they did not agree: the
            # window recorded start_mm_along_aorta and a diameter profile, the API
            # recorded s0_mm and nothing else. Everything downstream was written
            # against the window's shape, so a case re-cut in batch came out looking
            # like a case with no region at all -- arch_viewer printed no span in its
            # header and reported both boundaries as "not found", silently, on an
            # export that was perfectly good.
            #
            # Additive: s0_mm and the rest stay, since nub_session and the tests
            # restore a region from them.
            meta["region"].update(self._region_geometry(s0, s1, active, offs))
        # WHAT THE TOOL MOVED, and what it was asked for.
        #
        # The offsets below are the ones actually cut, so the mesh reproduces from
        # this file either way. What they cannot say is that a plane is not where the
        # caller put it: a vessel end trimmed shallower because that is where it comes
        # out flat, or a region boundary nudged clear of an ostium it was undercutting.
        # Both change where a measurement was taken, and someone reading this file
        # months later has no other way to know it happened.
        adjust = {}
        if getattr(self, "_trunk_moved", None):
            adjust["vessel_end"] = self._trunk_moved
        if getattr(self, "_bound_moved", None):
            adjust["region_bounds"] = self._bound_moved
        if float(getattr(self, "_region_moved_mm", 0.0) or 0.0):
            adjust["region_shifted_mm"] = round(
                float(self._region_moved_mm), 2)
        if adjust:
            meta["adjustments"] = adjust
        for t in active:
            m = self.matched.get(t.label) or {}
            # The plane as a plane, not only as a depth. The viewer draws a rim from
            # its centroid and normal, and an export missing them renders as a mesh
            # with no openings marked on it -- so a folder half from the window and
            # half from here reviewed inconsistently, on the same cases.
            point, outward = t.cut_plane(float(offs[t.label]))
            meta["openings"].append({
                "name": t.name,
                "spline": t.label,
                "anatomical_name": self.anat_names.get(t.label),
                "cut_mode": getattr(t, "mode", None),
                # A trim is one opening standing in for the group it removed.
                "trimmed_group": (sorted(t.subtree_segments)
                                  if getattr(t, "is_interior", False) else None),
                "offset_mm": offs[t.label],
                # The plane as it was actually set, so a cut can be reproduced
                # rather than approximated from its depth alone.
                "tilt_a_deg": round(float(t.tilt_a_deg), 2),
                "tilt_b_deg": round(float(t.tilt_b_deg), 2),
                "shift_a_mm": round(float(t.shift_a_mm), 2),
                "shift_b_mm": round(float(t.shift_b_mm), 2),
                "plane_size_mm": round(float(band_width(t)), 2),
                "plane_origin": [round(float(v), 4) for v in point],
                "plane_normal": [round(float(v), 6) for v in outward],
                "diameter_mm": (round(float(m["diameter"]), 2)
                                if m.get("diameter") is not None else None),
                "centroid": ([round(float(c), 4) for c in m["centroid"]]
                             if m.get("centroid") is not None else None),
                # HOW FLAT THIS OPENING ACTUALLY CAME OUT, in mm off its own plane.
                #
                # The cutter measures this and used to discard it unless it crossed
                # the refusal bar at 3 mm, so an opening anywhere below that was
                # exported looking identical to a perfect one. It is not identical:
                # between RIM_FLAW_MM and TORN_REFUSE_MM the rim is visibly torn and
                # the cut has taken a piece of a neighbouring vessel with it -- 31
                # mm2 of one on JK24_1's Branch 4, which is inside a cohort now.
                #
                # It belongs in the record because it is a property of the SURFACE a
                # downstream measurement is taken on. delta K is an area-weighted
                # variance of curvature; a torn rim is real area with real curvature
                # that the anatomy does not have. Someone comparing two numbers
                # months from now needs to be able to ask which meshes were clean,
                # and until this key existed there was nothing in the file to ask.
                #
                # 0.0 means measured and flat. null means the cheap test said the
                # boundary came entirely from clipping, so it is on the plane by
                # construction and was never measured -- a different statement from
                # "measured and flat", and worth keeping distinct.
                "rim_dev_mm": (round(float(t._rim_dev_mm), 3)
                               if getattr(t, "_rim_dev_mm", None) is not None
                               else None),
                "rim_flawed": bool(
                    getattr(t, "_rim_dev_mm", None) is not None
                    and float(t._rim_dev_mm) > RIM_FLAW_MM),
            })
        return meta

    def export(self, out_dir, name):
        """Write <name>.stl and <name>_cut_planes.json. Refuses an invalid cut."""
        if not self.ok:
            raise CutError(self.problem)
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        stl = out_dir / f"{name}.stl"
        side = out_dir / f"{name}_cut_planes.json"
        self.mesh.export(str(stl))
        side.write_text(json.dumps(self.sidecar(), indent=2))
        return stl, side
