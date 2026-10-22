# Interactive cut tool for vascular meshes (branched or unbranched)

# Load a mesh on its own to look at it, or a mesh + centerline to cut it. The
# tree, the trunk/branch split and the cut modes are all worked out from the
# geometry. Each opening gets a slider AND a mm box; clicking one opens its
# advanced controls (tilt, plane size), which stay hidden otherwise.

#   python branch_gui.py                      -> opens empty, pick files in the UI
#   python branch_gui.py <centerline.csv> <mesh.mat>
#   python branch_gui.py -- <mesh.mat>        -> mesh only, no cutting

import json
import re
import sys
from pathlib import Path

import numpy as np
import trimesh
from pyvistaqt import QtInteractor
from qtpy import QtCore, QtWidgets

import pyvista as pv

from branch_cut import (
    _component,
    band_width,
    boundary_loops,
    centerline_arc_angle,
    centerline_geodesic,
    classify_terminals,
    cut_terminal,
    keep_largest_component,
    load_mat_mesh,
    match_openings_to_terminals,
    measure_along_branch,
    anatomy_choices,
    aorta_route,
    min_triangle_quality,
    name_anatomy,
    set_manual_name,
    pre_bifurcation_candidates,
    region_plan,
    fit_region_reach,
    fit_region,
    region_presets,
    preflight,
    rim_sliver_count,
    section_diameter,
    segment_radius,
    snap_tolerance,
    suggest_pullback,
    flat_trunk_offset,
    clear_bound_offset,
    is_model_end,
    validate_cut,
    vessel_attribution,
    vessel_coords,
    ANATOMY_FRAME,
    MIN_CENTERLINE_COVERAGE,
    MIN_OFFSET_MM,
    PAIR_AREA_FRACTION,
    SIDES_ASSUMED_NOTE,
)
from branch_tree import load_tree, Terminal

DEFAULT_OFFSET_MM = 6.0

# CAMP consumes these cut meshes to compute delta K (k_fluct, the area-weighted
# variance of per-patch integrated Gaussian curvature). That number is EXTENSIVE
# and partition-dependent: cut extent, face count, trim fraction and topology all
# move it. So a case is only comparable against the reference cohort when all four
# match, and the sidecar records ours rather than leaving them to be re-derived.
#
# The branched thoracic reference ran at a fixed face count over a bounded area
# range. A surface outside that range cannot be placed on the master curve, the
# fit would be extrapolated past its data, not interpolated within it.
CAMP_TARGET_FACES = 100_000
CAMP_REF_AREA_MM2 = (26_252.0, 98_459.0)

# How far back a cut may travel, as a fraction of its vessel's length. The honest
# limit is per-vessel (81%-98% on the sample), so no single fraction is both
# generous and safe. Since a bad cut is caught reliably (validate_cut checks
# geometry as well as topology, and export is blocked), the slider opens up and
# the status bar says when a cut has gone too far.
MAX_FRAC_OF_SEGMENT = 0.95

MAX_TILT_DEG = 60.0

HOVER_NODE = "#7b2fbf"      # the branch node the distance is measured to
HOVER_PT = "#111111"        # where the cursor's pick lands on the centerline
HOVER_PATH = "#7b2fbf"      # the centerline stretch being measured
MESH_COLOR = "#c9424a"
TRUNK_LOOP = "#00b050"
BRANCH_LOOP = "#0070c0"
CL_COLOR = "#1b3b6f"
PLANE_COLOR = "#ffb000"
PLANE_SEL = "#ff5c00"

MARK_COLOR = "#0a7d2c"      # pinned distance marks
MARK_SPAN = "#0a7d2c"       # the along-centerline span between two marks

# Distinct, reasonably colour-blind-friendly hues for the "colour openings" option,
# assigned to openings in row order so a branch's surface, its plane disc, and its
# row label all read as the same colour.
OPENING_PALETTE = [
    "#0072b2", "#e69f00", "#009e73", "#cc79a7", "#d55e00",
    "#56b4e9", "#8a6d3b", "#7b2fbf", "#2ca02c", "#c0392b",
]


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) for i in (0, 2, 4)]


# A face is only tinted as a branch if it lies within this many of that branch's
# radii of the branch centerline. A real branch-wall face sits ~1 radius out; a
# fat-aorta-wall patch that the nearest-point rule mis-captured sits many branch
# radii out, so this keeps those patches red instead of painting them a branch.
RADIAL_TINT_FACTOR = 2.0


MODES = [("Z-plane (horizontal)", "z"), ("Centerline-normal", "normal")]


def to_pv(m):
    f = np.hstack([np.full((len(m.faces), 1), 3, np.int64), np.asarray(m.faces)]).ravel()
    return pv.PolyData(np.asarray(m.vertices), f)


def _pair(minv, maxv, val, step, suffix, decimals=1):
    """A slider + spin box bound to the same number."""
    s = QtWidgets.QSlider(QtCore.Qt.Horizontal)
    s.setRange(int(round(minv * 10)), int(round(maxv * 10)))
    s.setValue(int(round(val * 10)))
    b = QtWidgets.QDoubleSpinBox()
    b.setRange(minv, maxv)
    b.setDecimals(decimals)
    b.setSingleStep(step)
    b.setSuffix(suffix)
    b.setValue(val)
    b.setKeyboardTracking(False)
    return s, b


class TerminalRow(QtWidgets.QWidget):
    """One opening: name, cut mode, slider + mm box, and a collapsible panel of
    advanced plane controls (tilt, plane size).

    *bound* builds the row for one of an isolated REGION's two boundary cuts. Those
    are placed by the region's own Start/End sliders, so the offset control would be
    a second, conflicting way to move the same plane and is left out -- but the
    plane itself still has a tilt and a reach, and until this existed there was no
    way to touch either once a region was on. The panel starts open, because a row
    whose only content is hidden behind an arrow reads as an empty row.
    """

    preview = QtCore.Signal(object)   # cheap: move the plane (carries the row)
    changed = QtCore.Signal()         # commit: re-cut the mesh
    selected = QtCore.Signal(object)
    pullback = QtCore.Signal(object)  # "find me the deepest cut that works"
    renamed = QtCore.Signal(object)   # the artery was declared (carries the row)

    def __init__(self, terminal, limit_mm, auto_band_mm, parent=None, bound=False,
                 choices=()):
        super().__init__(parent)
        self.terminal = terminal
        self.limit = float(limit_mm)
        self.auto_band = float(auto_band_mm)
        self.bound = bool(bound)
        self.choices = tuple(choices or ())
        self._syncing = False

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(8, 6, 8, 6)
        outer.setSpacing(4)

        # ---- always-visible row ------------------------------------------------
        top = QtWidgets.QGridLayout()
        top.setHorizontalSpacing(8)
        top.setVerticalSpacing(3)

        self.expand = QtWidgets.QToolButton()
        self.expand.setCheckable(True)
        self.expand.setArrowType(QtCore.Qt.RightArrow)
        self.expand.setAutoRaise(True)
        self.expand.setToolTip("Tilt and plane size for this cut")
        self.expand.toggled.connect(self._toggle_advanced)

        # Editable rather than a plain box, because naming an opening has two jobs.
        # Free text labels it for your own reading; picking from the list DECLARES the
        # artery, and region_presets matches on the words in those names -- so a
        # hand-named renal is what brings back "Aortic neck" and "Infrarenal aorta"
        # on the 18 of 47 abdominal scans where the geometry cannot say which lateral
        # pair the renals are. Typing "left renal" would look identical and do nothing.
        self.name = QtWidgets.QComboBox()
        self.name.setEditable(True)
        self.name.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.name.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.name.addItem(terminal.name)
        for nm in (choices or ()):
            if nm != terminal.name:
                self.name.addItem(nm)
        self.name.setCurrentText(terminal.name)
        f = self.name.font()
        f.setBold(True)
        self.name.setFont(f)
        self.name.lineEdit().setFrame(False)
        self.name.setToolTip(
            "What this opening is. Type any label to rename it, or pick an artery\n"
            "from the list to declare it -- a declared artery drives the region\n"
            "presets, and overrides whatever the geometry guessed.")
        self.name.lineEdit().editingFinished.connect(self._rename)
        self.name.activated.connect(lambda _i: self._rename())

        self.mode = QtWidgets.QComboBox()
        for text, val in MODES:
            self.mode.addItem(text, val)
        self.mode.setCurrentIndex(0 if terminal.mode == "z" else 1)
        self.mode.setToolTip(
            "Z-plane: a horizontal cut, matches how the aorta ends were cut.\n"
            "Centerline-normal: perpendicular to the vessel, right for branch tips.")
        self.mode.currentIndexChanged.connect(self._mode_changed)

        self.slider, self.box = _pair(MIN_OFFSET_MM, self.limit, DEFAULT_OFFSET_MM,
                                      1.0, " mm")
        self.box.setToolTip(
            f"Exact trim distance ({MIN_OFFSET_MM:.0f} – {self.limit:.0f} mm)")
        self.slider.valueChanged.connect(self._slider_moved)   # live
        self.slider.sliderReleased.connect(self.changed.emit)  # commit
        self.box.valueChanged.connect(self._box_typed)

        # "Pull it back" with a number attached. A cut too close to the junction has no
        # clean cross-section to take, and the advice to move it was useless without
        # saying how far -- which the cutter can work out, by trying.
        self.back_btn = QtWidgets.QToolButton()
        self.back_btn.setText("↑")
        self.back_btn.setAutoRaise(True)
        self.back_btn.setToolTip(
            "Move this cut back to the deepest offset that still has a clean\n"
            "cross-section of its own. Only ever moves it away from the junction.")
        self.back_btn.clicked.connect(lambda: self.pullback.emit(self))

        # sliderReleased only follows a DRAG. Arrow keys, the scroll wheel and clicks
        # in the groove all change the value without it, so nothing committed them:
        # the plane moved on screen while the mesh stayed cut at the old value and the
        # status bar kept its previous (often green) verdict -- a stale cut, exportable.
        # Debounced, because a held arrow key would otherwise re-cut on every repeat.
        self._commit_timer = QtCore.QTimer(self)
        self._commit_timer.setSingleShot(True)
        self._commit_timer.setInterval(220)
        self._commit_timer.timeout.connect(self.changed.emit)
        # Starting a drag cancels a pending debounce: a groove click arms the timer
        # and the user may drag straight after, and a re-cut landing mid-drag is
        # seconds of stutter for a value that is about to change again anyway.
        for _s in (self.slider,):
            _s.sliderPressed.connect(self._commit_timer.stop)

        self.readout = QtWidgets.QLabel("-")
        self.readout.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.readout.setMinimumWidth(90)
        self.readout.setStyleSheet("color:#555;")

        top.addWidget(self.expand, 0, 0)
        top.addWidget(self.name, 0, 1)
        top.addWidget(self.readout, 0, 2)
        top.addWidget(self.slider, 1, 1)
        top.addWidget(self.box, 1, 2)
        top.addWidget(self.back_btn, 1, 3)
        top.addWidget(self.mode, 2, 1, 1, 2)
        top.setColumnStretch(1, 1)
        outer.addLayout(top)

        # ---- collapsible advanced panel ---------------------------------------
        self.adv = QtWidgets.QWidget()
        self.adv.setVisible(False)
        self.adv.setStyleSheet(
            "background:#f4f6f8; border:1px solid #e0e4e8; border-radius:4px;")
        a = QtWidgets.QGridLayout(self.adv)
        a.setContentsMargins(8, 6, 8, 6)
        a.setHorizontalSpacing(6)
        a.setVerticalSpacing(4)

        # Seeded FROM the terminal: rebuilding rows (adding or undoing a trim) must not
        # zero the widgets while the terminal quietly keeps applying its old tilt.
        _ta = float(np.clip(terminal.tilt_a_deg, -MAX_TILT_DEG, MAX_TILT_DEG))
        _tb = float(np.clip(terminal.tilt_b_deg, -MAX_TILT_DEG, MAX_TILT_DEG))
        self.tilt_a_s, self.tilt_a_b = _pair(-MAX_TILT_DEG, MAX_TILT_DEG, _ta, 1.0, "°")
        self.tilt_b_s, self.tilt_b_b = _pair(-MAX_TILT_DEG, MAX_TILT_DEG, _tb, 1.0, "°")
        for s, b, attr in [(self.tilt_a_s, self.tilt_a_b, "tilt_a_deg"),
                           (self.tilt_b_s, self.tilt_b_b, "tilt_b_deg")]:
            s.valueChanged.connect(lambda v, bb=b, ss=s: self._tilt_dragged(bb, v, ss))
            s.sliderReleased.connect(self.changed.emit)
            b.valueChanged.connect(lambda v, ss=s, at=attr: self._tilt_typed(ss, at, v))

        # In-plane shift: slide the cut's anchor sideways within its plane so its
        # footprint clears a neighbouring vessel. The band stops covering the cut
        # ring once the shift exceeds ~(reach - radius), so the slider is capped
        # there, past it the cut only tears, which is not worth exposing.
        _r = max(float(getattr(terminal.segment, "_radius_mm", None) or 5.0), 1.0)
        shift_range = max(self.auto_band - _r, 6.0)
        _sa = float(np.clip(terminal.shift_a_mm, -shift_range, shift_range))
        _sb = float(np.clip(terminal.shift_b_mm, -shift_range, shift_range))
        self.shift_a_s, self.shift_a_b = _pair(-shift_range, shift_range, _sa, 0.5, " mm")
        self.shift_b_s, self.shift_b_b = _pair(-shift_range, shift_range, _sb, 0.5, " mm")
        for s, b in [(self.shift_a_s, self.shift_a_b), (self.shift_b_s, self.shift_b_b)]:
            b.setToolTip(
                "Slide the cut sideways within its own plane, to move its reach off a\n"
                "neighbouring vessel a rough centerline would otherwise clip. Pairs with\n"
                "Tilt (angle) and Plane size (reach); it does not change the cut's depth.")
            s.valueChanged.connect(lambda v, bb=b, ss=s: self._shift_dragged(bb, v, ss))
            s.sliderReleased.connect(self.changed.emit)
            b.valueChanged.connect(lambda v, ss=s: self._shift_typed(ss, v))

        # Plane size: how far the cut reaches around its plane. Shrinking it is the
        # cure when a plane over-reaches on a tortuous branch and clips a stretch
        # that doubles back (which shows up as one cut making two holes).
        self.band_s, self.band_b = _pair(3.0, max(4.0 * self.auto_band, 30.0),
                                         self.auto_band, 1.0, " mm")
        self.band_b.setToolTip(
            "How far the cut reaches around its plane.\n"
            "SHRINK this if one cut is punching two holes, the plane is over-reaching\n"
            "and catching a part of the vessel that doubles back near it.\n"
            "Too small and there is nothing left to clip ('no cut band').")
        self.band_s.valueChanged.connect(self._band_dragged)
        self.band_s.sliderReleased.connect(self.changed.emit)
        self.band_b.valueChanged.connect(self._band_typed)

        self.reset_btn = QtWidgets.QPushButton("Reset plane")
        self.reset_btn.clicked.connect(self._reset_plane)

        a.addWidget(QtWidgets.QLabel("Tilt A"), 0, 0)
        a.addWidget(self.tilt_a_s, 0, 1)
        a.addWidget(self.tilt_a_b, 0, 2)
        a.addWidget(QtWidgets.QLabel("Tilt B"), 1, 0)
        a.addWidget(self.tilt_b_s, 1, 1)
        a.addWidget(self.tilt_b_b, 1, 2)
        a.addWidget(QtWidgets.QLabel("Shift A"), 2, 0)
        a.addWidget(self.shift_a_s, 2, 1)
        a.addWidget(self.shift_a_b, 2, 2)
        a.addWidget(QtWidgets.QLabel("Shift B"), 3, 0)
        a.addWidget(self.shift_b_s, 3, 1)
        a.addWidget(self.shift_b_b, 3, 2)
        a.addWidget(QtWidgets.QLabel("Plane size"), 4, 0)
        a.addWidget(self.band_s, 4, 1)
        a.addWidget(self.band_b, 4, 2)
        a.addWidget(self.reset_btn, 5, 2)
        a.setColumnStretch(1, 1)
        outer.addWidget(self.adv)

        if self.bound:
            # The region's sliders own where this plane sits, and its mode is fixed
            # by region_plan; leaving dead controls on screen invites the user to
            # move a plane by a route that does not drive it.
            for w in (self.slider, self.box, self.back_btn, self.mode):
                w.setVisible(False)
            self.expand.setChecked(True)          # nothing else is in this row
            self.name.setToolTip(
                "A boundary of the isolated region. Move it with the region's "
                "Start/End sliders above; tilt and reach are here.")

        self.setStyleSheet("QLineEdit { background: transparent; }")

    # -- offset -------------------------------------------------------------

    def value(self):
        return float(self.box.value())

    def set_value(self, mm):
        self._syncing = True
        mm = float(np.clip(mm, MIN_OFFSET_MM, self.limit))
        self.box.setValue(mm)
        self.slider.setValue(int(round(mm * 10)))
        self._syncing = False
        # Whoever set this value owns the commit; a debounce armed for the previous
        # value must not fire afterwards against a number that is no longer shown.
        self._commit_timer.stop()

    def _commit_if_not_dragging(self, slider):
        """Commit a slider change that did not come from dragging the handle."""
        if not slider.isSliderDown():
            self._commit_timer.start()

    def _slider_moved(self, v):
        if self._syncing:
            return
        self._syncing = True
        self.box.setValue(v / 10.0)
        self._syncing = False
        self.preview.emit(self)          # plane follows the drag; mesh waits
        self._commit_if_not_dragging(self.slider)

    def _box_typed(self, v):
        if self._syncing:
            return
        self._syncing = True
        self.slider.setValue(int(round(v * 10)))
        self._syncing = False
        self._commit_timer.stop()    # this commit supersedes any pending debounce
        self.changed.emit()          # typing a number IS a commit

    # -- advanced -----------------------------------------------------------

    def _toggle_advanced(self, on):
        self.adv.setVisible(on)
        self.expand.setArrowType(QtCore.Qt.DownArrow if on else QtCore.Qt.RightArrow)
        if on:
            self.selected.emit(self)
        self.preview.emit(self)

    def _tilt_dragged(self, box, v, slider=None):
        if self._syncing:
            return
        self._syncing = True
        box.setValue(v / 10.0)
        self._syncing = False
        self._push_tilt()
        self.preview.emit(self)
        if slider is not None:
            self._commit_if_not_dragging(slider)

    def _tilt_typed(self, slider, attr, v):
        if self._syncing:
            return
        self._syncing = True
        slider.setValue(int(round(v * 10)))
        self._syncing = False
        self._push_tilt()
        self.changed.emit()

    def _push_tilt(self):
        self.terminal.tilt_a_deg = float(self.tilt_a_b.value())
        self.terminal.tilt_b_deg = float(self.tilt_b_b.value())

    def _shift_dragged(self, box, v, slider=None):
        if self._syncing:
            return
        self._syncing = True
        box.setValue(v / 10.0)
        self._syncing = False
        self._push_shift()
        self.preview.emit(self)
        if slider is not None:
            self._commit_if_not_dragging(slider)

    def _shift_typed(self, slider, v):
        if self._syncing:
            return
        self._syncing = True
        slider.setValue(int(round(v * 10)))
        self._syncing = False
        self._push_shift()
        self.changed.emit()

    def _push_shift(self):
        self.terminal.shift_a_mm = float(self.shift_a_b.value())
        self.terminal.shift_b_mm = float(self.shift_b_b.value())

    def _band_dragged(self, v):
        if self._syncing:
            return
        self._syncing = True
        self.band_b.setValue(v / 10.0)
        self._syncing = False
        self.terminal.band_mm = float(self.band_b.value())
        self.preview.emit(self)
        self._commit_if_not_dragging(self.band_s)

    def _band_typed(self, v):
        if self._syncing:
            return
        self._syncing = True
        self.band_s.setValue(int(round(v * 10)))
        self._syncing = False
        self.terminal.band_mm = float(v)
        self.changed.emit()

    def _reset_plane(self):
        self.terminal.reset_plane()
        self._syncing = True
        for s, b in [(self.tilt_a_s, self.tilt_a_b), (self.tilt_b_s, self.tilt_b_b),
                     (self.shift_a_s, self.shift_a_b), (self.shift_b_s, self.shift_b_b)]:
            b.setValue(0.0)
            s.setValue(0)
        self.band_b.setValue(self.auto_band)
        self.band_s.setValue(int(round(self.auto_band * 10)))
        self._syncing = False
        self.terminal.band_mm = None      # back to automatic
        self._commit_timer.stop()         # this commit supersedes any pending one
        self.changed.emit()

    def _mode_changed(self, _):
        self.terminal.mode = self.mode.currentData()
        self._commit_timer.stop()
        self.changed.emit()

    def _rename(self):
        text = self.name.currentText().strip() or self.terminal.label
        was = getattr(self.terminal, "name", "")
        self.terminal.name = text
        if text != was:
            # Only the window knows whether this text is one of the canonical artery
            # names, and only it can rebuild the presets that follow from it.
            self.renamed.emit(self)

    def declared_artery(self):
        """The artery this row DECLARES, or "" if the text is just a label.

        An exact match against the offered list, not a fuzzy one. "left renal" and
        "Left renal artery" would both read as an assignment to a person, and only
        the second is a name region_presets can find -- so a near miss must not be
        silently promoted into a declaration that then fails to do anything.
        """
        text = self.name.currentText().strip()
        return text if text in self.choices else ""

    def set_trim(self, on_remove):
        """Mark this row as a trim (a cut before a split) and offer to undo it."""
        self.undo_trim = QtWidgets.QToolButton()
        self.undo_trim.setText("✕")
        self.undo_trim.setAutoRaise(True)
        self.undo_trim.setToolTip(
            "Remove this trim, the openings it replaced come back as separate cuts.")
        self.undo_trim.clicked.connect(on_remove)
        self.layout().itemAt(0).addWidget(self.undo_trim, 0, 3)
        self.name.setToolTip(
            "One cut on the vessel feeding this group, placed before it splits.\n"
            "Its distance is measured back from the split.")

    def set_colour(self, hexcol):
        """Tint the opening's name to match its branch colour, or clear it.

        Applied to the combo's line edit rather than the combo, so the drop-down
        arrow keeps its own styling -- a stylesheet on the box itself repaints the
        whole widget and the arrow disappears with it.
        """
        if hexcol:
            self.name.lineEdit().setStyleSheet(
                f"background:transparent; color:{hexcol}; font-weight:bold;")
        else:
            self.name.lineEdit().setStyleSheet("background:transparent;")

    def show_result(self, stats):
        self.readout.setText("-" if stats is None else f"⌀ {stats['diameter']:.1f} mm")
        if stats is not None:
            self.readout.setToolTip("Diameter of the cut opening.")

    def flag_pullback(self, on):
        """Make the pull-back button obvious on the row that is actually failing.

        It is a 22-pixel button at the end of a row and there is no reason anyone
        would hunt for it while reading a red status bar. When this row is the one
        being complained about, it says so.
        """
        self.back_btn.setStyleSheet(
            "QToolButton { background:#c62828; color:white; border-radius:3px;"
            " font-weight:bold; }" if on else "")
        self.back_btn.setToolTip(
            "This cut has no clean cross-section here.\n"
            "Move it back to the deepest offset that does."
            if on else
            "Move this cut back to the deepest offset that still has a clean\n"
            "cross-section of its own. Only ever moves it away from the junction.")

    def show_live_diameter(self, sd):
        """Live cross-section size at the plane while a slider is dragged.

        Headline is the equivalent-area diameter (2*sqrt(area/pi)); the min-max
        caliper range is appended when the section is noticeably non-circular so
        the ovality is visible, and the full breakdown is in the tooltip.
        """
        if sd is None:
            self.readout.setText("-")
            return
        ead, dmin, dmax = sd["ead_mm"], sd["dmin_mm"], sd["dmax_mm"]
        if dmax - dmin > 0.12 * max(ead, 1e-6):
            self.readout.setText(f"⌀ {ead:.1f} ({dmin:.0f}–{dmax:.0f})")
        else:
            self.readout.setText(f"⌀ {ead:.1f} mm")
        self.readout.setToolTip(
            f"At this cut plane:\n"
            f"  equivalent-area diameter {ead:.1f} mm\n"
            f"  width {dmin:.1f}–{dmax:.1f} mm (min–max caliper)\n"
            f"  circularity {sd['circularity']:.2f}  (1.0 = a circle)")


class Section(QtWidgets.QWidget):
    """A panel block that folds down to one clickable header line.

    The tool has grown enough controls that showing them all at once squeezed the
    openings list -- the surface the user actually works on -- down to a sliver.
    Setup blocks are touched once per model, so they collapse after loading and
    show a summary of what they are set to, giving the list the room instead.
    """

    def __init__(self, title, body, key="", settings=None, summary_fn=None):
        super().__init__()
        self._title, self._key, self._settings = title, key, settings
        self._summary = ""
        # Every control feeding a summary lives inside the block it summarises, so it
        # cannot change while folded -- computing the text at fold time is enough, and
        # saves sprinkling refresh calls through a dozen unrelated handlers.
        self._summary_fn = summary_fn

        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        self.header = QtWidgets.QToolButton()
        self.header.setCheckable(True)
        self.header.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        self.header.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                  QtWidgets.QSizePolicy.Fixed)
        self.header.setCursor(QtCore.Qt.PointingHandCursor)
        self.header.setStyleSheet(
            "QToolButton { border:none; background:#eef1f4; border-radius:3px;"
            " padding:4px 6px; font-weight:600; text-align:left; }"
            "QToolButton:hover { background:#e2e7ec; }")
        self.header.toggled.connect(self._apply)
        v.addWidget(self.header)

        self.body = body
        v.addWidget(self.body)

        want = True
        if settings is not None and key:
            want = str(settings.value(f"sec_{key}", "true")).lower() == "true"
        self.header.setChecked(want)
        self._apply(want)

    def _apply(self, on):
        self.body.setVisible(bool(on))
        self.header.setArrowType(QtCore.Qt.DownArrow if on else QtCore.Qt.RightArrow)
        if not on and self._summary_fn is not None:
            try:
                self._summary = self._summary_fn() or ""
            except Exception:
                self._summary = ""     # a summary is a nicety; never block the fold
        # The summary is the whole point of collapsing: folded away, the header still
        # has to say what the block is set to.
        self.header.setText(self._title if (on or not self._summary)
                            else f"{self._title}     {self._summary}")
        if self._settings is not None and self._key:
            self._settings.setValue(f"sec_{self._key}", "true" if on else "false")

    def set_summary(self, text):
        self._summary = text or ""
        self._apply(self.header.isChecked())

    def refresh(self):
        """Recompute the folded summary in place, without changing fold state."""
        self._apply(self.header.isChecked())

    def set_open(self, on):
        # setChecked emits nothing when the state is already correct, so the toggled
        # -> _apply connection would not fire and the summary (recomputed only in
        # _apply) would keep asserting a setting that has since changed underneath it.
        # Call _apply directly; it is idempotent.
        self.header.setChecked(bool(on))
        self._apply(bool(on))


class CutWindow(QtWidgets.QMainWindow):
    def __init__(self, csv=None, mat=None):
        super().__init__()
        self.setWindowTitle("Vascular cut tool")
        self.resize(1560, 940)

        self.settings = QtCore.QSettings("branches", "cuttool")
        self.csv_path = Path(csv) if csv else None
        self.mat_path = Path(mat) if mat else None
        self.out_dir = Path(self.settings.value("out_dir", str(Path.cwd() / "output")))

        self.tree = None
        self.base = None
        self.cut = None
        self.matched = {}
        self.problem = None
        self.rows = []
        self.bound_rows = []
        self.sel = None
        self.opening_colour = {}          # terminal label -> hex colour (row order)
        self.trims = []                   # active pre-bifurcation ("trim") cuts
        self.anat_names, self.anat_vessels, self.anat_notes = {}, {}, []
        # The measure plane: the stretch a measurement covered, where its disc sits,
        # and the widest station in it. Kept on the window rather than in the panel so
        # a re-cut can put the disc back after the scene is cleared.
        self._m_range = (0.0, 0.0)
        self._m_station = None
        self._m_widest = None
        self._m_pose = None
        self.region = None                # (s0, s1) along the aorta, or None
        self._region_name = ""            # which preset it came from, for the banner
        self.route = None
        self.mesh_actor = None

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)

        # Setup blocks fold away once a model is loaded; the openings list keeps the
        # room. See Section.
        self.sec_input = Section("Input", self._files_group(), "input", self.settings,
                                 self._sum_input)
        self.sec_cl = Section("Centerline", self._centerline_group(), "cl",
                              self.settings, self._sum_cl)
        self.sec_view = Section("View & measure", self._view_group(), "view",
                                self.settings, self._sum_view)
        self.sec_out = Section("Output", self._output_group(), "out", self.settings,
                               self._sum_out)

        self.openings_box = QtWidgets.QGroupBox("Openings")
        self.openings_lay = QtWidgets.QVBoxLayout(self.openings_box)
        self.openings_lay.setSpacing(2)
        self.hint = QtWidgets.QLabel(
            "Load a mesh to view it.\nAdd a centerline to cut it.")
        self.hint.setStyleSheet("color:#777; padding:12px;")
        self.openings_lay.addWidget(self.hint)

        # One scroll for the whole control column, rather than a scroll around the
        # openings alone. Nothing can then be pushed off the bottom of a short screen,
        # and the window is free to be much shorter than the controls it holds.
        stack = QtWidgets.QWidget()
        sl = QtWidgets.QVBoxLayout(stack)
        sl.setContentsMargins(10, 10, 6, 4)
        sl.setSpacing(6)
        for s in (self.sec_input, self.sec_cl, self.sec_view, self.sec_out):
            sl.addWidget(s)
        sl.addWidget(self.openings_box)
        sl.addWidget(self._branch_bar())
        sl.addStretch(1)

        self.panel_scroll = QtWidgets.QScrollArea()
        self.panel_scroll.setWidgetResizable(True)
        self.panel_scroll.setWidget(stack)
        self.panel_scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.panel_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)

        panel = QtWidgets.QWidget()
        # Wide enough for the branch bar ("Move all branches ▲ toward tips ▼ deeper"),
        # which sets the column's layout minimum -- below this the horizontal scrollbar
        # is off and the right-hand controls would simply be clipped.
        panel.setMinimumWidth(420)
        pl = QtWidgets.QVBoxLayout(panel)
        pl.setContentsMargins(0, 0, 0, 0)
        pl.setSpacing(6)
        pl.addWidget(self.panel_scroll, 1)

        # Status and the two actions are pinned below the scroll: whether the cut is
        # valid, and the button that exports it, must never scroll out of sight.
        foot = QtWidgets.QWidget()
        fl = QtWidgets.QVBoxLayout(foot)
        fl.setContentsMargins(10, 0, 10, 10)
        fl.setSpacing(6)
        self.status = QtWidgets.QLabel("No model loaded")
        self.status.setWordWrap(True)
        self.status.setAlignment(QtCore.Qt.AlignCenter)
        self.status.setMinimumHeight(46)
        fl.addWidget(self.status)

        btns = QtWidgets.QHBoxLayout()
        self.reset_btn = QtWidgets.QPushButton("Reset")
        self.reset_btn.setToolTip(
            "Put every cut plane back to its default depth and angle, and drop any\n"
            "isolated region. Trimmed branch groups are kept.")
        self.reset_btn.clicked.connect(self.reset)
        self.export_btn = QtWidgets.QPushButton("Export open mesh…")
        self.export_btn.setDefault(True)
        self.export_btn.clicked.connect(self.export)
        btns.addWidget(self.reset_btn)
        btns.addWidget(self.export_btn, 1)
        fl.addLayout(btns)
        pl.addWidget(foot)

        self.plotter = QtInteractor(central)
        self.plotter.set_background("white")
        # An anatomical orientation triad, labelled in ANATOMY_FRAME's own terms
        # rather than x/y/z. Every Left/Right name the tool prints follows from that
        # frame, and an anterior view puts the patient's LEFT on the viewer's RIGHT
        # -- so without a marker on screen there is nothing to check a side against,
        # and a user reading "left iliac" off the wrong half of the screen has no way
        # to notice. Reported from the cohort as iliacs coming out swapped; this does
        # not decide whether the tool or the eye was wrong, it makes the answer
        # visible either way. L/P/S are the +x/+y/+z of ANATOMY_FRAME.
        try:
            self.plotter.add_axes(xlabel="L", ylabel="P", zlabel="S",
                                  color="#333333")
        except Exception:                                       # noqa: BLE001
            pass          # an orientation widget is never worth failing a load over

        # A splitter instead of a fixed 450 px: long anatomical names and deep region
        # readouts need more width on some cases, less on others.
        self.split = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.split.addWidget(panel)
        self.split.addWidget(self.plotter.interactor)
        self.split.setStretchFactor(0, 0)
        self.split.setStretchFactor(1, 1)
        self.split.setChildrenCollapsible(False)
        # QSettings hands back strings on Windows, and a width saved on a 4K monitor
        # would swallow a laptop screen -- so parse defensively and clamp.
        try:
            pw = int(float(self.settings.value("panel_w", 450)))
        except (TypeError, ValueError):
            pw = 450
        self.split.setSizes([max(420, min(pw, 900)), 1100])
        self.split.splitterMoved.connect(
            lambda *_: self.settings.setValue("panel_w", self.split.sizes()[0]))
        root.addWidget(self.split)
        self._init_hover()

        self._set_enabled(False)
        if self.mat_path:
            self.load_model()
        else:
            # Nothing loaded: Browse/Load must be in front of the user whatever state
            # the sections were left in last session.
            self.sec_input.set_open(True)

    # ---------------- hover measuring ----------------

    def _init_hover(self):
        """Wire a mouse-move observer + surface picker for the hover readout.

        The observers are always attached but do nothing unless the checkbox is on
        and a centerline is loaded, so they cost a single attribute read otherwise.
        Picking is skipped while a mouse button is held (i.e. while rotating), so
        camera control stays smooth.
        """
        import vtk

        self._hover_picker = vtk.vtkCellPicker()
        self._hover_picker.SetTolerance(0.0008)
        self._hover_busy = False
        self._hover_dragging = False
        self._hover_shown = False
        self._hover_last = None
        self._hover_marker_r = 2.0
        self._marks = []                 # pinned measurements (survive re-cuts)
        self._lclick_xy = None

        iren = self.plotter.iren
        iren.add_observer("MouseMoveEvent", self._on_hover)
        # Left button gets its own handlers: a press+release with no drag is a
        # CLICK, which pins a mark. Right/middle only set the drag flag.
        iren.add_observer("LeftButtonPressEvent", self._left_press)
        iren.add_observer("LeftButtonReleaseEvent", self._left_release)
        for ev in ("RightButtonPressEvent", "MiddleButtonPressEvent"):
            iren.add_observer(ev, self._hover_press)
        for ev in ("RightButtonReleaseEvent", "MiddleButtonReleaseEvent"):
            iren.add_observer(ev, self._hover_release)

    def _hover_press(self, *args):
        self._hover_dragging = True

    def _hover_release(self, *args):
        self._hover_dragging = False

    def _left_press(self, *args):
        self._hover_dragging = True
        try:
            self._lclick_xy = self.plotter.iren.get_event_position()
        except Exception:
            self._lclick_xy = None

    def _left_release(self, *args):
        self._hover_dragging = False
        if not self.measure_cb.isChecked() or self.tree is None or self._lclick_xy is None:
            return
        try:
            xy = self.plotter.iren.get_event_position()
        except Exception:
            return
        moved = (abs(xy[0] - self._lclick_xy[0]) + abs(xy[1] - self._lclick_xy[1]))
        if moved <= 4:                   # a click, not a drag-to-rotate
            self._place_mark(xy)

    def _measure_toggled(self, on):
        self.settings.setValue("measure", "true" if on else "false")
        if not on:
            self._clear_hover()

    def _on_hover(self, *args):
        if (not self.measure_cb.isChecked() or self.tree is None
                or self.mesh_actor is None or self._hover_dragging
                or self._hover_busy):
            return
        self._hover_busy = True
        try:
            x, y = self.plotter.iren.get_event_position()
            picker = self._hover_picker
            hit = picker.Pick(x, y, 0, self.plotter.renderer)
            if not hit or picker.GetCellId() < 0:
                self._clear_hover()
                return
            world = np.asarray(picker.GetPickPosition(), float)
            if self._hover_last is not None and \
                    np.linalg.norm(world - self._hover_last) < 0.4:
                return                       # negligible move; leave the readout be
            self._hover_last = world
            res = measure_along_branch(self.tree, world)
            if res is None:
                self._clear_hover()
                return
            self._draw_hover(res)
        finally:
            self._hover_busy = False

    @staticmethod
    def _on_top(actor):
        """Let *actor* show through the mesh at any opacity.

        Every measurement marker is anchored on the centerline, i.e. inside the
        lumen. Below 100% the mesh renders in VTK's translucent pass and writes no
        depth, so the marker shows through; at exactly 100% the mesh is opaque,
        writes depth, and the marker is culled to nothing. A polygon offset biases
        the marker forward in depth only -- its true position is unchanged, and it
        stays a normal actor in renderer.actors so the usual removal path works.
        PyVista biases its own widget handles the same way.
        """
        try:
            m = actor.mapper
            m.SetResolveCoincidentTopologyToPolygonOffset()
            m.SetRelativeCoincidentTopologyPolygonOffsetParameters(0, -66000)
        except Exception:
            pass                 # a cosmetic nicety; never break drawing over it
        return actor

    def _draw_hover(self, res):
        p = self.plotter
        r = self._hover_marker_r
        node_word = "node" if res["node_kind"] == "junction" else "end"

        self._on_top(p.add_mesh(pv.Sphere(radius=1.15 * r, center=res["node_pos"]),
                                color=HOVER_NODE, name="hover_node", pickable=False))
        # The dot marks the point ON THE CENTERLINE that the reading refers to --
        # not the spot on the wall the ray happened to strike. Those differ by the
        # local vessel radius, which is ~20mm in the trunk: a wall dot reads as
        # belonging to neither the centerline nor the number beside it.
        self._on_top(p.add_mesh(pv.Sphere(radius=0.8 * r, center=res["proj"]),
                                color=HOVER_PT, name="hover_pt", pickable=False))
        if len(res["path"]) >= 2:
            self._on_top(p.add_mesh(pv.Spline(res["path"], max(len(res["path"]), 50)),
                                    color=HOVER_PATH, line_width=7, name="hover_path",
                                    pickable=False))
        p.add_text(
            f"{res['branch_name']}\n{res['distance_mm']:.1f} mm to {node_word}",
            position="upper_left", font_size=13, color=HOVER_NODE, name="hover_text")
        self.measure_lbl.setText(
            f"{res['branch_name']}: {res['distance_mm']:.1f} mm to {node_word}")
        self._hover_shown = True
        p.render()

    def _clear_hover(self):
        self._hover_last = None
        if not self._hover_shown:
            return
        for n in ("hover_node", "hover_pt", "hover_path", "hover_text"):
            try:
                self.plotter.remove_actor(n, render=False)
            except Exception:
                pass
        self.measure_lbl.setText("")
        self._hover_shown = False
        self.plotter.render()

    # ---------------- pinned marks ----------------

    def _place_mark(self, xy):
        """Pin the measurement at screen position *xy* as a persistent mark."""
        picker = self._hover_picker
        if not picker.Pick(xy[0], xy[1], 0, self.plotter.renderer) \
                or picker.GetCellId() < 0:
            return
        world = np.asarray(picker.GetPickPosition(), float)
        res = measure_along_branch(self.tree, world)
        if res is None:
            return
        self._marks.append({"world": world, "res": res})
        self._draw_marks()
        self._update_marks_label()

    def _draw_marks(self):
        """(Re)draw every pinned mark and the span between the last two."""
        p = self.plotter
        for n in list(p.renderer.actors):
            if n.startswith("mark_"):
                try:
                    p.remove_actor(n, render=False)
                except Exception:
                    pass
        if not self._marks:
            return

        r = self._hover_marker_r
        pts, labels = [], []
        for i, mk in enumerate(self._marks, 1):
            # Anchor on the centerline, not on the wall the pick ray hit: the label
            # is an along-centerline distance, so the dot has to sit where that
            # distance is actually measured from. Dot and label share the anchor --
            # move one without the other and they drift apart by the vessel radius.
            anchor = mk["res"]["proj"]
            self._on_top(p.add_mesh(pv.Sphere(radius=0.9 * r, center=anchor),
                                    color=MARK_COLOR, name=f"mark_pt_{i}",
                                    pickable=False))
            pts.append(anchor)
            labels.append(f"{i}: {mk['res']['distance_mm']:.1f} mm")
        p.add_point_labels(np.array(pts), labels, name="mark_labels",
                           font_size=12, text_color=MARK_COLOR, bold=True,
                           shape=None, show_points=False, always_visible=True,
                           pickable=False)
        # The along-centerline span between the last two marks is reported as a
        # number in the panel; we deliberately do NOT draw a connecting line,
        # which reads as visual clutter over the model.

    def _update_marks_label(self):
        n = len(self._marks)
        self.send_mark_btn.setEnabled(n > 0)
        if n == 0:
            self.marks_lbl.setText("")
            self.clear_marks_btn.setEnabled(False)
            return
        self.clear_marks_btn.setEnabled(True)
        if n == 1:
            self.marks_lbl.setText(
                f"1 mark, {self._marks[0]['res']['distance_mm']:.1f} mm to node")
            return
        g = centerline_geodesic(self.tree, self._marks[-2]["world"],
                                self._marks[-1]["world"])
        span = f"{g['distance_mm']:.1f} mm" if g else "-"
        self.marks_lbl.setText(
            f"{n} marks, last two: {span} apart along the centerline")

    def _clear_marks(self):
        self._marks = []
        self._draw_marks()
        self._update_marks_label()
        self.plotter.render()

    # ---------------- collapsed-section summaries ----------------
    # Each says what its block is set to, so folding one away never hides a setting
    # that would change the result.

    def _sum_input(self):
        mesh = Path(self.mat_label.text()).name if self.mat_label.text() else "no mesh"
        cl = Path(self.csv_label.text()).name if self.csv_label.text() else "mesh only"
        return f"{mesh}  ·  {cl}"

    def _sum_cl(self):
        if not self.smooth_cb.isChecked():
            return "as supplied"
        return f"smoothed {self.smooth_mm.value():.1f} mm"

    def _sum_view(self):
        bits = [f"opacity {self.opacity.value()}%"]
        if self.colour_cb.isChecked():
            bits.append("coloured")
        if self.measure_cb.isChecked():
            bits.append("measuring")
        return "  ·  ".join(bits)

    def _sum_out(self):
        return f"{Path(self.out_edit.text()).name or '?'}  ·  .{self.fmt.currentText()}"

    # ---------------- groups ----------------

    def _files_group(self):
        g = QtWidgets.QWidget()          # titled by its Section header
        lay = QtWidgets.QGridLayout(g)
        lay.setContentsMargins(6, 4, 2, 4)
        lay.setHorizontalSpacing(6)

        # Typed as well as browsed. The output folder has always been a text box and
        # these two were not, so a cohort of 54 scans in sibling folders meant 54
        # trips through a file dialog to change one path component. Paste-and-go is
        # also how anyone works from a spreadsheet of case paths.
        self.mat_label = QtWidgets.QLineEdit(str(self.mat_path or ""))
        self.mat_label.setPlaceholderText("surface mesh .mat / .stl / .ply  (required)")
        self.mat_label.setToolTip(
            "Type or paste a path, or use Browse…. Quotes and stray spaces are "
            "trimmed, so a path copied from Explorer can go straight in.")
        self.mat_label.editingFinished.connect(self._mat_typed)
        mat_btn = QtWidgets.QPushButton("Browse…")
        mat_btn.clicked.connect(self._pick_mat)

        self.csv_label = QtWidgets.QLineEdit(str(self.csv_path or ""))
        self.csv_label.setPlaceholderText(
            "centerline .csv / .vtp / .json  (optional, leave blank to just view)")
        self.csv_label.setToolTip(
            "Type or paste a path, or use Browse…. Leave it blank to view the mesh "
            "without cutting it.")
        self.csv_label.editingFinished.connect(self._csv_typed)
        csv_btn = QtWidgets.QPushButton("Browse…")
        csv_btn.clicked.connect(self._pick_csv)
        csv_clear = QtWidgets.QPushButton("✕")
        csv_clear.setFixedWidth(26)
        csv_clear.setToolTip("Clear the centerline, view the mesh only")
        csv_clear.clicked.connect(self._clear_csv)

        self.load_btn = QtWidgets.QPushButton("Load")
        self.load_btn.clicked.connect(self.load_model)

        lay.addWidget(QtWidgets.QLabel("Mesh"), 0, 0)
        lay.addWidget(self.mat_label, 0, 1)
        lay.addWidget(mat_btn, 0, 2)
        lay.addWidget(QtWidgets.QLabel("Centerline"), 1, 0)
        lay.addWidget(self.csv_label, 1, 1)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(2)
        row.addWidget(csv_btn)
        row.addWidget(csv_clear)
        w = QtWidgets.QWidget()
        w.setLayout(row)
        lay.addWidget(w, 1, 2)
        lay.addWidget(self.load_btn, 2, 1, 1, 2)
        lay.setColumnStretch(1, 1)
        return g

    def _centerline_group(self):
        """De-staircasing, off by default and never silent.

        A Simpleware spline export is already a fitted curve (0.6-2.8 deg/mm of
        turning) and smoothing it would only degrade it. A cl.py voxel-geodesic
        path is a staircase (75 deg/mm) and needs it. So: measure, report, and let
        the user decide, do not quietly reshape anyone's centerline.
        """
        g = QtWidgets.QWidget()          # titled by its Section header
        lay = QtWidgets.QGridLayout(g)
        lay.setContentsMargins(6, 4, 2, 4)
        lay.setHorizontalSpacing(6)

        self.quality_lbl = QtWidgets.QLabel("-")
        self.quality_lbl.setWordWrap(True)
        self.quality_lbl.setStyleSheet("color:#777; font-size:11px;")

        self.smooth_cb = QtWidgets.QCheckBox("Smooth (de-staircase)")
        self.smooth_cb.setToolTip(
            "OFF: the centerline is used exactly as supplied.\n\n"
            "ON: fits a smoothing spline through it, removing the voxel staircase\n"
            "left by a geodesic path (cl.py). On the sample this cut the tangent\n"
            "turning from 75 to 3 °/mm AND left the line deeper inside the lumen\n"
            "(min wall distance 3.67 → 4.74 mm), so it is not cutting corners.\n\n"
            "Leave it OFF for Simpleware/Scan-IP splines, they are already curves\n"
            "and smoothing can only degrade them.")
        self.smooth_cb.toggled.connect(self._smooth_toggled)

        self.smooth_mm = QtWidgets.QDoubleSpinBox()
        self.smooth_mm.setRange(0.2, 5.0)
        self.smooth_mm.setDecimals(1)
        self.smooth_mm.setSingleStep(0.5)
        self.smooth_mm.setValue(float(self.settings.value("smooth_mm", 1.0)))
        self.smooth_mm.setSuffix(" mm")
        self.smooth_mm.setKeyboardTracking(False)
        self.smooth_mm.setToolTip(
            "How far the smoothed curve may deviate from the supplied points (RMS).\n"
            "1 mm de-staircases without leaving the centre. Above ~4 mm it starts\n"
            "cutting corners and drifting toward the wall.")
        self.smooth_mm.valueChanged.connect(self._smooth_amount_changed)

        # Which anatomy this is decides how the openings are named and which regions
        # are offered, and it is inferred from one radius ratio. On 40 real aortas
        # that inference was wrong 11 times -- naming abdominal aortas "Ascending
        # aorta" and offering them Ishimaru zones. Declaring it fixed all 11, so the
        # guess stays as the default but stops being the last word.
        self.family = QtWidgets.QComboBox()
        self.family.addItem("Auto-detect", None)
        self.family.addItem("Abdominal, one inlet", 1)
        self.family.addItem("Thoracic, root + distal end", 2)
        self.family.setToolTip(
            "Whether this is an abdominal or a thoracic aorta.\n\n"
            "Auto-detect reads it from the relative size of the largest openings,\n"
            "which is usually right but cannot always tell an abdominal aorta with\n"
            "prominent iliacs from a thoracic one. Set it explicitly if the opening\n"
            "names or the offered regions look like the wrong anatomy.\n\n"
            "Changing this re-reads the centerline.")
        self.family.currentIndexChanged.connect(self._family_changed)

        lay.addWidget(self.smooth_cb, 0, 0)
        lay.addWidget(self.smooth_mm, 0, 1)
        lay.addWidget(self.quality_lbl, 1, 0, 1, 2)
        lay.addWidget(QtWidgets.QLabel("Anatomy"), 2, 0)
        lay.addWidget(self.family, 2, 1)
        lay.setColumnStretch(0, 1)
        return g

    def _view_group(self):
        g = QtWidgets.QWidget()          # titled by its Section header
        lay = QtWidgets.QGridLayout(g)
        lay.setContentsMargins(6, 4, 2, 4)
        lay.addWidget(QtWidgets.QLabel("Opacity"), 0, 0)
        self.opacity = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.opacity.setRange(5, 100)
        self.opacity.setValue(int(self.settings.value("opacity", 55)))
        self.opacity.valueChanged.connect(self._opacity_changed)
        self.opacity_lbl = QtWidgets.QLabel(f"{self.opacity.value()}%")
        self.opacity_lbl.setFixedWidth(38)
        lay.addWidget(self.opacity, 0, 1)
        lay.addWidget(self.opacity_lbl, 0, 2)

        self.colour_cb = QtWidgets.QCheckBox("Colour branches")
        self.colour_cb.setChecked(self.settings.value("colour", "false") == "true")
        self.colour_cb.setToolTip(
            "Tint each branch's vessel, its cut-plane disc, and its row label a\n"
            "distinct colour, so it is obvious which branch a plane belongs to.\n"
            "The aorta trunk stays red.\n\n"
            "Needs a centerline loaded.")
        self.colour_cb.toggled.connect(self._colour_toggled)
        lay.addWidget(self.colour_cb, 1, 0, 1, 3)

        self.measure_cb = QtWidgets.QCheckBox("Measure along branch (hover, click to pin)")
        self.measure_cb.setChecked(self.settings.value("measure", "false") == "true")
        self.measure_cb.setToolTip(
            "HOVER the cursor over the mesh to read the distance, in mm, along the\n"
            "centerline from that spot back to the node its branch springs from\n"
            "(for a trunk, to the nearest junction).\n\n"
            "CLICK to pin that spot as a mark that stays. With two or more marks, the\n"
            "distance ALONG THE CENTERLINE between the last two is shown too.\n\n"
            "Needs a centerline loaded.")
        self.measure_cb.toggled.connect(self._measure_toggled)
        lay.addWidget(self.measure_cb, 2, 0, 1, 3)

        self.measure_lbl = QtWidgets.QLabel("")
        self.measure_lbl.setStyleSheet("color:#7b2fbf; font-size:11px; font-weight:600;")
        lay.addWidget(self.measure_lbl, 3, 0, 1, 3)

        self.marks_lbl = QtWidgets.QLabel("")
        self.marks_lbl.setWordWrap(True)
        self.marks_lbl.setStyleSheet("color:#0a7d2c; font-size:11px; font-weight:600;")
        lay.addWidget(self.marks_lbl, 4, 0, 1, 3)

        self.send_mark_btn = QtWidgets.QPushButton("Send cut to last mark")
        self.send_mark_btn.setEnabled(False)
        self.send_mark_btn.setToolTip(
            "Move the cut plane of the opening your last mark sits on to that mark,\n"
            "instead of dragging its slider. Place a mark first (tick Measure, then\n"
            "click a spot on the vessel).")
        self.send_mark_btn.clicked.connect(self._send_cut_to_mark)
        lay.addWidget(self.send_mark_btn, 5, 0, 1, 2)

        self.clear_marks_btn = QtWidgets.QPushButton("Clear marks")
        self.clear_marks_btn.setEnabled(False)
        self.clear_marks_btn.clicked.connect(self._clear_marks)
        lay.addWidget(self.clear_marks_btn, 5, 2)

        # ── size and angulation over the region ──────────────────────────────
        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setStyleSheet("color:#e0e0e0;")
        lay.addWidget(line, 6, 0, 1, 3)

        self.size_btn = QtWidgets.QPushButton("Measure size and angulation")
        self.size_btn.setToolTip(
            "Equivalent-area diameter every 2 mm along the vessel, and how much the\n"
            "centerline turns over it.\n\n"
            "The angle is the arc a best-fit circle subtends: centerline length\n"
            "divided by that circle's radius. For a seal zone that is neck\n"
            "angulation, and unlike a hand-drawn pair of tangents it does not depend\n"
            "on where the tangents are drawn.\n\n"
            "Isolate a region first. A whole aorta is not one arc, and the fit error\n"
            "will say so.")
        self.size_btn.clicked.connect(self._measure_size)
        lay.addWidget(self.size_btn, 7, 0, 1, 3)

        self.size_lbl = QtWidgets.QLabel("")
        self.size_lbl.setWordWrap(True)
        self.size_lbl.setTextFormat(QtCore.Qt.RichText)
        self.size_lbl.setStyleSheet("color:#1c6b34; font-size:11px;")
        lay.addWidget(self.size_lbl, 8, 0, 1, 3)

        # ── the measure plane ────────────────────────────────────────────────
        #
        # A summary says "max 41.2 mm at 86 mm" and there is no way to see whether
        # 86 mm is where you meant. This is that station, drawn, and movable: the
        # measurement is only as good as its placement, and placement is the thing
        # the numbers cannot show you.
        self.mplane_box = QtWidgets.QWidget()
        mp = QtWidgets.QGridLayout(self.mplane_box)
        mp.setContentsMargins(0, 4, 0, 0)
        mp.setHorizontalSpacing(6)
        mp.setVerticalSpacing(2)

        lbl = QtWidgets.QLabel("Station")
        lbl.setStyleSheet("color:#666; font-size:11px;")
        mp.addWidget(lbl, 0, 0)
        self.mplane_slider, self.mplane_box_mm = _pair(0.0, 100.0, 0.0, 1.0, " mm")
        self.mplane_slider.valueChanged.connect(
            lambda v: self._measure_station(v / 10.0))
        self.mplane_box_mm.valueChanged.connect(self._measure_station)
        mp.addWidget(self.mplane_slider, 0, 1)
        mp.addWidget(self.mplane_box_mm, 0, 2)

        self.mplane_lbl = QtWidgets.QLabel("")
        self.mplane_lbl.setTextFormat(QtCore.Qt.RichText)
        self.mplane_lbl.setStyleSheet("color:#1c6b34; font-size:11px;")
        mp.addWidget(self.mplane_lbl, 1, 0, 1, 3)

        row = QtWidgets.QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        self.mplane_max = QtWidgets.QPushButton("Go to widest")
        self.mplane_max.setToolTip("Jump the plane to the widest station measured.")
        self.mplane_max.clicked.connect(self._measure_goto_max)
        row.addWidget(self.mplane_max, 1)
        # The disc sits inside the lumen and shows through the wall, which is the
        # point while placing it and in the way once it is placed -- particularly at
        # the arch, where it competes with the cut planes.
        self.mplane_show = QtWidgets.QCheckBox("Show plane")
        self.mplane_show.setChecked(True)
        self.mplane_show.setToolTip(
            "Hide the disc without losing the measurement. The station and its\n"
            "diameter stay, so it can be brought back where it was.")
        self.mplane_show.toggled.connect(self._measure_show_toggled)
        row.addWidget(self.mplane_show, 0)
        mp.addLayout(row, 2, 0, 1, 3)
        mp.setColumnStretch(1, 1)

        self.mplane_box.setVisible(False)      # nothing to place until it is measured
        lay.addWidget(self.mplane_box, 9, 0, 1, 3)

        lay.setColumnStretch(1, 1)
        return g

    # ── the measure plane: where the numbers came from ───────────────────────

    def _measure_station(self, s_mm):
        """Put the measure plane at *s_mm* along the route and section there."""
        if self.route is None or self.base is None:
            return
        s = float(np.clip(s_mm, self._m_range[0], self._m_range[1]))
        self._m_station = s
        for w, v in ((self.mplane_slider, int(round(s * 10))),
                     (self.mplane_box_mm, s)):
            w.blockSignals(True)
            w.setValue(v)
            w.blockSignals(False)

        nm, loc, _ = self.route.locate(s)
        pos, tan = self.tree.segments[nm].at_arclength(loc)
        sd = section_diameter(self.base, pos, tan, near_point=pos)
        dia = float(sd["ead_mm"]) if sd and sd.get("ead_mm") else 0.0
        if dia:
            txt = f"&empty; <b>{dia:.1f} mm</b> here"
            lo, hi = sd.get("min_mm"), sd.get("max_mm")
            if lo and hi and hi - lo > 0.05 * max(hi, 1.0):
                txt += f" &nbsp;·&nbsp; {lo:.1f}-{hi:.1f} across"
            self.mplane_lbl.setText(txt)
        else:
            self.mplane_lbl.setText("<i>no cross-section here</i>")
        self._draw_measure_plane(pos, tan, dia)

    def _draw_measure_plane(self, pos, tan, dia_mm=0.0):
        """Draw the disc, sized from the section it just measured.

        A little wider than the lumen so its edge is visible against the wall, and
        in its own colour so it is never mistaken for a cut plane -- this one
        removes nothing.
        """
        self._m_pose = (pos, tan, dia_mm)
        if not self.mplane_show.isChecked():
            return
        r = 0.8 * dia_mm if dia_mm > 0.1 else 8.0
        self.plotter.add_mesh(
            pv.Disc(center=pos, inner=0, outer=max(r, 2.0), normal=tan, c_res=64),
            color="#8e24aa", opacity=0.55, name="measure_plane", pickable=False)
        self.plotter.render()

    def _measure_show_toggled(self, on):
        """Hide or restore the disc, keeping the station it marks."""
        if on:
            if self._m_pose is not None:
                self._draw_measure_plane(*self._m_pose)
        else:
            try:
                self.plotter.remove_actor("measure_plane")
            except Exception:                                   # noqa: BLE001
                pass
            self.plotter.render()

    def _measure_goto_max(self):
        if getattr(self, "_m_widest", None) is not None:
            self._measure_station(self._m_widest)

    def _clear_measure_plane(self):
        try:
            self.plotter.remove_actor("measure_plane")
        except Exception:                                       # noqa: BLE001
            pass
        self.mplane_box.setVisible(False)

    def _measure_size(self):
        """Diameter along the vessel, and the angle the centerline turns through."""
        if self.tree is None or self.base is None or self.route is None:
            self._say("Load a mesh and a centerline first.", bad=True)
            return
        if self.region is None:
            s0, s1 = 0.0, float(self.route.length)
            where = "whole vessel"
        else:
            _a, _b = self.region
            s0 = float(_a) if _b is None else min(_a, _b)
            s1 = float(self.route.length) if _b is None else max(_a, _b)
            where = self._region_name or "region"

        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            prof = []
            for s in np.arange(s0, s1 + 1e-9, 2.0):
                nm, loc, _ = self.route.locate(float(s))
                pos, tan = self.tree.segments[nm].at_arclength(loc)
                sd = section_diameter(self.base, pos, tan, near_point=pos)
                if sd and sd.get("ead_mm"):
                    prof.append((float(s), float(sd["ead_mm"])))
            ang = centerline_arc_angle(self.tree, self.route, s0, s1)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

        if not prof:
            self.size_lbl.setText("")
            self._say("Could not section the mesh along this stretch.", bad=True)
            return

        ds = [d for _, d in prof]
        i = int(np.argmax(ds))

        # Hand the stretch that was measured to the plane, and open it at the widest
        # station -- the one the headline number is quoting, so the first thing shown
        # is where that number came from.
        self._m_range = (s0, s1)
        self._m_widest = float(prof[i][0])
        self.mplane_slider.setRange(int(round(s0 * 10)), int(round(s1 * 10)))
        self.mplane_box_mm.setRange(s0, s1)
        self.mplane_box.setVisible(True)
        self._measure_station(self._m_widest)

        bits = [f"<b>{where}</b>, {s1 - s0:.0f} mm",
                f"max &empty; <b>{ds[i]:.1f} mm</b> at {prof[i][0]:.0f} mm &nbsp;·&nbsp; "
                f"min {min(ds):.1f} &nbsp;·&nbsp; mean {float(np.mean(ds)):.1f}"]
        if ang is None:
            bits.append("too short to fit a circle")
        elif ang.get("straight"):
            bits.append("straight: no angulation")
        else:
            # The fit error decides whether the angle is worth quoting at all, so it
            # is shown next to it rather than buried.
            verdict = ("good fit" if ang["rms_mm"] <= 0.5 else
                       "fair fit" if ang["rms_mm"] <= 1.5 else
                       "<b>poor fit, treat the angle with caution</b>")
            bits.append(f"angle <b>{ang['angle_deg']:.0f}&deg;</b> "
                        f"(radius {ang['radius_mm']:.0f} mm) &nbsp;·&nbsp; "
                        f"fit error {ang['rms_mm']:.2f} mm, {verdict}")
            if getattr(self, "jagged_centerline", False) \
                    and not self.smooth_cb.isChecked():
                bits.append("centerline reads as jagged and is unsmoothed, which "
                            "inflates length and so the angle")
        self.size_lbl.setText("<br>".join(bits))
        self._say(f"Measured {where}: max &empty; {ds[i]:.1f} mm"
                  + (f", angle {ang['angle_deg']:.0f}&deg;"
                     if ang and not ang.get("straight") else ""))

    def _branch_bar(self):
        """A one-line toolbar that nudges ALL branch cuts together, by a set step."""
        self.branch_bar = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(self.branch_bar)
        h.setContentsMargins(2, 0, 2, 0)
        h.setSpacing(6)
        h.addWidget(QtWidgets.QLabel("Move all branches"))
        self.nudge_mm = QtWidgets.QDoubleSpinBox()
        self.nudge_mm.setRange(0.5, 50.0)
        self.nudge_mm.setDecimals(1)
        self.nudge_mm.setSingleStep(1.0)
        self.nudge_mm.setValue(5.0)
        self.nudge_mm.setSuffix(" mm")
        h.addWidget(self.nudge_mm)
        up = QtWidgets.QPushButton("▲ toward tips")
        up.setToolTip("Move every branch cut this far toward its tip (trim less).")
        up.clicked.connect(lambda: self._nudge_branches(-self.nudge_mm.value()))
        down = QtWidgets.QPushButton("▼ deeper")
        down.setToolTip("Move every branch cut this far away from its tip (trim more).")
        down.clicked.connect(lambda: self._nudge_branches(self.nudge_mm.value()))
        h.addWidget(up)
        h.addWidget(down)
        h.addStretch(1)

        wrap = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(wrap)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        v.addWidget(self.branch_bar)

        # Cut a vessel BEFORE it splits, e.g. a common iliac ahead of its internal /
        # external branches: one opening in place of the group.
        self.trim_btn = QtWidgets.QPushButton("Trim a branch group…")
        self.trim_btn.setMenu(QtWidgets.QMenu(self.trim_btn))
        self.trim_btn.setEnabled(False)
        self.trim_btn.setToolTip(
            "Cut a vessel BEFORE it splits, so the branches beyond it are removed\n"
            "together and leave ONE opening, a common iliac trimmed ahead of its\n"
            "internal/external split, or a celiac trunk ahead of its branches.")
        v.addWidget(self.trim_btn)

        # Anatomical names depend on which way round the volume is. A mirrored export
        # would swap left for right and nothing in the geometry could tell, so the
        # assumption is stated here (and recorded in the export) rather than hidden.
        self.anat_lbl = QtWidgets.QLabel("")
        self.anat_lbl.setWordWrap(True)
        self.anat_lbl.setStyleSheet("color:#8a6100; font-size:10px;")
        self.anat_lbl.setToolTip(
            f"Openings are named from the geometry, assuming {ANATOMY_FRAME}.\n"
            "Sides (left/right) follow from that; a mirrored export would swap them.\n"
            "Check one known vessel on a new case, and rename any opening by typing\n"
            "over its name. The assumption is recorded in the exported JSON.")
        v.addWidget(self.anat_lbl)

        # Isolate a stretch of the aorta -- an aortic neck, an Ishimaru landing zone,
        # the arch. Everything else is trimmed away; branches inside it stay.
        self.region_btn = QtWidgets.QPushButton("Isolate a region…")
        self.region_btn.setMenu(QtWidgets.QMenu(self.region_btn))
        self.region_btn.setEnabled(False)
        self.region_btn.setToolTip(
            "Keep only a stretch of the aorta, measured along its centerline:\n"
            "an aortic neck, an Ishimaru landing zone, the arch. Branches that\n"
            "leave from inside the region are kept and cut as usual.\n\n"
            "Only regions this anatomy can actually bound are offered.")
        v.addWidget(self.region_btn)

        # The nub-in preset. It was validated over 82 vessels in 30 thoracic cases and
        # then reachable only by running the cohort tool and reopening this window on
        # its output -- so the one cut with a measured answer behind it was the one cut
        # nobody could ask for. It sets every opening at once, which is why it is a
        # button here rather than an entry in the region menu: it is not a region.
        self.nub_btn = QtWidgets.QPushButton("Nub-in cut (preset)")
        self.nub_btn.setEnabled(False)
        self.nub_btn.setToolTip(
            "Cut every branch back to a short stub past the aortic wall, and square\n"
            "off the aorta's own ends. The stub length is 0.75 of each vessel's own\n"
            "radius -- the cohort's answer, not a preference: a branch leaves as a\n"
            "funnel and a plane through a funnel does not come back flat, so each\n"
            "vessel has a shortest nub it can take.\n\n"
            "A nub that will not cut cleanly is moved outward until it does, never\n"
            "inward, and every such rescue is reported. Sets all openings at once.")
        self.nub_btn.clicked.connect(self._apply_nub_preset)
        v.addWidget(self.nub_btn)

        # Commits a non-drag change to a region bound once the value settles, rather
        # than re-cutting on every key repeat. See _region_dragged.
        self._region_commit_timer = QtCore.QTimer(self)
        self._region_commit_timer.setSingleShot(True)
        self._region_commit_timer.setInterval(220)
        self._region_commit_timer.timeout.connect(self._region_commit)

        self.region_bar = QtWidgets.QWidget()
        rg = QtWidgets.QGridLayout(self.region_bar)
        rg.setContentsMargins(2, 0, 2, 0)
        rg.setHorizontalSpacing(6)
        rg.setVerticalSpacing(2)
        self.reg_s0_s, self.reg_s0_b = _pair(0.0, 1000.0, 0.0, 1.0, " mm")
        self.reg_s1_s, self.reg_s1_b = _pair(0.0, 1000.0, 0.0, 1.0, " mm")
        for s, b in ((self.reg_s0_s, self.reg_s0_b), (self.reg_s1_s, self.reg_s1_b)):
            s.valueChanged.connect(
                lambda vv, bb=b, ss=s: self._region_dragged(bb, vv, ss))
            s.sliderReleased.connect(self._region_commit)
            b.valueChanged.connect(lambda vv, ss=s: self._region_typed(ss, vv))
        rg.addWidget(QtWidgets.QLabel("Start"), 0, 0)
        rg.addWidget(self.reg_s0_s, 0, 1)
        rg.addWidget(self.reg_s0_b, 0, 2)
        rg.addWidget(QtWidgets.QLabel("End"), 1, 0)
        rg.addWidget(self.reg_s1_s, 1, 1)
        rg.addWidget(self.reg_s1_b, 1, 2)
        self.region_lbl = QtWidgets.QLabel("")
        self.region_lbl.setWordWrap(True)
        self.region_lbl.setStyleSheet("color:#0a5d2c; font-size:11px; font-weight:600;")
        rg.addWidget(self.region_lbl, 2, 0, 1, 2)
        self.region_clear = QtWidgets.QPushButton("Whole vessel")
        self.region_clear.setToolTip("Drop the region and go back to cutting the "
                                     "whole vessel.")
        self.region_clear.clicked.connect(lambda: self._set_region(None))
        rg.addWidget(self.region_clear, 2, 2)
        rg.setColumnStretch(1, 1)
        self.region_bar.setVisible(False)
        v.addWidget(self.region_bar)
        return wrap

    # ---------------- region isolation ----------------

    def _refresh_region_menu(self):
        if not hasattr(self, "region_btn"):
            return
        menu = self.region_btn.menu()
        menu.clear()
        presets = []
        if self.tree is not None and self.base is not None:
            try:
                # planar=True: this menu drives a planar cut, and region_presets no
                # longer filters on planar feasibility by default. See
                # vascular_api.CutSession.presets for the measurement.
                presets = region_presets(self.tree, self.base, self.anat_names,
                                         self.anat_vessels, planar=True)
            except Exception:
                presets = []
        # The way OUT belongs in the same menu as the way in: a user who turned a
        # region on comes back here to turn it off, and finding only presets reads
        # as "there is no way back".
        if self.region is not None:
            back = menu.addAction("← Back to the whole vessel")
            back.setToolTip("Drop the region and cut the whole vessel again.")
            back.triggered.connect(lambda: self._set_region(None))
            menu.addSeparator()
        for p in presets:
            # s1 None is an open distal end: one plane, everything below it kept,
            # so there is no length to quote.
            span = ("to the ends" if p["s1"] is None
                    else f"{p['s1'] - p['s0']:.0f} mm")
            act = menu.addAction(f"{p['name']}   ({span})")
            act.setToolTip(p["why"])
            act.triggered.connect(lambda _=False, pp=p: self._apply_preset(pp))
        # The nub preset needs a whole model, not a stretch of one: it cuts the aorta's
        # own ends, which an isolated region has already replaced with its bounds.
        if hasattr(self, "nub_btn"):
            self.nub_btn.setEnabled(self.tree is not None and self.base is not None
                                    and self.region is None)
        self.region_btn.setEnabled(bool(presets) or self.region is not None)
        self.region_btn.setText(
            f"Region: {self._region_name} ▾" if self.region is not None
            else "Isolate a region…")

    def _apply_preset(self, p):
        """Take a preset, and if it carries alternatives, the first that cuts.

        The aortic arch's distal end is a convention (last great vessel + 20 mm)
        that often lands in the hairpin, where no planar cut exists. Rather than
        hand the user a refusal, walk the bound out to the first one that
        validates, and say how far it had to go, so the region is never quietly
        redefined behind their back.
        """
        self._region_manual = False        # a preset picked these stations, not you
        options = p.get("s1_options") or [p["s1"]]
        for i, s1 in enumerate(options):
            self._set_region((p["s0"], s1), p["name"])
            if self.problem is None:
                if i:
                    self._say(f"{self.status.text()}  ,   {p['name']} extended to "
                              f"{s1 - p['s1'] + 20:.0f} mm past the last great vessel "
                              f"({s1 - p['s1']:+.0f} mm) to reach a clean cut.")
                return
        # nothing validated: leave it on the conventional bound with its own error,
        # which is the honest answer rather than the least-bad arbitrary span
        self._set_region((p["s0"], p["s1"]), p["name"])

    def _set_region(self, span, name=None):
        """Turn a region on (a station pair) or off (None)."""
        # Cancel any debounced bound-change still in flight. Without this, nudging a
        # bound and immediately leaving lets the timer fire afterwards and put the
        # region back -- with its bar hidden, which is the original trap made worse.
        self._region_commit_timer.stop()
        # A reach found for one region says nothing about the next: it was fitted to
        # one pair of planes at one place on the aorta. Cleared on EVERY change of
        # region, not only on leaving one, so moving from the neck to the infrarenal
        # segment does not inherit a scale chosen for a different pair of planes.
        self._region_reach = None
        self._region_moved_mm = 0.0
        if span is None:
            self._region_manual = False    # next region starts fresh, not "yours"
            # ...and so do its boundary planes. Carrying a tilt over into the NEXT
            # region would apply it to a plane somewhere else entirely.
            self._bound_planes = {}
        # span[1] None means an open distal end: one plane, everything below it kept.
        self.region = (None if span is None else
                       (float(span[0]),
                        None if span[1] is None else float(span[1])))
        self._region_name = name or ("custom" if self.region is not None else "")
        self.region_bar.setVisible(self.region is not None)
        self._sync_region_bar()
        self._build_rows()
        self._on_change()

    def _sync_region_bar(self):
        """Put the Start/End controls where self.region actually is.

        Factored out because the region can now move without the user touching it:
        a boundary with no cuttable cross-section is nudged along the aorta, and the
        controls have to follow or they would show a station the cut did not use.
        """
        if self.region is None or self.route is None:
            return
        L = self.route.length
        self._syncing_region = True
        pairs = [(self.reg_s0_s, self.reg_s0_b, self.region[0])]
        if self.region[1] is not None:
            pairs.append((self.reg_s1_s, self.reg_s1_b, self.region[1]))
        for s, b, val in pairs:
            b.setRange(0.0, L)
            s.setRange(0, int(round(L * 10)))
            b.setValue(val)
            s.setValue(int(round(val * 10)))
        self._syncing_region = False
        # There is no end station to drag when the region has no end.
        open_distal = self.region[1] is None
        for wdg in (self.reg_s1_s, self.reg_s1_b):
            wdg.setEnabled(not open_distal)
            wdg.setToolTip("This region has no end: everything below the start "
                           "is kept." if open_distal else "")

    def _span_from_controls(self):
        """The Start/End controls read back as a span, keeping an open end OPEN.

        The End control is DISABLED whenever the region has no end (_sync_region_bar
        above), so the number inside it is not the user's: it is whatever was left
        there -- 0.0 in a fresh session, from the widget's own construction, or the
        End of whichever two-plane region was used before.

        Reading it unconditionally is how "Infrarenal aorta with iliacs" lost its
        open end the moment its Start slider moved. That preset is the only region
        that ships s1 None, and deliberately: region_plan clamps a distal plane a
        millimetre short of the route end, and an abdominal aorta's iliac ostia sit
        AT the bifurcation, so a two-plane version drops them (branch_cut.py, the
        comment above "Infrarenal aorta with iliacs"). Once s1 was a number again the
        region either inverted (stale 0.0, so min/max in _region_terminals made the
        span run from the top of the aorta down to the drag point) or clamped short
        of the bifurcation -- and either way the iliacs fell outside it.

        The user could not undo it, because the control now driving that plane was
        still greyed out and still saying the region had no end. Reported from the
        cohort as "the distal cutting plane moves up and you cannot change it back".
        """
        s0 = float(self.reg_s0_b.value())
        # An open end can only be left by choosing a different region: there is no
        # control to fill it in, so a region that has one still has one.
        if self.region is not None and self.region[1] is None:
            return (s0, None)
        return (s0, float(self.reg_s1_b.value()))

    def _region_dragged(self, box, v, slider=None):
        if getattr(self, "_syncing_region", False):
            return
        self._syncing_region = True
        box.setValue(v / 10.0)
        self._syncing_region = False
        self.region = self._span_from_controls()
        # Rebuild the bounds before drawing them. _draw_planes reads cut_list(), whose
        # region bounds come from self._region_bounds -- and that list is only rebuilt
        # by _build_rows() on COMMIT. So dragging a region slider moved self.region and
        # then drew the planes from the old stations: the discs sat still while the
        # numbers changed, which reads as the control being dead. Replanning is cheap
        # (no mesh is cut), and the radius it needs is cached after the first call.
        self._region_manual = True     # a dragged boundary is the user's, do not nudge
        try:
            self._region_bounds, _inside, self._region_notes = self._region_terminals()
        except Exception:                                       # noqa: BLE001
            pass                       # a momentarily impossible span: keep the last
        self._draw_planes()
        self.plotter.render()
        # Arrow keys / wheel / groove clicks never emit sliderReleased, so without
        # this the bounds move (and the exported metadata with them) while the mesh
        # stays cut where it was, under an unchanged green status. Debounced.
        if slider is not None and not slider.isSliderDown():
            self._region_commit_timer.start()

    def _region_typed(self, slider, v):
        if getattr(self, "_syncing_region", False):
            return
        self._syncing_region = True
        slider.setValue(int(round(v * 10)))
        self._syncing_region = False
        self._region_commit()

    def _region_commit(self):
        # Reaching here means a slider or spin box moved: the bounds are the user's
        # from now on, so stop relocating them to clear ostia.
        self._region_manual = True
        # A dismissed region stays dismissed. _on_change runs a re-cut that takes far
        # longer than the 220 ms debounce, so a timer armed just before the user left
        # can still fire during it; re-reading the boxes here would resurrect the
        # region behind a hidden bar. Belt as well as braces with the stop() in
        # _set_region, because this is the failure that traps the user.
        if self.region is None:
            return
        self.region = self._span_from_controls()
        self._build_rows()
        self._on_change()

    def _region_terminals(self):
        """(bounds, inside, notes) for the active region, or ([], [], []) if none."""
        if self.region is None or self.tree is None:
            return [], [], []
        # s1 None is an open-ended region: one plane, everything below it kept.
        # sorted() cannot order None against a float, so unpack rather than sort.
        a, b = self.region
        s0, s1 = (float(a), None) if b is None else (min(a, b), max(a, b))
        # Pass the mesh so ostium clearances are measured, not guessed from a 2 mm
        # fallback, otherwise the region's bounds depend on what ran before it.
        # Once the user has moved a slider the boundary is THEIRS: clearing ostia
        # then relocates it behind their back, which is what made the slider snap.
        # A preset keeps the automatic clearance, since they picked no station.
        bounds, inside, _removed, notes = region_plan(
            self.tree, s0, s1, self.route, mesh=self.base,
            clear_ostia=not getattr(self, "_region_manual", False))
        # region_plan builds FRESH Terminal objects on every call, and it is called
        # again on every nudge of a region slider -- so a tilt or plane size the user
        # set would be discarded by the next drag, with no sign it had gone. Carry
        # those settings across. Keyed by the bound's role rather than its label:
        # the label carries the segment the plane happens to land on, which changes
        # under the user as they move the boundary along the aorta.
        saved = getattr(self, "_bound_planes", {})
        for t in bounds:
            for attr, val in saved.get(t.name, {}).items():
                setattr(t, attr, val)
        # ...and a nudge found for a boundary that was undercutting a branch's
        # ostium, carried the same way and for the same reason.
        nudge = getattr(self, "_bound_nudge", None) or {}
        for t in bounds:
            if t.name in nudge:
                t._region_offset = float(nudge[t.name])
        # A reach the fitter settled on for THIS region, reapplied to the fresh
        # objects. After the saved settings, so a plane size the user typed wins.
        self._apply_region_reach(bounds)
        return bounds, inside, notes

    def _apply_region_reach(self, bounds):
        """Scale the boundary planes by whatever reach was found to cut here."""
        scale = getattr(self, "_region_reach", None)
        if not scale:
            return
        for t in bounds:
            if not getattr(t, "band_mm", None):
                t.band_mm = float(band_width(t)) * float(scale)

    BOUND_PLANE_ATTRS = ("tilt_a_deg", "tilt_b_deg", "shift_a_mm", "shift_b_mm",
                         "band_mm")

    def _on_bound_change(self):
        """A region boundary's plane was tilted or resized: remember it, then re-cut."""
        self._bound_planes = getattr(self, "_bound_planes", {})
        for row in getattr(self, "bound_rows", []):
            t = row.terminal
            self._bound_planes[t.name] = {a: getattr(t, a, None)
                                          for a in self.BOUND_PLANE_ATTRS}
        self._on_change()

    def _apply_nub_preset(self):
        """Set every opening to the nub-in preset, and say what it had to do.

        The notes are not decoration. A rescued nub is longer than the rule asked
        for, and a trunk end that fell back to a perpendicular plane is not the
        horizontal end the paper wants -- both are cases for a person to look at, and
        neither is visible from the geometry alone once the cut is made.
        """
        if self.tree is None or self.base is None:
            self._say("Load a mesh and a centerline first.", bad=True)
            return
        try:
            import nub_session
        except Exception as e:                                      # noqa: BLE001
            self._say(f"The nub preset is unavailable: {e}", bad=True)
            return
        self.nub_btn.setEnabled(False)
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            notes = nub_session.apply_here(self)
        except Exception as e:                                      # noqa: BLE001
            self._say(f"The nub preset failed: {e}", bad=True)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
            self._refresh_region_menu()          # re-enables the button
        head = next((n for n in notes if n.startswith("nub-in cut")), "Nub preset applied.")
        rest = [n for n in notes if not n.startswith("nub-in cut")]
        self._say(head + (f"  {len(rest)} note(s) — see below." if rest else ""))
        if rest:
            self.anat_lbl.setText("Nub preset: " + "; ".join(rest[:4])
                                  + (f"; +{len(rest) - 4} more" if len(rest) > 4 else ""))

    def _on_renamed(self, row):
        """One opening's artery was declared, cleared, or just relabelled.

        Naming changes no geometry, so this must NOT re-cut and must NOT rebuild the
        rows: the row that emitted this is the one the user is still typing in, and
        replacing it mid-edit loses the keystroke and the focus. Only the things that
        FOLLOW from a name are rebuilt -- the region presets, the landmark stations
        they are placed against, and the caption saying how many openings are named.

        Declaring an artery can silently rename a DIFFERENT row, because a hand-named
        renal displaces an automatically named one (see _apply_manual). Those rows
        have to be re-read from the terminals or the panel would keep showing a name
        that no longer applies.
        """
        if self.tree is None:
            return
        set_manual_name(self.tree, row.terminal.label, row.declared_artery())
        try:
            self.anat_names, self.anat_vessels, self.anat_notes = name_anatomy(
                self.tree, self.base)
        except Exception as e:                                      # noqa: BLE001
            self._say(f"Could not re-read the anatomy: {e}", bad=True)
            return
        for t in self.tree.terminals:
            if t.label in self.anat_names:
                t.name = self.anat_names[t.label]
        self.route = aorta_route(self.tree)
        for r in self.all_rows():
            if r is row:
                continue
            want = getattr(r.terminal, "name", "")
            if want and r.name.currentText().strip() != want:
                r.name.blockSignals(True)
                r.name.setCurrentText(want)
                r.name.blockSignals(False)
        self._apply_colour_labels()
        self._refresh_region_menu()
        self._update_anat_label()
        what = row.declared_artery()
        self._say(f"{row.terminal.label} declared as {what}." if what else
                  f"{row.terminal.label} relabelled. Pick an artery from the list to "
                  f"drive the region presets with it.")

    def _update_anat_label(self):
        # Say which anatomy the names came from, and whether that was inferred or
        # declared. Reading "Ascending aorta" on an abdominal case is only puzzling
        # until you can see the tool decided this was a thoracic aorta.
        # A mesh that does not match the centerline makes every name below it
        # meaningless, so it pre-empts the whole caption rather than appending to it.
        # The cut refuses on the same evidence, but only once it is attempted, and by
        # then the wrong anatomy has been on screen for as long as the user has been
        # looking at it.
        if self.tree is not None and not getattr(self.tree, "mesh_matches", True):
            cov = float(getattr(self.tree, "coverage", 1.0))
            outside = list(getattr(self.tree, "radius_outside", None) or [])
            why = (f"only {100 * cov:.0f}% of the centerline lies within this mesh"
                   if cov < MIN_CENTERLINE_COVERAGE else
                   f"{len(outside)} of {len(self.tree.terminals)} centerlines run "
                   f"outside this mesh")
            self.anat_lbl.setText(
                f"⚠ This mesh and this centerline do not look like the same "
                f"export — {why}. Nothing below has been read from the anatomy. "
                f"The usual cause is a trunk-only mesh paired with a branched "
                f"centerline.")
            return

        fam = ""
        if self.tree is not None:
            n_trunk = sum(1 for t in self.tree.terminals
                          if getattr(t, "is_trunk", False))
            kind = "abdominal" if n_trunk == 1 else "thoracic"
            how = ("you set this"
                   if self.family.currentData() else "read from the geometry")
            fam = f"Treated as {kind}, {how}. "
            # A NEAR call is the one worth checking, and the user cannot tell a
            # near call from a certain one unless the margin is on screen. The
            # rule reads a trunk that ENDS BY DIVIDING as abdominal, and the cases
            # it still gets wrong are the ones close to the line -- so say so
            # there, and stay quiet when the answer is not in doubt.
            split = getattr(self.tree, "outflow_split", None)
            if split is not None and not self.family.currentData():
                if abs(float(split) - PAIR_AREA_FRACTION) < 0.08:
                    fam = (f"Treated as {kind}, and it was a CLOSE call "
                           f"({float(split):.2f} against {PAIR_AREA_FRACTION:.2f}) "
                           f"-- worth confirming. ")
        n = len(self.anat_names)
        # A hand-declared name is not a guess, so it does not carry the left/right
        # caveat and the count of them is worth stating separately -- it is the
        # difference between "the tool read this" and "you told it".
        manual = len(getattr(self.tree, "manual_names", None) or {}) if self.tree else 0
        by_hand = f" ({manual} by hand)" if manual else ""
        if not n:
            self.anat_lbl.setText(
                (fam + "No opening could be named unambiguously — name them yourself "
                       "from the list in each opening's name box, and the region "
                       "presets follow.") if fam else "")
            return
        # Left/Right is decided by ONE scalar for the whole case, so when the witness
        # that sets it is missing, EVERY side name is a guess from the frame -- and
        # until now the caption said the same thing either way. The note recording it
        # reached only the export JSON, which is not where a person looks.
        # 11 of 59 abdominal cases in the last cohort run carry this note.
        sides_assumed = any(SIDES_ASSUMED_NOTE in str(nt) for nt in self.anat_notes)
        if sides_assumed:
            self.anat_lbl.setText(
                f"⚠ {fam}{n} opening(s) named anatomically{by_hand}. LEFT AND RIGHT "
                f"ARE NOT CONFIRMED on this case: nothing in the anatomy witnessed "
                f"which side is which, so the sides follow only from {ANATOMY_FRAME}. "
                f"A mirrored export would swap every one of them and the geometry "
                f"could not tell. CHECK THEM before you rely on a side. Wrong "
                f"anatomy? Set it under Centerline. Any opening can be named by hand "
                f"from its own name box.")
            return
        self.anat_lbl.setText(
            f"⚠ {fam}{n} opening(s) named anatomically{by_hand}, assuming "
            f"{ANATOMY_FRAME}, confirm left/right on a new case. Wrong anatomy? Set "
            f"it under Centerline. Any opening can be named by hand from its own "
            f"name box.")

    def _output_group(self):
        g = QtWidgets.QWidget()          # titled by its Section header
        lay = QtWidgets.QGridLayout(g)
        lay.setContentsMargins(6, 4, 2, 4)
        lay.setHorizontalSpacing(6)

        self.out_edit = QtWidgets.QLineEdit(str(self.out_dir))
        self.out_edit.setPlaceholderText("folder for exported meshes")
        self.out_edit.setToolTip("Remembered between sessions.")
        self.out_edit.editingFinished.connect(self._out_typed)
        out_btn = QtWidgets.QPushButton("Browse…")
        out_btn.clicked.connect(self._pick_out)

        self.name_edit = QtWidgets.QLineEdit()
        self.name_edit.setPlaceholderText("file name (no extension)")
        self.name_edit.textChanged.connect(self._refresh_preview)

        self.preview_lbl = QtWidgets.QLabel("")
        self.preview_lbl.setStyleSheet("color:#777; font-size:11px;")

        self.fmt = QtWidgets.QComboBox()
        self.fmt.addItems(["stl", "ply", "obj"])
        self.fmt.setCurrentText(self.settings.value("format", "stl"))
        self.fmt.currentTextChanged.connect(self._fmt_changed)

        self.clean_cb = QtWidgets.QCheckBox("Sliver-free rim")
        self.clean_cb.setChecked(self.settings.value("snap", "false") == "true")
        self.clean_cb.setToolTip(
            "OFF (default): the cut is bit-identical to your mesh. A few hair-thin\n"
            "triangles are left on the rims.\n\n"
            "ON: vertices within ~5% of a mesh edge of a cut plane are nudged onto it\n"
            "before clipping, so no sliver is created. Openings stay exactly planar.\n"
            "The tolerance scales with YOUR mesh's density, not a fixed distance.\n"
            "Nothing is re-sampled or smoothed, this is not a remesh.")
        self.clean_cb.toggled.connect(self._snap_toggled)

        lay.addWidget(QtWidgets.QLabel("Folder"), 0, 0)
        lay.addWidget(self.out_edit, 0, 1)
        lay.addWidget(out_btn, 0, 2)
        lay.addWidget(QtWidgets.QLabel("File name"), 1, 0)
        lay.addWidget(self.name_edit, 1, 1)
        lay.addWidget(self.fmt, 1, 2)
        lay.addWidget(self.preview_lbl, 2, 1)
        lay.addWidget(self.clean_cb, 2, 2)
        lay.setColumnStretch(1, 1)
        return g

    # ---------------- file pickers ----------------

    def _pick_csv(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Centerline", str(Path.cwd()),
            "Centerlines (*.csv *.tsv *.txt *.dat *.vtp *.vtk *.vtu *.json);;"
            "Tables (*.csv *.tsv *.txt *.dat);;"
            "Polylines (*.vtp *.vtk *.vtu);;"
            "Markups (*.json);;All files (*)")
        if p:
            self.csv_path = Path(p)
            self._reset_family()      # a declared anatomy belongs to the old case
            self.csv_label.setText(p)

    def _clear_csv(self):
        self.csv_path = None
        self.csv_label.setText("")

    def _pick_mat(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Surface mesh", str(Path.cwd()),
            "Meshes (*.mat *.stl *.ply *.obj);;All files (*)")
        if p:
            self.mat_path = Path(p)
            self.mat_label.setText(p)
            self._reset_family()      # a declared anatomy belongs to the old case

    def _pick_out(self):
        p = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Folder for exported meshes", str(self.out_dir))
        if p:
            self.out_dir = Path(p)
            self.out_edit.setText(p)
            self.settings.setValue("out_dir", p)
            self._refresh_preview()

    def _out_typed(self):
        t = self.out_edit.text().strip()
        if t:
            self.out_dir = Path(t)
            self.settings.setValue("out_dir", t)
            self._refresh_preview()

    @staticmethod
    def _typed_path(text):
        """A pasted path, made usable. Returns None for nothing typed.

        Explorer's "Copy as path" wraps the whole thing in double quotes, and a path
        dragged from a terminal often arrives with a trailing space. Neither exists
        on disk, so without this a perfectly correct path is rejected as missing --
        the one failure that makes a person go back to the file dialog for good.
        """
        t = str(text or "").strip().strip('"').strip("'").strip()
        return Path(t) if t else None

    def _mat_typed(self):
        p = self._typed_path(self.mat_label.text())
        if p == self.mat_path:
            return                                # focus moved, nothing was edited
        self.mat_path = p
        self.mat_label.setText(str(p) if p else "")
        self._reset_family()                      # a declared anatomy is the old case's
        if p is not None and not p.exists():
            self._say(f"{p} does not exist. Check the path, or use Browse….", bad=True)
        elif p is not None:
            self._say(f"Mesh set to {p.name}. Press Load to open it.")

    def _csv_typed(self):
        p = self._typed_path(self.csv_label.text())
        if p == self.csv_path:
            return
        self.csv_path = p
        self.csv_label.setText(str(p) if p else "")
        self._reset_family()
        if p is None:
            self._say("Centerline cleared, the mesh will be shown without cutting.")
        elif not p.exists():
            self._say(f"{p} does not exist. Check the path, or use Browse….", bad=True)
        else:
            self._say(f"Centerline set to {p.name}. Press Load to open it.")

    def _fmt_changed(self, v):
        self.settings.setValue("format", v)
        self._refresh_preview()

    def export_stem(self):
        typed = self.name_edit.text().strip()
        if typed:
            return Path(typed).stem
        return f"{self.mat_path.stem}_open" if self.mat_path else "open_mesh"

    def _refresh_preview(self):
        stem = self.export_stem()
        self.preview_lbl.setText(
            f"→ {stem}.{self.fmt.currentText()}  +  {stem}_cut_planes.json")

    # ---------------- loading ----------------

    def load_model(self, keep_smooth_choice=False):
        if not self.mat_path:
            self._say("Pick a mesh first.", bad=True)
            return
        try:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            self.base = (load_mat_mesh(self.mat_path)
                         if self.mat_path.suffix.lower() == ".mat"
                         else trimesh.load(str(self.mat_path), process=False))
            self.tree = None
            if self.csv_path:
                # Load raw FIRST and measure it as supplied. Only then decide
                # whether to offer smoothing, and only suggest, never impose.
                probe = load_tree(self.csv_path, smooth_mm=0.0)
                if not keep_smooth_choice:
                    self.smooth_cb.blockSignals(True)
                    self.smooth_cb.setChecked(probe.jagged)
                    self.smooth_cb.blockSignals(False)
                # Kept, because an unsmoothed staircase inflates centerline length
                # and therefore any angle derived from it.
                self.jagged_centerline = bool(probe.jagged)
                self._show_quality(probe)

                smooth = self.smooth_mm.value() if self.smooth_cb.isChecked() else 0.0
                self.tree = (probe if smooth == 0.0
                             else load_tree(self.csv_path, smooth_mm=smooth))
                classify_terminals(self.tree, self.base,
                                   n_trunk_ends=self.family.currentData())
                # Name each opening after the artery it actually is, where the
                # geometry says so unambiguously; anything else keeps Branch N.
                self.anat_names, self.anat_vessels, self.anat_notes = name_anatomy(
                    self.tree, self.base)
                for t in self.tree.terminals:
                    if t.label in self.anat_names:
                        t.name = self.anat_names[t.label]
                self.route = aorta_route(self.tree)
        except Exception as e:
            # Do not leave the PREVIOUS model's openings on screen naming vessels that
            # are no longer loaded, clear the panel with the failure.
            self._say(f"Could not load: {e}", bad=True)
            self.tree = None
            self.trims = []
            self.region, self.route = None, None
            self._region_name = ""
            self.region_bar.setVisible(False)
            self._offset_memory = {}
            self._build_rows()
            self._set_enabled(False)
            # The load failed, so Browse/Load are the only useful controls left. A
            # previous success may have persisted Input folded, which would leave the
            # user staring at a failure with no visible way to pick another file.
            self.sec_input.set_open(True)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

        # Tuning the smoothing re-enters here for the SAME model, so neither the
        # export name nor the fold state may be reset on that path: it would discard
        # a typed filename and collapse the spin box out from under the cursor.
        if not keep_smooth_choice:
            self.name_edit.setText(f"{self.mat_path.stem}_open")
        self._refresh_preview()
        # Setup is done for this model: fold it away and give the openings the room.
        # View & measure stays as the user left it, it is used while cutting. With no
        # centerline the panel goes on to say "use Browse… next to Centerline, then
        # Load", so folding Input away there would hide the controls it names.
        if not keep_smooth_choice:
            if self.tree is not None:
                for s in (self.sec_input, self.sec_cl, self.sec_out):
                    s.set_open(False)
            else:
                # Mesh only: the openings panel is about to say "use Browse… next to
                # Centerline, then Load", and those live in Input, so a fold left
                # over from an earlier load would hide the controls it names.
                self.sec_input.set_open(True)
        # A load can change what a folded header is asserting (smoothing is auto-set
        # from the centerline's own jaggedness), so re-derive every summary whether or
        # not the fold state moved.
        for s in (self.sec_input, self.sec_cl, self.sec_view, self.sec_out):
            s.refresh()
        self.trims = []                  # a new model: old trims name old openings
        self.region = None
        # A measurement belongs to the geometry it was taken on, so a new model
        # retires the plane rather than leaving it pointing at a station that means
        # something else now.
        self.size_lbl.setText("")
        self._m_station = self._m_widest = self._m_pose = None
        self._clear_measure_plane()
        self._region_name = ""
        self.region_bar.setVisible(False)
        self._offset_memory = {}
        self._build_rows()
        self._set_enabled(self.tree is not None)

        if self.tree is None:
            self.cut = self.base
            self.matched = {}
            self.problem = None
            self._say(f"Mesh loaded, {len(self.base.faces):,} faces, "
                      f"watertight={self.base.is_watertight}. "
                      f"Add a centerline to cut it.")
        else:
            self.recut()
        self.draw(reset_camera=True)

    def _build_rows(self):
        # Adding or undoing a trim rebuilds every row, so carry the cut depths across
        # rather than dropping the user's work back to the default.
        self._offset_memory = getattr(self, "_offset_memory", {})
        self._offset_memory.update({r.terminal.label: r.value() for r in self.rows})
        self._rgb_cache_key = None       # the palette is about to be reassigned
        while self.openings_lay.count():
            item = self.openings_lay.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)      # deleteLater only SCHEDULES removal, and the
                w.deleteLater()        # old widget keeps painting until it happens
        self.rows = []
        self.bound_rows = []
        self.sel = None
        self.opening_colour = {}

        if self.tree is None:
            lbl = QtWidgets.QLabel(
                "No centerline loaded, viewing the mesh only.\n"
                "Use Browse… next to Centerline, then Load, to cut it.")
            lbl.setStyleSheet("color:#777; padding:12px;")
            lbl.setWordWrap(True)
            self.openings_lay.addWidget(lbl)
            self.openings_lay.addStretch(1)
            self.openings_box.setTitle("Openings")
            self._refresh_trim_menu()    # nothing is trimmable without a centerline
            self._refresh_region_menu()  # ...nor regionable: drop the old model's
            return                       # presets and the "Region: X" caption

        # With a region active its two bounds are driven by the region sliders, not by
        # rows of their own -- only the branches inside it get rows.
        self._region_bounds, self._region_notes = [], []
        rowed = self.active_terminals()
        if self.region is not None:
            bounds, inside, notes = self._region_terminals()
            self._region_bounds, self._region_notes = bounds, notes
            rowed = self.region_rows()
            # A region silently empties this list when no branch falls inside it. Say
            # so, and keep the way back where the user is already looking.
            self.openings_lay.addWidget(self._region_banner())
            # The region's own two cuts, so their planes can be tilted and resized.
            # They are NOT in self.rows: that list is the offset-bearing cut list,
            # and cut_list() adds the bounds separately from _region_bounds.
            for t in self._region_bounds:
                row = TerminalRow(t, self._row_limit(t), band_width(t), bound=True)
                row.changed.connect(self._on_bound_change)
                row.preview.connect(self._on_preview)
                row.selected.connect(self._on_select)
                self.openings_lay.addWidget(row)
                self.bound_rows.append(row)
                line = QtWidgets.QFrame()
                line.setFrameShape(QtWidgets.QFrame.HLine)
                line.setStyleSheet("color:#e0e0e0;")
                self.openings_lay.addWidget(line)
            if not rowed:
                empty = QtWidgets.QLabel(
                    "No branches leave the vessel inside this region, only its two "
                    "ends are cut. Widen the region, or go back to the whole vessel."
                    if self._region_bounds else
                    "This region collapses to nothing once its branch ostia are "
                    "cleared, so nothing is being cut. Widen it, or go back to the "
                    "whole vessel.")
                empty.setWordWrap(True)
                empty.setStyleSheet("color:#777; padding:8px 4px;")
                self.openings_lay.addWidget(empty)

        choices = anatomy_choices(self.tree)
        for t in rowed:
            row = TerminalRow(t, self._row_limit(t), band_width(t), choices=choices)
            if t.label in self._offset_memory:
                row.set_value(self._offset_memory[t.label])   # clamps to the new limit
            row.changed.connect(self._on_change)
            row.preview.connect(self._on_preview)
            row.selected.connect(self._on_select)
            row.pullback.connect(self._on_pullback)
            row.renamed.connect(self._on_renamed)
            if getattr(t, "is_interior", False):
                row.set_trim(lambda _=None, term=t: self._remove_trim(term))
            self.openings_lay.addWidget(row)
            self.rows.append(row)
            line = QtWidgets.QFrame()
            line.setFrameShape(QtWidgets.QFrame.HLine)
            line.setStyleSheet("color:#e0e0e0;")
            self.openings_lay.addWidget(line)
        self.openings_lay.addStretch(1)

        # Only branches get a distinct colour; the aorta trunk stays red. Index by
        # branch order so Branch 1/2/3 map to stable, well-separated hues.
        self.opening_colour = {}
        bi = 0
        for r in self.rows:
            if not r.terminal.is_trunk:
                self.opening_colour[r.terminal.label] = OPENING_PALETTE[bi % len(OPENING_PALETTE)]
                bi += 1
        self._apply_colour_labels()

        n_br = sum(1 for r in self.rows if not r.terminal.is_trunk)
        kind = "unbranched trunk" if n_br == 0 else f"trunk + {n_br} branch(es)"
        trimmed = (f", {len(self.trims)} trimmed group(s)" if self.trims else "")
        if self.region is not None:
            n = len(self.rows) + len(self._region_bounds)
            self.openings_box.setTitle(
                f"Openings, {n} (isolated region, {len(self.rows)} branch(es) inside)")
        else:
            self.openings_box.setTitle(
                f"Openings, {len(self.rows)} ({kind}{trimmed})")
        self._refresh_trim_menu()
        self._refresh_region_menu()
        self._update_anat_label()
        self._update_region_label()

    def _region_banner(self):
        """A loud 'you are looking at a region, here is the way out' strip.

        The old exit was a short button in the corner of the region's slider grid,
        below a wrapped label -- easy to miss, and the region menu offered no way
        back at all. This puts the escape at the top of the list it emptied.
        """
        bar = QtWidgets.QFrame()
        bar.setStyleSheet(
            "QFrame { background:#e8f4ec; border:1px solid #7fb894; border-radius:4px; }")
        h = QtWidgets.QHBoxLayout(bar)
        h.setContentsMargins(8, 5, 5, 5)
        h.setSpacing(6)
        name = self._region_name or "custom"
        lbl = QtWidgets.QLabel(f"Isolated region, {name}")
        lbl.setStyleSheet("color:#0a5d2c; font-weight:600; border:none;")
        lbl.setToolTip("Only this stretch of the aorta is kept. Branches that leave "
                       "from inside it are cut as usual.")
        h.addWidget(lbl)
        h.addStretch(1)
        back = QtWidgets.QPushButton("✕  Whole vessel")
        back.setToolTip("Drop the region and go back to cutting the whole vessel.")
        back.setStyleSheet("border:1px solid #7fb894; border-radius:3px; padding:2px 8px;")
        back.clicked.connect(lambda: self._set_region(None))
        h.addWidget(back)
        return bar

    def _update_region_label(self):
        if self.region is None or self.route is None:
            self.region_lbl.setText("")
            return
        _a, _b = self.region
        s0, s1 = ((float(_a), None) if _b is None
                  else (min(_a, _b), max(_a, _b)))
        if s1 is None:
            bits = [f"from {s0:.0f} mm to the vessel ends"]
            stations = (("start", s0),)
        else:
            bits = [f"{s1 - s0:.0f} mm of aorta"]
            stations = (("start", s0), ("end", s1))
        if self.base is not None:
            for tag, s in stations:
                nm, loc, _ = self.route.locate(s)
                pos, tan = self.tree.segments[nm].at_arclength(loc)
                sd = section_diameter(self.base, pos, tan, near_point=pos)
                if sd:
                    bits.append(f"{tag} ⌀ {sd['ead_mm']:.1f}")
        self.region_lbl.setText(" · ".join(bits) +
                                ("\n" + "; ".join(self._region_notes)
                                 if self._region_notes else ""))

    # ---------------- trims (cut before a split) ----------------

    def trims_within(self, terminals):
        """The trims whose whole group lies inside `terminals`.

        A trim only means something where every opening it swallows is present. An
        iliac trim offered while you are looking at the visceral segment would
        remove openings that are not on screen.
        """
        labels = {t.label for t in terminals}
        return [tr for tr in self.trims
                if {t.label for t in tr["cand"]["terminals"]} <= labels and
                {t.label for t in tr["cand"]["terminals"]}]

    def _apply_trims(self, terminals, trims):
        """Replace each trim's group of openings with its single interior cut."""
        gone = set()
        for tr in trims:
            gone |= {t.label for t in tr["cand"]["terminals"]}
        return ([t for t in terminals if t.label not in gone]
                + [tr["terminal"] for tr in trims])

    def region_rows(self):
        """The openings inside the active region, with its trims applied."""
        _bounds, inside, _ = self._region_terminals()
        return self._apply_trims(inside, self.trims_within(inside))

    def active_terminals(self):
        """The openings actually being cut: real terminals, minus the ones each trim
        swallows, plus one interior cut per trim.

        With a region active the answer is different in kind -- the region's own two
        bounds plus whatever branches leave from inside it -- because everything
        outside is discarded rather than cut. Trims still apply inside it: a region
        used to drop them entirely, so trimming the iliacs while looking at the
        infrarenal aorta recorded the trim and changed nothing.
        """
        if self.tree is None:
            return []
        if self.region is not None:
            bounds, _inside, _ = self._region_terminals()
            return self.region_rows() + bounds
        return self._apply_trims(self.tree.terminals, self.trims)

    def _row_limit(self, t):
        """How far this cut's slider may travel back from its tip.

        An interior (trim) cut is hemmed in at BOTH ends -- run it too far back and
        the plane arrives at the junction feeding the vessel, where it tears. So
        stop a couple of radii short of that, rather than at the segment's end.
        """
        L = t.segment.length
        if getattr(t, "is_interior", False):
            r = max(float(getattr(t.segment, "_radius_mm", None) or 5.0), 1.0)
            # ...but never past the vessel itself: on a segment shorter than 2 radii
            # the "safe" limit would otherwise sit beyond its far end.
            return float(np.clip(L - 2.0 * r, MIN_OFFSET_MM + 1.0,
                                 MAX_FRAC_OF_SEGMENT * L))
        return MAX_FRAC_OF_SEGMENT * L

    def _trim_label(self, cand):
        names = [t.name for t in cand["terminals"]]
        if len(names) <= 3:
            return " + ".join(names)
        return f"{names[0]} + {len(names) - 1} more"

    def _trim_name(self, cand):
        """What the trimmed opening IS: the feeding vessel, named if we know it."""
        vessel = self.anat_vessels.get(cand["segment"])
        return vessel if vessel else f"Trim before {self._trim_label(cand)}"

    def _refresh_trim_menu(self):
        """Rebuild the 'trim a branch group' menu from what is still trimmable."""
        if not hasattr(self, "trim_btn"):
            return
        menu = self.trim_btn.menu()
        menu.clear()
        cands = (pre_bifurcation_candidates(self.tree, self.base)
                 if self.tree is not None else [])
        # A trim already inside another trim's removed subtree is not offered again.
        taken = set()
        for tr in self.trims:
            taken |= set(tr["cand"]["subtree"]) | {tr["cand"]["segment"]}
        # Inside a region, only offer trims whose whole group is in view. The menu
        # used to list every candidate in the tree, so the visceral segment offered
        # to trim the iliacs, which are not part of what is being cut.
        visible = None
        if self.region is not None:
            _b, inside, _ = self._region_terminals()
            visible = {t.label for t in inside}
        n = 0
        for c in cands:
            if c["segment"] in taken:
                continue
            if visible is not None and not {t.label for t in c["terminals"]} <= visible:
                continue
            seg = self.tree.segments[c["segment"]]
            vessel = self.anat_vessels.get(c["segment"], c["segment"])
            act = menu.addAction(
                f"{vessel}  ,  before {self._trim_label(c)} ({seg.length:.0f} mm)")
            act.setToolTip(f"Cut {c['segment']} before it splits, removing "
                           f"{len(c['terminals'])} openings as one.")
            act.triggered.connect(lambda _=False, cc=c: self._add_trim(cc))
            n += 1
        self.trim_btn.setEnabled(n > 0)
        self.trim_btn.setText(f"Trim a branch group…  ({n})" if n
                              else "Trim a branch group…")

    def _add_trim(self, cand):
        """Replace a group of openings with ONE cut on the vessel feeding them."""
        seg = self.tree.segments[cand["segment"]]
        if getattr(seg, "_radius_mm", None) is None:
            seg._radius_mm = segment_radius(self.base, seg)   # band width needs it
        t = Terminal(seg, cand["at_end"], label=f"{seg.name}.trim", mode="normal",
                     name=self._trim_name(cand))
        t.is_interior = True
        t.subtree_segments = set(cand["subtree"])
        # Open it where it actually cuts. A trim sits on a short interior segment and
        # the usual 6 mm default is often already past the point where its parent
        # takes over -- KY160's brachiocephalic trunk is clean 2 mm before its split
        # and torn by 6 -- so the first thing the user would see is a red status bar
        # on a cut they were told they could make. The search is the same one the
        # pull-back button uses, so this lands on the deepest offset that works.
        start = min(DEFAULT_OFFSET_MM, max(seg.length * 0.5, MIN_OFFSET_MM))
        try:
            found = suggest_pullback(self.base, self.tree, t, start)
        except Exception:                                       # noqa: BLE001
            found = None
        if found is not None:
            self._offset_memory[t.label] = found[0]

        # A nested trim is redundant once its parent is trimmed away.
        self.trims = [tr for tr in self.trims
                      if tr["cand"]["segment"] not in cand["subtree"]]
        self.trims.append({"cand": cand, "terminal": t})
        self._rebuild_after_trim(f"Trimming {self._trim_label(cand)} as one opening "
                                 f"on {seg.name}.")

    def _remove_trim(self, terminal):
        """Undo a trim: its openings come back as individual cuts."""
        self.trims = [tr for tr in self.trims if tr["terminal"] is not terminal]
        self._rebuild_after_trim("Trim removed, those openings are cut individually "
                                 "again.")

    def _rebuild_after_trim(self, message):
        self._build_rows()
        self._on_change()
        if self.problem is None:
            self._say(message)

    def _apply_colour_labels(self):
        """Match each row's name colour to its branch, or clear when colouring is off."""
        on = self.colour_cb.isChecked()
        for r in self.rows:
            r.set_colour(self.opening_colour.get(r.terminal.label) if on else None)

    def _colour_toggled(self, on):
        self.settings.setValue("colour", "true" if on else "false")
        self._apply_colour_labels()
        self.draw()

    def _nudge_branches(self, delta_mm):
        """Shift every branch cut by *delta_mm* (negative = toward the tips)."""
        if self.tree is None:
            return
        branch_rows = [r for r in self.rows if not r.terminal.is_trunk]
        if not branch_rows:
            self._say("No branches to move, this vessel has only trunk ends.", bad=True)
            return
        before = {r.terminal.label: r.value() for r in branch_rows}
        for r in branch_rows:
            r.set_value(r.value() + delta_mm)   # set_value clamps to the row's range
        moved = [r for r in branch_rows
                 if abs(r.value() - before[r.terminal.label]) > 1e-6]
        if not moved:
            self._say("Branches are already at the end of their range, "
                      "nothing to move.", bad=True)
            return
        clamped = len(branch_rows) - len(moved)
        self._on_change()
        if clamped and self.problem is None:    # note the ones that hit their limit
            self._say(f"Moved {len(moved)} branch(es); {clamped} already at their "
                      f"limit.")

    def _send_cut_to_mark(self):
        """Move the cut of the opening the last mark sits on to that marked point."""
        if not self._marks:
            self._say("Place a mark first (tick Measure, then click a spot on a vessel).",
                      bad=True)
            return
        mk = self._marks[-1]
        seg_name = mk["res"]["segment"]
        world = np.asarray(mk["world"], float)
        cands = [r for r in self.rows if r.terminal.segment.name == seg_name]
        if not cands:
            self._say("That mark isn't on a vessel that has its own cut.", bad=True)
            return
        # A segment can carry two openings (an unbranched trunk); send the nearer one.
        row = min(cands, key=lambda r: float(
            np.linalg.norm(np.asarray(r.terminal.tip, float) - world)))
        t = row.terminal
        s = float(mk["res"]["s_mm"])
        requested = (t.segment.length - s) if t.at_end else s
        row.set_value(requested)             # clamps to the row's allowed range
        placed = row.value()
        self._on_change()
        if self.problem is not None:
            return                     # recut already showed why the cut was refused
        if abs(placed - requested) > 0.1:
            self._say(f"{t.name}'s cut placed at {placed:.1f} mm from its tip, "
                      f"the mark is beyond this cut's allowed range.")
        else:
            self._say(f"Moved {t.name}'s cut to the mark "
                      f"({placed:.1f} mm from its tip).")

    def _set_enabled(self, on):
        self.reset_btn.setEnabled(on)
        self.export_btn.setEnabled(on)
        self.colour_cb.setEnabled(on)
        has_branch = on and any(not t.is_trunk for t in self.tree.terminals)
        self.branch_bar.setEnabled(has_branch)

    def _on_select(self, row):
        """Only one advanced panel open at a time, that is what keeps it tidy."""
        for r in self.rows:
            if r is not row and r.expand.isChecked():
                r.expand.setChecked(False)
        self.sel = row
        self.draw()

    # ---------------- cutting ----------------

    def _on_pullback(self, row):
        """Move one cut back to the deepest offset that still cuts cleanly.

        The search is the cutter's own -- it probes candidate offsets with the real
        rule and accepts one only if the rim closes on the plane and the cut takes
        no more than a square millimetre of any vessel it is not cutting. So the
        button cannot offer something the cutter will then refuse, and it never
        moves a cut CLOSER to the junction, which would be advice to repeat the
        mistake.
        """
        if self.tree is None or self.base is None:
            return
        cur = row.value()
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            found = suggest_pullback(self.base, self.tree, row.terminal, cur)
        except Exception as e:                                  # noqa: BLE001
            found = None
            self._say(f"{row.terminal.name}: could not search for a clean cut ({e})",
                      bad=True)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        if found is None:
            self._say(f"{row.terminal.name}: no offset along this vessel has a clean "
                      f"cross-section of its own, so there is nowhere to pull back to.",
                      bad=True)
            return
        off, frac = found
        if off >= cur - 0.05:
            self._say(f"{row.terminal.name} is already at the deepest clean cut "
                      f"({cur:.0f} mm).")
            return
        row.set_value(off)
        self._say(f"{row.terminal.name} pulled back to {off:.0f} mm "
                  f"({frac * 100:.0f}% along the vessel), the deepest cut here that "
                  f"leaves a flat opening.")
        self._on_change()

    def all_rows(self):
        """Every row showing a cut: the branch rows, plus a region's own two bounds.

        Kept separate from self.rows on purpose -- that list is what carries an
        OFFSET, and a region bound's offset comes from the region sliders instead
        (see cut_list). Anything that reports a result belongs here; anything that
        reads or writes a depth belongs on self.rows.
        """
        return list(self.rows) + list(getattr(self, "bound_rows", []))

    def cut_list(self):
        """(terminals to cut, offsets), the single source of truth for a re-cut."""
        act = [r.terminal for r in self.rows]
        offs = {r.terminal.label: r.value() for r in self.rows}
        for t in getattr(self, "_region_bounds", []):
            act.append(t)
            offs[t.label] = float(getattr(t, "_region_offset", DEFAULT_OFFSET_MM))
        return act, offs

    def offsets(self):
        return self.cut_list()[1]

    def _retry_bound_clearance(self):
        """A region boundary undercutting an ostium gets nudged clear of it.

        Same rule as CutSession, here for the same reason the trunk backoff is: a
        fix in one path and not the other is how the window missed the region-reach
        retry while the benchmark called it shipped. See clear_bound_offset.

        A BOUNDARY THE USER DRAGGED IS THEIRS. region_plan already refuses to clear
        ostia once _region_manual is set -- "a dragged boundary is the user's, do not
        nudge" -- and this retry has to obey the same rule or it undoes it: moving a
        hand-placed bound up to 8 mm to make it cut is exactly the slider-snapping
        this repository fixed once already. Reported from the tool as the start
        jumping when it is moved a millimetre.

        So it fires only on a boundary a PRESET placed. If a station someone chose by
        hand will not cut, the honest answer is the refusal and the reason.
        """
        if (self.region is None or not self.problem
                or getattr(self, "_region_manual", False)
                or getattr(self, "_bound_clear_tried", False)):
            return
        m = re.search(r"(\d+) openings, expected (\d+)", str(self.problem))
        if not m or int(m.group(1)) <= int(m.group(2)):
            return                        # too FEW openings is a different fault
        act, offs = self.cut_list()
        bounds = [t for t in act if getattr(t, "_region_offset", None) is not None]
        if not bounds:
            return
        attr = vessel_attribution(self.base, self.tree)
        snap = snap_tolerance(self.base) if self.clean_cb.isChecked() else 0.0
        moved = {}
        for t in bounds:
            was = float(offs[t.label])
            found = clear_bound_offset(self.base, self.tree, t, was,
                                       snap_mm=snap, attr=attr)
            if found is not None:
                moved[t.name] = (was, found)
        if not moved:
            return
        keep = dict(getattr(self, "_bound_nudge", None) or {})
        self._bound_clear_tried = True
        try:
            self._bound_nudge = dict(keep)
            self._bound_nudge.update({n: f for n, (_w, f) in moved.items()})
            self._cut_steps = None        # the planes moved; the cache has not
            self._recut_once()
            if self.problem is None:
                self._bounds_moved = [{"name": n, "from_mm": round(w, 2),
                                       "to_mm": round(f, 2)}
                                      for n, (w, f) in moved.items()]
            else:
                self._bound_nudge = keep  # no better: put them back as planned
                self._cut_steps = None
                self._recut_once()
        finally:
            self._bound_clear_tried = False

    def _retry_trunk_offset(self):
        """A trunk end that will not cut flat gets one retry at an offset that does.

        Same rule as CutSession, and it has to be here too rather than only in the
        API: a fix that lives in one path and not the other is how the desktop window
        quietly missed the region-reach retry while the benchmark reported it as
        shipped. See flat_trunk_offset for why a MODEL END's offset may be moved when
        no other cut's may.

        The row's own number is updated, so the window shows where the cut actually
        went rather than where it was asked to go.
        """
        if getattr(self, "_trunk_backoff_tried", False):
            return
        act, offs = self.cut_list()
        t = next((x for x in act
                  if str(self.problem).startswith(str(x.name))
                  and is_model_end(self.tree, x)), None)
        # Or an opening count that names nobody: FIND the cut at fault rather than
        # guess at it. Same rule as CutSession._double_opening_cut. A plane running
        # through a neighbouring vessel opens that vessel too, and the message then
        # names no terminal -- JK14_2's inlet crosses the celiac trunk and JK30_4's
        # Branch 4 crosses a neighbour, so keying this to model ends alone misses it.
        count = re.search(r"(\d+) openings, expected (\d+)", str(self.problem))
        if t is None and count and int(count.group(1)) > int(count.group(2)):
            attr = vessel_attribution(self.base, self.tree)
            snap = snap_tolerance(self.base) if self.clean_cb.isChecked() else 0.0
            for x in act:
                if getattr(x, "_region_offset", None) is not None:
                    continue        # a bound has its own remedy, along the aorta
                try:
                    out = cut_terminal(self.base, self.tree, x, float(offs[x.label]),
                                       snap_mm=snap, _attr=attr, _suggest=False)
                except Exception:                                   # noqa: BLE001
                    continue
                if len(boundary_loops(out)) > 1:
                    t = x
                    break
        row = next((r for r in self.rows if r.terminal.label == t.label), None) \
            if t is not None else None
        if row is None:
            return
        want = float(offs[t.label])
        found = flat_trunk_offset(
            self.base, self.tree, t, want,
            snap_mm=snap_tolerance(self.base) if self.clean_cb.isChecked() else 0.0)
        if found is None:
            return
        self._trunk_backoff_tried = True
        try:
            row.set_value(found)              # set_value guards its own signals
            self._cut_steps = None            # the plane moved; the cache has not
            self._recut_once()
            self._trunk_moved = ({"name": str(t.name), "from_mm": round(want, 2),
                                  "to_mm": round(found, 2)}
                                 if self.problem is None else None)
            if self.problem is not None:      # no better: put it back as asked
                row.set_value(want)
                self._cut_steps = None
                self._recut_once()
        finally:
            self._trunk_backoff_tried = False

    def recut(self):
        """Re-apply every cut, and give a failing REGION one retry at a wider or
        narrower boundary plane.

        Only on failure and only once. Probing up front costs two extra cuts on every
        rebuild, and a region slider rebuilds on every nudge -- doing it there made
        the whole suite stop finishing. Same shape as the pull-back button, which
        also only looks for a working offset once a cut has refused.
        """
        if self.tree is None:
            return
        self._recut_once()
        if self.problem is not None:
            self._retry_trunk_offset()
        if self.problem is not None:
            self._retry_bound_clearance()
        if (self.problem is None or self.region is None
                or getattr(self, "_reach_tried", False)):
            return
        bounds = list(getattr(self, "_region_bounds", []))
        if any(getattr(t, "band_mm", None) for t in bounds):
            return                              # a hand-set plane size is the answer
        self._reach_tried = True
        try:
            got = fit_region(self.base, self.tree, self.route, *self.region)
            if got is None:
                return
            scale = got["scale"]
            if abs(got["moved_mm"]) > 1e-9:
                self._region_moved_mm = got["moved_mm"]
                self.region = (got["s0"], got["s1"])
                self._sync_region_bar()
            elif abs(scale - 1.0) < 1e-9:
                return
            # Kept as a SCALE on the window, not as a band_mm on these Terminal
            # objects. _region_terminals rebuilds them on every slider nudge and
            # saves their plane attributes into _bound_planes -- so a fitted band
            # written here would come back looking hand-set, and the fitter would
            # then refuse to touch it for the rest of the session.
            self._region_reach = scale
            self._apply_region_reach(bounds)
            self._cut_steps = None          # the planes moved; the cache has not
            self._build_rows()               # the bounds may have moved
            self._recut_once()
            if self.problem is None:
                bits = []
                if abs(scale - 1.0) > 1e-9:
                    bits.append(f"plane reach {scale:.2f}x the automatic size")
                if abs(got["moved_mm"]) > 1e-9:
                    bits.append(f"boundary moved {got['moved_mm']:+.0f} mm along "
                                f"the aorta")
                self._say("Ready — " + ", ".join(bits) + ", which is what cuts here.")
        finally:
            self._reach_tried = False

    def _recut_once(self):
        # The rows ARE the cut list (plus a region's own bounds, which the region
        # sliders drive): a trim replaces the openings it swallows, so cutting
        # tree.terminals here would re-cut branches that no longer exist.
        active, offs = self.cut_list()

        # A region can collapse to nothing once its branch ostia are cleared. With no
        # cuts to make, every check below passes trivially on the untouched watertight
        # mesh -- "Ready, 0 openings, euler 2" -- and export would write the CLOSED
        # original under an _open name. Refuse instead.
        if not active:
            self.cut = self.base.copy()
            self.matched = {}
            self.problem = (
                "Nothing is being cut, this region collapses to nothing once its "
                "branch ostia are cleared. Widen it, or go back to the whole vessel."
                if self.region is not None else
                "Nothing is being cut, there are no openings selected.")
            for r in self.all_rows():
                r.show_result(None)
            self._say(self.problem, bad=True)
            return

        pre = preflight(self.tree, offs, terminals=active)
        if pre is not None:
            self.problem = pre
            if self.cut is None:
                self.cut = self.base.copy()
            # A drag can leave a live section diameter in the readout; the cut it
            # belonged to was just refused, so put each row back to its last
            # committed opening diameter rather than let it contradict the (red)
            # status bar.
            for r in self.all_rows():
                r.show_result(self.matched.get(r.terminal.label))
            self._say(pre, bad=True)
            return

        snap = snap_tolerance(self.base) if self.clean_cb.isChecked() else 0.0

        # RESUME rather than redo. Moving one slider used to re-cut every opening from
        # the base mesh, so a five-opening model paid for five cuts to change one. Each
        # cut's result depends only on the cuts before it, so keeping the mesh after
        # each one lets a re-cut restart at the first opening whose offset actually
        # moved. Dragging the last row then costs one cut instead of five.
        #
        # Keyed on the whole step -- which terminal, at what offset, with what snap --
        # so any change to any of them invalidates that step and everything after it.
        # `_cut_steps` is dropped wholesale whenever the model or the row list changes,
        # since the mesh it starts from would no longer be the same mesh.
        # ...and CUT THE ONE THAT MOVED LAST, so the prefix stays reusable. Openings
        # are cut in whatever order the rows happen to be in, and dragging the first
        # row invalidated everything after it -- a full rebuild every time, while
        # dragging the last row cost one cut. The order is free to change: cutting
        # these in any sequence gives a byte-identical mesh, checked over eight
        # permutations each on a 5-opening and an 11-opening model, because each cut
        # removes a different vessel. So a changed opening is moved to the end, and
        # the next drag of that same slider resumes from everything before it.
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

        # meshes[k] is the mesh after the first k cuts, so meshes[0] is the base.
        self.problem = None
        meshes = (cache["meshes"][:reuse] if reuse else []) + [m]
        for t in active[reuse:]:
            try:
                m = cut_terminal(m, self.tree, t, offs[t.label], snap_mm=snap)
            except Exception as e:
                self.problem = f"{t.name}: {e}"
            meshes.append(m)
        # Only a run of cuts that all succeeded can be resumed from: a step that threw
        # left the mesh untouched, and replaying from there would silently skip it.
        self._cut_steps = ({"base_id": id(self.base), "order": order,
                            "steps": want, "meshes": meshes}
                           if self.problem is None else None)
        self.cut = keep_largest_component(m)
        self.matched = match_openings_to_terminals(self.cut, self.tree, offs,
                                                   terminals=active)
        if self.problem is None:
            # Geometry as well as topology: a cut that escapes into the parent
            # vessel still leaves one loop and the right euler number, so topology
            # alone cannot see it. Only the opening's position can.
            self.problem = validate_cut(self.cut, self.tree, offs, terminals=active)

        for r in self.all_rows():
            r.show_result(self.matched.get(r.terminal.label))
            # Every message that ends in "pull it back" opens with the opening's name.
            r.flag_pullback(bool(self.problem)
                            and self.problem.startswith(r.terminal.name))

        if self.problem is None:
            n = rim_sliver_count(self.cut)
            # A cut the tool MOVED is said out loud. The slider shows the new number
            # too, but the status line is what people read.
            moved = getattr(self, "_trunk_moved", None)
            self._say(f"Ready, {len(boundary_loops(self.cut))} openings, "
                      f"euler {int(self.cut.euler_number)}"
                      f"{f', {n} rim sliver(s)' if n else ', no slivers'}"
                      + (f". {moved['name']} trimmed at {moved['to_mm']:g} mm "
                         f"instead of {moved['from_mm']:g} mm, which is where it "
                         f"comes out flat" if moved else "")
                      + ("". join(
                          f". {b['name']} moved "
                          f"{abs(b['to_mm'] - b['from_mm']):g} mm to clear a branch "
                          f"it was undercutting"
                          for b in (getattr(self, "_bounds_moved", None) or []))))
        else:
            self._say(self.problem, bad=True)

    def _on_change(self):
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            self.recut()
            self.draw()
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

    def _on_preview(self, row=None):
        """Live while a slider is dragged: move the planes and show the diameter at
        the plane, but leave the (expensive) re-cut for release.

        A full re-cut is ~0.5 s, too slow per drag frame; moving the discs and
        sectioning the mesh once (~3 ms) is cheap, so the plane and its diameter
        track the drag and the mesh catches up on release.
        """
        if self.tree is None:
            return
        # One render for the whole frame. Redrawing the discs rendered the window
        # once per disc, which on a five-opening model was most of a 217 ms drag
        # step -- the plane lagged the cursor for frames nobody asked to see.
        p = self.plotter
        p.suppress_rendering = True
        try:
            self._draw_planes()
            if row is not None and self.base is not None:
                t = row.terminal
                point, outward = t.cut_plane(row.value())
                row.show_live_diameter(
                    section_diameter(self.base, point, outward, near_point=point))
        finally:
            p.suppress_rendering = False
        p.render()

    def _show_quality(self, probe):
        """Say what the supplied centerline actually is, in numbers."""
        turn = probe.worst_turning
        if probe.jagged:
            self.quality_lbl.setText(
                f"⚠ Voxel staircase: {turn:.0f}°/mm of tangent turning. "
                f"Cut planes would swing wildly. Smoothing suggested.")
            self.quality_lbl.setStyleSheet("color:#a3231a; font-size:11px;")
        else:
            self.quality_lbl.setText(
                f"✓ Smooth curve: {turn:.1f}°/mm of tangent turning. "
                f"No smoothing needed, leave it off.")
            self.quality_lbl.setStyleSheet("color:#1c6b34; font-size:11px;")

    def _smooth_toggled(self, on):
        self.settings.setValue("smooth_on", "true" if on else "false")
        if self.csv_path and self.base is not None:
            self.load_model(keep_smooth_choice=True)

    def _smooth_amount_changed(self, v):
        self.settings.setValue("smooth_mm", v)
        if self.smooth_cb.isChecked() and self.csv_path and self.base is not None:
            self.load_model(keep_smooth_choice=True)

    def _reset_family(self):
        """Back to Auto-detect, a declaration is about one case, not the next."""
        if self.family.currentIndex() != 0:
            self.family.blockSignals(True)   # no reload: the file has not loaded yet
            self.family.setCurrentIndex(0)
            self.family.blockSignals(False)

    def _family_changed(self, _):
        """Re-read the centerline under the declared anatomy.

        The whole classification hangs off this -- which openings are trunk ends,
        what they are named, which regions are offered -- so the tree has to be
        rebuilt rather than relabelled. keep_smooth_choice keeps the user's
        smoothing decision and their export name across the reload.
        """
        if self.csv_path and self.base is not None:
            self.load_model(keep_smooth_choice=True)

    def _snap_toggled(self, on):
        self.settings.setValue("snap", "true" if on else "false")
        if self.base is not None and self.tree is not None:
            self._on_change()

    def _opacity_changed(self, v):
        self.opacity_lbl.setText(f"{v}%")
        self.settings.setValue("opacity", v)
        if self.mesh_actor is not None:
            self.mesh_actor.GetProperty().SetOpacity(v / 100.0)
            self.plotter.render()

    def _say(self, text, bad=False):
        self.status.setText(text)
        self.status.setStyleSheet(
            "background:#fdecea; color:#a3231a; border:1px solid #f5c2bd;"
            "border-radius:4px; padding:6px; font-weight:600;" if bad else
            "background:#eaf7ee; color:#1c6b34; border:1px solid #bfe4cb;"
            "border-radius:4px; padding:6px; font-weight:600;")

    # ---------------- drawing ----------------

    def _opening_face_rgb(self, shown):
        """Per-face RGB tint by which opening's vessel owns each face; None if off.

        Reuses vessel_coords (the same face-to-segment ownership the cutter uses),
        so a branch's surface, its plane disc, and its row label all share a colour.
        """
        if (self.tree is None or not self.colour_cb.isChecked()
                or not self.opening_colour):
            return None
        # The whole per-face tint depends only on the cut mesh (the palette is fixed
        # per model), so cache it: vessel_coords is a ~50 ms KD query and the
        # per-branch connectivity walk adds more, but a row select / opacity change
        # redraws the same cut mesh and can reuse the result.
        key = (id(shown), len(shown.faces))
        if getattr(self, "_rgb_cache_key", None) == key:
            return self._rgb_cache

        # vessel_attribution, not vessel_coords: same answer, but it carries the
        # model's attribution forward instead of re-deriving it for the cut mesh, so
        # a redraw costs a lookup rather than 130 ms of re-deciding ownership.
        owner, _, radial, _r = vessel_attribution(shown, self.tree)
        centers = np.asarray(shown.triangles_center)
        by_seg = {}
        for r in self.rows:
            if r.terminal.label in self.opening_colour:
                by_seg.setdefault(r.terminal.segment.name, []).append(r.terminal)

        out = np.tile(np.array(_hex_to_rgb(MESH_COLOR), np.uint8), (len(owner), 1))
        for seg, terms in by_seg.items():
            if len(terms) != 1:
                continue                      # shared segment: ambiguous -> leave red
            t = terms[0]
            R = max(float(getattr(t.segment, "_radius_mm", None) or 5.0), 1.0)
            # Two gates together. RADIAL: nearest-owner AND within a couple of the
            # branch's own radii of its centerline, which rejects the fat-aorta-wall
            # patches the KD rule mis-captures near an ostium. CONTINUITY: keep only
            # the ONE surface-connected piece grown from the branch's own tip, so a
            # patch that survives the radial gate but is detached from the branch
            # (a scrap sitting over on another branch) is not painted this colour.
            cand = (owner == seg) & (radial <= RADIAL_TINT_FACTOR * R)
            if not cand.any():
                continue
            ci = np.flatnonzero(cand)
            tip_face = int(ci[np.argmin(
                np.linalg.norm(centers[ci] - np.asarray(t.tip, float), axis=1))])
            keep = _component(shown, cand, [tip_face])
            out[keep] = _hex_to_rgb(self.opening_colour[t.label])

        self._rgb_cache_key, self._rgb_cache = key, out
        return out

    def _draw_planes(self):
        """Just the cut-plane discs. Named actors, so re-adding replaces in place."""
        active, offs = self.cut_list()
        colour_on = self.colour_cb.isChecked()
        # Only rebuild the discs that moved. Dragging one slider re-made all of them,
        # every frame, and a disc is a 64-segment polygon plus an actor -- the four
        # that did not move cost as much as the one that did.
        keys = getattr(self, "_plane_keys", None)
        if keys is None:
            keys = self._plane_keys = {}
        live = self.plotter.renderer.actors
        for t in active:
            r = next((x for x in self.all_rows() if x.terminal is t), None)
            point, outward = t.cut_anchor(offs[t.label])   # disc shows the band
            radius = band_width(t)                 # the disc IS the cut's reach
            is_sel = (self.sel is r)
            if colour_on and t.label in self.opening_colour:
                colour = self.opening_colour[t.label]       # a branch: its own hue
            else:
                colour = PLANE_SEL if is_sel else PLANE_COLOR   # trunk / colour off
            name = f"plane_{t.label}"
            key = (tuple(np.round(point, 4)), tuple(np.round(outward, 4)),
                   round(float(radius), 3), str(colour), bool(is_sel))
            if keys.get(name) == key and name in live:
                continue                 # unchanged, and still on screen
            keys[name] = key
            self.plotter.add_mesh(
                pv.Disc(center=point, inner=0, outer=radius, normal=outward, c_res=64),
                color=colour,
                opacity=0.60 if is_sel else 0.30,   # the selected disc reads darker
                name=f"plane_{t.label}",
                pickable=False,          # a translucent disc must never catch the
            )                            # hover picker instead of the vessel wall

    def draw(self, reset_camera=False):
        # Build the scene with the render window switched off, and show it once at the
        # end. Every add_mesh renders by default, so drawing eighteen actors rendered
        # the window eighteen times -- 19 renders for one re-cut, and the frames
        # nobody saw were most of what made a slider release feel slow.
        p = self.plotter
        p.suppress_rendering = True
        try:
            self._draw_scene(p, reset_camera)
        finally:
            p.suppress_rendering = False
        p.render()

    def _draw_scene(self, p, reset_camera=False):
        cam = None if reset_camera else p.camera_position
        p.clear()                        # also removes hover markers and pinned marks
        self._plane_keys = {}            # the discs went with it, so none are current
        self.mesh_actor = None
        self._hover_shown = False
        self._hover_last = None
        if reset_camera:                 # a fresh model, old marks no longer apply
            self._marks = []
            self._update_marks_label()

        shown = self.cut if self.cut is not None else self.base
        if shown is None:
            return

        ext = np.asarray(shown.bounds[1]) - np.asarray(shown.bounds[0])
        self._hover_marker_r = max(0.006 * float(np.linalg.norm(ext)), 1.2)

        face_rgb = self._opening_face_rgb(shown)   # None unless colouring is on
        if face_rgb is not None:
            self.mesh_actor = p.add_mesh(
                to_pv(shown), scalars=face_rgb, rgb=True, preference="cell",
                smooth_shading=True, opacity=self.opacity.value() / 100.0,
                specular=0.2, name="mesh")
        else:
            self.mesh_actor = p.add_mesh(
                to_pv(shown), color=MESH_COLOR, smooth_shading=True,
                opacity=self.opacity.value() / 100.0, specular=0.2, name="mesh")

        if self.tree is not None:
            for seg in self.tree.segments.values():
                p.add_mesh(pv.Spline(seg.pts, 300), color=CL_COLOR, line_width=4,
                           name=f"cl_{seg.name}")
            self._draw_planes()

            offs = self.offsets()
            active = [r.terminal for r in self.rows] or self.tree.terminals
            for i, loop in enumerate(boundary_loops(shown)):
                pts = shown.vertices[loop]
                c = pts.mean(axis=0)
                near = min(active,
                           key=lambda t: np.linalg.norm(
                               t.cut_plane(offs.get(t.label, DEFAULT_OFFSET_MM))[0] - c))
                p.add_points(pts,
                             color=TRUNK_LOOP if near.mode == "z" else BRANCH_LOOP,
                             point_size=10, render_points_as_spheres=True,
                             name=f"loop_{i}")

        # Pinned marks are data, not tied to the mesh, re-render them after the
        # clear() above so they survive a re-cut.
        self._draw_marks()

        # So is the measure plane. It marks a station on the centerline, which a cut
        # does not move, and having it vanish on the next slider release would make
        # it useless for the one thing it is for: looking at a cut and its
        # measurement together.
        if (self.mplane_box.isVisible() and self._m_pose is not None
                and self.mplane_show.isChecked()):
            pos, tan, dia = self._m_pose
            r = 0.8 * dia if dia > 0.1 else 8.0
            p.add_mesh(pv.Disc(center=pos, inner=0, outer=max(r, 2.0),
                               normal=tan, c_res=64),
                       color="#8e24aa", opacity=0.55, name="measure_plane",
                       pickable=False)

        # Only the mesh may be picked, so the hover cell-picker can't latch onto a
        # cut-plane disc, a centerline, or a mark in front of the surface.
        for name, act in p.renderer.actors.items():
            try:
                act.SetPickable(name == "mesh")
            except Exception:
                pass

        if reset_camera:
            p.reset_camera()
        elif cam is not None:
            p.camera_position = cam
        p.render()

    # ---------------- actions ----------------

    def reset(self):
        # Reset is where a stuck user goes, so it puts EVERYTHING back: the region,
        # the trims, the hand-declared arteries and every plane. Trims used to
        # survive, on the reasoning that a trimmed group is a structural choice
        # rather than a viewing mode. That reasoning is wrong from the chair: someone
        # who presses Reset wants the model they loaded, and leaving three branches
        # collapsed into one makes the button look broken in precisely the situation
        # it exists for.
        if self.region is not None:
            self._set_region(None)
        if self.trims:
            self.trims = []
        if self.tree is not None and getattr(self.tree, "manual_names", None):
            self.tree.manual_names = {}
            try:
                self.anat_names, self.anat_vessels, self.anat_notes = name_anatomy(
                    self.tree, self.base)
                for t in self.tree.terminals:
                    if t.label in self.anat_names:
                        t.name = self.anat_names[t.label]
            except Exception:                                   # noqa: BLE001
                pass
        self._offset_memory = {}
        self._bound_planes = {}
        self._region_reach = None
        self._build_rows()          # the trims are gone: the rows are a different set
        # Silence the rows while they are put back. Each one emits `changed` as it is
        # set, and every emission was re-cutting the whole model, validating it and
        # redrawing the scene -- six full passes for one button press, 3.9 seconds of
        # them, five of which nobody ever saw because the next one overwrote it.
        for r in self.rows:
            r.blockSignals(True)
            try:
                r.terminal.reset_plane()
                r._reset_plane()
                r.set_value(DEFAULT_OFFSET_MM)
            finally:
                r.blockSignals(False)
        self._on_change()

    def export(self):
        if self.tree is None:
            self._say("Nothing to export, load a centerline to cut the mesh.", bad=True)
            return
        if self.problem is not None:
            QtWidgets.QMessageBox.warning(
                self, "Cut is not valid",
                f"{self.problem}\n\nAdjust the offending opening until the status "
                f"bar turns green, then export.")
            return

        try:
            self.out_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            self._say(f"Cannot write to {self.out_dir}: {e}", bad=True)
            return

        stem = self.export_stem()
        ext = self.fmt.currentText()
        out = self.out_dir / f"{stem}.{ext}"
        side = self.out_dir / f"{stem}_cut_planes.json"

        # Re-exporting the same case is the normal way to work here -- adjust a cut,
        # export, look at it, adjust again -- so overwriting is not an accident to
        # guard against with a modal. The export is derived data, rebuilt in seconds
        # from the mesh, the centerline and the JSON written beside it, and the default
        # file name follows the mesh, so two cases do not collide. It is reported in
        # the status bar instead, which says what happened without interrupting.
        replaced = out.exists()

        self.cut.export(out)

        offs = self.offsets()
        n_loops = len(boundary_loops(self.cut))
        area_mm2 = float(self.cut.area)
        lo, hi = CAMP_REF_AREA_MM2
        meta = {
            # File NAMES, not paths. The name identifies which case this came from,
            # which is all a downstream reader needs; the full path additionally
            # published wherever the file happened to live, a lab share layout on
            # a desktop, and a server's storage layout once this is hosted.
            "mesh": Path(self.mat_path).name if self.mat_path else None,
            "centerline": Path(self.csv_path).name if self.csv_path else None,
            "euler_number": int(self.cut.euler_number),
            "n_openings": n_loops,
            "n_faces": len(self.cut.faces),
            # Total surface area of the open cut. This is the x-axis a curvature
            # comparison is plotted against, so it belongs next to the topology
            # rather than being recomputed by every reader from the STL.
            "area_mm2": round(area_mm2, 2),
            "snap_mm": round(snap_tolerance(self.base), 4)
                       if self.clean_cb.isChecked() else 0.0,
            "rim_slivers": rim_sliver_count(self.cut),
            "min_triangle_quality": round(min_triangle_quality(self.cut), 5),
            # Anatomical names are inferred from the geometry under this frame. It is
            # recorded so a downstream reader can check the left/right convention
            # rather than trusting it, a mirrored volume would swap the sides.
            "anatomy": {
                "frame": ANATOMY_FRAME,
                "left_right_verified": False,
                "named_openings": len(self.anat_names),
                # Which names were DECLARED rather than inferred. A reader checking
                # the left/right convention needs to know that these ones did not
                # come from the frame above and are not subject to it -- and a
                # measurement placed against a hand-named renal should say so.
                "named_by_hand": dict(getattr(self.tree, "manual_names", None) or {}),
                "names": dict(self.anat_names),
                "notes": list(self.anat_notes),
            },
            # For the curvature pipeline downstream. It remeshes before measuring,
            # and a remesh is allowed to move vertices but not to open or close a
            # hole, so it records what the topology was BEFORE remeshing and
            # asserts against it afterwards. This tool is the only stage that knows
            # the intended opening count, which makes the check cheap here and
            # guesswork anywhere else.
            "camp": {
                "target_faces": CAMP_TARGET_FACES,
                "pre_remesh_euler": int(self.cut.euler_number),
                "pre_remesh_loops": n_loops,
                "reference_area_mm2": list(CAMP_REF_AREA_MM2),
                "area_in_reference_range": bool(lo <= area_mm2 <= hi),
            },
            "openings": [],
        }
        if self.region is not None and self.route is not None:
            # What was isolated, and where along the vessel, so the measurement can
            # be reproduced or compared across cases without re-deriving it.
            _a, _b = self.region
            s0, s1 = ((float(_a), None) if _b is None
                      else (min(_a, _b), max(_a, _b)))
            ends = {}
            # An open-ended region has no distal station to measure at.
            stations = (("start", s0),) if s1 is None else (("start", s0), ("end", s1))
            for tag, s in stations:
                nm, loc, _ = self.route.locate(s)
                pos, tan = self.tree.segments[nm].at_arclength(loc)
                sd = section_diameter(self.base, pos, tan, near_point=pos)
                ends[tag] = round(sd["ead_mm"], 2) if sd else None
            prof = []
            for s in np.arange(s0, s1, 2.0):
                nm, loc, _ = self.route.locate(s)
                pos, tan = self.tree.segments[nm].at_arclength(loc)
                sd = section_diameter(self.base, pos, tan, near_point=pos)
                if sd:
                    prof.append((round(float(s - s0), 2), round(sd["ead_mm"], 2)))
            meta["region"] = {
                "start_mm_along_aorta": round(s0, 2),
                "end_mm_along_aorta": round(s1, 2),
                "length_mm": round(s1 - s0, 2),
                "aorta_route_length_mm": round(self.route.length, 2),
                "start_diameter_mm": ends["start"],
                "end_diameter_mm": ends["end"],
                "max_diameter_mm": max((d for _, d in prof), default=None),
                "branches_inside": [r.terminal.name for r in self.rows],
                "boundary_notes": list(self._region_notes),
                "diameter_profile_mm": prof,
            }
        for t in [r.terminal for r in self.rows]:
            s = self.matched.get(t.label)
            point, outward = t.cut_plane(offs[t.label])
            meta["openings"].append({
                "name": t.name,
                "spline": t.label,
                # What the tool inferred, kept separate from `name` so a hand-typed
                # rename is distinguishable from the automatic label.
                "anatomical_name": self.anat_names.get(t.label)
                                   or self.anat_vessels.get(t.segment.name),
                "cut_mode": t.mode,
                # A trim is one opening standing in for the group it removed.
                "trimmed_group": sorted(t.subtree_segments)
                                 if getattr(t, "is_interior", False) else None,
                "offset_mm": round(offs[t.label], 3),
                "tilt_a_deg": round(t.tilt_a_deg, 2),
                "tilt_b_deg": round(t.tilt_b_deg, 2),
                "shift_a_mm": round(t.shift_a_mm, 2),
                "shift_b_mm": round(t.shift_b_mm, 2),
                "plane_size_mm": round(band_width(t), 2),
                "plane_origin": [round(float(v), 4) for v in point],
                "plane_normal": [round(float(v), 6) for v in outward],
                "diameter_mm": round(s["diameter"], 2) if s else None,
                "centroid": [round(float(c), 4) for c in s["centroid"]] if s else None,
            })
        side.write_text(json.dumps(meta, indent=2))
        self._say(f"{'Replaced' if replaced else 'Exported'} {out.name} + "
                  f"{side.name} → {self.out_dir}")


def _opt(args, flag):
    """Pull "--flag value" out of *args*, returning the value or None."""
    if flag not in args:
        return None
    i = args.index(flag)
    val = args[i + 1] if i + 1 < len(args) else None
    del args[i:i + 2]
    return val


def main():
    args = [a for a in sys.argv[1:] if a != "--"]
    # Optional and entirely additive: --nubs opens the window on a case with the
    # nub-in batch cut already made, so a pass over the cohort starts from what the
    # preset produced instead of a default offset on every opening. Nothing about
    # the window changes when they are absent.
    nubs = _opt(args, "--nubs")
    out_dir = _opt(args, "--out")
    name = _opt(args, "--name")

    csv = mat = None
    if len(args) == 1:
        mat = args[0]                      # mesh only
    elif len(args) >= 2:
        csv, mat = args[0], args[1]

    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")
    w = CutWindow(csv, mat)
    w.show()
    if nubs or out_dir or name:
        # Inside the try: nub_session is part of the cohort tooling and is NOT in the
        # collaborator distribution, so on that copy this import is the failure, not
        # anything it would have done. A missing extra must not stop the window
        # opening on a case somebody wants to cut by hand.
        try:
            import nub_session
            if nubs:
                for note in nub_session.apply_from_file(w, nubs, out_dir, name):
                    print(f"  {note}")
            else:
                nub_session.apply(w, {"cuts": []}, out_dir, name)
        except Exception as exc:                                   # noqa: BLE001
            # A seed that fails must not cost the window. The case is still loaded
            # and still cuttable by hand, which is the whole point of opening it.
            print(f"  could not restore the nub cut: {exc}")
    sys.exit(app.exec_() if hasattr(app, "exec_") else app.exec())


if __name__ == "__main__":
    main()
