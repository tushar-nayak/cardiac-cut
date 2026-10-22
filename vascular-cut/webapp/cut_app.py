# The cut tool in a browser, built to match the desktop window.
#
# Same VTK pipeline as branch_gui.py: a PyVista plotter, rendered on the server and
# streamed to the browser by trame. That is what lets the cut planes stay live and
# draggable instead of degrading to a static preview.
#
#   python webapp/cut_app.py            -> http://localhost:8080
#
# The geometry, the ranges and the defaults are the desktop's exactly, so a cut made
# here is the cut made there; anything that CUTS differently is a bug. The page itself
# is organised for someone who has never seen the desktop tool: panels in the order
# the work goes, plain labels, and the desktop's own name for a control in its
# tooltip (Axial is its Z-plane, Keep more its "toward tips", Remove rim slivers its
# "Sliver-free rim").

import asyncio
import copy
import io
import json
import os
import re
import secrets
import sys
import time
import zipfile
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "webapp") not in sys.path:
    sys.path.insert(0, str(ROOT / "webapp"))

pv.OFF_SCREEN = True

from trame.app import get_server, asynchronous             # noqa: E402
from trame.ui.vuetify3 import SinglePageWithDrawerLayout   # noqa: E402
from trame.widgets import vuetify3 as v, vtk as vtkw, html, client  # noqa: E402

from vascular_api import (                                 # noqa: E402
    CutSession, CutError, CUT_MODES, DEFAULT_OFFSET_MM, MIN_OFFSET_MM,
)
from branch_cut import SIDES_ASSUMED_NOTE                  # noqa: E402
from scene import Scene, planes_for                        # noqa: E402
from shape import (                                        # noqa: E402
    ShapeJob, ShapeError, FACE_TARGETS, TRIM_OPTIONS, SEED_OPTIONS,
)
import shape_metrics as SM                                 # noqa: E402
# The Cardiac section's geometry: a closed ventricle split into inner and outer walls
# and cut at the base. The library layer, like branch_cut for vessels.
from heart_split import HeartSession, HeartError           # noqa: E402
# The zone labelling, from the LIBRARY layer. It used to live in
# tools/zone_stub_bench.py, which this app cannot import -- nothing that ships
# may depend on a bench directory, and tools/ is never on this app's sys.path.
from region_station import (                               # noqa: E402
    final_labels, zone_face_colours, ZONE_COLOURS, ZONE_NAMES,
)
# At module top, deliberately. Importing the geometry stack LAZILY, inside a handler,
# after trame and vtk were already loaded, killed the server outright on Windows --
# a stack overflow in pandas' import, exit 139, no traceback, and a page whose socket
# simply closed. Nothing here may be imported from inside a trigger.
import wslink.chunking as _wsc                              # noqa: E402
from aiohttp import web                                     # noqa: E402
from submissions import (                                   # noqa: E402
    SubmissionStore, UploadRejected, StorageError, KINDS,
    NOTICE_HEADING, NOTICE_PARAGRAPHS, notice_placeholders,
)
from visit_gate import VisitGate                            # noqa: E402

SAMPLE_CSV = ROOT / "KY160_1_CentrelinePoints.csv"
SAMPLE_MAT = ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"

server = get_server(client_type="vue3")
state, ctrl = server.state, server.controller

# The stylesheet lives in its own file so it can be edited without touching Python,
# and is injected with client.Style -- NOT html.Style, which renders nothing: the
# class names reached the browser and the rules did not, so the panel kept Vuetify's
# stock look while appearing to be styled.
# Serve the stylesheet as a static file and let trame link it. Every in-band route
# failed: html.Style renders a <style> Vue never puts in the DOM, client.Style binds
# to a state key it never fills, and a mounted script cannot reach `document` --
# trame evaluates those expressions in a sandbox without it, and the exception took
# the whole toolbar down with it.
STATIC = Path(__file__).resolve().parent / "static"
# Vuetify's own module FIRST, so its stylesheet is linked before app.css and the rules
# here win on order rather than on !important. The layout would enable it anyway, but
# later -- after app.css -- which is why the old stylesheet needed !important everywhere.
from trame_vuetify.module import vue3 as _vuetify_module    # noqa: E402
server.enable_module(_vuetify_module)
server.enable_module({
    "styles": ["vc_static/app.css"],
    # Uploads and downloads over HTTP beside the websocket -- see vc_upload.js.
    # vc_boot.js: the copy-ID helper, and the drawer starting closed on a phone.
    "scripts": ["vc_static/vc_upload.js", "vc_static/vc_boot.js"],
    "serve": {"vc_static": str(STATIC)},
})

# Shown in the footer. If this string is not on screen, the browser is running
# cached code and nothing else in this file is being seen either -- which is worth
# being able to tell at a glance rather than inferring.
BUILD = "ui-24"

# The page's colours, in one place. Vuetify turns each into a theme colour
# (--v-theme-<name>, and .bg-<name> / .text-<name>), and the refused page below is
# formatted from the same dict; app.css carries the full scale of each (iris-25 to
# iris-800, mist-25 to ink). One accent, iris, marks what you can act on or have
# chosen, and is never drawn as a small dot or square, so it cannot be read as a branch
# swatch. The app bar is night. `sel` is the selected cut disc's orange (scene.py
# plane_active) and appears only where it means that disc. Red, amber and green are
# status only. Light only for now.
UI_COLOURS = {
    "chrome": "#15133A", "on-chrome": "#F4F3FF",
    "background": "#F2F3FA", "on-background": "#161A36",
    "surface": "#FFFFFF", "on-surface": "#161A36", "surface-2": "#EBEDF6",
    "ink2": "#454D6D", "ink3": "#5F6788", "line": "#D9DCE9",
    "primary": "#5A48E6", "on-primary": "#FFFFFF", "accent-soft": "#EFEEFF",
    "sel": "#D55E00", "sel-soft": "#FBEDE3",     # the selected disc only (key, HUD)
    "error": "#B42318", "danger": "#D92D20",
    "warning": "#B54708", "warning-soft": "#FFFAEB",
    "success": "#067647", "info": "#5A48E6",
}
VUETIFY_CONFIG = {
    "theme": {"defaultTheme": "light",
              "themes": {"light": {"dark": False, "colors": UI_COLOURS}}},
    "defaults": {"VTooltip": {"location": "top", "maxWidth": 320, "openDelay": 400}},
    # Below 960 px the drawer becomes an overlay over the view rather than beside it.
    "display": {"mobileBreakpoint": "md"},
}

# classes=, never class_=.
#
# trame maps ("classes" -> "class") and knows nothing about class_, so class_ is
# accepted, ignored, and reaches no element. Every rule aimed at a Vuetify widget was
# therefore dead: .note, .statusbar, .op-name, .mode-tog, .adv-reset and every mr-2
# and mt-2 in the file. Only html.Div carried its classes, which is why .op-card
# looked styled and nothing around it did -- and why this went unnoticed.

# ── who is allowed to reach this, and what it may touch ──────────────────────
#
# The Mesh and Centerline boxes take a PATH and hand it to trimesh.load(). On
# loopback that is not a vulnerability, it is the feature: it is your machine, your
# files, and confining it to the project directory would break the ordinary case of
# pointing at a geometry on a lab share. The moment the listener is not loopback the
# same box is a remote file-read primitive against the server, and the error text
# reports whether a path exists.
#
# So the boundary enforced here is the BINDING, not a path filter. A non-loopback
# host requires VC_PUBLIC=1, which is the configuration that has no path input at
# all. Getting this wrong should be impossible by accident rather than merely
# documented, because the failure is silent and the README is not in the room.
LOOPBACK = {"127.0.0.1", "::1", "localhost"}
PUBLIC = os.environ.get("VC_PUBLIC", "").strip().lower() in {"1", "true", "yes", "on"}

# ── uploads: where they go, and how big ──────────────────────────────────────
# Everything uploaded is KEPT (tools/HOSTING.md job 7), under a generated id, with the
# original filename recorded as metadata and never used as a path. Public mode refuses
# to start without VC_STORE_DIR, because a container's own disk is thrown away with
# the container and "we kept your files" would then be untrue.
STORE = SubmissionStore.from_env(ROOT, public=PUBLIC)
print(f"[uploads] {STORE.root}", flush=True)
MAX_UPLOAD_BYTES = max(STORE.limits.values())

# THE UPLOAD CAP HAS TO BE ENFORCED BEFORE THE APP SEES ANYTHING. After the handshake
# wslink raises every client's reassembly limit to 4 GB, and its un-chunker allocates
# the full claimed size from the FIRST chunk's header, before any content arrives --
# so one 13-byte message could make this server reserve 4 GB (its own source says
# "This un-chunker is vulnerable to DOS"). A size check in an upload handler runs far
# too late: the file is already in memory. So the reassembly limit is clamped to the
# largest upload this instance accepts, plus room for the message around it. An
# honest oversize file never gets this far -- the page refuses it before sending --
# so what this stops is a client that ignored the page.
_wsc_set_max = _wsc.UnChunker.set_max_message_size


# Uploads no longer travel the websocket at all (they go over HTTP, below), so
# nothing legitimate a page sends over it is more than a few kilobytes. The clamp can
# therefore be tight rather than sized to the largest upload.
WS_MESSAGE_CAP = 8 * 2 ** 20


def _clamped_max(self, size):
    _wsc_set_max(self, min(int(size), WS_MESSAGE_CAP))


_wsc.UnChunker.set_max_message_size = _clamped_max


def safe_stem(name, fallback="cut"):
    """A file NAME, never a path.

    `export_name` is joined onto the output directory, so "../../x" walks out of it.
    That is a bug even on loopback -- it writes somewhere the person did not choose
    and the status bar still says it exported -- so this applies in every mode
    rather than only in the hosted one.
    """
    stem = Path(str(name or "")).name                 # drops any directory part
    stem = "".join(c for c in stem if c.isalnum() or c in "-_. ").strip(" .")
    return stem[:120] or fallback

SCENE = Scene(off_screen=True)
# Where the camera starts. Kept because an empty scene's reset_camera() does NOT go back
# here -- it keeps the last person's view direction and zoom -- and the wipe between
# visitors restores it rather than hand the next one a view framed on someone's aorta.
INITIAL_CAMERA = SCENE.plotter.camera_position
SESSION = {"s": None, "sub": None}   # not in state: a CutSession is not JSON;
                                     # sub is the stored submission it came from

# The desktop window's empty states, verbatim. Each says what to do next, which is
# the whole reason there are four of them rather than one blank list.
EMPTY_PRELOAD = "Load a mesh to view it.\nAdd a centerline to cut it."
EMPTY_NO_CL = ("No centerline loaded, viewing the mesh only.\n"
               "Add its centerline under Input, then load again, to cut it.")
EMPTY_REGION_ENDS = ("No branches leave the vessel inside this region. Only its "
                     "two ends are cut. Widen the region, or go back to the whole "
                     "vessel.")
EMPTY_REGION_COLLAPSED = ("This region collapses to nothing once its branch ostia "
                          "are cleared, so nothing is being cut. Widen it, or go "
                          "back to the whole vessel.")

# Named, not inline: this dict and _SHAPE_STATE below are the only complete definition
# of a blank page, and the wipe between visitors replays them (VISIT_DEFAULTS).
_CUT_STATE = {
    "trame__title": "VascularCut",
    # Never in public mode. These went to every visitor in the first state message
    # and sat in the (disabled) boxes on screen, and they are this machine's absolute
    # paths -- the account name included.
    "csv_path": "" if PUBLIC else (str(SAMPLE_CSV) if SAMPLE_CSV.exists() else ""),
    "mat_path": "" if PUBLIC else (str(SAMPLE_MAT) if SAMPLE_MAT.exists() else ""),
    "up_mesh": None,
    "up_cl": None,
    "up_mesh_pct": -1,
    "up_cl_pct": -1,
    "vc_token": "",
    "up_mesh_note": "",
    "up_cl_note": "",
    "up_mesh_bad": False,
    "up_cl_bad": False,
    "submission_id": "",
    "max_mesh_bytes": STORE.limits["surface"],
    "max_cl_bytes": STORE.limits["centerline"],
    "loaded": False,
    "family": "auto",
    "smooth": True,
    "clean_rims": False,
    "show_centerline": True,
    "show_planes": True,
    "show_loops": True,
    "colour_branches": True,
    "colour_zones": False,
    "zones_busy": False,
    "zones_note": "",
    "zone_legend": [],
    "marks": [],
    "marks_note": "",
    "opacity": 55,
    "openings": [],
    "selected": None,
    "openings_title": "Openings",
    "empty_text": EMPTY_PRELOAD,
    "status_text": "No case loaded",
    "status_bad": False,
    "region_options": ["Whole vessel"],
    "region_pick": "Whole vessel",
    "region_s0": 0.0,
    "region_s1": 0.0,
    # True while the region has NO distal end -- one plane, everything below it kept.
    # The End controls are disabled and their number is not the user's, so nothing may
    # read region_s1 without checking this first. See _sync_region_boxes.
    "region_open": False,
    "route_length": 0.0,
    "in_region": False,
    "region_name": "",
    "anat_text": "",
    "quality_text": "",
    "busy": False,
    "area_mm2": 0.0,
    "n_openings": 0,
    "euler": 0,
    "n_faces": 0,
    "export_name": "cut",
    "export_msg": "",
    "export_bad": False,
    "angle": None,
    "nudge_mm": 2.0,
    "trim_options": [],
    "trim_pick": None,
    "can_trim": False,
    "artery_options": [],
    "can_nub": False,
    "can_nudge": False,
    "measure": None,
    "measuring": False,
    # The measure plane: where it sits, what the vessel is there, and the range it
    # may be moved over. Nothing to place until something has been measured.
    "m_station": 0.0,
    "m_lo": 0.0,
    "m_hi": 0.0,
    "m_note": "",
    "m_show": True,
    "cut_modes": [{"title": t, "value": k} for k, t in CUT_MODES.items()],
    # ── what the page shows about the case at a glance ───────────────────────
    # The longer story behind the one-line status (e.g. each nub that had to move).
    "status_notes": [],
    # The anatomy call, split out so the page can badge it and gate on it: nub-in
    # and the Ishimaru zones exist for thoracic aortas only.
    "is_thoracic": False,
    "anat_kind": "",
    "anat_how": "",
    "anat_level": "info",       # info | warn | bad -- how much to trust the call
    "case_name": "",
    "cut_ok": False,
    # What the case was LOADED with, so the page can tell when the Anatomy or Smooth
    # boxes no longer match it and offer to apply them.
    "loaded_family": "auto",
    "loaded_smooth": True,
    "load_seq": 0,              # bumped on every load; the page collapses Input on it
    "region_items": [{"title": "Whole vessel", "value": "Whole vessel",
                      "subtitle": ""}],
    "region_why": "",
}
state.update(_CUT_STATE)

# ── shape metrics ────────────────────────────────────────────────────────────
#
# The second half of the site, and deliberately separable from the first: delta K
# takes any open triangular surface, so someone who already has geometries never has
# to touch the cut tool. What links them is a handoff, not a dependency.
JOB = ShapeJob()

# Which surface the camera has already framed. Framing on every redraw would yank the
# view back to a default every time the colouring changed; never framing leaves a
# newly loaded surface off screen. Frame each source once.
_FRAMED = {"src": None}

def _repeat_text(band, at):
    """How far two runs of the same surface can disagree at a trim. Neutral wording:
    it is a property of the setting, not a warning about this surface."""
    med, worst = band
    return (f"Repeatability at trim {at:g}: two runs of the same surface differ by "
            f"{med:.1f}× typically, {worst:.1f}× at worst.")


# Numbers as strings, because a VBtnToggle compares its buttons' value ATTRIBUTES to
# the model. An int in state never equals the string the attribute carries, so the
# reference setting would come up with nothing selected.
_SHAPE_STATE = {
    "mode": "cut",
    "sm_path": "",
    "sm_faces": str(SM.REFERENCE["target_faces"]),
    "sm_trim": f"{SM.REFERENCE['trim_frac']:g}",
    # One partition is one draw. ONE is the default for EVERY session of this app,
    # local runs included -- there is no hosted-only path, and an earlier comment
    # here wrongly implied there was. The reason is hosting (a small CPU-only
    # container, tools/HOSTING.md, where four partitions is four times the runtime of
    # the expensive half of the site), but the consequence is not hosting-specific.
    # The error that buys is not hidden: +/-16% is printed beside every appearance of
    # the value. Raise it to 4 or 8 whenever you are measuring rather than looking.
    "sm_seeds": "1",
    "sm_outlines": False,
    "sm_dim_trimmed": True,
    "sm_range": "",
    "sm_scale": "robust",     # robust | full | manual
    "sm_clip": 5.0,           # per cent clipped off each end in robust mode
    "sm_cmin": 0.0,
    "sm_cmax": 1.0,
    "sm_full": "",            # the data range, so a manual bound has something to aim at
    "sm_info": None,
    "sm_result": None,
    "sm_notes": [],
    "sm_headline": "",
    "sm_busy": False,
    "sm_phase": "",
    "sm_elapsed": 0.0,
    "sm_msg": "",
    "sm_bad": False,
    "sm_export_name": "delta_k",
    "sm_export_msg": "",
    "sm_export_bad": False,
    "sm_camp": SM.campsrc_origin(),
    "sm_ref_partition_cv": f"{100 * SM.PARTITION_CV:.0f}",
    "sm_ref_faces": f"{SM.REFERENCE['target_faces']:,}",
    "sm_ref_trim": f"{SM.REFERENCE['trim_frac']:g}",
    "sm_ref_n": SM.REFERENCE["n_cases"],
    # How much two runs of the same surface can disagree at the CURRENTLY SELECTED
    # trim. Kept beside the value rather than in a footnote: without it there is
    # nothing to stop someone reading a 40% difference as a finding. Seeded here
    # rather than left to the change handler, which does not fire until the trim moves.
    "sm_repeat": _repeat_text(*SM.repeatability_band()),
    # The reference face count as the toggle's own string, so the page can say whether
    # the chosen settings match the cohort's without comparing a number to a string.
    "sm_ref_faces_n": str(SM.REFERENCE["target_faces"]),
}
state.update(_SHAPE_STATE)

# ── cardiac ──────────────────────────────────────────────────────────────────
#
# The third section: a closed left-ventricle myocardium, split into its inner and
# outer walls and cut at the base by a plane (heart_split.py). Like Cut, it hands its
# surfaces to Shape rather than measuring them itself.
HEART = {"s": None, "seq": 0}
# The cardiac meshes: the project's dataset/, beside this folder.
CARDIAC_DATASET = ROOT.parent / "dataset"
# Where the section exports walls: the project's results/, one folder per ventricle.
CARDIAC_WALLS = ROOT.parent / "results" / "cardiac_walls"
HEART_DEPTH_MM = 3.0          # how far below the fitted rim the plane starts


def _hs_case_items():
    if PUBLIC or not CARDIAC_DATASET.is_dir():
        return []
    return [{"title": p.stem, "value": str(p)}
            for p in sorted(CARDIAC_DATASET.glob("*.stl")) if not p.name.startswith("._")]


_HEART_STATE = {
    "hs_cases": _hs_case_items(),
    "hs_path": (_hs_case_items() or [{"value": ""}])[0]["value"],
    "hs_loaded": False,
    "hs_busy": False,
    "hs_name": "",
    "hs_depth": HEART_DEPTH_MM,
    "hs_tilt_a": 0.0,
    "hs_tilt_b": 0.0,
    # Holes are left open by default: CAMP does not need a closed surface, and a
    # patch adds curvature the wall does not have (CAMP_INTEGRATION_GUIDE.md, Q1).
    "hs_fill": False,
    "hs_fill_mm": 60.0,
    "hs_show_inner": True,
    "hs_show_outer": True,
    "hs_show_removed": True,
    "hs_show_plane": True,
    "hs_ok": False,
    "hs_problem": "",
    "hs_info": None,
    "hs_notes": [],
    "hs_msg": "",
    "hs_bad": False,
    "hs_export_name": "heart",
    "hs_export_msg": "",
    "hs_export_bad": False,
}
state.update(_HEART_STATE)

# A blank page: what a browser gets from a fresh process. Every key the app ever sets
# is in one of the dicts above, so replaying this resets all of it. Deep-copied at
# each use: the lists in it are mutable, and state holds on to what it is given.
VISIT_DEFAULTS = copy.deepcopy({**_CUT_STATE, **_SHAPE_STATE, **_HEART_STATE})


def _say(text, bad=False, notes=None):
    """One line in the status bar; anything longer goes in `notes`, which the page
    lists under Details. Every call clears the notes of the call before."""
    state.status_text = text
    state.status_bad = bool(bad)
    state.status_notes = list(notes or [])


# ── state <-> session ────────────────────────────────────────────────────────
def _sync_openings():
    s = SESSION["s"]
    if s is None:
        state.openings = []
        state.openings_title = "Openings"
        state.is_thoracic = False
        state.cut_ok = False
        return
    rows = []
    declared = s.declared_arteries()
    # Offered per row rather than globally: the list depends on the anatomy, and an
    # abdominal aorta must not be offered a subclavian.
    state.artery_options = [{"title": "Auto (from the geometry)", "value": ""}] + [
        {"title": n, "value": n} for n in s.anatomy_choices()]
    # The colour each opening's branch is drawn in, so its card can carry the same
    # swatch. The aorta itself is drawn in the mesh colour.
    try:
        colours = s.opening_colours() or {}
    except Exception:                                           # noqa: BLE001
        colours = {}
    for o in s.openings():
        rows.append({
            "label": o["label"],
            "name": o["name"],
            # Rounded for the box: float arithmetic otherwise shows 179.20000000000002.
            "offset": round(float(o["offset_mm"]), 1),
            "colour": colours.get(o["label"]) or ("#c9424a" if o.get("is_trunk")
                                                  else ""),
            "limit": float(o["limit_mm"]),
            "diameter": o["diameter_mm"],
            "mode": o["cut_mode"] or "z",
            "bound": bool(o["is_region_bound"]),
            "tilt_a": o["tilt_a_deg"],
            "tilt_b": o["tilt_b_deg"],
            "shift_a": o["shift_a_mm"],
            "shift_b": o["shift_b_mm"],
            "band": o["band_mm"],
            "band_auto": o["band_auto"],
            # Ranges are per-opening and move with the plane size, so they travel
            # with the row rather than being global.
            "lim_tilt": o["limit_tilt"],
            "lim_shift": o["limit_shift"],
            "lim_band_min": o["limit_band_min"],
            "lim_band_max": o["limit_band_max"],
            "is_trim": bool(o.get("is_trim")),
            "trim_segment": o.get("trim_segment"),
            "is_trunk": bool(o.get("is_trunk")),
            # The artery this opening has been DECLARED to be, if any. Distinct from
            # the name: a name is a label, a declaration drives the region presets.
            "artery": declared.get(o["label"], ""),
            # Set below, once the status is known: every message that ends in "pull
            # it back" opens with the opening's name, so the row it belongs to can
            # show its own way out instead of leaving the reader to find it.
            "bad": False,
        })
    prob = None if s.ok else (s.problem or "")
    if prob:
        for r in rows:
            r["bad"] = prob.startswith(r["name"])
    state.openings = rows

    # What can still be trimmed. Rebuilt with the rows, because trimming one group
    # can remove another candidate: a trim inside a trim's removed subtree is not
    # offered again.
    cands = s.trim_candidates()
    state.trim_options = [
        {"title": f"{c['vessel']}  ·  before {c['group']}  ({c['length_mm']:.0f} mm)",
         "value": c["segment"]}
        for c in cands]
    state.can_trim = bool(cands)
    # Branches, not openings: a trunk-only model has nothing to nudge.
    state.can_nudge = any(not r["is_trunk"] and not r["bound"] for r in rows)
    # The nub preset squares off the aorta's own ends, which a region has already
    # replaced with its bounds -- so it is a whole-model action or nothing. And it is
    # for thoracic aortas only, where it was measured; the API refuses the rest.
    state.is_thoracic = bool(s.is_thoracic)
    state.can_nub = s.region is None and s.is_thoracic

    # The group title says what is being cut, so a region is never mistaken for the
    # whole vessel just because the list got shorter.
    if s.region is not None:
        inside = sum(1 for r in rows if not r["bound"])
        state.openings_title = (f"Openings, {len(rows)}: {inside} "
                                f"branch{'' if inside == 1 else 'es'} inside the region")
        state.empty_text = (EMPTY_REGION_ENDS if any(r["bound"] for r in rows)
                            else EMPTY_REGION_COLLAPSED)
    else:
        state.openings_title = f"Openings, {len(rows)}"
        state.empty_text = EMPTY_NO_CL


def _anat_caption(s):
    """What the panel says about the anatomy, and whether to believe it.

    A centerline running outside the mesh makes every radius meaningless, and the
    family call, the vessel names and the region presets are all measured off those
    radii. The cut refuses on the same evidence, but only once it is attempted -- by
    which time the wrong anatomy has been on screen the whole time. So a mismatch
    REPLACES the caption rather than qualifying it.
    """
    a = s.anatomy_summary()
    if not a.get("mesh_matches", True):
        why = (f"only {100 * a['coverage']:.0f}% of the centerline lies within it"
               if a["coverage"] < 0.98 else
               f"{len(a['outside_mesh'])} of {a['n_terminals']} centerlines run "
               f"outside it")
        return (f"This mesh and this centerline do not look like the same export: "
                f"{why}. Nothing here has been read from the anatomy. The "
                f"usual cause is a trunk-only mesh paired with a branched centerline.")
    hand = len(s.declared_arteries())
    # Appended to whichever caption is returned below, never returned on its own:
    # the close-call branch returns early, so writing this only into the ordinary
    # branch lost it on exactly the cases carrying BOTH problems -- a marginal
    # family call AND unwitnessed sides. Those are the least trustworthy cases in
    # the cohort, and they were the ones told least.
    sides = ""
    if any(SIDES_ASSUMED_NOTE in str(nt) for nt in a.get("notes", [])):
        sides = (" LEFT AND RIGHT ARE NOT CONFIRMED here: nothing in the anatomy "
                 "witnessed which side is which, so the sides follow only from the "
                 "LPS frame. A mirrored export would swap every one of them. Check "
                 "against the L/P/S marker in the corner of the view before relying "
                 "on a side.")
    # A CLOSE family call leads, because it changes every name below it. On the
    # centerlines this tool actually cuts the two families overlap, so a marginal
    # reading is not a curiosity -- it is the case most likely to be wrong.
    if a.get("close_call"):
        other = "abdominal" if a["kind"] == "thoracic" else "thoracic"
        return (f"Read as {a['kind']}, but it was a CLOSE call "
                f"({a['outflow_split']:.2f} against {a['family_cutoff']:.2f}): "
                f"{other} is arguable, and every name below follows from it. Set the "
                f"anatomy by hand if this is wrong. {a['named_openings']} opening(s) "
                f"named{f', {hand} by hand' if hand else ''}.{sides}")
    # No opening could be named. The desktop says what to do about it and the web
    # said nothing, which is the JK25_1 situation exactly: no renal named, so the
    # infrarenal presets are never offered and the region list collapses to one
    # entry, with nothing on screen explaining why or how to get out of it.
    if not a.get("named_openings"):
        return (f"Treated as {a['kind']} ({a['how']}). No opening could be named "
                f"unambiguously: name them yourself (open an opening’s details with "
                f"its arrow and set Artery), and the region presets follow. A named "
                f"renal is what "
                f"brings back the aortic neck and the infrarenal regions.{sides}")
    # anatomy_summary() has carried "notes" all along and this caption dropped them.
    # The one that matters says the tool found no witness in the anatomy for which
    # side is which -- no brachiocephalic on an arch, no splenic/hepatic pair on an
    # abdomen -- so every Left/Right below rests on the frame alone. Left unsaid, the
    # caption reads identically whether handedness was READ or GUESSED.
    if sides:
        return (f"Treated as {a['kind']} ({a['how']}). {a['named_openings']} "
                f"opening(s) named{f', {hand} by hand' if hand else ''}.{sides}")
    return (f"Treated as {a['kind']} ({a['how']}). "
            f"{a['named_openings']} opening(s) named anatomically"
            f"{f', {hand} by hand' if hand else ''}, assuming an LPS frame. "
            f"Confirm left/right on a new case.")


def _anat_facts(s):
    """(kind, how, level) for the page's anatomy badge.

    level is how much to trust the call, and so how loudly the caption speaks: "bad"
    when the mesh and centerline do not match, "warn" for a close family call, no
    opening named, or sides not confirmed, and "info" otherwise -- the routine caption
    used to be a warning on every case, so it stopped meaning anything.
    """
    a = s.anatomy_summary()
    how = "set by hand" if a.get("declared") else "detected"
    if not a.get("mesh_matches", True):
        level = "bad"
    elif (a.get("close_call") or not a.get("named_openings")
          or any(SIDES_ASSUMED_NOTE in str(nt) for nt in a.get("notes", []))):
        level = "warn"
    else:
        level = "info"
    return a["kind"], how, level


def _set_anatomy_state(s):
    """Everything on the page that follows from the names and the family call."""
    presets = s.presets(planar=True)
    state.region_options = ["Whole vessel"] + [p["name"] for p in presets]
    state.region_items = ([{"title": "Whole vessel", "value": "Whole vessel",
                            "subtitle": "The whole model, every opening cut"}]
                          + [{"title": p["name"], "value": p["name"],
                              "subtitle": p.get("why", "") or ""} for p in presets])
    state.anat_text = _anat_caption(s)
    state.anat_kind, state.anat_how, state.anat_level = _anat_facts(s)


def _sync_anatomy():
    """The region list and the anatomy caption, after the names have changed.

    Both are derived from the names rather than from the geometry, so declaring an
    artery has to rebuild them -- that is the entire point of declaring one. The
    current region choice is preserved where it survives; a preset that no longer
    exists falls back to the whole vessel rather than pointing at nothing.
    """
    s = SESSION["s"]
    if s is None:
        return
    _set_anatomy_state(s)
    if state.region_pick not in state.region_options:
        state.region_pick = "Whole vessel"


def _refresh_status():
    s = SESSION["s"]
    if s is None:
        state.cut_ok = False
        return
    st = s.status()
    _say(st["message"], bad=not st["ok"])
    state.cut_ok = bool(st["ok"])
    state.n_openings = st.get("n_openings", 0)
    if st["ok"]:
        state.area_mm2 = st["area_mm2"]
        state.euler = st["euler_number"]
        state.n_faces = st["n_faces"]


def _redraw(full=True):
    s = SESSION["s"]
    if s is None:
        SCENE.clear()
    else:
        SCENE.style["mesh_opacity"] = max(0.05, float(state.opacity) / 100.0)
        if full:
            rgb = None
            if (state.colour_zones and ZONES["rgb"] is not None
                    and ZONES["for_cut"] == _zone_cut_key(s)):
                rgb = ZONES["rgb"]
            elif state.colour_branches:
                try:
                    rgb = s.face_colours()
                except Exception:                               # noqa: BLE001
                    rgb = None      # a tint is never worth losing the surface over
            SCENE.set_mesh(s.mesh, valid=s.ok, face_rgb=rgb)
            SCENE.set_centerline(s.tree, show=state.show_centerline)
            SCENE.set_loops(s.mesh if state.show_loops else None)
        SCENE.set_planes(planes_for(s, state.selected) if state.show_planes else [])
    if ctrl.view_update:
        ctrl.view_update()


# ── Ishimaru zones ───────────────────────────────────────────────────────────
# Computed on request, never on redraw: 0.7 s on a 4k-vertex mesh and 4.2 s on a 28k
# one, all of it GIL-held, so doing it automatically would stall the whole server on
# every slider release. The result is cached against the cut it came from and dropped
# the moment that cut changes -- a zone map of the PREVIOUS cut drawn over the current
# surface is worse than none, because it still looks authoritative.
ZONES = {"rgb": None, "legend": [], "for_cut": None}


def _zone_cut_key(s):
    """What the labelling depends on. Any change here invalidates the cache."""
    if s is None:
        return None
    return (len(s.mesh.vertices), len(s.mesh.faces), bool(s.ok),
            getattr(s, "region_name", None))


def _zones_clear(note=""):
    ZONES.update(rgb=None, legend=[], for_cut=None)
    state.zone_legend = []
    state.zones_note = note
    # With no map there is nothing to colour by; leaving the switch on would colour
    # nothing and still read as "zones shown".
    state.colour_zones = False


def _zone_colours_for(s):
    """(face rgb, legend, zones present). Pure -- it runs off the event loop."""
    mesh, route, zones, lab, isb, adj, flat, reg, oz = final_labels(s)
    present = sorted({int(z) for z in lab if z >= 0})
    if not present:
        raise RuntimeError("no zone could be placed on this anatomy")
    # zones are the presets the anatomy offers, keyed by NAME. Parse the number out
    # rather than trusting position: on a bovine arch the list is 0, 2, 3, 4 and the
    # third entry is zone 3, not zone 2.
    spans = {}
    for z in zones or []:
        m = re.search(r"zone\s+(\d+)", str(z.get("name", "")), re.I)
        if m:
            spans[int(m.group(1))] = (z.get("s0"), z.get("s1"))
    legend = []
    for z in present:
        s0, s1 = spans.get(z, (None, None))
        span = ("" if s0 is None
                else f"{s0:.0f}-{'end' if s1 is None else format(s1, '.0f')} mm")
        legend.append({"zone": z, "colour": ZONE_COLOURS[z % len(ZONE_COLOURS)],
                       "name": ZONE_NAMES.get(z, ""), "span": span})
    return zone_face_colours(mesh, lab), legend, present


def _zones_outdated(gen, s):
    """True when a zone result must not be installed. Handles the second case itself."""
    if gen != _gen():
        return True     # the visitor left: the next person's page must not get this
    if SESSION["s"] is not s:
        with state:
            state.zones_busy = False     # a different case is on screen now
        return True
    return False


async def _zones_job(gen, s):
    # gen and s are taken when the job is SCHEDULED, not when it first runs: the
    # visitor can leave, and the page be wiped, in between, and a job that read them
    # then would work on the blank page and write its failure into it.
    if gen != _gen() or s is None:
        return
    key = _zone_cut_key(s)
    try:
        loop = asyncio.get_event_loop()
        rgb, legend, present = await loop.run_in_executor(None, _zone_colours_for, s)
    except Exception as e:                                       # noqa: BLE001
        if _zones_outdated(gen, s):
            return
        # `with state:` -- this runs as a task after an await, outside any trame
        # callback, and state changed there is not sent to the page until something
        # else flushes it. The delta K watcher does the same, for the same reason.
        with state:
            _zones_clear(f"Zones unavailable: {e}")
            state.colour_zones = False
            state.zones_busy = False
            _redraw(full=True)
        return
    if _zones_outdated(gen, s):
        return
    with state:
        _zones_finish(rgb, legend, present, key)


def _zones_finish(rgb, legend, present, key):
    ZONES.update(rgb=rgb, legend=legend, for_cut=key)
    state.zone_legend = legend
    # SAY WHY A ZONE IS MISSING, AND SAY THE TRUE REASON. There are two, and they read
    # very differently to someone checking the anatomy. Zones 1 and 2 both end at the
    # left common carotid, so if no vessel was named that, neither has a landmark and
    # neither is drawn. If it WAS named but comes off the brachiocephalic -- a bovine
    # arch -- zone 1 has zero length. This used to claim a shared origin for every
    # missing zone, and on the repo's KY160 sample, which has no left carotid named at
    # all, that was simply false.
    missing = [z for z in range(5) if z not in present and z < max(present)]
    s = SESSION["s"]
    names = {v for v in ((s.anat_names or {}).values() if s else ()) if v}
    if not missing:
        why = ""
    elif "Left common carotid" not in names and set(missing) & {1, 2}:
        why = (": no left common carotid was found on this model, and it is the "
               "landmark those zones end at.")
    else:
        why = ": a shared origin, so it has no length."
    state.zones_note = (
        "Zones " + ", ".join(str(z) for z in present) + "."
        + (" No zone " + ", ".join(str(z) for z in missing) + why if missing else ""))
    state.zones_busy = False
    _redraw(full=True)


@ctrl.set("run_zones")
def run_zones():
    s = SESSION["s"]
    if s is None:
        _zones_clear("Load a model first.")
        return
    # Ishimaru zones are an arch classification. The page does not offer them on an
    # abdominal case; this holds when a client fires the trigger anyway.
    if not s.is_thoracic:
        _zones_clear("Ishimaru zones are defined on the thoracic aorta only, and this "
                     "case reads as abdominal.")
        return
    if not s.ok:
        _zones_clear("Cut the vessel first -- zones are read off the cut surface.")
        return
    if state.zones_busy:
        return
    state.zones_busy = True
    state.colour_zones = True
    state.zones_note = "Finding the ostia and placing the boundaries..."
    asynchronous.create_task(_zones_job(_gen(), s))


def _recut_and_draw():
    s = SESSION["s"]
    if s is None:
        return
    state.busy = True
    try:
        s.recut()
        # The zone map belonged to the cut that has just been replaced.
        if ZONES["rgb"] is not None:
            _zones_clear("The cut changed. Find the zones again.")
        _sync_openings()
        _refresh_status()
        _redraw(full=True)
        # A measurement of the previous cut is worse than none, so drop it. It is
        # recomputed on request rather than automatically: sectioning the mesh at
        # every station takes long enough to be felt on each slider release.
        state.measure = None
    finally:
        state.busy = False


@server.trigger("run_measure")
def run_measure():
    """Equivalent-area diameter along the vessel, or along the isolated region."""
    s = SESSION["s"]
    if s is None:
        return
    state.measuring = True
    try:
        m = s.measure(step_mm=2.0)
        if m is None:
            state.measure = None
            _say("Nothing to measure. Load a centerline first.", bad=True)
            return
        # The profile is hundreds of points and nothing on screen reads it; sending
        # it would just make every state sync bigger.
        state.measure = {k: v for k, v in m.items() if k != "profile"}
        # ...and the angle, which is the other half of what the desktop's button
        # says it does. The fit error travels with it: the arc a best-fit circle
        # subtends is only worth quoting when the centerline is actually an arc, and
        # a whole aorta is not one.
        try:
            a = s.angulation()
            # "Poor fit" decided here, not in the template: a bare > in a v-if closes
            # the tag and the element silently vanishes.
            state.angle = (dict(a, poor_fit=bool(not a.get("straight")
                                                 and float(a.get("rms_mm") or 0) > 1.5))
                           if a else None)
        except Exception:                                       # noqa: BLE001
            state.angle = None
        # Open the plane at the widest station, which is what the headline quotes,
        # so the first thing shown is where that number came from.
        lo, hi = s.measured_range()
        state.m_lo, state.m_hi = round(lo, 1), round(hi, 1)
        _place_measure_plane(float(m["max_at_mm"]))
    finally:
        state.measuring = False


def _measure_clear():
    """No measurement, no measure plane and no marks: all of them belong to the case
    they were taken on."""
    state.measure = None
    state.angle = None
    state.m_note = ""
    state.m_station = state.m_lo = state.m_hi = 0.0
    SCENE.set_measure_plane(None)
    state.marks = []
    _sync_marks([])
    SCENE.set_marks([])


def _place_measure_plane(s_mm):
    """Put the measure plane at a station and say how wide the vessel is there."""
    s = SESSION["s"]
    if s is None:
        return
    st = s.measure_station(s_mm)
    if st is None:
        state.m_note = ""
        SCENE.set_measure_plane(None)
    else:
        state.m_station = st["s_mm"]
        if st["ead_mm"]:
            note = f"⌀ {st['ead_mm']:.1f} mm here"
            lo, hi = st["min_mm"], st["max_mm"]
            if lo and hi and hi - lo > 0.05 * max(hi, 1.0):
                note += f"  ·  {lo:.1f}-{hi:.1f} across"
            state.m_note = note
        else:
            state.m_note = "no cross-section here"
        SCENE.set_measure_plane(st if state.m_show else None)
    if ctrl.view_update:
        ctrl.view_update()


@server.trigger("measure_station")
def measure_station(value=None):
    if value is not None:
        _place_measure_plane(float(value))


@server.trigger("measure_show")
def measure_show(on=True):
    """Hide the disc without losing the measurement, as on the desktop."""
    state.m_show = bool(on)
    if not state.m_show:
        SCENE.set_measure_plane(None)
        if ctrl.view_update:
            ctrl.view_update()
    elif state.m_station is not None:
        _place_measure_plane(float(state.m_station))


def _idx(index):
    """Row index as an int, or None.

    It does not always arrive as one. Iterating the DeepReactive proxy yields the
    v-for index as a STRING, so every handler was raising
    "'<' not supported between instances of 'str' and 'int'" on each event -- the
    thumb still moved, because that is pure client-side, while nothing reached the
    scene. A silent server-side exception is invisible from the browser, so coerce
    here rather than trusting the wire.
    """
    try:
        return int(index)
    except (TypeError, ValueError):
        return None


def _rows():
    """The openings as a list, whatever shape they came back in.

    DeepReactive syncs its proxy back to the server, and the list returns as an
    OBJECT keyed by stringified indices -- so `rows[0]` raises KeyError and the
    handler dies silently mid-drag. Normalise on read; the alternative is every
    caller guessing.
    """
    rows = state.openings
    if isinstance(rows, dict):
        return [rows[k] for k in sorted(rows, key=lambda k: int(k))]
    return list(rows or [])


def _row(index):
    rows = _rows()
    i = _idx(index)
    if i is None or i < 0 or i >= len(rows):
        return None
    return rows[i]


def _put(index, **fields):
    """Write fields back into a row and mark the list dirty.

    Every handler takes its value from the event rather than reading state.openings,
    because a v-model bound to a field INSIDE a v-for does not reliably propagate
    back: trame dirty-checks top-level keys by reference, and mutating an object
    inside the list leaves the reference unchanged. The server then reads the value
    the row had when it was built -- which is why dragging appeared to do nothing.
    Rebuilding the list is what makes the change visible in both directions.
    """
    rows = _rows()
    i = _idx(index)
    if i is None or i < 0 or i >= len(rows):
        return None
    rows[i] = dict(rows[i], **fields)
    state.openings = rows
    state.dirty("openings")
    return rows[i]


# ── loading ──────────────────────────────────────────────────────────────────
@state.change("csv_path", "mat_path")
def on_source_changed(**_):
    """A declared anatomy belongs to the old case, not the next one.

    branch_gui calls _reset_family from every path picker, with the same one-line
    reason: "a declared anatomy belongs to the old case". The web kept it, so
    someone who declared a thoracic aorta by hand and then typed a new path got the
    next case forced to thoracic silently -- and the family call decides which
    openings are trunk ends, what they are named and which regions are offered, so
    it is the single most consequential piece of state to carry over by accident.

    Only the family. Smoothing and rim cleaning are preferences about how the user
    likes to work; the family is a claim about one particular scan.
    """
    if state.family and state.family != "auto":
        state.family = "auto"


def _scrub(text):
    """Error text with the store's location taken out.

    A loader failing on an uploaded file names the file it was reading, and that name
    is a path inside the store -- this server's disk layout, shown to a stranger.
    """
    msg = str(text)
    root = Path(STORE.root)
    forms = {str(root), str(root.resolve())}
    forms |= {f.replace("\\", "/") for f in forms}
    for form in sorted(forms, key=len, reverse=True):
        msg = msg.replace(form, "<upload>").replace(form.lower(), "<upload>")
    return msg


@ctrl.set("load_model")
def load_model():
    if PUBLIC:
        # A hosted instance has no path input, so refuse here as well as not
        # rendering the boxes. The check that matters is the server-side one: any
        # client can write state.csv_path and fire this trigger from the page's own
        # JavaScript, whatever the widgets look like.
        _say("This instance takes uploaded files, not paths. Choose a mesh and a "
             "centerline above; the case loads once both have arrived.", bad=True)
        return
    csv, mat = state.csv_path.strip(), state.mat_path.strip()
    if not mat:
        _say("Pick a mesh first.", bad=True)
        return
    if not csv:
        # The old placeholder said "or leave blank to just view". It never worked:
        # CutSession("") reads Path("."), and the error was a permission-denied on
        # the current folder.
        _say("Pick a centerline too. The tool needs both: it cannot derive a "
             "centerline from the surface.", bad=True)
        return
    # A typed path is not a submission, so nothing it produces is retained.
    SESSION["sub"], state.submission_id = None, ""
    _open_session(csv, mat)


def _session_args():
    """What building a session needs from the page, read ON the event loop."""
    return dict(family=state.family, smooth_mm=None if state.smooth else 0.0,
                clean_rims=state.clean_rims)


def _used(args):
    """What a session was built with, in the page's own terms."""
    return {"family": args["family"], "smooth": args["smooth_mm"] is None}


# The files a case on THIS machine was opened from, so "Apply to this case" can read
# them again with a changed Anatomy or Smooth. Not in SESSION: the wipe between
# visitors checks SESSION's shape exactly, and this is local-mode only anyway.
LOADED_FROM = {"csv": None, "mat": None}


def _open_session(csv, mat):
    """Build the session and draw it, in one go. The typed path, on this machine."""
    state.busy = True
    try:
        args = _session_args()
        _install_session(CutSession(csv, mat, **args), _used(args),
                         name=Path(str(mat)).stem)
        LOADED_FROM.update(csv=csv, mat=mat)
    except CutError as e:
        SESSION["s"], state.loaded, state.cut_ok = None, False, False
        _say(_scrub(e), bad=True)
    except Exception as e:
        SESSION["s"], state.loaded, state.cut_ok = None, False, False
        _say(_scrub(f"{type(e).__name__}: {e}"), bad=True)
    finally:
        state.busy = False


def _install_session(s, used=None, name=None):
    """Put a built session on screen. Touches the scene, so the event loop only.

    `used` is what it was built with ({family, smooth}), so the page can tell when
    the boxes no longer match; `name` is the case's name as its owner gave it.
    """
    SESSION["s"] = s
    state.loaded = True
    state.selected = None
    state.route_length = round(float(s.route.length), 1) if s.route else 0.0
    state.region_pick = "Whole vessel"
    state.in_region = False
    _set_anatomy_state(s)
    if used is not None:
        state.loaded_family = used["family"]
        state.loaded_smooth = bool(used["smooth"])
    if name is not None:
        state.case_name = str(name)
    # A zone map and its legend belong to the previous case; they used to stay on
    # screen over the next one. So did a measurement: its readout, its plane and the
    # pinned marks, all drawn on the new case with the old case's numbers.
    _zones_clear()
    _measure_clear()
    state.load_seq = int(state.load_seq or 0) + 1
    state.quality_text = (
        f"Centerline read as jagged; smoothed {s.smooth_mm:.1f} mm."
        if s.jagged else "Centerline used as supplied.")
    _sync_openings()
    _refresh_status()
    _redraw(full=True)
    SCENE.reset_camera(s.mesh)
    if ctrl.view_reset_camera:
        ctrl.view_reset_camera()


# ── upload ───────────────────────────────────────────────────────────────────
# The two files arrive separately, each from its own control -- .vtp and .vtk are valid
# for both roles, so one multi-picker could not tell which file is which. Each is
# checked the moment it arrives and held in memory; nothing is written until both are
# here, so a half-finished upload leaves nothing on disk. Either can then be replaced
# on its own, and each replacement is a new submission.
PENDING = {"surface": None, "centerline": None}
_UP = {"surface": ("up_mesh_note", "up_mesh_bad"),
       "centerline": ("up_cl_note", "up_cl_bad")}

# The page's side of the HTTP routes: a random token the page sends with each upload.
# It is not authentication -- this app has none -- and it stops another site's page
# from posting files into this one. Every connected browser receives it in state (an
# earlier comment here said only one page could know it; trame gives every page the
# whole state), so in public mode it is minted per VISIT, on admit, and withdrawn when
# the visitor leaves. A tab left over from the last person, or anyone who copied the
# token, then gets a 403 rather than posting into the next person's session. On this
# machine there is no gate, and one token serves the process as before.
HTTP_TOKEN = None if PUBLIC else secrets.token_urlsafe(24)
state.vc_token = HTTP_TOKEN or ""


def _set_token(tok):
    global HTTP_TOKEN
    HTTP_TOKEN = tok
    with state:
        state.vc_token = tok or ""


def _token_ok(sent):
    """Constant-time, and False whenever no token is current -- a bare `!=` against
    None would let a request with no header through."""
    tok = HTTP_TOKEN
    return bool(tok) and secrets.compare_digest(str(sent or "").encode("utf-8"),
                                                tok.encode("utf-8"))
# Downloads waiting to be fetched: one-time link -> (filename, bytes, mime). Only the
# newest few are kept, since each can be tens of MB and an unclicked one is abandoned.
HTTP_DL = {}
HTTP_DL_KEEP = 3
# The byte cap binds before the face limit for ASCII STL, which is about five times the
# size of the same surface in binary: 50 MB of ASCII is roughly 190k faces, 50 MB of
# binary about a million. So an oversize mesh is told the cheap way round it.
_STL_HINT = (" If it is an ASCII STL, saving it as binary STL makes it about five "
             "times smaller.")


def _check_upload(kind, name, data):
    """Check one arriving file. Returns (checked upload, None) or (None, refusal).

    Pure, so it can run off the event loop: face counting a large ASCII STL takes a
    second or two, and that must not freeze the page either. It does NOT hold the
    file. PENDING used to be written here, from the worker thread, where nothing can
    tell whether the person who sent the file is still the one on the page; it is now
    written on the loop, after that has been checked.
    """
    try:
        return STORE.check_upload(kind, name, data), None
    except UploadRejected as e:
        return None, str(e)


def _upload_note(kind):
    c = PENDING[kind]
    size = f"{len(c.data) / 2 ** 20:.1f} MB"
    if c.faces:
        size += f", {c.faces:,} faces"
    return f"{c.name}  ({size})"


async def _http_upload(request):
    """POST /vc/upload/{kind}: one file, streamed, capped before it is read."""
    if not _token_ok(request.headers.get("X-VC-Token")):
        return _out_of_date(403)
    kind = request.match_info.get("kind")
    if kind not in KINDS:
        return web.json_response({"error": "Unknown upload."}, status=404)
    note, bad = _UP[kind]
    label = KINDS[kind][1]
    cap = STORE.limits[kind]
    name = unquote(request.headers.get("X-VC-Name", ""))
    # The visit this upload belongs to. If it ends while the bytes are still arriving,
    # or while they are being checked, what they produce belongs to nobody on the page
    # any more, so it is dropped -- never held for, or shown to, whoever comes next.
    gen = _gen()

    def refuse(msg, status):
        if gen != _gen():
            return _out_of_date(409)
        with state:
            PENDING[kind] = None
            state[note], state[bad] = msg, True
        return web.json_response({"error": msg}, status=status)

    too_big = (f"The {label} file is over the {cap / 2 ** 20:.0f} MB limit here. "
               f"Nothing was stored." + (_STL_HINT if kind == "surface" else ""))
    # Refused from the declared length BEFORE a byte is read, and again while
    # streaming in case the declaration lied. aiohttp's read() would cap at 1 MiB, so
    # the body is streamed rather than read.
    if request.content_length is not None and request.content_length > cap:
        return refuse(too_big, 413)
    buf = bytearray()
    async for chunk in request.content.iter_chunked(1 << 16):
        if gen != _gen():
            return _out_of_date(409)
        if len(buf) + len(chunk) > cap:
            return refuse(too_big, 413)
        buf += chunk
    loop = asyncio.get_event_loop()
    checked, why = await loop.run_in_executor(None, _check_upload, kind, name,
                                              bytes(buf))
    del buf
    if gen != _gen():
        return _out_of_date(409)
    if why:
        return refuse(why, 422)
    # On the loop, with no await from here on: the check above still holds.
    PENDING[kind] = checked
    with state:
        state[note], state[bad] = _upload_note(kind), False
        if not (PENDING["surface"] and PENDING["centerline"]):
            other = "centerline" if kind == "surface" else "surface mesh"
            _say(f"Received the {label}. Now add the {other}.")
    if PENDING["surface"] and PENDING["centerline"]:
        # The pair is handed over as it is NOW, so the load cannot pick up a file
        # that arrives -- or is wiped -- while it is storing the first one. And it is
        # numbered: a second pair completed while the first is still building makes
        # the first stale, so only the newest can reach the screen.
        LOAD["seq"] += 1
        asynchronous.create_task(_load_uploads(
            gen, LOAD["seq"], PENDING["surface"], PENDING["centerline"]))
    return web.json_response({"ok": True})


# The newest load of this visit. Two can be in flight at once -- a file replaced while
# the last pair is still building -- and whichever finished last used to be installed,
# with the other's submission id beside it, so a download filed one case's cut under
# the other case's submission.
LOAD = {"seq": 0}


def _out_of_date(status):
    return web.json_response(
        {"error": "This page is out of date. Reload it and try again."}, status=status)


async def _http_download(request):
    """GET /vc/download/{tok}: a one-time link, so a large zip never blocks the page."""
    item = HTTP_DL.pop(request.match_info.get("tok", ""), None)
    if item is None:
        return web.Response(status=404, text="This download link has already been "
                                             "used. Press Download again.")
    name, data, mime = item
    return web.Response(body=data, content_type=mime, headers={
        "Content-Disposition": f'attachment; filename="{name}"',
        "Cache-Control": "no-store"})


@ctrl.add("on_server_bind")
def _bind_http(wslink_server):
    # trame's aiohttp app, reached through its bind hook. aiohttp only -- which is
    # the backend this app runs on.
    wslink_server.app.router.add_post("/vc/upload/{kind}", _http_upload)
    wslink_server.app.router.add_get("/vc/download/{tok}", _http_download)
    if GATE is not None:
        GATE.attach(wslink_server)


async def _load_uploads(gen, seq, surf, cl):
    """Store the pair, build the session OFF the event loop, then draw it on it.

    Building is the slow part -- 18-30 s for a 600k-face mesh, measured -- and it is
    pure numpy, scipy and trimesh, so it runs in a worker
    thread and the page stays live. Drawing touches the VTK scene and stays on the
    loop: VTK work beside the live render window, off its thread, has crashed this
    process before (see shape_worker). A .vtp/.vtk/.vtu CENTERLINE is read through
    VTK, so that one case is built on the loop as before.

    Every await is followed by a check that this load still matters: that the visitor
    who sent the pair is still the one on the page (gen), and that no newer pair has
    replaced it (seq). A worker thread cannot be stopped, so a session built for
    someone who has left, or for a pair since replaced, is thrown away, never shown.
    """
    def outdated():
        return gen != _gen() or seq != LOAD["seq"]

    if outdated():
        return
    loop = asyncio.get_event_loop()
    with state:
        state.busy = True
        _say(f"Loading {surf.name}. The page stays usable while it does.")
    try:
        sub = await loop.run_in_executor(
            None, STORE.create_submission, (surf.name, surf.data), (cl.name, cl.data))
    except UploadRejected as e:
        if outdated():
            return
        with state:
            state.busy = False
            _say(str(e), bad=True)
        return
    except StorageError as e:
        print(f"[uploads] {e}", flush=True)
        if outdated():
            return
        with state:
            state.busy = False
            _say("The upload could not be stored, so nothing was loaded.", bad=True)
        return
    if outdated():
        # Stored, but its id never reached the page -- the visitor left, or replaced
        # the pair first -- so nobody can quote it to have it removed. This line is the
        # only record tying the two together. The id only, never a filename.
        why = "its visitor had left" if gen != _gen() else "a newer upload replaced it"
        print(f"[uploads] {sub.id}: stored after {why}", flush=True)
        return
    with state:
        # Shown now, because it is how the person finds their files and a build can
        # take a minute. SESSION["sub"] -- which a download files its cut under --
        # changes only at install, together with the session it belongs to.
        state.submission_id = sub.id
        # Same reason on_source_changed gives: a declared anatomy belongs to the
        # previous case, and the family call decides which openings are trunk ends.
        if state.family and state.family != "auto":
            state.family = "auto"
        args = _session_args()
    await _build_and_install(
        outdated, str(sub.centerline.path), str(sub.surface.path),
        sub.centerline.extension in (".vtp", ".vtk", ".vtu"), args, _used(args),
        sub_id=sub.id, name=Path(surf.name).stem)


async def _build_and_install(outdated, csv, mat, vtk_centerline, args, used,
                             sub_id=None, name=None, keep_sub=False):
    """Build a session in a worker thread (on the loop for a VTK centerline, which is
    read through VTK) and install it, unless `outdated()` says it no longer matters.

    The tail of an upload's load, and all of "Apply to this case". A rebuild keeps the
    submission it came from even when it fails, so it can be applied again.
    """
    loop = asyncio.get_event_loop()
    try:
        if vtk_centerline:
            s = CutSession(csv, mat, **args)
        else:
            s = await loop.run_in_executor(None, lambda: CutSession(csv, mat, **args))
    except Exception as e:                                        # noqa: BLE001
        if outdated():
            return
        with state:
            SESSION["s"], state.loaded, state.cut_ok = None, False, False
            if not keep_sub:
                SESSION["sub"] = None
            state.busy = False
            _say(_scrub(e if isinstance(e, CutError) else f"{type(e).__name__}: {e}"),
                 bad=True)
        return
    if outdated():
        return
    with state:
        try:
            if sub_id is not None:
                SESSION["sub"] = sub_id
            _install_session(s, used, name=name)
        finally:
            state.busy = False


@ctrl.set("reload_case")
def reload_case():
    """Build the open case again with the Anatomy and Smooth boxes as they are now.

    Both are read only when a case loads, so changing them afterwards did nothing --
    and on the hosted site nothing could make them count, because every upload starts
    from auto. A thoracic case misread as abdominal could never reach nub-in or the
    zones there. This rebuilds from the files the case came from: the stored upload,
    or on this machine the paths it was opened from. Never from a path in state.
    """
    if state.busy:
        return
    args = _session_args()
    sub = SESSION.get("sub")
    if sub:
        try:
            p = STORE.upload_paths(sub)
        except Exception:                                           # noqa: BLE001
            _say("The stored upload could not be read again, so nothing changed.",
                 bad=True)
            return
        csv, mat = str(p["centerline"]), str(p["surface"])
        LOAD["seq"] += 1
        gen, seq = _gen(), LOAD["seq"]
        state.busy = True
        _say("Applying the new settings to this case.")
        asynchronous.create_task(_build_and_install(
            lambda: gen != _gen() or seq != LOAD["seq"], csv, mat,
            Path(csv).suffix.lower() in (".vtp", ".vtk", ".vtu"), args, _used(args),
            sub_id=sub, keep_sub=True))
    elif not PUBLIC and LOADED_FROM["mat"]:
        _open_session(LOADED_FROM["csv"], LOADED_FROM["mat"])
    else:
        _say("There is nothing to read again. Load a case first.", bad=True)


# ── one visitor at a time, and nothing of theirs left for the next ───────────
# Public mode only. See visit_gate.py for why trame needs this and how it works: one
# browser is admitted per process, everyone else is refused before they receive any
# state, and the moment the admitted browser's socket closes, everything it could have
# left behind is wiped here. On this machine nothing changes -- tabs share one state,
# as they always have, because they are all yours.
#
# Shown, by the refused page, in place of the app. Static on purpose: it goes in with
# innerHTML, so nothing a visitor supplied may ever be put into it.
BUSY_HTML = (
    "<style>"
    "@font-face{{font-family:\"Figtree\";font-style:normal;font-weight:300 900;"
    "font-display:swap;src:url(\"vc_static/fonts/figtree-latin-wght-normal.woff2\") "
    "format(\"woff2\")}}"
    "html,body{{margin:0;background:{background}}}"
    "</style>"
    "<div style='font:15px/1.55 Figtree,\"Segoe UI\",system-ui,-apple-system,sans-serif;"
    "max-width:34em;margin:14vh auto;padding:0 20px;color:{on-surface}'>"
    "<div style='font-weight:700;font-size:17px;letter-spacing:-.01em;color:{on-surface};"
    "margin:0 0 22px'>Vascular<span style='color:{primary}'>Cut</span></div>"
    "<h2 style='font-size:20px;font-weight:650;margin:0 0 8px'>This tool is in use</h2>"
    "<p style='margin:0 0 20px;color:{ink2}'>One person uses it at a time, so that nobody "
    "can see anyone else&#8217;s files. Please try again in a few minutes.</p>"
    "<a href='' style='display:inline-block;height:36px;line-height:36px;padding:0 16px;"
    "border-radius:10px;background:linear-gradient(180deg,#6C5DF2,#5644E4);"
    "color:{on-primary};text-decoration:none;font-size:14px;font-weight:600;"
    "box-shadow:0 8px 18px -6px rgba(90,72,230,.55)'>Try again</a></div>"
).format_map(UI_COLOURS)


def _gen():
    """Which visit this is. Changes on every admit and every leave; 0 on this machine."""
    return GATE.gen if GATE is not None else 0


# The state keys a page may write: the app's own, plus the layout's drawer toggle.
# Never trame's own trame__* keys -- the page's template, its script and style lists,
# the Vue plugins it evaluates -- which is what the next visitor's page is built from.
PAGE_KEYS = frozenset(VISIT_DEFAULTS) | {"main_drawer"}


def _page_may_write(key):
    return not key.startswith("trame__") and key in PAGE_KEYS


# The whole state as the first visitor is about to receive it -- trame's keys
# included -- taken once, at the first admit, when no page has been able to write
# anything yet. The wipe restores all of it, not just the app's keys, so a write that
# got past the filter above still does not outlive its visit.
_SNAPSHOT = {"full": None}


def _state_now():
    """Every key trame holds, with its value, pushed to a page or not yet.

    Not state.to_dict(): that returns only what has already been PUSHED, which until
    the server is serving is nothing at all -- a snapshot taken through it was empty
    and the wipe restored nothing, which is how tests/test_visit_wipe.py caught it.
    Reads trame-server's two dicts directly; it is pinned in requirements.txt.
    """
    keys = set(state._pushed_state) | set(state._pending_update)
    return {k: state[k] for k in keys}


def _wipe_visit():
    """Everything a visitor leaves behind, gone before anyone else is let in.

    Order matters. The globals go first, so the @state.change listeners that fire when
    the state is reset find no session and do nothing. The scene is emptied IN PLACE --
    the browser's view is bound to this render window -- and rendered, so no frame of
    the old surface is left in its buffer, and the camera goes back to where it
    started. Then the state goes back to a fresh process's, inside `with state:` so it
    is flushed: trame caches what it last pushed, and a stale cache would suppress any
    update for the next visitor that happened to equal the last one's value.

    Nothing is deleted from disk. Stored submissions are kept by policy; a delta K
    run's working files stay in its temp folder, which nothing on the page can reach.
    """
    PENDING.update(surface=None, centerline=None)
    HTTP_DL.clear()
    SESSION.update(s=None, sub=None)
    LOADED_FROM.update(csv=None, mat=None)
    ZONES.update(rgb=None, legend=[], for_cut=None)
    _FRAMED["src"] = None
    JOB.reset()
    HEART.update(s=None)
    SCENE.clear()
    SCENE.plotter.mesh = None           # pyvista keeps the last dataset it was given
    SCENE.plotter.mapper = None
    SCENE.plotter.camera_position = INITIAL_CAMERA
    SCENE.render()
    full = _SNAPSHOT["full"]
    with state:
        if full is not None:
            # Keys no fresh page has, blanked; then every key as it first was.
            state.update({k: None for k in _state_now() if k not in full})
            state.update(copy.deepcopy(full))
        else:
            state.update(copy.deepcopy(VISIT_DEFAULTS))
    if ctrl.view_update:
        try:
            ctrl.view_update()
        except Exception:                                           # noqa: BLE001
            pass                        # nobody is connected to push a frame to


def _admit_visit():
    if _SNAPSHOT["full"] is None:
        _SNAPSHOT["full"] = copy.deepcopy(_state_now())
    # Wiped again on the way in, not only on the way out: if the last leave's wipe
    # failed part-way, this one either succeeds or raises, and the gate then refuses
    # the visitor rather than show them what was left.
    _wipe_visit()
    _set_token(secrets.token_urlsafe(24))


def _end_visit():
    _set_token(None)
    _wipe_visit()


GATE = None
if PUBLIC:
    GATE = VisitGate(
        busy_html=BUSY_HTML, on_admit=_admit_visit, on_leave=_end_visit,
        may_write=_page_may_write, log=lambda m: print(m, flush=True),
        # Only for a deployment under a path prefix: the prefix, in lower case, e.g.
        # ws://USE_HOST/vc/ws. Untested behind a proxy; see tools/HOSTING.md.
        session_url=os.environ.get("VC_SESSION_URL") or "ws://USE_HOST/ws",
        # An open tab that sends nothing for this long is closed, so a forgotten tab
        # cannot hold the instance for ever. 0 turns it off.
        idle_limit=60.0 * float(os.environ.get("VC_VISIT_IDLE_MIN") or 30))
    server.add_protocol_to_configure(lambda p: p.registerLinkProtocol(GATE.protocol()))


# ── depth ────────────────────────────────────────────────────────────────────
_LAST_PREVIEW = {"t": 0.0}


def _throttled(min_interval=0.04):
    """Drop preview frames that arrive faster than the renderer can serve them.

    Every preview pushes a freshly rendered image down the socket. Un-throttled, a
    drag queues far more frames than can be drawn, the channel backs up, and the
    whole UI -- including the slider thumb -- goes sticky.
    """
    now = time.monotonic()
    if now - _LAST_PREVIEW["t"] < min_interval:
        return False
    _LAST_PREVIEW["t"] = now
    return True


@server.trigger("preview_offset")
def preview_offset(index, value=None):
    """While dragging: move the disc only. A full re-cut per frame is far too slow;
    the desktop window makes the same split and it is why dragging feels live."""
    s = SESSION["s"]
    row = _row(index)
    if s is None or row is None or not _throttled():
        return
    depth = float(row["offset"] if value is None else value)
    try:
        t = s.terminal(row["label"])
        origin, normal = t.cut_anchor(depth)
    except Exception:
        return
    planes = planes_for(s, row["label"])
    for p in planes:
        if p["label"] == row["label"]:
            p["origin"] = np.asarray(origin, float)
            p["normal"] = np.asarray(normal, float)
    SCENE.set_planes(planes)
    if ctrl.view_update:
        ctrl.view_update()


def _sync_marks(marks=None):
    """One line describing the pins, composed here rather than in a template.

    Takes the list rather than reading state back: a value just written is not
    reliably the value that comes out of the reactive proxy in the same handler,
    and the note silently stayed empty while the pin itself worked.
    """
    marks = list(state.marks or []) if marks is None else list(marks)
    if not marks:
        state.marks_note = ""
        return
    last = marks[-1]
    state.marks_note = (f"{len(marks)} mark(s), last: {last['vessel']} "
                        f"{last['mm_to_node']:.0f} mm to node")


def _pick_world(x, y):
    """World point under a click in the 3-D view, or None if it missed the mesh.

    Only the mesh is pickable -- the discs, the centerline and the marks are all
    set unpickable -- so a click either lands on the wall or on nothing.
    """
    import vtk
    ren = SCENE.plotter.renderer
    picker = vtk.vtkCellPicker()
    picker.SetTolerance(0.002)
    if not picker.Pick(float(x), float(y), 0.0, ren):
        return None
    return list(picker.GetPickPosition())


@server.trigger("view_click")
def view_click(pos=None):
    """Pin a mark where the user clicked on the surface.

    The desktop measures under the cursor and pins on click. Here it is click only:
    every hover would be a round trip to the server and back, which would feel worse
    than not having it. The measurement itself is identical.
    """
    s = SESSION["s"]
    if s is None or not pos:
        return
    size = SCENE.plotter.ren_win.GetSize()
    # vtk.js reports the position in display pixels with the origin at bottom-left,
    # which is what vtkCellPicker wants -- but it arrives scaled to the client's
    # canvas, so anything outside the window is a miss rather than a wild pick.
    x, y = float(pos[0]), float(pos[1])
    if not (0 <= x <= size[0] and 0 <= y <= size[1]):
        return
    world = _pick_world(x, y)
    if world is None:
        return
    m = s.measure_at(world)
    if not m:
        return
    marks = list(state.marks or [])
    marks.append(m)
    state.marks = marks
    _sync_marks(marks)
    SCENE.set_marks([mk["point"] for mk in marks])
    _say(f"{m['vessel']}: {m['mm_to_node']:.0f} mm to node"
         f"{' (from the vessel end)' if m['from_vessel_end'] else ''}.")
    if ctrl.view_update:
        ctrl.view_update()


@server.trigger("clear_marks")
def clear_marks():
    state.marks = []
    _sync_marks([])
    SCENE.set_marks([])
    if ctrl.view_update:
        ctrl.view_update()


@server.trigger("send_cut_to_mark")
def send_cut_to_mark():
    """Put a cut where the last mark is.

    Pinning a point and then reading its number off the screen to type into a
    slider is the long way round; this is what the pin was for.
    """
    s = SESSION["s"]
    if s is None or not state.marks:
        return
    plan = s.offset_for_mark(state.marks[-1])
    if plan is None:
        _say("That mark is not on a vessel being cut, so there is no cut to move.",
             bad=True)
        return
    s.set_offset(plan["label"], plan["offset_mm"])
    state.selected = plan["label"]
    _recut_and_draw()
    if s.ok:
        _say(f"{plan['name']} cut moved to {plan['offset_mm']:.0f} mm, "
             f"where the mark is.")


@server.trigger("nudge_branches")
def nudge_branches(delta=None):
    """Move every branch cut together, without touching the trunk ends.

    Setting eight branches one at a time to try a slightly shallower trim is the
    kind of thing that makes a tool feel like work. Trunk ends are excluded on
    purpose: they are the aorta's own two cuts and have nothing to do with how far
    a branch is trimmed.
    """
    s = SESSION["s"]
    if s is None or delta is None:
        return
    d = float(delta)
    rows = [o for o in s.openings()
            if not o["is_trunk"] and not o["is_region_bound"]]
    if not rows:
        _say("No branches to move, this vessel has only trunk ends.", bad=True)
        return
    moved = 0
    for o in rows:
        before = float(o["offset_mm"])
        after = min(max(before + d, MIN_OFFSET_MM), float(o["limit_mm"]))
        if abs(after - before) > 1e-6:
            s.set_offset(o["label"], after)
            moved += 1
    if not moved:
        _say("Branches are already at the end of their range, nothing to move.",
             bad=True)
        return
    _recut_and_draw()
    if s.ok:
        _say(f"Moved {moved} branch cut(s) "
             f"{abs(d):.0f} mm {'deeper' if d > 0 else 'toward the tips'}.")


@server.trigger("declare_artery")
def declare_artery(index=None, artery=None):
    """Say which artery an opening IS -- which is what drives the region presets.

    Renaming an opening is a label and nothing more. Declaring it is the fix for the
    half of the abdominal cohort where no rule can tell which lateral pair the renals
    are, and where the presets therefore collapse to "Whole abdominal aorta".
    """
    s = SESSION["s"]
    if s is None or index is None:
        return
    try:
        row = state.openings[int(index)]
    except (IndexError, TypeError, ValueError):
        return
    try:
        s.declare_artery(row["label"], artery or None)
    except Exception as e:                                      # noqa: BLE001
        _say(str(e), bad=True)
        return
    _sync_openings()
    _sync_anatomy()
    _say(f"{row['label']} declared as {artery}." if artery else
         f"{row['label']} handed back to the geometry.")


@server.trigger("nub_preset")
def nub_preset():
    """Cut every branch back to a short stub and square off the aorta's ends.

    The notes are shown rather than swallowed: a nub that would not cut cleanly is
    moved outward until it does, and a trunk end too oblique for a horizontal plane
    falls back to perpendicular. Both are cases for a person to look at, and neither
    is visible in the finished geometry.
    """
    s = SESSION["s"]
    if s is None:
        return
    try:
        notes = s.nub_preset()
    except Exception as e:                                      # noqa: BLE001
        _say(str(e), bad=True)
        return
    _sync_openings()
    _refresh_status()
    _redraw(full=True)
    # The headline on the status line; every rescue listed under Details, where it
    # can be read, instead of three of them run together and the rest cut off.
    head, rest = notes[0], notes[1:]
    _say(head + (f"  {len(rest)} note(s) under Details." if rest else ""), notes=rest)


@server.trigger("add_trim")
def add_trim(segment=None):
    """One cut on the vessel feeding a group, instead of one per opening."""
    s = SESSION["s"]
    if s is None or not segment:
        return
    try:
        s.add_trim(str(segment))
    except Exception as e:                                      # noqa: BLE001
        _say(str(e), bad=True)
        return
    state.trim_pick = None
    _recut_and_draw()


@server.trigger("remove_trim")
def remove_trim(segment=None):
    """Undo a trim: the openings it swallowed come back as individual cuts."""
    s = SESSION["s"]
    if s is None or not segment:
        return
    try:
        s.remove_trim(str(segment))
    except Exception as e:                                      # noqa: BLE001
        _say(str(e), bad=True)
        return
    _recut_and_draw()


@server.trigger("pull_back")
def pull_back(index):
    """Move one opening to the deepest offset that still cuts cleanly."""
    s, row = SESSION["s"], _row(index)
    if s is None or row is None:
        return
    state.busy = True
    try:
        off, note = s.pull_back(row["label"])
    except Exception as e:                                      # noqa: BLE001
        state.busy = False
        _say(str(e), bad=True)
        return
    finally:
        state.busy = False
    _put(index, offset=float(off))
    state.selected = row["label"]
    _recut_and_draw()
    if s.ok:
        _say(note)          # the recut's own "Ready" would otherwise bury the number


@server.trigger("commit_offset")
def commit_offset(index, value=None):
    s, row = SESSION["s"], _row(index)
    if s is None or row is None:
        return
    depth = float(row["offset"] if value is None else value)
    _put(index, offset=depth)
    s.set_offset(row["label"], depth)
    state.selected = row["label"]
    _recut_and_draw()


# ── the plane: tilt, shift, size ─────────────────────────────────────────────
FIELD_TO_KW = {
    "tilt_a": "tilt_a_deg", "tilt_b": "tilt_b_deg",
    "shift_a": "shift_a_mm", "shift_b": "shift_b_mm", "band": "band_mm",
}


def _apply_plane(index, field, value):
    s, row = SESSION["s"], _row(index)
    if s is None or row is None or field not in FIELD_TO_KW:
        return None, None
    row = _put(index, **{field: float(value)}) if value is not None else row
    s.set_plane(row["label"],
                tilt_a_deg=row["tilt_a"], tilt_b_deg=row["tilt_b"],
                shift_a_mm=row["shift_a"], shift_b_mm=row["shift_b"],
                band_mm=row["band"])
    return s, row


@server.trigger("preview_plane")
def preview_plane(index, field, value=None):
    """Tilt and shift move the disc live too -- they are plane adjustments, and
    watching them is the entire point of having them."""
    if not _throttled():
        return
    s, row = _apply_plane(index, field, value)
    if s is None:
        return
    SCENE.set_planes(planes_for(s, row["label"]))
    if ctrl.view_update:
        ctrl.view_update()


@server.trigger("commit_plane")
def commit_plane(index, field, value=None):
    s, row = _apply_plane(index, field, value)
    if s is None:
        return
    state.selected = row["label"]
    _recut_and_draw()


@server.trigger("reset_plane")
def reset_plane(index):
    s, row = SESSION["s"], _row(index)
    if s is None or row is None:
        return
    s.reset_plane(row["label"])
    state.selected = row["label"]
    _recut_and_draw()


@server.trigger("set_mode")
def set_mode(index, value=None):
    s, row = SESSION["s"], _row(index)
    if s is None or row is None:
        return
    mode = value if value is not None else row["mode"]
    try:
        s.set_mode(row["label"], mode)
    except CutError as e:
        _say(str(e), bad=True)
        return
    _put(index, mode=mode)
    state.selected = row["label"]
    _recut_and_draw()


@server.trigger("rename_opening")
def rename_opening(index, value=None):
    s, row = SESSION["s"], _row(index)
    if s is None or row is None:
        return
    name = value if value is not None else row["name"]
    s.rename(row["label"], name)
    _sync_openings()


@server.trigger("select_row")
def select_row(label):
    """Expanding a row also selects it, and the selected row's disc goes darker --
    that is how you tell which disc you are about to move."""
    state.selected = None if state.selected == label else label
    _redraw(full=False)


# ── region ───────────────────────────────────────────────────────────────────
@state.change("region_pick")
def on_region(region_pick, **_):
    s = SESSION["s"]
    if s is None or not state.loaded:
        return
    try:
        if region_pick == "Whole vessel":
            s.clear_region()
            state.in_region = False
            state.region_name = ""
            state.region_why = ""
        else:
            p = s.apply_preset(region_pick)
            state.in_region = True
            state.region_name = p["name"]
            state.region_why = p.get("why", "") or ""
            _sync_region_boxes(s)
        _recut_and_draw()
    except Exception as e:
        _say(f"{type(e).__name__}: {e}", bad=True)


def _sync_region_boxes(s):
    """Put the Start/End boxes where the region ACTUALLY is. Two traps here.

    s1 None is an OPEN distal end -- one plane, everything below it kept. It is what
    "Infrarenal aorta with iliacs" ships (branch_cut.region_presets), because a
    two-plane region's distal plane is clamped a millimetre short of the route end
    while an abdominal aorta's iliac ostia sit AT the bifurcation, so they would fall
    outside it. `round(float(None), 1)` raises TypeError, which is why picking that
    preset used to fail outright on the web: the raise happened BEFORE the redraw, so
    the banner said the region was isolated over a scene still showing the old cut.

    And apply_preset walks a ladder when the station it was asked for will not cut,
    so `s.region` is not always the preset's own s0/s1. Reading the preset dict --
    which is what this used to do -- showed the user a station the cut did not use.
    Read the session, which is the only thing that knows where the planes went.
    """
    reg = getattr(s, "region", None)
    if not reg:
        state.region_open = False
        return
    s0, s1 = reg
    state.region_s0 = round(float(s0), 1)
    state.region_open = s1 is None
    # Left at its previous value when the end is open. It is disabled and unread;
    # writing route_length here would look like a real station sitting in a live box.
    if s1 is not None:
        state.region_s1 = round(float(s1), 1)


@ctrl.set("commit_region")
def commit_region():
    s = SESSION["s"]
    if s is None or not state.in_region:
        return
    # An open end stays open. The End box is disabled for such a region, so its number
    # is whatever was last left there -- 0.0 by default, or the End of whichever
    # two-plane region ran before. Reading it unconditionally is how the DESKTOP lost
    # this region's open end on the first drag of the Start slider, reported from a
    # 49-case cohort and fixed in branch_gui._span_from_controls; the same read was
    # here, unguarded, and the End control was never even disabled.
    s1 = None if state.region_open else float(state.region_s1)
    s.set_region(float(state.region_s0), s1,
                 state.region_name or state.region_pick, manual=True)
    _recut_and_draw()


@ctrl.set("clear_region")
def clear_region():
    state.region_pick = "Whole vessel"


# ── view ─────────────────────────────────────────────────────────────────────
@state.change("show_centerline", "show_planes", "show_loops", "opacity",
              "colour_branches", "colour_zones")
def on_view_toggle(**_):
    if SESSION["s"] is not None:
        _redraw(full=True)


@state.change("clean_rims")
def on_clean(clean_rims, **_):
    s = SESSION["s"]
    if s is None:
        return
    s.clean_rims = bool(clean_rims)
    _recut_and_draw()


# ── reset / export ───────────────────────────────────────────────────────────
@ctrl.set("reset_all")
def reset_all():
    """Back to the model as loaded: region, trims, hand names and planes all dropped.

    This used to clear only the offsets, the planes and the region, so a trimmed
    branch group and a hand-declared artery both survived a Reset -- while the
    desktop's Reset dropped them. The same button did two different things depending
    on which front end you were in, and on the web it left the model visibly not
    reset in exactly the situation the button exists for. The scope now lives in
    CutSession.reset_all so the two cannot drift again.
    """
    s = SESSION["s"]
    if s is None:
        return
    s.reset_all()
    state.region_pick = "Whole vessel"
    state.in_region = False
    state.region_open = False
    state.selected = None
    state.trim_pick = None
    _recut_and_draw()
    SCENE.reset_camera(s.mesh)
    if ctrl.view_reset_camera:
        ctrl.view_reset_camera()


@server.trigger("download_cut")
def download_cut():
    """The cut and its sidecar, zipped, handed to the browser as a download.

    One zip rather than two files: the sidecar is what makes the STL reproducible, a
    single file cannot be half-downloaded, it halves the transfer, and two downloads
    from one click make the browser ask permission. Returns None when there is nothing
    to give -- the page only downloads when it gets bytes back, so a refusal never
    becomes a file containing the word null.
    """
    s = SESSION["s"]
    if s is None:
        state.export_bad, state.export_msg = True, "Nothing is loaded yet."
        return None
    if not s.ok:
        state.export_bad, state.export_msg = True, f"Cut is not valid. {s.problem}"
        return None
    stem = safe_stem(state.export_name)
    stl = s.mesh.export(file_type="stl")
    side = s.sidecar()
    kept = ""
    if SESSION.get("sub"):
        # Retained with the upload it came from, under that submission's id (job 7).
        # The stored sidecar names the files by role -- surface.mat, not the name
        # somebody gave it -- because filenames are where identifiers leak.
        side["submission"] = SESSION["sub"]
        try:
            STORE.add_cut(SESSION["sub"], stl, side)
            kept = " A copy is kept with your upload."
        except StorageError as e:
            print(f"[uploads] {e}", flush=True)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{stem}.stl", stl)
        z.writestr(f"{stem}_cut_planes.json", json.dumps(side, indent=2))
    tok = secrets.token_urlsafe(24)
    HTTP_DL[tok] = (f"{stem}.zip", buf.getvalue(), "application/zip")
    while len(HTTP_DL) > HTTP_DL_KEEP:
        HTTP_DL.pop(next(iter(HTTP_DL)))
    state.export_bad = False
    state.export_msg = (f"Downloaded {stem}.zip: {stem}.stl and "
                        f"{stem}_cut_planes.json.{kept}")
    # A link, not the bytes: a large zip sent back over the websocket would freeze
    # the page for as long as it took to arrive.
    return f"vc/download/{tok}"


@ctrl.set("export_cut")
def export_cut():
    s = SESSION["s"]
    if s is None:
        return
    if PUBLIC:
        # Writing into the server's own output/ helps nobody who is not at this
        # machine, and every visitor would share -- and overwrite -- that one folder.
        state.export_bad = True
        state.export_msg = "This instance does not write files for you. Use Download."
        return
    if not s.ok:
        state.export_bad = True
        state.export_msg = f"Cut is not valid. {s.problem}"
        return
    try:
        stl, side = s.export(ROOT / "output", safe_stem(state.export_name))
        state.export_bad = False
        # Say WHOSE disk. "Wrote ... to output/" reads as "it is on your machine",
        # which is true only when the server is your machine. There is no download
        # yet -- it pairs with the upload work in tools/HOSTING.md, neither built --
        # so a remote user needs to know the file is sitting on the server rather
        # than waiting for a save dialog that is never coming.
        state.export_msg = (
            f"Wrote {stl.name} and {side.name} to output/ ON THE SERVER"
            f"{' (this machine)' if not PUBLIC else ''}."
            + ("" if not PUBLIC else
               " There is no download yet, so ask whoever runs this instance for "
               "the file."))
    except Exception as e:
        state.export_bad = True
        state.export_msg = f"{type(e).__name__}: {e}"


def _upload_field(key, note, bad, label, kind, max_key):
    """One file picker. It sends over HTTP, not the websocket, with its progress shown.

    The size check is done IN THE BROWSER first, so an oversize file is never sent;
    the server checks again from the declared length before reading, and while it
    streams. `accept` goes through __properties because trame's VFileInput wrapper does
    not list it, and an unlisted keyword is dropped without a word.

    The field is bound to a key that is never set, so it clears itself once a file is
    chosen and the same file can be chosen again; what arrived is shown in the note.
    """
    pct = f"{key}_pct"
    route = f"vc/upload/{kind}"
    hint = _STL_HINT.replace("'", "") if kind == "surface" else ""
    # What it takes and how big, said before the file is chosen rather than after it
    # is refused.
    takes = (f"{' '.join(STORE.accepted_extensions(kind))}"
             f" · up to {STORE.limit_mb(kind):.0f} MB")
    test_id = "upload-mesh" if kind == "surface" else "upload-centerline"
    # A drop zone: the whole field opens the file dialog, and a file dropped on it is
    # taken. The label element stays (visually hidden) because it is what names the
    # input for assistive tech and for the upload checks.
    with html.Div(classes="fld vc-drop-fld"):
        html.Label(label, classes="fld-lbl vc-sr")
        with v.VFileInput(
            model_value=(key, None),
            density="compact", variant="outlined", hide_details=True,
            prepend_icon="", show_size=True, classes="vc-drop",
            accept=STORE.accept_attribute(kind), placeholder="Choose a file",
            disabled=("busy",),
            raw_attrs=[f'data-test="{test_id}"'],
            __properties=["accept", "placeholder"],
            __events=["update:modelValue"],
            update_modelValue=(
                f"$event && (($event.size > {max_key}) ? "
                f"({note} = 'That file is ' + ($event.size / 1048576).toFixed(1)"
                f" + ' MB; the limit here is ' + ({max_key} / 1048576).toFixed(0)"
                f" + ' MB. Nothing was sent.{hint}', {bad} = true) : "
                f"(({note} = 'Sending ' + $event.name + '...'), ({bad} = false), "
                f"utils.get('vcUpload')($event, '{route}', vc_token, '{pct}')"
                # Anything but 200 is shown. The app's own refusals (413, 422) also
                # arrive through state, in the same words; the rest -- a stale
                # token's 403 after a restart, a proxy's body-size limit -- used to
                # leave "Sending..." on screen for ever. vc_upload.js hands back an
                # error only when the body was this app's JSON, never a proxy's page.
                f".then(r => {{ if (r && r.status !== 200) {{ {note} = "
                f"(r.body && r.body.error) || ('The upload was refused (HTTP ' + "
                f"r.status + '). Reload the page and try again.'); {bad} = true; }} }})))")):
            # In the field, because a file input never shows its placeholder. Always
            # true: the field clears itself once a file is chosen (see above).
            with html.Template(v_slot_prepend_inner=True):
                with html.Span(classes="vc-drop-ic"):
                    v.VIcon("mdi-tray-arrow-up", size="21")
                with html.Span(classes="vc-drop-tx", raw_attrs=['aria-hidden="true"']):
                    html.Span(label, classes="vc-drop-t")
                    html.Span("Choose a file…", classes="vc-drop-a")
                    with html.Span(classes="vc-drop-m"):
                        html.Span(f"{' '.join(STORE.accepted_extensions(kind))} · ")
                        html.Span(f"up to {STORE.limit_mb(kind):.0f} MB",
                                  style="white-space:nowrap")
        # !== -1, not >= 0: a bare > in a bound attribute closes the tag when the
        # template is built, and the element silently vanishes.
        v.VProgressLinear(model_value=(pct, -1), v_if=(f"{pct} !== -1",),
                          height=3, color="primary", classes="mt-1")
        with html.Div(v_if=(note,), classes=(f"{bad} ? 'upnote bad' : 'upnote'",)):
            v.VIcon("mdi-alert-circle-outline", v_if=(bad,), size="14",
                    classes="mr-1")
            html.Span(f"{{{{ {note} }}}}")


# ── shape metrics: delta K over the CAMP partition ───────────────────────────
def _sm_say(text, bad=False):
    state.sm_msg = text
    state.sm_bad = bool(bad)


def _draw_shape():
    """Whatever the shape section should be showing right now.

    Three states, and each of them draws something: no surface, a surface waiting to
    be measured, and a partition. The middle one matters -- a blank viewport while you
    pick the settings reads as a page that failed to load.
    """
    verts, faces = JOB.polydata_arrays()
    if verts is None:
        SCENE.clear_patches()
        state.sm_range = state.sm_full = ""
    elif not JOB.done:
        SCENE.set_patches(verts, faces)
        state.sm_range = state.sm_full = ""
    else:
        fv = JOB.face_values(scale=state.sm_scale, lo=state.sm_cmin,
                             hi=state.sm_cmax, clip_pct=state.sm_clip)
        SCENE.set_patches(
            verts, faces, scalars=fv["values"], clim=fv["clim"], cmap=fv["cmap"],
            # The colour bar is the page's (the key over the view), not the scene's:
            # drawn in the render as well, it sat beside it and said the same thing.
            n_colors=fv["n_colors"], title=fv["title"], show_bar=False,
            kept_mask=JOB.kept_face_mask() if state.sm_dim_trimmed else None,
            edges=JOB.patch_edges() if state.sm_outlines else None)
        lo, hi = fv["range"]
        state.sm_full = f"{lo:.4g} to {hi:.4g}" if fv["bar"] else ""
        # Keep the manual boxes showing what is actually on the bar, so switching to
        # Manual starts from the current picture instead of from whatever was typed
        # against a different colouring.
        if state.sm_scale != "manual":
            state.sm_cmin = round(float(fv["clim"][0]), 6)
            state.sm_cmax = round(float(fv["clim"][1]), 6)
        if not fv["bar"]:
            state.sm_range = ""
        elif fv["clipped"]:
            hidden = float((fv["values"] < fv["clim"][0]).mean()
                           + (fv["values"] > fv["clim"][1]).mean())
            state.sm_range = (f"Bar {fv['clim'][0]:.4g} to {fv['clim'][1]:.4g}. "
                              f"{100 * hidden:.1f}% of the surface is past an end and "
                              f"drawn at it. Full range {lo:.4g} to {hi:.4g}.")
        else:
            state.sm_range = (f"Bar spans the whole range, {lo:.4g} to {hi:.4g}, "
                              f"which is how CAMP scales it.")

    src = str(JOB.source) if JOB.source is not None else None
    if src != _FRAMED["src"]:
        _FRAMED["src"] = src
        SCENE.reset_camera(mesh=JOB.source_mesh)
        if ctrl.view_reset_camera:
            ctrl.view_reset_camera()
    if ctrl.view_update:
        ctrl.view_update()


def _sm_after_source():
    state.sm_info = JOB.source_info()
    state.sm_result = None
    state.sm_notes = []
    state.sm_headline = ""
    state.sm_export_name = safe_stem(Path(JOB.source_name).stem, "delta_k")
    state.sm_export_msg = ""
    _sm_say(f"{JOB.source_name}: {JOB.source_faces:,} faces, "
            f"{2 - JOB.source_euler} opening(s). Ready to measure.")
    if state.mode == "shape":
        _draw_shape()


@ctrl.set("sm_use_cut")
def sm_use_cut():
    """Measure the surface this session just cut, without a round trip through disk.

    Refuses an invalid cut rather than measuring it. delta K on a surface with a hole
    where a cut punched through two walls is a number, and it is a number about the
    defect, not about the anatomy.
    """
    s = SESSION["s"]
    if s is None or s.mesh is None:
        _sm_say("Nothing has been cut yet. Load a case in the cut section first.",
                bad=True)
        return False
    if not s.ok:
        _sm_say(f"The cut is not valid, so measuring it would not mean anything. "
                f"{s.problem}", bad=True)
        return False
    try:
        JOB.adopt_mesh(s.mesh, safe_stem(state.export_name),
                       expect_euler=int(s.mesh.euler_number))
    except ShapeError as e:
        _sm_say(str(e), bad=True)
        return False
    except Exception as e:                                  # noqa: BLE001
        _sm_say(f"{type(e).__name__}: {e}", bad=True)
        return False
    _sm_after_source()
    return True


@ctrl.set("sm_handoff")
def sm_handoff():
    """The cut section's way into the shape section."""
    if sm_use_cut():
        state.mode = "shape"


@ctrl.set("sm_load_path")
def sm_load_path():
    if PUBLIC:
        _sm_say("This instance takes uploaded files, not paths. Upload is not built "
                "yet, so loading is unavailable here.", bad=True)
        return
    try:
        JOB.adopt_path(state.sm_path)
    except ShapeError as e:
        _sm_say(str(e), bad=True)
        return
    _sm_after_source()


def _sm_band_pct(s):
    """How far this delta K could move if the partition were drawn again, per cent.

    Over one partition that is not an estimate from this run -- a single draw cannot
    report its own spread -- it is the cohort figure, PARTITION_CV, measured over 30
    aortas at 8 seeds each. Over several it is this run's own standard error, which
    is the better number because it is measured on the surface in front of you.

    Returned so the UI can print it in the SAME string as the value. A delta K and
    its band get separated the moment they are in different elements: the value ends
    up quoted and the band does not travel with it.
    """
    n_seeds = int(s.get("n_seeds", 1) or 1)
    se, k = s.get("k_fluct_se"), float(s.get("k_fluct") or 0.0)
    if n_seeds > 1 and se and k:
        return round(100.0 * float(se) / k)
    return round(100.0 * SM.partition_band(n_seeds))


def _sm_finished():
    s = JOB.summary
    state.sm_result = {
        "k_fluct": s["k_fluct"],
        "band_pct": _sm_band_pct(s),
        "n_seeds": int(s.get("n_seeds", 1) or 1),
        "A_mm2": s["A_mm2"],
        "n_patches": s["n_patches"],
        "n_kept": s["n_kept"],
        "pct_floor": s["pct_at_min_size"],
        "pct_failed": s["pct_failed_quality"],
        "faces": s["faces"],
        "euler": s["euler_number"],
        "seconds": s.get("seconds", 0.0),
    }
    state.sm_notes = [{"level": lv, "text": t} for lv, t in JOB.notes()]
    state.sm_headline = JOB.headline()
    _sm_say(f"delta K = {s['k_fluct']:.5g} +/- {_sm_band_pct(s)}% over "
            f"{s['n_patches']:,} patches, in {s.get('seconds', 0):.0f} s.")
    _draw_shape()


async def _sm_watch(gen, proc):
    """Poll the worker from the event loop.

    Polling, not a callback, because the work is in another process. It is in another
    process because ACVD and CAMP are VTK underneath and running them on a thread
    beside this app's live render window segfaults the interpreter. See shape_worker.

    It watches ONE run for ONE visit, both fixed when it is scheduled. Once the
    visitor has left, or this run has been cancelled or replaced, it stops without
    writing anything -- otherwise it would report a result, or "idle", onto the next
    person's page, or attach itself to their run.
    """
    while True:
        await asyncio.sleep(0.35)
        if gen != _gen() or JOB.proc is not proc:
            return
        with state:
            st = JOB.poll()
            state.sm_elapsed = round(float(st.get("elapsed", 0.0)), 1)
            if st["state"] == "running":
                state.sm_phase = st["message"]
                continue
            state.sm_busy = False
            state.sm_phase = ""
            if st["state"] == "done":
                _sm_finished()
            elif st["state"] == "failed":
                _sm_say(st["message"], bad=True)
            # "idle" means it was cancelled from the panel, which already said so.
            return


@ctrl.set("sm_run")
def sm_run():
    if JOB.source is None:
        _sm_say("Choose a surface first.", bad=True)
        return
    if state.sm_busy:
        return
    try:
        JOB.start(int(state.sm_faces), float(state.sm_trim),
                  n_seeds=int(state.sm_seeds or 1))
    except ShapeError as e:
        _sm_say(str(e), bad=True)
        return
    state.sm_busy = True
    state.sm_result = None
    state.sm_notes = []
    state.sm_headline = ""
    state.sm_elapsed = 0.0
    state.sm_phase = "Starting"
    _sm_say(f"Measuring {JOB.source_name} at {int(state.sm_faces):,} faces.")
    # The previous partition is not this one. Leaving it up while a new run goes would
    # let someone read a colour off the old surface and a number off the new panel.
    _draw_shape()
    asynchronous.create_task(_sm_watch(_gen(), JOB.proc))


@ctrl.set("sm_cancel")
def sm_cancel():
    JOB.cancel()
    state.sm_busy = False
    state.sm_phase = ""
    _sm_say("Cancelled.", bad=True)


@ctrl.set("sm_export")
def sm_export():
    if PUBLIC:
        # As export_cut: the server's output/ is shared by every visitor and read by
        # none of them, and a table named after one person's cut would sit there for
        # the next to overwrite. The button is not rendered in public mode; this is
        # the check that holds when a client fires the trigger anyway.
        state.sm_export_bad = True
        state.sm_export_msg = "This instance does not write files for you."
        return
    if not JOB.done:
        state.sm_export_bad = True
        state.sm_export_msg = "There is no result to export yet."
        return
    try:
        table, side = JOB.export(ROOT / "output",
                                 safe_stem(state.sm_export_name, "delta_k"))
        state.sm_export_bad = False
        state.sm_export_msg = f"Wrote {table.name} and {side.name} to output/"
    except Exception as e:                                  # noqa: BLE001
        state.sm_export_bad = True
        state.sm_export_msg = f"{type(e).__name__}: {e}"


@state.change("sm_outlines", "sm_dim_trimmed", "sm_scale", "sm_clip")
def on_sm_view(**_):
    if state.mode == "shape" and JOB.done:
        _draw_shape()


@state.change("sm_trim")
def on_sm_trim(sm_trim=None, **_):
    """Restate the repeatability the moment the trim moves.

    The band belongs to the SETTING, not to the result, so it has to be right before
    anything is computed -- that is when someone decides whether trim 0 is what they
    want. Waiting until the run finished would show the number only once it was too
    late to act on.
    """
    try:
        band, at = SM.repeatability_band(float(sm_trim))
    except (TypeError, ValueError):
        return
    state.sm_repeat = _repeat_text(band, at)


@ctrl.set("sm_apply_bounds")
def sm_apply_bounds():
    """Manual colour bounds, applied when a box is committed rather than per keypress."""
    if state.sm_scale != "manual":
        state.sm_scale = "manual"      # the change handler redraws
        return
    if state.mode == "shape" and JOB.done:
        _draw_shape()


# ── cardiac: split a ventricle, cut its base, hand the walls to Shape ────────
def _hs_say(text, bad=False):
    state.hs_msg = text
    state.hs_bad = bool(bad)


def _hs_draw(frame=False):
    """The cardiac view: both walls, the removed cap, the plane."""
    s = HEART["s"]
    if s is None:
        SCENE.clear_heart()
    else:
        SCENE.set_heart(s.inner, s.outer, s.removed,
                        show={"inner": state.hs_show_inner,
                              "outer": state.hs_show_outer,
                              "removed": state.hs_show_removed})
        SCENE.set_heart_plane(s.disc() if state.hs_show_plane else None)
        if frame or HEART.get("framed") is not s:
            HEART["framed"] = s
            SCENE.reset_camera(mesh=s.main)
            if ctrl.view_reset_camera:
                ctrl.view_reset_camera()
    if ctrl.view_update:
        ctrl.view_update()


def _hs_publish(s):
    """The session's verdict and readout into state, and one line for the footer."""
    state.hs_ok = bool(s.ok)
    state.hs_problem = s.problem or ""
    state.hs_info = {k: v for k, v in s.info.items() if k not in ("notes", "problems")}
    state.hs_notes = list(s.info.get("notes", []))
    if s.ok:
        i = s.info
        how = "plane alone" if i["separated"] else "plane + automatic seam"
        _hs_say(f"Cut valid ({how}): {i['removed_pct']:.0f}% of the outer wall removed, "
                f"{i['holes_open']} hole(s) open, {i['holes_filled']} patched.")
    else:
        _hs_say(f"Cut not valid. {s.problem}", bad=True)


async def _hs_build(outdated, path):
    loop = asyncio.get_event_loop()
    fill = float(state.hs_fill_mm) if state.hs_fill else 0.0
    try:
        # Seeds and the geodesic seam: seconds on a 300k-face shell, so off the loop.
        # No VTK is touched here (trimesh and scipy only), so a thread is safe.
        s = await loop.run_in_executor(
            None, lambda: HeartSession(path, depth_mm=HEART_DEPTH_MM, fill_holes_mm=fill))
    except Exception as e:                                          # noqa: BLE001
        if outdated():
            return
        with state:
            HEART["s"] = None
            state.hs_loaded = state.hs_busy = False
            _hs_say(str(e) if isinstance(e, HeartError) else f"{type(e).__name__}: {e}",
                    bad=True)
            _hs_draw()
        return
    if outdated():
        return
    with state:
        HEART["s"] = s
        state.hs_loaded = True
        state.hs_name = s.path.stem
        state.hs_export_name = s.path.stem
        state.hs_depth, state.hs_tilt_a, state.hs_tilt_b = HEART_DEPTH_MM, 0.0, 0.0
        state.hs_export_msg = ""
        _hs_publish(s)
        state.hs_busy = False
        if state.mode == "heart":
            _hs_draw(frame=True)


@ctrl.set("hs_load")
def hs_load():
    if PUBLIC:
        _hs_say("The Cardiac section reads meshes from this machine, so it is "
                "available only in the local version.", bad=True)
        return
    path = str(state.hs_path or "").strip().strip('"')
    if not path:
        _hs_say("Choose a ventricle mesh first.", bad=True)
        return
    HEART["seq"] += 1
    gen, seq = _gen(), HEART["seq"]
    state.hs_busy = True
    _hs_say(f"Splitting {Path(path).name} into inner and outer walls…")
    asynchronous.create_task(_hs_build(lambda: gen != _gen() or seq != HEART["seq"],
                                       path))


def _hs_wanted():
    fill = float(state.hs_fill_mm) if state.hs_fill else 0.0
    return (float(state.hs_depth), float(state.hs_tilt_a), float(state.hs_tilt_b), fill)


async def _hs_recut(outdated, wanted):
    loop = asyncio.get_event_loop()
    s = HEART["s"]
    depth, ta, tb, fill = wanted
    s.set_plane(depth_mm=depth, tilt_a=ta, tilt_b=tb)
    s.fill_holes_mm = fill
    try:
        await loop.run_in_executor(None, s.recut)
    except Exception as e:                                          # noqa: BLE001
        if outdated():
            return
        with state:
            state.hs_busy = False
            _hs_say(f"{type(e).__name__}: {e}", bad=True)
        return
    if outdated():
        return
    with state:
        state.hs_busy = False
        _hs_publish(s)
        _hs_draw()
        # The plane moved again while this cut ran: cut once more, to where it is now.
        if _hs_wanted() != wanted:
            hs_commit()


@ctrl.set("hs_commit")
def hs_commit(*_):
    """Re-cut at the plane and hole settings now in state. On slider release and on a
    committed number box, never per drag frame -- the disc alone follows the drag."""
    if HEART["s"] is None or state.hs_busy:
        return
    HEART["seq"] += 1
    gen, seq = _gen(), HEART["seq"]
    state.hs_busy = True
    _hs_say("Cutting…")
    asynchronous.create_task(_hs_recut(lambda: gen != _gen() or seq != HEART["seq"],
                                       _hs_wanted()))


@state.change("hs_depth", "hs_tilt_a", "hs_tilt_b")
def on_hs_plane(**_):
    """While dragging: move the disc only."""
    s = HEART["s"]
    if s is None or state.mode != "heart" or not state.hs_show_plane or not _throttled():
        return
    try:
        SCENE.set_heart_plane(s.disc(depth_mm=float(state.hs_depth),
                                     tilt_a=float(state.hs_tilt_a),
                                     tilt_b=float(state.hs_tilt_b)))
    except (TypeError, ValueError):
        return
    if ctrl.view_update:
        ctrl.view_update()


@state.change("hs_show_inner", "hs_show_outer", "hs_show_removed", "hs_show_plane")
def on_hs_view(**_):
    if state.mode == "heart":
        _hs_draw()


@state.change("hs_fill")
def on_hs_fill(**_):
    if HEART["s"] is not None:
        hs_commit()


@ctrl.set("hs_reset")
def hs_reset():
    """Back to the fitted rim plane, 3 mm toward the apex, untilted."""
    if HEART["s"] is None:
        return
    state.hs_depth, state.hs_tilt_a, state.hs_tilt_b = HEART_DEPTH_MM, 0.0, 0.0
    hs_commit()


@ctrl.set("hs_export")
def hs_export():
    s = HEART["s"]
    if PUBLIC or s is None:
        state.hs_export_bad = True
        state.hs_export_msg = "Nothing to export." if s is None else \
            "This instance does not write files for you."
        return
    try:
        # One folder per ventricle under the project's results/cardiac_walls/: the
        # layout run_cardiac_delta_k.py --reviewed-splits reads, so walls cut here
        # can be batch-measured as they are.
        stem = safe_stem(state.hs_export_name, "heart")
        paths = s.export(CARDIAC_WALLS / stem, stem)
        state.hs_export_bad = False
        state.hs_export_msg = (f"Wrote {paths['inner'].name}, {paths['outer'].name} and "
                               f"{paths['sidecar'].name} to results/cardiac_walls/{stem}/")
    except Exception as e:                                          # noqa: BLE001
        state.hs_export_bad = True
        state.hs_export_msg = str(e) if isinstance(e, HeartError) else \
            f"{type(e).__name__}: {e}"


@server.trigger("hs_measure")
def hs_measure(side):
    """Hand one wall to Shape, as Cut's "Measure δK on this cut" does. Refuses an
    invalid cut: δK on a wall that still carries part of the base measures the cut."""
    s = HEART["s"]
    if s is None or not s.ok:
        _hs_say("The cut is not valid, so measuring a wall would not mean anything."
                + (f" {s.problem}" if s is not None and s.problem else ""), bad=True)
        return
    wall = s.inner if side == "inner" else s.outer
    try:
        JOB.adopt_mesh(wall, f"{safe_stem(state.hs_export_name, 'heart')}_{side}",
                       expect_euler=int(wall.euler_number), kind="heart")
    except ShapeError as e:
        _hs_say(str(e), bad=True)
        return
    _sm_after_source()
    state.mode = "shape"


@state.change("mode")
def on_mode(mode, **_):
    """One render window, three views, so entering one takes the others down."""
    if mode == "shape":
        SCENE.clear_cut()
        SCENE.clear_heart()
        _draw_shape()
    elif mode == "heart":
        SCENE.clear_cut()
        SCENE.clear_patches()
        _hs_draw()
    else:
        SCENE.clear_patches()
        SCENE.clear_heart()
        _redraw(full=True)
        # Put back what clear_cut took down with the cut: the marks, and the measure
        # plane if the measurement is still on screen.
        SCENE.set_marks([mk["point"] for mk in (state.marks or []) if "point" in mk])
        if state.measure and state.m_show and SESSION["s"] is not None:
            _place_measure_plane(float(state.m_station))
    if ctrl.view_update:
        ctrl.view_update()


# ── layout ───────────────────────────────────────────────────────────────────
#
# The drawer holds the work, in the order it goes: Input, Cuts, Region, View & measure,
# each a white card. The case being worked on is pinned above them, and the verdict
# with Download is pinned below them, so neither scrolls out of reach. Everything that
# explains rather than controls lives behind an (i), so what stays on screen is the
# controls and their state. Colours, type and spacing are in webapp/static/app.css;
# the design and its reasons are recorded in the ui-23 build spec.
#
# Copy that reaches the template is written here, in Python, and placed as a plain
# text child. An apostrophe inside a bound expression takes the whole panel down
# without a word (see the marks note below), so the copy uses the typographic ’.
# Nothing here may carry "load" in a button's name ahead of the Load button, which the
# checks find by that word: "Upload" and "Download" both count.
CARD = "op-card"      # styled in webapp/static/app.css

TXT_TOOLS_ON = "Thoracic tools on: the Nub-in cut and Ishimaru zones."
TXT_TOOLS_OFF = ("Nub-in and Ishimaru zones are thoracic-only, so they are off for "
                 "this case.")
TXT_NUB_RULE = ("Cuts every branch to a stub 0.75 × the aorta’s local radius past the "
                "wall, and squares off the aorta’s own ends. A stub that will not cut "
                "cleanly is moved outward until it does; each move is listed under "
                "Details. Disabled while a region is isolated.")
TXT_NUB_REGION = ("Nub-in works on the whole vessel: it squares off the aorta’s own "
                  "ends, which this region has replaced.")
TXT_NUB_OFF = "Nub-in is a thoracic preset, and this case reads as abdominal."
TXT_ZONES_OFF = "Ishimaru zones: thoracic cases only."
TXT_APPLY_HINT = "Anatomy and smoothing are read when a case loads."
TXT_DL_REFUSE = "Download will refuse until the opening marked in red cuts."
TXT_AXIAL = "A horizontal slice, as the aorta’s ends are cut. The desktop’s Z-plane."
TXT_PERP = ("Square to the vessel, right for branch tips. The desktop’s "
            "centerline-normal.")
TXT_AXES = ("A and B are two directions across the cut plane, at right angles. "
            "Shrink Plane size if one cut makes two holes.")
TXT_EULER_TIP = ("Euler characteristic of the cut surface. A clean tube tree with N "
                 "openings has 2 − N; any other value means a hole or handle "
                 "somewhere.")
TXT_RESET_CONFIRM = ("Put every cut back to its default, drop trims, declared "
                     "arteries and the region, and reframe the view?")
TXT_FROM_TIP = ("How far back from the open end of the vessel this cut sits. The "
                "desktop calls it the trim distance.")
TXT_EMPTY_LOCAL = ("Choose a surface mesh and its centerline under Input, or open the "
                   "files already named there.")
TXT_SM_PUBLIC = "Here, δK is measured on the cut from the Cut section."
TXT_SM_DEFAULT = "Choose an open surface under Surface, then Compute δK."
TXT_PIN_HINT = ("Click the surface to pin a point and read how far it sits from where "
                "its vessel branches off.")
TXT_DIAMETER = ("Equivalent-area diameter, 2*sqrt(A/pi), sampled every 2 mm. A section "
                "is rarely a circle, so this is the diameter of the circle with the "
                "same area.")
TXT_SLIVERS = ("Tidies the cut edges. Changes the exported mesh. The desktop’s "
               "“Sliver-free rim”.")
TXT_ZONES = "Colours the arch by landing zone and lists where each zone runs."
TXT_PLASMA = ("Patches are coloured by integrated K, the quantity delta K is the "
              "area-weighted variance of, on CAMP’s own plasma map.")
TXT_SIDECAR = ("A row per patch, and a sidecar recording the face count, the trim and "
               "which CAMP build answered. A delta K without those cannot be "
               "reproduced or compared.")
# The "Change anatomy" link and the region chip both open things elsewhere in the
# drawer; they bring the drawer back first if it was put away (phone).
OPEN_REGION = ("cut_panels = cut_panels.includes('region') ? cut_panels "
               ": cut_panels.concat(['region']); main_drawer = true; "
               "utils.get('vcReveal')('[data-test=sec-region]')")
OPEN_DIAMETER = ("cut_panels = cut_panels.includes('view') ? cut_panels "
                 ": cut_panels.concat(['view']); main_drawer = true; "
                 "utils.get('vcReveal')('[data-test=vc-diameter]')")
DOWNLOAD_JS = ("trigger('download_cut').then(h => h && utils.get('vcDownload')(h))")
# The first failing opening's label, and the server's reason for it on its card.
# `openings` reaches the page as an object keyed "0", "1", ... (DeepReactive's copy),
# not an array: v-for iterates it, but .find/.some/.filter do not exist on it and
# .length is undefined. Everything outside the card loop reads it through
# Object.values().
BAD_LABEL = "(Object.values(openings || {}).find(o => o.bad) || {}).label"
WHY = "(utils.get('vcWhy') || (() => ''))(status_text, o.name)"
# The notice's lead-ins. The words are not changed: the split pieces rejoin to the
# paragraph exactly, and notice_version hashes the constants anyway.
NOTICE_LEADS = re.compile(r"(What is kept:|What they are used for:|How long they are "
                          r"kept:|Who can see them:|\[PI TO SUPPLY: [^\]]+\])")

# Page-only state: which panels are open, whether a click pins, whether the ID was
# noted, the colour key, the anatomy menu. Never sent to the server, so the gate never
# sees it and a wipe cannot leave it behind; each new page starts from these.
state.update({
    "cut_panels": ["input"],
    "out_panels": [],
    "shape_a": ["surface", "partition"],
    "shape_b": [],
    "heart_panels": ["hs_input", "hs_cut", "hs_holes"],
    # Off: a press is also how a rotate starts, so with Pin on every turn of the
    # model would drop a mark. Switched on when a point is to be measured.
    "pin_mode": False,
    "id_noted": "",
    "id_copied": False,
    "show_key": True,
    "fine_open": False,
    "anat_menu": False,
})
state.client_only("cut_panels", "out_panels", "shape_a", "shape_b", "heart_panels",
                  "pin_mode",
                  "id_noted", "id_copied", "main_drawer", "show_key", "fine_open",
                  "anat_menu")


def _info(tip, aria):
    """An (i) whose tooltip opens beside the drawer, over the canvas."""
    with v.VTooltip(text=tip, location="end", max_width="320"):
        with html.Template(v_slot_activator=("{ props }",)):
            v.VIcon("mdi-information-outline", v_bind=("props",), size="16",
                    classes="vc-info",
                    raw_attrs=['tabindex="0"', f'aria-label="{aria}"'])


def _tip_btn(tip, label=None, **kw):
    """A button with a tooltip. Only for ENABLED controls: a tooltip never shows on a
    disabled Vuetify button, so a disabled control's reason is visible text instead."""
    with v.VTooltip(text=tip, location="top"):
        with html.Template(v_slot_activator=("{ props }",)):
            if label is None:
                v.VBtn(v_bind=("props",), **kw)
            else:
                v.VBtn(label, v_bind=("props",), **kw)


def _fld(label):
    """A field with its label ABOVE it, never floating in its border: Vuetify's
    outlined label sits on the top edge, so stacked fields put each label through the
    box above it."""
    d = html.Div(classes="fld")
    with d:
        html.Label(label, classes="fld-lbl")
    return d


def _switch(model, label, tip=None, aria=None, test=None, meta=None, **kw):
    """A switch row: the label (its own, so clicking the words toggles it), an optional
    (i), optional state text at the right, and the switch."""
    attrs = [f'data-test="{test}"'] if test else []
    with v.VSwitch(v_model=(model,), inset=True, color="primary", density="compact",
                   hide_details=True, classes="vc-sw", raw_attrs=attrs, **kw):
        with html.Template(v_slot_label=True):
            html.Span(label, classes="vc-sw-t")
            if tip:
                _info(tip, aria or label)
            if meta:
                html.Span(f"{{{{ {meta} }}}}", v_if=meta, classes="vc-sw-m")


def _seg(model, items, test=None, update=None, full=True, disabled=None, dots=None):
    """A segmented control. `items` is [(label, value, v_if or None)]."""
    kw = dict(v_model=(model,)) if update is None else dict(
        model_value=(model,), __events=["update:modelValue"], update_modelValue=update)
    if disabled:
        kw["disabled"] = (disabled,)
    with v.VBtnToggle(mandatory=True, density="compact", variant="text",
                      classes="mode-tog vc-seg" + (" vc-seg-full" if full else ""),
                      raw_attrs=[f'data-test="{test}"'] if test else [], **kw):
        for label, value, cond in items:
            extra = dict(v_if=cond) if cond else {}
            with v.VBtn(value=value, size="small",
                        raw_attrs=[f'data-value="{value}"'], **extra):
                html.Span(label)
                if dots and value == dots:
                    html.I(classes="vc-ref")


def _ptitle(icon, title, done=None):
    """A panel title: an icon tile, then the exact title text. Pills follow."""
    if done:
        with html.Span(classes=(f"{done} ? 'vc-ptile done' : 'vc-ptile'",)):
            v.VIcon(icon=(f"{done} ? 'mdi-check' : '{icon}'",), size="17")
    else:
        with html.Span(classes="vc-ptile"):
            v.VIcon(icon, size="17")
    html.Span(title, classes="vc-ptitle")


def _notice(where):
    """The upload notice, whole, words unchanged, as a definition list: each lead-in
    in the left column, its text in the right. Hosted only (HOSTING.md job 7)."""
    with html.Div(classes=f"vc-notice vc-notice-{where}",
                  raw_attrs=['data-test="upload-notice"']):
        if notice_placeholders():
            with html.Div(classes="vc-draft"):
                v.VIcon("mdi-pencil-outline", size="15")
                html.Span("Draft: the bracketed parts have not been decided yet.")
        with html.Div(classes="vc-notice-h"):
            with html.Span(classes="vc-notice-ic"):
                v.VIcon("mdi-shield-check-outline", size="22")
            html.Span(NOTICE_HEADING, classes="vc-notice-title")
        last = len(NOTICE_PARAGRAPHS) - 1
        with html.Div(classes="vc-nt"):
            for k, para in enumerate(NOTICE_PARAGRAPHS):
                if k == last:
                    continue
                if para.startswith("Do not put names"):
                    with html.Div(classes="vc-callout-div warn"):
                        v.VIcon("mdi-shield-alert-outline", size="18")
                        html.Div(para)
                    continue
                # One row per lead-in; text before the first lead-in gets a row with an
                # empty key. The pieces are the paragraph's own, in order.
                rows, cur = [], None
                for part in NOTICE_LEADS.split(para):
                    if not part:
                        continue
                    if NOTICE_LEADS.fullmatch(part) and not part.startswith("[PI"):
                        cur = [part, []]
                        rows.append(cur)
                    else:
                        if cur is None:
                            cur = ["", []]
                            rows.append(cur)
                        cur[1].append(part)
                for key, parts in rows:
                    with html.Div(classes="vc-nt-row"):
                        html.Div(key, classes="vc-nt-k")
                        with html.Div(classes="vc-nt-v"):
                            for part in parts:
                                if part.startswith("[PI TO SUPPLY"):
                                    html.Span(part, classes="vc-ph")
                                else:
                                    html.Span(part)
        with html.Div(classes="vc-agree"):
            v.VIcon("mdi-check-circle-outline", size="18")
            html.Span(NOTICE_PARAGRAPHS[last])


def _pickers():
    # A surface AND a centerline: the tool cannot derive one from the other. The case
    # loads by itself once both are in.
    with html.Div(classes="vc-pickers"):
        _upload_field("up_mesh", "up_mesh_note", "up_mesh_bad",
                      "Surface mesh", "surface", "max_mesh_bytes")
        _upload_field("up_cl", "up_cl_note", "up_cl_bad",
                      "Centerline", "centerline", "max_cl_bytes")


def _anatomy_seg():
    with _fld("Anatomy"):
        _seg("family", [("Detect", "auto", None), ("Thoracic", "thoracic", None),
                        ("Abdominal", "abdominal", None)], test="family")


def _apply_block():
    # Anatomy and Smooth are read when a case loads, so changing either afterwards did
    # nothing until the next load -- and on the hosted site there was no next load.
    with html.Div(v_if="loaded && (family !== loaded_family || smooth !== loaded_smooth)",
                  classes="vc-callout-div warn vc-apply"):
        v.VIcon("mdi-alert", size="18")
        with html.Div():
            html.Div(TXT_APPLY_HINT)
            v.VBtn("Apply to this case", variant="outlined", size="x-small",
                   prepend_icon="mdi-restore", classes="vc-btn-secondary mt-2",
                   click=ctrl.reload_case, disabled=("busy",),
                   raw_attrs=['data-test="apply-anatomy"'])


def _readout(rows, cls="vc-rrow"):
    """Readout rows: label, value, unit, in that order -- the checks read them so."""
    for key, lbl, expr, unit, cond in rows:
        extra = dict(v_if=cond) if cond else {}
        with html.Div(classes=f"ro-row {cls}", raw_attrs=[f'data-key="{key}"'], **extra):
            html.Span(lbl, classes="ro-lbl")
            html.Span(f"{{{{ {expr} }}}}", classes="ro-val")
            html.Span(unit, classes="ro-u")


def _num(model, unit, change=None, mini=False, **kw):
    """A number box: the value, right-aligned, with its unit beside it."""
    ev = dict(__events=["change"], change=change) if change else {}
    with html.Div(classes="num mini" if mini else "num"):
        v.VTextField(type="number", density="compact", variant="plain",
                     hide_details=True, **model, **ev, **kw)
        html.Span(unit, classes="u")


with SinglePageWithDrawerLayout(server, vuetify_config=VUETIFY_CONFIG) as layout:
    layout.title.classes = "d-none"
    layout.icon.classes = "d-md-none vc-navicon"
    layout.toolbar.color = "chrome"
    layout.toolbar.height = 56
    layout.toolbar.flat = True
    layout.toolbar.classes = "vc-appbar"

    # ── app bar ──────────────────────────────────────────────────────────────
    with layout.toolbar:
        with html.Div(classes="vc-brand"):
            with html.Span(classes="vc-logo"):
                # An aortic arch and its three branches, every vessel end cut flat.
                html.Span(classes="vc-logo-glyph")
            with html.Span(classes="vc-wordmark d-none d-sm-inline-flex"):
                html.Span("Vascular")
                html.B("Cut")
        # Two sections, not two pages. They share a render window and a session, which
        # is what makes cut-then-measure one gesture instead of a download and an
        # upload -- and either one is usable without ever touching the other.
        with v.VBtnToggle(v_model=("mode",), mandatory=True, variant="text",
                          classes="top-tog vc-seg-top",
                          raw_attrs=['data-test="section-switch"']):
            v.VBtn("Cut", value="cut", raw_attrs=['data-test="section-cut"'])
            v.VBtn("Cardiac", value="heart", raw_attrs=['data-test="section-heart"'])
            with v.VBtn(value="shape", raw_attrs=['data-test="section-shape"']):
                html.Span("Shape δK")
        html.Span(classes="vc-vsep d-none d-md-inline-block",
                  v_if="(mode === 'cut' && loaded) || (mode === 'shape' && sm_result) "
                       "|| (mode === 'heart' && hs_loaded)")
        v.VChip("{{ hs_info ? hs_info.removed_pct.toFixed(0) : '–' }}% of outer removed",
                v_if="mode === 'heart' && hs_loaded && hs_info && $vuetify.display.mdAndUp",
                variant="flat", size="small", classes="vc-bchip",
                raw_attrs=['data-test="stat-heart-removed"'])
        v.VChip("{{ hs_info ? hs_info.holes_open + hs_info.holes_filled : 0 }} hole(s)",
                v_if="mode === 'heart' && hs_loaded && hs_info && $vuetify.display.mdAndUp",
                variant="flat", size="small",
                classes=("hs_info && hs_info.holes_open + hs_info.holes_filled ? "
                         "'vc-bchip bad' : 'vc-bchip ok'",),
                raw_attrs=['data-test="stat-heart-holes"'])
        v.VChip("{{ n_openings }} openings",
                v_if="mode === 'cut' && loaded && $vuetify.display.mdAndUp",
                variant="flat", size="small", classes="vc-bchip",
                raw_attrs=['data-test="stat-openings"', ':data-value="n_openings"'])
        with v.VTooltip(text=TXT_EULER_TIP, location="bottom",
                        v_if="mode === 'cut' && loaded && $vuetify.display.mdAndUp"):
            with html.Template(v_slot_activator=("{ props }",)):
                v.VChip("Euler {{ euler }}", v_bind=("props",), variant="flat",
                        size="small",
                        classes=("euler === 2 - n_openings ? 'vc-bchip ok' "
                                 ": 'vc-bchip bad'",),
                        prepend_icon=("euler === 2 - n_openings ? 'mdi-check' "
                                      ": 'mdi-alert-outline'",),
                        raw_attrs=['data-test="stat-euler"', ':data-value="euler"'])
        with v.VTooltip(location="bottom",
                        v_if="mode === 'cut' && loaded && $vuetify.display.mdAndUp"):
            with html.Template(v_slot_activator=("{ props }",)):
                v.VChip("{{ Math.round(area_mm2).toLocaleString() }} mm²",
                        v_bind=("props",), variant="flat", size="small",
                        classes="vc-bchip",
                        raw_attrs=['data-test="stat-area"', ':data-value="area_mm2"'])
            html.Span("Surface area of the cut mesh, "
                      "{{ n_faces.toLocaleString() }} triangles")
        # The band is in the SAME chip as the value, not a neighbouring one, so that
        # a screenshot of the header cannot separate them.
        v.VChip("δK {{ sm_result.k_fluct.toPrecision(4) }} "
                "± {{ sm_result.band_pct }}%",
                v_if="mode === 'shape' && sm_result", variant="flat", size="small",
                classes="vc-bchip accent")
        v.VChip("{{ sm_result.n_patches.toLocaleString() }} patches",
                v_if="mode === 'shape' && sm_result && $vuetify.display.mdAndUp",
                variant="flat", size="small", classes="vc-bchip")
        v.VSpacer()
        v.VProgressCircular(v_if="busy || sm_busy || hs_busy", indeterminate=True, size="20",
                            width="2", classes="mx-2 vc-spin")
        with v.VMenu(location="bottom end", close_on_content_click=False):
            with html.Template(v_slot_activator=("{ props }",)):
                v.VBtn(icon="mdi-information-outline", variant="text",
                       v_bind=("props",), classes="vc-bib",
                       raw_attrs=['aria-label="About"', 'data-test="about"'])
            with v.VCard(classes="vc-menu vc-about", width="340"):
                with html.Div(classes="vc-about-h"):
                    with html.Span(classes="vc-logo sm"):
                        html.Span(classes="vc-logo-glyph")
                    with html.Div():
                        html.Div("VascularCut", classes="vc-about-name")
                        html.Div(f"Build {BUILD} · CAMP: {{{{ sm_camp }}}}",
                                 classes="vc-about-meta")
                html.P("Cuts the same geometry as the desktop tool (branch_gui).")
                with html.Ol(classes="vc-about-steps"):
                    html.Li("Input: a surface mesh and its centerline.")
                    html.Li("Cuts: a preset, then single openings.")
                    html.Li("Region and View & measure, if needed.")
                    html.Li("Download.")
        # Delayed by CSS, so a slider drag does not flicker it. trame__busy, because
        # busy is set and cleared inside one synchronous trigger and never reaches
        # the page.
        v.VProgressLinear(indeterminate=True, absolute=True, location="bottom",
                          height="2", color="primary", classes="vc-busyline",
                          v_if=("trame__busy || busy || sm_busy || zones_busy || hs_busy",))

    # ── drawer ───────────────────────────────────────────────────────────────
    with layout.drawer as drawer:
        # Beside the view from 960 px up; below that an overlay over it, which
        # vc_boot.js starts closed so a phone opens on the model.
        drawer.permanent = None
        drawer.disable_resize_watcher = None
        drawer.mobile_breakpoint = "md"
        drawer.width = ("$vuetify.display.lgAndUp ? 456 : ($vuetify.display.mdAndUp "
                        "? 392 : 420)",)
        drawer.classes = "vc-drawer"

        # ── case header: pinned above the panels ─────────────────────────────
        with html.Template(v_slot_prepend=True):
            with html.Div(v_if="mode === 'cut'", classes="vc-case",
                          raw_attrs=['data-test="context"']):
                with html.Div(classes="vc-case-row"):
                    html.Span("{{ loaded ? (case_name || 'This case') : 'No case open' }}",
                              classes=("loaded ? 'vc-case-name' : 'vc-case-name none'",))
                    # The anatomy call, and the way to change it.
                    with v.VMenu(v_if="loaded", v_model=("anat_menu",),
                                 location="bottom start", close_on_content_click=False):
                        with html.Template(v_slot_activator=("{ props }",)):
                            with v.VBtn(v_bind=("props",), variant="flat",
                                        append_icon="mdi-chevron-down",
                                        classes=("'vc-hchip vc-anat ' + anat_level",),
                                        raw_attrs=['data-test="anat-chip"',
                                                   ':data-value="anat_kind"',
                                                   ':title="anat_text"']):
                                html.Span(classes="vc-anat-dot",
                                          v_if="anat_level === 'info'")
                                v.VIcon("mdi-alert-outline", size="15",
                                        v_if="anat_level !== 'info'")
                                html.Span("{{ anat_kind === 'thoracic' ? 'Thoracic' "
                                          ": 'Abdominal' }}")
                                html.Span("· {{ anat_how }}", classes="vc-anat-how")
                        with v.VCard(classes="vc-menu", width="340"):
                            html.Div("Anatomy for {{ case_name || 'this case' }}",
                                     classes="vc-menu-h")
                            html.Div("{{ anat_text }}", v_if="anat_text",
                                     classes="vc-menu-cap")
                            html.Div(TXT_TOOLS_ON, v_if="is_thoracic",
                                     classes="vc-menu-cap",
                                     raw_attrs=['data-test="thoracic-note"',
                                                ':data-thoracic="is_thoracic"'])
                            html.Div(TXT_TOOLS_OFF, v_if="!is_thoracic",
                                     classes="vc-menu-cap",
                                     raw_attrs=['data-test="thoracic-note"',
                                                ':data-thoracic="is_thoracic"'])
                            for val, title, sub in [
                                    ("auto", "Detect from the geometry",
                                     "'Currently reads as ' + anat_kind"),
                                    ("thoracic", "Thoracic",
                                     "'Nub-in cut and Ishimaru zones on'"),
                                    ("abdominal", "Abdominal",
                                     "'Nub-in cut and Ishimaru zones off'")]:
                                with html.Div(classes=(f"family === '{val}' ? 'vc-opt on' "
                                                       f": 'vc-opt'",),
                                              click=f"family = '{val}'",
                                              raw_attrs=['role="menuitemradio"',
                                                         'tabindex="0"',
                                                         f':aria-checked="family === \'{val}\'"']):
                                    with html.Div(classes="vc-opt-t"):
                                        html.Span(title)
                                        v.VIcon("mdi-check", size="16",
                                                v_if=f"family === '{val}'")
                                    html.Div(f"{{{{ {sub} }}}}", classes="vc-opt-s")
                            with html.Div(classes="vc-menu-sw"):
                                _switch("smooth", "Smooth the centerline if it is jagged")
                            _apply_block()
                    # The region, and the one-click way back to the whole vessel.
                    with v.VBtn(v_if="loaded", variant="flat", click=OPEN_REGION,
                                prepend_icon="mdi-arrow-expand-vertical",
                                classes=("in_region ? 'vc-hchip vc-regchip on' "
                                         ": 'vc-hchip vc-regchip'",),
                                raw_attrs=['data-test="region-chip"']):
                        html.Span("{{ in_region ? region_name : 'Whole vessel' }}",
                                  classes="vc-regchip-t")
                        v.VIcon("mdi-chevron-down", size="16", v_if="!in_region")
                    v.VBtn(icon="mdi-close", v_if="loaded && in_region", variant="text",
                           size="x-small", classes="vc-regclose",
                           click=ctrl.clear_region,
                           raw_attrs=['aria-label="Whole vessel"', 'title="Whole vessel"'])
                html.Div("{{ n_openings }} openings · Euler {{ euler }} · "
                         "{{ Math.round(area_mm2).toLocaleString() }} mm²",
                         v_if="loaded && $vuetify.display.smAndDown",
                         classes="vc-stat")
                # The submission ID: the only way to have uploaded files found and
                # removed, so it stays in view for as long as the case is open.
                with html.Div(v_if="submission_id", classes="vc-subid",
                              raw_attrs=['data-test="submission-id"']):
                    with html.Div(classes="vc-subid-h"):
                        html.Span("Submission ID")
                        v.VSpacer()
                        v.VBtn("{{ id_copied ? 'Copied' : 'Copy' }}", variant="text",
                               size="small", color="primary", classes="vc-btn-text",
                               prepend_icon=("id_copied ? 'mdi-check' "
                                             ": 'mdi-content-copy'",),
                               click=("utils.get('vcCopy')(submission_id); "
                                      "id_copied = true"),
                               raw_attrs=['aria-label="Copy the submission ID"'])
                    html.Div("{{ submission_id }}", classes="vc-id")
                    # Hosted, the page is wiped the moment it closes, so the only
                    # moment to note the ID down is now.
                    if PUBLIC:
                        with html.Div(v_if="id_noted !== submission_id",
                                      classes="vc-callout-div warn mt-2"):
                            v.VIcon("mdi-alert", size="18")
                            with html.Div():
                                html.Div("Note this ID now: closing or reloading this "
                                         "page clears it, and it is how your files are "
                                         "found.")
                                v.VBtn("Got it", variant="outlined", size="x-small",
                                       classes="vc-btn-secondary mt-2",
                                       click="id_noted = submission_id")
                # How far to trust the anatomy call, when that is in doubt. A routine
                # call says nothing here; its words are in the chip's menu.
                v.VAlert("{{ anat_text }}",
                         v_if="loaded && anat_text && anat_level !== 'info'",
                         density="compact", variant="tonal",
                         classes=("anat_level === 'bad' ? 'vc-callout bad mt-2' "
                                  ": 'vc-callout warn mt-2'",),
                         type=("anat_level === 'bad' ? 'error' : 'warning'",),
                         icon=("anat_level === 'bad' ? '$error' : 'mdi-alert'",),
                         raw_attrs=['data-test="anat-caption"',
                                    ':data-level="anat_level"'])
            with html.Div(v_if="mode === 'heart'", classes="vc-case",
                          raw_attrs=['data-test="context-heart"']):
                with html.Div(classes="vc-case-row"):
                    html.Span("{{ hs_loaded ? hs_name : 'No ventricle open' }}",
                              classes=("hs_loaded ? 'vc-case-name' : 'vc-case-name none'",))
                    with html.Span(v_if="hs_loaded && hs_info", classes="vc-hchip vc-srcchip"):
                        v.VIcon("mdi-heart-outline", size="15")
                        html.Span("{{ hs_info.separated ? 'Plane split' "
                                  ": 'Plane + seam' }}")
            with html.Div(v_if="mode === 'shape'", classes="vc-case",
                          raw_attrs=['data-test="context-shape"']):
                with html.Div(classes="vc-case-row"):
                    html.Span("{{ case_name || 'Shape δK' }}", classes="vc-case-name")
                    with html.Span(v_if="sm_info", classes="vc-hchip vc-srcchip"):
                        v.VIcon("mdi-cube-outline", size="15")
                        html.Span("{{ sm_info.name }}")

        # ── dock: pinned below the panels ────────────────────────────────────
        with html.Template(v_slot_append=True):
            with html.Div(v_if="mode === 'cut' && loaded", classes="vc-dock",
                          raw_attrs=['data-test="dock"']):
                # The verdict, beside the Output title, right above Download.
                v.VChip("{{ busy ? 'Working' : (cut_ok ? 'Cut valid' : 'Cut not valid') }}",
                        v_if="$vuetify.display.mdAndUp", variant="flat",
                        classes=("busy ? 'vc-verdict idle' : (cut_ok ? 'vc-verdict ok' "
                                 ": 'vc-verdict bad')",),
                        prepend_icon=("busy ? 'mdi-progress-clock' : (cut_ok ? "
                                      "'mdi-check-circle' : 'mdi-alert-circle')",),
                        raw_attrs=['data-test="verdict"',
                                   ':data-value="cut_ok ? \'ok\' : \'bad\'"'])
                with v.VExpansionPanels(v_model=("out_panels",), multiple=True,
                                        flat=True, classes="vc-out"):
                    with v.VExpansionPanel(value="out",
                                           raw_attrs=['data-test="sec-output"']):
                        with v.VExpansionPanelTitle():
                            v.VIcon("mdi-tune-variant", size="17", classes="mr-2")
                            html.Span("Output")
                        with v.VExpansionPanelText():
                            with _fld("Name"):
                                v.VTextField(v_model=("export_name",), suffix=".zip",
                                             density="compact", variant="outlined",
                                             hide_details=True,
                                             raw_attrs=['data-test="export-name"'])
                            if not PUBLIC:
                                v.VBtn("Export open mesh…", click=ctrl.export_cut,
                                       block=True, variant="outlined",
                                       classes="vc-btn-secondary mt-3",
                                       prepend_icon="mdi-content-save-outline",
                                       disabled=("!loaded",),
                                       raw_attrs=['data-test="export"'])
                            # The handoff. Nothing is written to disk on the way, so
                            # the surface measured is the one on screen rather than
                            # whichever export happened to be lying in output/.
                            with v.VBtn(click=ctrl.sm_handoff, block=True,
                                        variant="tonal", color="primary",
                                        classes="vc-btn-feature mt-2",
                                        append_icon="mdi-arrow-right",
                                        disabled=("!loaded",),
                                        raw_attrs=['data-test="handoff"']):
                                html.Span("Measure δK on this cut")
                            if not PUBLIC:
                                html.Div("Export writes the mesh and its settings to "
                                         "output/ on this machine.",
                                         classes="vc-cap mt-2")
                # r && ...: the trigger returns None when there is nothing to give,
                # and utils.download would otherwise save that as a file.
                v.VBtn("Download open mesh (.zip)", click=DOWNLOAD_JS, block=True,
                       size="large", variant="flat",
                       classes=("cut_ok ? 'vc-btn-primary vc-dl' : 'vc-btn-bad vc-dl'",),
                       prepend_icon=("cut_ok ? 'mdi-tray-arrow-down' "
                                     ": 'mdi-alert-circle-outline'",),
                       raw_attrs=['data-test="download"'])
                html.Div("Saves {{ export_name || 'cut' }}.zip: the open mesh (.stl) "
                         "and its cut settings (.json).", v_if="cut_ok",
                         classes="vc-dock-note")
                html.Div(TXT_DL_REFUSE, v_if="!cut_ok", classes="vc-dock-note bad")
                v.VAlert("{{ export_msg }}", v_if="export_msg",
                         density="compact", variant="tonal",
                         classes=("export_bad ? 'vc-callout bad mt-2' "
                                  ": 'vc-callout ok mt-2'",),
                         type=("export_bad ? 'error' : 'success'",),
                         raw_attrs=['data-test="export-msg"',
                                    ':data-bad="export_bad"'])
            with html.Div(v_if="mode === 'heart' && hs_loaded", classes="vc-dock",
                          raw_attrs=['data-test="dock-heart"']):
                # Compact on purpose: the panels above hold the work, so the dock is
                # only the verdict and the way into Shape. Export is under Output.
                with html.Div(classes="d-flex align-center flex-wrap ga-2"):
                    v.VChip("{{ hs_busy ? 'Working' : (hs_ok ? 'Cut valid' "
                            ": 'Cut not valid') }}",
                            variant="flat",
                            classes=("hs_busy ? 'vc-verdict idle' : (hs_ok ? "
                                     "'vc-verdict ok' : 'vc-verdict bad')",),
                            prepend_icon=("hs_busy ? 'mdi-progress-clock' : (hs_ok ? "
                                          "'mdi-check-circle' : 'mdi-alert-circle')",),
                            raw_attrs=['data-test="verdict-heart"'])
                html.Div("{{ hs_problem }}", v_if="!hs_ok && hs_problem",
                         classes="vc-dock-note bad")
                # The handoff, one wall at a time: Shape measures one surface.
                with html.Div(classes="d-flex ga-2 mt-2"):
                    for side, label in (("inner", "δK: inner"), ("outer", "δK: outer")):
                        with v.VBtn(click=f"trigger('hs_measure', ['{side}'])",
                                    variant="flat", classes="vc-btn-primary flex-grow-1",
                                    append_icon="mdi-arrow-right",
                                    disabled=("!hs_ok || hs_busy",),
                                    raw_attrs=[f'data-test="hs-measure-{side}"',
                                               f'title="Measure delta K on the {side} wall"']):
                            html.Span(label)
            with html.Div(v_if="mode === 'shape'", classes="vc-dock vc-dock-shape"):
                with v.VBtn(click=ctrl.sm_run, block=True, size="large",
                            variant="flat", classes="vc-btn-primary",
                            prepend_icon="mdi-chart-bell-curve",
                            disabled=("!sm_info || sm_busy",),
                            raw_attrs=['data-test="sm-compute"']):
                    html.Span("Compute δK")
                html.Div("Choose a surface first, under Surface.", v_if="!sm_info",
                         classes="vc-dock-note")

        # The Cardiac Ventricle panel folds once a ventricle is split, so Top cut sits
        # high, as Cut's Input does on load.
        client.ClientStateChange(
            value="hs_loaded",
            change="hs_loaded && (heart_panels = heart_panels.filter(k => k !== 'hs_input'))")
        # Input folds away and Cuts opens on every load, including a reload of the
        # case already on screen -- which is why this watches a counter and not
        # `loaded`.
        client.ClientStateChange(value="load_seq",
                                 change="load_seq && (cut_panels = ['cuts'], out_panels = [])")
        # The opening that will not cut, brought into view when the cut turns bad.
        client.ClientStateChange(
            value="cut_ok",
            change=(f"cut_ok || (cut_panels.includes('cuts') || (cut_panels = "
                    f"cut_panels.concat(['cuts'])), utils.get('vcRevealOp')({BAD_LABEL}))"))
        # δK arriving folds Surface and Partition, so the result sits high; their
        # summaries stay in the panels' titles.
        client.ClientStateChange(
            value="sm_result",
            change=("sm_result && (shape_a = shape_a.filter(k => k !== 'partition' "
                    "&& k !== 'surface'))"))

        # ═══ cut ═════════════════════════════════════════════════════════════
        with html.Div(v_if="mode === 'cut'", classes="vc-scroll"):
            with v.VExpansionPanels(v_model=("cut_panels",), multiple=True, flat=True,
                                    classes="vc-panels",
                                    raw_attrs=['data-test="cut-panels"']):

                # ── Input ────────────────────────────────────────────────────
                with v.VExpansionPanel(value="input", raw_attrs=['data-test="sec-input"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-folder-open-outline", "Input", done="loaded")
                        html.Span(classes="vc-dot bad", v_if="up_mesh_bad || up_cl_bad")
                    with v.VExpansionPanelText():
                        if PUBLIC:
                            # Before a case: what to do, in order. The notice and the
                            # pickers are beside this panel, over the empty view.
                            with html.Div(v_if="!loaded", classes="vc-steps"):
                                for n, (t, s) in enumerate([
                                        ("Read the data notice",
                                         "Beside this panel. Your files are kept."),
                                        ("Add the surface mesh",
                                         f"{' or '.join(STORE.accepted_extensions('surface')[:2])}, "
                                         f"up to {STORE.limit_mb('surface'):.0f} MB"),
                                        ("Add the centerline",
                                         "A point table or polydata file, up to "
                                         f"{STORE.limit_mb('centerline'):.0f} MB")], 1):
                                    with html.Div(classes="vc-step now" if n == 1
                                                  else "vc-step"):
                                        html.Span(str(n), classes="vc-step-n")
                                        with html.Div():
                                            html.Div(t, classes="vc-step-t")
                                            html.Div(s, classes="vc-step-s")
                                with html.Div(classes="vc-callout-div info mt-2"):
                                    v.VIcon("mdi-information-outline", size="18")
                                    html.Div("The case opens by itself once both have "
                                             "arrived. The anatomy is read from the "
                                             "geometry, and you can change it "
                                             "afterwards.")
                            # After a case: read it again with other options, or send
                            # another. The notice comes first, in full, every time.
                            with html.Div(v_if="loaded"):
                                _anatomy_seg()
                                _switch("smooth", "Smooth the centerline if it is jagged",
                                        tip="De-staircase: smooths a stepped centerline "
                                            "before cutting.")
                                html.Div("{{ quality_text }}", v_if="quality_text",
                                         classes="vc-cap")
                                _apply_block()
                                html.Div(classes="vc-divider")
                                _notice("drawer")
                                _pickers()
                                html.Div("Submission ID {{ submission_id }}",
                                         v_if="submission_id", classes="upnote sub-id")
                        else:
                            # THE TYPED PATH, on this machine only. It is how a
                            # researcher points at a lab share, and in public mode it is
                            # an arbitrary-file-read primitive, so it is not rendered
                            # at all there -- and load_model refuses it server-side
                            # anyway, since any client can fire that trigger with a
                            # path in state.
                            with _fld("Surface mesh path"):
                                v.VTextField(
                                    v_model=("mat_path",), density="compact",
                                    variant="outlined", hide_details=True,
                                    prepend_inner_icon="mdi-file-outline",
                                    placeholder=".mat / .stl / .ply",
                                    __events=[("keyup_enter", "keyup.enter")],
                                    keyup_enter=ctrl.load_model,
                                    raw_attrs=['data-test="path-mesh"'])
                                html.Div("{{ (utils.get('vcBase') || (x => x))(mat_path) }}",
                                         v_if="mat_path && $vuetify.display.smAndDown",
                                         classes="vc-cap vc-basename")
                            with _fld("Centerline path"):
                                v.VTextField(
                                    v_model=("csv_path",), density="compact",
                                    variant="outlined", hide_details=True,
                                    prepend_inner_icon="mdi-file-outline",
                                    placeholder=".csv / .vtp / .json",
                                    __events=[("keyup_enter", "keyup.enter")],
                                    keyup_enter=ctrl.load_model,
                                    raw_attrs=['data-test="path-centerline"'])
                                html.Div("{{ (utils.get('vcBase') || (x => x))(csv_path) }}",
                                         v_if="csv_path && $vuetify.display.smAndDown",
                                         classes="vc-cap vc-basename")
                            _anatomy_seg()
                            _switch("smooth", "Smooth the centerline if it is jagged",
                                    tip="De-staircase: smooths a stepped centerline "
                                        "before cutting.")
                            html.Div("{{ quality_text }}", v_if="quality_text",
                                     classes="vc-cap")
                            _apply_block()
                            v.VBtn("Load", click=ctrl.load_model, block=True,
                                   variant="flat", classes="vc-btn-primary vc-load mt-3",
                                   raw_attrs=['data-test="load"'])
                            with html.Details(classes="vc-disclose"):
                                with html.Summary():
                                    v.VIcon("mdi-tray-arrow-up", size="18")
                                    html.Span("Or upload files")
                                    v.VIcon("mdi-chevron-down", size="18",
                                            classes="vc-disclose-chev")
                                _pickers()
                                html.Div("Submission ID {{ submission_id }}",
                                         v_if="submission_id", classes="upnote sub-id")

                # ── Cuts ─────────────────────────────────────────────────────
                with v.VExpansionPanel(value="cuts", raw_attrs=['data-test="sec-cuts"'],
                                       classes=("loaded ? '' : 'vc-dim'",)):
                    with v.VExpansionPanelTitle(classes="vc-sticky"):
                        _ptitle("mdi-content-cut", "Cuts")
                        html.Span("{{ n_openings }} openings", v_if="loaded",
                                  classes="vc-pill")
                        with html.Span(v_if="loaded && !cut_ok && "
                                            "Object.values(openings || {}).some(o => o.bad)",
                                       classes="vc-pill bad"):
                            v.VIcon("mdi-alert-circle", size="14")
                            html.Span("{{ Object.values(openings || {}).filter(o => o.bad).length }} "
                                      "won’t cut")
                        html.Span("After a case opens", v_if="!loaded",
                                  classes="vc-pill")
                    with v.VExpansionPanelText():
                        # A failing opening, named, where it cannot be missed.
                        with html.Div(v_if="loaded && !cut_ok && Object.values(openings || {}).some(o => o.bad)",
                                      classes="vc-banner-wrap",
                                      raw_attrs=['data-test="banner"']):
                            with html.Div(classes="vc-banner"):
                                v.VIcon("mdi-alert-circle", size="20")
                                with html.Div(classes="vc-banner-tx"):
                                    html.Div("{{ (Object.values(openings || {}).find(o => o.bad) || {}).name }} "
                                             "won’t cut", classes="vc-banner-t")
                                    html.Div("{{ Object.values(openings || {}).filter(o => o.bad).length === 1 ? "
                                             "'The other ' + (Object.keys(openings || {}).length - 1) + "
                                             "' openings cut cleanly.' : "
                                             "Object.values(openings || {}).filter(o => o.bad).length + "
                                             "' openings will not cut.' }}",
                                             classes="vc-banner-s")
                                v.VBtn("Show", variant="outlined", size="x-small",
                                       classes="vc-btn-bad",
                                       prepend_icon="mdi-crosshairs-gps",
                                       click=(f"selected === {BAD_LABEL} || "
                                              f"trigger('select_row', [{BAD_LABEL}]); "
                                              f"utils.get('vcRevealOp')({BAD_LABEL})"))

                        # The nub-in preset: every opening at once, to the cohort's
                        # own rule, measured over 82 vessels in 30 THORACIC cases --
                        # so it exists here only for a thoracic case, and anything
                        # else says why it is absent rather than greying it out.
                        with html.Div(v_if="loaded && is_thoracic",
                                      classes=("in_region ? 'vc-preset off' "
                                               ": 'vc-preset'",)):
                            with html.Span(classes="vc-preset-ic"):
                                v.VIcon("mdi-content-cut", size="19")
                            with html.Div(classes="vc-preset-tx"):
                                with html.Div(classes="vc-preset-t"):
                                    html.Span("Nub-in cut")
                                    _info(TXT_NUB_RULE, "Nub-in cut, the thoracic preset")
                                html.Div("Thoracic preset · every opening at once",
                                         v_if="!in_region", classes="vc-preset-s")
                                with html.Div(v_if="in_region", classes="vc-preset-s"):
                                    html.Span(TXT_NUB_REGION)
                                    v.VBtn("Whole vessel", variant="text", size="x-small",
                                           color="primary", classes="vc-btn-text",
                                           click=ctrl.clear_region)
                            v.VBtn("Apply", variant="flat", size="small",
                                   classes="vc-btn-primary",
                                   disabled=("!can_nub || busy",),
                                   click="trigger('nub_preset')",
                                   raw_attrs=['data-test="nub-in"',
                                              'aria-label="Apply the Nub-in cut"'])
                        with html.Div(v_if="loaded && !is_thoracic",
                                      classes="vc-callout-div info mb-3",
                                      raw_attrs=['data-test="nub-off"']):
                            v.VIcon("mdi-information-outline", size="18")
                            with html.Div():
                                html.Div(TXT_NUB_OFF)
                                v.VBtn("Change anatomy", variant="text", size="x-small",
                                       color="primary", classes="vc-btn-text mt-1",
                                       prepend_icon="mdi-pencil-outline",
                                       click="anat_menu = true; main_drawer = true")

                        # Move every branch together. Setting eight of them one at a
                        # time to try a slightly shallower trim is the kind of thing
                        # that makes a tool feel like work.
                        # NB: the condition is a flag set in Python, not
                        # "openings.length > 1". A bare > inside a v-if attribute
                        # closes the tag when the template is rendered, and the whole
                        # element vanishes from the DOM without a word.
                        with html.Div(v_if="loaded && can_nudge", classes="vc-nudge"):
                            with html.Span(classes="vc-lbl"):
                                html.Span("All branch cuts")
                                _info("Moves every branch cut together.",
                                      "All branch cuts")
                            with html.Div(classes="vc-stepper"):
                                _tip_btn("Every branch cut moves deeper, away from its "
                                         "tip. The desktop’s ▼ deeper.", "Keep less",
                                         variant="text", prepend_icon="mdi-arrow-down",
                                         click="trigger('nudge_branches', [nudge_mm])",
                                         raw_attrs=['data-test="nudge-less"'])
                                _num(dict(v_model_number=("nudge_mm",)), "mm", step="1")
                                _tip_btn("Every branch cut moves toward its tip, keeping "
                                         "more of each branch. The desktop’s ▲ toward "
                                         "tips.", "Keep more", variant="text",
                                         append_icon="mdi-arrow-up",
                                         click="trigger('nudge_branches', [-nudge_mm])",
                                         raw_attrs=['data-test="nudge-more"'])

                        html.Pre("Available once a case is loaded." if PUBLIC
                                 else EMPTY_PRELOAD, v_if="!loaded", classes="empty")
                        html.Pre("{{ empty_text }}",
                                 v_if="loaded && Object.keys(openings || {}).length === 0",
                                 classes="empty")

                        # ── one card per opening ─────────────────────────────
                        #
                        # DeepReactive, or the thumb does not move until you let go.
                        #
                        # trame's state is shallow: assigning o.offset inside a v-for
                        # mutates an object in the array without changing the array
                        # reference, so Vue never re-renders. The control sits still
                        # while the 3-D view tracks the drag perfectly. Measured it --
                        # the thumb stayed at x=41 through an entire drag and jumped to
                        # x=181 on release. This wraps the list in a deeply reactive
                        # proxy so the assignment is seen, and every control then
                        # assigns locally AND fires a trigger.
                        with client.DeepReactive("openings"):
                          with html.Div(classes="vc-ops"):
                            # A plain div, not a VCard: the card brought Vuetify styling
                            # that fought the stylesheet.
                            with html.Div(
                                    v_for="(o, i) in openings", key="o.label",
                                    classes=("'op-card' + (selected === o.label ? ' is-sel' "
                                             ": '') + (o.bad ? ' is-bad' : '')",),
                                    raw_attrs=['data-test="opening"',
                                               ':data-label="o.label"',
                                               ':data-bad="o.bad"']):
                                # Row 1: which opening, how its plane lies, how wide.
                                # A click anywhere on it but the name and the plane
                                # selects the opening and opens its details.
                                with html.Div(classes="op-head",
                                              click="trigger('select_row', [o.label])"):
                                    html.Span(
                                        v_if="o.colour",
                                        classes=("colour_branches && !colour_zones "
                                                 "? 'op-sw' : 'op-sw ring'",),
                                        style=("'--c:' + o.colour",),
                                        title="Branch colour")
                                    v.VTextField(
                                        v_model=("o.name",), density="compact",
                                        variant="plain", hide_details=True,
                                        classes="op-name",
                                        __events=["change", ("click_stop", "click.stop")],
                                        change="trigger('rename_opening', "
                                               "[i, $event.target.value])",
                                        click_stop="0",
                                        title=("o.name + ' · click to rename'",),
                                        raw_attrs=['data-test="opening-name"',
                                                   'aria-label="Opening name"'])
                                    # Compact tags, the word in the tooltip, so the name
                                    # keeps its room.
                                    with html.Span(v_if="o.bound", classes="op-tag",
                                                   title="Set by the region: this end is "
                                                         "where the region starts or stops.",
                                                   raw_attrs=['aria-label="Region end"']):
                                        v.VIcon("mdi-arrow-collapse-vertical", size="13")
                                    with html.Span(v_if="o.is_trim", classes="op-tag",
                                                   title="A trimmed branch group: one cut "
                                                         "in place of several openings.",
                                                   raw_attrs=['aria-label="Trim"']):
                                        v.VIcon("mdi-source-branch", size="13")
                                    with v.VMenu(location="bottom end"):
                                        with html.Template(v_slot_activator=("{ props }",)):
                                            v.VBtn("{{ o.mode === 'normal' ? "
                                                   "'Perpendicular' : 'Axial' }}",
                                                   v_bind=("props",), variant="text",
                                                   size="x-small", classes="vc-plane",
                                                   append_icon="mdi-chevron-down",
                                                   disabled=("o.bound",),
                                                   __events=[("click_stop", "click.stop")],
                                                   click_stop="0",
                                                   raw_attrs=['data-test="opening-mode"'])
                                        with v.VCard(classes="vc-menu", width="300"):
                                            for val, title, sub in [
                                                    ("z", "Axial", TXT_AXIAL),
                                                    ("normal", "Perpendicular", TXT_PERP)]:
                                                with html.Div(
                                                        classes=(f"o.mode === '{val}' ? "
                                                                 f"'vc-opt on' : 'vc-opt'",),
                                                        click=(f"o.mode = '{val}'; "
                                                               f"trigger('set_mode', "
                                                               f"[i, '{val}'])"),
                                                        raw_attrs=['role="menuitemradio"',
                                                                   'tabindex="0"']):
                                                    with html.Div(classes="vc-opt-t"):
                                                        html.Span(title)
                                                        v.VIcon("mdi-check", size="16",
                                                                v_if=f"o.mode === '{val}'")
                                                    html.Div(sub, classes="vc-opt-s")
                                    with html.Span(classes="op-dia"):
                                        v.VIcon("mdi-diameter-variant", size="14")
                                        html.Span("{{ o.diameter ? o.diameter.toFixed(1) "
                                                  ": '-' }}")
                                        html.Span("mm", v_if="o.diameter", classes="u")
                                    # A trim stands for several openings, so it needs a
                                    # way back.
                                    v.VBtn(icon="mdi-undo-variant", v_if="o.is_trim",
                                           size="x-small", variant="text",
                                           classes="trim-undo",
                                           __events=[("click_stop", "click.stop")],
                                           click_stop=("trigger('remove_trim', "
                                                       "[o.trim_segment])"),
                                           raw_attrs=['aria-label="Undo this trim"',
                                                      'data-test="opening-undo-trim"'])
                                    v.VIcon(icon=("selected === o.label ? 'mdi-chevron-up' "
                                                  ": 'mdi-chevron-down'",),
                                            size="18", classes="op-exp",
                                            raw_attrs=['aria-label="Plane details"',
                                                       'data-test="opening-expand"',
                                                       'role="button"'])

                                # Row 2: the depth, as slider and number, and the
                                # pull back with the number worked out.
                                with html.Div(classes="op-row"):
                                    html.Span("From tip", classes="vc-lbl",
                                              title=TXT_FROM_TIP)
                                    v.VSlider(model_value=("o.offset",), min="1",
                                              max=("o.limit",), step="0.1",
                                              density="compact", hide_details=True,
                                              disabled=("o.bound",), color="primary",
                                              track_size="4", thumb_size="16",
                                              __events=["end", "update:modelValue"],
                                              update_modelValue=(
                                                  "o.offset = $event; "
                                                  "trigger('preview_offset', [i, $event])"),
                                              end="trigger('commit_offset', [i, $event])",
                                              raw_attrs=['data-test="opening-depth"'])
                                    _num(dict(v_model_number=("o.offset",)), "mm",
                                         change="trigger('commit_offset', "
                                                "[i, Number($event.target.value)])",
                                         disabled=("o.bound",), step="1",
                                         raw_attrs=['data-test="opening-depth-num"'])
                                    v.VBtn(icon="mdi-arrow-collapse-up", v_if="!o.bad",
                                           size="small", variant="text",
                                           disabled=("o.bound",), classes="pull-btn",
                                           click="trigger('pull_back', [i])",
                                           title="Pull back to a clean cross-section",
                                           raw_attrs=['aria-label="Pull back to a clean '
                                                      'cross-section"',
                                                      'data-test="opening-pull"'])
                                    v.VBtn("Pull back", v_if="o.bad", size="small",
                                           variant="flat", classes="vc-btn-danger pull-btn",
                                           prepend_icon="mdi-arrow-collapse-up",
                                           disabled=("o.bound",),
                                           click="trigger('pull_back', [i])",
                                           raw_attrs=['data-test="opening-pull"'])
                                with html.Div(v_if="selected === o.label",
                                              classes="op-range"):
                                    html.Span("1 mm")
                                    html.Span("{{ (+o.limit).toFixed(1) }} mm")
                                # Row 3, on a failing opening: the server's own reason.
                                with html.Div(v_if=f"o.bad && {WHY}", classes="op-why"):
                                    v.VIcon("mdi-alert-circle", size="16")
                                    html.Span(f"{{{{ {WHY} }}}}")

                                # ── details, revealed by selecting the opening ─────
                                with html.Div(v_if="selected === o.label", classes="op-adv"):
                                    # WHICH ARTERY this is, as opposed to what it is
                                    # called. The name is a label; this drives the region
                                    # presets, and on the half of the abdominal cohort
                                    # where no rule can tell which lateral pair the renals
                                    # are, it is the only thing that brings the aortic
                                    # neck back.
                                    with html.Div(classes="op-adv-row"):
                                        with html.Span(classes="vc-lbl"):
                                            html.Span("Artery")
                                            _info("Declaring the artery drives the region "
                                                  "presets.", "Artery")
                                        v.VSelect(model_value=("o.artery",),
                                                  menu_icon="mdi-chevron-down",
                                                  items=("artery_options",),
                                                  density="compact", variant="outlined",
                                                  hide_details=True,
                                                  classes="artery-pick vc-field",
                                                  clearable=True,
                                                  raw_attrs=['data-test="opening-artery"'],
                                                  __events=["update:modelValue"],
                                                  update_modelValue=(
                                                      "o.artery = $event || ''; "
                                                      "trigger('declare_artery', "
                                                      "[i, $event || ''])"))
                                    with html.Div(classes="vc-fine",
                                                  raw_attrs=['data-test="fine-tune"']):
                                        with html.Div(classes="vc-fine-h",
                                                      click="fine_open = !fine_open",
                                                      raw_attrs=['role="button"',
                                                                 'tabindex="0"',
                                                                 ':aria-expanded="fine_open"']):
                                            v.VIcon("mdi-tune-variant", size="17")
                                            with html.Div(classes="vc-fine-tx"):
                                                with html.Div(classes="vc-fine-t"):
                                                    html.Span("Fine-tune the plane")
                                                    _info(TXT_AXES, "Plane axes")
                                                html.Div(
                                                    "Tilt {{ (+o.tilt_a).toFixed(0) }}°, "
                                                    "{{ (+o.tilt_b).toFixed(0) }}° · Shift "
                                                    "{{ (+o.shift_a).toFixed(0) }}, "
                                                    "{{ (+o.shift_b).toFixed(0) }} mm · "
                                                    "{{ (+o.band).toFixed(0) }} mm"
                                                    "{{ o.band_auto ? ' auto' : '' }}",
                                                    v_if="!fine_open", classes="vc-fine-s")
                                            v.VBtn("Reset plane", v_if="fine_open",
                                                   variant="text", size="x-small",
                                                   color="primary",
                                                   classes="vc-btn-text",
                                                   __events=[("click_stop", "click.stop")],
                                                   click_stop="trigger('reset_plane', [i])")
                                            v.VIcon(icon=("fine_open ? 'mdi-chevron-up' "
                                                          ": 'mdi-chevron-down'",),
                                                    size="18", classes="vc-fine-chev")
                                        with html.Div(v_if="fine_open", classes="vc-fine-b"):
                                            with html.Div(classes="vc-grid2"):
                                                for lbl, key, lim, unit in [
                                                        ("Tilt A", "tilt_a", "o.lim_tilt", "°"),
                                                        ("Tilt B", "tilt_b", "o.lim_tilt", "°"),
                                                        ("Shift A", "shift_a", "o.lim_shift", "mm"),
                                                        ("Shift B", "shift_b", "o.lim_shift", "mm")]:
                                                    with html.Div(classes="vc-cell"):
                                                        with html.Div(classes="vc-cell-h"):
                                                            html.Span(lbl, classes="vc-lbl")
                                                            _num(dict(model_value=(
                                                                f"o.{key} == null ? '' : "
                                                                f"(+o.{key}).toFixed(1)",)),
                                                                unit, mini=True,
                                                                change=(f"trigger('commit_plane', "
                                                                        f"[i, '{key}', Number("
                                                                        f"$event.target.value)])"))
                                                        v.VSlider(
                                                            model_value=(f"o.{key}",),
                                                            min=(f"-{lim}",), max=(lim,),
                                                            step="0.1", density="compact",
                                                            hide_details=True, track_size="4",
                                                            thumb_size="16", color="primary",
                                                            classes="vc-bip",
                                                            __events=["end",
                                                                      "update:modelValue"],
                                                            update_modelValue=(
                                                                f"o.{key} = $event; "
                                                                f"trigger('preview_plane', "
                                                                f"[i, '{key}', $event])"),
                                                            end=(f"trigger('commit_plane', "
                                                                 f"[i, '{key}', $event])"))
                                            with html.Div(classes="op-adv-row size"):
                                                html.Span("Plane size", classes="vc-lbl")
                                                v.VSlider(
                                                    model_value=("o.band",),
                                                    min=("o.lim_band_min",),
                                                    max=("o.lim_band_max",), step="0.1",
                                                    density="compact", hide_details=True,
                                                    track_size="4", thumb_size="16",
                                                    color="primary",
                                                    __events=["end", "update:modelValue"],
                                                    update_modelValue=(
                                                        "o.band = $event; "
                                                        "trigger('preview_plane', "
                                                        "[i, 'band', $event])"),
                                                    end="trigger('commit_plane', "
                                                        "[i, 'band', $event])")
                                                with html.Div(classes="num"):
                                                    v.VTextField(
                                                        model_value=("o.band == null ? '' "
                                                                     ": (+o.band).toFixed(1)",),
                                                        type="number", density="compact",
                                                        variant="plain", hide_details=True,
                                                        __events=["change"],
                                                        change=("trigger('commit_plane', "
                                                                "[i, 'band', Number("
                                                                "$event.target.value)])"))
                                                    html.Span("{{ o.band_auto ? 'auto' "
                                                              ": 'mm' }}", classes="u")

                        # Less often needed, so under the list.
                        with html.Div(v_if="loaded", classes="vc-tools"):
                            # Trim a branch group: one cut on the vessel FEEDING several
                            # openings, placed before it splits. Which group is
                            # collapsed is the whole choice, so it is a menu of them.
                            with v.VMenu(v_if="can_trim", location="end top",
                                         max_height="420"):
                                with html.Template(v_slot_activator=("{ props }",)):
                                    with html.Div(v_bind=("props",), classes="vc-tool",
                                                  raw_attrs=['data-test="trim-pick"',
                                                             'role="button"',
                                                             'tabindex="0"']):
                                        with html.Span(classes="vc-tool-ic"):
                                            v.VIcon("mdi-source-branch", size="16")
                                        html.Span("Trim a branch group…",
                                                  classes="vc-tool-t")
                                        html.Span("{{ trim_options.length }} "
                                                  "{{ trim_options.length === 1 ? "
                                                  "'group' : 'groups' }}",
                                                  classes="vc-tool-m")
                                        v.VIcon("mdi-chevron-right", size="18",
                                                classes="vc-tool-chev")
                                with v.VCard(classes="vc-menu", width="340"):
                                    html.Div("One cut before the vessel splits, in place "
                                             "of its openings.", classes="vc-menu-cap")
                                    with html.Div(v_for="t in trim_options",
                                                  key="t.value", classes="vc-opt",
                                                  click="trigger('add_trim', [t.value])",
                                                  raw_attrs=['role="menuitem"',
                                                             'tabindex="0"']):
                                        html.Div("{{ String(t.title).split(' · ')[0] }}",
                                                 classes="vc-opt-t")
                                        html.Div("{{ String(t.title).split(' · ')"
                                                 ".slice(1).join(' · ') }}",
                                                 classes="vc-opt-s")
                            with html.Div(classes="vc-tool"):
                                with html.Span(classes="vc-tool-ic"):
                                    v.VIcon("mdi-blur", size="16")
                                with html.Div(classes="vc-tool-sw"):
                                    _switch("clean_rims", "Remove rim slivers",
                                            tip=TXT_SLIVERS, test="clean-rims",
                                            meta="!clean_rims && (utils.get('vcSlivers') "
                                                 "|| (() => ''))(status_text)")

                # ── Region ───────────────────────────────────────────────────
                with v.VExpansionPanel(value="region", raw_attrs=['data-test="sec-region"'],
                                       classes=("loaded ? '' : 'vc-dim'",)):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-arrow-expand-vertical", "Region")
                        html.Span("{{ region_name }}", v_if="loaded && in_region",
                                  classes="vc-pill accent")
                        html.Span("After a case opens", v_if="!loaded",
                                  classes="vc-pill")
                    with v.VExpansionPanelText():
                        html.Div("Available once a case is loaded.", v_if="!loaded",
                                 classes="empty")
                        with html.Div(v_if="loaded"):
                            with _fld("Isolate a region"):
                                v.VSelect(v_model=("region_pick",),
                                          menu_icon="mdi-chevron-down",
                                          items=("region_items",), item_props=True,
                                          density="compact", variant="outlined",
                                          hide_details=True, disabled=("!loaded",),
                                          classes="vc-field vc-region-pick",
                                          raw_attrs=['data-test="region-pick"'])
                            html.Div("{{ region_why }}", v_if="in_region && region_why",
                                     classes="vc-cap")
                            with html.Div(v_if="in_region", classes="vc-bounds"):
                                for lbl, key, dis in [("Start", "region_s0", None),
                                                      ("End", "region_s1", "region_open")]:
                                    extra = dict(disabled=(dis,)) if dis else {}
                                    with html.Div(classes="op-row"):
                                        html.Span(lbl, classes="vc-lbl")
                                        v.VSlider(v_model=(key,), min="0",
                                                  max=("route_length",), step="0.1",
                                                  density="compact", hide_details=True,
                                                  track_size="4", thumb_size="16",
                                                  color="primary", __events=["end"],
                                                  end=ctrl.commit_region, **extra)
                                        _num(dict(v_model_number=(key,)), "mm",
                                             change=ctrl.commit_region, **extra)
                                html.Div("Distances are measured along the centerline.",
                                         classes="vc-cap")
                                # Disabled when the region has no end, as the desktop
                                # does (branch_gui._sync_region_bar).
                                html.Div("This region has no end: everything below the "
                                         "start is kept.", v_if="region_open",
                                         classes="vc-cap")

                # ── View & measure ───────────────────────────────────────────
                with v.VExpansionPanel(value="view", raw_attrs=['data-test="sec-view"'],
                                       classes=("loaded ? '' : 'vc-dim'",)):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-ruler-square", "View & measure")
                        html.Span("Max {{ measure.max_mm.toFixed(1) }} mm",
                                  v_if="measure", classes="vc-pill")
                        html.Span("After a case opens", v_if="!loaded",
                                  classes="vc-pill")
                    with v.VExpansionPanelText():
                        # ── measurement, first: it is what this panel is for ────
                        with html.Div(classes="vc-sub",
                                      raw_attrs=['data-test="vc-diameter"']):
                            html.Span("Diameter along the vessel")
                            _info(TXT_DIAMETER, "Diameter along the vessel")
                        v.VBtn("Measure", click="trigger('run_measure')",
                               block=True, variant="tonal", color="primary",
                               classes="vc-btn-feature", prepend_icon="mdi-ruler",
                               disabled=("!loaded || measuring",),
                               loading=("measuring",),
                               raw_attrs=['data-test="measure"'])
                        with html.Div(v_if="measure", classes="readout",
                                      raw_attrs=['data-test="readout-measure"']):
                            with html.Div(classes="vc-hero2"):
                                with html.Div(classes="ro-row vc-htile",
                                              raw_attrs=['data-key="max"']):
                                    html.Span("Max", classes="ro-lbl")
                                    html.Span("{{ measure.max_mm.toFixed(1) }}",
                                              classes="ro-val")
                                    html.Span("mm", classes="ro-u")
                                # The angle, and how well a circle actually fitted.
                                # Quoting one without the other invites a number from a
                                # fit that never described the vessel.
                                with html.Div(v_if="angle", classes="ro-row vc-htile",
                                              raw_attrs=['data-key="angle"',
                                                         'data-test="readout-angle"']):
                                    with html.Span(classes="ro-lbl"):
                                        html.Span("Angle")
                                        v.VIcon("mdi-alert", size="15",
                                                classes="vc-warn-ic",
                                                v_if="!angle.straight && angle.poor_fit")
                                    html.Span("{{ angle.straight ? '-' : "
                                              "angle.angle_deg.toFixed(0) }}",
                                              classes="ro-val")
                                    html.Span("°", classes="ro-u deg")
                            with html.Div(classes="vc-rgrid"):
                                _readout([
                                    ("at", "Widest at", "measure.max_at_mm.toFixed(0)",
                                     "mm along", None),
                                    ("min", "Min", "measure.min_mm.toFixed(1)", "mm", None),
                                    ("mean", "Mean", "measure.mean_mm.toFixed(1)", "mm",
                                     None),
                                    ("length", "Length", "measure.length_mm.toFixed(1)",
                                     "mm", None),
                                    ("start", "Start", "measure.start_mm.toFixed(1)", "mm",
                                     None),
                                    ("end", "End", "measure.end_mm.toFixed(1)", "mm", None),
                                    ("radius", "Radius", "angle.radius_mm.toFixed(0)",
                                     "mm", "angle && !angle.straight"),
                                    ("fit", "Fit", "angle.rms_mm.toFixed(2)", "mm rms",
                                     "angle && !angle.straight"),
                                    ("stations", "Stations", "measure.n_stations", "",
                                     None),
                                ])
                            html.Div("Region: {{ measure.region }}",
                                     v_if="measure.region", classes="vc-cap")
                            html.Div("Straight: no angulation to report.",
                                     v_if="angle && angle.straight", classes="vc-cap")
                            with html.Div(v_if="angle && !angle.straight && angle.poor_fit",
                                          classes="vc-callout-div warn mt-2"):
                                v.VIcon("mdi-alert", size="18")
                                html.Div("Poor fit, treat this angle with caution. A "
                                         "whole aorta is not one arc, isolate the zone "
                                         "first.")
                            # The station those numbers came from, drawn and movable. A
                            # summary can say "widest at 86 mm" and still be measuring
                            # the wrong place; this is the only way to see that.
                            with html.Div(classes="op-row vc-mplane"):
                                html.Span("Plane at", classes="vc-lbl")
                                v.VSlider(v_model=("m_station",),
                                          min=("m_lo",), max=("m_hi",), step="0.5",
                                          density="compact", hide_details=True,
                                          track_size="4", thumb_size="16",
                                          color="primary",
                                          __events=["end", "update:modelValue"],
                                          update_modelValue=("m_station = $event",),
                                          end="trigger('measure_station', [$event])")
                                _num(dict(v_model_number=("m_station",)), "mm",
                                     change="trigger('measure_station', "
                                            "[Number($event.target.value)])", step="1")
                            html.Div("{{ m_note }}", v_if="m_note", classes="vc-cap")
                            with html.Div(classes="vc-row-between"):
                                v.VBtn("Go to widest", size="small", variant="text",
                                       color="primary", classes="vc-btn-text",
                                       prepend_icon="mdi-target",
                                       click="trigger('measure_station', "
                                             "[measure.max_at_mm])")
                                v.VSwitch(v_model=("m_show",), label="Show plane",
                                          inset=True, color="primary",
                                          density="compact", hide_details=True,
                                          classes="vc-sw vc-sw-inline",
                                          __events=["update:modelValue"],
                                          update_modelValue=("m_show = $event; "
                                                             "trigger('measure_show', "
                                                             "[$event])"))

                        # ── pin a point ─────────────────────────────────────────
                        html.Div(classes="vc-divider")
                        with html.Div(classes="vc-sub"):
                            html.Span("Pin a point")
                            _info(TXT_PIN_HINT, "Pin a point")
                        _switch("pin_mode", "Pin a point on click", test="pin-switch")
                        # The text is composed in Python, not in the template. A Vue
                        # expression with nested quotes -- and an apostrophe in a
                        # neighbouring title -- took the WHOLE panel body down: every
                        # control in it vanished, with nothing in the console to say
                        # why.
                        html.Div("{{ marks_note || 'No marks yet.' }}", classes="vc-cap")
                        with html.Div(classes="vc-btnrow"):
                            _tip_btn("Move that cut to the pinned point.",
                                     "Send cut to last mark", size="small",
                                     variant="outlined", classes="vc-btn-secondary",
                                     prepend_icon="mdi-map-marker-outline",
                                     disabled=("!marks_note",),
                                     click="trigger('send_cut_to_mark')",
                                     raw_attrs=['data-test="send-to-mark"'])
                            v.VBtn("Clear marks", size="small", variant="text",
                                   classes="vc-btn-quiet", disabled=("!marks_note",),
                                   click="trigger('clear_marks')",
                                   raw_attrs=['data-test="clear-marks"'])

                        # ── what the view shows ─────────────────────────────────
                        html.Div(classes="vc-divider")
                        html.Div("Show", classes="vc-sub")
                        with html.Div(classes="vc-tchips"):
                            for key, label, glyph in [
                                    ("show_centerline", "Centerline", "cl"),
                                    ("show_planes", "Cut discs", "disc"),
                                    ("show_loops", "Opening rims", "rim")]:
                                with v.VBtn(variant="flat", size="small",
                                            classes=(f"{key} ? 'vc-tchip on' "
                                                     f": 'vc-tchip'",),
                                            click=f"{key} = !{key}",
                                            raw_attrs=[f':aria-pressed="{key}"']):
                                    html.Span(classes=f"vc-glyph {glyph}")
                                    html.Span(label)
                        with _fld("Colour the surface by"):
                            _seg("colour_zones ? 'zone' : (colour_branches ? 'branch' "
                                 ": 'none')",
                                 [("Branch", "branch", None),
                                  ("Zone", "zone", "is_thoracic"),
                                  ("None", "none", None)],
                                 test="colour-by",
                                 update=("colour_zones = $event === 'zone'; "
                                         "$event === 'zone' || "
                                         "(colour_branches = $event === 'branch')"))
                        with html.Div(classes="op-row vc-opacity"):
                            html.Span("Surface opacity", classes="vc-lbl")
                            v.VSlider(v_model=("opacity",), min="5", max="100",
                                      step="5", density="compact", hide_details=True,
                                      track_size="4", thumb_size="16", color="primary")
                            _num(dict(v_model_number=("opacity",)), "%")
                        # ISHIMARU ZONES, thoracic only, right under the colouring that
                        # uses them. A button rather than a checkbox because finding
                        # the ostia takes seconds, not milliseconds.
                        with html.Div(v_if="loaded && is_thoracic",
                                      raw_attrs=['data-test="zones"']):
                            with html.Div(classes="vc-sub vc-sub-row"):
                                html.Span("Ishimaru zones")
                                _info(TXT_ZONES, "Ishimaru zones")
                                v.VSpacer()
                                v.VBtn("Find zones", click=ctrl.run_zones,
                                       loading=("zones_busy",),
                                       disabled=("zones_busy",),
                                       variant="tonal", color="primary", size="small",
                                       classes="vc-btn-feature",
                                       prepend_icon="mdi-map-legend",
                                       raw_attrs=['data-test="find-zones"'])
                            html.Div("{{ zones_note }}", v_show=("zones_note",),
                                     classes="vc-cap")
                            with html.Div(v_for="z in zone_legend", key="z.zone",
                                          classes="zone-row"):
                                html.Span(style=("`background:${z.colour}`",),
                                          classes="zone-sw")
                                html.Span("Zone {{ z.zone }}", classes="zone-n")
                                html.Span("{{ z.name }}", classes="zone-name")
                                html.Span("{{ z.span }}", classes="zone-span")
                        html.Div(TXT_ZONES_OFF, v_if="loaded && !is_thoracic",
                                 classes="vc-cap vc-off",
                                 raw_attrs=['data-test="zones-off"'])

        # ═══ cardiac ═════════════════════════════════════════════════════════
        with html.Div(v_if="mode === 'heart'", classes="vc-scroll"):
            with v.VExpansionPanels(v_model=("heart_panels",), multiple=True, flat=True,
                                    classes="vc-panels",
                                    raw_attrs=['data-test="heart-panels"']):
                # ── Ventricle ────────────────────────────────────────────────
                with v.VExpansionPanel(value="hs_input",
                                       raw_attrs=['data-test="sec-hs-input"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-heart-outline", "Ventricle", done="hs_loaded")
                    with v.VExpansionPanelText():
                        if PUBLIC:
                            html.Div("The Cardiac section reads meshes from this "
                                     "machine, so it is available only in the local "
                                     "version.", classes="vc-cap")
                        else:
                            html.Div("A closed left-ventricle myocardium. It is split "
                                     "into inner and outer walls, then cut at the base.",
                                     classes="vc-cap")
                            with _fld("Mesh"):
                                v.VCombobox(v_model=("hs_path",), items=("hs_cases",),
                                            item_title="title", item_value="value",
                                            return_object=False, density="compact",
                                            variant="outlined", hide_details=True,
                                            prepend_inner_icon="mdi-file-outline",
                                            placeholder="dataset/…stl, or a path",
                                            raw_attrs=['data-test="hs-path"'])
                            v.VBtn("{{ hs_busy ? 'Splitting…' : 'Split this ventricle' }}",
                                   block=True, click=ctrl.hs_load, variant="flat",
                                   classes="vc-btn-primary mt-2",
                                   prepend_icon="mdi-call-split",
                                   disabled=("hs_busy",),
                                   raw_attrs=['data-test="hs-load"'])
                            html.Div("Takes a few seconds: the inner and outer seeds and "
                                     "the seam between them are found once, here.",
                                     classes="vc-cap mt-1")

                # ── Top cut ──────────────────────────────────────────────────
                with v.VExpansionPanel(value="hs_cut", disabled=("!hs_loaded",),
                                       raw_attrs=['data-test="sec-hs-cut"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-content-cut", "Top cut")
                        html.Span("{{ (+hs_depth).toFixed(1) }} mm below rim",
                                  v_if="hs_loaded", classes="vc-pill")
                    with v.VExpansionPanelText():
                        _info_txt = ("The plane starts fitted to the basal rim, where "
                                     "the outer wall folds over into the cavity, so a "
                                     "tilted base is followed. Depth moves it toward "
                                     "the apex along its own normal; the tilts turn it "
                                     "about two axes across it. Everything on the base "
                                     "side is removed from both walls.")
                        with html.Div(classes="vc-cap vc-lbl"):
                            html.Span("The disc follows the drag; the walls re-cut "
                                      "when you let go.")
                            _info(_info_txt, "How the top cut works")
                        for lbl, key, lo, hi, step, unit in [
                                ("Depth", "hs_depth", "-20", "60", "0.5", "mm"),
                                ("Tilt A", "hs_tilt_a", "-30", "30", "1", "°"),
                                ("Tilt B", "hs_tilt_b", "-30", "30", "1", "°")]:
                            with html.Div(classes="op-row"):
                                html.Span(lbl, classes="vc-lbl")
                                v.VSlider(v_model=(key,), min=lo, max=hi, step=step,
                                          density="compact", hide_details=True,
                                          track_size="4", thumb_size="16",
                                          color="primary", disabled=("hs_busy",),
                                          __events=["end"], end=ctrl.hs_commit,
                                          raw_attrs=[f'data-test="{key}"'])
                                _num(dict(v_model_number=(key,)), unit,
                                     change=ctrl.hs_commit, mini=True, step=step)
                        html.Div("Positive depth cuts deeper, toward the apex. "
                                 "Negative moves the plane up past the rim.",
                                 classes="vc-cap")
                        v.VBtn("Back to the fitted rim plane", block=True,
                               variant="outlined", classes="vc-btn-secondary mt-2",
                               prepend_icon="mdi-restore", click=ctrl.hs_reset,
                               disabled=("!hs_loaded || hs_busy",),
                               raw_attrs=['data-test="hs-reset"'])

                # ── Holes ────────────────────────────────────────────────────
                with v.VExpansionPanel(value="hs_holes", disabled=("!hs_loaded",),
                                       raw_attrs=['data-test="sec-hs-holes"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-circle-off-outline", "Holes")
                        html.Span("{{ hs_info.holes_open }} open · "
                                  "{{ hs_info.holes_filled }} patched",
                                  v_if="hs_info", classes="vc-pill")
                    with v.VExpansionPanelText():
                        html.Div("Where the segmentation punched through the wall, the "
                                 "inner and outer walls are joined by a small tube. "
                                 "The automatic seam cuts each one, leaving a hole in "
                                 "each wall, drawn as a red ring.", classes="vc-cap")
                        html.Div("Left open by default. CAMP does not need a closed "
                                 "surface, and a patch adds curvature the wall does "
                                 "not have — the reference method never fills them.",
                                 classes="vc-cap")
                        _switch("hs_fill", "Patch holes instead",
                                tip=("Closes each hole up to the size below with a fan "
                                     "of triangles. The sidecar records it. Use it to "
                                     "compare, not as the default."),
                                test="hs-fill", disabled=("hs_busy",))
                        with html.Div(v_if="hs_fill", classes="op-row"):
                            html.Span("Up to", classes="vc-lbl")
                            _num(dict(v_model_number=("hs_fill_mm",)), "mm around",
                                 change=ctrl.hs_commit, step="5")
                        html.Div("Hole sizes (mm around): {{ hs_info.hole_perimeters_mm"
                                 ".join(', ') }}",
                                 v_if="hs_info && hs_info.hole_perimeters_mm.length",
                                 classes="vc-cap mt-2")

                # ── View ─────────────────────────────────────────────────────
                with v.VExpansionPanel(value="hs_view", disabled=("!hs_loaded",),
                                       raw_attrs=['data-test="sec-hs-view"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-eye-outline", "View")
                    with v.VExpansionPanelText():
                        _switch("hs_show_inner", "Inner wall (blue)")
                        _switch("hs_show_outer", "Outer wall (tan)")
                        _switch("hs_show_removed", "What the cut removes (grey)")
                        _switch("hs_show_plane", "Cut plane")

                # ── Output ───────────────────────────────────────────────────
                with v.VExpansionPanel(value="hs_out", disabled=("!hs_loaded",),
                                       raw_attrs=['data-test="sec-hs-output"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-export-variant", "Output")
                    with v.VExpansionPanelText():
                        if PUBLIC:
                            html.Div("The walls are written only by the local version "
                                     "of this tool.", classes="vc-cap")
                        else:
                            with _fld("Name"):
                                v.VTextField(v_model=("hs_export_name",),
                                             density="compact", variant="outlined",
                                             hide_details=True,
                                             raw_attrs=['data-test="hs-export-name"'])
                            v.VBtn("Export walls…", click=ctrl.hs_export, block=True,
                                   variant="outlined", classes="vc-btn-secondary mt-2",
                                   prepend_icon="mdi-content-save-outline",
                                   disabled=("!hs_ok || hs_busy",),
                                   raw_attrs=['data-test="hs-export"'])
                            html.Div("Writes the inner and outer walls (.stl) and the "
                                     "cut settings (.json) to a folder of this name "
                                     "in results/cardiac_walls/, ready for the batch "
                                     "δK run.",
                                     classes="vc-cap mt-1")
                        v.VAlert("{{ hs_export_msg }}", v_if="hs_export_msg",
                                 density="compact", variant="tonal",
                                 classes=("hs_export_bad ? 'vc-callout bad mt-2' "
                                          ": 'vc-callout ok mt-2'",),
                                 type=("hs_export_bad ? 'error' : 'success'",))

            # ── the readout: not a panel, so it is never folded away ─────────
            with html.Div(v_if="hs_info", classes="readout vc-rgrid mx-3 mt-2",
                          raw_attrs=['data-test="readout-heart"']):
                _readout([
                    ("split", "Split", "hs_info.separated ? 'Plane alone' "
                     ": 'Plane + automatic seam'", "", None),
                    ("removed", "Outer removed", "hs_info.removed_pct.toFixed(1)", "%",
                     None),
                    ("inner_area", "Inner area",
                     "Math.round(hs_info.inner_area_mm2).toLocaleString()", "mm²", None),
                    ("outer_area", "Outer area",
                     "Math.round(hs_info.outer_area_mm2).toLocaleString()", "mm²", None),
                    ("inner_open", "Inner openings", "hs_info.inner_loops", "", None),
                    ("outer_open", "Outer openings", "hs_info.outer_loops", "", None),
                    ("inner_handles", "Inner bridges", "hs_info.inner_handles", "", None),
                    ("outer_handles", "Outer bridges", "hs_info.outer_handles", "", None),
                ])
            with html.Ul(v_if="hs_notes.length", classes="vc-details-notes mx-3"):
                html.Li("{{ n }}", v_for="n in hs_notes", key="n")

        # ═══ shape metrics ═══════════════════════════════════════════════════
        with html.Div(v_if="mode === 'shape'", classes="vc-scroll"):
            with v.VExpansionPanels(v_model=("shape_a",), multiple=True, flat=True,
                                    classes="vc-panels"):
                # ── Surface ──────────────────────────────────────────────────
                with v.VExpansionPanel(value="surface",
                                       raw_attrs=['data-test="sec-surface"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-cube-outline", "Surface", done="sm_info")
                        html.Span("{{ sm_info.faces.toLocaleString() }} faces · Euler "
                                  "{{ sm_info.euler }}", v_if="sm_info",
                                  classes="vc-pill")
                    with v.VExpansionPanelText():
                        _tip_btn("Measures the surface currently on screen in the "
                                 "Cut section.", "Use the cut from the Cut section",
                                 block=True, click=ctrl.sm_use_cut, variant="flat",
                                 classes="vc-btn-primary",
                                 prepend_icon="mdi-content-cut",
                                 disabled=("!loaded",),
                                 raw_attrs=['data-test="sm-use-cut"'])
                        html.Div("Open a case in the Cut section first.",
                                 v_if="!loaded", classes="vc-cap")
                        html.Div("For a ventricle, use δK: inner or δK: outer in the "
                                 "Cardiac section.", classes="vc-cap")
                        if PUBLIC:
                            html.Div(TXT_SM_PUBLIC, classes="vc-cap")
                        else:
                            html.Div("Or open a surface on this machine: any open "
                                     "triangular mesh, from this tool or not.",
                                     classes="vc-cap mt-3")
                            with _fld("Surface"):
                                v.VTextField(v_model=("sm_path",), density="compact",
                                             variant="outlined", hide_details=True,
                                             prepend_inner_icon="mdi-file-outline",
                                             placeholder=".stl / .ply / .obj",
                                             raw_attrs=['data-test="sm-path"'])
                            v.VBtn("Load surface", block=True,
                                   click=ctrl.sm_load_path, variant="outlined",
                                   classes="vc-btn-secondary",
                                   prepend_icon="mdi-tray-arrow-up",
                                   raw_attrs=['data-test="sm-load"'])
                        with html.Div(v_if="sm_info", classes="readout vc-rgrid",
                                      raw_attrs=['data-test="readout-surface"']):
                            _readout([
                                ("file", "File", "sm_info.name", "", None),
                                ("faces", "Faces", "sm_info.faces.toLocaleString()", "",
                                 None),
                                ("openings", "Openings", "sm_info.openings", "", None),
                                ("euler", "Euler", "sm_info.euler", "", None),
                                ("area", "Area",
                                 "Math.round(sm_info.area_mm2).toLocaleString()", "mm²",
                                 None),
                                ("folds", "Folds", "sm_info.folds === null ? "
                                 "'not checked' : sm_info.folds", "", None),
                            ])
                        # CAMP was built for open geometries. A closed one still
                        # partitions, so nothing here would fail; it just answers a
                        # different question, which is worse.
                        v.VAlert("This surface is closed. CAMP expects an open "
                                 "geometry, so cut it first.", icon="mdi-alert",
                                 v_if="sm_info && sm_info.closed", type="warning",
                                 density="compact", variant="tonal",
                                 classes="vc-callout warn mt-3")
                        # Known before the run, so said before the run.
                        v.VAlert("{{ sm_info.folds }} folded edge(s), worst "
                                 "dihedral {{ Math.round(sm_info.max_dihedral) }}"
                                 "°. The surface doubles back through itself "
                                 "there: segmentation geometry, not anatomy. δK "
                                 "is still computed; the folds are not corrected "
                                 "for, because across the reference cohort their "
                                 "count does not predict how much δK moves.",
                                 v_if="sm_info && sm_info.folds", icon="mdi-alert",
                                 type="warning", density="compact",
                                 variant="tonal", classes="vc-callout warn mt-3")

                # ── Partition ────────────────────────────────────────────────
                with v.VExpansionPanel(value="partition",
                                       raw_attrs=['data-test="sec-partition"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-vector-polyline", "Partition")
                        with html.Span(classes=("sm_faces === sm_ref_faces_n && sm_trim "
                                                "=== sm_ref_trim ? 'vc-pill ok' "
                                                ": 'vc-pill warn'",)):
                            v.VIcon(icon=("sm_faces === sm_ref_faces_n && sm_trim === "
                                          "sm_ref_trim ? 'mdi-check' : 'mdi-alert-outline'",),
                                    size="14")
                            html.Span("{{ Math.round(+sm_faces / 1000) }}k · trim "
                                      "{{ sm_trim }} · {{ sm_seeds }} "
                                      "{{ sm_seeds === '1' ? 'partition' : 'partitions' }}")
                    with v.VExpansionPanelText():
                        for lbl, key, opts, tid, ref in [
                            ("Faces", "sm_faces",
                             [(f"{n // 1000}k", str(n), None) for n in FACE_TARGETS],
                             "sm-faces", str(SM.REFERENCE["target_faces"])),
                            ("Trim", "sm_trim",
                             [(f"{t:g}", f"{t:g}", None) for t in TRIM_OPTIONS],
                             "sm-trim", f"{SM.REFERENCE['trim_frac']:g}"),
                            ("Partitions", "sm_seeds",
                             [(str(n), str(n), None) for n in SEED_OPTIONS],
                             "sm-seeds", "1"),
                        ]:
                            with html.Div(classes="vc-row-lbl"):
                                html.Span(lbl, classes="vc-lbl")
                                _seg(key, opts, test=tid, disabled="sm_busy", dots=ref)
                        # A ventricle wall has no reference cohort: the one here is 30
                        # aortas, so matching its settings would say nothing.
                        with html.Div(v_if="sm_info && sm_info.kind === 'heart'",
                                      classes="vc-callout-div warn mt-3",
                                      raw_attrs=['data-test="sm-ref-heart"']):
                            v.VIcon("mdi-heart-outline", size="18")
                            html.Div("A ventricle wall. The reference cohort is "
                                     "aortic, so there is none for this surface: "
                                     "compare only with other ventricle walls run at "
                                     "the same face count and trim.")
                        # Whether a number from these settings can be set beside the
                        # reference cohort's at all. Face count especially.
                        with html.Div(v_if="!(sm_info && sm_info.kind === 'heart')",
                                      classes=("sm_faces === sm_ref_faces_n && sm_trim "
                                               "=== sm_ref_trim ? 'vc-callout-div ok mt-3' "
                                               ": 'vc-callout-div warn mt-3'",),
                                      raw_attrs=['data-test="sm-ref"',
                                                 ':data-match="sm_faces === sm_ref_faces_n '
                                                 '&& sm_trim === sm_ref_trim"']):
                            v.VIcon(icon=("sm_faces === sm_ref_faces_n && sm_trim === "
                                          "sm_ref_trim ? 'mdi-check-circle' "
                                          ": 'mdi-alert'",), size="18")
                            html.Div("Matches the reference cohort ({{ sm_ref_n }} "
                                     "cases, {{ sm_ref_faces }} faces, trim "
                                     "{{ sm_ref_trim }}): values are comparable.",
                                     v_if="sm_faces === sm_ref_faces_n && "
                                          "sm_trim === sm_ref_trim")
                            html.Div("Differs from the reference cohort "
                                     "({{ sm_ref_faces }} faces, trim "
                                     "{{ sm_ref_trim }}): compare only with runs "
                                     "at the same settings.",
                                     v_if="sm_faces !== sm_ref_faces_n || "
                                          "sm_trim !== sm_ref_trim")
                        html.Div("{{ sm_repeat }}", v_if="sm_repeat",
                                 classes="vc-cap mt-2",
                                 raw_attrs=['data-test="sm-repeat"'])
                        with html.Details(classes="vc-disclose"):
                            with html.Summary():
                                v.VIcon("mdi-help-circle-outline", size="18")
                                html.Span("What these settings do")
                                v.VIcon("mdi-chevron-down", size="18",
                                        classes="vc-disclose-chev")
                            with html.Div(classes="vc-more"):
                                html.P("CAMP splits the surface until every patch "
                                       "takes a converged Monge fit, so there is no "
                                       "patch count to choose. These two are what "
                                       "delta K moves with.")
                                # The single most consequential control on this panel,
                                # and the least obvious. delta K is carried by ~22
                                # patches out of ~850 and k-means places their
                                # boundaries arbitrarily, so one partition is a draw.
                                html.P("One partition is one draw: re-seeding the same "
                                       "mesh moves delta K by about "
                                       "{{ sm_ref_partition_cv }}%. Averaging n "
                                       "independent partitions divides that by root n, "
                                       "and costs n times the runtime.")
                                html.P("The trim discards the highest and lowest "
                                       "integrated-K patches -- both tails, so the "
                                       "number typed is per tail.")
                                html.P("Reference cohort: {{ sm_ref_n }} cases at "
                                       "{{ sm_ref_faces }} faces, trim "
                                       "{{ sm_ref_trim }}. delta K is extensive and "
                                       "partition dependent, so two values are only "
                                       "comparable when both settings match. Face "
                                       "count especially: delta K drifts with it "
                                       "because min_patch_size is a vertex count, and "
                                       "averaging partitions does not remove that.")

            # ── the result: not a panel, so it is never folded away ──────────
            with html.Div(classes="vc-result", raw_attrs=['data-test="sec-result"']):
                with html.Div(classes="vc-result-h"):
                    with html.Span(classes="vc-ptile accent"):
                        v.VIcon("mdi-chart-bell-curve", size="17")
                    html.Span("Result", classes="vc-ptitle")
                    v.VSpacer()
                    v.VChip("{{ sm_busy ? 'Computing ' + sm_elapsed.toFixed(0) + ' s' : "
                            "(sm_result ? 'δK computed' : (sm_bad ? 'Failed' : "
                            "(sm_info ? 'Ready' : 'No surface'))) }}",
                            variant="flat", size="small",
                            classes=("sm_busy ? 'vc-verdict idle' : (sm_result ? "
                                     "'vc-verdict ok' : (sm_bad ? 'vc-verdict bad' "
                                     ": 'vc-verdict idle'))",),
                            prepend_icon=("sm_busy ? 'mdi-progress-clock' : (sm_result "
                                          "? 'mdi-check-circle' : (sm_bad ? "
                                          "'mdi-alert-circle' : undefined))",),
                            raw_attrs=['data-test="verdict"'])
                html.Pre("No surface yet. Use the cut from the Cut section"
                         + ("." if PUBLIC else ", or open a surface on this machine."),
                         v_if="!sm_info", classes="empty")
                html.Pre("Ready. Compute δK to partition this surface.",
                         v_if="sm_info && !sm_result && !sm_busy", classes="empty")

                # While it runs. The phase matters more than a spinner: CAMP is silent
                # for most of a run, so without it the panel looks stuck at ten seconds.
                with html.Div(v_if="sm_busy", classes="readout vc-running"):
                    v.VProgressLinear(indeterminate=True, height="4",
                                      color="primary", bg_color="#E2E0FF",
                                      rounded=True, classes="mb-2")
                    with html.Div(classes="ro-row vc-rrow"):
                        html.Span("Elapsed", classes="ro-lbl")
                        html.Span("{{ sm_elapsed.toFixed(0) }}", classes="ro-val")
                        html.Span("s", classes="ro-u")
                    with html.Div(classes="vc-res-line"):
                        v.VIcon("mdi-progress-clock", size="15")
                        html.Span("{{ sm_phase }}")
                    with html.Div(classes="vc-row-between"):
                        v.VSpacer()
                        v.VBtn("Cancel", click=ctrl.sm_cancel, size="small",
                               variant="text", color="primary", classes="vc-btn-text",
                               raw_attrs=['data-test="sm-cancel"'])

                with html.Div(v_if="sm_result", classes="readout",
                              raw_attrs=['data-test="readout-result"']):
                    # The band rides in the same row as the value, at the same size. It
                    # is not a caption.
                    with html.Div(classes="ro-row vc-res-v", raw_attrs=['data-key="dk"']):
                        html.Span("δK", classes="ro-lbl")
                        html.Span("{{ sm_result.k_fluct.toPrecision(4) }}",
                                  classes="ro-val")
                        html.Span("± {{ sm_result.band_pct }}%", classes="ro-u vc-band")
                    with html.Div(classes="vc-stats"):
                        _readout([
                            ("area", "Area",
                             "Math.round(sm_result.A_mm2).toLocaleString()", "mm²", None),
                            ("patches", "Patches", "sm_result.n_patches.toLocaleString()",
                             "", None),
                            ("kept", "Kept", "sm_result.n_kept.toLocaleString()",
                             "after trim", None),
                            ("floor", "At floor", "sm_result.pct_floor.toFixed(1)",
                             "% unsplittable", None),
                            ("fail", "Fit fail", "sm_result.pct_failed.toFixed(1)",
                             "% below R2/eta", None),
                            ("faces", "Faces", "sm_result.faces.toLocaleString()", "",
                             None),
                            ("euler", "Euler", "sm_result.euler", "", None),
                            ("took", "Took", "sm_result.seconds.toFixed(0)", "s", None),
                        ], cls="vc-stat-cell")
                with html.Div(v_if="sm_result", classes="vc-res-notes"):
                    with html.Div(classes="vc-res-line"):
                        v.VIcon("mdi-information-outline", size="15")
                        html.Span("{{ sm_headline }}")
                    # Directly under the value, because the value on its own invites a
                    # comparison this number cannot support.
                    with html.Div(v_if="sm_repeat", classes="vc-res-line"):
                        v.VIcon("mdi-swap-vertical", size="15")
                        html.Span("{{ sm_repeat }}")
                    # Everything that would stop this number being comparable to
                    # another. In full: the conditions are half of the value.
                    v.VAlert("{{ n.text }}", v_for="n in sm_notes", key="n.text",
                             density="compact", variant="tonal",
                             classes=("n.level === 'warn' ? 'vc-callout warn' "
                                      ": 'vc-callout info'",),
                             type=("n.level === 'warn' ? 'warning' : 'info'",),
                             icon=("n.level === 'warn' ? 'mdi-alert' "
                                   ": 'mdi-information-outline'",),
                             raw_attrs=['data-test="sm-note"', ':data-level="n.level"'])

            with v.VExpansionPanels(v_model=("shape_b",), multiple=True, flat=True,
                                    classes="vc-panels mt-2"):
                # ── View ─────────────────────────────────────────────────────
                with v.VExpansionPanel(value="view",
                                       raw_attrs=['data-test="sec-shape-view"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-palette-outline", "View")
                        html.Span("{{ sm_scale === 'robust' ? 'Robust · clip ' + "
                                  "sm_clip + '%' : (sm_scale === 'full' ? 'Full range' "
                                  ": 'Manual') }}", classes="vc-pill")
                    with v.VExpansionPanelText():
                        _switch("sm_outlines", "Patch outlines", tip=TXT_PLASMA,
                                aria="How patches are coloured",
                                disabled=("!sm_result",))
                        _switch("sm_dim_trimmed", "Grey out trimmed patches",
                                disabled=("!sm_result",))
                        # Integrated K is long-tailed: a handful of patches carry
                        # values many times the rest, so the full range puts almost
                        # every patch in the top few per cent of the map and the vessel
                        # reads as one flat colour. Robust is the default for that
                        # reason; Full is CAMP's own auto-scale, kept because a figure
                        # meant to sit beside theirs needs it.
                        with html.Div(classes="vc-row-lbl mt-2"):
                            html.Span("Scale", classes="vc-lbl")
                            _seg("sm_scale", [("Robust", "robust", None),
                                              ("Full", "full", None),
                                              ("Manual", "manual", None)],
                                 test="sm-scale", disabled="!sm_result")
                        with html.Div(classes="op-row vc-clip",
                                      v_if="sm_scale === 'robust'"):
                            html.Span("Clip", classes="vc-lbl")
                            v.VSlider(v_model=("sm_clip",), min="0", max="20",
                                      step="0.5", density="compact",
                                      hide_details=True, track_size="4",
                                      thumb_size="16", color="primary",
                                      disabled=("!sm_result",))
                            _num(dict(v_model_number=("sm_clip",)), "%",
                                 disabled=("!sm_result",))
                        with html.Div(classes="vc-grid2 mt-2",
                                      v_if="sm_scale === 'manual'"):
                            for lbl, key in (("Min", "sm_cmin"), ("Max", "sm_cmax")):
                                with html.Div(classes="vc-cell-h"):
                                    html.Span(lbl, classes="vc-lbl")
                                    _num(dict(v_model_number=(key,)), "",
                                         change=ctrl.sm_apply_bounds,
                                         disabled=("!sm_result",))
                        html.Div("Data range {{ sm_full }}",
                                 v_if="sm_full && sm_result", classes="vc-cap mt-2")
                        html.Div("{{ sm_range }}", v_if="sm_range && sm_result",
                                 classes="vc-cap")

                # ── Output ───────────────────────────────────────────────────
                with v.VExpansionPanel(value="out",
                                       raw_attrs=['data-test="sec-shape-output"']):
                    with v.VExpansionPanelTitle():
                        _ptitle("mdi-export-variant", "Output")
                    with v.VExpansionPanelText():
                        with _fld("Name"):
                            v.VTextField(v_model=("sm_export_name",),
                                         density="compact", variant="outlined",
                                         hide_details=True)
                        # It writes to the server's own disk, so it is not offered
                        # when hosted -- sm_export refuses it there as well.
                        if PUBLIC:
                            html.Div("The patch table is written only by the "
                                     "desktop and local versions of this tool.",
                                     classes="vc-cap")
                        else:
                            with html.Div(classes="vc-btnrow"):
                                v.VBtn("Export patch table…", click=ctrl.sm_export,
                                       block=True, variant="outlined",
                                       classes="vc-btn-secondary",
                                       prepend_icon="mdi-content-save-outline",
                                       disabled=("!sm_result",),
                                       raw_attrs=['data-test="sm-export"'])
                            with html.Div(classes="vc-cap vc-lbl"):
                                html.Span("A row per patch, and a sidecar")
                                _info(TXT_SIDECAR, "What the export writes")
                        v.VAlert("{{ sm_export_msg }}", v_if="sm_export_msg",
                                 density="compact", variant="tonal",
                                 classes=("sm_export_bad ? 'vc-callout bad mt-2' "
                                          ": 'vc-callout ok mt-2'",),
                                 type=("sm_export_bad ? 'error' : 'success'",))

    # ── the 3-D view ─────────────────────────────────────────────────────────
    with layout.content:
        # vc-stage makes the container the overlays' positioning box. The view stays
        # its direct child, as it always was: wrapped in a div of its own it would
        # take a percentage height from a parent that has only a min-height, and
        # collapse.
        with v.VContainer(fluid=True, classes="fill-height pa-0 ma-0 vc-stage"):
            # LeftButtonPress comes back with the click position, which is what
            # turns the view from a picture into something you can measure on.
            # Gated by the Pin switch, so a click meant to turn the model does not
            # also drop a mark.
            view = vtkw.VtkRemoteView(
                SCENE.plotter.ren_win, interactive_ratio=1,
                interactor_events=("view_events", ["LeftButtonPress"]),
                LeftButtonPress="(e) => pin_mode && trigger('view_click', "
                                "[e.position ? [e.position.x, e.position.y] "
                                ": null])",
            )
            ctrl.view_update = view.update
            ctrl.view_reset_camera = view.reset_camera

            # Everything over the canvas, in one wrapper the pixel checks hide
            # before they sample it. Nothing goes in the bottom-left corner, where
            # the L/P/S triad is.
            with html.Div(classes="vc-overlays",
                          raw_attrs=['data-test="view-overlay"']):
                # ── top left: what is selected, and what is wrong ───────────
                with html.Div(v_if="mode === 'cut' && loaded", classes="vc-tl"):
                    with html.Div(v_if="!selected && $vuetify.display.mdAndUp",
                                  classes="vc-glass vc-hint"):
                        v.VIcon("mdi-cursor-default-click-outline", size="18")
                        html.Span("Select an opening to inspect its cut")
                    with html.Div(v_if="selected && $vuetify.display.mdAndUp && "
                                       "Object.values(openings || {}).find(o => o.label === selected)",
                                  classes="vc-glass vc-hud"):
                        html.Span(classes="vc-glyph disc sel")
                        html.B("{{ Object.values(openings || {}).find(o => o.label === selected).name }}")
                        html.I(classes="vc-hud-sep d-none d-md-inline-block")
                        with html.Span(classes="d-none d-md-inline-flex"):
                            v.VIcon("mdi-diameter-variant", size="15")
                            html.Span("{{ (+Object.values(openings || {}).find(o => o.label === selected)"
                                      ".diameter || 0).toFixed(1) }} mm")
                        html.I(classes="vc-hud-sep d-none d-md-inline-block")
                        html.Span("{{ (+Object.values(openings || {}).find(o => o.label === selected)"
                                  ".offset).toFixed(1) }} mm from tip",
                                  classes="d-none d-md-inline-flex")
                        html.I(classes="vc-hud-sep d-none d-md-inline-block")
                        html.Span("{{ Object.values(openings || {}).find(o => o.label === selected).mode === "
                                  "'normal' ? 'Perpendicular' : 'Axial' }}",
                                  classes="d-none d-md-inline-flex")
                    with html.Div(v_if="!cut_ok", classes="vc-statechip"):
                        v.VIcon("mdi-alert-circle-outline", size="17")
                        html.Span("Surface greyed until the cut is valid")
                    with html.Div(v_if="measure && m_show && $vuetify.display.mdAndUp",
                                  classes="vc-glass vc-hint"):
                        html.Span(classes="vc-glyph mplane")
                        html.Span("Measure plane at {{ (+m_station).toFixed(0) }} mm")
                        v.VBtn("Hide plane", variant="text", size="x-small",
                               color="primary", classes="vc-btn-text",
                               click="m_show = false; trigger('measure_show', [false])")

                # ── top right: the view's own tools ─────────────────────────
                with html.Div(v_if="(mode === 'cut' && loaded) || "
                                   "(mode === 'shape' && sm_info) || "
                                   "(mode === 'heart' && hs_loaded)",
                              classes="vc-glass vc-viewbar",
                              raw_attrs=['data-test="view-toolbar"']):
                    _tip_btn("Fit view", icon="mdi-fit-to-screen-outline",
                             variant="text", size="small", classes="vc-vt-ic",
                             click=ctrl.view_reset_camera,
                             raw_attrs=['aria-label="Fit view"'])
                    html.Span(classes="vc-vt-sep d-none d-md-inline-block")
                    with v.VBtn(v_if="mode === 'cut'", variant="text",
                                classes=("pin_mode ? 'vc-vt on' : 'vc-vt'",),
                                click="pin_mode = !pin_mode",
                                raw_attrs=['aria-label="Pin points on click"',
                                           ':aria-pressed="pin_mode"',
                                           'title="Pin a point on click"']):
                        v.VIcon("mdi-map-marker-plus-outline", size="18")
                        html.Span("Pin", classes="vc-vt-t")
                        html.Span(classes=("pin_mode ? 'vc-minisw on' "
                                           ": 'vc-minisw'",))
                    with v.VBtn(v_if="mode === 'cut'", variant="text", classes="vc-vt",
                                click=OPEN_DIAMETER,
                                raw_attrs=['aria-label="Diameter along the vessel"',
                                           'title="Diameter along the vessel"']):
                        v.VIcon("mdi-ruler", size="18")
                        html.Span("Diameter", classes="vc-vt-t")
                    with v.VMenu(v_if="mode === 'cut' && $vuetify.display.smAndDown",
                                 location="start"):
                        with html.Template(v_slot_activator=("{ props }",)):
                            v.VBtn(icon="mdi-palette-outline", v_bind=("props",),
                                   variant="text", classes="vc-vt-ic",
                                   raw_attrs=['aria-label="Colour key"'])
                        with v.VCard(classes="vc-menu vc-keymenu"):
                            html.Div("Colour key", classes="vc-menu-h")
                            for glyph, label in [("disc", "Cut disc"),
                                                 ("disc sel", "Selected cut"),
                                                 ("rim", "Opening rim"),
                                                 ("cl", "Centerline")]:
                                with html.Div(classes="vc-ki"):
                                    html.Span(classes=f"vc-glyph {glyph}")
                                    html.Span(label)
                    html.Span(classes="vc-vt-sep d-none d-md-inline-block",
                              v_if="mode === 'cut'")
                    with v.VBtn(variant="text", v_if="$vuetify.display.mdAndUp",
                                classes=("show_key ? 'vc-vt on' : 'vc-vt'",),
                                click="show_key = !show_key",
                                raw_attrs=['aria-label="Colour key"',
                                           ':aria-pressed="show_key"',
                                           'data-test="colour-key"']):
                        v.VIcon("mdi-palette-outline", size="18")
                        html.Span("Key", classes="vc-vt-t")

                # ── bottom centre: the colour key ────────────────────────────
                with html.Div(v_if="mode === 'cut' && loaded && show_key && "
                                   "$vuetify.display.mdAndUp",
                              classes="vc-glass vc-legend"):
                    html.Span("Key", classes="vc-ki-h")
                    for glyph, label, cond in [
                            ("disc", "Cut disc", None),
                            ("disc sel", "Selected cut", None),
                            ("rim", "Opening rim", None),
                            ("cl", "Centerline", None),
                            ("mplane", "Measure plane", "measure && m_show"),
                            ("mark", "Pinned point", "marks.length")]:
                        extra = dict(v_if=cond) if cond else {}
                        with html.Span(classes="vc-ki", **extra):
                            html.Span(classes=f"vc-glyph {glyph}")
                            html.Span(label)
                    with html.Span(v_for="z in zone_legend", key="z.zone",
                                   v_if="colour_zones", classes="vc-ki"):
                        html.Span(classes="vc-glyph zone",
                                  style=("`background:${z.colour}`",))
                        html.Span("Zone {{ z.zone }}")
                with html.Div(v_if="mode === 'shape' && sm_result && sm_full && show_key",
                              classes="vc-glass vc-cbar"):
                    with html.Div(classes="vc-cbar-h"):
                        html.Span("Integrated K per patch")
                        html.Span("{{ sm_scale === 'robust' ? 'Robust scale · clip ' + "
                                  "sm_clip + '%' : (sm_scale === 'full' ? 'Full range' "
                                  ": 'Manual scale') }}", classes="vc-cbar-m")
                    html.Div(classes="vc-cbar-g")
                    with html.Div(classes="vc-cbar-t"):
                        html.Span("{{ (+sm_cmin).toPrecision(3) }}")
                        with html.Span(v_if="sm_dim_trimmed", classes="vc-cbar-gr"):
                            html.I()
                            html.Span("Trimmed")
                        html.Span("{{ (+sm_cmax).toPrecision(3) }}")

                # ── pins: what a click will do, and what it did ──────────────
                with html.Div(v_if="mode === 'cut' && loaded && marks_note && "
                                   "$vuetify.display.mdAndUp",
                              classes="vc-glass vc-pinchip"):
                    v.VIcon("mdi-map-marker-plus-outline", size="17")
                    html.Span("{{ marks_note || 'Click the surface to pin a point' }}")
                    v.VBtn("Send cut here", v_if="marks_note", variant="text",
                           size="x-small", color="primary", classes="vc-btn-text",
                           click="trigger('send_cut_to_mark')")
                    v.VBtn("Clear", v_if="marks_note", variant="text", size="x-small",
                           classes="vc-btn-quiet", click="trigger('clear_marks')")

                # ── phone: the counts, and the way to the controls ──────────
                with html.Div(v_if="mode === 'cut' && loaded && "
                                   "$vuetify.display.smAndDown",
                              classes="vc-glass vc-statpill"):
                    html.Span("{{ n_openings }} openings")
                    html.I(classes="vc-hud-sep")
                    with html.Span(classes=("euler === 2 - n_openings ? 'ok' : 'bad'",)):
                        v.VIcon(icon=("euler === 2 - n_openings ? 'mdi-check' "
                                      ": 'mdi-alert-outline'",), size="15")
                        html.Span("Euler {{ euler }}")
                with v.VBtn(v_if="mode === 'cut' && loaded && $vuetify.display.smAndDown",
                            variant="flat", classes="vc-glass vc-adjust",
                            prepend_icon="mdi-tune-variant",
                            click="main_drawer = true; cut_panels = ['cuts']"):
                    html.Span("Adjust cuts")

                # ── nothing to look at yet ───────────────────────────────────
                if PUBLIC:
                    # The notice and the pickers, over the empty view, before the
                    # first case. Words unchanged; the pickers below it, never above.
                    with html.Div(v_if="mode === 'cut' && !loaded",
                                  classes="vc-consent-wrap"):
                        with html.Div(classes="vc-consent"):
                            _notice("card")
                            with html.Div(classes="vc-consent-pickers"):
                                _pickers()
                else:
                    with html.Div(v_if="mode === 'cut' && !loaded", classes="vc-empty",
                                  raw_attrs=['data-test="view-empty"']):
                        with html.Span(classes="vc-empty-ic"):
                            v.VIcon(icon=("busy ? 'mdi-progress-clock' "
                                          ": 'mdi-cube-outline'",), size="24")
                        html.Div("{{ busy ? 'Loading the case' : 'No case open' }}",
                                 classes="vc-empty-t")
                        html.P("The page stays usable while it does.", v_if="busy")
                        html.P(TXT_EMPTY_LOCAL, v_if="!busy")
                        v.VBtn("Open the named files", v_if="!busy && mat_path && csv_path",
                               variant="flat", classes="vc-btn-primary",
                               prepend_icon="mdi-folder-open-outline",
                               click=ctrl.load_model,
                               raw_attrs=['data-test="open-named"'])
                        v.VBtn("Open controls", v_if="$vuetify.display.smAndDown",
                               variant="tonal", color="primary",
                               classes="vc-btn-feature mt-2",
                               click="main_drawer = true")
                with html.Div(v_if="mode === 'heart' && hs_loaded && show_key && "
                                   "$vuetify.display.mdAndUp",
                              classes="vc-glass vc-legend"):
                    html.Span("Key", classes="vc-ki-h")
                    for colour, label in [("#4c78a8", "Inner wall"),
                                          ("#d8b7a0", "Outer wall"),
                                          ("#9e9e9e", "Removed"),
                                          ("#d55e00", "Cut plane"),
                                          ("#0070c0", "Basal cut edge"),
                                          ("#d62728", "Hole")]:
                        with html.Span(classes="vc-ki"):
                            html.Span(classes="vc-glyph zone",
                                      style=f"background:{colour}")
                            html.Span(label)
                with html.Div(v_if="mode === 'heart' && !hs_loaded", classes="vc-empty",
                              raw_attrs=['data-test="view-empty-heart"']):
                    with html.Span(classes="vc-empty-ic"):
                        v.VIcon(icon=("hs_busy ? 'mdi-progress-clock' "
                                      ": 'mdi-heart-outline'",), size="24")
                    html.Div("{{ hs_busy ? 'Splitting the ventricle' "
                             ": 'No ventricle open' }}", classes="vc-empty-t")
                    html.P("Choose a closed myocardium mesh under Ventricle, then "
                           "Split this ventricle.", v_if="!hs_busy")
                    v.VBtn("Split the chosen mesh", v_if="!hs_busy && hs_path",
                           variant="flat", classes="vc-btn-primary",
                           prepend_icon="mdi-call-split", click=ctrl.hs_load,
                           raw_attrs=['data-test="hs-open-named"'])
                with html.Div(v_if="mode === 'shape' && !sm_info", classes="vc-empty",
                              raw_attrs=['data-test="view-empty"']):
                    with html.Span(classes="vc-empty-ic"):
                        v.VIcon("mdi-cube-outline", size="24")
                    html.Div("No surface yet", classes="vc-empty-t")
                    html.P("Use the cut from the Cut section, under Surface."
                           + ("" if PUBLIC else " Or open a surface on this machine."))
                    v.VBtn("Open controls", v_if="$vuetify.display.smAndDown",
                           variant="tonal", color="primary",
                           classes="vc-btn-feature mt-2", click="main_drawer = true")

    # ── footer: the status line, Details, Reset ──────────────────────────────
    with layout.footer as footer:
        footer.clear()
        footer.classes = "vc-footer"
        with html.Div(classes="vc-fstatus"):
            html.Span(classes=("mode === 'shape' ? (sm_bad ? 'vc-sdot bad' : (sm_result "
                               "? 'vc-sdot' : 'vc-sdot idle')) : (mode === 'heart' ? "
                               "(hs_bad ? 'vc-sdot bad' : (hs_loaded ? 'vc-sdot' : "
                               "'vc-sdot idle')) : (status_bad ? "
                               "'vc-sdot bad' : (loaded ? 'vc-sdot' : 'vc-sdot idle')))",))
            # .statusbar holds the server's words and nothing else: the checks read it
            # whole, and look for text-error on it.
            html.Span("{{ status_text }}", v_if="mode === 'cut'",
                      classes=("status_bad ? 'statusbar text-error' : 'statusbar'",),
                      raw_attrs=['data-test="status"',
                                 ':data-status="status_bad ? \'bad\' : \'ok\'"'])
            html.Span("{{ sm_msg || '" + TXT_SM_DEFAULT + "' }}", v_if="mode === 'shape'",
                      classes=("sm_bad ? 'statusbar text-error' : 'statusbar'",),
                      raw_attrs=['data-test="status"'])
            html.Span("{{ hs_msg || 'Choose a ventricle under Ventricle, then split it.' }}",
                      v_if="mode === 'heart'",
                      classes=("hs_bad ? 'statusbar text-error' : 'statusbar'",),
                      raw_attrs=['data-test="status"'])
            with v.VMenu(v_if="mode === 'cut' && loaded && $vuetify.display.mdAndUp",
                         location="top end", close_on_content_click=False):
                with html.Template(v_slot_activator=("{ props }",)):
                    with v.VBtn(v_bind=("props",), variant="text", size="small",
                                classes="vc-btn-quiet",
                                prepend_icon="mdi-text-box-outline",
                                raw_attrs=['aria-label="Status details"',
                                           'data-test="status-details"']):
                        html.Span("Details", classes="d-none d-sm-inline")
                        html.Span("{{ status_notes.length }}",
                                  v_if="status_notes.length",
                                  classes=("status_bad ? 'vc-count bad' : 'vc-count'",))
                with v.VCard(classes="vc-menu vc-details", max_width="460"):
                    html.Div("Status", classes="vc-menu-h")
                    html.P("{{ status_text }}",
                           classes=("status_bad ? 'vc-details-main bad' "
                                    ": 'vc-details-main'",))
                    with html.Ul(v_if="status_notes.length", classes="vc-details-notes"):
                        html.Li("{{ n }}", v_for="n in status_notes", key="n")
        with html.Div(classes="vc-factions"):
            # On a phone the verdict lives here, the only one on screen.
            v.VChip("{{ busy ? 'Working' : (cut_ok ? 'Cut valid' : 'Cut not valid') }}",
                    v_if="mode === 'cut' && loaded && $vuetify.display.smAndDown",
                    variant="flat",
                    classes=("busy ? 'vc-verdict idle' : (cut_ok ? 'vc-verdict ok' "
                             ": 'vc-verdict bad')",),
                    prepend_icon=("busy ? 'mdi-progress-clock' : (cut_ok ? "
                                  "'mdi-check-circle' : 'mdi-alert-circle')",),
                    raw_attrs=['data-test="verdict"'])
            v.VSpacer(classes="d-md-none")
            with v.VMenu(v_if="mode === 'cut' && loaded && $vuetify.display.mdAndUp",
                         location="top end"):
                with html.Template(v_slot_activator=("{ props }",)):
                    with v.VBtn(v_bind=("props",), variant="text", size="small",
                                classes="vc-btn-quiet", prepend_icon="mdi-restore",
                                disabled=("!loaded",), raw_attrs=['data-test="reset"']):
                        html.Span("Reset cuts…", classes="d-none d-sm-inline")
                with v.VCard(classes="vc-menu vc-confirm", max_width="340"):
                    html.P(TXT_RESET_CONFIRM)
                    with html.Div(classes="vc-confirm-row"):
                        v.VBtn("Keep", variant="text", size="small",
                               classes="vc-btn-quiet")
                        v.VBtn("Reset", variant="flat", size="small",
                               classes="vc-btn-danger", click=ctrl.reset_all,
                               raw_attrs=['data-test="reset-confirm"'])
            # On a phone: Details and Reset in one menu, and Download labelled, in the
            # same state as the dock's.
            with v.VMenu(v_if="mode === 'cut' && $vuetify.display.smAndDown",
                         location="top end", close_on_content_click=False):
                with html.Template(v_slot_activator=("{ props }",)):
                    v.VBtn(icon="mdi-dots-horizontal", v_bind=("props",), variant="flat",
                           classes="vc-btn-secondary vc-phone-more",
                           raw_attrs=['aria-label="Details and Reset cuts"'])
                with v.VCard(classes="vc-menu vc-details", max_width="340"):
                    html.Div("Status", classes="vc-menu-h")
                    html.P("{{ status_text }}",
                           classes=("status_bad ? 'vc-details-main bad' "
                                    ": 'vc-details-main'",))
                    with html.Ul(v_if="status_notes.length", classes="vc-details-notes"):
                        html.Li("{{ n }}", v_for="n in status_notes", key="n")
                    with html.Div(v_if="loaded", classes="vc-menu-sw"):
                        html.P(TXT_RESET_CONFIRM, classes="vc-details-main mt-2")
                        with html.Div(classes="vc-confirm-row"):
                            v.VBtn("Reset cuts", variant="flat", size="small",
                                   prepend_icon="mdi-restore", classes="vc-btn-danger",
                                   click=ctrl.reset_all)
            with v.VBtn(variant="flat",
                        v_if="mode === 'cut' && loaded && $vuetify.display.smAndDown",
                        classes=("cut_ok ? 'vc-btn-primary vc-phone-dl' "
                                 ": 'vc-btn-bad vc-phone-dl'",),
                        prepend_icon=("cut_ok ? 'mdi-tray-arrow-down' "
                                      ": 'mdi-alert-circle-outline'",),
                        click=DOWNLOAD_JS,
                        raw_attrs=['aria-label="Download open mesh (.zip)"']):
                html.Span("Download")
            html.Span(f"build {BUILD}", classes="vc-build",
                      raw_attrs=['data-test="build"'])


if __name__ == "__main__":
    # The port is overridable because trame's built-in server holds ONE shared state
    # per process: every tab connected to it sees and changes the same case. Hosted,
    # the gate above admits one visitor per process, so a second process is the only
    # way to serve two people at once. Set VC_CUT_PORT to run one beside another
    # rather than having the second fail to bind, which surfaces as the page never
    # mounting rather than as an error anyone can read.
    port = int(os.environ.get("VC_CUT_PORT") or 8080)
    host = (os.environ.get("VC_CUT_HOST") or "127.0.0.1").strip()

    # Refuse, rather than warn. Binding a routable address without hosted mode puts
    # the path boxes on the network, and they read any file the server account can
    # read. There is no authentication anywhere in this app. A warning printed to a
    # terminal nobody is watching would not be a control.
    if host not in LOOPBACK and not PUBLIC:
        sys.exit(
            f"Refusing to bind {host}. That address is reachable from the network, "
            f"this app has no authentication, and its input is a file PATH, so "
            f"anything the server account can read would be readable remotely.\n"
            f"Use 127.0.0.1, or set VC_PUBLIC=1 for the hosted configuration.")
    # Uploads are kept, and the page has to say for how long, who sees them and how
    # to have them removed before anyone uploads. Until the PI has decided those, the
    # notice carries bracketed placeholders -- fine on this machine for testing, not
    # on an address a stranger can reach, where it would read as policy.
    if host not in LOOPBACK and notice_placeholders():
        sys.exit(
            f"Refusing to bind {host}: the upload notice still has undecided parts "
            f"({len(notice_placeholders())}), and uploads are kept. Fill them in "
            f"webapp/submissions.py first (tools/HOSTING.md job 7).")

    # timeout=0 disables trame's idle shutdown. Without it the server exits once no
    # browser is attached, so starting it and opening the tab a minute later finds
    # nothing listening -- which reads as a crash rather than as a timeout.
    #
    # Hosted, that same setting is a resource leak: an abandoned tab keeps its VTK
    # render window and mesh alive, measured at ~0.4 GB resident per session, and
    # concurrency is what this costs money by. So reclaim idle sessions there.
    idle = int(os.environ.get("VC_IDLE_TIMEOUT") or 300) if PUBLIC else 0
    server.start(port=port, open_browser=False, host=host, timeout=idle)
