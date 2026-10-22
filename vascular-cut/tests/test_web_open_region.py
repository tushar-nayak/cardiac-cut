"""The web tool must handle a region with no distal end, as the desktop does.

Found by a desktop/web parity audit after the same bug was fixed in the desktop
window (branch_gui._span_from_controls). The web was worse in two ways:

  1. it could not select "Infrarenal aorta with iliacs" AT ALL. That preset ships
     s1 = None, and on_region did `round(float(p["s1"]), 1)`. float(None) raises,
     the raise happened BEFORE _recut_and_draw, and the exception was swallowed by a
     blanket handler into a raw "TypeError: float() argument..." in the status line.
     state.in_region had already been set, so the "Isolated region" banner and the
     Start/End sliders appeared over a scene still showing the previous whole-vessel
     cut. The preset is in the web's own menu -- nothing filters it out.

  2. commit_region read the End box unconditionally and nothing ever disabled that
     box, so it carried both halves of the defect the desktop only had one of.

vascular_api is CLEAN and always was: set_region and apply_preset both honour
s1 None. The bug was entirely in the web layer on top of it.

This drives the real ctrl handlers against a real CutSession, with no browser: the
defect is in Python state, not in the DOM, and tests/web/* need playwright.

    .venv\\Scripts\\python tests\\test_web_open_region.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
sys.path.insert(0, str(ROOT))

import cut_app as APP                                              # noqa: E402
# on_region is a @state.change handler, not a ctrl method, so it is called
# directly off the module; commit_region IS a ctrl and is called through it.
from vascular_api import CutSession                                # noqa: E402

fails = []
state, ctrl = APP.state, APP.ctrl

s = CutSession(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
               str(ROOT / "JK5_1_M10.mat"))
APP.SESSION["s"] = s
state.loaded = True
state.route_length = float(s.route.length)

presets = s.presets()
OPEN = "Infrarenal aorta with iliacs"
TWO_PLANE = "Infrarenal aorta"
print("=" * 76)
print("presets offered to a web user")
print("=" * 76)
for p in presets:
    print(f"  {p['name']:<34} s0={p['s0']:7.2f}  s1={p['s1']}")
if not any(p["name"] == OPEN for p in presets):
    print(f"\nFAILURES: {OPEN} not offered on JK5_1")
    raise SystemExit(1)


def iliacs_inside():
    """Names of iliac openings the current region keeps."""
    try:
        return [o["name"] for o in s.openings()
                if "iliac" in o["name"].lower() and not o.get("is_region_bound")]
    except Exception as e:                                          # noqa: BLE001
        return [f"<openings raised {type(e).__name__}: {e}>"]


# ---- 1. selecting the open preset must not raise ---------------------------
print()
print("=" * 76)
print(f"1. select '{OPEN}' (the one that used to raise TypeError)")
print("=" * 76)
state.status_bad = False
state.region_pick = OPEN
APP.on_region(OPEN)
print(f"  status      : {str(state.status)[:80]}")
print(f"  in_region   : {state.in_region}   region_open: {state.region_open}")
print(f"  session     : {s.region}")
print(f"  iliacs kept : {iliacs_inside()}")
if state.status_bad:
    fails.append(f"selecting the open preset reported an error: {state.status}")
if "TypeError" in str(state.status):
    fails.append("the float(None) TypeError is still raised")
if s.region is None or s.region[1] is not None:
    fails.append(f"the session did not take an open end: {s.region}")
if not state.region_open:
    fails.append("region_open was not set, so the End box stays live")
opening_iliacs = iliacs_inside()
if len(opening_iliacs) < 2:
    fails.append(f"the iliacs were not kept: {opening_iliacs}")

# ---- 2. moving Start must not close it -------------------------------------
print()
print("=" * 76)
print("2. drag Start +5 mm and commit -- the open end must survive")
print("=" * 76)
moved_to = round(float(state.region_s0) + 5.0, 1)
state.region_s0 = moved_to
ctrl.commit_region()
print(f"  session     : {s.region}")
print(f"  region_open : {state.region_open}")
print(f"  iliacs kept : {iliacs_inside()}")
if s.region[1] is not None:
    fails.append(f"the commit closed the open end: s1 became {s.region[1]} "
                 f"(End box held {state.region_s1})")
if abs(float(s.region[0]) - moved_to) > 0.06:
    fails.append(f"the Start plane did not move: {s.region[0]} != {moved_to}")
if len(iliacs_inside()) < 2:
    fails.append(f"the commit cut the iliacs out: {iliacs_inside()}")

# ---- 3. the carry-over case: a stale End from a two-plane region -----------
# The nastier of the two failure modes. Pick a two-plane region first so the End box
# holds a real number, then pick the open one: the stale number must not be used.
print()
print("=" * 76)
print(f"3. '{TWO_PLANE}' first (fills the End box), then the open one")
print("=" * 76)
if any(p["name"] == TWO_PLANE for p in presets):
    state.region_pick = TWO_PLANE
    APP.on_region(TWO_PLANE)
    stale = float(state.region_s1)
    print(f"  after '{TWO_PLANE}': region={s.region}  End box={stale:.1f}")
    if state.region_open:
        fails.append("a two-plane region was marked open")
    state.region_pick = OPEN
    APP.on_region(OPEN)
    state.region_s0 = round(float(state.region_s0) + 3.0, 1)
    ctrl.commit_region()
    print(f"  then '{OPEN}' + drag: region={s.region}  (stale End was {stale:.1f})")
    print(f"  iliacs kept : {iliacs_inside()}")
    if s.region[1] is not None:
        fails.append(f"the stale End ({stale:.1f}) closed the open region: "
                     f"{s.region[1]}")
    if len(iliacs_inside()) < 2:
        fails.append("the carry-over case cut the iliacs out")
else:
    fails.append(f"'{TWO_PLANE}' not offered, so the carry-over case is untested")

# ---- 4. and a two-plane region must still track BOTH boxes ----------------
print()
print("=" * 76)
print("4. a two-plane region still moves its End")
print("=" * 76)
state.region_pick = TWO_PLANE
APP.on_region(TWO_PLANE)
end_to = round(float(state.region_s1) - 8.0, 1)
state.region_s1 = end_to
ctrl.commit_region()
print(f"  region      : {s.region}   region_open: {state.region_open}")
if s.region[1] is None:
    fails.append("a two-plane region lost its End")
elif abs(float(s.region[1]) - end_to) > 0.06:
    fails.append(f"the End did not move: {s.region[1]} != {end_to}")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
