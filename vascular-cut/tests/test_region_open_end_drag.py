"""Dragging the Start plane must not close an open-ended region.

Reported from the abdominal cohort, by a user who cut 48 of 49 meshes:

    "'Infrarenal aorta with iliacs' function sometimes cuts quite low, and if you
     try to move the proximal cutting plane it messes it up and cuts off the iliacs
     (distal cutting plane moves up and you cannot change it back)"

"Infrarenal aorta with iliacs" is the only region that ships ``s1 = None`` -- one
plane, everything below it kept -- and it does so deliberately, because region_plan
clamps a distal plane a millimetre short of the route end while an abdominal aorta's
iliac ostia sit AT the bifurcation. A two-plane version drops them.

_sync_region_bar DISABLES the End control for such a region, so the number inside it
is never the user's. It is whatever was left there: 0.0 from the widget's own
construction in a fresh session, or the End of whichever two-plane region ran before.
_region_dragged and _region_commit read it back unconditionally, which replaced None
with that stale number on the first touch of the Start slider.

Both outcomes match the report. With a stale 0.0 the span INVERTS -- _region_terminals
does min/max -- so the region becomes the aorta ABOVE the drag point and the distal
plane appears at the top, which is the "distal cutting plane moves up". With a stale
route length it clamps short of the bifurcation and the iliacs fall outside. Neither
could be undone, because the End control driving that plane was still greyed out and
still carrying the tooltip saying the region had no end.

WHAT THIS TEST WOULD SHOW WITHOUT THE FIX
Revert _span_from_controls in branch_gui.py (make _region_dragged and _region_commit
read ``(reg_s0_b.value(), reg_s1_b.value())`` again) and the checks below fail at
"open end survived the drag", with two bounds instead of one and the iliacs gone.

    .venv\\Scripts\\python tests\\test_region_open_end_drag.py
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import region_presets

app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []


def pump(n=8):
    for _ in range(n):
        app.processEvents()


w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
pump(10)

presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
target = next((p for p in presets if p["name"] == "Infrarenal aorta with iliacs"),
              None)
if target is None:
    print("FAILURES: the open-ended preset was not offered on JK5_1")
    raise SystemExit(1)

w._set_region((target["s0"], target["s1"]), target["name"])
pump()

# The staleness itself, before anything is dragged. This is the value the old code
# read back. Printed rather than asserted: it is the CAUSE, and pinning its exact
# number would make this test brittle for no gain.
stale = float(w.reg_s1_b.value())
bounds0, inside0, _ = w._region_terminals()
iliacs0 = [t for t in inside0 if "iliac" in t.name.lower()]

print("=" * 76)
print(f"{target['name']}: s0 = {target['s0']:.1f} mm, s1 = {w.region[1]}")
print("=" * 76)
print(f"  End control enabled : {w.reg_s1_b.isEnabled()}  (must be False)")
print(f"  End control value   : {stale:.1f} mm  <- stale, not the user's")
print(f"  bounds              : {len(bounds0)}   iliacs inside: {len(iliacs0)}")
if w.region[1] is not None:
    fails.append("the preset did not arrive with an open distal end")
if w.reg_s1_b.isEnabled():
    fails.append("the End control is editable on a region that has no end")
if len(iliacs0) < 2:
    fails.append(f"the iliacs were not inside the region to begin with "
                 f"({len(iliacs0)})")

# ---- the report: move the PROXIMAL plane -----------------------------------
# Down the aorta by 5 mm, the way a user nudges a boundary clear of something.
moved_to = float(target["s0"]) + 5.0
print()
print("=" * 76)
print(f"drag Start {target['s0']:.1f} -> {moved_to:.1f} mm")
print("=" * 76)
w._region_dragged(w.reg_s0_b, int(round(moved_to * 10)), slider=w.reg_s0_s)
pump()

bounds1, inside1, _ = w._region_terminals()
iliacs1 = [t for t in inside1 if "iliac" in t.name.lower()]
print(f"  region              : {w.region}")
print(f"  bounds              : {len(bounds1)}   (one plane, still)")
print(f"  iliacs inside       : {[t.name for t in iliacs1]}")
print(f"  End control enabled : {w.reg_s1_b.isEnabled()}")

if w.region[1] is not None:
    fails.append(f"open end survived the drag: s1 became {w.region[1]} "
                 f"(stale End box held {stale:.1f})")
if len(bounds1) != 1:
    fails.append(f"the drag added a distal plane: {len(bounds1)} bounds, expected 1")
if len(iliacs1) < 2:
    fails.append(f"the drag cut the iliacs out of the region "
                 f"({len(iliacs0)} -> {len(iliacs1)})")
# The spin box carries one decimal, so a station read back through it is quantised
# to 0.1 mm. Half a step is the tightest this can honestly be.
QUANTUM = 0.05 + 1e-9
if abs(float(w.region[0]) - moved_to) > QUANTUM:
    fails.append(f"the Start plane did not move: {w.region[0]} != {moved_to}")

# ---- and the committed path, which is a second copy of the same read --------
print()
print("=" * 76)
print("commit (the debounced path, and what a typed value takes)")
print("=" * 76)
w._region_commit()
pump()
bounds2, inside2, _ = w._region_terminals()
iliacs2 = [t for t in inside2 if "iliac" in t.name.lower()]
print(f"  region        : {w.region}")
print(f"  bounds        : {len(bounds2)}   iliacs inside: {len(iliacs2)}")
if w.region[1] is not None:
    fails.append(f"commit closed the open end: s1 became {w.region[1]}")
if len(iliacs2) < 2:
    fails.append(f"commit cut the iliacs out ({len(iliacs2)})")

# ---- a two-plane region must still read BOTH boxes -------------------------
# The fix keeps an open end open; it must not freeze the End of a normal region.
print()
print("=" * 76)
print("a two-plane region still tracks its End control")
print("=" * 76)
plain = next((p for p in presets if p["name"] == "Infrarenal aorta"), None)
if plain is None:
    fails.append("'Infrarenal aorta' was not offered, so the guard is untested")
else:
    w._set_region((plain["s0"], plain["s1"]), plain["name"])
    pump()
    end_to = float(plain["s1"]) - 8.0
    w._region_dragged(w.reg_s1_b, int(round(end_to * 10)), slider=w.reg_s1_s)
    pump()
    print(f"  region  : {w.region}")
    print(f"  End enabled: {w.reg_s1_b.isEnabled()}")
    if w.region[1] is None:
        fails.append("a two-plane region lost its End")
    elif abs(float(w.region[1]) - end_to) > QUANTUM:
        fails.append(f"the End plane did not move: {w.region[1]} != {end_to}")

w.close()
print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
