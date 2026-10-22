"""Several trims at once on the whole vessel, which is where they can apply today.

JK5_1 offers four candidates: each iliac pair, the celiac group, and the aorta
above the bifurcation. Applying two independent ones must leave two interior cuts,
not one.
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import boundary_loops

app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []

w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
for _ in range(10):
    app.processEvents()

print("=" * 74)
print("whole vessel, trims applied one after another")
print("=" * 74)
print(f"  start: {len(w.rows)} rows, {len(boundary_loops(w.cut))} openings")
print(f"         {w.status.text()[:66]}")

applied = 0
for step in range(3):
    acts = [a for a in w.trim_btn.menu().actions()]
    if not acts:
        print(f"  no more trims offered after {applied}")
        break
    label = acts[0].text()[:46]
    acts[0].trigger()
    for _ in range(6):
        app.processEvents()
    applied += 1
    interior = [t for t in w.active_terminals() if getattr(t, "is_interior", False)]
    print(f"\n  trim {applied}: {label}")
    print(f"    rows={len(w.rows)}  recorded={len(w.trims)}  "
          f"interior cuts={len(interior)}  openings={len(boundary_loops(w.cut))}")
    print(f"    {w.status.text()[:66]}")
    if len(interior) != len(w.trims):
        fails.append(f"after {applied} trims: {len(w.trims)} recorded but "
                     f"{len(interior)} interior cuts in the active list")

if applied >= 2 and len(w.trims) >= 2:
    print(f"\n  {len(w.trims)} trims coexist, each producing its own interior cut")
elif applied >= 2:
    fails.append("a second trim did not stick")

# undo one and check the swallowed openings come back
if w.trims:
    n_before = len(w.rows)
    w._remove_trim(w.trims[0]["terminal"])
    for _ in range(6):
        app.processEvents()
    print(f"\n  removed one trim: rows {n_before} -> {len(w.rows)}, "
          f"recorded={len(w.trims)}")
    if len(w.rows) <= n_before:
        fails.append("removing a trim did not restore its openings")

w.close()
print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
