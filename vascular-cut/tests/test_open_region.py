"""Cut above the renals, keep the iliacs, then trim them as a group.

The workflow that could not be expressed before: an open-ended region is one plane
with everything below it kept, so the iliacs stay as ordinary openings and become
trimmable.
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import boundary_loops, region_presets

app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []

w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
for _ in range(10):
    app.processEvents()

presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
print("=" * 76)
print("presets offered")
print("=" * 76)
for p in presets:
    span = "open" if p["s1"] is None else f"{p['s1']:.1f}"
    print(f"  {p['name']:<32} {p['s0']:6.1f} -> {span}")

target = next((p for p in presets if p["name"] == "Infrarenal aorta with iliacs"),
              None)
if target is None:
    print("\nFAILURES: the open-ended preset was not offered")
    raise SystemExit(1)

w._set_region((target["s0"], target["s1"]), target["name"])
for _ in range(8):
    app.processEvents()

bounds, inside, _ = w._region_terminals()
print()
print("=" * 76)
print(f"{target['name']}, s0 = {target['s0']:.1f} mm, open distal end")
print("=" * 76)
print(f"  status : {w.status.text()[:70]}")
print(f"  bounds : {len(bounds)}  (one plane, not two)")
print(f"  inside : {[t.name for t in inside]}")
print(f"  rows   : {len(w.rows)}   openings in the mesh: "
      f"{len(boundary_loops(w.cut)) if w.cut is not None else -1}")
if len(bounds) != 1:
    fails.append(f"expected a single bound, got {len(bounds)}")
iliacs = [t for t in inside if "iliac" in t.name.lower()]
print(f"  iliacs kept: {[t.name for t in iliacs]}")
if len(iliacs) < 2:
    fails.append(f"the iliacs were not kept ({len(iliacs)} found)")
if w.problem is not None:
    fails.append(f"the region did not cut: {w.problem[:60]}")

print()
print("=" * 76)
print("now trim the iliacs as a group, inside that region")
print("=" * 76)
acts = list(w.trim_btn.menu().actions())
print(f"  trims offered here: {len(acts)}")
for a in acts:
    print(f"    {a.text()[:66]}")
if not acts:
    fails.append("no trim offered inside the open region")
else:
    before = len(w.rows)
    acts[0].trigger()
    for _ in range(8):
        app.processEvents()
    interior = [t for t in w.active_terminals() if getattr(t, "is_interior", False)]
    print(f"\n  after the trim: rows {before} -> {len(w.rows)}, "
          f"recorded={len(w.trims)}, interior cuts={len(interior)}")
    print(f"  openings: {len(boundary_loops(w.cut)) if w.cut is not None else -1}")
    print(f"  status: {w.status.text()[:70]}")
    if len(w.rows) >= before:
        fails.append("the trim did not collapse any openings")

    acts2 = list(w.trim_btn.menu().actions())
    print(f"\n  still offered: {len(acts2)}")
    if acts2:
        acts2[0].trigger()
        for _ in range(8):
            app.processEvents()
        print(f"  after a SECOND trim: rows={len(w.rows)}, "
              f"recorded={len(w.trims)}, "
              f"openings={len(boundary_loops(w.cut)) if w.cut is not None else -1}")
        print(f"  status: {w.status.text()[:70]}")
        if len(w.trims) < 2:
            fails.append("a second trim inside the region did not stick")

# and back out cleanly
w._set_region(None)
for _ in range(6):
    app.processEvents()
print(f"\n  cleared: {len(w.rows)} rows, {w.status.text()[:52]}")

w.close()
print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
