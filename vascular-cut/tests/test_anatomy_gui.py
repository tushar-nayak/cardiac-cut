"""Anatomical naming through the GUI, on both geometries. argv[1] = jk5|ky160"""
import sys
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import pre_bifurcation_candidates

MODELS = {
    "jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
    "ky160": ("KY160_1_CentrelinePoints.csv", "Geometries & Meshes/mesh/KY160_1.mat"),
}
which = sys.argv[1]
csv, mat = MODELS[which]

app = QtWidgets.QApplication(sys.argv[:1])
# a modal warning from export() would block a headless run forever
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []
w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(10): app.processEvents()

print("=" * 80)
print(f"{which.upper()}, rows as the user sees them")
print("=" * 80)
for r in w.rows:
    print(f"    {r.terminal.name:34s} {r.terminal.segment.name:8s} {r.readout.text()}")
print(f"\n  anatomy note: {w.anat_lbl.text()[:96]}")
if not w.anat_lbl.text():
    fails.append("no-anat-label")
if w.problem is not None:
    fails.append(("invalid", w.problem[:50]))

print()
print("  trim menu:")
for a in w.trim_btn.menu().actions():
    print(f"    - {a.text()}")

# trim one and check the row takes the vessel's anatomical name
cands = pre_bifurcation_candidates(w.tree, w.base)
if cands:
    c = cands[0]
    w._add_trim(c)
    for _ in range(6): app.processEvents()
    tr = next(r for r in w.rows if getattr(r.terminal, "is_interior", False))
    print(f"\n  after trimming {c['segment']}: row is named '{tr.terminal.name}'")
    print(f"  status: {w.status.text()[:70]}")
    if w.problem is not None:
        fails.append(("trim-invalid", w.problem[:50]))
    vessel = w.anat_vessels.get(c["segment"])
    if vessel and tr.terminal.name != vessel:
        fails.append(("trim-name", tr.terminal.name))

# export metadata
w.out_dir = Path(__file__).resolve().parent / "_out"
w.name_edit.setText(f"anat_{which}")
w.export()
for _ in range(4): app.processEvents()
side = w.out_dir / f"anat_{which}_cut_planes.json"
if side.exists():
    meta = json.loads(side.read_text())
    print(f"\n  export anatomy block: {json.dumps(meta['anatomy'])[:120]}")
    named = [(o["name"], o["anatomical_name"]) for o in meta["openings"]]
    print("  export openings (name / inferred):")
    for n, a in named:
        print(f"    {n:34s} {a}")
    if "frame" not in meta["anatomy"]:
        fails.append("no-frame-in-export")
else:
    fails.append("no-export")

w.close()
print()
print("=" * 80)
print("ALL PASS" if not fails else f"FAILURES: {fails}")
print("=" * 80)

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
