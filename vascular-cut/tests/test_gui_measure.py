"""The size-and-angulation readout in the desktop window.

The maths already shipped inside branch_cut, but only the web layer could reach it,
so a collaborator with the zip had the code and no way to run it. Drive the button.
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

w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
for _ in range(10):
    app.processEvents()


def strip(html):
    import re
    return re.sub(r"<[^>]+>", "", html).replace("&empty;", "D").replace("&deg;", "deg")


print("=" * 76)
print("whole vessel")
print("=" * 76)
w.size_btn.click()
for _ in range(6):
    app.processEvents()
print("  " + strip(w.size_lbl.text()).replace("\n", "\n  "))
if not w.size_lbl.text():
    fails.append("no readout on the whole vessel")

for name in ["Aortic neck", "Infrarenal aorta with iliacs"]:
    p = next((x for x in region_presets(w.tree, w.base, w.anat_names,
                                        w.anat_vessels) if x["name"] == name), None)
    if p is None:
        continue
    w._set_region((p["s0"], p["s1"]), name)
    for _ in range(6):
        app.processEvents()
    w.size_btn.click()
    for _ in range(6):
        app.processEvents()
    print()
    print("=" * 76)
    print(name)
    print("=" * 76)
    txt = strip(w.size_lbl.text())
    print("  " + txt.replace("\n", "\n  "))
    if "angle" not in txt.lower() and "straight" not in txt.lower():
        fails.append(f"{name}: no angle reported")
    if "max D" not in txt:
        fails.append(f"{name}: no max diameter reported")

print()
print("=" * 76)
print("the button is safe before anything is loaded")
print("=" * 76)
w2 = CutWindow(None, None)
w2.size_btn.click()
for _ in range(4):
    app.processEvents()
print(f"  status: {w2.status.text()[:60]}")
if w2.size_lbl.text():
    fails.append("produced a readout with nothing loaded")
w2.close()

w.close()
print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
