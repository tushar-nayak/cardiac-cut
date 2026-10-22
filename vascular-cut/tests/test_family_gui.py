"""The Anatomy control: declaring the family must re-read the case end to end.

argv[1] = jk5|ky160
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]

app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(15):
    app.processEvents()


def pump(n=12):
    for _ in range(n):
        app.processEvents()


def state():
    n_trunk = sum(1 for t in w.tree.terminals if getattr(t, "is_trunk", False))
    return {"rows": len(w.rows), "trunk": n_trunk,
            "named": len(w.anat_names),
            "presets": tuple(a.text().split("   ")[0]
                             for a in w.region_btn.menu().actions()
                             if not a.text().startswith("←") and a.text()),
            "problem": w.problem is not None}


fails = []
print("=" * 78)
print(f"{which.upper()} ,  Anatomy declaration")
print("=" * 78)
print(f"  control present: {w.family.count()} options "
      f"{[w.family.itemText(i) for i in range(w.family.count())]}")
if w.family.count() != 3:
    fails.append("combo missing options")

auto = state()
print(f"\n  Auto-detect : {auto['trunk']} trunk end(s), {auto['rows']} rows, "
      f"{auto['named']} named")
print(f"     presets: {auto['presets']}")
print(f"     label  : {w.anat_lbl.text()[:96]}")

for idx, label in ((1, "Abdominal"), (2, "Thoracic")):
    w.family.setCurrentIndex(idx)
    pump(25)
    s = state()
    print(f"\n  {label:<11}: {s['trunk']} trunk end(s), {s['rows']} rows, "
          f"{s['named']} named, cut {'INVALID' if s['problem'] else 'valid'}")
    print(f"     presets: {s['presets']}")
    print(f"     label  : {w.anat_lbl.text()[:96]}")
    want = 1 if idx == 1 else 2
    if s["trunk"] != want:
        fails.append(f"{label}: {s['trunk']} trunk ends, expected {want}")
    if "you set" not in w.anat_lbl.text() and s["named"]:
        fails.append(f"{label}: label does not say the anatomy was declared")

# back to auto must restore exactly
w.family.setCurrentIndex(0)
pump(25)
back = state()
print(f"\n  back to Auto: {back['trunk']} trunk end(s), {back['rows']} rows, "
      f"{back['named']} named")
if back != auto:
    for k in auto:
        if auto[k] != back[k]:
            fails.append(f"auto not restored: {k} {auto[k]} -> {back[k]}")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
