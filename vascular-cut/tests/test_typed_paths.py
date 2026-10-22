"""Mesh and centerline paths can be typed, not only browsed.

The output folder has always been a text box and these two were not, so a cohort of
54 scans in sibling folders meant 54 trips through a file dialog to change one path
component. Requested by the lab in those words: "possibility of manually writing the
location into the input brackets, like at the output".

    python tests/test_typed_paths.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qtpy import QtWidgets                                          # noqa: E402

MAT = ROOT / "Geometries & Meshes/mesh/KY160_1.mat"
CSV = ROOT / "KY160_1_CentrelinePoints.csv"

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


app = QtWidgets.QApplication(sys.argv)
from branch_gui import CutWindow                                    # noqa: E402

w = CutWindow(None, None)

print("=" * 78)
print("the boxes accept typing at all")
print("=" * 78)
check(not w.mat_label.isReadOnly(), "the mesh path can be typed into")
check(not w.csv_label.isReadOnly(), "the centerline path can be typed into")

print()
print("=" * 78)
print("a path pasted from Explorer works as pasted")
print("=" * 78)
# "Copy as path" wraps the whole thing in double quotes, and a path dragged from a
# terminal often arrives with a trailing space. Neither exists on disk, so without
# trimming, a perfectly correct path is rejected as missing -- which is the one
# failure that sends someone back to the file dialog for good.
w.mat_label.setText(f'  "{MAT}" ')
w._mat_typed()
print(f"  typed '  \"{MAT.name}\" ' -> {w.mat_path}")
check(w.mat_path == MAT, "quotes and surrounding spaces are stripped")
check(w.mat_label.text() == str(MAT), "and the box is rewritten with the clean path")

w.csv_label.setText(f"'{CSV}'")
w._csv_typed()
check(w.csv_path == CSV, "single quotes too")

print()
print("=" * 78)
print("Load opens what was typed")
print("=" * 78)
w.load_model()
w._on_change()
print(f"  {len(w.rows)} opening(s), status: {w.status.text()[:60]}")
check(w.tree is not None, "the centerline loaded")
check(w.cut is not None, "and the mesh was cut")

print()
print("=" * 78)
print("a path that is not there says so, rather than failing at Load")
print("=" * 78)
w.mat_label.setText("Z:/nope/missing.mat")
w._mat_typed()
print(f"  {w.status.text()[:70]}")
check("does not exist" in w.status.text(), "a missing file is reported when typed")

print()
print("=" * 78)
print("an emptied centerline box means 'view the mesh only'")
print("=" * 78)
w.csv_label.setText("")
w._csv_typed()
check(w.csv_path is None, "clearing the box clears the centerline")
check("without cutting" in w.status.text(),
      f"and says what that means: {w.status.text()[:60]}")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
