"""Exporting the same case twice must not raise a dialog and must report clearly."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow

app = QtWidgets.QApplication(sys.argv)
# a modal warning from export() would block a headless run forever
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []

# Any modal opened during export is a failure: trip a flag if one appears.
opened = []
for nm in ("question", "warning", "information", "critical"):
    orig = getattr(QtWidgets.QMessageBox, nm)
    def spy(*a, _n=nm, **k):
        opened.append(_n)
        return QtWidgets.QMessageBox.No
    setattr(QtWidgets.QMessageBox, nm, staticmethod(spy))

w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
for _ in range(10): app.processEvents()

out = Path(__file__).parent / "out_overwrite"
w.out_dir = out
w.name_edit.setText("repeat_case")
stem = out / f"repeat_case.{w.fmt.currentText()}"

for i in (1, 2, 3):
    w.export()
    for _ in range(3): app.processEvents()
    print(f"  export #{i}: dialogs={opened or 'none'}   file={'yes' if stem.exists() else 'NO'}")
    print(f"    status: {w.status.text()[:78]}")
    if opened:
        fails.append(("dialog-on-export", i, list(opened)))
        opened.clear()

if not stem.exists():
    fails.append("no-file-written")
# second and third runs should say Replaced, the first Exported
if "Replaced" not in w.status.text():
    fails.append(("no-replaced-wording", w.status.text()[:60]))

w.close()
print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
