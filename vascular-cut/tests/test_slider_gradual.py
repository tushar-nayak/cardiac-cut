"""The region sliders must move the boundary GRADUALLY, not snap to a few spots."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numpy as np
from qtpy import QtWidgets
from branch_gui import CutWindow

app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
w = CutWindow(str(ROOT / "KY160_1_CentrelinePoints.csv"),
              str(ROOT / "Geometries & Meshes/mesh/KY160_1.mat"))
w.show()
for _ in range(15):
    app.processEvents()

acts = {a.text().split("   ")[0]: a for a in w.region_btn.menu().actions() if a.text()}
acts["Aortic arch"].trigger()
for _ in range(25):
    app.processEvents()
print(f"arch preset -> {w.region[0]:.1f} .. {w.region[1]:.1f} mm, "
      f"{'valid' if w.problem is None else 'INVALID'}")

# now drag the start across the range that used to be dead (56 -> 118 mm)
print(f"\n{'asked':>8}{'actual bound':>14}   note")
seen = []
for ask in np.arange(56.0, 120.0, 4.0):
    w.reg_s0_b.setValue(float(ask))
    for _ in range(6):
        app.processEvents()
    bounds, _inside, notes = w._region_terminals()
    got = float(w.region[0])
    seen.append(round(got, 1))
    warn = next((n for n in notes if "sits" in n), "")
    print(f"{ask:8.1f}{got:14.1f}   {warn[:64]}")

distinct = len(set(seen))
print(f"\n  {len(seen)} slider positions produced {distinct} distinct boundaries")
ok = distinct >= 0.9 * len(seen)
print(f"  slider tracks the boundary: {ok}")
print()
print("ALL PASS" if ok else "FAIL: the slider still snaps")
