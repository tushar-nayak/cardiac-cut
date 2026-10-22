"""Selecting 'Aortic arch' must now land on a bound that actually cuts."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import boundary_loops

app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
w = CutWindow(str(ROOT / "KY160_1_CentrelinePoints.csv"),
              str(ROOT / "Geometries & Meshes/mesh/KY160_1.mat"))
w.show()
for _ in range(15):
    app.processEvents()

fails = []


def actions():
    # re-fetch every time: _refresh_region_menu clears the menu and deletes the
    # old QAction objects, so a cached handle is a dangling pointer
    return {a.text().split("   ")[0]: a
            for a in w.region_btn.menu().actions() if a.text()}


print("presets offered:", list(actions()))

for name in ("Aortic arch", "Ishimaru zone 0", "Ishimaru zone 4"):
    a = actions().get(name)
    if a is None:
        print(f"\n{name}: not offered")
        continue
    a.trigger()
    for _ in range(25):
        app.processEvents()
    ok = w.problem is None
    print(f"\n{name}:")
    print(f"  region {w.region[0]:.1f} -> {w.region[1]:.1f} mm "
          f"({w.region[1]-w.region[0]:.0f} mm), {len(boundary_loops(w.cut))} openings")
    print(f"  {'VALID' if ok else 'INVALID'}   {w.status.text()[:120]}")
    if name == "Aortic arch" and not ok:
        fails.append("arch still invalid after the search")
    w._set_region(None)
    for _ in range(10):
        app.processEvents()

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
