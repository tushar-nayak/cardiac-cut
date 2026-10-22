"""Do keyboard / wheel / groove slider changes actually re-cut the mesh?

valueChanged drives the live preview; sliderReleased commits. But sliderReleased
fires ONLY after dragging the handle -- arrow keys, the scroll wheel and clicking
in the groove all change the value without it. If nothing commits, the plane moves
on screen while self.cut stays stale, and export writes a mesh that disagrees with
the sliders (and with the JSON metadata written beside it).

Programmatic slider.setValue() takes exactly the same path as a keypress: it emits
valueChanged and never sliderReleased.  argv[1] = jk5|ky160
"""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import region_presets

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


def pump(n=10):
    for _ in range(n):
        app.processEvents()


def settle(seconds=1.5):
    """Pump the event loop in REAL time so the debounced commit timer can fire."""
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        app.processEvents()
        time.sleep(0.01)


# Face count is a poor proxy: a region bound can be nudged back onto the same
# ostium clearance, giving an identical mesh from a genuinely re-run cut. Count the
# re-cuts directly, "did the change get committed" is the actual invariant.
_recuts = [0]
_real_recut = w.recut


def _counted_recut():
    _recuts[0] += 1
    return _real_recut()


w.recut = _counted_recut

bad = []
print("=" * 78)
print(f"{which.upper()} ,  does a non-drag slider change re-cut?")
print("=" * 78)

# ---- 1. a branch cut-depth slider ----
row = next((r for r in w.rows if not r.terminal.is_trunk), w.rows[0])
faces0 = len(w.cut.faces)
val0 = row.value()
print(f"\n1. cut depth ({row.terminal.name})")
print(f"   before: value={val0:.1f} mm, cut has {faces0} faces")
n0 = _recuts[0]
row.slider.setValue(row.slider.value() + 60)      # == pressing Up/PageUp
settle()                                          # let the debounced commit fire
val1 = row.value()
faces1 = len(w.cut.faces)
print(f"   after slider change: value={val1:.1f} mm, cut has {faces1} faces")
committed = _recuts[0] > n0
print(f"   re-cuts triggered: {_recuts[0] - n0}   (faces changed: {faces1 != faces0})")
if abs(val1 - val0) < 1e-9:
    print("   (slider did not move, inconclusive)")
elif not committed:
    bad.append("cut-depth-not-committed")
    print("   !! the slider says one depth, the mesh is cut at another")

# ---- 2. region bound sliders ----
presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
if presets:
    p = presets[0]
    w._set_region((p["s0"], p["s1"]), p["name"])
    pump(12)
    faces0 = len(w.cut.faces)
    prob0 = w.problem
    s_before = w.reg_s1_b.value()
    print(f"\n2. region end bound ('{p['name']}')")
    print(f"   before: end={s_before:.1f} mm, {faces0} faces, problem={prob0 is not None}")
    n0 = _recuts[0]
    w.reg_s1_s.setValue(max(w.reg_s1_s.minimum(),
                            w.reg_s1_s.value() - 300))   # 30 mm shorter, via keyboard
    settle()                                        # let the debounced commit fire
    s_after = w.reg_s1_b.value()
    faces2 = len(w.cut.faces)
    print(f"   after slider change: end={s_after:.1f} mm, {faces2} faces, "
          f"problem={w.problem is not None}")
    print(f"   self.region now {tuple(round(x, 1) for x in w.region)}")
    committed = _recuts[0] > n0
    print(f"   re-cuts triggered: {_recuts[0] - n0}   (faces changed: {faces2 != faces0})")
    if abs(s_after - s_before) < 1e-9:
        print("   (slider did not move, inconclusive)")
    elif not committed:
        bad.append("region-bound-not-committed")
        print("   !! self.region moved but the mesh was never re-cut")
        # and would it export that mismatch?
        print(f"   status bar still reads: {w.status.text()[:64]!r}")
    w._set_region(None)
    pump()

print()
print("ALL PASS" if not bad else f"BUG CONFIRMED: {bad}")
sys.stdout.flush()
