"""Regressions for the five UI defects the adversarial review confirmed.

argv[1] = jk5|ky160
"""
import sys
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
# no modal dialog may block a headless run
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []
w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(15):
    app.processEvents()


def pump(n=8):
    for _ in range(n):
        app.processEvents()


print("=" * 78)
print(f"{which.upper()} ,  UI defect regressions")
print("=" * 78)

# --- 1. set_open(False) on an ALREADY folded section must recompute the summary ---
print("\n1. folded summary cannot go stale")
w.sec_cl.set_open(False)
pump()
before = w.sec_cl.header.text()
w.smooth_cb.blockSignals(True)
w.smooth_cb.setChecked(not w.smooth_cb.isChecked())   # change it from outside
w.smooth_cb.blockSignals(False)
w.sec_cl.set_open(False)                              # already folded -> was a no-op
pump()
after = w.sec_cl.header.text()
print(f"   before={before!r}")
print(f"   after ={after!r}")
ok = before != after
print(f"   summary re-derived while folded: {ok}")
if not ok:
    fails.append("stale-summary")
w.smooth_cb.blockSignals(True)
w.smooth_cb.setChecked(not w.smooth_cb.isChecked())   # put it back
w.smooth_cb.blockSignals(False)
w.sec_cl.refresh()
pump()

# --- 2. tuning the smoothing must not fold sections or clobber the export name ---
print("\n2. smoothing re-load keeps the user's context")
w.sec_cl.set_open(True)
w.name_edit.setText("my_export_name")
pump()
open_before = [s.header.isChecked() for s in
               (w.sec_input, w.sec_cl, w.sec_view, w.sec_out)]
w.smooth_mm.setValue(min(w.smooth_mm.value() + 0.5, 5.0))   # triggers a re-load
pump(12)
open_after = [s.header.isChecked() for s in
              (w.sec_input, w.sec_cl, w.sec_view, w.sec_out)]
name_kept = w.name_edit.text() == "my_export_name"
print(f"   section fold state {open_before} -> {open_after}")
print(f"   export name kept: {name_kept} ({w.name_edit.text()!r})")
if open_before != open_after:
    fails.append("smoothing-folded-sections")
if not name_kept:
    fails.append("smoothing-clobbered-name")

# --- 3. an empty region must be refused, not exported ---
print("\n3. collapsed region is refused")
L = w.route.length
w._set_region((0.0, 0.2), "degenerate")
pump(10)
n_active = len(w.cut_list()[0])
refused = w.problem is not None
print(f"   cut_list={n_active} terminals, problem={str(w.problem)[:56]!r}")
print(f"   refused: {refused}")
if n_active == 0 and not refused:
    fails.append("empty-region-accepted")
# the openings list must say WHY, not claim its two ends are cut
texts = [w.openings_lay.itemAt(i).widget().text()
         for i in range(w.openings_lay.count())
         if isinstance(w.openings_lay.itemAt(i).widget(), QtWidgets.QLabel)]
said = any("collapses to nothing" in t for t in texts)
print(f"   list explains the collapse: {said}")
if n_active == 0 and not said:
    fails.append("collapse-unexplained")
w._set_region(None)
pump()

# --- 4. region button must reset when the centerline goes away ---
print("\n4. region button resets on a mesh-only load")
presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
if presets:
    w._set_region((presets[0]["s0"], presets[0]["s1"]), presets[0]["name"])
    pump()
    print(f"   with region: {w.region_btn.text()!r} enabled={w.region_btn.isEnabled()}")
    w.csv_path = None
    w.csv_label.setText("")
    w.load_model()
    pump(12)
    txt, en = w.region_btn.text(), w.region_btn.isEnabled()
    print(f"   after mesh-only load: {txt!r} enabled={en}  region={w.region}")
    if txt != "Isolate a region…" or en:
        fails.append(("stale-region-btn", txt, en))
    # and Input must NOT be folded when the panel tells you to use Browse/Load
    print(f"   Input section open (panel says 'use Browse…'): "
          f"{w.sec_input.header.isChecked()}")
    if not w.sec_input.header.isChecked():
        fails.append("input-folded-mesh-only")

# --- 5. a failed load must leave Input reachable ---
print("\n5. failed load reopens Input")
w.sec_input.set_open(False)
pump()
w.mat_path = ROOT / "definitely_not_a_real_mesh.mat"
w.load_model()
pump(10)
print(f"   status: {w.status.text()[:60]!r}")
print(f"   Input section open: {w.sec_input.header.isChecked()}")
if not w.sec_input.header.isChecked():
    fails.append("failed-load-traps-user")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
sys.stdout.flush()

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
