"""Can a pending debounce timer resurrect a region the user just dismissed?

The debounced commit fires 220 ms after a non-drag slider change. If the user
nudges a region bound with the wheel and then immediately clicks "Whole vessel",
the timer fires afterwards and _region_commit re-reads the spin boxes -- putting
self.region back while region_bar is hidden. That would recreate the exact trap
the region-escape work set out to remove, and leave it HARDER to escape.
argv[1] = jk5|ky160
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


def settle(seconds=1.5):
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        app.processEvents()
        time.sleep(0.01)


bad = []
base_rows = len(w.rows)
presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
print("=" * 78)
print(f"{which.upper()} ,  can a pending timer resurrect a dismissed region?")
print("=" * 78)

if presets:
    p = presets[0]
    for route_name, escape in (("banner / corner button", lambda: w.region_clear.click()),
                               ("Reset", lambda: w.reset_btn.click())):
        w._set_region((p["s0"], p["s1"]), p["name"])
        settle(0.6)
        print(f"\n  region on ({p['name']}), leaving via {route_name}")
        # wheel/arrow nudge, then dismiss the region INSIDE the debounce window
        w.reg_s1_s.setValue(max(w.reg_s1_s.minimum(), w.reg_s1_s.value() - 100))
        escape()                     # user leaves immediately, timer still pending
        for _ in range(6):
            app.processEvents()
        print(f"    right after escaping: region={w.region}, "
              f"bar visible={w.region_bar.isVisible()}, rows={len(w.rows)}")
        settle(1.5)                  # let any pending timer fire
        print(f"    after the timer window: region={w.region}, "
              f"bar visible={w.region_bar.isVisible()}, rows={len(w.rows)}, "
              f"btn={w.region_btn.text()!r}")
        if w.region is not None:
            bad.append((route_name, "region resurrected"))
            print("    !! the region came back after the user dismissed it")
            if not w.region_bar.isVisible():
                print("    !! and the region bar is HIDDEN, no sliders, worse than before")
        if len(w.rows) != base_rows:
            bad.append((route_name, f"rows {len(w.rows)} != {base_rows}"))
        w._set_region(None)
        settle(0.5)

# a row's own timer must not survive its row being rebuilt by a trim/region change
print("\n  row commit timer vs a rebuild")
row = next((r for r in w.rows if not r.terminal.is_trunk), w.rows[0])
row.slider.setValue(row.slider.value() + 40)
w._build_rows()                       # rebuild while the timer is pending
for _ in range(6):
    app.processEvents()
settle(1.2)
print(f"    survived rebuild without error, rows={len(w.rows)}")

print()
print("ALL PASS" if not bad else f"BUG CONFIRMED: {bad}")
sys.stdout.flush()
