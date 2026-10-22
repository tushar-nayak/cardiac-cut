"""Every way a stuck user might try to leave a region must actually work.

One model per process: two QtInteractor windows in one process segfaults VTK.
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
fails = []
w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(15):
    app.processEvents()
presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
base_rows = len(w.rows)
print("=" * 78)
print(f"{which.upper()} ,  {base_rows} openings, {len(presets)} preset region(s)")
print("=" * 78)


def enter(p):
    w._set_region((p["s0"], p["s1"]), p["name"])
    for _ in range(8):
        app.processEvents()


def escaped(route):
    for _ in range(8):
        app.processEvents()
    ok = (w.region is None and len(w.rows) == base_rows
          and not w.region_bar.isVisible()
          and w.region_btn.text() == "Isolate a region…")
    print(f"    {route:34s} -> region={w.region}, {len(w.rows)} rows, "
          f"btn={w.region_btn.text()!r}   {'OK' if ok else 'STUCK'}")
    if not ok:
        fails.append((which, route))


if presets:
    p = presets[0]
    print(f"\n  entering '{p['name']}'")

    # route 1: banner button at the top of the openings list
    enter(p)
    banner = w.openings_lay.itemAt(0).widget()
    btn = banner.findChild(QtWidgets.QPushButton) if banner else None
    btn_txt = repr(btn.text()) if btn else None
    print(f"  banner: {banner.__class__.__name__}, button={btn_txt}")
    if btn is None:
        fails.append((which, "no-banner"))
    else:
        btn.click()
        escaped("banner exit")

    # route 2: the region menu itself
    enter(p)
    print(f"  button renamed while active: {w.region_btn.text()!r}")
    acts = [a for a in w.region_btn.menu().actions()
            if "whole vessel" in a.text().lower() and a.text().startswith("←")]
    print(f"  menu exit entries: {[a.text() for a in acts]}")
    if not acts:
        fails.append((which, "no-menu-exit"))
    else:
        acts[0].trigger()
        escaped("menu exit")

    # route 3: Reset
    enter(p)
    w.reset_btn.click()
    escaped("Reset button")

    # route 4: the original corner button
    enter(p)
    w.region_clear.click()
    escaped("original corner button")

    # a region with nothing inside must explain itself, not just go blank
    for pp in presets:
        enter(pp)
        if w.rows:
            continue
        texts = [w.openings_lay.itemAt(i).widget().text()
                 for i in range(w.openings_lay.count())
                 if isinstance(w.openings_lay.itemAt(i).widget(), QtWidgets.QLabel)]
        said = any("No branches leave" in t for t in texts)
        print(f"\n  empty region '{pp['name']}' explains itself: {said}")
        if not said:
            fails.append((which, f"empty-silent:{pp['name']}"))
        break
    w._set_region(None)

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
sys.stdout.flush()

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
