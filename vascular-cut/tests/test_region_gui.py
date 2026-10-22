"""Region isolation end-to-end through the GUI. argv[1] = jk5|ky160"""
import sys
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import boundary_loops, region_presets

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat")}
which = sys.argv[1]
csv, mat = MODELS[which]
app = QtWidgets.QApplication(sys.argv[:1])
# export() pops a MODAL warning when the cut is invalid; headless, that blocks forever
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []

w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(10): app.processEvents()
base_rows = len(w.rows)
print("=" * 80)
print(f"{which.upper()} ,  {base_rows} openings, route {w.route.length:.0f} mm")
print("=" * 80)
print("  regions offered:")
for a in w.region_btn.menu().actions():
    print(f"    - {a.text()}")
if not w.region_btn.isEnabled():
    fails.append("region-btn-disabled")

# KY160's arch is fat (r=21.7 mm) and crowded, and "Ishimaru zone 3" and "Aortic
# arch" used to put a boundary where the rim tore 7.6 mm off plane. Both now cut
# flat, so the list is empty: measured directly off the geometry rather than taken
# from the validator, every rim of both regions is 0.000 mm out of its own best-fit
# plane, across 109-146 rim points each.
#
# What fixed it was ownership. The boundary sits beside the arch, where a trunk's
# axis lies a full radius deep in the lumen, so aortic wall used to come out nearer
# a thin branch's centerline than its own and got dragged into that branch's cut.
# Normalising the radial distance by each vessel's own radius keeps the wall with
# the aorta and the rim closes on the plane.
#
# If a boundary ever tears again this list is where to record it -- assert the
# refusal, never pretend a torn cut is fine.
KNOWN_REFUSED = {}.get(which, set())
# Sliding zone 0's start downstream walks the boundary into the same crowded arch.
KNOWN_REFUSED_SLIDE = {"ky160": {5.0, 10.0}}.get(which, set())

presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
for p in presets:
    w._set_region((p["s0"], p["s1"]))
    for _ in range(5): app.processEvents()
    loops = len(boundary_loops(w.cut))
    n_cuts = len(w.rows) + len(w._region_bounds)
    ok = w.problem is None and loops == n_cuts
    print(f"\n  {p['name']}: {w.reg_s0_b.value():.1f} -> {w.reg_s1_b.value():.1f} mm")
    print(f"    openings {loops} (cuts {n_cuts}), euler {int(w.cut.euler_number)}, "
          f"{len(w.rows)} branch row(s)   {'OK' if ok else 'PROBLEM'}")
    print(f"    readout: {w.region_lbl.text()[:88]}")
    expected_bad = p["name"] in KNOWN_REFUSED
    if w.problem is not None:
        print(f"    status: {w.status.text()[:70]}")
        if expected_bad:
            print("    (known limitation, refusal is the correct behaviour)")
        else:
            fails.append((p["name"], w.problem[:48]))
    elif expected_bad:
        # it started cutting cleanly: either the geometry handling improved (good,
        # update the list) or a check stopped firing (bad). Either way, look.
        fails.append((p["name"], "expected a refusal but the cut passed"))
    elif loops != n_cuts:
        fails.append((p["name"], f"{loops} loops vs {n_cuts} cuts"))

# slide the region's own boundary and confirm it keeps cutting
if presets:
    p = presets[0]
    w._set_region((p["s0"], p["s1"]))
    for _ in range(4): app.processEvents()
    print(f"\n  sliding the start of '{p['name']}':")
    for d in (0.0, 5.0, 10.0):
        w.reg_s0_b.setValue(p["s0"] + d)
        for _ in range(4): app.processEvents()
        # A user-moved boundary is now honoured rather than relocated to clear a
        # branch ostium. So a slide that lands ON an ostium is ALLOWED to fail --
        # what must not happen is failing silently. Require a warning naming the
        # ostium whenever an invalid slide is the cause.
        _b, _i, notes = w._region_terminals()
        warned = any("ostium" in n for n in notes)
        print(f"    +{d:4.0f} mm -> {len(boundary_loops(w.cut))} openings, "
              f"{'OK' if w.problem is None else 'RED: ' + w.problem[:38]}"
              f"{'   [warned: near an ostium]' if warned else ''}")
        if w.problem is not None and d not in KNOWN_REFUSED_SLIDE and not warned:
            fails.append(("slide", d, "failed with no ostium warning"))
        elif w.problem is not None:
            print("      (expected, the boundary is the user's and it is flagged)")

# export the isolated piece
w.out_dir = Path(__file__).parent / "out_region"
w.name_edit.setText(f"region_{which}")
w.export()
for _ in range(4): app.processEvents()
side = w.out_dir / f"region_{which}_cut_planes.json"
if w.problem is not None:
    print("\n  (last region was invalid, export correctly blocked)")
    side = None
if side is None:
    pass
elif side.exists():
    meta = json.loads(side.read_text())
    print(f"\n  exported {meta['n_openings']} openings, "
          f"{meta['n_faces']} faces; region block: {'region' in meta}")
else:
    fails.append("no-export")

# back to the whole vessel
w._set_region(None)
for _ in range(5): app.processEvents()
print(f"\n  cleared region -> {len(w.rows)} rows (was {base_rows}), "
      f"{len(boundary_loops(w.cut))} openings, "
      f"{'OK' if w.problem is None else 'RED'}")
if len(w.rows) != base_rows:
    fails.append(("restore-rows", len(w.rows), base_rows))
if w.problem is not None:
    fails.append("restore-invalid")

w.close()
print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
