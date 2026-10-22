"""A region's own two boundary cuts must be adjustable. argv[1] = jk5|ky160

Reported from the lab: "Is there a way to manually adjust the cutting plane size
and tilt once you enable the 'isolate region' function, I can't find it?" There
was not -- an isolated region gave rows only to the branches INSIDE it, and its two
boundary planes, the ones actually doing the isolating, had no controls at all.

The trap this guards is not that the controls exist but that they LAST.
region_plan builds fresh Terminal objects on every call and is called again on
every nudge of a region slider, so a tilt set by hand is discarded by the next drag
unless it is carried across -- silently, with the plane snapping back while the
user is looking at the number they typed.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets                                          # noqa: E402
from branch_gui import CutWindow                                    # noqa: E402
from branch_cut import band_width                                   # noqa: E402

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]
app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []

w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(10):
    app.processEvents()

print("=" * 78)
print(f"{which.upper()} ,  {len(w.rows)} openings, route {w.route.length:.0f} mm")
print("=" * 78)

# ---- no region: no bound rows, and nothing else changed ---------------------
if w.bound_rows:
    fails.append(f"bound rows before any region: {len(w.bound_rows)}")
print(f"  whole vessel : {len(w.rows)} rows, {len(w.bound_rows)} bound row(s)")

# ---- enter a two-plane region ----------------------------------------------
preset = next((p for p in w.region_btn.menu().actions()
               if "Back to" not in p.text()), None)
if preset is None:
    print("no region presets offered")
    raise SystemExit(1)
preset.trigger()
for _ in range(10):
    app.processEvents()

print(f"  region '{w._region_name}' : {len(w.rows)} branch row(s), "
      f"{len(w.bound_rows)} bound row(s), {len(w._region_bounds)} bound(s)")
if len(w.bound_rows) != len(w._region_bounds):
    fails.append(f"{len(w._region_bounds)} bounds but {len(w.bound_rows)} rows")
if not w.bound_rows:
    fails.append("region has no adjustable boundary planes")

for r in w.bound_rows:
    # The controls Matija could not find, and the offset slider that would be a
    # second way to move a plane the region sliders already own.
    if not r.adv.isVisible():
        fails.append(f"{r.terminal.name}: advanced panel starts collapsed")
    if r.slider.isVisible() or r.box.isVisible():
        fails.append(f"{r.terminal.name}: offset control shown on a region bound")
    print(f"    {r.terminal.name:14s} tilt=({r.tilt_a_b.value():.0f}, "
          f"{r.tilt_b_b.value():.0f})  plane size={r.band_b.value():.1f} mm  "
          f"panel open={r.adv.isVisible()}")

# ---- set a tilt and a plane size, and check they reach the cut --------------
row = w.bound_rows[0]
row.tilt_a_b.setValue(7.0)
row.band_b.setValue(round(row.band_b.value() + 4.0, 1))
for _ in range(10):
    app.processEvents()
want_tilt, want_band = 7.0, row.band_b.value()
print(f"\n  set {row.terminal.name}: tilt A = {want_tilt}°, "
      f"plane size = {want_band:.1f} mm")

acted = [t for t in w.cut_list()[0] if t.name == row.terminal.name]
if not acted:
    fails.append("the adjusted bound is not in the cut list")
elif abs(float(acted[0].tilt_a_deg) - want_tilt) > 1e-6:
    fails.append(f"cut list carries tilt {acted[0].tilt_a_deg}, set {want_tilt}")
else:
    print(f"  cut list carries tilt {acted[0].tilt_a_deg:.0f}°, "
          f"band {acted[0].band_mm}")

# Reaching the cut list is not the same as reaching the CUTTER: a control wired to
# an attribute the cutter never reads would pass every check above and still do
# nothing. band_width() is what cut_terminal asks, so ask it the same way.
#
# What is NOT asserted, and why. An earlier version of this demanded that tilting
# change the resulting mesh. Both halves of that are wrong:
#
#   * Tilting a boundary on a 15 mm-radius aorta tears the rim past the 0.5 mm
#     flatness bar, and the cutter REFUSES -- "the opening is not flat, its rim
#     strays 0.6 mm off the cut plane". That is the cutter working, not failing.
#     JK5's aortic neck refuses from 3 deg outward.
#   * Shrinking the plane size legitimately changes nothing here. band_w only gates
#     `near_cut` (branch_cut.py), which limits how far a cut may reach onto
#     ambiguous faces to close its rim; everything past the cut by arc length goes
#     regardless. On a vessel that is not reach-limited the mesh is identical at
#     37.9 mm and at 6 mm -- and identically so for an ordinary branch row.
#
# So the honest assertion is that the setting arrives at the cutter, and that a
# refusal (when one comes) names this boundary rather than something else.
live = next((x for x in w.cut_list()[0] if x.name == row.terminal.name), None)
print(f"  band_width(cutter sees) = {band_width(live):.2f} mm  (set {want_band:.1f})")
if abs(band_width(live) - want_band) > 1e-6:
    fails.append(f"cutter sees band {band_width(live)}, set {want_band}")

for probe_tilt in (0.0, 12.0):
    row.tilt_a_b.setValue(probe_tilt)
    for _ in range(10):
        app.processEvents()
    live = next((x for x in w.cut_list()[0] if x.name == row.terminal.name), None)
    verdict = w.problem or f"ready, {len(w.cut.faces)} faces"
    print(f"  tilt {probe_tilt:5.1f}° -> cutter sees {live.tilt_a_deg:.0f}°, "
          f"{verdict[:60]}")
    if abs(float(live.tilt_a_deg) - probe_tilt) > 1e-6:
        fails.append(f"cutter sees tilt {live.tilt_a_deg}, set {probe_tilt}")
    # A refusal must not blame a DIFFERENT vessel. It may perfectly well be a
    # whole-model verdict -- "3 openings, expected 2" names nothing because nothing
    # is what it is about, and a steep tilt on a boundary really can sever a
    # neighbour. Requiring the terminal's name in every message treated those as
    # bugs.
    #
    # This assertion only ever passed because tilts were not reaching the cut at all:
    # the re-cut cache keyed on the offset alone, so changing a tilt reused the old
    # mesh and nothing could break. Fixing the cache is what exposed it.
    other = [x.name for x in w.cut_list()[0] if x.name != row.terminal.name]
    if w.problem and any(w.problem.startswith(nm) for nm in other):
        fails.append(f"tilting {row.terminal.name} was blamed on another vessel: "
                     f"{w.problem[:60]}")
row.tilt_a_b.setValue(want_tilt)
for _ in range(10):
    app.processEvents()

# ---- THE TRAP: nudge a region slider, which rebuilds the bounds -------------
before = (w.reg_s0_b.value(), w.reg_s1_b.value())
w.reg_s0_b.setValue(before[0] + 2.0)
for _ in range(15):
    app.processEvents()
w._region_commit()
for _ in range(15):
    app.processEvents()

kept = [t for t in w.cut_list()[0] if t.name == row.terminal.name]
print(f"\n  after nudging Start {before[0]:.1f} -> {w.reg_s0_b.value():.1f} mm")
if not kept:
    fails.append("the bound vanished from the cut list after a region nudge")
else:
    got_t, got_b = float(kept[0].tilt_a_deg), kept[0].band_mm
    print(f"    tilt A survives : {got_t:.0f}°  (set {want_tilt:.0f}°)")
    print(f"    plane size      : {got_b}  (set {want_band:.1f})")
    if abs(got_t - want_tilt) > 1e-6:
        fails.append(f"tilt lost on region rebuild: {got_t} not {want_tilt}")
    if got_b is None or abs(float(got_b) - want_band) > 1e-6:
        fails.append(f"plane size lost on region rebuild: {got_b} not {want_band}")

# ---- leaving the region must NOT carry the tilt onto the whole vessel -------
w._set_region(None)
for _ in range(10):
    app.processEvents()
if w.bound_rows:
    fails.append(f"{len(w.bound_rows)} bound row(s) left after leaving the region")
if getattr(w, "_bound_planes", {}):
    fails.append("boundary plane settings outlived their region")
print(f"\n  back to whole vessel: {len(w.rows)} rows, "
      f"{len(w.bound_rows)} bound row(s), status {w.problem or 'ready'}")

print()
if fails:
    print("FAILURES:")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
