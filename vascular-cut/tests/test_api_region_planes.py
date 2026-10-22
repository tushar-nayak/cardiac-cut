"""A region boundary's plane settings must survive a re-cut in the API too.

tests/test_region_planes.py guards this for the desktop window and says why:

    "region_plan builds fresh Terminal objects on every call and is called again on
     every nudge of a region slider, so a tilt set by hand is discarded by the next
     drag unless it is carried across -- silently, with the plane snapping back
     while the user is looking at the number they typed."

But that is a CutWindow test, so it never covered vascular_api, and the API had the
trap unfixed. `terminal(label)` resolves a boundary through `_cut_list()`, which calls
region_plan, which manufactures the Terminal for that cut only -- so `set_plane` on a
boundary wrote to an object that was thrown away, and returned normally. The web draws
Tilt / Shift / Plane size / Reset for exactly those boundaries, so the controls were
live and inert: found by the 28 Aug desktop/web parity audit.

The API already carried `_region_reach` and `_bound_nudge` across the same rebuild,
and the comment beside them names this very trap -- so the mechanism was understood
and hand-set planes had simply been missed.

    .venv\\Scripts\\python tests\\test_api_region_planes.py
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from vascular_api import CutSession                                  # noqa: E402

fails = []
s = CutSession(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
               str(ROOT / "JK5_1_M10.mat"))

preset = next((p for p in s.presets() if p["name"] == "Infrarenal aorta"), None)
if preset is None:
    print("FAILURES: 'Infrarenal aorta' not offered on JK5_1")
    raise SystemExit(1)
s.apply_preset("Infrarenal aorta")


def bounds():
    return s._region_terminals()[0]


def by_role(attr):
    return {t.name: float(getattr(t, attr, 0.0) or 0.0) for t in bounds()}


b = bounds()[0]
role, label = b.name, b.label
print("=" * 76)
print(f"region 'Infrarenal aorta', boundary {label!r}  role {role!r}")
print("=" * 76)
print(f"  boundaries: {[t.name for t in bounds()]}")
print(f"  before    : tilt_a {by_role('tilt_a_deg')}")

# ---- 1. set a tilt and a plane size on the boundary -------------------------
s.set_plane(label, tilt_a_deg=7.0)
lim = s.plane_limits(label)
band_to = round(min(max(lim["band_min"] + 1.0, 9.0), lim["band_max"]), 3)
s.set_plane(label, band_mm=band_to)
print(f"  set       : tilt_a 7.0, band {band_to}")
if abs(by_role("tilt_a_deg").get(role, 0.0) - 7.0) > 1e-9:
    fails.append("the tilt did not take at all")

# ---- 2. THE TRAP: move the region, which rebuilds every boundary Terminal ---
s.set_region(preset["s0"] + 4.0, preset["s1"], "Infrarenal aorta", manual=True)
after = by_role("tilt_a_deg")
band_after = {t.name: getattr(t, "band_mm", None) for t in bounds()}
print()
print("  after moving the region (region_plan rebuilds the Terminals):")
print(f"    tilt_a  : {after}")
print(f"    band_mm : {band_after}")
if abs(after.get(role, 0.0) - 7.0) > 1e-9:
    fails.append(f"the boundary tilt was discarded by the re-cut: "
                 f"{after.get(role)} != 7.0")
if band_after.get(role) is None or abs(float(band_after[role]) - band_to) > 1e-6:
    fails.append(f"the boundary plane size was discarded: {band_after.get(role)}")

# ---- 3. reset must actually reset, not be re-applied from the record --------
now = next((t for t in bounds() if t.name == role), None)
if now is None:
    fails.append(f"the {role!r} boundary vanished after the move")
else:
    s.reset_plane(now.label)
    back = by_role("tilt_a_deg")
    band_back = {t.name: getattr(t, "band_mm", None) for t in bounds()}
    print()
    print(f"  after reset_plane: tilt_a {back}  band {band_back}")
    if abs(back.get(role, 0.0)) > 1e-9:
        fails.append("reset_plane did not clear the tilt -- the saved value was "
                     "re-applied over it")
    if band_back.get(role) is not None:
        fails.append("reset_plane did not restore the automatic plane size")

# ---- 4. leaving the region entirely drops the record ------------------------
# Carrying a tilt into the NEXT region would apply it to a plane somewhere else.
s.set_region(preset["s0"], preset["s1"], "Infrarenal aorta", manual=True)
s.set_plane(bounds()[0].label, tilt_a_deg=5.0)
s.clear_region()
s.apply_preset("Infrarenal aorta")
fresh = by_role("tilt_a_deg")
print()
print(f"  after clear_region and re-picking the preset: {fresh}")
if any(abs(v) > 1e-9 for v in fresh.values()):
    fails.append(f"a boundary tilt survived leaving the region: {fresh}")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
