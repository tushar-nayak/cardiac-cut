"""Reset must put back the model as loaded, in the API as in the window.

branch_gui.CutWindow.reset says why, and the reasoning is not desktop-specific:

    "Reset is where a stuck user goes, so it puts EVERYTHING back: the region, the
     trims, the hand-declared arteries and every plane. Trims used to survive, on
     the reasoning that a trimmed group is a structural choice rather than a viewing
     mode. That reasoning is wrong from the chair: someone who presses Reset wants
     the model they loaded, and leaving three branches collapsed into one makes the
     button look broken in precisely the situation it exists for."

The web's reset_all cleared offsets, planes and the region, and left trims and
hand-declared arteries standing -- so the same button did two different things
depending on which front end you were in. Found by the 28 Aug parity audit.

The scope now lives in CutSession.reset_all, called by both, so they cannot drift
again.

    .venv\\Scripts\\python tests\\test_api_reset_all.py
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from vascular_api import CutSession                                  # noqa: E402

fails = []
s = CutSession(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
               str(ROOT / "JK5_1_M10.mat"))
names_at_load = dict(s.anat_names)

# ---- build up every kind of state a Reset is supposed to undo ---------------
cands = s.trim_candidates()
if not cands:
    print("FAILURES: nothing trimmable on JK5_1, cannot test")
    raise SystemExit(1)
s.add_trim(cands[0]["segment"])

lbl = next(o["label"] for o in s.openings()
           if not o["is_region_bound"] and not o["is_trim"])
choices = s.anatomy_choices()
pick = next((c for c in choices if "renal" in c.lower()), choices[0])
s.declare_artery(lbl, pick)

s.apply_preset("Infrarenal aorta")
bound = s._region_terminals()[0][0]
s.set_plane(bound.label, tilt_a_deg=6.0)
s.set_offset(lbl, 12.0)
s.set_plane(lbl, tilt_a_deg=4.0)

print("=" * 76)
print("before reset")
print("=" * 76)
print(f"  trims            : {len(s._trims)}")
print(f"  declared         : {s.declared_arteries()}")
print(f"  region           : {s.region_name}  {s.region}")
print(f"  offsets set      : {len(s._offsets)}")
print(f"  boundary planes  : {s._bound_planes}")
before = dict(trims=len(s._trims), declared=len(s.declared_arteries()),
              offsets=len(s._offsets))
if not all(before.values()):
    fails.append(f"the test failed to set up state to clear: {before}")

# ---- reset ------------------------------------------------------------------
s.reset_all()
print()
print("=" * 76)
print("after reset")
print("=" * 76)
print(f"  trims            : {len(s._trims)}")
print(f"  declared         : {s.declared_arteries()}")
print(f"  region           : {s.region}")
print(f"  offsets set      : {len(s._offsets)}")
print(f"  boundary planes  : {s._bound_planes}")

if s._trims:
    fails.append(f"trims survived Reset: {[t['segment'] for t in s._trims]}")
if s.declared_arteries():
    fails.append(f"hand-declared arteries survived Reset: {s.declared_arteries()}")
if s.region is not None:
    fails.append(f"the region survived Reset: {s.region}")
if s._offsets:
    fails.append(f"offsets survived Reset: {s._offsets}")
if s._bound_planes:
    fails.append(f"boundary planes survived Reset: {s._bound_planes}")

tilts = {t.label: float(getattr(t, "tilt_a_deg", 0.0) or 0.0)
         for t in s.tree.terminals}
hot = {k: v for k, v in tilts.items() if abs(v) > 1e-9}
print(f"  planes still tilted: {hot}")
if hot:
    fails.append(f"a hand-set tilt survived Reset: {hot}")

# ---- and the names must be re-derived, not left displaced -------------------
# declare_artery displaces the automatic name of the same artery, so a Reset that
# drops the declaration without re-deriving would leave that opening unnamed.
print()
print(f"  names at load {len(names_at_load)}, after reset {len(s.anat_names)}")
if s.anat_names != names_at_load:
    only_before = set(names_at_load) - set(s.anat_names)
    only_after = set(s.anat_names) - set(names_at_load)
    fails.append(f"the automatic names were not restored: "
                 f"lost {sorted(only_before)}, gained {sorted(only_after)}")

# ---- the cut must still be usable afterwards --------------------------------
print(f"  cut ok after reset : {s.ok}"
      + ("" if s.ok else f"  ({str(s.problem)[:60]})"))

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
