"""A region boundary that undercuts a branch is nudged clear of it, and says so.

THE PROBLEM. Twelve region presets refused with "N openings, expected N-1" -- one
opening too many. The extra one is not a torn rim and not a severed neighbour. On
JK17_2's infrarenal bound the aortic cross-section is 163 points, 15.46 mm across and
0.00 mm off its plane, which is flawless, and beside it sits a 19-point hole 1.79 mm
across, 35 mm away. That second loop is a BRANCH'S OSTIUM: the boundary sits a
millimetre or two from where a branch leaves, so the removed side takes aortic wall
from around the base of a vessel the region KEEPS, and the hole it leaves is an
opening nobody asked for.

WHY NOTHING IN THE CUTTER FIXED IT. The faces taken are the cut's own -- measured,
zero foreign surface -- so no ownership rule applies. Four things were tried and
ruled out before the cause was found, and they are recorded so they are not retried:

    the plane's reach        swept 55 mm down to 6 mm on six presets, no effect
    ownership of the faces   zero foreign in every case
    subtree removed in part  taking the subtree whole changes nothing
    subtree membership       no subtree segment attaches inside the region

WHAT FIXES IT. One to three millimetres of clearance:

    shift      0    +1   +2   +3   +5
    JK17_2     2     1    1    1    1      openings left by the bound's own cut
    JK30_5     2     2    1    1    1
    JK4_2      2     2    2    1    1

Abdominal regions went 190 of 210 to 206 of 215 end to end, refusals 20 to 9.

WHAT IS ASSERTED HERE. The two things that make it safe rather than merely effective:
a nudge must SURVIVE region_plan rebuilding its Terminals -- it does not live on the
objects, which are thrown away and rebuilt on every call, and mutating them instead
was the first attempt and did nothing at all -- and a boundary that MOVED must say so,
because it changes which stretch of aorta was measured.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from vascular_api import CutSession                                # noqa: E402
import branch_cut as bc                                            # noqa: E402

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]
fails = []
print("=" * 78)
print(f"region boundary clearance , {which.upper()}")
print("=" * 78)

s = CutSession(str(ROOT / csv), str(ROOT / mat))
presets = [p["name"] for p in s.presets()]
print(f"  {len(presets)} preset(s): {', '.join(p[:22] for p in presets)}")

# ---- 1. a nudge survives region_plan rebuilding the bounds -------------------
#
# THE TRAP, and the first version fell straight into it. region_plan builds fresh
# Terminal objects every time _cut_list() is called, so a nudge written onto a bound
# is gone by the next re-cut with nothing to show it ever applied. It has to be
# carried on the session and reapplied, keyed by the bound's ROLE -- the label names
# whichever segment the plane happens to land on, which changes as the boundary moves.
if not presets:
    print("  (no region presets on this model)")
else:
    s.apply_preset(presets[0])
    bounds = [t for t in s._cut_list()[0]
              if getattr(t, "_region_offset", None) is not None]
    if not bounds:
        fails.append("a region preset produced no boundary terminals")
    else:
        role = bounds[0].name
        was = float(bounds[0]._region_offset)
        s._bound_nudge = {role: was + 2.0}
        again = [t for t in s._cut_list()[0]
                 if getattr(t, "_region_offset", None) is not None]
        got = next((float(t._region_offset) for t in again if t.name == role), None)
        print(f"\n  {role}: planned at {was:.2f} mm, nudged to {was + 2.0:.2f}, "
              f"after rebuild {got}")
        if got is None or abs(got - (was + 2.0)) > 1e-6:
            fails.append(f"the nudge did not survive region_plan rebuilding the "
                         f"bounds: {got} not {was + 2.0}")
        s._bound_nudge = {}

        # ...and it must not touch a bound it was not asked about
        other = [t for t in again if t.name != role]
        if other:
            print(f"  {other[0].name}: untouched at "
                  f"{float(other[0]._region_offset):.2f} mm")

# ---- 2. the ladder: planned offset first, both ways, and bounded -------------
#
# THREE properties, each of which was wrong in a version of this.
#
#   0.0 FIRST, so a boundary that already cuts cleanly is left exactly where it was
#     planned. Without it the search moves every bound it is handed, including ones
#     with nothing wrong -- JK30_5's region end was being shifted a millimetre for
#     no reason.
#   BOTH DIRECTIONS. Outward-only left rescuable bounds refusing, because which side
#     is clear depends on which side the branch is on: JK25_3's aortic neck needs
#     -6 mm and refuses outright at +1 and +2; JK6_1's visceral bound needs -8.
#   BOUNDED. A boundary is a station on the aorta, so an unbounded search would
#     relocate the region rather than clear an ostium beside it.
if bc.BOUND_NUDGE_MM[0] != 0.0:
    fails.append(f"BOUND_NUDGE_MM does not try the planned offset first: "
                 f"{bc.BOUND_NUDGE_MM[:3]} -- a bound with nothing wrong gets moved")
if not (any(d > 0 for d in bc.BOUND_NUDGE_MM)
        and any(d < 0 for d in bc.BOUND_NUDGE_MM)):
    fails.append("BOUND_NUDGE_MM searches only one direction; which side is clear "
                 "depends on which side the branch is on")
if max(abs(d) for d in bc.BOUND_NUDGE_MM) > 10.0:
    fails.append(f"BOUND_NUDGE_MM reaches "
                 f"{max(abs(d) for d in bc.BOUND_NUDGE_MM)} mm, far enough to move "
                 f"the region rather than clear an ostium")
steps = sorted({abs(d) for d in bc.BOUND_NUDGE_MM})
print(f"\n  BOUND_NUDGE_MM: planned offset first, both ways, up to "
      f"{max(steps):g} mm ({len(bc.BOUND_NUDGE_MM)} steps)")

# ---- 3. a moved boundary is REPORTED ----------------------------------------
st = s.recut()
moved = st.get("bounds_moved")
print(f"  recut: ok={st['ok']}  bounds_moved={moved}")
if "bounds_moved" not in st:
    fails.append("status() does not carry bounds_moved at all")
if moved:
    for b in moved:
        if b["name"] not in st["message"]:
            fails.append(f"{b['name']} moved and the status line does not name it")
    print(f"    reported: {st['message'][:100]}")
else:
    print("    (nothing needed moving on this model, which is the common case)")
s.clear_region()

print()
if fails:
    print("FAILURES:")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
