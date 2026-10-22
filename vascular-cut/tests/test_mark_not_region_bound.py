""""Send cut to mark" must not land on a region boundary and claim success.

Found by the 28 Aug desktop/web parity audit.

vascular_api.offset_for_mark already tried to exclude region boundaries:

    if x.segment.name == seg and not getattr(x, "is_region_bound", False)

but NOTHING EVER SET that attribute. It existed only as a dict key that openings()
computes on the fly (`id(t) not in real`), so the getattr always returned False and
the guard was inert.

The consequence is a control that lies. A region boundary is not a vessel end: its
depth comes from region_plan via `_region_offset`, and `_cut_list` overwrites the
boundary's offset with that value on every cut. So sending a cut to a mark on the
boundary's segment set an offset that was immediately discarded, re-cut, and reported
"moved to N mm, where the mark is" -- having moved nothing. The same class of fault
the cohort reported against the region sliders.

Fixed at the source: region_plan now sets `is_region_bound = True` on the boundary
Terminals it manufactures, because that is the one place a boundary is made.

    .venv\\Scripts\\python tests\\test_mark_not_region_bound.py
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from vascular_api import CutSession                                  # noqa: E402

fails = []
s = CutSession(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
               str(ROOT / "JK5_1_M10.mat"))

# The open-ended preset, so there are real openings INSIDE as well as a boundary.
OPEN = "Infrarenal aorta with iliacs"
if not any(p["name"] == OPEN for p in s.presets()):
    print(f"FAILURES: {OPEN} not offered on JK5_1")
    raise SystemExit(1)
s.apply_preset(OPEN)
bounds, inside, _ = s._region_terminals()

print("=" * 76)
print(f"region {OPEN!r}")
print("=" * 76)
print(f"  boundaries : {[(t.name, getattr(t, 'is_region_bound', 'MISSING')) for t in bounds]}")
print(f"  inside     : {[t.name for t in inside]}")

# ---- 1. the flag exists on boundaries and not on ordinary openings ----------
for t in bounds:
    if getattr(t, "is_region_bound", False) is not True:
        fails.append(f"boundary {t.name!r} does not carry is_region_bound")
for t in inside:
    if getattr(t, "is_region_bound", False):
        fails.append(f"ordinary opening {t.name!r} is flagged as a region bound")

# ---- 2. a mark on a boundary's segment must not move the boundary -----------
print()
print("=" * 76)
print("a mark on a BOUNDARY's segment")
print("=" * 76)
def looks_like_a_bound(label):
    """Independent of the flag under test, on purpose.

    region_plan labels its boundaries "<segment>.regionA" / ".regionB" (branch_cut,
    one construction site). Checking is_region_bound here instead would make the
    test circular: without the fix that attribute is missing, so a boundary would
    read as an ordinary opening and the assertion could never fire -- which is
    exactly how the inert guard survived in the first place.
    """
    return str(label).endswith((".regionA", ".regionB"))


for b in bounds:
    mark = {"segment": b.segment.name, "s_mm": float(b.segment.length) * 0.5}
    plan = s.offset_for_mark(mark)
    print(f"  {b.name:<14} segment {b.segment.name:<10} -> {plan}")
    if plan is None:
        continue                       # refused outright, which is correct
    if looks_like_a_bound(plan["label"]):
        fails.append(f"{b.name}: send-to-mark targeted the region boundary "
                     f"{plan['label']!r}, whose offset _cut_list overwrites -- "
                     f"it would report success and move nothing")
    elif next((x for x in s._cut_list()[0] if x.label == plan["label"]), None) is None:
        fails.append(f"{b.name}: planned a cut on {plan['label']}, which is not "
                     f"in the cut list at all")

# ---- 3. and a mark on a real opening inside the region STILL works ----------
# The guard must not be so broad that the feature stops working in a region.
print()
print("=" * 76)
print("a mark on an ORDINARY opening inside the region")
print("=" * 76)
worked = 0
for t in inside:
    mark = {"segment": t.segment.name, "s_mm": float(t.segment.length) * 0.5}
    plan = s.offset_for_mark(mark)
    if plan is not None:
        worked += 1
        if looks_like_a_bound(plan["label"]):
            fails.append(f"a mark on {t.name!r} resolved to the boundary "
                         f"{plan['label']!r}")
print(f"  {worked} of {len(inside)} openings inside accepted a mark")
if inside and not worked:
    fails.append("no opening inside the region accepts a mark -- the guard is "
                 "too broad and the feature is dead in a region")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
