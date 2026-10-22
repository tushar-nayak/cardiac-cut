"""A marginal abdominal/thoracic call has to say it is marginal.

The desktop window has reported close calls since the family rule went in. The API
did not, so the web app asserted a marginal reading as a fact.

That matters more than a missing caption. Measured on the centerlines the tool
ACTUALLY cuts -- smoothed when jagged, which is 52 of 79 cases -- the two families
overlap: an abdominal aorta reaches 0.225 and a thoracic one 0.249, either side of a
0.25 cut-off. No threshold separates them, so 77 of 79 is this rule's ceiling rather
than a tuning failure, and the two it misses are exactly the marginal ones. Saying
"close" is the remedy, and setting the family by hand is the fix.

JK25_1 is the abdominal aorta at 0.225 that reads as thoracic.

    python tests/test_family_close_call.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "webapp"))

from vascular_api import (CutSession, FAMILY_CLOSE_BAND)             # noqa: E402
from branch_cut import PAIR_AREA_FRACTION                            # noqa: E402

CSV = ROOT / "JK5_1_CentrelinePoints_branches.csv"
MAT = ROOT / "JK5_1_M10.mat"

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


print("=" * 78)
print("a confident call is not flagged")
print("=" * 78)
s = CutSession(str(CSV), str(MAT))
a = s.anatomy_summary()
print(f"  {a['kind']}, split {a['outflow_split']} against {a['family_cutoff']}")
check(a["kind"] == "abdominal", "JK5 reads as abdominal")
check(a["outflow_split"] is not None, "the evidence behind the call is reported")
check(not a["close_call"],
      f"and a split of {a['outflow_split']} is not called close")

print()
print("=" * 78)
print("a call near the cut-off IS flagged")
print("=" * 78)
# Forced rather than fished for: the band is what decides, so drive it directly.
s.tree.outflow_split = PAIR_AREA_FRACTION - FAMILY_CLOSE_BAND / 2.0
near = s.anatomy_summary()
print(f"  split {near['outflow_split']} -> close_call={near['close_call']}")
check(near["close_call"], "a split half a band from the cut-off is called close")
s.tree.outflow_split = PAIR_AREA_FRACTION - FAMILY_CLOSE_BAND * 2.0
far = s.anatomy_summary()
print(f"  split {far['outflow_split']} -> close_call={far['close_call']}")
check(not far["close_call"], "and one two bands away is not")

print()
print("=" * 78)
print("a family the user DECLARED is never called close")
print("=" * 78)
# They said so. Second-guessing a declaration with "but it was marginal" is noise:
# the margin describes a reading nobody is relying on any more.
d = CutSession(str(CSV), str(MAT), family="abdominal")
d.tree.outflow_split = PAIR_AREA_FRACTION - FAMILY_CLOSE_BAND / 2.0
da = d.anatomy_summary()
print(f"  declared {da['kind']} ({da['how']}), close_call={da['close_call']}")
check(da["declared"], "the declaration is recorded")
check(not da["close_call"], "and no close-call warning is raised against it")

print()
print("=" * 78)
print("the web caption leads with it")
print("=" * 78)
import cut_app                                                       # noqa: E402
s.tree.outflow_split = PAIR_AREA_FRACTION - FAMILY_CLOSE_BAND / 2.0
text = cut_app._anat_caption(s)
print(f"  {text[:104]}")
check("CLOSE call" in text, "the caption says the call was close")
check("by hand" in text, "and points at the way to fix it")
s.tree.outflow_split = 0.995
plain = cut_app._anat_caption(s)
check("CLOSE call" not in plain,
      f"a confident call gets the ordinary caption: {plain[:56]}")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
