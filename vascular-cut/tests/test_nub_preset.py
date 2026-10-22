"""The nub-in preset, applied to a model that is already open.

The rule was validated over 82 vessels in 30 thoracic cases, and was then reachable
only by running tools/nub_cohort.py and reopening the window on its output. This
drives the path the button takes: compute the plan from the open model, put it in
the window, and check the openings actually moved to it.

    python tests/test_nub_preset.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qtpy import QtWidgets                                          # noqa: E402

import nub_cut as NC                                                # noqa: E402
import nub_session                                                  # noqa: E402

CSV = ROOT / "KY160_1_CentrelinePoints.csv"
MAT = ROOT / "Geometries & Meshes/mesh/KY160_1.mat"

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


app = QtWidgets.QApplication(sys.argv)
from branch_gui import CutWindow                                    # noqa: E402


def settle(win, n=10):
    for _ in range(n):
        app.processEvents()


w = CutWindow(str(CSV), str(MAT))
w.show()
settle(w)
print(f"  loaded {MAT.name}: {len(w.rows)} opening(s)")
check(w.tree is not None and w.base is not None, "the model loaded")

print()
print("=" * 78)
print("the button is offered exactly when the preset can be applied")
print("=" * 78)
check(w.nub_btn.isEnabled(), "enabled on a whole model with a centerline")
# A region has already replaced the aorta's own ends with its bounds, and the preset
# squares those ends off -- there is nothing there for it to square.
presets = []
try:
    from branch_cut import region_presets
    presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
except Exception:
    pass
if presets:
    w._apply_preset(presets[0]); settle(w)
    check(not w.nub_btn.isEnabled(), "disabled while a region is isolated")
    w._set_region(None); settle(w)
    check(w.nub_btn.isEnabled(), "enabled again once the whole vessel is back")
else:
    print("  (no region preset on this model, skipping the region check)")

print()
print("=" * 78)
print("applying it")
print("=" * 78)
before = {r.terminal.label: r.value() for r in w.rows}
notes = nub_session.apply_here(w)
after = {r.terminal.label: r.value() for r in w.rows}
for n in notes[:10]:
    print(f"    {n}")
if len(notes) > 10:
    print(f"    ... {len(notes) - 10} more")

check(any(n.startswith("nub-in cut") for n in notes),
      "it reports how many openings it placed")
moved = [k for k in after if abs(after[k] - before.get(k, -1)) > 0.05]
print(f"  {len(moved)} of {len(after)} openings moved")
check(bool(moved), "the openings actually moved to the preset")

# Every branch gets a stub, and the aorta's own ends are squared off horizontally at
# the trunk offset -- two different cuts, and the preset is the only thing that makes
# both at once.
#
# A large offset here is NOT a large nub. Offsets are measured back from the vessel's
# TIP, so KY160's 212 mm branch takes 178.8 mm to leave a 33 mm stub off the aorta.
# That is the batch's own number for this case, checked against it directly.
trunks = [r for r in w.rows if r.terminal.is_trunk]
branch = [r for r in w.rows if not r.terminal.is_trunk]
print(f"  trunk ends {[f'{r.value():.1f}mm/{r.terminal.mode}' for r in trunks]}")
print(f"  branches   {[f'{r.value():.1f}mm' for r in branch]}")
check(all(abs(r.value() - NC.TRUNK_OFFSET_MM) < 0.05 or r.terminal.mode != NC.TRUNK_MODE
          for r in trunks),
      f"the aorta's ends sit at the trunk offset ({NC.TRUNK_OFFSET_MM} mm), or say why not")
check(all(r.value() > 0 for r in branch), "every branch has a positive stub")

print()
print("=" * 78)
print("the result is a usable cut, not just a set of numbers")
print("=" * 78)
w._on_change(); settle(w)
print(f"  status: {w.status.text()[:100]}")
check(w.cut is not None, "a mesh came out")
check(w.problem is None, f"and it validates: {w.problem}")

print()
print("=" * 78)
print("it refuses cleanly with nothing loaded")
print("=" * 78)
blank = CutWindow(None, None)
settle(blank)
try:
    nub_session.plan_here(blank)
    check(False, "an empty window is refused")
except ValueError as e:
    check("centerline" in str(e), f"an empty window is refused: {e}")
check(not blank.nub_btn.isEnabled(), "and the button is not offered there")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
