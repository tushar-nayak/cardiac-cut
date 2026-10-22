"""A model end that will not cut flat is trimmed shallower, and says so.

THE PROBLEM. Four abdominal inlets refused at the default 6 mm with "the opening is
not flat, its rim strays N mm off the cut plane", and every one of them cuts EXACTLY
flat a few millimetres shallower. The deviation grows monotonically with depth:

    offset      1     2     3     4     5   6 (default)  7
    JK14_3   0.00  0.00  0.00  0.98  1.15    1.40       2.38
    JK23_3   0.00  0.00  0.00  0.70  1.26    2.25        --
    JK30_3   0.00  0.00  0.00  0.91  1.73    2.54        --
    JK6_2    0.00  0.00  0.00  0.00  0.00    0.00       1.02

The cause is geometric, and four separate face-selection rules were built against it
before this was understood. An inlet segment is 15 to 22 mm long, and the deeper the
plane travels the more likely it crosses a NEIGHBOURING branch that runs upward. That
branch's collar borders the stub being removed, so taking it cuts into the branch and
leaving it puts its above-plane corners into the rim -- one to six faces, 0.3 to
1.6 mm proud, which is exactly the deviation reported.

WHY MOVING THE OFFSET IS ALLOWED HERE AND NOWHERE ELSE. At a model end the offset
says how much of the end stub to trim, and no anatomy depends on the answer -- a
shallower cut keeps MORE vessel. Every other opening sits at a real anatomical
station, where moving the plane would change what is measured.

WHAT IS ASSERTED, and what is NOT. The load-bearing half is is_model_end: a branch
tip is a free end of the model too, and letting this rule reach branches is exactly
what cost 840,623 mm2 of foreign surface against 19,375 in an earlier attempt. That
discrimination -- trunk yes; branch, trim, region bound, junction no -- is checked
exhaustively here, on a synthetic tree, deterministically.

The RESCUE itself needs geometry neither shipped model has: both JK5_1 and KK2_5 cut
flat at every depth their trunk allows, so there is nothing for the backoff to find.
Section 2 sweeps for a torn depth and says so when it finds none, rather than passing
on a function that did nothing. The rescue is exercised by tools/smoke_bench.py over
the cohort, where it takes four refusing inlets to clean cuts.

The last assertion is the one worth keeping honest: when a cut IS moved, it must say
so. A cut placed somewhere other than where the caller asked is not reproducible, and
a status line reading "Ready" either way is how that goes unnoticed.
"""
import sys
from pathlib import Path
import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from branch_tree import Segment, Terminal, Tree                    # noqa: E402
import branch_cut as bc                                            # noqa: E402
from vascular_api import CutSession                                # noqa: E402

fails = []
print("=" * 78)
print("trunk-end offset backoff")
print("=" * 78)

# ---- 1. is_model_end tells a trunk end from everything that looks like one ---
z = np.linspace(0.0, 40.0, 41)
pts = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
seg = Segment("A", pts, np.tile([0.0, 0.0, 1.0], (41, 1)))
tree = Tree({"A": seg}, [], [])
t = Terminal(seg, True, label="A.end")
tree.terminals = [t]

t.is_trunk = True
if not bc.is_model_end(tree, t):
    fails.append("a free trunk end is not recognised as a model end")
t.is_trunk = False
if bc.is_model_end(tree, t):
    fails.append("a BRANCH tip is treated as a model end -- it is a free end too, "
                 "and letting the rule reach branches is what cost 840,623 mm2")
t.is_trunk = True
t.subtree_segments = {"B"}
if bc.is_model_end(tree, t):
    fails.append("a TRIM is treated as a model end; it is meant to swallow openings")
t.subtree_segments = None
t._region_offset = 4.0
if bc.is_model_end(tree, t):
    fails.append("a REGION BOUND is treated as a model end; it is a plane inside "
                 "the vessel, not at its end")
t._region_offset = None
# ...and an end that is a junction is not a model end: something continues there.
tree.junctions = [(pts[-1], [("A", True)])]
tree._segnodes_cache = None
if bc.is_model_end(tree, t):
    fails.append("an end at a JUNCTION is treated as a model end")
print("  is_model_end: trunk free end yes; branch, trim, region bound, junction no")

# ---- 2. it finds a flat offset where a deep one is torn ----------------------
#
# Found on the model rather than asserted from a table: sweep the trunk end until an
# offset comes out over the flatness bar, then require the backoff to return one
# below it that does not. Skips honestly if this model has no torn depth to find,
# rather than passing on a function that did nothing.
MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "kk2": ("KK2_5.csv", "KK2_5.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]
s = CutSession(str(ROOT / csv), str(ROOT / mat))
trunk = next((x for x in s.tree.terminals if bc.is_model_end(s.tree, x)), None)
print(f"\n  {which.upper()}: model end = {trunk.name if trunk else None}")
if trunk is None:
    fails.append("no model end found on this model at all")
else:
    attr = bc.vessel_attribution(s.base, s.tree)

    def rim_dev(off):
        try:
            out = bc.cut_terminal(s.base, s.tree, trunk, off, _attr=attr,
                                  _suggest=False)
        except ValueError:
            return None
        loops = bc.boundary_loops(out)
        if not loops:
            return None
        pt, nrm = trunk.cut_plane(off)
        v = np.asarray(out.vertices)
        rim = min(loops, key=lambda L: float(
            np.linalg.norm(v[np.asarray(L)].mean(axis=0) - pt)))
        return float(np.abs((v[np.asarray(rim)] - pt) @ nrm).max())

    torn = None
    for off in (6.0, 9.0, 12.0, 16.0, 20.0, 26.0, 32.0):
        if off > 0.9 * trunk.segment.length:
            break
        d = rim_dev(off)
        if d is not None and d > bc.RIM_FLATNESS_MM:
            torn = (off, d)
            break
    if torn is None:
        print("    no torn depth on this model -- backoff has nothing to rescue here")
    else:
        got = bc.flat_trunk_offset(s.base, s.tree, trunk, torn[0], attr=attr)
        after = rim_dev(got) if got is not None else None
        print(f"    {torn[0]:.0f} mm gives {torn[1]:.2f} mm off plane; "
              f"backoff -> {got}  giving {after}")
        if got is None:
            fails.append(f"no flat offset found below {torn[0]:.0f} mm, though the "
                         f"rim there is {torn[1]:.2f} mm off plane")
        else:
            if got > torn[0] + 1e-9:
                fails.append(f"backoff went DEEPER ({got}) than the offset that tore "
                             f"({torn[0]}); it must keep more vessel, not less")
            if after is None or after > bc.RIM_FLATNESS_MM:
                fails.append(f"backoff returned {got} whose rim is {after}, still "
                             f"over the {bc.RIM_FLATNESS_MM} mm bar")

# ---- 3. a moved cut is REPORTED ---------------------------------------------
st = s.recut()
moved = st.get("trunk_offset_moved")
print(f"\n  {which.upper()}: ok={st['ok']}  moved={moved}")
if moved is not None:
    if str(moved.get("name", "")) not in st["message"]:
        fails.append("a cut was moved and the status line does not name it")
    if f"{moved['to_mm']:g}" not in st["message"]:
        fails.append("a cut was moved and the status line does not give the offset")
    print(f"    reported: {st['message'][:96]}")
else:
    print("    (nothing needed moving on this model, which is the common case)")
if st["ok"] and "trunk_offset_moved" not in st:
    fails.append("status() does not carry trunk_offset_moved at all")

print()
if fails:
    print("FAILURES:")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
