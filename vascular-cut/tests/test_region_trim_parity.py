"""The API and the window must cut the SAME LIST inside a region. argv[1]=jk5|ky160

THE DIVERGENCE. branch_gui.region_rows() replaces each trimmed group with its single
interior cut before a region's rows are built, and says so -- it even records that a
region used to drop trims entirely, "so trimming the iliacs while looking at the
infrarenal aorta recorded the trim and changed nothing". vascular_api._cut_list simply
did not make that call. With a region active the API therefore cut a DIFFERENT MESH
from the window for the same session: every trimmed group came back as separate
openings instead of one.

WHY IT MATTERS BEYOND TIDINESS. vascular_api's own header sets the contract -- "the
window remains the reference implementation: where the two could drift, this file
follows branch_gui.CutWindow" -- and the website drives the API, so the two would have
shipped different geometry for the same choices. It is also what an arch-plus-nub
cohort needs: a nub on a vessel that splits IS a trim, so a batch runner using the API
would have produced arch meshes with no nubs on them and no complaint.

WHAT IS ASSERTED. Not that the meshes are identical -- that would drag the whole Qt
render path in -- but that both agree on WHICH cuts are being made and where, with a
region active and a trim inside it, which is the thing that had drifted.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets                                          # noqa: E402
from branch_gui import CutWindow                                    # noqa: E402
from vascular_api import CutSession                                 # noqa: E402

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]
fails = []


def _is_trim(t):
    """A trim, not a region bound. Both set is_interior; only a bound carries
    _region_offset, so that is what tells them apart."""
    return (bool(getattr(t, "is_interior", False))
            and getattr(t, "_region_offset", None) is None)


print("=" * 78)
print(f"region + trim parity , {which.upper()}")
print("=" * 78)

s = CutSession(str(ROOT / csv), str(ROOT / mat))
cands = s.trim_candidates()
presets = [p["name"] for p in s.presets()]
print(f"  {len(cands)} trim candidate(s), {len(presets)} region preset(s)")
if not cands or not presets:
    print("  (this model offers no trim or no region; nothing to compare)")
    print("\nALL PASS")
    raise SystemExit(0)

# The trim whose group actually falls inside a region, so the comparison is about the
# case that drifted rather than one where the trim is out of view either way.
app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
win = CutWindow(str(ROOT / csv), str(ROOT / mat))
win.show()
for _ in range(10):
    app.processEvents()

# CHOSEN WITHOUT ASKING _cut_list, which is the thing under test. Selecting the case
# by "does a trim survive the cut list" made this pass VACUOUSLY against the bug: with
# trims dropped inside a region there was never a case to compare, so it reported "no
# trim falls inside any region" and exited green. The region's own membership --
# _region_terminals, which the bug does not touch -- decides instead.
found = None
for cand in cands:
    for name in presets:
        s.apply_preset(name)
        inside_segs = {t.segment.name for t in s._region_terminals()[1]}
        s.clear_region()
        if set(cand["subtree"] if "subtree" in cand else []) <= inside_segs and                 {c["segment"] for c in [cand]} <= inside_segs | {cand["segment"]}:
            pass
        # the group this trim swallows, by segment, has to lie inside the region
        grp = {t.segment.name for t in (cand.get("terminals") or [])} or None
        if grp is None:
            import branch_cut as _bc
            full = next((x for x in _bc.pre_bifurcation_candidates(s.tree, s.base)
                         if x["segment"] == cand["segment"]), None)
            grp = {t.segment.name for t in (full or {}).get("terminals", [])}
        if grp and grp <= inside_segs:
            found = (cand, name, None)
            break
    if found:
        break

if found is None:
    print("  no trim falls inside any region on this model")
    print("\nALL PASS")
    raise SystemExit(0)

cand, region_name, api_labels = found
print(f"  comparing: trim on {cand['segment']} ({cand['group'][:36]}) "
      f"inside region '{region_name}'")

s.add_trim(cand["segment"])
s.apply_preset(region_name)
api_active, api_offs = s._cut_list()

gui_cand = next((c for c in __import__("branch_cut").pre_bifurcation_candidates(
    win.tree, win.base) if c["segment"] == cand["segment"]), None)
if gui_cand is None:
    fails.append(f"the window does not offer a trim on {cand['segment']} at all")
else:
    win._add_trim(gui_cand)
    span = next((p for p in win.region_presets_list()
                 if p["name"] == region_name), None) \
        if hasattr(win, "region_presets_list") else None
    win._set_region((s.region[0], s.region[1]), name=region_name)
    for _ in range(10):
        app.processEvents()
    gui_active, gui_offs = win.cut_list()

    api_names = sorted(str(t.name) for t in api_active)
    gui_names = sorted(str(t.name) for t in gui_active)
    print(f"\n    API    ({len(api_names):2d}): {', '.join(n[:22] for n in api_names)}")
    print(f"    window ({len(gui_names):2d}): {', '.join(n[:22] for n in gui_names)}")
    if api_names != gui_names:
        only_api = sorted(set(api_names) - set(gui_names))
        only_gui = sorted(set(gui_names) - set(api_names))
        fails.append(f"the API and the window cut different lists inside a region: "
                     f"only the API cuts {only_api}, only the window cuts {only_gui}")

    n_api_trim = sum(1 for t in api_active if _is_trim(t))
    n_gui_trim = sum(1 for t in gui_active if _is_trim(t))
    print(f"    interior (trim) cuts: API {n_api_trim}, window {n_gui_trim}")
    if n_api_trim != n_gui_trim:
        fails.append(f"trims inside the region: API {n_api_trim}, window {n_gui_trim} "
                     f"-- a trim that vanishes in one path changes the mesh silently")
    if n_api_trim == 0:
        fails.append("the trim was dropped inside the region by BOTH paths, so this "
                     "test is no longer exercising the case it exists for")

print()
if fails:
    print("FAILURES:")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
