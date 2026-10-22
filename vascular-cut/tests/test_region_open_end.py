"""A region boundary on the vessel's own end gets no plane -- it is already an opening.

"Whole abdominal aorta" runs inlet to bifurcation, so its proximal boundary sits
exactly on the inlet. Ishimaru zone 4 runs to the distal end. region_terminals clamps
those to MIN_OFFSET_MM, so the tool tried to slice a 1 mm wafer off the end of the
model and preflight refused it: "a cut of 1.0 mm sits on the vessel tip".

Across both cohorts 101 presets place a bound within 2 mm of a vessel end. Dropping
that plane fixes 20 of them and changes nothing on the other 79 -- a plane cutting a
wafer off an already-open end was never doing useful work.

    python tests/test_region_open_end.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from branch_tree import load_tree                                   # noqa: E402
from branch_cut import (load_mat_mesh, classify_terminals,          # noqa: E402
                        name_anatomy, region_presets, region_plan,
                        aorta_route, _junction_graph, MIN_OFFSET_MM)

CSV = ROOT / "JK5_1_CentrelinePoints_branches.csv"
MAT = ROOT / "JK5_1_M10.mat"

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


tree = load_tree(str(CSV))
mesh = load_mat_mesh(str(MAT))
classify_terminals(tree, mesh, verbose=False)
names, vessels, _ = name_anatomy(tree, mesh)
jg = _junction_graph(tree)
route = aorta_route(tree, jg)
presets = region_presets(tree, mesh, names, vessels)

print("=" * 78)
print("a region running to the vessel's end does not place a plane there")
print("=" * 78)
whole = next((p for p in presets if p["name"] == "Whole abdominal aorta"), None)
check(whole is not None,
      f"'Whole abdominal aorta' is still offered ({[p['name'] for p in presets]})")
if whole:
    s1 = route.length if whole["s1"] is None else whole["s1"]
    bounds, inside, removed, notes = region_plan(tree, whole["s0"], s1, route, jg,
                                                 mesh=mesh, clear_ostia=True)
    offs = [round(float(getattr(t, "_region_offset", -1)), 2) for t in bounds]
    print(f"  s0={whole['s0']:.1f} s1={s1:.1f} (route {route.length:.1f})")
    print(f"  {len(bounds)} bound(s), offsets {offs}")
    for n in notes:
        print(f"    note: {n}")
    check(len(bounds) < 2,
          f"the boundary on the inlet is dropped ({len(bounds)} plane(s), not 2)")
    check(all(o > MIN_OFFSET_MM for o in offs),
          f"no plane is left sitting on a vessel tip: {offs}")
    check(any("already an opening" in n for n in notes),
          "and it says why, rather than silently cutting one fewer plane")
    # Dropping the proximal plane means everything above it stays, so the branches
    # up there need their own cuts -- membership has to widen with the boundary.
    print(f"  {len(inside)} branch(es) inside the region")
    check(len(inside) > 0, "branches above the dropped boundary are still cut")

print()
print("=" * 78)
print("a region with real boundaries is untouched")
print("=" * 78)
neck = next((p for p in presets if p["name"] == "Aortic neck"), None)
check(neck is not None, "'Aortic neck' is offered")
if neck:
    s1 = route.length if neck["s1"] is None else neck["s1"]
    b, _i, _r, _n = region_plan(tree, neck["s0"], s1, route, jg, mesh=mesh,
                                clear_ostia=True)
    offs = [round(float(getattr(t, "_region_offset", -1)), 2) for t in b]
    print(f"  {len(b)} bound(s), offsets {offs}")
    check(len(b) == 2, "both planes are still placed on a region inside the vessel")

print()
print("=" * 78)
print("and the region still cuts, end to end")
print("=" * 78)
from vascular_api import CutSession                                 # noqa: E402
s = CutSession(str(CSV), str(MAT))
for nm in ("Whole abdominal aorta", "Aortic neck"):
    if nm not in [p["name"] for p in s.presets()]:
        check(False, f"{nm} is offered by the session")
        continue
    s.apply_preset(nm)
    st = s.recut()
    s.clear_region()
    print(f"  {nm:24s} ok={st['ok']}  {st['message'][:44]}")
    check(st["ok"], f"{nm} cuts: {st['message']}")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
