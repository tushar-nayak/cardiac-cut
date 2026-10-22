"""A region boundary's plane size must be measured, and must not depend on call order.

A region bound sits on an INTERIOR aortic segment. classify_terminals only caches a
radius for segments that carry terminals, so band_width found none there and fell
through to its own r = 5.0 default -- handing every region boundary in the abdominal
cohort a 12.50 mm reach whatever the aorta measured. True radii on those bounds run
14.1 to 25.6 mm.

Worse than wrong, it was order dependent. pre_bifurcation_candidates(tree, mesh)
fills that same cache as a side effect, so a session that happened to open the trim
menu first cut the identical region with a 35 mm band instead of a 12.5 mm one. Same
input, same request, two different cuts, and nothing said so.

    python tests/test_region_band.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from branch_tree import load_tree                                   # noqa: E402
from branch_cut import (load_mat_mesh, classify_terminals,          # noqa: E402
                        name_anatomy, region_presets, region_plan,
                        aorta_route, band_width, segment_radius,
                        pre_bifurcation_candidates, _junction_graph)

CSV = ROOT / "JK5_1_CentrelinePoints_branches.csv"
MAT = ROOT / "JK5_1_M10.mat"

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


def build():
    tree = load_tree(str(CSV))
    mesh = load_mat_mesh(str(MAT))
    classify_terminals(tree, mesh, verbose=False)
    names, vessels, _ = name_anatomy(tree, mesh)
    return tree, mesh, names, vessels


def bounds_for(tree, mesh, names, vessels, preset="Aortic neck"):
    jg = _junction_graph(tree)
    route = aorta_route(tree, jg)
    p = next(x for x in region_presets(tree, mesh, names, vessels)
             if x["name"] == preset)
    s1 = route.length if p["s1"] is None else p["s1"]
    b, _i, _r, _n = region_plan(tree, p["s0"], s1, route, jg, mesh=mesh,
                                clear_ostia=True)
    return b


print("=" * 78)
print("the band is measured from the aorta, not from a default")
print("=" * 78)
tree, mesh, names, vessels = build()
bounds = bounds_for(tree, mesh, names, vessels)
check(len(bounds) == 2, f"the region has two bounds ({len(bounds)})")
for t in bounds:
    true_r = segment_radius(mesh, t.segment)
    band = band_width(t)
    print(f"  {t.name:14s} seg={t.segment.name:10s} radius={true_r:6.2f} mm  "
          f"band={band:6.2f} mm  (2.5x radius = {2.5 * true_r:.2f})")
    check(abs(band - 2.5 * true_r) < 0.5,
          f"{t.name}: the band is 2.5x the measured radius, not 2.5x a 5 mm default")
    # 12.50 is exactly 2.5 * 5.0, the fallback. Anything landing there on a 15 mm
    # aorta is the bug, not a coincidence.
    check(abs(band - 12.5) > 0.01 or abs(true_r - 5.0) < 0.5,
          f"{t.name}: not the 12.50 mm fallback")

print()
print("=" * 78)
print("and it does not depend on what ran first")
print("=" * 78)
# Fresh tree, but with the trim menu opened first -- which fills the same radius
# cache as a side effect. Before the fix this alone changed the plane size.
t2, m2, n2, v2 = build()
pre_bifurcation_candidates(t2, m2)          # the side effect
b2 = bounds_for(t2, m2, n2, v2)
cold = [round(band_width(t), 3) for t in bounds]
warm = [round(band_width(t), 3) for t in b2]
print(f"  cold (region asked for first): {cold}")
print(f"  warm (trim menu opened first): {warm}")
check(cold == warm,
      f"the same region gets the same plane size either way ({cold} vs {warm})")

print()
print("=" * 78)
print("the region still cuts")
print("=" * 78)
from branch_cut import cut_terminal, snap_tolerance, boundary_loops   # noqa: E402
m = mesh.copy()
snap = snap_tolerance(mesh)
err = ""
for t in bounds:
    try:
        m = cut_terminal(m, tree, t, getattr(t, "_region_offset", 6.0), snap_mm=snap)
    except Exception as e:                                          # noqa: BLE001
        err = str(e)[:70]
        break
loops = len(boundary_loops(m))
print(f"  loops={loops} euler={int(m.euler_number)} {err}")
check(not err, f"both boundary planes cut: {err}")
# Two planes through a closed surface leave exactly two openings. More means one of
# them reached sideways through a vessel it was not aimed at.
check(loops == 2, f"and leave exactly two openings, not {loops}")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
