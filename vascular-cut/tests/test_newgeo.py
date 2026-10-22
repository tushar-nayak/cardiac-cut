"""Test section_diameter and centerline_geodesic."""
import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from branch_tree import load_tree
from branch_cut import (load_mat_mesh, classify_terminals, section_diameter,
                        centerline_geodesic, cut_terminal, keep_largest_component,
                        boundary_loops, match_openings_to_terminals)

tree = load_tree(ROOT / "KY160_1_CentrelinePoints.csv")
mesh = load_mat_mesh(ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat")
classify_terminals(tree, mesh, verbose=False)
fails = []

print("=" * 78)
print("A. section_diameter, sanity + comparison to the actual cut opening")
print("=" * 78)
print(f"  {'terminal':<22}{'off':>5}{'EAD':>8}{'dmin':>7}{'dmax':>7}{'circ':>7}"
      f"{'  cut-opening dia':>18}")
for t in tree.terminals:
    for off in [6.0, 20.0]:
        p, n = t.cut_plane(off)
        sd = section_diameter(mesh, p, n, near_point=p)
        # ground-truth-ish: the actual cut opening diameter at the same plane
        m = mesh.copy()
        try:
            for x in tree.terminals:
                m = cut_terminal(m, tree, x, off if x is t else 6.0)
            m = keep_largest_component(m)
            mt = match_openings_to_terminals(m, tree,
                    {x.label: (off if x is t else 6.0) for x in tree.terminals})
            open_dia = mt[t.label]["diameter"] if mt.get(t.label) else float("nan")
        except Exception:
            open_dia = float("nan")
        if sd is None:
            print(f"  {t.name:<22}{off:5.0f}   (no section)")
            continue
        print(f"  {t.name:<22}{off:5.0f}{sd['ead_mm']:8.1f}{sd['dmin_mm']:7.1f}"
              f"{sd['dmax_mm']:7.1f}{sd['circularity']:7.2f}{open_dia:15.1f} mm")
        # EAD must lie between dmin and dmax, and be positive
        if not (sd['dmin_mm'] - 0.5 <= sd['ead_mm'] <= sd['dmax_mm'] + 0.5):
            fails.append(("ead-bounds", t.name, off, sd['ead_mm'], sd['dmin_mm'], sd['dmax_mm']))

print()
print("  A circle check: section a synthetic cylinder -> EAD == dmin == dmax")
import trimesh
cyl = trimesh.creation.cylinder(radius=10, height=100, sections=200)
sd = section_diameter(cyl, [0, 0, 0], [0, 0, 1])
print(f"    cylinder r=10: EAD={sd['ead_mm']:.2f} dmin={sd['dmin_mm']:.2f} "
      f"dmax={sd['dmax_mm']:.2f} circ={sd['circularity']:.3f}  (expect ~20, ~1.0)")
if abs(sd['ead_mm'] - 20) > 0.5 or abs(sd['circularity'] - 1.0) > 0.02:
    fails.append(("cylinder", sd))

print()
print("=" * 78)
print("B. centerline_geodesic")
print("=" * 78)
print("  same-segment: distance between two points on Line7 == |s1 - s2|")
seg = tree.segments["Line7"]
for s1, s2 in [(10, 50), (5, 300), (100, 100.0001)]:
    pa = seg.at_arclength(s1)[0]
    pb = seg.at_arclength(s2)[0]
    g = centerline_geodesic(tree, pa, pb)
    expect = abs(s1 - s2)
    ok = g["same_segment"] and abs(g["distance_mm"] - expect) < 0.5
    print(f"    s1={s1:6.1f} s2={s2:6.1f} -> {g['distance_mm']:7.2f} mm "
          f"(expect {expect:6.2f})  same_seg={g['same_segment']}  {'OK' if ok else 'FAIL'}")
    if not ok:
        fails.append(("same-seg", s1, s2, g["distance_mm"], expect))

print()
print("  cross-segment: a point on branch Line4 to a point on the descending aorta")
print("    must go via the junctions and be >= straight-line distance")
pa = tree.segments["Line4"].at_arclength(50)[0]     # on an arch branch
pb = tree.segments["Line7"].at_arclength(200)[0]    # on the descending aorta
g = centerline_geodesic(tree, pa, pb)
straight = np.linalg.norm(pa - pb)
# lower bound: Line4 residual (50 to J1) + ... just check it's plausible & >= straight
print(f"    geodesic={g['distance_mm']:.1f} mm   straight-line={straight:.1f} mm   "
      f"same_seg={g['same_segment']}")
ok = (not g["same_segment"]) and g["distance_mm"] >= straight - 1
print(f"    {'OK' if ok else 'FAIL'} (geodesic >= straight-line)")
if not ok:
    fails.append(("cross-seg", g["distance_mm"], straight))

# cross-seg ground truth: Line4(50->J1) + Line1(J1->J2) + Line2(J2->J3) + Line7(J3->200)
# J1 is Line4.start so residual=50; Line1 full=29.31; Line2 full=38.55; Line7 200 from start(tip) to J3(end): L-200
L7 = tree.segments["Line7"].length
gt = 50 + tree.segments["Line1"].length + tree.segments["Line2"].length + (L7 - 200)
print(f"    hand-computed path length = 50 + Line1 + Line2 + (L7-200) = {gt:.1f} mm")
print(f"    match: {abs(g['distance_mm'] - gt) < 1.0}")
if abs(g["distance_mm"] - gt) > 1.0:
    fails.append(("cross-seg-gt", g["distance_mm"], gt))

print()
print("=" * 78)
print("ALL PASS" if not fails else f"FAILURES: {fails}")
print("=" * 78)

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
