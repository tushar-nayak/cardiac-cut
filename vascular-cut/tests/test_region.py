"""Can we isolate a stretch of aorta -- the neck, a landing zone, the arch?"""
import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from branch_tree import load_tree
from branch_cut import (load_mat_mesh, classify_terminals, name_anatomy, aorta_route,
                        aorta_landmarks, region_terminals, cut_terminal,
                        keep_largest_component, boundary_loops, snap_tolerance,
                        section_diameter, _junction_graph)

fails = []


def isolate(tree, mesh, s0, s1, branch_offset=6.0, verbose=True):
    """Cut the model down to the stretch [s0, s1] of the aorta, branches kept."""
    from branch_cut import region_plan
    bounds, inside, removed, notes = region_plan(tree, s0, s1)
    if verbose and notes:
        for n in notes: print(f"      note: {n}")
    if len(bounds) != 2:
        return None, None
    active = inside + bounds
    m = mesh.copy()
    snap = snap_tolerance(mesh)
    for t in active:
        off = getattr(t, "_region_offset", branch_offset)
        try:
            m = cut_terminal(m, tree, t, off, snap_mm=snap)
        except Exception as e:
            if verbose:
                print(f"      cut {t.name} failed: {str(e)[:52]}")
            return None, None
    return keep_largest_component(m), active


# ---------------------------------------------------------------- abdominal neck
print("=" * 82)
print("JK5, isolate the AORTIC NECK (lowest renal -> where the sac begins)")
print("=" * 82)
tree = load_tree(ROOT / "JK5_1_CentrelinePoints_branches.csv")
mesh = load_mat_mesh(ROOT / "JK5_1_M10.mat")
classify_terminals(tree, mesh, verbose=False)
names, vessels, _ = name_anatomy(tree, mesh)
route = aorta_route(tree)
land = aorta_landmarks(tree, route, names)
print(f"  route length {route.length:.1f} mm over {len(route.chain)} segments")
for k, v in sorted(land.items(), key=lambda kv: kv[1]):
    print(f"    {v:7.1f} mm  {k}")

lo_renal = max(v for k, v in land.items() if "renal" in k.lower())
# distal end = where the aorta first exceeds the reference calibre by 20%
ref = section_diameter(mesh, route.point(lo_renal + 1),
                       tree.segments[route.locate(lo_renal + 1)[0]].at_arclength(
                           route.locate(lo_renal + 1)[1])[1])
sac = None
for s in np.arange(lo_renal + 2, route.length, 2.0):
    nm, loc, _ = route.locate(s)
    pos, tan = tree.segments[nm].at_arclength(loc)
    sd = section_diameter(mesh, pos, tan, near_point=pos)
    if sd and sd["ead_mm"] > 1.2 * ref["ead_mm"]:
        sac = s
        break
print(f"\n  lowest renal at {lo_renal:.1f} mm (D {ref['ead_mm']:.1f} mm); "
      f"dilates past +20% at {sac:.1f} mm  ->  neck {sac - lo_renal:.1f} mm")

# The neck begins just BELOW the lowest renal: a plane sitting exactly on the ostium
# clips the renal's own mouth and leaves a partial hole in the aorta wall.
CLEAR_MM = 3.0
out, active = isolate(tree, mesh, lo_renal, sac)
pass
if out is None:
    fails.append("neck-cut-failed")
else:
    loops = boundary_loops(out)
    print(f"  isolated neck: {len(out.faces)} faces, {len(loops)} openings, "
          f"euler {int(out.euler_number)}, {len(active)} cuts "
          f"({len(active)-2} branch(es) inside)")
    ext = np.asarray(out.bounds[1]) - np.asarray(out.bounds[0])
    print(f"  bounding box {np.round(ext,1)} mm")
    if len(loops) != len(active):
        fails.append(("neck-openings", len(loops), len(active)))

# ---------------------------------------------------------------- thoracic arch
print()
print("=" * 82)
print("KY160, isolate the AORTIC ARCH (innominate -> 20 mm past the last vessel)")
print("=" * 82)
t2 = load_tree(ROOT / "KY160_1_CentrelinePoints.csv")
m2 = load_mat_mesh(ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat")
classify_terminals(t2, m2, verbose=False)
n2, v2, _ = name_anatomy(t2, m2)
r2 = aorta_route(t2)
l2 = aorta_landmarks(t2, r2, n2)
print(f"  route length {r2.length:.1f} mm over {len(r2.chain)} segments")
for k, v in sorted(l2.items(), key=lambda kv: kv[1]):
    print(f"    {v:7.1f} mm  {k}")
# the arch: from just before the first great vessel to 20 mm past the last
first, last = min(l2.values()), max(l2.values())
s0, s1 = max(first - 10.0, 1.0), min(last + 20.0, r2.length - 1.0)
print(f"\n  arch region: {s0:.1f} -> {s1:.1f} mm ({s1-s0:.1f} mm of aorta)")
out2, act2 = isolate(t2, m2, s0, s1)
if out2 is None:
    fails.append("arch-cut-failed")
else:
    loops = boundary_loops(out2)
    print(f"  isolated arch: {len(out2.faces)} faces, {len(loops)} openings, "
          f"euler {int(out2.euler_number)}, {len(act2)} cuts "
          f"({len(act2)-2} branch(es) kept)")
    ext = np.asarray(out2.bounds[1]) - np.asarray(out2.bounds[0])
    print(f"  bounding box {np.round(ext,1)} mm")
    if len(loops) != len(act2):
        fails.append(("arch-openings", len(loops), len(act2)))

print()
print("=" * 82)
print("ALL PASS" if not fails else f"FAILURES: {fails}")
print("=" * 82)

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
