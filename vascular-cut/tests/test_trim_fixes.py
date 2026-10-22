"""Verify the review fixes: no more false-green trims, offsets survive rebuilds,
smoothing keeps clearance, candidates are viable."""
import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import (load_mat_mesh, classify_terminals, pre_bifurcation_candidates,
                        cut_terminal, keep_largest_component, boundary_loops,
                        validate_cut, snap_tolerance, vessel_coords, segment_radius)
from branch_tree import load_tree, Terminal

app = QtWidgets.QApplication(sys.argv)
fails = []

# ---------------------------------------------------------------- 1. false-green
print("=" * 80)
print("1. FALSE-GREEN: a trim must not leave subtree attached while reporting Ready")
print("=" * 80)
tree = load_tree(ROOT / "JK5_1_CentrelinePoints_branches.csv", smooth_mm=1.0)
mesh = load_mat_mesh(ROOT / "JK5_1_M10.mat")
classify_terminals(tree, mesh, verbose=False)
snap = snap_tolerance(mesh)
cands = pre_bifurcation_candidates(tree, mesh)
print(f"  candidates offered: {[c['segment'] for c in cands]}")

print(f"  {'trim':8s}{'off':>6}{'valid':>8}{'subtree left (mm2)':>20}   verdict")
for c in cands:
    seg = tree.segments[c["segment"]]
    for off in [2.0, 8.0, 15.0]:
        lim = max(seg.length - 2.0 * (seg._radius_mm or 5.0), 2.0)
        if off > lim:
            continue
        t = Terminal(seg, c["at_end"], label=f"{seg.name}.trim", mode="normal", name="trim")
        t.is_interior = True; t.subtree_segments = set(c["subtree"])
        gone = {x.label for x in c["terminals"]}
        active = [x for x in tree.terminals if x.label not in gone] + [t]
        offs = {x.label: 6.0 for x in active}; offs[t.label] = off
        try:
            m = mesh.copy()
            for x in active:
                m = cut_terminal(m, tree, x, offs[x.label], snap_mm=snap)
            m = keep_largest_component(m)
        except Exception as e:
            print(f"  {seg.name:8s}{off:6.1f}{'error':>8}{'-':>20}   {str(e)[:28]}")
            continue
        prob = validate_cut(m, tree, offs, terminals=active)
        # What matters is not how much subtree material remains (a collar right at the
        # ostium is legitimately owned by the child vessels) but how FAR it reaches:
        # a surviving blind-ended branch sticks tens of mm past the plane.
        own, _ = vessel_coords(m, tree)
        keep = np.isin(own, list(c["subtree"]))
        pt = t.cut_plane(off)[0]
        reach = (float(np.linalg.norm(
            np.asarray(m.triangles_center)[keep] - pt, axis=1).max()) if keep.any() else 0.0)
        rr = max(float(seg._radius_mm or 5.0), 1.0)
        green = prob is None
        bad = green and reach > 3.0 * rr      # a real stub, not an ostium collar
        print(f"  {seg.name:8s}{off:6.1f}{('GREEN' if green else 'red'):>8}"
              f"{reach:12.1f} (3r={3*rr:.0f})   "
              f"{'FALSE-GREEN!' if bad else ('ok' if green else 'refused')}")
        if bad:
            fails.append(("false-green", seg.name, off, round(reach, 1)))

# ---------------------------------------------------------------- 2. GUI state
print()
print("=" * 80)
print("2. GUI: offsets + tilt survive a trim rebuild; colour cache refreshes")
print("=" * 80)
w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
for _ in range(10): app.processEvents()

r0 = w.rows[1]
r0.set_value(17.0); r0.terminal.tilt_a_deg = 12.0
w._on_change()
for _ in range(3): app.processEvents()
label0, before = r0.terminal.label, r0.value()
c15 = next(c for c in pre_bifurcation_candidates(w.tree, w.base) if c["segment"] == "Line15")
w._add_trim(c15)
for _ in range(6): app.processEvents()
after_row = next(r for r in w.rows if r.terminal.label == label0)
print(f"  {label0}: offset {before:.1f} -> {after_row.value():.1f} mm after adding a trim")
print(f"  tilt widget {after_row.tilt_a_b.value():.1f}° vs terminal {after_row.terminal.tilt_a_deg:.1f}°")
if abs(after_row.value() - before) > 0.05:
    fails.append(("offset-reset", before, after_row.value()))
if abs(after_row.tilt_a_b.value() - after_row.terminal.tilt_a_deg) > 0.05:
    fails.append(("tilt-desync",))

w.colour_cb.setChecked(True)
for _ in range(3): app.processEvents()
rgb_a = w._opening_face_rgb(w.cut)
w._remove_trim(next(r.terminal for r in w.rows if getattr(r.terminal, "is_interior", False)))
for _ in range(6): app.processEvents()
rgb_b = w._opening_face_rgb(w.cut)
changed = rgb_a is None or rgb_b is None or not np.array_equal(rgb_a, rgb_b)
print(f"  colour map refreshed after undoing the trim: {changed}")
if not changed:
    fails.append("colour-cache-stale")
print(f"  status: {w.status.text()[:64]}")
if w.problem is not None:
    fails.append(("gui-invalid", w.problem[:50]))
w.close()

# ---------------------------------------------------------------- 3. smoothing
print()
print("=" * 80)
print("3. SMOOTHING: topology intact AND the line keeps its clearance from the wall")
print("=" * 80)
import trimesh
for name, csv, mat in [("JK5", "JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
                       ("KY160", "KY160_1_CentrelinePoints.csv",
                        "Geometries & Meshes/mesh/KY160_1.mat")]:
    raw = load_tree(ROOT / csv, smooth_mm=0.0)
    msh = load_mat_mesh(ROOT / mat)
    sm = load_tree(ROOT / csv, smooth_mm=1.0)
    topo = (len(sm.junctions) == len(raw.junctions)
            and len(sm.terminals) == len(raw.terminals))
    def clearance(tr):
        worst = np.inf
        for s in tr.segments.values():
            d = trimesh.proximity.signed_distance(msh, s.pts[::5])
            worst = min(worst, float(d.min()))
        return worst
    c_raw, c_sm = clearance(raw), clearance(sm)
    turn_raw = raw.worst_turning
    turn_sm = max(s.turning for s in sm.segments.values())
    # Judge against the RAW line, not against zero: these centerlines already graze
    # the wall at their tips as supplied, so the question is whether smoothing makes
    # it worse (the un-windowed pin cost ~1.2 mm on JK5).
    print(f"  {name:6s} topology {'kept' if topo else 'BROKEN'}   "
          f"min wall clearance raw {c_raw:6.2f} -> smoothed {c_sm:6.2f} mm "
          f"(Δ {c_sm - c_raw:+.2f})   turning {turn_raw:.1f} -> {turn_sm:.1f} deg/mm")
    if not topo:
        fails.append((name, "topology"))
    if c_sm < c_raw - 0.25:
        fails.append((name, "smoothing-erodes-clearance", round(c_sm - c_raw, 2)))

print()
print("=" * 80)
print("ALL PASS" if not fails else f"FAILURES: {fails}")
print("=" * 80)

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
