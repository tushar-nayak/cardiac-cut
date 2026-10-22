"""A segment may not own wall past the end of its own centerline. argv[1]=jk5|ky160|kk2

THE BUG. vessel_coords settles a doubtful face by asking which vessel's WALL it sits
nearest -- |distance to that centerline - that vessel's radius there|. The distance
comes from a KD-tree over the centerline's POINTS, and a KD-tree clamps: ask it about
a point past the end of the line and it answers with the distance to the last point,
which is mostly a distance ALONG the vessel rather than across it. Compared against a
radius, that asks whether the face sits on a wall the segment does not have there,
because the segment stopped.

WHAT IT COST. On JK1_2 the aorta is split at z=517.6, Line09 below and Line06 above.
All 86 faces of Line06's inlet stub sit 16 mm ABOVE Line09's last centerline point and
clamp to it: 18.74 mm away against a 15.67 mm radius, a miss of 3.06, which beats the
5.11 mm miss against Line06 whose stub they actually are. So Line09 took a quarter of
the surface the inlet cut has to remove, the cutter would not remove another vessel's
wall, and the leftover tongue held the rim 25-35% off plane at EVERY offset, plane size
and mode -- the tool refusing a cut that a plain plane makes flat to 0.000 mm. It is
not a JK1_2 quirk: every one of the 49 abdominal scans has surface claimed past some
endpoint, a median of 5.9% and up to 23.8%.

WHAT IS ASSERTED. Not the JK1_2 numbers, which need the cohort. The invariant behind
them, on a model that ships with the repo: a face taken by the wall contest must lie
where the winning segment actually RUNS.

Two exemptions, both load bearing, both found by breaking them.

  The PRIMARY assignment. A model's end cap genuinely belongs to the segment that
    ends there, and nothing else is nearer, so nearest-centerline is right about it.
    Only the contest -- the step that takes a face AWAY from its nearest centerline
    -- is held to this.

  A FREE end. An endpoint means two unrelated things. At a junction the vessel stops
    and another carries on. At a free end it does not stop at all: the centerline
    ran out before the mesh did, and the cap past it is this vessel's own to claim.
    Gating both alike is byte-identical on cut_bench over 4,145 cuts -- past a free
    end there is only the cap, and nearest-centerline has already awarded it -- so
    this is asserted as the correct statement of when a segment has ended, not as a
    measured improvement.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from vascular_api import CutSession                                # noqa: E402
import branch_cut as bc                                            # noqa: E402

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat"),
          "kk2": ("KK2_5.csv", "KK2_5.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]
fails = []

s = CutSession(str(ROOT / csv), str(ROOT / mat))
mesh, tree = s.mesh, s.tree
print("=" * 78)
print(f"{which.upper()} , {len(mesh.faces)} faces, {len(tree.segments)} segments")
print("=" * 78)

# ---- 1. the unit: distance across the vessel, and how far past its end -------
#
# A straight 40 mm segment up +z. A point beside it is exactly its radial distance
# away with no overshoot; the same point lifted 15 mm past the tip is still 10 mm
# across, and 15 mm past.
from branch_tree import Segment                                    # noqa: E402
pts = np.array([[0.0, 0.0, z] for z in np.linspace(0.0, 40.0, 41)])
straight = Segment("S", pts, np.tile([0.0, 0.0, 1.0], (41, 1)))
probe = np.array([[10.0, 0.0, 20.0],        # beside it
                  [10.0, 0.0, 55.0],        # 15 mm past the tip
                  [10.0, 0.0, -15.0]])      # 15 mm before the start
if not hasattr(bc, "_wall_distance"):
    fails.append("branch_cut has no _wall_distance: the clamped KD distance is "
                 "still being compared against a radius")
else:
    dd, ii = bc._segment_kd(straight).query(probe)
    across, over = bc._wall_distance(straight, probe, dd, ii)
    print(f"  spacing {bc._point_spacing(straight):.2f} mm")
    for p, d, a, o in zip(probe, dd, across, over):
        print(f"    z={p[2]:6.1f}  clamped distance {d:6.2f}  across {a:6.2f}  "
              f"past the end {o:6.2f}")
    for want_a, want_o, got_a, got_o in ((10, 0, across[0], over[0]),
                                         (10, 15, across[1], over[1]),
                                         (10, 15, across[2], over[2])):
        if abs(got_a - want_a) > 1e-6 or abs(got_o - want_o) > 1e-6:
            fails.append(f"_wall_distance gave across={got_a:.2f} past={got_o:.2f}, "
                         f"expected {want_a} / {want_o}")

    # A FREE end must come back untouched: the vessel runs on past where its
    # centerline stopped, and the cap there is its own to claim.
    free_a, free_o = bc._wall_distance(straight, probe, dd, ii, (False, False))
    if not np.allclose(free_a, dd) or float(np.abs(free_o).max()) != 0.0:
        fails.append(f"a free end was gated: across={free_a} past={free_o}")
    print(f"    free end -> distance unchanged {np.allclose(free_a, dd)}, "
          f"no overshoot recorded {float(np.abs(free_o).max()) == 0.0}")

# ---- 2. the invariant, on real geometry -------------------------------------
kd, names, arcs = bc._centerline_kd(tree)
centers = np.asarray(mesh.triangles_center)
_d, idx = kd.query(centers, workers=-1)
primary = names[idx]                        # nearest centerline point, before any contest
owner, _s, _r, _ro = bc.vessel_attribution(mesh, tree)

stolen = np.flatnonzero(owner != primary)
print(f"\n  the wall contest moved {len(stolen)} of {len(centers)} faces "
      f"({100*len(stolen)/len(centers):.1f}%)")

# Measured HERE rather than through branch_cut's own helper, so this section still
# reports the real overshoots when the helper does not exist yet.
#
# Only JUNCTION ends are asserted on. At a free end the vessel does not stop -- the
# centerline runs out before the mesh does, and the cap past it is that vessel's own
# surface to reclaim, so a claim past a free end is not a violation of anything.
nodes = bc._segment_nodes(tree)


def junction_ends(name):
    nd = nodes.get(name) or {}
    return (nd.get("start") is not None, nd.get("end") is not None)


def overshoot(seg, q, ii):
    P = np.asarray(seg.pts, float)
    out = np.zeros(len(q))
    stops = junction_ends(seg.name)
    for idx, nb, stop in ((0, 1, stops[0]),
                          (len(P) - 1, len(P) - 2, stops[1])):
        if not stop:
            continue
        at = np.flatnonzero(ii == idx)
        if not len(at):
            continue
        tang = P[idx] - P[nb]
        tang = tang / max(float(np.linalg.norm(tang)), 1e-9)
        out[at] = np.maximum((q[at] - P[idx]) @ tang, 0.0)
    return out


def spacing(seg):
    P = np.asarray(seg.pts, float)
    return float(np.median(np.linalg.norm(np.diff(P, axis=0), axis=1)))


worst = 0.0
by_seg = {}
for name in np.unique(owner[stolen]):
    seg = tree.segments[name]
    at = stolen[owner[stolen] == name]
    dd, ii = bc._segment_kd(seg).query(centers[at])
    over = overshoot(seg, centers[at], ii)
    tol = spacing(seg)
    by_seg[name] = (len(at), float(over.max()), tol)
    worst = max(worst, float(over.max()) - tol)
    if over.max() > tol + 1e-6:
        fails.append(f"{name} won {(over > tol).sum()} face(s) up to "
                     f"{over.max():.2f} mm past its own junction end (spacing {tol:.2f} mm)")

for name, (n, mx, tol) in sorted(by_seg.items(), key=lambda kv: -kv[1][1])[:6]:
    flag = "  <-- PAST ITS END" if mx > tol + 1e-6 else ""
    print(f"    {name:10s} took {n:5d} face(s), furthest {mx:5.2f} mm past a junction end "
          f"(tolerance {tol:.2f}){flag}")
print(f"\n  worst overshoot past a junction end, beyond tolerance: {worst:+.2f} mm")

# ---- 3. the contest must still DO something ---------------------------------
#
# A gate that refuses everything would pass the invariant above and silently undo the
# fix the contest exists for: a fat trunk's wall claimed by a thin branch running
# alongside it, which is what tore openings nowhere near a junction.
if len(stolen) == 0:
    fails.append("the wall contest no longer moves any face -- the gate refuses all")

print()
if fails:
    print("FAILURES:")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
