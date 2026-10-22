"""segment_radius must measure the LUMEN, not the model's end cap.

THE BUG. segment_radius takes the median distance from the centerline to the surface
over the middle 20-80% of a segment. Its docstring says the middle is sampled to dodge
the cap at a terminal -- but on a SHORT segment the middle is still beside the cap, so
the dodge fails and the number returned is the distance to the end of the model.

A wide, short vessel is the worst case and it is exactly the shape of an aortic inlet
segment. A 15 mm-radius tube only 20 mm long: every station in the middle is nearer to
one of its two caps than to its own wall, so the answer comes back at roughly a third
of the truth.

WHAT IT COST. On JK1_2's Line06 -- 17.8 mm long, aorta ~15 mm across -- it read
8.586 mm. The same centerline, the same 25 stations, measured again after the inlet cap
had been cut away, reads 14.673 mm. 8.586 is the value the whole tool used: it feeds
trunk-versus-branch ranking, anatomical naming, band_width, and the normalisation that
decides which vessel owns a face. Across the abdominal cohort it is 39 of 891 segments,
always about one per case, always the short segment carrying the inlet.

Everything here is synthetic on purpose. The fault is a property of the measurement,
not of any scan, so it should be provable without a cohort.
"""
import sys
from pathlib import Path
import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from branch_tree import Segment                                    # noqa: E402
import branch_cut as bc                                            # noqa: E402

fails = []
print("=" * 78)
print("segment_radius against a capped tube")
print("=" * 78)


def axis_segment(name, z0, z1, npts=41):
    z = np.linspace(z0, z1, npts)
    pts = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    return Segment(name, pts, np.tile([0.0, 0.0, 1.0], (npts, 1)))


# ---- 1. short and wide: the shape of an inlet segment -----------------------
#
# trimesh's cylinder is centred on the origin and closed at both ends, so this is a
# vessel whose every mid-station is closer to a cap than to its own wall.
R, H = 15.0, 20.0
tube = trimesh.creation.cylinder(radius=R, height=H, sections=64)
seg = axis_segment("short", -H / 2, H / 2)
got = bc.segment_radius(tube, seg)
print(f"  tube r={R} h={H}  (every mid-station is nearer a cap than the wall)")
print(f"    segment_radius -> {got:6.3f}   truth {R:.3f}")
if abs(got - R) > 0.15 * R:
    fails.append(f"short wide tube: segment_radius {got:.3f}, truth {R:.3f}")

# ---- 2. long and thin: nothing may change -----------------------------------
#
# The ordinary case, where the middle really is far from both caps. This is the guard
# against "fixing" the short case by changing the answer everywhere.
r2, h2 = 3.0, 60.0
thin = trimesh.creation.cylinder(radius=r2, height=h2, sections=48)
seg2 = axis_segment("long", -h2 / 2, h2 / 2)
got2 = bc.segment_radius(thin, seg2)
print(f"\n  tube r={r2} h={h2}  (the ordinary case, must be untouched)")
print(f"    segment_radius -> {got2:6.3f}   truth {r2:.3f}")
if abs(got2 - r2) > 0.1 * r2:
    fails.append(f"long thin tube: segment_radius {got2:.3f}, truth {r2:.3f}")

# ---- 3. a centerline OUTSIDE the mesh must still read negative --------------
#
# A negative radius is not a failure to measure, it is the signal that the centerline
# does not lie inside this surface -- mismatched() in cut_bench and _no_radius_reason()
# both key off it, and cut_bench REJECTS a case when a third of its segments read
# negative. Rays fired from outside a surface come back with a confident positive
# number, so the fallback must never run here.
outside = axis_segment("outside", -h2 / 4, h2 / 4)
outside.pts[:, 0] += r2 * 4.0            # shift the whole line clear of the tube
outside._kd_own = None
got3 = bc.segment_radius(thin, outside)
print(f"\n  the same tube, centerline moved {r2*4:.0f} mm outside it")
print(f"    segment_radius -> {got3:6.3f}   (must stay negative)")
if not (got3 < 0.0):
    fails.append(f"a centerline outside the mesh read {got3:.3f}, not negative: "
                 f"the mesh/centerline mismatch signal is gone")

print()
if fails:
    print("FAILURES:")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
