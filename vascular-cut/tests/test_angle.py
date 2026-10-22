"""Neck angulation from an osculating circle: theta = L / R.

Checked against shapes whose answer is known before it is pointed at an aorta.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from branch_cut import centerline_arc_angle

fails = []


class FakeSeg:
    def __init__(self, pts):
        self.pts = np.asarray(pts, float)
        d = np.linalg.norm(np.diff(self.pts, axis=0), axis=1)
        self.cum = np.concatenate([[0.0], np.cumsum(d)])
        self.length = float(self.cum[-1])

    def at_arclength(self, s):
        i = int(np.clip(np.searchsorted(self.cum, s) - 1, 0, len(self.pts) - 2))
        f = ((s - self.cum[i]) / max(self.cum[i + 1] - self.cum[i], 1e-9))
        p = self.pts[i] + f * (self.pts[i + 1] - self.pts[i])
        t = self.pts[i + 1] - self.pts[i]
        return p, t / max(np.linalg.norm(t), 1e-9)


class FakeTree:
    def __init__(self, seg):
        self.segments = {"S": seg}


class FakeRoute:
    def __init__(self, seg):
        self.length = seg.length

    def locate(self, s):
        return "S", s, True


def arc(radius, sweep_deg, n=200, noise=0.0, tilt=0.0):
    a = np.linspace(0, np.radians(sweep_deg), n)
    p = np.column_stack([radius * np.cos(a), radius * np.sin(a), np.zeros_like(a)])
    if tilt:
        p[:, 2] = tilt * np.sin(np.linspace(0, 3 * np.pi, n))
    if noise:
        rng = np.random.default_rng(0)
        p = p + rng.normal(0, noise, p.shape)
    return p


print("=" * 74)
print("known arcs: a circle of radius R over sweep t should return t")
print("=" * 74)
for R, sweep in [(30.0, 60.0), (30.0, 90.0), (50.0, 45.0), (120.0, 20.0)]:
    seg = FakeSeg(arc(R, sweep))
    res = centerline_arc_angle(FakeTree(seg), FakeRoute(seg), 0.0, seg.length)
    err = abs(res["angle_deg"] - sweep)
    print(f"  R={R:5.0f}  sweep={sweep:5.1f}  ->  angle={res['angle_deg']:6.1f}  "
          f"R_fit={res['radius_mm']:6.1f}  rms={res['rms_mm']:.3f}  "
          f"planarity={res['planarity']:.3f}   err {err:.2f} deg")
    if err > 0.6:
        fails.append(f"R={R} sweep={sweep}: off by {err:.2f} deg")
    if abs(res["radius_mm"] - R) > 0.6:
        fails.append(f"R={R}: fitted {res['radius_mm']}")

print()
print("=" * 74)
print("a straight line has no angulation and an enormous radius")
print("=" * 74)
line = np.column_stack([np.linspace(0, 60, 120), np.zeros(120), np.zeros(120)])
seg = FakeSeg(line)
res = centerline_arc_angle(FakeTree(seg), FakeRoute(seg), 0.0, seg.length)
print(f"  angle={res['angle_deg']}  radius={res['radius_mm']}  rms={res['rms_mm']}")
if res["angle_deg"] > 2.0:
    fails.append(f"a straight line returned {res['angle_deg']} deg")

print()
print("=" * 74)
print("the diagnostics have to flag a curve that is not an arc")
print("=" * 74)
seg = FakeSeg(arc(30.0, 90.0, noise=1.2))
res = centerline_arc_angle(FakeTree(seg), FakeRoute(seg), 0.0, seg.length)
print(f"  noisy arc     : angle={res['angle_deg']:6.1f}  rms={res['rms_mm']:.2f}")
seg = FakeSeg(arc(30.0, 90.0, tilt=6.0))
res = centerline_arc_angle(FakeTree(seg), FakeRoute(seg), 0.0, seg.length)
print(f"  out of plane  : angle={res['angle_deg']:6.1f}  rms={res['rms_mm']:.2f}  "
      f"planarity={res['planarity']:.3f}")
if res["planarity"] < 0.05:
    fails.append("a non-planar curve reported planarity near zero")

print()
print("=" * 74)
print("too short to fit")
print("=" * 74)
seg = FakeSeg(arc(30.0, 1.0, n=6))
print(f"  1 degree of arc -> {centerline_arc_angle(FakeTree(seg), FakeRoute(seg), 0.0, seg.length)}")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
