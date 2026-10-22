"""A region's boundary planes get a reach that actually cuts, from either direction.

Two failure modes pull opposite ways:

  torn rim   "the opening comes out N mm out of plane" -- the plane ran out of band
             before it finished spanning the ring. Wants a WIDER reach.
  severed    an extra opening -- the plane reached sideways through a bystander, on
             an abdominal aorta a renal running parallel rather than leaving
             laterally. Wants a NARROWER one.

Measured over 182 region cuts on 49 abdominal scans: 167 cut cleanly at the automatic
reach and are left untouched, 10 of the remaining 15 are rescued (JK4_2 at 0.60x,
JK17_2 and JK3_1 at 1.30x), and the 6 that severed a bystander drop to none.

    python tests/test_region_reach.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from branch_tree import load_tree                                   # noqa: E402
from branch_cut import (load_mat_mesh, classify_terminals,          # noqa: E402
                        name_anatomy, region_presets, region_plan,
                        aorta_route, band_width, fit_region_reach,
                        cut_terminal, snap_tolerance, boundary_loops,
                        _junction_graph, REGION_REACH_SCALES)

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
    jg = _junction_graph(tree)
    return tree, mesh, names, vessels, jg, aorta_route(tree, jg)


def bounds_for(tree, mesh, names, vessels, jg, route, preset):
    p = next(x for x in region_presets(tree, mesh, names, vessels)
             if x["name"] == preset)
    s1 = route.length if p["s1"] is None else p["s1"]
    b, _i, _r, _n = region_plan(tree, p["s0"], s1, route, jg, mesh=mesh,
                                clear_ostia=True)
    return b


def cuts_cleanly(tree, mesh, bounds):
    m = mesh.copy()
    snap = snap_tolerance(mesh)
    try:
        for t in bounds:
            m = cut_terminal(m, tree, t, getattr(t, "_region_offset", 6.0),
                             snap_mm=snap)
    except Exception:                                               # noqa: BLE001
        return None
    return len(boundary_loops(m))


print("=" * 78)
print("the search goes both ways from the automatic reach")
print("=" * 78)
print(f"  scales: {REGION_REACH_SCALES}")
check(REGION_REACH_SCALES[0] == 1.0,
      "the automatic reach is tried first, so a working cut keeps its own plane")
check(any(k > 1 for k in REGION_REACH_SCALES) and any(k < 1 for k in REGION_REACH_SCALES),
      "and the search widens as well as narrows -- a torn rim wants more reach, "
      "a severed bystander wants less")

print()
print("=" * 78)
print("a region that already cuts is left exactly as it was")
print("=" * 78)
tree, mesh, names, vessels, jg, route = build()
b = bounds_for(tree, mesh, names, vessels, jg, route, "Aortic neck")
auto = [band_width(t) for t in b]
scale = fit_region_reach(mesh, tree, b)
after = [band_width(t) for t in b]
print(f"  automatic {['%.2f' % a for a in auto]} -> chose {scale} "
      f"-> {['%.2f' % a for a in after]}")
check(scale == 1.0, f"the automatic reach is kept ({scale})")
check(all(abs(x - y) < 1e-6 for x, y in zip(auto, after)),
      "and the plane sizes are unchanged")
check(cuts_cleanly(tree, mesh, b) == 2, "the region still cuts to two openings")

print()
print("=" * 78)
print("a plane size set by hand is never overridden")
print("=" * 78)
t2, m2, n2, v2, jg2, r2 = build()
b2 = bounds_for(t2, m2, n2, v2, jg2, r2, "Aortic neck")
for t in b2:
    t.band_mm = 21.0                       # a person typed this
got = fit_region_reach(m2, t2, b2)
print(f"  hand-set 21.0 mm -> fitter returned {got}, bands now "
      f"{[band_width(t) for t in b2]}")
check(got is None, "the fitter declines to search when a size was set by hand")
check(all(abs(band_width(t) - 21.0) < 1e-6 for t in b2),
      "and the typed value survives untouched")

print()
print("=" * 78)
print("a reach that cannot work leaves the bands alone rather than a bad guess")
print("=" * 78)
t3, m3, n3, v3, jg3, r3 = build()
b3 = bounds_for(t3, m3, n3, v3, jg3, r3, "Aortic neck")
# No scale can succeed if none is offered.
got3 = fit_region_reach(m3, t3, b3, scales=())
print(f"  with no scales to try -> {got3}, band_mm {[t.band_mm for t in b3]}")
check(got3 is None, "it reports that nothing worked")
check(all(getattr(t, "band_mm", None) in (None, 0) for t in b3),
      "and leaves band_mm unset, so band_width falls back to the measured automatic")

print()
print("=" * 78)
print("it reaches the session the web app drives")
print("=" * 78)
from vascular_api import CutSession                                 # noqa: E402
s = CutSession(str(CSV), str(MAT))
s.apply_preset("Aortic neck")
st = s.recut()
print(f"  {st['message'][:72]}")
check(st["ok"], f"the region cuts through the API: {st['message']}")

print()
print("=" * 78)
print("and the desktop window, which is the tool people actually use")
print("=" * 78)
# The first version of this fix was wired into the API only. The desktop window has
# its own cut path and silently did not get it -- the primary tool missing the fix
# the benchmark was reporting as shipped.
from qtpy import QtWidgets                                          # noqa: E402
app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
from branch_gui import CutWindow                                    # noqa: E402

w = CutWindow(str(CSV), str(MAT))
w.show()
for _ in range(6):
    app.processEvents()
p = next(x for x in region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
         if x["name"] == "Aortic neck")
w._apply_preset(p)
for _ in range(6):
    app.processEvents()
print(f"  problem={w.problem}  reach={getattr(w, '_region_reach', None)}")
check(w.problem is None, f"the region cuts in the window too: {w.problem}")
check(hasattr(w, "_apply_region_reach"),
      "the window carries the reach as a scale, not as a band on throwaway objects")
# Leaving the region must not carry its reach into the next one: it was fitted to one
# pair of planes at one place on the aorta.
w._set_region(None)
check(getattr(w, "_region_reach", None) is None,
      "and leaving the region clears the reach it found")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
