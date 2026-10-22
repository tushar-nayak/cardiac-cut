"""Does an EMPTY region export the original closed mesh under a green 'Ready'?

The adversarial review claims: when region_plan yields no bounds, cut_list() is
empty, recut() cuts nothing, self.cut becomes the untouched watertight base, and
validate_cut passes because it expects 0 openings -- so export() writes a CLOSED
mesh named '<stem>_open.stl'. That would be silently wrong geometry.
argv[1] = jk5|ky160
"""
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import boundary_loops, region_plan

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]

app = QtWidgets.QApplication(sys.argv[:1])
# export() pops a MODAL warning when the cut is invalid, which blocks forever in a
# headless run. Swallow it and record that it fired.
_warned = []
QtWidgets.QMessageBox.warning = staticmethod(
    lambda *a, **k: _warned.append(a[2] if len(a) > 2 else ""))

w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(15):
    app.processEvents()

L = w.route.length
print("=" * 78)
print(f"{which.upper()} ,  route {L:.1f} mm; hunting for a region with NO bounds")
print("=" * 78)

# find a span where region_plan returns fewer than 2 bounds
victim = None
for s0 in [x * L / 40.0 for x in range(40)]:
    for width in (0.2, 0.5, 1.0):
        s1 = s0 + width
        if s1 >= L:
            continue
        try:
            bounds, inside, _removed, notes = region_plan(w.tree, s0, s1, w.route)
        except Exception:
            continue
        if len(bounds) != 2:
            victim = (s0, s1, len(bounds), notes)
            break
    if victim:
        break

if victim is None:
    print("  no empty-bound span found in the sweep, claim not reproducible here")
    print("\nNO REPRO")
    sys.exit(0)

s0, s1, nb, notes = victim
print(f"  span {s0:.2f} -> {s1:.2f} mm gives {nb} bound(s); notes={notes}")

w._set_region((s0, s1), "degenerate")
for _ in range(10):
    app.processEvents()

n_loops = len(boundary_loops(w.cut))
print(f"\n  after _set_region:")
print(f"    region bounds : {len(w._region_bounds)}")
print(f"    rows          : {len(w.rows)}")
print(f"    cut_list()    : {len(w.cut_list()[0])} terminals")
print(f"    problem       : {w.problem!r}")
print(f"    status text   : {w.status.text()[:72]!r}")
print(f"    openings      : {n_loops}")
print(f"    watertight    : {w.cut.is_watertight}")
print(f"    faces vs base : {len(w.cut.faces)} vs {len(w.base.faces)}")

identical = len(w.cut.faces) == len(w.base.faces) and w.cut.is_watertight
green = w.problem is None
print(f"\n  VERDICT: closed mesh presented as a cut = {identical and green}")

# would it actually export?
out = Path(__file__).parent / "out_empty"
w.out_dir = out
w.name_edit.setText(f"empty_{which}")
w.export()
for _ in range(6):
    app.processEvents()
stl = out / f"empty_{which}.stl"
side = out / f"empty_{which}_cut_planes.json"
print(f"    export blocked : {not stl.exists()}  "
      f"(warned user: {bool(_warned)})")
if stl.exists():
    meta = json.loads(side.read_text()) if side.exists() else {}
    print(f"    !! WROTE {stl.name} with n_openings={meta.get('n_openings')} "
          f"faces={meta.get('n_faces')}")
    stl.unlink()
    if side.exists():
        side.unlink()

bad = identical and green
print()
print("BUG CONFIRMED" if bad else "NO BUG (empty region is refused)")
sys.stdout.flush()
