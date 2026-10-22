"""The trim feature end-to-end through the GUI, on the JK5 abdominal aorta."""
import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import boundary_loops, pre_bifurcation_candidates

app = QtWidgets.QApplication(sys.argv)
fails = []

w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
for _ in range(10): app.processEvents()

print("=" * 78)
print("BASELINE, abdominal aorta, no trims")
print("=" * 78)
print(f"  rows: {len(w.rows)}   status: {w.status.text()[:60]}")
print(f"  openings: {len(boundary_loops(w.cut))}  euler {int(w.cut.euler_number)}")
print(f"  trim menu offers: {len(w.trim_btn.menu().actions())} candidates")
for a in w.trim_btn.menu().actions():
    print(f"    - {a.text()}")
base_rows = len(w.rows)
if w.problem is not None:
    fails.append(("baseline-invalid", w.problem[:60]))
if not w.trim_btn.isEnabled():
    fails.append("trim-btn-disabled")

# ---- trim the RIGHT common iliac (Line15 -> removes its 2 distal branches) ----
print()
print("=" * 78)
print("TRIM 1, right common iliac (Line15), before its internal/external split")
print("=" * 78)
cands = pre_bifurcation_candidates(w.tree)
c15 = next(c for c in cands if c["segment"] == "Line15")
w._add_trim(c15)
for _ in range(6): app.processEvents()
row15 = next(r for r in w.rows if getattr(r.terminal, "is_interior", False))
print(f"  rows: {base_rows} -> {len(w.rows)}   (2 branch rows replaced by 1 trim row)")
print(f"  trim row: '{row15.terminal.name}'  offset={row15.value():.1f}mm  "
      f"limit={row15.limit:.1f}mm  readout={row15.readout.text()!r}")
print(f"  openings: {len(boundary_loops(w.cut))}  euler {int(w.cut.euler_number)}")
print(f"  status: {w.status.text()[:70]}")
if len(w.rows) != base_rows - 1:
    fails.append(("trim1-rowcount", len(w.rows)))
if w.problem is not None:
    fails.append(("trim1-invalid", w.problem[:60]))
if len(boundary_loops(w.cut)) != len(w.rows):
    fails.append(("trim1-openings", len(boundary_loops(w.cut)), len(w.rows)))

# ---- trim the LEFT common iliac too ----
print()
print("=" * 78)
print("TRIM 2, left common iliac (Line16) as well")
print("=" * 78)
c16 = next(c for c in pre_bifurcation_candidates(w.tree) if c["segment"] == "Line16")
w._add_trim(c16)
for _ in range(6): app.processEvents()
print(f"  rows: {len(w.rows)}   trims: {len(w.trims)}")
print(f"  openings: {len(boundary_loops(w.cut))}  euler {int(w.cut.euler_number)}")
print(f"  status: {w.status.text()[:70]}")
for r in w.rows:
    tag = "TRIM" if getattr(r.terminal, "is_interior", False) else "    "
    print(f"    {tag} {r.terminal.name:34s} {r.value():5.1f}mm  {r.readout.text()}")
if w.problem is not None:
    fails.append(("trim2-invalid", w.problem[:60]))
if len(boundary_loops(w.cut)) != len(w.rows):
    fails.append(("trim2-openings",))

# ---- move the trim's plane and confirm it stays valid ----
print()
print("=" * 78)
print("SLIDE the trim cut back along the common iliac")
print("=" * 78)
r15 = next(r for r in w.rows if r.terminal.label == "Line15.trim")
for off in [6.0, 12.0, 20.0, 28.0]:
    r15.set_value(off)
    w._on_change()
    for _ in range(3): app.processEvents()
    print(f"  offset {off:5.1f}mm -> openings {len(boundary_loops(w.cut))}, "
          f"{r15.readout.text():>12}, {'OK' if w.problem is None else 'RED: ' + w.problem[:40]}")
    if w.problem is not None:
        fails.append(("slide", off, w.problem[:40]))
r15.set_value(10.0); w._on_change()
for _ in range(3): app.processEvents()

# ---- undo one trim ----
print()
print("=" * 78)
print("UNDO the right-iliac trim")
print("=" * 78)
w._remove_trim(r15.terminal)
for _ in range(6): app.processEvents()
print(f"  rows: {len(w.rows)}   trims: {len(w.trims)}   "
      f"openings {len(boundary_loops(w.cut))}")
print(f"  status: {w.status.text()[:70]}")
if len(w.trims) != 1:
    fails.append(("undo-trims", len(w.trims)))
if w.problem is not None:
    fails.append(("undo-invalid", w.problem[:60]))

# ---- export metadata records the trim ----
print()
print("=" * 78)
print("EXPORT metadata")
print("=" * 78)
trim_rows = [r for r in w.rows if getattr(r.terminal, "is_interior", False)]
print(f"  trim rows remaining: {[r.terminal.name for r in trim_rows]}")
print(f"  subtree recorded: {sorted(trim_rows[0].terminal.subtree_segments)}")

w.close()
print()
print("=" * 78)
print("ALL PASS" if not fails else f"FAILURES: {fails}")
print("=" * 78)

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
