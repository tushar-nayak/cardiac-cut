"""Exercise the four quick-win features: in-plane shift, move-all-branches,
send-to-mark, and per-branch colour."""
import sys
import math
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import (open_terminals, validate_cut, snap_tolerance,
                        measure_along_branch, vessel_coords)

app = QtWidgets.QApplication(sys.argv)
fails = []


class FakePicker:
    def __init__(self): self.world = None
    def Pick(self, x, y, z, ren): return 1
    def GetCellId(self): return 3
    def GetPickPosition(self): return tuple(self.world)


w = CutWindow(str(ROOT / "KY160_1_CentrelinePoints.csv"),
              str(ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"))
w.show()
for _ in range(8): app.processEvents()

branch_rows = [r for r in w.rows if not r.terminal.is_trunk]
trunk_rows = [r for r in w.rows if r.terminal.is_trunk]
print(f"rows: {len(w.rows)}  branches={len(branch_rows)}  trunks={len(trunk_rows)}")

# ---------------------------------------------------------------- 1. in-plane shift
print("\n" + "=" * 74)
print("1. IN-PLANE SHIFT, cut_plane stays TRUE; cut_anchor carries the shift")
print("=" * 74)
t = branch_rows[0].terminal
off = 6.0
p_true, n0 = t.cut_plane(off)
t.shift_a_mm, t.shift_b_mm = 5.0, -3.0
p_true2, n2 = t.cut_plane(off)          # opening location: must NOT move with shift
p_anch, n_anch = t.cut_anchor(off)      # band centre: carries the shift
d = p_anch - p_true
print(f"  cut_plane unchanged by shift : {np.allclose(p_true, p_true2) and np.allclose(n0, n2)}")
print(f"  anchor keeps the plane normal: {np.allclose(n_anch, n0)}")
print(f"  anchor shift ⟂ normal        : |d·n|={abs(d @ n0):.2e} mm")
print(f"  |shift| = √(5²+3²)           : got {np.linalg.norm(d):.4f}, expect {math.hypot(5,3):.4f}")
if not (np.allclose(p_true, p_true2) and np.allclose(n0, n2)):
    fails.append("cutplane-moved-by-shift")
if not np.allclose(n_anch, n0): fails.append("anchor-changed-normal")
if abs(d @ n0) > 1e-9: fails.append("anchor-not-in-plane")
if abs(np.linalg.norm(d) - math.hypot(5, 3)) > 1e-6: fails.append("anchor-magnitude")

# a modest shift must still cut cleanly
t.shift_a_mm, t.shift_b_mm = 2.0, 0.0
offs = {x.terminal.label: 6.0 for x in w.rows}
try:
    out = open_terminals(w.base, w.tree, offs, 6.0,
                         snap_mm=snap_tolerance(w.base), verbose=False)
    prob = validate_cut(out, w.tree, offs)
    print(f"  cut with a 2mm shift  : {'VALID' if prob is None else 'INVALID: ' + prob[:50]}")
    if prob is not None: fails.append(("shift-cut-invalid", prob[:40]))
except Exception as e:
    print(f"  cut with a 2mm shift  : EXCEPTION {e}")
    fails.append(("shift-cut-exception", str(e)[:40]))
t.shift_a_mm = t.shift_b_mm = 0.0        # reset
branch_rows[0].terminal.reset_plane()

# ---------------------------------------------------------------- 2. move all branches
print("\n" + "=" * 74)
print("2. MOVE ALL BRANCHES, nudges every branch, leaves the trunk alone")
print("=" * 74)
before_b = {r.terminal.label: r.value() for r in branch_rows}
before_t = {r.terminal.label: r.value() for r in trunk_rows}
w._nudge_branches(-5.0)                   # 5 mm toward the tips
for _ in range(3): app.processEvents()
ok_b = all(abs(r.value() - max(before_b[r.terminal.label] - 5.0, 1.0)) < 0.05
           for r in branch_rows)
ok_t = all(abs(r.value() - before_t[r.terminal.label]) < 1e-6 for r in trunk_rows)
for r in branch_rows:
    print(f"  {r.terminal.name:12s}: {before_b[r.terminal.label]:5.1f} -> {r.value():5.1f} mm")
for r in trunk_rows:
    print(f"  {r.terminal.name:22s}: {before_t[r.terminal.label]:5.1f} -> {r.value():5.1f} mm (unchanged)")
print(f"  branches moved -5mm: {ok_b}   trunk untouched: {ok_t}")
if not ok_b: fails.append("nudge-branches")
if not ok_t: fails.append("nudge-touched-trunk")
w._nudge_branches(5.0)                    # put them back
for _ in range(3): app.processEvents()
# no-op feedback when the branches are already pinned at the range floor
for r in branch_rows: r.set_value(1.0)
w._nudge_branches(-5.0)                    # cannot go below the floor
for _ in range(2): app.processEvents()
noop_ok = "nothing to move" in w.status.text().lower()
print(f"  nudge at the floor gives feedback (not silent): {noop_ok}  -> {w.status.text()[:46]!r}")
if not noop_ok: fails.append("nudge-noop-silent")
for r in branch_rows: r.set_value(6.0)
w._on_change()
for _ in range(2): app.processEvents()

# ---------------------------------------------------------------- 3. send cut to mark
print("\n" + "=" * 74)
print("3. SEND CUT TO MARK, plane jumps to a marked centerline point")
print("=" * 74)
fake = FakePicker(); w._hover_picker = fake
w.measure_cb.setChecked(True)
for _ in range(2): app.processEvents()
br = branch_rows[0]
seg = br.terminal.segment
s_target = 0.55 * seg.length
axis = seg.at_arclength(s_target)[0]
vtx = w.base.vertices[np.argmin(np.linalg.norm(w.base.vertices - axis, axis=1))]
fake.world = vtx
w._place_mark((400, 300))
res = measure_along_branch(w.tree, vtx)
print(f"  marked on {res['segment']} at s={res['s_mm']:.1f} mm (branch = {br.terminal.name})")
w._send_cut_to_mark()
for _ in range(3): app.processEvents()
t = br.terminal
expected = (t.segment.length - res["s_mm"]) if t.at_end else res["s_mm"]
expected = float(np.clip(expected, 1.0, br.limit))
print(f"  {t.name}'s offset -> {br.value():.1f} mm   (expected {expected:.1f})")
if abs(br.value() - expected) > 0.2:
    fails.append(("send-to-mark", br.value(), expected))
else:
    print("  OK, plane moved to the mark")

# ---------------------------------------------------------------- 4. per-branch colour
print("\n" + "=" * 74)
print("4. COLOUR OPENINGS, per-face tint, disc, and row label share a colour")
print("=" * 74)
w.colour_cb.setChecked(True)
for _ in range(3): app.processEvents()
rgb = w._opening_face_rgb(w.cut)
uniq = set(map(tuple, rgb)) if rgb is not None else set()
n_branch = sum(1 for r in w.rows if not r.terminal.is_trunk)
print(f"  face-rgb array: {None if rgb is None else rgb.shape}   distinct colours: {len(uniq)}"
      f"   (expect {n_branch} branches + red base = {n_branch + 1})")
# ONLY branch rows carry a hex colour; trunk rows must be red (untinted).
# Read the combo's LINE EDIT, not the combo. set_colour styles the line edit on
# purpose -- a stylesheet on the box itself repaints the whole widget and takes the
# drop-down arrow with it (branch_gui.py, OpeningRow.set_colour). This test read the
# combo, where the stylesheet is always "", so branch_tinted was always False and
# trunk_plain always passed VACUOUSLY: it was asserting nothing in either direction.
branch_tinted = all(
    w.opening_colour[r.terminal.label].lower()
    in r.name.lineEdit().styleSheet().lower()
    for r in w.rows if not r.terminal.is_trunk)
trunk_plain = all(
    "color:#" not in r.name.lineEdit().styleSheet()
    for r in w.rows if r.terminal.is_trunk)
print(f"  branch labels tinted: {branch_tinted}   trunk labels stay red/plain: {trunk_plain}")
label_ok = branch_tinted and trunk_plain
# the trunk must NOT be assigned a colour
trunk_uncoloured = all(r.terminal.label not in w.opening_colour
                       for r in w.rows if r.terminal.is_trunk)
print(f"  trunk excluded from the palette: {trunk_uncoloured}")
if not trunk_uncoloured: fails.append("trunk-coloured")
if len(uniq) != n_branch + 1: fails.append(("wrong-colour-count", len(uniq)))
# plane discs got the opening colours (actors exist)
disc_actors = [n for n in w.plotter.renderer.actors if n.startswith("plane_")]
print(f"  plane disc actors: {len(disc_actors)} (expect {len(w.rows)})")
if rgb is None or len(uniq) < 2: fails.append("colour-facergb")
if not label_ok: fails.append("colour-labels")
if len(disc_actors) != len(w.rows): fails.append("colour-discs")
# toggling off clears the label tint
w.colour_cb.setChecked(False)
for _ in range(2): app.processEvents()
off_clear = all("color:#" not in r.name.lineEdit().styleSheet() for r in w.rows)
print(f"  toggling off clears the tint: {off_clear}   (face-rgb now {w._opening_face_rgb(w.cut)})")
if not off_clear: fails.append("colour-toggle-off")

w.close()
print("\n" + "=" * 74)
print("ALL PASS" if not fails else f"FAILURES: {fails}")
print("=" * 74)

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
