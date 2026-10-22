"""Feature B (live diameter) + Feature A (click-to-mark + span) through the widget."""
import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import section_diameter, centerline_geodesic

app = QtWidgets.QApplication(sys.argv)
fails = []


class FakePicker:
    def __init__(self): self.world = None; self.cell = 3
    def Pick(self, x, y, z, ren): return 1 if self.cell >= 0 else 0
    def GetCellId(self): return self.cell
    def GetPickPosition(self): return tuple(self.world)


w = CutWindow(str(ROOT / "KY160_1_CentrelinePoints.csv"),
              str(ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"))
w.show()
for _ in range(6): app.processEvents()

print("=" * 76)
print("FEATURE B: live diameter at the plane while a slider is dragged")
print("=" * 76)
# the proximal trunk end, now named anatomically ("Ascending aorta"), so pick it
# by role rather than by a hard-coded label
r = max((x for x in w.rows if x.terminal.is_trunk),
        key=lambda x: float(x.terminal.tip[2]))
for off in [6.0, 20.0, 40.0]:
    r.box.blockSignals(True); r.box.setValue(off); r.box.blockSignals(False)
    w._on_preview(r)                        # simulate a drag frame
    txt = r.readout.text()
    sd = section_diameter(w.base, *r.terminal.cut_plane(off), near_point=r.terminal.cut_plane(off)[0])
    ok = "⌀" in txt and f"{sd['ead_mm']:.1f}" in txt
    print(f"  offset {off:5.1f} -> readout {txt!r}   (EAD {sd['ead_mm']:.1f})  {'OK' if ok else 'FAIL'}")
    if not ok: fails.append(("live-dia", off, txt))
# tooltip carries min/max + circularity
print(f"  tooltip: {r.readout.toolTip().splitlines()[1:]!r}")

print()
print("=" * 76)
print("FEATURE A: click to pin marks")
print("=" * 76)
fake = FakePicker()
w._hover_picker = fake
w.measure_cb.setChecked(True)
for _ in range(2): app.processEvents()

# place mark 1
v1 = w.base.vertices[500]; fake.world = v1
w._place_mark((400, 300))
print(f"  after 1 click: marks={len(w._marks)}  label={w.marks_lbl.text()!r}  "
      f"clear-enabled={w.clear_marks_btn.isEnabled()}")
mark_actors = [n for n in w.plotter.renderer.actors if n.startswith("mark_")]
print(f"    mark actors: {mark_actors}")
if len(w._marks) != 1 or not w.clear_marks_btn.isEnabled():
    fails.append(("mark1",))

# place mark 2 -> span appears
v2 = w.base.vertices[9000]; fake.world = v2
w._place_mark((420, 320))
g = centerline_geodesic(w.tree, v1, v2)
print(f"  after 2 clicks: marks={len(w._marks)}  label={w.marks_lbl.text()!r}")
print(f"    expected span (geodesic) = {g['distance_mm']:.1f} mm")
# span is reported as a NUMBER only; the connecting line was intentionally removed
span_ok = len(w._marks) == 2 and f"{g['distance_mm']:.1f}" in w.marks_lbl.text() \
          and "mark_span" not in w.plotter.renderer.actors
print(f"    span number shown, no line drawn: {'OK' if span_ok else 'FAIL'}")
if not span_ok: fails.append(("span", w.marks_lbl.text()))

print()
print("=" * 76)
print("MARKS SURVIVE A RE-CUT (draw() re-renders them)")
print("=" * 76)
r.box.setValue(30.0)            # triggers commit -> recut -> draw()
for _ in range(4): app.processEvents()
still = [n for n in w.plotter.renderer.actors if n.startswith("mark_pt_")]
print(f"  marks in data after re-cut : {len(w._marks)}")
print(f"  mark actors after re-cut   : {len(still)}  (must be re-rendered)")
if len(w._marks) != 2 or len(still) != 2:
    fails.append(("survive-recut", len(w._marks), len(still)))

print()
print("=" * 76)
print("CLICK DETECTION: press+release without drag = mark; with drag = no mark")
print("=" * 76)
w._clear_marks()
w.plotter.iren.get_event_position = lambda: (500, 400)
fake.world = w.base.vertices[3000]
w._lclick_xy = (500, 400)
w._left_release()                # same pos as press -> click
print(f"  no-drag release -> marks={len(w._marks)}  (expect 1)")
if len(w._marks) != 1: fails.append(("click-detect",))
w._clear_marks()
w._lclick_xy = (500, 400)
w.plotter.iren.get_event_position = lambda: (560, 470)   # moved a lot -> a rotate
w._left_release()
print(f"  drag release    -> marks={len(w._marks)}  (expect 0)")
if len(w._marks) != 0: fails.append(("drag-no-mark",))

print()
print("=" * 76)
print("CLEAR MARKS")
print("=" * 76)
fake.world = w.base.vertices[100]; w._place_mark((1, 1))
fake.world = w.base.vertices[200]; w._place_mark((2, 2))
print(f"  before clear: {len(w._marks)} marks")
w._clear_marks()
left = [n for n in w.plotter.renderer.actors if n.startswith("mark_")]
print(f"  after clear : {len(w._marks)} marks, {len(left)} actors, "
      f"label={w.marks_lbl.text()!r}, clear-enabled={w.clear_marks_btn.isEnabled()}")
if len(w._marks) != 0 or len(left) != 0 or w.clear_marks_btn.isEnabled():
    fails.append(("clear",))

w.close()
print()
print("=" * 76)
print("ALL PASS" if not fails else f"FAILURES: {fails}")
print("=" * 76)

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
