"""Say so on screen when left and right were never confirmed, and draw the frame.

Reported from the abdominal cohort: "Sometimes it mixes left and right iliacs."

Side is decided by ONE scalar for a whole case -- left_sign in name_anatomy -- read
from a witness in the anatomy: the splenic/hepatic pair on an abdomen, the
brachiocephalic on an arch. When no witness is found, left_sign falls back to +1 and
a note is appended saying the sides now rest on ANATOMY_FRAME alone. 11 of 59 cases
in the last abdominal cohort run carry that note.

The note reached only the exported JSON. On screen the caption read exactly the same
whether handedness had been READ from the anatomy or GUESSED from the frame, so the
one case a user needed to check looked identical to the 48 that were sound.

Two things are asserted here, and neither claims to know which way round any
particular case went -- that is the point. They make the difference visible:

  1. the caption distinguishes the two, and reverts when the note goes away
  2. the view carries an L/P/S orientation triad, so a side can be checked against
     something. An anterior view puts patient-LEFT on the viewer's RIGHT, which is
     its own way to misread an iliac and has precedent in this repo:
     tests/test_trim_gui.py calls Line15 "the RIGHT common iliac" on JK5_1, where
     the tool says Left and the geometry agrees with the tool.

    .venv\\Scripts\\python tests\\test_sides_unconfirmed.py
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from qtpy import QtWidgets
from branch_gui import CutWindow
from branch_cut import SIDES_ASSUMED_NOTE

app = QtWidgets.QApplication(sys.argv[:1])
QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: None)
QtWidgets.QMessageBox.critical = staticmethod(lambda *a, **k: None)
fails = []

w = CutWindow(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
              str(ROOT / "JK5_1_M10.mat"))
w.show()
for _ in range(10):
    app.processEvents()

print("=" * 76)
print("1. JK5_1 finds its celiac, so handedness IS witnessed")
print("=" * 76)
print(f"  notes   : {w.anat_notes}")
confirmed_caption = w.anat_lbl.text()
print(f"  caption : {confirmed_caption[:120]}")
if any(SIDES_ASSUMED_NOTE in str(n) for n in w.anat_notes):
    fails.append("JK5_1 unexpectedly reports unconfirmed sides")
if "NOT CONFIRMED" in confirmed_caption:
    fails.append("a witnessed case is being warned about as if it were not")
if "named anatomically" not in confirmed_caption:
    fails.append("the ordinary caption is missing")

print()
print("=" * 76)
print("2. the same case WITH the note: the caption must change")
print("=" * 76)
# Injected rather than found: no local geometry reaches the fallback, and the point
# under test is the display, not the naming. The exact sentence comes from
# branch_cut so this cannot pass against a marker that drifted.
w.anat_notes = ["no splenic/hepatic pair found to confirm left from right, "
                + SIDES_ASSUMED_NOTE]
w._update_anat_label()
warned = w.anat_lbl.text()
print(f"  caption : {warned[:200]}")
if "NOT CONFIRMED" not in warned:
    fails.append("the unconfirmed-sides warning never reached the caption")
if warned == confirmed_caption:
    fails.append("the caption is identical with and without the note")
if "mirrored export" not in warned.lower():
    fails.append("the caption does not say what would go wrong")

print()
print("=" * 76)
print("3. and it reverts when the note goes")
print("=" * 76)
w.anat_notes = []
w._update_anat_label()
back = w.anat_lbl.text()
print(f"  caption : {back[:120]}")
if "NOT CONFIRMED" in back:
    fails.append("the warning stuck after the note was cleared")
if back != confirmed_caption:
    fails.append("the caption did not return to what it was")

print()
print("=" * 76)
print("4. the view carries an anatomical orientation triad")
print("=" * 76)
ren = w.plotter.renderer
widget = (getattr(w.plotter, "axes_widget", None)
          or getattr(ren, "axes_widget", None))
actor = getattr(ren, "axes_actor", None) or getattr(ren, "marker_actor", None)
print(f"  axes widget : {type(widget).__name__ if widget is not None else None}")
print(f"  axes actor  : {type(actor).__name__ if actor is not None else None}")
if widget is None and actor is None:
    fails.append("no orientation marker in the desktop view")
else:
    labels = []
    if actor is not None:
        for get in ("GetXAxisLabelText", "GetYAxisLabelText", "GetZAxisLabelText"):
            fn = getattr(actor, get, None)
            if fn is not None:
                labels.append(fn())
    print(f"  labels      : {labels}")
    # L/P/S, not x/y/z: the frame is the thing being communicated.
    if labels and [s.upper() for s in labels] != ["L", "P", "S"]:
        fails.append(f"the triad is not labelled L/P/S: {labels}")

w.close()
print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
