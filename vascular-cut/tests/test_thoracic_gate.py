"""Nub-in and the Ishimaru zones exist on the web page for thoracic aortas only.

Nub-in is the cohort's rule measured over 82 vessels in 30 THORACIC cases, and the
Ishimaru zones are an arch classification; neither means anything on an abdominal
aorta. So an abdominal case must not be offered either one, must be told why, and must
be refused on the server when a page fires the trigger anyway. A wrong anatomy call is
corrected by setting Anatomy and applying it to the case ("Apply to this case").

This drives the real app with no browser: the page's conditions are read from the
built template, and the handlers are called against real sessions (local mode).

WHAT IT PINS
  1. JK5 (abdominal): is_thoracic False; run_zones refuses with a reason and leaves no
     zone state behind; the nub_preset trigger refuses and leaves the cut unchanged.
  2. KY160 (thoracic): is_thoracic True and nub-in offered; isolating a region turns
     it off again (its own reason, shown on the page).
  3. A thoracic case's zone colouring does not carry over to the next case.
  4. In the template, the Nub-in button, Find zones and the Zone colour button are
     only rendered under `is_thoracic`, and the reasons under `!is_thoracic`.
  5. "Apply to this case" on this machine rebuilds from the paths the case came from,
     never from the path boxes, and a thoracic call brings nub-in back.

    .venv\\Scripts\\python tests\\test_thoracic_gate.py
"""
import os
import sys
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
sys.path.insert(0, str(ROOT))
os.environ.pop("VC_PUBLIC", None)

import cut_app as APP                                              # noqa: E402

state = APP.state
KY_CSV = ROOT / "KY160_1_CentrelinePoints.csv"
KY_MAT = ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"
JK_CSV = ROOT / "JK5_1_CentrelinePoints_branches.csv"
JK_MAT = ROOT / "JK5_1_M10.mat"
fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        fails.append(name)


print("=" * 78)
print("1. an abdominal case")
print("=" * 78)
APP._open_session(str(JK_CSV), str(JK_MAT))
s = APP.SESSION["s"]
check("control: JK5 loaded", s is not None and state.loaded, state.status_text)
check("is_thoracic is False", state.is_thoracic is False)
check("the page is told it is abdominal", state.anat_kind == "abdominal", state.anat_kind)
check("can_nub is False", state.can_nub is False)

APP.run_zones()
check("run_zones refuses, saying why",
      "thoracic aorta only" in state.zones_note, state.zones_note)
check("  ...and leaves nothing running or coloured",
      not state.zones_busy and not state.colour_zones and not state.zone_legend
      and APP.ZONES["rgb"] is None)

before = [(o["label"], o["offset_mm"]) for o in s.openings()]
APP.nub_preset()
check("the nub_preset trigger refuses", state.status_bad
      and "thoracic" in state.status_text, state.status_text)
check("  ...and the cut is unchanged",
      [(o["label"], o["offset_mm"]) for o in s.openings()] == before)

print("\n" + "=" * 78)
print("2. a thoracic case")
print("=" * 78)
APP._open_session(str(KY_CSV), str(KY_MAT))
s = APP.SESSION["s"]
check("control: KY160 loaded", s is not None and state.loaded, state.status_text)
check("is_thoracic is True", state.is_thoracic is True)
check("anat_kind is thoracic", state.anat_kind == "thoracic", state.anat_kind)
check("nub-in is offered", state.can_nub is True)
check("the case is named", state.case_name == "KY160_1", state.case_name)
check("the load counter moved", state.load_seq >= 2, state.load_seq)
pick = next(p for p in state.region_options if p != "Whole vessel")
state.region_pick = pick
APP.on_region(pick)
check(f"isolating a region ({pick}) turns nub-in off",
      state.in_region and state.can_nub is False)
APP.clear_region()
# clear_region only sets the picker; in a browser its change listener does the rest,
# and listeners do not fire without one.
APP.on_region(state.region_pick)
check("  ...and the whole vessel turns it back on", state.can_nub is True)

print("\n" + "=" * 78)
print("3. nothing of one case carries over to the next")
print("=" * 78)
state.colour_zones = True
state.zone_legend = [{"zone": 0, "name": "x", "colour": "#000", "span": ""}]
state.zones_note = "zone note from KY160"
# ...and nor does a measurement: its readout, its plane, or the pinned marks. They
# used to be drawn on the next case, with the previous case's numbers beside them.
APP.run_measure()
check("control: KY160 measured, its plane drawn",
      state.measure is not None and "measure" in APP.SCENE._actors)
state.marks = [{"label": "Mark 1", "vessel": "Left subclavian", "mm_to_node": 12.0,
                "point": [0.0, 0.0, 0.0]}]
APP.SCENE.set_marks([[0.0, 0.0, 0.0]])
# The delta K view shares the render window: it must not show the cut view's
# measure plane or marks, and going back to Cut must put them back.
APP.on_mode("shape")
check("the delta K view shows no measure plane and no marks",
      "measure" not in APP.SCENE._actors
      and not any(k.startswith("mark:") for k in APP.SCENE._actors))
APP.on_mode("cut")
check("  ...and back in Cut both are drawn again", "measure" in APP.SCENE._actors
      and any(k.startswith("mark:") for k in APP.SCENE._actors))
APP._open_session(str(JK_CSV), str(JK_MAT))
check("colour_zones is off after loading another case", state.colour_zones is False)
check("the legend and note are gone", not state.zone_legend and not state.zones_note)
check("the measurement is gone", state.measure is None and state.angle is None
      and not state.m_note)
check("so are its plane and the marks", "measure" not in APP.SCENE._actors
      and not state.marks and not any(k.startswith("mark:") for k in APP.SCENE._actors))

print("\n" + "=" * 78)
print("4. the page renders the thoracic tools under is_thoracic only")
print("=" * 78)


class Conds(HTMLParser):
    """For each piece of text in the template, the v-if conditions above it."""

    def __init__(self):
        super().__init__()
        self.stack, self.hits = [], []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self.stack.append((tag, a.get("v-if") or a.get("v-show") or ""))

    def handle_startendtag(self, tag, attrs):
        pass

    def handle_endtag(self, tag):
        for k in range(len(self.stack) - 1, -1, -1):
            if self.stack[k][0] == tag:
                del self.stack[k:]
                break

    def handle_data(self, data):
        t = data.strip()
        if t:
            self.hits.append((t, [c for _, c in self.stack if c]))


p = Conds()
p.feed(APP.layout.html)


def conds_of(text):
    return [c for t, c in p.hits if t == text]


def under(text, cond):
    """Every occurrence of text sits under a v-if that requires cond."""
    found = conds_of(text)
    return bool(found) and all(any(cond in c.replace("!" + cond, "")
                                   for c in cs) for cs in found)


def under_not(text, cond):
    found = conds_of(text)
    return bool(found) and all(any("!" + cond in c for c in cs) for cs in found)


check("the Nub-in card", under("Nub-in cut", "is_thoracic"), conds_of("Nub-in cut"))
check("Find zones", under("Find zones", "is_thoracic"), conds_of("Find zones"))
check("the Zone colour button", under("Zone", "is_thoracic"), conds_of("Zone"))
check("the Nub-in reason, under !is_thoracic", under_not(APP.TXT_NUB_OFF, "is_thoracic"))
check("the zones reason, under !is_thoracic", under_not(APP.TXT_ZONES_OFF, "is_thoracic"))
check("the Nub-in text states the aorta's radius, not the vessel's",
      "aorta’s local radius" in APP.TXT_NUB_RULE
      and "own radius" not in APP.layout.html)

print("\n" + "=" * 78)
print("5. Apply to this case, on this machine")
print("=" * 78)
check("control: JK5 is loaded, read as abdominal", state.anat_kind == "abdominal"
      and state.loaded_family == "auto")
state.mat_path = str(KY_MAT)       # the boxes now name another case entirely
state.csv_path = str(KY_CSV)
state.flush()                      # (a path change resets the anatomy to auto)
state.family = "thoracic"
state.flush()
APP.reload_case()
check("the case was rebuilt as thoracic", state.is_thoracic is True
      and state.anat_kind == "thoracic" and state.anat_how == "set by hand",
      f"{state.anat_kind} / {state.anat_how}")
check("from the files it came from, not the path boxes",
      state.case_name == "JK5_1_M10" and APP.LOADED_FROM["mat"] == str(JK_MAT),
      state.case_name)
check("what it was loaded with is recorded", state.loaded_family == "thoracic")
check("nub-in is now offered", state.can_nub is True)
APP._wipe_visit()
check("the wipe forgets where the case came from",
      APP.LOADED_FROM == {"csv": None, "mat": None})
APP.reload_case()
check("with nothing loaded, Apply refuses and says so",
      state.status_bad and "nothing to read again" in state.status_text.lower(),
      state.status_text)

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
