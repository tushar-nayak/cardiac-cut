"""What the web app can now do: declare an artery, and apply the nub preset.

The trame app drives CutSession and nothing else, so checking it here checks the
website without a browser -- which is the only practical way to test a cut that takes
seconds and a preset that re-cuts the whole model.

    python tests/test_api_naming_nub.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import nub_cut as NC                                                # noqa: E402
from vascular_api import CutSession, CutError                       # noqa: E402

TH_CSV = ROOT / "KY160_1_CentrelinePoints.csv"
TH_MAT = ROOT / "Geometries & Meshes/mesh/KY160_1.mat"
AB_CSV = ROOT / "JK5_1_CentrelinePoints_branches.csv"
AB_MAT = ROOT / "JK5_1_M10.mat"

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


print("=" * 78)
print("declaring an artery drives the region presets")
print("=" * 78)
ab = CutSession(str(AB_CSV), str(AB_MAT))
choices = ab.anatomy_choices()
print(f"  {len(choices)} names offered, e.g. {choices[:3]}")
check("Left renal artery" in choices, "the renals are offered on an abdominal aorta")
check("Left subclavian artery" not in choices, "and arch vessels are not")

base = [p["name"] for p in ab.presets()]
renals = [k for k, v in ab.anat_names.items() if "renal" in v.lower()]
print(f"  as loaded: {len(base)} presets, renals at {renals}")

# Reproduce the 18-scan situation: no renal named, so no renal landmark.
for lbl in renals:
    ab.declare_artery(lbl, "Lumbar artery")
blind = [p["name"] for p in ab.presets()]
print(f"  renals hidden: {blind}")
check("Aortic neck" not in blind, "with no renal, the aortic neck is not offered")

ab.declare_artery(renals[0], "Left renal artery")
ab.declare_artery(renals[1], "Right renal artery")
hand = [p["name"] for p in ab.presets()]
print(f"  declared     : {hand}")
check("Aortic neck" in hand, "declaring them brings the aortic neck back")
check("Infrarenal aorta" in hand, "and the infrarenal aorta")
check(sorted(hand) == sorted(base), "the full preset list is restored")
check(ab.declared_arteries() == {renals[0]: "Left renal artery",
                                 renals[1]: "Right renal artery"},
      "and the declarations are reported as such")

# The presets have to be usable, not just listed -- a name that produces a region
# nobody can cut would be worse than no name.
p = next(x for x in ab.presets() if x["name"] == "Aortic neck")
ab.apply_preset("Aortic neck")
st = ab.recut()
print(f"  aortic neck  : {p['s0']:.0f}-{p['s1']:.0f} mm -> {st['message']}")
check(st["ok"], f"the region declared by hand actually cuts: {st['message']}")
ab.clear_region()
ab.recut()

for lbl in list(ab.declared_arteries()):
    ab.declare_artery(lbl, None)
check([q["name"] for q in ab.presets()] == base,
      "clearing every declaration returns the presets to the geometry's answer")
try:
    ab.declare_artery(renals[0], "Left subclavian artery")
    check(False, "a name this anatomy cannot carry is refused")
except CutError as e:
    check("not an artery" in str(e), f"a name this anatomy cannot carry is refused")

print()
print("=" * 78)
print("the nub-in preset, applied through the same API the web app uses")
print("=" * 78)
th = CutSession(str(TH_CSV), str(TH_MAT))
before = {o["label"]: o["offset_mm"] for o in th.openings()}
# Nub-in is for thoracic aortas only: refused on the abdominal case, with a reason,
# and the case left exactly as it was.
from vascular_api import CutError                                  # noqa: E402
ab_before = (len(ab.mesh.faces), dict(ab._offsets))
try:
    ab.nub_preset()
    refused = None
except CutError as e:
    refused = str(e)
check(refused is not None and "thoracic" in refused,
      "nub-in is refused on an abdominal aorta, saying why")
check((len(ab.mesh.faces), dict(ab._offsets)) == ab_before,
      "  ...and the abdominal case is left untouched")
check(th.is_thoracic and not ab.is_thoracic, "is_thoracic agrees with the anatomy call")

notes = th.nub_preset()
for n in notes[:8]:
    print(f"    {n}")
after = {o["label"]: o["offset_mm"] for o in th.openings()}
st = th.status()
print(f"  {st['message']}")
check(notes and notes[0].startswith("nub-in cut applied"),
      "it says how many openings it placed")
check(st["ok"], f"the result is a valid cut: {st['message']}")
check(set(after) != set(before) or any(abs(after[k] - before[k]) > 0.05
                                       for k in after if k in before),
      "the openings moved to the preset")

# The two things the preset does that no single control does: a stub on every branch,
# and the aorta's own ends squared off horizontally.
trunks = [o for o in th.openings() if o.get("is_trunk")]
print(f"  trunk ends: {[(o['name'], o['offset_mm'], o['cut_mode']) for o in trunks]}")
check(all(abs(o["offset_mm"] - NC.TRUNK_OFFSET_MM) < 0.05
          or o["cut_mode"] != NC.TRUNK_MODE for o in trunks),
      f"the aorta's ends sit at {NC.TRUNK_OFFSET_MM} mm, or fell back and said so")

# Refused inside a region, because the preset squares off ends a region has replaced.
th.set_region(10.0, 60.0, "test")
th.recut()
try:
    th.nub_preset()
    check(False, "refused while a region is isolated")
except CutError as e:
    check("region" in str(e).lower(), f"refused while a region is isolated: {e}")
th.clear_region()
th.recut()

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
