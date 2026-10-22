"""Naming an opening by hand, and the regions that follow from it.

The point of the feature is not the label. Over 47 abdominal scans name_anatomy
names no renal on 18, and on 14 of those because several lateral pairs fit the renal
pattern equally well -- so it correctly refuses to choose, and region_presets then
offers "Whole abdominal aorta" and nothing else. A person can see which pair is the
renals. This checks that saying so actually rebuilds the presets, rather than just
changing a caption.

    python tests/test_manual_naming.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from branch_tree import load_tree                                    # noqa: E402
from branch_cut import (load_mat_mesh, classify_terminals, name_anatomy,   # noqa: E402
                        region_presets, anatomy_choices, set_manual_name,
                        ANATOMY_NAMES_ABDOMINAL, ANATOMY_NAMES_THORACIC)

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


def load(csv, mat):
    tree = load_tree(ROOT / csv)
    mesh = load_mat_mesh(ROOT / mat)
    classify_terminals(tree, mesh, verbose=False)
    return tree, mesh


def presets_of(tree, mesh):
    names, vessels, _ = name_anatomy(tree, mesh)
    return [p["name"] for p in region_presets(tree, mesh, names, vessels)], names


print("=" * 78)
print("the list offered depends on the anatomy")
print("=" * 78)
ab_tree, ab_mesh = load("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat")
th_tree, th_mesh = load("KY160_1_CentrelinePoints.csv",
                        "Geometries & Meshes/mesh/KY160_1.mat")
ab_ch, th_ch = anatomy_choices(ab_tree), anatomy_choices(th_tree)
print(f"  abdominal: {len(ab_ch)} names, first {ab_ch[0]!r}")
print(f"  thoracic : {len(th_ch)} names, first {th_ch[0]!r}")
check(ab_ch == ANATOMY_NAMES_ABDOMINAL, "an abdominal aorta is offered abdominal names")
check(th_ch == ANATOMY_NAMES_THORACIC, "a thoracic aorta is offered thoracic names")
check(not (set(ab_ch) & {"Left subclavian artery", "Ascending aorta"}),
      "and not arch vessels, which cannot be on an abdominal aorta")
# region_presets matches on words, so the words have to be there to match.
check(sum(1 for n in ab_ch if "renal" in n.lower()) >= 4,
      "accessory and polar renals carry the word 'renal', so they count toward "
      "the lowest-renal boundary")

print()
print("=" * 78)
print("declaring a renal on a scan where the geometry cannot name one")
print("=" * 78)
# Strip the automatic renals to reproduce the 18-scan situation on a model that
# ships with the repository: JK5 names its renals fine, so take them away first.
base_presets, base_names = presets_of(ab_tree, ab_mesh)
auto_renals = sorted(k for k, v in base_names.items() if "renal" in v.lower())
print(f"  as loaded: {len(base_presets)} presets, renals named {auto_renals}")
check("Aortic neck" in base_presets, "JK5 names its renals automatically")

# Rename both renals to something with no anatomical meaning, which is exactly what
# a Branch N looks like to region_presets.
renamed = {}
for t in ab_tree.terminals:
    if t.label in auto_renals:
        renamed[t.label] = t.name
        set_manual_name(ab_tree, t.label, "Lumbar artery")
blind_presets, blind_names = presets_of(ab_tree, ab_mesh)
print(f"  renals called lumbar: {len(blind_presets)} presets -> {blind_presets}")
check("Aortic neck" not in blind_presets and "Infrarenal aorta" not in blind_presets,
      "with no renal named, the renal-dependent presets are gone")
check("Whole abdominal aorta" in blind_presets,
      "and the one that needs no landmark survives")

# Now declare them by hand, which is the whole feature.
labels = sorted(renamed)
set_manual_name(ab_tree, labels[0], "Left renal artery")
set_manual_name(ab_tree, labels[1], "Right renal artery")
hand_presets, hand_names = presets_of(ab_tree, ab_mesh)
print(f"  declared by hand    : {len(hand_presets)} presets -> {hand_presets}")
check("Aortic neck" in hand_presets, "declaring the renals brings back the aortic neck")
check("Infrarenal aorta" in hand_presets, "and the infrarenal aorta")
check("Infrarenal aorta with iliacs" in hand_presets, "and the version with iliacs")
check(len(hand_presets) >= len(base_presets),
      f"a hand-named case is no worse off than an auto-named one "
      f"({len(hand_presets)} vs {len(base_presets)})")

print()
print("=" * 78)
print("a hand-named artery displaces the automatic one of the same name")
print("=" * 78)
# Otherwise region_presets sees two "renal" landmarks and places the aortic neck at
# whichever is lower -- possibly the very vessel being corrected.
set_manual_name(ab_tree, labels[0], None)          # let the left be automatic again
other = next(t.label for t in ab_tree.terminals
             if t.label not in labels and not getattr(t, "is_trunk", False))
set_manual_name(ab_tree, other, "Right renal artery")
dup_names = name_anatomy(ab_tree, ab_mesh)[0]
right = sorted(k for k, v in dup_names.items() if v == "Right renal artery")
print(f"  openings called 'Right renal artery': {right}")
check(right == [other], "exactly one opening carries the declared name")
check(dup_names.get(other) == "Right renal artery", "and it is the declared one")

print()
print("=" * 78)
print("clearing a declaration returns the opening to the geometry's answer")
print("=" * 78)
for lbl in list(getattr(ab_tree, "manual_names", {})):
    set_manual_name(ab_tree, lbl, None)
back_presets, back_names = presets_of(ab_tree, ab_mesh)
print(f"  {len(back_presets)} presets, renals {sorted(k for k, v in back_names.items() if 'renal' in v.lower())}")
check(back_presets == base_presets,
      "everything is exactly as it was before any name was declared")
check(not getattr(ab_tree, "manual_names", {}),
      "and no declaration is left behind")

print()
print("=" * 78)
print("a declaration for an opening that is not in this model is ignored")
print("=" * 78)
set_manual_name(ab_tree, "NotALabel.end", "Left renal artery")
stale_names = name_anatomy(ab_tree, ab_mesh)[0]
check("NotALabel.end" not in stale_names,
      "a label from another model does not appear in the names")
check(stale_names == base_names,
      "and it changes nothing else")
set_manual_name(ab_tree, "NotALabel.end", None)

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
