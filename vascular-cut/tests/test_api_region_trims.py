"""Inside a region, only offer trims whose whole group is in view.

branch_gui._refresh_trim_menu has done this since the note in Whats-New: "The menu
is also scoped to what is in view: it used to offer trims for branches the region
does not contain." vascular_api.trim_candidates never got the same rule, so the web
listed every candidate in the tree while a region was active.

That is not cosmetic. _cut_list applies trims only through _trims_within, which keeps
a trim only when EVERY opening it swallows is inside the region. So a trim taken on a
group outside the region was recorded in the session and then silently dropped: the
caller got a success, a trim in the list, and no change to the mesh at all. Same fault
the cohort reported against the region sliders -- a control that reports success and
does nothing.

    .venv\\Scripts\\python tests\\test_api_region_trims.py
"""
import sys
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from vascular_api import CutSession                                  # noqa: E402

fails = []
s = CutSession(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
               str(ROOT / "JK5_1_M10.mat"))

whole = {c["segment"] for c in s.trim_candidates()}
print("=" * 76)
print("whole vessel: every trimmable group is offered")
print("=" * 76)
print(f"  {sorted(c['vessel'] for c in s.trim_candidates())}")
if not whole:
    print("FAILURES: nothing trimmable on JK5_1")
    raise SystemExit(1)

CASES = [
    # region, a vessel that must NOT be offered because its group is outside it
    ("Visceral segment", "iliac"),
    ("Infrarenal aorta with iliacs", "Celiac"),
]
names = {p["name"] for p in s.presets()}
for region, forbidden in CASES:
    if region not in names:
        fails.append(f"{region!r} not offered on JK5_1, so it is untested")
        continue
    s.clear_region()
    s.apply_preset(region)
    _b, inside, _n = s._region_terminals()
    visible = {t.label for t in inside}
    cands = s.trim_candidates()
    print()
    print("=" * 76)
    print(f"region {region!r}")
    print("=" * 76)
    print(f"  openings inside : {sorted(t.name for t in inside)}")
    print(f"  trims offered   : {sorted(c['vessel'] for c in cands)}")

    # Every offered trim's whole group must be inside the region -- that is the rule.
    from branch_cut import pre_bifurcation_candidates                # noqa: PLC0415
    all_c = {c["segment"]: c for c in pre_bifurcation_candidates(s.tree, s.base)}
    for c in cands:
        grp = {t.label for t in all_c[c["segment"]]["terminals"]}
        if not grp <= visible:
            fails.append(f"{region}: offered {c['vessel']!r} whose group is not "
                         f"wholly inside ({sorted(grp - visible)} outside)")
    if any(forbidden.lower() in c["vessel"].lower() for c in cands):
        fails.append(f"{region}: still offers a {forbidden!r} trim")
    # ...and it must still offer SOMETHING, or the scoping is just switched off.
    if not cands:
        fails.append(f"{region}: nothing offered at all, the scoping is too strict")

    # The point of the rule: anything it now hides would have been a silent no-op.
    hidden = whole - {c["segment"] for c in cands}
    for seg in sorted(hidden):
        grp = {t.label for t in all_c[seg]["terminals"]}
        if grp <= visible:
            fails.append(f"{region}: hid {seg}, whose group IS inside the region")

print()
print("=" * 76)
print("and leaving the region restores the full list")
print("=" * 76)
s.clear_region()
back = {c["segment"] for c in s.trim_candidates()}
print(f"  {len(back)} offered, was {len(whole)}")
if back != whole:
    fails.append(f"clearing the region did not restore the trim list: "
                 f"missing {sorted(whole - back)}, extra {sorted(back - whole)}")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
