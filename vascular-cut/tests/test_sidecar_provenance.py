"""The exported record must say which names a person declared, not just how many.

Found by the 28 Aug desktop/web parity audit. Both front ends write a sidecar with
the SAME FILENAME (<name>_cut_planes.json), so a downstream reader cannot tell which
one produced the file -- but they did not write the same thing. The window recorded

    anatomy.names          which opening got which name
    anatomy.named_by_hand  which of those a person DECLARED

and vascular_api.sidecar() recorded neither, only `named_openings`, a count.

Why it matters, in the window's own words: every side name inferred from the geometry
rests on the coordinate frame recorded beside it, and a mirrored export would swap all
of them -- but a hand-declared name does NOT rest on that frame and is not subject to
the doubt. This became sharper when the tool started warning that left and right were
never confirmed on a case (11 of 59 in the last abdominal cohort run): on such a case
the inferred sides are a frame assumption and the declared ones are not, and a reader
placing a measurement against a hand-named renal has to be able to tell which it is.

So this is a provenance defect in exported data, not a UI difference: any headless or
web export made so far records anatomy that cannot be audited.

    .venv\\Scripts\\python tests\\test_sidecar_provenance.py
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))

from vascular_api import CutSession                                  # noqa: E402

fails = []
s = CutSession(str(ROOT / "JK5_1_CentrelinePoints_branches.csv"),
               str(ROOT / "JK5_1_M10.mat"))

print("=" * 76)
print("1. the sidecar carries the name MAP, not just a count")
print("=" * 76)
side = s.sidecar()
anat = side.get("anatomy", {})
print(f"  keys           : {sorted(anat)}")
print(f"  named_openings : {anat.get('named_openings')}")
print(f"  names          : {len(anat.get('names', {}))} entries")
for k in ("names", "named_by_hand"):
    if k not in anat:
        fails.append(f"anatomy.{k} missing from the API sidecar")
if len(anat.get("names", {})) != anat.get("named_openings"):
    fails.append(f"names ({len(anat.get('names', {}))}) disagrees with "
                 f"named_openings ({anat.get('named_openings')})")
if anat.get("named_by_hand"):
    fails.append(f"nothing was declared yet, but named_by_hand is "
                 f"{anat['named_by_hand']}")

print()
print("=" * 76)
print("2. declare one by hand -- it must be distinguishable from the inferred")
print("=" * 76)
choices = s.anatomy_choices()
target = None
for o in s.openings():
    lbl = o.get("label")
    if lbl and not o.get("is_region_bound"):
        target = lbl
        break
pick = next((c for c in choices if "renal" in c.lower()), None) or choices[0]
print(f"  declaring {target} as {pick!r}")
s.declare_artery(target, pick)

side = s.sidecar()
anat = side["anatomy"]
by_hand = anat.get("named_by_hand", {})
names = anat.get("names", {})
print(f"  named_by_hand  : {by_hand}")
print(f"  is it in names : {target in names}  -> {names.get(target)!r}")
if target not in by_hand:
    fails.append(f"{target} was declared by hand but is not in named_by_hand")
elif by_hand[target] != pick:
    fails.append(f"named_by_hand records {by_hand[target]!r}, expected {pick!r}")
if names.get(target) != pick:
    fails.append(f"names records {names.get(target)!r} for the declared opening")
# The whole point: inferred and declared must be separable.
inferred = {k: v for k, v in names.items() if k not in by_hand}
print(f"  inferred       : {len(inferred)}   declared: {len(by_hand)}")
if not inferred:
    fails.append("every name reads as hand-declared, so the two cannot be told apart")

print()
print("=" * 76)
print("3. and it survives a real export to disk")
print("=" * 76)
tmp = Path(tempfile.mkdtemp(prefix="sidecar_"))
try:
    if s.ok:
        _stl, sfile = s.export(tmp, "provenance_check")
        on_disk = json.loads(Path(sfile).read_text(encoding="utf-8"))
        d = on_disk["anatomy"]
        print(f"  file           : {Path(sfile).name}")
        print(f"  named_by_hand  : {d.get('named_by_hand')}")
        if target not in (d.get("named_by_hand") or {}):
            fails.append("the declaration did not survive export to disk")
        if not d.get("names"):
            fails.append("the name map did not survive export to disk")
    else:
        # Not a failure of THIS behaviour -- sidecar() is what is under test and it
        # was checked above. Say so rather than passing silently.
        print(f"  cut is not exportable on this case ({str(s.problem)[:60]}), "
              f"so the on-disk half is untested")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
raise SystemExit(1 if fails else 0)
