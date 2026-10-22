"""A left carotid may not leave the aorta distal to the left subclavian.

    python tests/test_carotid_order.py

WHAT IS BEING PINNED

The naming rule separates the two by SHAPE -- a carotid ascends to the neck, a
subclavian reaches out to the arm, measured from takeoff to centerline tip. That is
sound when both vessels were traced to a comparable extent and wrong when they were
not: a subclavian whose centerline stops early rises more than it spreads, and reads
as a carotid. Nothing in the shape test can notice that it has produced an arch that
does not exist.

Order along the aorta can, and this is the invariant: on the LEFT the carotid leaves
before the subclavian, always. The consequences of getting it wrong are not confined
to a label. "Aortic arch" ends a fixed distance past the subclavian, so a swap put
the arch's distal bound 36 to 61 mm past it instead of 20; and Ishimaru zone 2 IS
carotid-to-subclavian, so with the names reversed its span was negative and the
preset was silently dropped as impossible.

TWO DIRECTIONS, because a guard that fires on everything is worse than none.

    KY68_1    wrong without the fix: the vessel named carotid left 21.9 mm DISTAL to
              the one named subclavian, and is the fatter of the two (r 4.58 against
              2.40), which a subclavian is.

    DJ18_2    right without it, and must stay right. The first version of this guard
              compared attach_idx -- an index into tree.junctions, whose order says
              nothing about position -- and swapped the cases that were already
              correct while leaving the wrong ones alone. Both names must survive,
              not merely the ordering: a guard that withdrew them would pass a test
              that only checked which came first.

Both live in the thoracic cohort rather than in this repository, so both are SKIPPED
when it is not mounted. KY160_1, which is in the repository, carries no left common
carotid at all -- the guard needs both vessels to have been named before it can have
an opinion -- so it is reported and not scored.

RIGHT SIDE DELIBERATELY UNTESTED, because the invariant does not hold there: an
aberrant right subclavian arises from the distal arch, past the carotid, and the
same guard applied to it would "correct" real anatomy into nonsense.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import branch_cut as BC                                              # noqa: E402
from vascular_api import CutSession                                  # noqa: E402

COHORT = Path(r"Z:/Martino/branches_CAMP/GT")
LOCAL = (ROOT / "KY160_1_CentrelinePoints.csv",
         ROOT / "Geometries & Meshes/mesh/KY160_1.mat")

fails = []


def stations(s):
    """Left carotid and left subclavian, as mm along the aorta. None if unnamed.

    Walked out to the route, because a great vessel may leave the aorta THROUGH
    another one -- on three of the four cases the guard fires on, the mislabelled
    subclavian hangs off the brachiocephalic, so its own junction is not an aortic
    junction and a plain lookup finds nothing.
    """
    jg = BC._junction_graph(s.tree)
    jst = {}
    for i, name in enumerate(s.route.chain):
        ns, ne = jg[name]
        seg = s.tree.segments[name]
        a, b = (ns, ne) if s.route.forward[i] else (ne, ns)
        if a is not None:
            jst.setdefault(a, s.route.starts[i])
        if b is not None:
            jst.setdefault(b, s.route.starts[i] + seg.length)

    def of(want):
        t = next((x for x in s.tree.terminals
                  if str(x.name).strip().lower() == want), None)
        if t is None:
            return None
        ns, ne = jg.get(t.segment.name, (None, None))
        for j in (ns, ne):
            st = BC._station_on_route(jg, jst, j)
            if st is not None:
                return float(st)
        return None

    return of("left common carotid"), of("left subclavian")


def check(tag, csv, mat, must_be_named=True):
    s = CutSession(str(csv), str(mat), family="thoracic")
    car, sub = stations(s)
    if car is None or sub is None:
        note = (f"{tag}: the left pair is not both named "
                f"(carotid {car}, subclavian {sub})")
        if must_be_named:
            fails.append(note + " -- the guard must not withdraw a name to "
                                "satisfy the ordering")
        else:
            print(f"  {note} -- the rule has no opinion here, not scored")
        return None
    print(f"  {tag}: left carotid at {car:.1f} mm, left subclavian at {sub:.1f} mm "
          f"({sub - car:+.1f})")
    if car >= sub:
        fails.append(f"{tag}: the left carotid leaves {car - sub:.1f} mm DISTAL to "
                     f"the left subclavian, which no arch does")
    return car, sub


print("a left carotid is proximal to the left subclavian\n")

csv, mat = LOCAL
if csv.exists() and mat.exists():
    check("KY160_1 (in the repository)", csv, mat, must_be_named=False)
else:
    print("  KY160_1: not in the repository, skipped")

WANT = {"KY68_1": "wrong without the guard",
        "DJ18_2": "right without it, and must stay right"}
seen = set()
if COHORT.exists():
    from tools.anat_bench import find_cases
    for cid, cl, mp in find_cases(str(COHORT)):
        short = cid.split("/")[-1]
        if short not in WANT:
            continue
        seen.add(short)
        check(f"{short} ({WANT[short]})", cl, mp)
for miss in sorted(set(WANT) - seen):
    print(f"  {miss}: thoracic cohort not mounted, skipped")

print()
if fails:
    for f in fails:
        print(f"FAIL  {f}")
    raise SystemExit(1)
print("PASS")
