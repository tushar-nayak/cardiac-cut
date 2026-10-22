"""The aortic root is found from where the branches leave, even through a stem.

    .venv\\Scripts\\python tests\\test_proximal_end.py
    .venv\\Scripts\\python tests\\test_proximal_end.py --prove

WHAT THIS GUARDS. _proximal_end decides which trunk end the blood enters by, and its
first and best signal is that the great vessels leave from the ASCENDING end. It used to
count a branch only when the branch's OWN segment touched the aorta. Most arches put the
right subclavian and right carotid on the brachiocephalic, so those two never count; and
a centerline that draws the left carotid and left subclavian off a short common stem
leaves no branch touching the aorta at all. With no takeoffs it fell through to the
taper rule -- and on ZF_Example's KK16_6 (_2 centerline, 2 Oct) the descending end is the
fatter one (r 21.6 against 14.6 mm), so the descending aorta was named the root, zone 0
ran 316 mm down it, and zone 2 vanished.

The synthetic trees below are that centerline's topology and segment lengths, plus the
same case's earlier centerline (whose carotid and subclavian touch the aorta directly,
and which was always read correctly), so the check needs no patient data. When the
ZF_Example folder is present the real files are checked too.

Two more guard the fix itself. A branch group with three tips low on the descending
aorta must not outvote the arch: the count is one vote per TAKEOFF, not per tip (the
first version of this fix counted tips, and a review showed that flips a root the old
rule had right). And a branch on a loop that reaches the aorta at two junctions is
skipped, so the answer does not depend on the order of segments in the file.

--prove puts each earlier count back (direct-only, and one vote per tip), in memory,
and checks the suite notices both.
"""
import sys
from pathlib import Path
from types import SimpleNamespace as NS

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import branch_cut as BC                                      # noqa: E402

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print(f"    FAIL  {msg}")
    else:
        print(f"    ok    {msg}")
    return cond


def tree_from(lengths, jg, terminals):
    """A stand-in tree with just what _proximal_end reads."""
    segs = {n: NS(name=n, length=L) for n, L in lengths.items()}
    terms = [NS(label=lbl, segment=segs[seg], tip=(0.0, 0.0, z)) for lbl, seg, z in terminals]
    return NS(segments=segs, terminals=terms), {t.label: t for t in terms}


# ZF_Example KK16_6, _2 centerline: Line13 is the ascending aorta. The brachiocephalic
# (Line16) carries the right subclavian and right carotid; a 10 mm stem (Line04) carries
# the left carotid and left subclavian. No branch segment touches the aorta.
STEM = tree_from(
    {"Line19": 316.5, "Line12": 34.8, "Line13": 53.7, "Line16": 79.2, "Line18": 131.9,
     "Line15": 54.3, "Line04": 10.0, "Line08": 50.1, "Line17": 94.9},
    None, [("Line19.start", "Line19", -326.1), ("Line13.start", "Line13", -181.1),
           ("Line15.end", "Line15", -17.4), ("Line18.end", "Line18", -19.3),
           ("Line17.end", "Line17", -19.4), ("Line08.end", "Line08", -63.3)])
STEM_JG = {"Line16": (0, 1), "Line04": (2, 3), "Line12": (0, 2), "Line13": (None, 0),
           "Line08": (3, None), "Line18": (1, None), "Line19": (None, 2),
           "Line17": (3, None), "Line15": (1, None)}
STEM_RAD = {"Line19.start": 21.55, "Line13.start": 14.63}

# The same patient's earlier centerline: the left carotid (Line67) and left subclavian
# (Line66) leave the aorta directly. Line29 is the ascending aorta.
DIRECT = tree_from(
    {"Line76": 302.4, "Line29": 57.7, "Line52": 26.0, "Line66": 55.1, "Line72": 78.1,
     "Line73": 54.3, "Line75": 138.5, "Line67": 112.4, "Line33": 22.0},
    None, [("Line29.end", "Line29", -181.1), ("Line76.start", "Line76", -326.1),
           ("Line73.end", "Line73", -16.7), ("Line67.end", "Line67", -18.0),
           ("Line75.end", "Line75", -20.2), ("Line66.end", "Line66", -63.3)])
DIRECT_JG = {"Line76": (None, 0), "Line29": (1, None), "Line52": (1, 2), "Line66": (0, None),
             "Line72": (3, 1), "Line73": (3, None), "Line75": (3, None), "Line67": (2, None),
             "Line33": (2, 0)}
DIRECT_RAD = {"Line29.end": 14.87, "Line76.start": 21.57}


# The earlier centerline again, with a celiac-like group low on the descending aorta: an
# interior trunk (Cel) 30 mm above the distal end, ending in three tips. One vote per TIP
# puts the mean at the middle of the aorta, the taper rule then decides, and the fatter
# descending end wins; one vote per TAKEOFF keeps the ascending end.
LOW = tree_from(
    {"Line76a": 30.0, "Line76b": 272.4, "Line29": 57.7, "Line52": 26.0, "Line66": 55.1,
     "Line72": 78.1, "Line73": 54.3, "Line75": 138.5, "Line67": 112.4, "Line33": 22.0,
     "Cel": 15.0, "C1": 20.0, "C2": 20.0, "C3": 20.0},
    None, [("Line29.end", "Line29", -181.1), ("Line76a.start", "Line76a", -326.1),
           ("Line73.end", "Line73", -16.7), ("Line67.end", "Line67", -18.0),
           ("Line75.end", "Line75", -20.2), ("Line66.end", "Line66", -63.3),
           ("C1.end", "C1", -300.0), ("C2.end", "C2", -301.0), ("C3.end", "C3", -302.0)])
LOW_JG = {"Line76a": (None, 9), "Line76b": (9, 0), "Line29": (1, None), "Line52": (1, 2),
          "Line66": (0, None), "Line72": (3, 1), "Line73": (3, None), "Line75": (3, None),
          "Line67": (2, None), "Line33": (2, 0), "Cel": (9, 10), "C1": (10, None),
          "C2": (10, None), "C3": (10, None)}
LOW_RAD = {"Line29.end": 14.87, "Line76a.start": 21.57}


def run_synthetic():
    (tree, t), jg = STEM, STEM_JG
    for order in (("Line19.start", "Line13.start"), ("Line13.start", "Line19.start")):
        trunk = [t[order[0]], t[order[1]]]
        got, how = BC._proximal_end(tree, trunk, jg, STEM_RAD)
        check(got.label == "Line13.start" and how == "branches",
              f"stem centerline, trunk order {order[0][:6]} first: root is the ascending "
              f"end Line13.start by branches (got {got.label} by {how})")

    (tree, t), jg = DIRECT, DIRECT_JG
    for order in (("Line29.end", "Line76.start"), ("Line76.start", "Line29.end")):
        trunk = [t[order[0]], t[order[1]]]
        got, how = BC._proximal_end(tree, trunk, jg, DIRECT_RAD)
        check(got.label == "Line29.end" and how == "branches",
              f"direct centerline, trunk order {order[0][:6]} first: root unchanged, "
              f"Line29.end by branches (got {got.label} by {how})")

    (tree, t), jg = LOW, LOW_JG
    for order in (("Line29.end", "Line76a.start"), ("Line76a.start", "Line29.end")):
        trunk = [t[order[0]], t[order[1]]]
        got, how = BC._proximal_end(tree, trunk, jg, LOW_RAD)
        check(got.label == "Line29.end" and how == "branches",
              f"three-tip group low on the descending aorta, trunk order {order[0][:7]} first: "
              f"one vote per takeoff keeps the ascending end (got {got.label} by {how})")

    # A loop: branch X sits on a path that joins the aorta at BOTH arch junctions. It has
    # no single takeoff, so it is skipped -- and listing the loop's segments in either
    # order must give the same answer.
    (tree, t) = STEM
    base = dict(STEM_JG)
    trunk = [t["Line19.start"], t["Line13.start"]]
    chain = BC._aorta_chain(tree, base, "Line19", "Line13")
    stations = BC._aorta_stations(tree, base, "Line19", "Line13")
    plain = BC._branch_takeoffs(tree, base, trunk, chain, stations)
    for first, second in (("P1", "P2"), ("P2", "P1")):
        loop_tree, lt = tree_from(
            {**{n: s.length for n, s in tree.segments.items()}, "P1": 20.0, "P2": 20.0, "X": 15.0},
            None, [(x.label, x.segment.name, x.tip[2]) for x in tree.terminals]
            + [("X.end", "X", -40.0)])
        ends = {"P1": (0, 20), "P2": (20, 2)}
        jg = dict(base)
        jg[first], jg[second], jg["X"] = ends[first], ends[second], (20, None)
        ltrunk = [lt["Line19.start"], lt["Line13.start"]]
        got = BC._branch_takeoffs(loop_tree, jg, ltrunk, chain, stations)
        check(got == plain,
              f"loop listed {first} first: the looping branch is skipped, takeoffs unchanged "
              f"({got} against {plain})")

    # No branches at all: nothing to count, so the taper rule still decides.
    tree, t = tree_from({"A": 100.0}, None, [("A.start", "A", 0.0), ("A.end", "A", -50.0)])
    got, how = BC._proximal_end(tree, [t["A.start"], t["A.end"]], {"A": (None, None)},
                                {"A.start": 10.0, "A.end": 14.0})
    check(got.label == "A.end" and how == "taper",
          f"unbranched vessel: the fatter end by taper, as before (got {got.label} by {how})")


def run_real():
    ex = ROOT / "ZF_Example"
    cl = ex / "Zihan_branches_KK16_6_CentrelinePoints_2.csv"
    mp = ex / "Zihan_branches_KK16_6_2.mat"
    if not (cl.exists() and mp.exists()):
        print("  ZF_Example KK16_6 _2 files not present; real-data check skipped")
        return
    from vascular_api import CutSession
    s = CutSession(str(cl), str(mp))
    names = s.anat_names or {}
    check(names.get("Line13.start") == "Ascending aorta",
          f"ZF_Example KK16_6 _2: Line13.start named the ascending aorta "
          f"(root by {s.tree.proximal_from})")
    zone0 = next((p for p in s.presets() if p["name"] == "Ishimaru zone 0"), None)
    check(zone0 is not None and zone0["s1"] < 80.0,
          f"ZF_Example KK16_6 _2: zone 0 ends at the brachiocephalic, not 316 mm down "
          f"the descending aorta (s1={None if zone0 is None else round(zone0['s1'], 1)})")


def run():
    print("synthetic trees")
    run_synthetic()
    print("real data")
    run_real()


def direct_only(tree, jg, trunk, chain, stations):
    """The original count: a branch only when its own segment touches the aorta."""
    out = []
    for t in tree.terminals:
        if t in trunk:
            continue
        ns, ne = jg[t.segment.name]
        j = ns if ns is not None else ne
        if j in stations:
            out.append(stations[j])
    return out


def per_tip(tree, jg, trunk, chain, stations):
    """The first version of the fix: through stems, but one vote per tip, nearest by hops."""
    link = {}
    for name, (ns, ne) in jg.items():
        if name in set(chain) or ns is None or ne is None:
            continue
        link.setdefault(ns, []).append(ne)
        link.setdefault(ne, []).append(ns)
    out = []
    for t in tree.terminals:
        if t in trunk:
            continue
        ns, ne = jg[t.segment.name]
        start = ns if ns is not None else ne
        seen, frontier, hit = {start}, [start], None
        while frontier and hit is None:
            nxt = []
            for j in frontier:
                if j in stations:
                    hit = stations[j]
                    break
                for o in link.get(j, ()):
                    if o not in seen:
                        seen.add(o)
                        nxt.append(o)
            frontier = nxt
        if hit is not None:
            out.append(hit)
    return out


def prove():
    """Put each earlier takeoff count back, in memory, and check the suite notices."""
    if not hasattr(BC, "_branch_takeoffs"):
        print("  nothing to reinstate: _branch_takeoffs does not exist (the fix is not in)")
        fails.append("prove: _branch_takeoffs missing")
        return
    real = BC._branch_takeoffs
    for name, impostor in (("direct-only count", direct_only), ("one vote per tip", per_tip)):
        BC._branch_takeoffs = impostor
        try:
            del fails[:]
            print(f"{'=' * 78}\nPROVE: {name} reinstated")
            run_synthetic()
        finally:
            BC._branch_takeoffs = real
        if fails:
            print(f"  caught it: {len(fails)} assertion(s) fired")
            del fails[:]
        else:
            print(f"  NOT CAUGHT -- the test does not guard against the {name}")
            fails.append(f"prove: {name} went undetected")


if __name__ == "__main__":
    if "--prove" in sys.argv:
        prove()
    else:
        run()
    print("\n" + "=" * 78)
    if fails:
        print(f"{len(fails)} FAILURE(S)")
        for f in fails:
            print(f"  {f}")
        sys.exit(1)
    print("all checks passed")
