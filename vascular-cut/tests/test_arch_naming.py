"""The brachiocephalic's own daughters keep their names.

    .venv\\Scripts\\python tests\\test_arch_naming.py
    .venv\\Scripts\\python tests\\test_arch_naming.py --prove

WHAT THIS GUARDS. Arch vessels are named two different ways, and only one of them is a
guess. A vessel is called a carotid or a subclavian because it climbs more than it
spreads, or spreads more than it climbs -- a pattern match, and a body carries plenty of
vessels that resemble a great one. But a vessel arising FROM THE BRACHIOCEPHALIC TRUNK
is named by where it branches, which is a fact about the tree.

When those two collided, both names were withdrawn as undistinctive. That cost three
cases their zone 1: on DJ8_1, KK20_2 and KY7_1 the trunk's own carotid daughter lost
"Right common carotid", a four-vessel elimination pass then handed it "Left common
carotid", and the left carotid came out sharing the right subclavian's station -- so
zone 1 spanned nothing and was never offered.

The invariant, and it needs no cohort to state: THE TWO CAROTIDS CANNOT SHARE A STATION
WITH THE TWO SUBCLAVIANS UNLESS THE ARCH IS BOVINE. On a bovine arch the left carotid
really does share the brachiocephalic origin, and then all three of the right carotid,
left carotid and right subclavian sit together. Anything between those two pictures --
the LEFT carotid sharing the RIGHT subclavian's station while the RIGHT carotid sits
alone further along -- is the swap, and describes no arch.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import branch_cut as BC                                      # noqa: E402

LOCAL = [("KY160_1", "KY160_1_CentrelinePoints.csv",
          "Geometries & Meshes/mesh/KY160_1.mat")]
COHORT = r"Z:/Martino/branches_CAMP/GT"
# The three the swap actually cost a zone, plus one normal and one truly bovine arch
# so the fix is shown not to have simply renamed everything.
COHORT_CASES = ("DJ8_1", "KK20_2", "KY7_1", "KR52_1", "KK5_1")
TOL_MM = 1.0

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print(f"    FAIL  {msg}")
    return cond


def cases():
    out = [(cid, ROOT / csv, ROOT / mat) for cid, csv, mat in LOCAL]
    try:
        from anat_bench import find_cases
        for cid, cl, mp in find_cases(COHORT):
            if cid in COHORT_CASES:
                out.append((cid, Path(cl), Path(mp)))
    except Exception:                                        # noqa: BLE001
        pass
    return out


def run():
    for cid, cl, mp in cases():
        s = CutSession(str(cl), str(mp))
        s.nub_preset()
        s.recut()
        route = BC.aorta_route(s.tree)
        jg = BC._junction_graph(s.tree)
        land = BC.aorta_landmark_stations(s.tree, route, jg)
        need = ("Right common carotid", "Left common carotid",
                "Right subclavian", "Left subclavian")
        got = {k: v for k, v in land.items() if k in need}
        print(f"\n{cid}: " + "  ".join(f"{k}={v:.1f}" for k, v in
                                       sorted(got.items(), key=lambda kv: kv[1])))
        if len(got) < 4:
            # Not this test's business. Whether all four get named at all is a
            # coverage question -- the repo's own KY160_1 names three -- and the
            # invariants below only say anything once all four exist.
            print(f"  only {len(got)} of the four named, nothing to check here")
            continue
        rcc, lcc = got["Right common carotid"], got["Left common carotid"]
        rsc, lsc = got["Right subclavian"], got["Left subclavian"]
        bovine = abs(lcc - rcc) <= TOL_MM and abs(lcc - rsc) <= TOL_MM

        # THE RIGHT PAIR SHARE AN ORIGIN. They are the brachiocephalic's daughters, so
        # on any arch, bovine or not, they leave at the same station.
        check(abs(rcc - rsc) <= TOL_MM,
              f"{cid}: the right carotid ({rcc:.1f}) and right subclavian "
              f"({rsc:.1f}) do not share the brachiocephalic origin")
        # THE LEFT SUBCLAVIAN IS LAST. It is the distal bound of zone 2 and the whole
        # zone map hangs off it.
        check(lsc >= lcc - TOL_MM,
              f"{cid}: the left subclavian ({lsc:.1f}) is proximal to the left "
              f"carotid ({lcc:.1f}), which no arch is")
        # AND THE LEFT CAROTID IS EITHER ON THE TRUNK OR PAST IT -- never sharing the
        # right subclavian's station while the right carotid sits further along, which
        # is the swap this guards.
        check(bovine or lcc > rsc + TOL_MM,
              f"{cid}: the left carotid shares the right subclavian's station "
              f"({lcc:.1f}) but the arch is not bovine -- the carotids are swapped, "
              f"and zone 1 will span nothing")
        print(f"  {'bovine arch' if bovine else 'three-vessel arch'}: "
              f"zone 1 spans {max(lcc - rcc, 0.0):.1f} mm")


def prove():
    """Put the blanket withdrawal back and check the suite notices."""
    import re
    src = (ROOT / "branch_cut.py").read_text(encoding="utf-8")
    marker = "            keep = [(src, lbl) for src, lbl in hits if lbl in prot]"
    if marker not in src:
        print("  cannot reinstate: _unique no longer protects anything")
        fails.append("prove: could not reinstate the blanket withdrawal")
        return
    patched = src.replace(marker,
                          "            keep = []  # PROVE: blanket withdrawal", 1)
    mod = ROOT / "branch_cut.py"
    original = src
    try:
        mod.write_text(patched, encoding="utf-8")
        for m in [k for k in list(sys.modules) if k in ("branch_cut", "vascular_api")]:
            del sys.modules[m]
        globals()["BC"] = __import__("branch_cut")
        globals()["CutSession"] = __import__("vascular_api").CutSession
        del fails[:]
        print(f"\n{'=' * 78}\nPROVE: blanket withdrawal reinstated")
        run()
    finally:
        mod.write_text(original, encoding="utf-8")
    if fails:
        print(f"  caught it: {len(fails)} assertion(s) fired")
        del fails[:]
    else:
        print("  NOT CAUGHT -- the test does not guard the withdrawal")
        fails.append("prove: blanket withdrawal went undetected")


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
