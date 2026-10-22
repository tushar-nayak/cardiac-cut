"""Does the SESSION recover the cuts that raw cut_terminal lost against main?

tools/cut_bench.py probes `cut_terminal` directly, one terminal at a time, at fixed
fractions of its segment. Against `main` this branch scores ALL clean 296 -> 276 and
lists 20 cuts clean on main and not now -- and every one of the 20 is a TRUNK cut
("Aorta root (proximal)" / "Aorta distal") on a GT case.

That pattern has an obvious candidate explanation: this branch added
`flat_trunk_offset`, which main has no equivalent of. A session detects a trunk cut
that comes out badly and backs it off to where the vessel is flat, recording the move
in the sidecar as `adjustments.vessel_end`. cut_bench cannot see that, because it
never builds a session.

So the benchmark delta is only a regression if the SESSION also fails. This asks
exactly that, on exactly those 20 cuts:

    raw   BC.cut_terminal at the probed offset            (what cut_bench measures)
    sess  CutSession with that terminal at the same offset (what a person gets)

and reports, for each, whether the session produced a valid cut, whether it had to
move the plane to do it, and how far.

    .venv\\Scripts\\python tools\\session_vs_raw.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import branch_cut as BC                                             # noqa: E402
import cut_bench as CB                                              # noqa: E402
from vascular_api import CutSession                                 # noqa: E402

# The 20, verbatim from the vs-main run: (case, terminal name, fraction)
PROBES = [
    ("GT/DJ8_1", "Aorta root (proximal)", 0.70),
    ("GT/DJ8_1", "Aorta root (proximal)", 0.80),
    ("GT/DJ8_1", "Aorta root (proximal)", 0.85),
    ("GT/DJ8_1", "Aorta distal", 0.90),
    ("GT/KK15_1", "Aorta distal", 0.80),
    ("GT/KK15_1", "Aorta distal", 0.85),
    ("GT/KK15_1", "Aorta distal", 0.90),
    ("GT/KK16_1", "Aorta distal", 0.85),
    ("GT/KK16_1", "Aorta distal", 0.90),
    ("GT/KK20_2", "Aorta root (proximal)", 0.80),
    ("GT/KK20_2", "Aorta root (proximal)", 0.85),
    ("GT/KK20_2", "Aorta root (proximal)", 0.90),
    ("GT/KK32_0", "Aorta root (proximal)", 0.50),
    ("GT/KK32_0", "Aorta root (proximal)", 0.70),
    ("GT/KK32_0", "Aorta distal", 0.20),
    ("GT/KK32_0", "Aorta distal", 0.50),
    ("GT/KK32_0", "Aorta distal", 0.80),
    ("GT/KK32_0", "Aorta distal", 0.85),
    ("GT/KK32_0", "Aorta distal", 0.90),
    ("GT/KK4_2", "Aorta distal", 0.20),
]


def main():
    paths = {n: (c, m) for n, c, m in CB.find_cases()}
    missing = sorted({c for c, _, _ in PROBES} - set(paths))
    if missing:
        print(f"cases not found: {missing}")

    rows = []
    cache = {}
    for case, term, frac in PROBES:
        if case not in paths:
            continue
        csv, mat = paths[case]
        if case not in cache:
            cache[case] = CB.build(BC, csv, mat)
        mesh, tree = cache[case]
        t = next((x for x in tree.terminals if x.name == term), None)
        if t is None:
            rows.append((case, term, frac, "no such terminal", "", "", ""))
            continue
        off = round(frac * float(t.segment.length), 1)

        # RAW: what cut_bench measures.
        try:
            BC.cut_terminal(mesh, tree, t, off, snap_mm=BC.snap_tolerance(mesh))
            dev = getattr(t, "_rim_dev_mm", None)
            raw = "clean" if (dev is not None and dev <= BC.RIM_FLAW_MM) else \
                  (f"flawed {dev:.1f}mm" if dev is not None else "cut")
        except Exception as e:                                      # noqa: BLE001
            raw = f"REFUSED ({str(e).split(',')[0][:28]})"

        # SESSION: what a person gets, flat_trunk_offset and all.
        #
        # Matched by LABEL, never by name. A CutSession overwrites t.name with the
        # ANATOMICAL name ("Ascending aorta"), while cut_bench's tree keeps the
        # classifier's ("Aorta root (proximal)") -- so a name lookup finds nothing on
        # every single row. It did, and because a miss was counted the same as a
        # failure the whole run reported 0 of 20 recovered and a verdict of "real
        # regression". Labels are the stable identity; cut_bench itself matches on
        # them for the same reason, and says so.
        s = CutSession(str(csv), str(mat))
        st_t = next((x for x in s.tree.terminals if x.label == t.label), None)
        if st_t is None:
            # Not a failure of the cutter -- a failure to compare. Marked so it
            # cannot be silently counted as either outcome.
            rows.append((case, term, frac, raw, "UNCOMPARABLE (no such label)",
                         "", ""))
            continue
        s.set_offset(st_t.label, off)
        s.recut()
        moved = getattr(s, "_trunk_moved", None)
        final = None
        for o in s.openings():
            if o["label"] == st_t.label:
                final = o["offset_mm"]
        sess = "VALID" if s.ok else f"invalid ({str(s.problem)[:34]})"
        movetxt = ("-" if not moved or moved.get("label") != st_t.label
                   else f"{moved['from_mm']:.1f}->{moved['to_mm']:.1f}")
        rows.append((case, term, frac, raw, sess, movetxt,
                     f"{final:.1f}" if final is not None else "?"))

    print("=" * 104)
    print(f"{'case':<12}{'terminal':<24}{'frac':<6}{'RAW cut_terminal':<26}"
          f"{'SESSION':<22}{'moved':<14}{'final'}")
    print("=" * 104)
    for r in rows:
        print(f"{r[0]:<12}{r[1]:<24}{r[2]:<6}{r[3]:<26}{r[4]:<22}{r[5]:<14}{r[6]}")

    # Anything that could not be compared is excluded from the denominator rather
    # than counted against the cutter. Counting it as a failure is what turned a
    # broken lookup into a "real regression" verdict on the first run.
    bad = [r for r in rows if str(r[4]).startswith("UNCOMPARABLE")]
    rows = [r for r in rows if not str(r[4]).startswith("UNCOMPARABLE")]
    n = len(rows)
    ok = sum(1 for r in rows if str(r[4]).startswith("VALID"))
    mv = sum(1 for r in rows if r[5] not in ("-", ""))
    print()
    if bad:
        print(f"  {len(bad)} probe(s) could not be compared and are EXCLUDED, not "
              f"counted as failures:")
        for r in bad:
            print(f"    {r[0]} {r[1]} {r[2]}")
        print()
    if not n:
        print("  nothing comparable -- no verdict.")
        return
    print(f"  of {n} cuts that raw cut_terminal loses against main:")
    print(f"    the session produces a VALID cut for   {ok}")
    print(f"    ...of which it moved the plane to do it {mv}")
    print()
    if ok == n:
        print("  VERDICT: not a regression in what a person gets. Every one is")
        print("  recovered by the session layer, which is exactly what")
        print("  flat_trunk_offset was added for and what cut_bench cannot see.")
    elif ok:
        print(f"  VERDICT: PARTIAL. {n - ok} are not recovered by the session either,")
        print("  and those are real losses against main. Listed above as 'invalid'.")
    else:
        print("  VERDICT: NOT an artefact. The session recovers none of them, so the")
        print("  benchmark delta is a real regression against main. Investigate before")
        print("  merging.")


if __name__ == "__main__":
    main()
