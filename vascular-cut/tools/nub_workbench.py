"""Work through the thoracic cohort case by case, starting from the nub-in cut.

    python tools/nub_workbench.py                 # every case, one after another
    python tools/nub_workbench.py --list          # what is done and what is not
    python tools/nub_workbench.py --case KY160_1  # one case
    python tools/nub_workbench.py --flagged       # only the ones the batch flagged
    python tools/nub_workbench.py --redo KY160_1,KY215_1   # reopen ones already done

Each case opens in the ordinary desktop tool with the batch's nub cut ALREADY MADE
-- the trims added, the sliders where the preset put them, the aorta's ends squared
off. Adjust whatever needs adjusting, press Export, close the window, and the next
case opens. Everything the tool can normally do is available; nothing about it
changes.

Exports land in output/nubs_final/ as <case>_nub.stl and <case>_nub_cut_planes.json.
A case counts as DONE when its file is there, so the flow survives being interrupted
-- close the window without exporting and that case simply comes round again.

EVERY case is offered, including the ones the batch did not flag. "Nothing flagged"
means the three things the run can measure came out right -- the rim lies on its
plane, one opening was made, no other vessel was touched -- and says nothing about
whether the nub is in the right place on that artery. That is the judgement being
made here, so no case is skipped on its behalf. --flagged narrows it if wanted; the
list marks which is which either way.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

HERE = Path(__file__).resolve().parent
DEFAULT_RUN = HERE / "nub_k075.jsonl"
DEFAULT_NUBS = ROOT / "output" / "nubs_k0.75"
DEFAULT_FINAL = ROOT / "output" / "nubs_final"
PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"


def _needs_a_look(rec):
    """The same test the reviewer uses, so the two agree on which cases are flagged."""
    sys.path.insert(0, str(ROOT))
    import nub_viewer
    return nub_viewer._needs_a_look(rec)


def load_cases(run, nubs_dir):
    out = []
    for line in Path(run).read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        cid = rec["case"].split("/")[-1]
        side = Path(nubs_dir) / f"{cid}_nub_k{rec['k']:g}.json"
        if not side.exists():
            continue
        out.append((cid, rec, side))
    return out


def why(rec):
    """One line on what this case wants looking at."""
    bits = []
    grown = [v for v in rec["vessels"] if v.get("rescued_mm") or v.get("retried_mm")]
    clamped = [v for v in rec["vessels"] if v.get("clamped")]
    marginal = [v for v in rec["vessels"] if v.get("marginal")]
    torn = [v for v in rec["vessels"] if v.get("ok") and not v.get("cut_ok")]
    uncut = [v for v in rec["vessels"]
             if not v.get("ok") and "no offset" in (v.get("why") or "")]
    handed = [v for v in rec["vessels"] if v.get("handed_to")]
    closed = [n for n in rec.get("notes", []) if "left closed" in n]
    perp = [n for n in rec.get("notes", []) if "perpendicular to the" in n]
    if uncut:
        bits.append(f"{len(uncut)} vessel(s) UNCUT")
    if torn:
        bits.append(f"{len(torn)} torn")
    if marginal:
        bits.append(f"{len(marginal)} barely a nub")
    if clamped:
        bits.append(f"{len(clamped)} shortened")
    if grown:
        worst = max(v.get("rescued_mm") or 0 for v in grown)
        bits.append(f"{len(grown)} grown (up to +{worst:.1f} mm)")
    if handed:
        bits.append(f"{len(handed)} trunk(s) handed to daughters")
    if closed:
        bits.append(f"{len(closed)} aortic end(s) left closed")
    if perp:
        bits.append(f"{len(perp)} end(s) cut perpendicular, too oblique for z")
    return "; ".join(bits) or "nothing flagged"


def done(cid, final):
    return any((Path(final) / f"{cid}_nub.{e}").exists() for e in ("stl", "ply", "obj"))


def cmd_list(cases, final):
    print(f"{len(cases)} cases   exports -> {final}\n")
    n_done = n_flag = 0
    for cid, rec, _ in cases:
        d = done(cid, final)
        f = _needs_a_look(rec)
        n_done += d
        n_flag += f and not d
        mark = "DONE" if d else "todo"
        flag = " *" if f and not d else "  "
        print(f"  [{mark}]{flag}{cid:12s} {rec['openings'] or 0} openings   {why(rec)}")
    print(f"\n{n_done} exported, {len(cases) - n_done} still to do "
          f"({n_flag} of them flagged, marked *)")
    return 0


def open_case(cid, rec, side, final, wait=True):
    csv, mat = rec.get("csv"), rec.get("mat")
    if not csv or not mat:
        print(f"  {cid}: the run record has no source paths, re-run nub_cohort.py")
        return False
    exe = str(PYTHON) if PYTHON.exists() else sys.executable
    cmd = [exe, str(ROOT / "branch_gui.py"), csv, mat,
           "--nubs", str(side), "--out", str(final), "--name", f"{cid}_nub"]
    print(f"\n=== {cid} ===  {why(rec)}")
    print(f"  export as {cid}_nub  ->  {final}")
    proc = subprocess.Popen(cmd, cwd=str(ROOT))
    if not wait:
        return True
    proc.wait()
    return done(cid, final)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=str(DEFAULT_RUN))
    ap.add_argument("--nubs-dir", default=str(DEFAULT_NUBS))
    ap.add_argument("--final", default=str(DEFAULT_FINAL))
    ap.add_argument("--case")
    ap.add_argument("--redo")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--flagged", action="store_true",
                    help="only the cases the batch flagged, instead of all of them")
    ap.add_argument("--all", action="store_true",
                    help="every case (the default; accepted so it is never wrong)")
    args = ap.parse_args()

    final = Path(args.final)
    final.mkdir(parents=True, exist_ok=True)
    cases = load_cases(args.run, args.nubs_dir)
    if not cases:
        raise SystemExit(f"no cases in {args.run} with sidecars in {args.nubs_dir}")

    if args.list:
        return cmd_list(cases, final)

    one = args.case or args.redo
    if one:
        # A comma-separated list, because redoing a rule change means redoing the
        # handful of cases it moved, and doing that one command at a time is a way to
        # miss one.
        want = [w.strip() for w in str(one).split(",") if w.strip()]
        picked = []
        for w in want:
            hit = next((c for c in cases if c[0] == w), None)
            if hit is None:
                raise SystemExit(f"no case {w}; --list shows what there is")
            picked.append(hit)
        for i, (cid, rec, side) in enumerate(picked, 1):
            print(f"--- {i} of {len(picked)} ---")
            open_case(cid, rec, side, final)
        return cmd_list(cases, final) if len(picked) > 1 else 0

    queue = [c for c in cases
             if not done(c[0], final) and (not args.flagged or _needs_a_look(c[1]))]
    if not queue:
        print("Nothing left to work through. --redo <case> reopens one.")
        return cmd_list(cases, final)

    print(f"{len(queue)} case(s) to work through. Export and close each window to "
          f"move on; close without exporting to leave a case for later, and "
          f"Ctrl-C here to stop.\n")
    for i, (cid, rec, side) in enumerate(queue, 1):
        print(f"--- {i} of {len(queue)} ---")
        try:
            ok = open_case(cid, rec, side, final)
        except KeyboardInterrupt:
            print("\nstopped")
            break
        print(f"  {cid}: {'exported' if ok else 'left for later'}")
    return cmd_list(cases, final)


if __name__ == "__main__":
    sys.exit(main())
