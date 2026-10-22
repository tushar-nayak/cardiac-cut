"""Work through the thoracic cohort case by case, starting from the arch-and-nub cut.

    python tools/arch_workbench.py                 # every case, one after another
    python tools/arch_workbench.py --list          # what is done and what is not
    python tools/arch_workbench.py --case KY160_1  # one case
    python tools/arch_workbench.py --flagged       # only the ones the batch flagged
    python tools/arch_workbench.py --redo KY160_1,KK20_2   # reopen ones already done

Each case opens in the ordinary desktop tool with the batch's cut ALREADY MADE -- the
arch isolated between its two stations, the great vessels left as nubs, the trims
added and the sliders where the preset put them. Adjust whatever needs adjusting,
press Export, close the window, and the next case opens. Everything the tool can
normally do is available; nothing about it changes.

Exports land in output/arch_final/ as <case>_arch.stl and <case>_arch_cut_planes.json.
A case counts as DONE when its file is there, so the flow survives being interrupted --
close the window without exporting and that case simply comes round again.

EVERY case is offered, including the ones the batch did not flag. "Nothing flagged"
means the things the run can measure came out right -- the rims lie on their planes,
the expected number of openings was made, the tool did not have to move a plane -- and
says nothing about whether the arch is bounded where an anatomist would bound it, or
whether a nub sits where it should on its artery. That is the judgement being made
here, so no case is skipped on its behalf. --flagged narrows it if wanted; the list
marks which is which either way.

One window per process, always. Two QtInteractor windows in one interpreter segfault
VTK, which is why this launches a child and waits rather than looping in-process.
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
DEFAULT_RUN = HERE / "arch_run.jsonl"
DEFAULT_CUTS = ROOT / "output" / "arch_nubs"
DEFAULT_FINAL = ROOT / "output" / "arch_final"
PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"


def load_cases(run, cuts_dir):
    """Every case in the run that has a sidecar to restore from."""
    out = []
    for line in Path(run).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        cid = rec["case"].split("/")[-1]
        side = Path(cuts_dir) / f"{cid}_arch_k{rec.get('k', 0.75):g}.json"
        if not side.exists():
            # A case the batch could not cut has no sidecar, so there is nothing to
            # restore -- but it is still worth opening, from the file as it loads,
            # because a person can often do what the batch could not.
            out.append((cid, rec, None))
            continue
        out.append((cid, rec, side))
    return out


def why(rec):
    """One line on what this case wants looking at."""
    if not rec.get("ok"):
        return f"BATCH REFUSED: {str(rec.get('error') or rec.get('message'))[:60]}"
    bits = []
    moved = rec.get("moved") or []
    if moved:
        names = ", ".join(str(m.get("name", "?"))[:18] for m in moved[:2])
        bits.append(f"{len(moved)} plane(s) moved by the tool ({names})")
    flat = rec.get("flat_mm")
    if flat is not None and flat > 0.5:
        bits.append(f"worst rim {flat:.2f} mm off plane")
    adjusted = [n for n in (rec.get("notes") or [])
                if "rescue" in n or "grown" in n or "clamped" in n]
    if adjusted:
        bits.append(f"{len(adjusted)} nub(s) adjusted")
    if rec.get("openings") != rec.get("expected"):
        bits.append(f"{rec.get('openings')} openings, {rec.get('expected')} asked for")
    return "; ".join(bits) or "nothing flagged"


def done(cid, final):
    return any((Path(final) / f"{cid}_arch.{e}").exists()
               for e in ("stl", "ply", "obj"))


def cmd_list(cases, final):
    print(f"{len(cases)} cases   exports -> {final}\n")
    n_done = n_flag = 0
    for cid, rec, side in cases:
        d = done(cid, final)
        f = why(rec) != "nothing flagged"
        n_done += d
        n_flag += f and not d
        mark = "DONE" if d else "todo"
        flag = " *" if f and not d else "  "
        miss = "" if side else "   (no batch cut, opens as it loads)"
        print(f"  [{mark}]{flag}{cid:12s} {why(rec)}{miss}")
    print(f"\n{n_done} exported, {len(cases) - n_done} still to do "
          f"({n_flag} of them flagged, marked *)")
    return 0


def open_case(cid, rec, side, final, wait=True):
    csv, mat = rec.get("csv"), rec.get("mat")
    if not csv or not mat:
        print(f"  {cid}: the run record has no source paths, re-run arch_cohort.py")
        return False
    exe = str(PYTHON) if PYTHON.exists() else sys.executable
    cmd = [exe, str(ROOT / "branch_gui.py"), csv, mat,
           "--out", str(final), "--name", f"{cid}_arch"]
    if side is not None:
        # The same switch the nub workbench uses. nub_session.apply restores the
        # trims, the sliders, the plane modes AND -- when the record carries one --
        # the region, which is what makes this the arch cut rather than a whole aorta
        # with nubs on it.
        cmd += ["--nubs", str(side)]
    print(f"\n=== {cid} ===  {why(rec)}")
    print(f"  export as {cid}_arch  ->  {final}")
    proc = subprocess.Popen(cmd, cwd=str(ROOT))
    if not wait:
        return True
    proc.wait()
    return done(cid, final)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=str(DEFAULT_RUN))
    ap.add_argument("--cuts-dir", default=str(DEFAULT_CUTS))
    ap.add_argument("--final", default=str(DEFAULT_FINAL))
    ap.add_argument("--case")
    ap.add_argument("--redo")
    ap.add_argument("--flagged", action="store_true")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    run = Path(args.run)
    if not run.exists():
        print(f"no run at {run}. Do this first:\n"
              f"    python tools/arch_cohort.py run")
        return 1
    final = Path(args.final)
    final.mkdir(parents=True, exist_ok=True)
    cases = load_cases(run, args.cuts_dir)
    if args.list:
        return cmd_list(cases, final)

    redo = {x.strip() for x in (args.redo or "").split(",") if x.strip()}
    if args.case:
        cases = [c for c in cases if c[0] == args.case]
        if not cases:
            print(f"no case {args.case} in {run}")
            return 1
    todo = [c for c in cases
            if (c[0] in redo) or (not done(c[0], final) and not redo)]
    if args.flagged:
        todo = [c for c in todo if why(c[1]) != "nothing flagged"]
    if not todo:
        print(f"nothing to do: all {len(cases)} exported to {final}")
        return 0

    print(f"{len(todo)} case(s) to review   exports -> {final}")
    print("Close each window when you are finished with it; the next one opens.\n")
    for i, (cid, rec, side) in enumerate(todo, 1):
        print(f"[{i}/{len(todo)}]", end="")
        ok = open_case(cid, rec, side, final)
        verdict = ("exported" if ok else
                   "closed without exporting, it will come round again")
        print(f"  -> {verdict}")
    left = [c for c in cases if not done(c[0], final)]
    print(f"\n{len(cases) - len(left)} of {len(cases)} exported to {final}")
    if left:
        print(f"still to do: {', '.join(c[0] for c in left)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
