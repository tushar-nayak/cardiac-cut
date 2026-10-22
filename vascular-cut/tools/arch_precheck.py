"""Open every case the way the workbench will, and say which ones arrive exportable.

    python tools/arch_precheck.py               # every case in the run
    python tools/arch_precheck.py --only KY160_1,KR6_1

The window refuses to export while its status bar is red, so a case that opens with a
complaint cannot be signed off until someone fixes it by hand. Worth knowing before
sitting down to a pass rather than one case at a time.

This is the only check that exercises what the reviewer will actually see: a real
CutWindow, loaded from the real files, with the batch's cut restored into it through
nub_session -- the same path, including the region restore that makes this the ARCH
cut rather than a whole aorta with nubs on it. The batch scoring its own mesh cannot
see any of that, because it never opens a window.

It also checks the one thing that would be silent and ruinous: that the region really
came back. A record whose region failed to restore still opens, still cuts, still
exports -- as the whole thoracic aorta. Nobody would notice until the meshes were
measured.

One case per process. Two QtInteractor windows in one interpreter segfault VTK, which
is why this shells out rather than looping.
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
PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"

CHILD = r'''
import json, sys
from pathlib import Path
ROOT = Path(sys.argv[1])
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tools"))
from qtpy import QtWidgets
import branch_cut as BC, nub_session
from branch_gui import CutWindow
csv, mat, side = sys.argv[2], sys.argv[3], sys.argv[4]
app = QtWidgets.QApplication([])
w = CutWindow(csv, mat); w.show()
for _ in range(6): app.processEvents()
notes = nub_session.apply_from_file(w, side)
for _ in range(12): app.processEvents()
print("@@" + json.dumps({
    "openings": len(BC.boundary_loops(w.cut)),
    "euler": int(w.cut.euler_number),
    "rows": len(w.rows),
    "trims": len(w.trims),
    "bounds": len(getattr(w, "_region_bounds", []) or []),
    "region": None if w.region is None else [
        None if x is None else round(float(x), 2) for x in w.region],
    "region_name": getattr(w, "_region_name", "") or "",
    "problem": None if w.problem is None else str(w.problem),
    "notes": list(notes),
}))
'''


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=str(HERE / "arch_run.jsonl"))
    ap.add_argument("--cuts-dir", default=str(ROOT / "output" / "arch_nubs"))
    ap.add_argument("--only")
    args = ap.parse_args()

    run = Path(args.run)
    if not run.exists():
        print(f"no run at {run}. Do this first:\n    python tools/arch_cohort.py run")
        return 1
    recs = [json.loads(l) for l in run.read_text(encoding="utf-8").splitlines()
            if l.strip()]
    if args.only:
        want = {o.strip() for o in args.only.split(",")}
        recs = [r for r in recs if r["case"].split("/")[-1] in want]

    exe = str(PYTHON) if PYTHON.exists() else sys.executable
    script = Path(args.cuts_dir) / "_arch_precheck_child.py"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(CHILD, encoding="utf-8")

    good, bad, lost = [], [], []
    for i, rec in enumerate(recs, 1):
        cid = rec["case"].split("/")[-1]
        side = Path(args.cuts_dir) / f"{cid}_arch_k{rec.get('k', 0.75):g}.json"
        if not (side.exists() and rec.get("csv") and rec.get("mat")):
            print(f"  [{i}/{len(recs)}] {cid:12s} SKIPPED, the batch wrote no sidecar")
            continue
        p = subprocess.run([exe, str(script), str(ROOT), rec["csv"], rec["mat"],
                            str(side)], capture_output=True, text=True, timeout=900)
        line = next((l for l in (p.stdout or "").splitlines()
                     if l.startswith("@@")), None)
        if line is None:
            bad.append(cid)
            tail = (p.stderr or p.stdout or "").strip().splitlines()[-1:] or [""]
            print(f"  [{i}/{len(recs)}] {cid:12s} CHILD FAILED  {tail[0][:60]}")
            continue
        d = json.loads(line[2:])
        # The region is the thing that would go missing quietly.
        has_region = d["region"] is not None and d["bounds"] > 0
        if not has_region:
            lost.append(cid)
        state = ("exportable" if d["problem"] is None else "REFUSES")
        print(f"  [{i}/{len(recs)}] {cid:12s} {state:11s} "
              f"{d['openings']} openings, {d['rows']} rows, {d['trims']} trim(s), "
              f"{d['bounds']} bound(s)"
              + ("" if has_region else "   <<< REGION DID NOT RESTORE")
              + ("" if d["problem"] is None else f"   {d['problem'][:44]}"))
        (good if d["problem"] is None else bad).append(cid)

    print(f"\n{'=' * 70}\n  {len(good)} open exportable, {len(bad)} arrive with a "
          f"complaint\n{'=' * 70}")
    if lost:
        print(f"  REGION MISSING on {len(lost)}: {', '.join(lost)}")
        print("  Those would export the whole aorta with nubs, not the arch.")
    if bad:
        print(f"  need a hand: {', '.join(bad)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
