"""Open every case the way the workbench will, and say which ones arrive exportable.

    python tools/nub_precheck.py            # all 30
    python tools/nub_precheck.py --only DJ18_2,KY238_1

The window refuses to export while its status bar is red, so a case that opens with a
complaint cannot be signed off until someone fixes it by hand. Worth knowing before
sitting down to a pass rather than one case at a time.

This is the only check that exercises what the reviewer will actually see: a real
CutWindow, loaded from the real files, with the batch cut restored into it through
nub_session -- the same path, including the split-line repair and the cache it has to
invalidate. The batch scoring its own mesh cannot see any of that.

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
for _ in range(10): app.processEvents()
print("@@" + json.dumps({
    "openings": len(BC.boundary_loops(w.cut)),
    "euler": int(w.cut.euler_number),
    "rows": len(w.rows),
    "trims": len(w.trims),
    "problem": None if w.problem is None else str(w.problem),
    "status": w.status.text(),
    "notes": [n for n in notes],
}))
'''


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=str(HERE / "nub_k075.jsonl"))
    ap.add_argument("--nubs-dir", default=str(ROOT / "output" / "nubs_k0.75"))
    ap.add_argument("--only")
    args = ap.parse_args()

    recs = [json.loads(l) for l in Path(args.run).read_text().splitlines() if l.strip()]
    if args.only:
        want = {o.strip() for o in args.only.split(",")}
        recs = [r for r in recs if r["case"].split("/")[-1] in want]

    exe = str(PYTHON) if PYTHON.exists() else sys.executable
    script = Path(args.nubs_dir) / "_precheck_child.py"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(CHILD)

    good, bad = [], []
    for i, rec in enumerate(recs, 1):
        cid = rec["case"].split("/")[-1]
        side = Path(args.nubs_dir) / f"{cid}_nub_k{rec['k']:g}.json"
        if not (side.exists() and rec.get("csv") and rec.get("mat")):
            print(f"  [{i}/{len(recs)}] {cid:12s} SKIPPED, no sidecar or source paths")
            continue
        out = subprocess.run([exe, str(script), str(ROOT), rec["csv"], rec["mat"],
                              str(side)], capture_output=True, text=True,
                             cwd=str(ROOT))
        line = next((l for l in out.stdout.splitlines() if l.startswith("@@")), None)
        if line is None:
            print(f"  [{i}/{len(recs)}] {cid:12s} FAILED TO OPEN")
            bad.append((cid, (out.stderr or out.stdout)[-160:].strip()))
            continue
        d = json.loads(line[2:])
        if d["problem"] is None:
            good.append(cid)
            print(f"  [{i}/{len(recs)}] {cid:12s} ready   {d['openings']} openings, "
                  f"euler {d['euler']}, {d['rows']} rows, {d['trims']} trim(s)")
        else:
            bad.append((cid, d["problem"]))
            print(f"  [{i}/{len(recs)}] {cid:12s} NOT EXPORTABLE  {d['problem'][:96]}")
    script.unlink(missing_ok=True)

    print(f"\n{len(good)} of {len(good) + len(bad)} open ready to export")
    if bad:
        print("\nthese need a hand before they can be exported:")
        for cid, why in bad:
            print(f"  {cid:12s} {why[:110]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
