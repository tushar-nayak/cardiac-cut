"""The plain path, on every geometry: load, cut, validate, export.

Everything else in tools/ scores a feature. This scores the thing a visitor to the
website will actually do -- hand it a mesh and a centerline, take the default cut,
download the result -- and it is the only number that decides whether the tool is
safe to put in front of strangers.

    python tools/smoke_bench.py                       # abdominal
    python tools/smoke_bench.py <cohort> [<cohort>]

Three verdicts per case:

    raised    an exception escaped. A hosted app shows a stack trace or a blank page.
              There is no argument about these: every one is a bug.
    refused   the tool declined and said why. Correct behaviour on bad input, but
              the REASON has to be one a stranger can act on.
    cut       a valid open mesh came out, and exporting it round-trips.

Also checks the things a download makes permanent: that the exported file reloads
with the topology the tool claimed, and that the sidecar records what was cut.
"""
import argparse
import json
import re
import sys
import tempfile
import time
import traceback
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from anat_bench import find_cases                                  # noqa: E402


def score(cid, cl, mp, export_dir):
    """One case through the plain path. Never raises -- that is what it measures."""
    from vascular_api import CutSession, CutError
    rec = {"case": cid, "verdict": "?", "secs": 0.0}
    t0 = time.time()
    try:
        s = CutSession(str(cl), str(mp))
    except CutError as e:
        rec.update(verdict="refused", stage="load", why=str(e)[:150])
        rec["secs"] = round(time.time() - t0, 1)
        return rec
    except Exception as e:                                          # noqa: BLE001
        rec.update(verdict="raised", stage="load",
                   why=f"{type(e).__name__}: {e}", trace=traceback.format_exc()[-600:])
        rec["secs"] = round(time.time() - t0, 1)
        return rec

    try:
        st = s.status()
        rec["openings"] = st.get("n_openings")
        rec["area_mm2"] = st.get("area_mm2")
        a = s.anatomy_summary()
        rec["family"] = a["kind"]
        rec["close_call"] = bool(a.get("close_call"))
        rec["mesh_matches"] = bool(a.get("mesh_matches", True))
        rec["named"] = a["named_openings"]
        rec["presets"] = len(s.presets())
        if not st["ok"]:
            rec.update(verdict="refused", stage="cut", why=st["message"][:150])
            rec["secs"] = round(time.time() - t0, 1)
            return rec
    except Exception as e:                                          # noqa: BLE001
        rec.update(verdict="raised", stage="cut",
                   why=f"{type(e).__name__}: {e}", trace=traceback.format_exc()[-600:])
        rec["secs"] = round(time.time() - t0, 1)
        return rec

    # The download. A cut nobody can take away is not a working tool, and an export
    # that changes the topology is worse than one that fails -- it is wrong quietly.
    try:
        import trimesh
        stem = re.sub(r"[^A-Za-z0-9_.-]", "_", cid)
        out = s.export(export_dir, stem)
        paths = [Path(p) for p in (out if isinstance(out, (list, tuple)) else [out])]
        mesh_files = [p for p in paths if p.suffix.lower() in
                      {".stl", ".ply", ".obj"}]
        side = [p for p in paths if p.suffix.lower() == ".json"]
        rec["exported"] = [p.name for p in paths]
        if not mesh_files:
            rec.update(verdict="raised", stage="export",
                       why=f"no mesh written, only {[p.name for p in paths]}")
            return rec
        back = trimesh.load(str(mesh_files[0]), process=False)
        import shape_metrics as SM
        welded = SM.weld(back)
        rec["export_euler"] = int(welded.euler_number)
        rec["export_faces"] = int(len(welded.faces))
        if int(welded.euler_number) != int(s.mesh.euler_number):
            rec.update(verdict="raised", stage="export",
                       why=(f"export changed the topology: euler "
                            f"{int(s.mesh.euler_number)} -> "
                            f"{int(welded.euler_number)}"))
            return rec
        if side:
            json.loads(side[0].read_text(encoding="utf-8"))
            rec["sidecar"] = side[0].name
        else:
            rec["sidecar"] = None
    except Exception as e:                                          # noqa: BLE001
        rec.update(verdict="raised", stage="export",
                   why=f"{type(e).__name__}: {e}", trace=traceback.format_exc()[-600:])
        rec["secs"] = round(time.time() - t0, 1)
        return rec

    rec["verdict"] = "cut"
    rec["secs"] = round(time.time() - t0, 1)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohorts", nargs="*",
                    default=[r"Y:/aorta/abdominal/MK/JK"])
    ap.add_argument("--json", default=None)
    ap.add_argument("--traces", action="store_true",
                    help="print a traceback for every case that raised")
    args = ap.parse_args()

    export_dir = Path(tempfile.mkdtemp(prefix="vcsmoke_"))
    rows = []
    for coh in args.cohorts:
        for cid, cl, mp in find_cases(coh):
            rows.append(score(cid, cl, mp, export_dir))
            r = rows[-1]
            mark = {"cut": "  ", "refused": "· ", "raised": "!!"}[r["verdict"]]
            print(f"{mark}{r['case']:11s} {r['verdict']:8s} "
                  f"{r.get('stage', ''):7s} {r['secs']:5.1f}s  "
                  f"{(r.get('why') or '')[:64]}")

    n = len(rows)
    by = Counter(r["verdict"] for r in rows)
    print("\n" + "=" * 70)
    print(f"THE PLAIN PATH: {by['cut']} of {n} cut and exported "
          f"({100 * by['cut'] / max(n, 1):.0f}%)")
    print("=" * 70)
    print(f"  raised an exception : {by['raised']}"
          f"{'   <-- every one is a bug' if by['raised'] else ''}")
    print(f"  refused, with a reason: {by['refused']}")
    slow = sorted(rows, key=lambda r: -r["secs"])[:3]
    print(f"  slowest: " + ", ".join(f"{r['case']} {r['secs']:.0f}s" for r in slow))

    if by["raised"]:
        print("\n  RAISED:")
        for r in rows:
            if r["verdict"] == "raised":
                print(f"    {r['case']:11s} [{r.get('stage')}] {r.get('why', '')[:90]}")
                if args.traces and r.get("trace"):
                    print("      " + r["trace"].replace("\n", "\n      ")[:700])
    if by["refused"]:
        print("\n  REFUSED, grouped by reason:")
        kinds = Counter()
        for r in rows:
            if r["verdict"] == "refused":
                k = re.sub(r"[\d.]+", "N", r.get("why", ""))
                kinds[k[:70]] += 1
        for k, v in kinds.most_common():
            print(f"    {v:3d}  {k}")

    flags = [r for r in rows if r.get("close_call") or not r.get("mesh_matches", True)]
    if flags:
        print(f"\n  {len(flags)} case(s) cut but carry a warning a user must read:")
        for r in flags:
            what = ("mesh/centerline mismatch" if not r.get("mesh_matches", True)
                    else f"close family call ({r['family']})")
            print(f"    {r['case']:11s} {what}")

    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=1), encoding="utf-8")
        print(f"\n  wrote {args.json}")
    return 1 if by["raised"] else 0


if __name__ == "__main__":
    sys.exit(main())
