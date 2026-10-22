"""Score the ANATOMY half of the tool over a whole cohort, and show what moved.

tools/cut_bench.py scores the cutter -- clean openings, foreign surface, refusals.
This scores everything upstream of it, which is where the abdominal complaints
actually live:

    family      abdominal or thoracic, against ground-truth labels if given
    route       the arc-length ruler every region and name is placed against
    naming      how many openings get an artery, and whether the renals do
    presets     the region tools that follow -- which is what a user notices

Usage:

    python tools/anat_bench.py --indir <cohort>              # score, save, compare
    python tools/anat_bench.py --indir <cohort> --labels labels.csv
    python tools/anat_bench.py --compare a.jsonl b.jsonl     # two saved runs
    python tools/anat_bench.py --indir <cohort> --case JK4_1 --verbose

A cohort directory is scanned for cases: any folder holding a centerline CSV and a
mesh, in the layout <patient>/<scan>/{centerline,mesh}/ or flat. Point --indir at
whatever you have; nothing here assumes a particular cohort.

The labels CSV needs a column of scan IDs and a column saying abdominal/thoracic
(or 1/2 trunk ends). Column names are sniffed, so the lab's own sheet works as-is.

Every run appends to tools/anat_bench.jsonl and prints the delta against the run
before it. That is the point: a change is worth keeping only if this moves.
"""
import argparse
import json
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from branch_tree import load_tree                                  # noqa: E402
import branch_cut as BC                                            # noqa: E402
from vascular_api import DEFAULT_SMOOTH_MM                         # noqa: E402

OUT = ROOT / "tools" / "anat_bench.jsonl"

# What "the region tools work" means, in one number. An abdominal aorta offers three
# renal-dependent presets, a visceral one and the whole vessel; four or more means a
# person can isolate a neck or an infrarenal segment, which is the actual complaint.
GOOD_PRESETS = 4


# ── finding cases ────────────────────────────────────────────────────────────
# Order matters and is not cosmetic. A case folder often holds BOTH a branched
# centerline and a single-line one, and the single line has no branches at all: every
# case then reads as two trunk ends with nothing leaving them, i.e. a thoracic aorta
# with no vessels. The first run of this benchmark scored 49 of 49 that way and
# reported zero names, zero presets and zero renals across the entire cohort. The
# branched file is what the tool is for, so it is looked for first.
CL_PAT = ("*_branches.csv", "*Branches*.csv",
          "*CentrelinePoints*.csv", "*CenterlinePoints*.csv", "*centerline*.csv",
          "*.vtp", "*_cl.csv")
MESH_PAT = ("*.mat", "*.stl", "*.ply", "*.obj")


def _first(folder, patterns):
    for p in patterns:
        hits = sorted(folder.glob(p))
        if hits:
            # An explicitly remeshed copy is not the source geometry; prefer the plain
            # one so the benchmark scores what a user would open.
            plain = [h for h in hits if not re.search(r"_o\d+K", h.name)]
            return (plain or hits)[0]
    return None


# Where the BRANCHED data lives, most specific first. Cohorts are laid out
# differently -- one keeps the branched centerline beside the plain one in
# centerline/, another separates them into branch_centerline/ and centerline/ -- and
# picking the plain one means every case reads as a vessel with no branches at all.
# That is not a subtle failure: it scored a 29-case thoracic cohort at zero names,
# zero presets and zero renals, on geometry that names its arch vessels perfectly.
CL_DIRS = ("branch_centerline", "branches", "centerline", "")
# branch_mesh and branch_model are the same idea under two names, and the thoracic
# cohort uses both -- KY160_1 has branch_mesh/, KK40_3 has branch_model/, and the
# rest have an EMPTY branch_mesh/ beside a populated branch_model/. Falling through
# to mesh/ pairs a branched centerline with a mesh of a different extent, and the
# anatomy call then goes wrong silently: KK40_3, KR47_1 and KR52_1 all read as
# abdominal that way, and read correctly the moment the right mesh is used.
MESH_DIRS = ("branch_mesh", "branch_model", "mesh", "")

# Never descend into these. Source imaging is not this tool's business, it is large,
# and it is the one thing that must not be opened by a sweep that walks a folder.
SKIP_DIRS = {"dicom", "dicom_sip", "sip", "mhd", "parquet", "raw", "__pycache__"}


def find_cases(indir):
    """Every case under *indir*: (id, centerline, mesh). Layout-agnostic."""
    indir = Path(indir)
    if not indir.exists():
        raise SystemExit(f"{indir} is not reachable")
    out, seen = [], set()
    folders = [indir]
    for p in sorted(indir.rglob("*")):
        if not p.is_dir():
            continue
        if any(part.lower() in SKIP_DIRS for part in p.relative_to(indir).parts):
            continue
        folders.append(p)
    for folder in folders:
        cl = next((f for d in CL_DIRS
                   for f in [_first(folder / d if d else folder, CL_PAT)] if f), None)
        mesh = next((f for d in MESH_DIRS
                     for f in [_first(folder / d if d else folder, MESH_PAT)] if f),
                    None)
        if cl is None or mesh is None:
            continue
        cid = folder.name if folder != indir else indir.name
        if cid in seen:
            continue
        seen.add(cid)
        out.append((cid, cl, mesh))
    return sorted(out)


def read_labels(path):
    """{scan id: 1 or 2 trunk ends}, sniffing the columns rather than demanding them."""
    if not path:
        return {}
    import pandas as pd
    df = pd.read_csv(path) if str(path).lower().endswith(".csv") else \
        pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]
    idc = next((c for c in df.columns
                if re.search(r"scan|case|id\b|name", c, re.I)), df.columns[0])
    labc = next((c for c in df.columns
                 if re.search(r"famil|anatom|abdom|thorac|type|trunk|label", c, re.I)),
                None)
    if labc is None:
        raise SystemExit(f"{path}: no column looks like the abdominal/thoracic label. "
                         f"Columns are {list(df.columns)}")
    out = {}
    for _, r in df.iterrows():
        sid = str(r[idc]).strip()
        raw = str(r[labc]).strip().lower()
        if not sid or sid == "nan" or not raw or raw == "nan":
            continue
        if raw.startswith("abd") or raw in ("1", "1.0"):
            out[sid] = 1
        elif raw.startswith("thor") or raw in ("2", "2.0"):
            out[sid] = 2
    print(f"  labels: {len(out)} from {Path(path).name} "
          f"(id column {idc!r}, label column {labc!r})")
    return out


# ── scoring one case ─────────────────────────────────────────────────────────
def score_case(cid, cl, mesh_path, want=None, verbose=False,
               cut_presets=False):
    rec = {"case": cid, "centerline": str(cl), "mesh": str(mesh_path)}
    t0 = time.time()
    # Smooth the way the TOOL does, or this scores a centerline nobody ever cuts.
    # CutSession probes the line as supplied and, when its turning exceeds the jagged
    # threshold, reloads it smoothed -- 6 of 30 thoracic cases and several abdominals.
    # Loading raw here made JK30_4 offer two region presets through the session and
    # none through the benchmark, on the same files.
    probe = load_tree(str(cl), smooth_mm=0.0)
    rec["jagged"] = bool(probe.jagged)
    tree = (load_tree(str(cl), smooth_mm=DEFAULT_SMOOTH_MM) if probe.jagged
            else probe)
    mesh = BC.load_mat_mesh(str(mesh_path)) if str(mesh_path).endswith(".mat") else None
    if mesh is None:
        import trimesh
        mesh = trimesh.load(str(mesh_path), process=False)

    BC.classify_terminals(tree, mesh, verbose=verbose)
    n_trunk = sum(1 for t in tree.terminals if getattr(t, "is_trunk", False))
    rec["family"] = n_trunk
    rec["family_want"] = want
    rec["family_ok"] = (None if want is None else int(n_trunk == want))
    rec["outflow_split"] = (None if getattr(tree, "outflow_split", None) is None
                            else round(float(tree.outflow_split), 3))
    # A call close to the cut-off is the one worth a human eye, so count them.
    rec["family_close"] = (None if rec["outflow_split"] is None else
                           int(abs(rec["outflow_split"] - BC.PAIR_AREA_FRACTION) < 0.08))

    names, vessels, notes = BC.name_anatomy(tree, mesh)
    jg = BC._junction_graph(tree)
    route = BC.aorta_route(tree, jg)
    rec["n_terminals"] = len(tree.terminals)
    rec["route_mm"] = None if route is None else round(float(route.length), 1)
    rec["n_segments"] = 0 if route is None else len(route.chain)
    # The ruler must actually span the vessel. A route far shorter than the model is
    # the failure mode that hides every landmark below it.
    tips = np.array([t.tip for t in tree.terminals], float)
    rec["model_span_mm"] = round(float(np.linalg.norm(tips.max(0) - tips.min(0))), 1)
    rec["route_frac"] = (None if not rec["route_mm"] or not rec["model_span_mm"]
                         else round(rec["route_mm"] / rec["model_span_mm"], 2))

    rec["n_named"] = len(names)
    rec["names"] = sorted(names.values())
    rec["renals"] = sorted(v for v in names.values() if "renal" in v.lower())
    rec["n_renals"] = len(rec["renals"])
    rec["iliacs"] = sorted(v for v in list(names.values()) + list(vessels.values())
                           if "iliac" in v.lower())
    rec["notes"] = list(notes)
    rec["ambiguous_renals"] = int(any("could be the renals" in n for n in notes))

    try:
        presets = BC.region_presets(tree, mesh, names, vessels)
    except Exception as e:                                          # noqa: BLE001
        presets = []
        rec["preset_error"] = f"{type(e).__name__}: {e}"
    rec["n_presets"] = len(presets)
    rec["presets"] = [p["name"] for p in presets]
    rec["presets_ok"] = int(len(presets) >= GOOD_PRESETS)

    # Offering a region is not the same as being able to cut one. Measured over 64
    # region cuts on 36 abdominal scans, 4 leave an EXTRA opening -- a boundary plane
    # reaching sideways through a vessel it was not aimed at, which on an abdominal
    # aorta is a renal running parallel rather than leaving laterally. That is the
    # lab's "too big cutting plane as the renal goes straight down", and it is
    # enriched five-fold in the scans they reported: 2 of 10 against 2 of 52.
    #
    # The tool refuses those cuts rather than exporting them, which is right. What it
    # does not do is offer a way forward, the way a branch cut offers pull-back --
    # and a smaller plane demonstrably fixes both cases. Counted here so that
    # remains visible.
    if cut_presets and presets:
        rec["cut"] = _score_presets(tree, mesh, names, vessels, presets)
        rec["presets_cut_ok"] = sum(1 for c in rec["cut"] if c.get("ok"))
        rec["presets_cut_extra"] = sum(1 for c in rec["cut"]
                                       if c.get("extra_openings"))
    rec["secs"] = round(time.time() - t0, 1)
    return rec


def _score_presets(tree, mesh, names, vessels, presets):
    """Apply each region preset and say whether its two boundary planes cut cleanly.

    Two planes through a closed surface leave exactly two openings. More than that
    means one of them severed a bystander, and counting REMOVED AREA cannot see it:
    these presets begin at the lowest renal, so taking renal surface above that
    boundary is the region doing its job, not clipping it.

    THIS IS A PROXY, and it flatters. It cuts only the two boundary planes; a real
    region cut also cuts every branch inside the region and runs validate_cut over
    the lot. On the abdominal cohort the proxy says 177 of 182 and the tool actually
    accepts 146 -- a gap of 31 cuts that steering by this number alone would never
    show. tools/region_bench.py measures the real thing; use this to find out WHY a
    region fails, not WHETHER it does.
    """
    jg = BC._junction_graph(tree)
    route = BC.aorta_route(tree, jg)
    if route is None:
        return []
    snap = BC.snap_tolerance(mesh)
    out = []
    for p in presets:
        s1 = route.length if p["s1"] is None else p["s1"]
        rec = {"preset": p["name"]}
        try:
            bounds, _inside, _rm, notes = BC.region_plan(
                tree, p["s0"], s1, route, mesh=mesh, clear_ostia=True)
            if not bounds:
                rec["err"] = "no bounds to place"
                out.append(rec)
                continue
            scale = BC.fit_region_reach(mesh, tree, bounds, snap_mm=snap)
            rec["reach_scale"] = scale
            m = mesh.copy()
            for t in bounds:
                m = BC.cut_terminal(m, tree, t, getattr(t, "_region_offset", 6.0),
                                    snap_mm=snap)
            # N planes through a closed surface leave N openings -- not always two.
            # A region whose boundary IS the vessel's own end places one plane, and
            # hard-coding 2 scored 49 perfectly good cuts as failures.
            loops = len(BC.boundary_loops(m))
            want = len(bounds)
            rec.update(loops=loops, extra_openings=max(0, loops - want),
                       euler=int(m.euler_number), ok=(loops == want),
                       n_bounds=want, notes=list(notes)[:2])
        except Exception as e:                                      # noqa: BLE001
            rec["err"] = f"{type(e).__name__}: {str(e)[:60]}"
        out.append(rec)
    return out


# ── the run ──────────────────────────────────────────────────────────────────
def summarise(rows):
    ok = [r for r in rows if "error" not in r]
    ab = [r for r in ok if r["family"] == 1]
    labelled = [r for r in ok if r.get("family_ok") is not None]
    return {
        "cases": len(rows),
        "failed": len(rows) - len(ok),
        "family_correct": sum(r["family_ok"] for r in labelled) if labelled else None,
        "family_labelled": len(labelled),
        "family_close_calls": sum(r.get("family_close") or 0 for r in ok),
        "abdominal": len(ab),
        "no_route": sum(1 for r in ok if r["route_mm"] is None),
        "short_route": sum(1 for r in ok
                           if r.get("route_frac") is not None and r["route_frac"] < 0.5),
        "renals_named": sum(1 for r in ab if r["n_renals"] >= 2),
        "renals_ambiguous": sum(r.get("ambiguous_renals") or 0 for r in ab),
        "presets_ok": sum(r["presets_ok"] for r in ok),
        "presets_total": sum(r["n_presets"] for r in ok),
        "named_total": sum(r["n_named"] for r in ok),
        "presets_cut_ok": (sum(r.get("presets_cut_ok", 0) for r in ok)
                           if any("cut" in r for r in ok) else None),
        "presets_cut_extra": (sum(r.get("presets_cut_extra", 0) for r in ok)
                              if any("cut" in r for r in ok) else None),
    }


KEYS = [("cases", "cases scored", 0), ("failed", "raised an exception", -1),
        ("family_correct", "family called correctly", 1),
        ("family_close_calls", "close family calls", 0),
        ("no_route", "no aortic route at all", -1),
        ("short_route", "route under half the model", -1),
        ("renals_named", "abdominals with renals named", 1),
        ("renals_ambiguous", "...of which were a choice", 0),
        ("presets_ok", f"cases with {GOOD_PRESETS}+ region presets", 1),
        ("presets_total", "region presets offered", 1),
        ("named_total", "openings named", 1),
        ("presets_cut_ok", "region BOUNDS that cut cleanly", 1),
        ("presets_cut_extra", "region bounds severing a bystander", -1)]


def show(summary, prev=None):
    print("\n" + "=" * 68)
    print("ANATOMY BENCHMARK")
    print("=" * 68)
    for key, label, better in KEYS:
        now = summary.get(key)
        if now is None:
            continue
        line = f"  {label:38s} {now:6}"
        if prev is not None and prev.get(key) is not None:
            d = now - prev[key]
            if d:
                mark = ("BETTER" if (d > 0) == (better > 0) and better
                        else "WORSE" if better else "")
                line += f"   {d:+d}  {mark}"
            else:
                line += "    ="
        if key == "family_correct" and summary.get("family_labelled"):
            line += f"   of {summary['family_labelled']} labelled"
        print(line)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=None, help="cohort folder to scan")
    ap.add_argument("--labels", default=None,
                    help="CSV/XLSX of ground-truth abdominal/thoracic labels")
    # Cohort membership IS the label. An abdominal cohort is abdominal; a thoracic
    # one is thoracic. That is ground truth for the family call over every case,
    # without anyone maintaining a file, and it is how this benchmark scores the one
    # decision that gates everything downstream.
    ap.add_argument("--family-label", type=int, default=None, choices=(1, 2),
                    help="every case here has this many trunk ends "
                         "(1 abdominal, 2 thoracic). Overridden per case by --labels.")
    ap.add_argument("--except-case", action="append", default=None,
                    metavar="ID=N",
                    help="a known exception to --family-label, e.g. KP35_1=2")
    ap.add_argument("--case", action="append", default=None,
                    help="score only these case ids (repeatable)")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--compare", nargs=2, default=None,
                    help="two saved runs, show the delta and exit")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--cut", action="store_true",
                    help="also CUT each region preset and check it leaves "
                         "exactly two openings (slow, ~2s per preset)")
    ap.add_argument("--worst", type=int, default=12,
                    help="how many problem cases to list")
    args = ap.parse_args()

    if args.compare:
        a, b = (json.loads(Path(p).read_text().splitlines()[-1]) for p in args.compare)
        show(b["summary"], a["summary"])
        return 0
    if not args.indir:
        ap.error("--indir is required unless --compare is used")

    labels = read_labels(args.labels)
    if args.family_label:
        for cid, _cl, _m in find_cases(args.indir):
            labels.setdefault(cid, int(args.family_label))
        for spec in (args.except_case or []):
            cid, _, n = spec.partition("=")
            labels[cid.strip()] = int(n)
        print(f"  family label: {args.family_label} for this cohort"
              + (f", except {args.except_case}" if args.except_case else ""))
    cases = find_cases(args.indir)
    if args.case:
        want = set(args.case)
        cases = [c for c in cases if c[0] in want]
    print(f"  {len(cases)} case(s) under {args.indir}")

    rows = []
    for cid, cl, mesh in cases:
        try:
            rows.append(score_case(cid, cl, mesh, labels.get(cid), args.verbose,
                                   cut_presets=args.cut))
        except Exception as e:                                      # noqa: BLE001
            print(f"  {cid:12s} RAISED {type(e).__name__}: {str(e)[:60]}")
            if args.verbose:
                traceback.print_exc()
            rows.append({"case": cid, "error": f"{type(e).__name__}: {e}"})

    print(f"\n{'case':12s} {'fam':>4} {'route':>7} {'frac':>5} {'namd':>5} "
          f"{'ren':>4} {'pre':>4}  flags")
    for r in rows:
        if "error" in r:
            print(f"{r['case']:12s}  --  RAISED {r['error'][:44]}")
            continue
        flags = []
        if r.get("family_ok") == 0:
            flags.append(f"FAMILY: called {r['family']}, labelled {r['family_want']}")
        if r.get("family_close"):
            flags.append("close family call")
        if r["route_mm"] is None:
            flags.append("NO ROUTE")
        elif r.get("route_frac") is not None and r["route_frac"] < 0.5:
            flags.append("short route")
        if r["family"] == 1 and r["n_renals"] < 2:
            flags.append("no renals")
        if r.get("ambiguous_renals"):
            flags.append("renals were a choice")
        if not r["presets_ok"]:
            flags.append(f"only {r['n_presets']} preset(s)")
        route = "   none" if r["route_mm"] is None else f"{r['route_mm']:7.0f}"
        frac = "    -" if r.get("route_frac") is None else f"{r['route_frac']:5.2f}"
        print(f"{r['case']:12s} {r['family']:4d} {route} {frac} "
              f"{r['n_named']:5d} {r['n_renals']:4d} {r['n_presets']:4d}  "
              f"{'; '.join(flags)}")

    summary = summarise(rows)
    prev = None
    out = Path(args.out)
    if out.exists():
        lines = [ln for ln in out.read_text().splitlines() if ln.strip()]
        if lines:
            prev = json.loads(lines[-1]).get("summary")
    show(summary, prev)

    if args.cut:
        # Said every time, because this is the trap the loop keeps walking into: a
        # change can move the BOUNDS count and move the number a user feels by
        # nothing at all. Extending the plane-reach ladder to 0.2x did exactly that,
        # 169 -> 170 bounds and 146 -> 146 end to end.
        print("\n  NB: 'region BOUNDS' counts two planes on a pristine mesh. For "
              "what\n      the tool actually accepts, run:  "
              f"python tools/region_bench.py {args.indir}")

    rec = {"when": time.strftime("%Y-%m-%d %H:%M:%S"), "indir": str(args.indir),
           "summary": summary, "cases": rows}
    with out.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"\n  appended to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
