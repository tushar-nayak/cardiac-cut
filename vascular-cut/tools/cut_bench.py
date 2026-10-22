"""Score the cutter against the lab cohorts, and against its own past self.

This is the benchmark the cutting work is steered by. It exists because two sample
geometries are not enough to settle a geometric rule: every threshold tried during the
tearing investigation won on one of them and lost on the other, and the trunk-versus-
branch split only became visible at ~60 cases.

    python tools/cut_bench.py run                 # both cutters, every case
    python tools/cut_bench.py run --limit 8       # a quick pass while iterating
    python tools/cut_bench.py score               # read the last run

WHAT IT MEASURES, and why only these three things

    flat      the opening's rim sits on the plane that made it (<= 0.5 mm).
              A torn rim is the failure the whole investigation started from.
    euler     the result has exactly one more opening than it started with.
              Catches a cut that punched a second hole somewhere else.
    foreign   surface area removed from a vessel this cut was not cutting.
              Nothing else detects this. A cut can leak into the parent, still
              leave one connected disk, still leave one boundary loop, and keep
              the euler number exactly right. It was measured at 206,162 mm2
              across the cohort before this was tracked, and nobody could see it.

CLEAN means all three at once. That is the number to move.

Judged on OUTPUT ONLY, never on which faces either implementation selected, so the
comparison stays honest across a rewrite.

WHERE THE CASES COME FROM

Cohort roots are not committed. Set VC_COHORTS to a semicolon-separated list of
directories, or drop a cohorts.local.json beside this file holding a list of paths.
Both cohorts are laid out differently and neither is guessable, so the mesh and
centerline are found by trying known names in order -- see MESH_PREF / CL_PREF. Get
that wrong and the failure is silent and total: pairing a branch centerline with a
trunk-only mesh puts every branch outside the surface, and 96 of 150 segments came
back with a negative radius before it was caught.
"""
import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import traceback
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "bench_last.jsonl"

# Depths as a fraction of vessel length. The usable range is where a cut has to work;
# 85-90% is where every failure lives, and is scored separately so a fix that helps
# one cannot hide behind the other.
FRACS = (0.20, 0.50, 0.70, 0.80, 0.85, 0.90)
USABLE_MAX = 0.80
FLAT_MM = 0.5
FOREIGN_OK_MM2 = 1.0          # ostium fuzz, not damage

# The branched mesh, under the names the two cohorts happen to use. branch_model/
# holds it in 25 of 27 GT cases and branch_mesh/ in the other two; mesh/<id>.mat is
# the trunk-only surface in GT and the branched one in JK, so it comes last.
MESH_PREF = ("branch_model/{id}.mat", "branch_mesh/{id}.mat",
             "mesh/{id}_branches.mat", "mesh/{id}.mat")
CL_PREF = ("branch_centerline/{id}_CentrelinePoints_branches.csv",
           "centerline/{id}_CentrelinePoints_branches.csv",
           "branch_centerline/{id}_CentrelinePoints.csv",
           "centerline/{id}_CentrelinePoints.csv")


def cohort_roots():
    """Where the cases live. Loud on every way this can be wrong.

    A silently empty cohort scores zero cuts and reads as "nothing regressed", so
    a bad config has to stop the run rather than shrink it. Write the paths with
    forward slashes: a lone backslash is an invalid JSON escape.
    """
    env = os.environ.get("VC_COHORTS", "").strip()
    if env:
        roots = [Path(p.strip()) for p in env.split(";") if p.strip()]
    else:
        cfg = HERE / "cohorts.local.json"
        if not cfg.exists():
            raise SystemExit(f"no cohort roots: set VC_COHORTS or write {cfg}")
        try:
            roots = [Path(p) for p in json.loads(cfg.read_text(encoding="utf-8"))]
        except json.JSONDecodeError as e:
            raise SystemExit(
                f"{cfg} is not valid JSON ({e}). Use forward slashes in the paths: "
                'a single backslash is an invalid escape.') from None
    if not roots:
        raise SystemExit("no cohort roots configured")
    return roots


def _first(root, patterns):
    for pat in patterns:
        p = root / pat.format(id=root.name)
        if p.exists():
            return p
    return None


def find_cases():
    """(name, csv, mat) for every case with both files.

    Handles a cohort laid out as <root>/<case>/ and as <root>/<case>/<sub>/, because
    GT is the first and JK is the second.
    """
    out = []
    for root in cohort_roots():
        if not root.exists():
            print(f"  cohort root missing: {root}", file=sys.stderr)
            continue
        for case in sorted(p for p in root.iterdir() if p.is_dir()):
            mat, csv = _first(case, MESH_PREF), _first(case, CL_PREF)
            if mat and csv:
                out.append((f"{root.name}/{case.name}", csv, mat))
                continue
            for sub in sorted(p for p in case.iterdir() if p.is_dir()):
                mat, csv = _first(sub, MESH_PREF), _first(sub, CL_PREF)
                if mat and csv:
                    out.append((f"{root.name}/{sub.name}", csv, mat))
    return out


def load_rev(alias, rev):
    """branch_cut.py at a git revision, importable alongside the working copy."""
    src = subprocess.run(["git", "show", f"{rev}:branch_cut.py"], cwd=ROOT,
                         capture_output=True, check=True).stdout
    tmp = Path(tempfile.mkdtemp()) / f"{alias}.py"
    tmp.write_bytes(src)                       # bytes: the file is not pure ASCII
    spec = importlib.util.spec_from_file_location(alias, str(tmp))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[alias] = mod
    spec.loader.exec_module(mod)
    return mod


def build(BC, csv, mat, smooth_mm=0.0, n_trunk_ends=None):
    from branch_tree import load_tree
    mesh = BC.load_mat_mesh(str(mat))
    tree = load_tree(str(csv), smooth_mm=smooth_mm)
    BC.classify_terminals(tree, mesh, n_trunk_ends=n_trunk_ends, verbose=False)
    BC.name_anatomy(tree, mesh)
    # EVERY segment, not just the terminals'. An interior segment otherwise keeps a
    # radius of None, which the mismatch check below reads as zero.
    for seg in tree.segments.values():
        if getattr(seg, "_radius_mm", None) is None:
            seg._radius_mm = BC.segment_radius(mesh, seg)
    return mesh, tree


def mismatched(tree):
    """How many segments have a centerline lying outside the mesh.

    A negative segment_radius means that stretch of centerline is outside the surface.
    A few is a local defect; half means the two files are not the same anatomy and
    nothing measured on them means anything.
    """
    rs = [float(getattr(s, "_radius_mm", 0.0) or 0.0) for s in tree.segments.values()]
    return sum(1 for r in rs if r <= 0.0), len(rs)


def survived(mesh, out):
    """Which original faces came through intact, by exact centroid.

    A face the cut did not touch keeps its three vertices to the bit. Implementation
    independent on purpose: it measures what a cutter DID, not how it decided.
    """
    a = np.round(np.asarray(mesh.triangles_center, float), 6)
    b = set(map(tuple, np.round(np.asarray(out.triangles_center, float), 6).tolist()))
    return np.fromiter((tuple(r) in b for r in a.tolist()), dtype=bool, count=len(a))


def judge(BC, mesh, tree, t, off, owner):
    point, normal = t.cut_plane(off)
    out = BC.cut_terminal(mesh, tree, t, off, snap_mm=BC.snap_tolerance(mesh))
    loops = BC.boundary_loops(out)
    verts = np.asarray(out.vertices)
    dev = None
    if loops:
        best, bd = None, np.inf
        for loop in loops:
            c = verts[np.asarray(loop)].mean(axis=0)
            d = float(np.linalg.norm(c - point))
            if d < bd:
                best, bd = loop, d
        dev = float(np.abs((verts[np.asarray(best)] - point) @ normal).max())
    gone = ~survived(mesh, out)
    foreign = gone & (owner != t.segment.name)
    sub = getattr(t, "subtree_segments", None)
    if sub:
        foreign = foreign & ~np.isin(owner, list(sub))
    return {"dev": dev, "euler": int(out.euler_number), "loops": len(loops),
            "foreign_mm2": float(mesh.area_faces[foreign].sum())}


def run_case(name, csv, mat, variants):
    built = {k: build(BC, csv, mat) for k, BC in variants.items()}
    ref = list(variants)[-1]
    mesh_r, tree_r = built[ref]
    nbad, ntot = mismatched(tree_r)
    if ntot and nbad > 0.3 * ntot:
        raise ValueError(f"mesh and centerline do not match: {nbad} of {ntot} "
                         f"segments have a centerline outside the surface")
    owner, _ = variants[ref].vessel_coords(mesh_r, tree_r)

    # Matched by LABEL, never by position. classify_terminals REORDERS the list it
    # returns -- `list(trunk) + branches`, branches sorted by height -- so the moment
    # two variants disagree about how many trunk ends a vessel has, terminals[i] is a
    # different artery in each and the bench silently compares one vessel's cut
    # against another's, at an offset taken from a third's length. Every cutter change
    # up to now classified identically, so index matching happened to hold; the first
    # change that touches classify_terminals turns the whole comparison into noise
    # (KK40_3 headed the "got worse" list on exactly this).
    by_label = {k: {t.label: t for t in tree.terminals}
                for k, (_m, tree) in built.items()}
    rows = []
    for t_ref in tree_r.terminals:
        L = t_ref.segment.length
        for f in FRACS:
            off = max(f * L, 1.0)
            rec = {"terminal": t_ref.name, "seg": t_ref.segment.name,
                   "label": t_ref.label,
                   "L": round(L, 1), "frac": f, "off": round(off, 1),
                   "r": round(float(t_ref.segment._radius_mm), 1),
                   "at_end": bool(t_ref.at_end), "mode": t_ref.mode,
                   "trunk": bool(getattr(t_ref, "is_trunk", False))}
            for k, BC in variants.items():
                mesh, tree = built[k]
                t_k = by_label[k].get(t_ref.label)
                if t_k is None:
                    # The variants do not even agree this opening exists. Say so
                    # rather than fall back to a positional guess.
                    rec[k] = {"err": f"no terminal {t_ref.label} in this variant"}
                    continue
                try:
                    rec[k] = judge(BC, mesh, tree, t_k, off, owner)
                except Exception as e:                      # noqa: BLE001
                    rec[k] = {"err": f"{type(e).__name__}: {str(e)[:140]}"}
                # A row where the variants disagree on trunk-vs-branch is a real
                # comparison of the same artery, but of two different DECISIONS about
                # it -- different name, different cut mode. Flagged so scoring can
                # separate "this rule cuts better" from "this rule cuts elsewhere".
                rec[k]["is_trunk"] = bool(getattr(t_k, "is_trunk", False))
                rec[k]["mode"] = t_k.mode
            rows.append(rec)
    reclassified = sum(1 for r in rows
                       if len({r[k].get("is_trunk") for k in variants
                               if isinstance(r.get(k), dict)
                               and "is_trunk" in r[k]}) > 1)
    return {"case": name, "faces": int(len(mesh_r.faces)),
            "terminals": len(tree_r.terminals), "rows": rows,
            "reclassified_rows": reclassified}


def cmd_run(args):
    import branch_cut as BC_new
    variants = {"base": load_rev("bc_base", args.baseline), "new": BC_new}
    cases = find_cases()
    if not cases:
        raise SystemExit("no cases found; set VC_COHORTS or tools/cohorts.local.json")
    cases = cases[:args.limit]
    out_path = Path(args.out)

    seen = set()
    if out_path.exists() and args.resume:
        seen = {json.loads(ln)["case"] for ln in
                out_path.read_text(encoding="utf-8").splitlines() if ln.strip()}
    elif out_path.exists():
        # Keep the run being replaced. A full sweep is 35 minutes and the record of
        # what a rule did is the only reason to have run it; a --limit 12 pass
        # otherwise silently overwrites the 60-case reference it is being judged
        # against, and the numbers are gone until someone spends the 35 minutes again.
        prev = out_path.with_suffix(".prev.jsonl")
        prev.unlink(missing_ok=True)
        out_path.rename(prev)
        print(f"previous run kept as {prev.name}", flush=True)
    todo = [c for c in cases if c[0] not in seen]
    print(f"{len(cases)} cases, {len(todo)} to run, baseline {args.baseline}",
          flush=True)

    with out_path.open("a", encoding="utf-8") as fh:
        for n, (name, csv, mat) in enumerate(todo, 1):
            t0 = time.time()
            try:
                rec = run_case(name, csv, mat, variants)
            except Exception as e:                          # noqa: BLE001
                rec = {"case": name, "error": f"{type(e).__name__}: {e}",
                       "trace": traceback.format_exc()[-300:], "rows": []}
            rec["seconds"] = round(time.time() - t0, 1)
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print(f"  [{n}/{len(todo)}] {name}: "
                  f"{rec.get('error', str(len(rec['rows'])) + ' cuts')} "
                  f"{rec['seconds']:.0f}s", flush=True)
    return cmd_score(args)


def _bucket(row):
    return ("trunk" if row.get("trunk") else "branch",
            "usable" if row["frac"] <= USABLE_MAX + 1e-9 else "junction")


def cmd_score(args):
    path = Path(args.out)
    if not path.exists():
        raise SystemExit(f"no results at {path}; run first")
    recs = [json.loads(ln) for ln in
            path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    ok = [r for r in recs if not r.get("error")]
    rejected = [r for r in recs if r.get("error")]

    agg = defaultdict(lambda: defaultdict(float))
    skipped = 0
    reclassified = 0
    variants = ("base", "new")
    for r in ok:
        for row in r["rows"]:
            # A cut neither cutter can attempt says nothing about either. "no distal
            # region" means the centerline runs past the end of the mesh.
            if any("no distal region" in row.get(v, {}).get("err", "")
                   for v in variants):
                skipped += 1
                continue
            # Nor does a cut the two variants make on DIFFERENT TERMS. If they
            # disagree about whether this artery is a trunk end, they give it a
            # different cut mode and a different plane, and the row measures that
            # disagreement rather than the cutting. Counted, reported, and kept out
            # of the table -- otherwise a change to classify_terminals reads as a
            # change to the cutter (a family fix showed CLEAN 2056 -> 2029 this way,
            # entirely from 14 cases being cut as a different anatomy).
            trunkness = {row[v].get("is_trunk") for v in variants
                         if isinstance(row.get(v), dict) and "is_trunk" in row[v]}
            if len(trunkness) > 1:
                reclassified += 1
                continue
            key = _bucket(row)
            for v in variants:
                d = row.get(v, {})
                a = agg[key + (v,)]
                a["n"] += 1
                if "err" in d:
                    a["refused"] += 1
                    continue
                flat = d["dev"] is not None and d["dev"] <= FLAT_MM
                a["flat"] += flat
                a["euler"] += d["euler"] == 1
                a["clean"] += (flat and d["euler"] == 1
                               and d["foreign_mm2"] <= FOREIGN_OK_MM2)
                a["foreign"] += d["foreign_mm2"]

    print(f"\ncases {len(recs)}: {len(ok)} scored, {len(rejected)} rejected"
          f"   cuts skipped {skipped}"
          + (f", reclassified {reclassified}" if reclassified else ""))
    if reclassified:
        print(f"    {reclassified} cut(s) held out: the variants disagree on "
              f"trunk-vs-branch there, so they are not cutting on the same terms. "
              f"The table below compares the CUTTER; the family call is judged "
              f"separately (see tools/ANATOMY.md).")
    for r in rejected[:5]:
        print(f"    {r['case']}: {r['error'][:88]}")
    print()
    hdr = f"{'':<22}{'cuts':>6}{'flat':>12}{'euler':>12}{'CLEAN':>12}" \
          f"{'refused':>9}{'foreign mm2':>13}"
    print(hdr)
    print("-" * len(hdr))
    totals = defaultdict(lambda: defaultdict(float))
    for kind in ("branch", "trunk"):
        for band in ("usable", "junction"):
            print(f"{kind} {band}")
            for v in variants:
                a = agg[(kind, band, v)]
                n = a["n"] or 1
                for f in ("n", "flat", "euler", "clean", "refused", "foreign"):
                    totals[v][f] += a[f]
                print(f"  {v:<20}{int(a['n']):6d}"
                      f"{int(a['flat']):6d}{100 * a['flat'] / n:5.0f}%"
                      f"{int(a['euler']):6d}{100 * a['euler'] / n:5.0f}%"
                      f"{int(a['clean']):6d}{100 * a['clean'] / n:5.0f}%"
                      f"{int(a['refused']):9d}{a['foreign']:13.0f}")
    print("-" * len(hdr))
    print("ALL")
    for v in variants:
        a = totals[v]
        n = a["n"] or 1
        print(f"  {v:<20}{int(a['n']):6d}"
              f"{int(a['flat']):6d}{100 * a['flat'] / n:5.0f}%"
              f"{int(a['euler']):6d}{100 * a['euler'] / n:5.0f}%"
              f"{int(a['clean']):6d}{100 * a['clean'] / n:5.0f}%"
              f"{int(a['refused']):9d}{a['foreign']:13.0f}")

    # Every cut the change makes worse, named. A summary that only moves totals hides
    # a rule that fixes forty cases and breaks five in a way nobody would accept.
    worse, worse_reclass = [], []
    for r in ok:
        for row in r["rows"]:
            b, n = row.get("base", {}), row.get("new", {})
            if "no distal region" in b.get("err", ""):
                continue
            ok_b = ("err" not in b and b["dev"] is not None and b["dev"] <= FLAT_MM
                    and b["euler"] == 1 and b["foreign_mm2"] <= FOREIGN_OK_MM2)
            ok_n = ("err" not in n and n["dev"] is not None and n["dev"] <= FLAT_MM
                    and n["euler"] == 1 and n["foreign_mm2"] <= FOREIGN_OK_MM2)
            if ok_b and not ok_n:
                # Same split as the table: a cut the two variants made on different
                # terms is not the cutter getting worse.
                tk = {row[v].get("is_trunk") for v in variants
                      if isinstance(row.get(v), dict) and "is_trunk" in row[v]}
                (worse_reclass if len(tk) > 1 else worse).append(
                    (r["case"], row, b, n))

    if worse_reclass:
        print(f"\n{len(worse_reclass)} cut(s) clean before and not now ONLY because "
              f"the artery is classified differently -- a different cut mode on a "
              f"different plane. Judge those on the family call, not here.")
    print(f"\nclean under the baseline, not clean now: {len(worse)}")
    for case, row, b, n in worse[:args.show]:
        # .get's default is built EAGERLY, so a refused cut -- which carries "err" and
        # no "dev" -- used to raise KeyError here and take the whole report with it,
        # after the 35-minute run had already finished.
        why = n.get("err") or (f"dev {n['dev']:.2f} euler {n['euler']} "
                               f"foreign {n['foreign_mm2']:.0f}")
        print(f"  {case:<16} {row['terminal']:<22} {100 * row['frac']:3.0f}%  "
              f"r={row['r']:5.1f}  {'trunk' if row.get('trunk') else 'branch'}  "
              f"-> {why[:70]}")
    return 0


def cmd_cases(args):
    """Show what each case resolved to, so a wrong pairing is visible before a run.

    Pairing a branch centerline with a trunk-only mesh is silent at load time and
    ruins every number downstream, so the resolution gets printed rather than
    trusted.
    """
    cases = find_cases()
    if not cases:
        raise SystemExit("no cases found under the configured roots")
    def shape(p, case):
        return f"{p.parent.name}/{p.name.replace(case, '<id>')}"

    tally = {}
    for name, csv, mat in cases:
        case = name.split("/", 1)[1]
        key = (shape(mat, case), shape(csv, case))
        tally.setdefault(key, []).append(name)
    print(f"{len(cases)} cases\n")
    for (m, c), names in sorted(tally.items(), key=lambda kv: -len(kv[1])):
        print(f"{len(names):4d}  mesh {m:28s} centerline {c:44s}")
        if args.show:
            print(f"      {', '.join(names[:args.show])}"
                  f"{' ...' if len(names) > args.show else ''}")

    dupes = [n for n, k in Counter(c[0] for c in cases).items() if k > 1]
    if dupes:
        print(f"\nduplicate case names: {', '.join(dupes)}")

    if not args.verify:
        print("\n--verify loads each pair and checks the centerline is inside the mesh")
        return 0

    # The only test that settles a pairing. A file name is a hint; a branch
    # centerline sitting outside the surface is proof the two are different exports.
    import branch_cut as BC
    print()
    bad_cases = []
    for name, csv, mat in cases:
        try:
            _, tree = build(BC, csv, mat)
            bad, tot = mismatched(tree)
        except Exception as e:
            print(f"  {name:16s} FAILED TO LOAD: {type(e).__name__}: {e}")
            bad_cases.append(name)
            continue
        if tot and bad / tot > 0.25:
            print(f"  {name:16s} {bad}/{tot} segments outside the mesh  <-- MISMATCH")
            bad_cases.append(name)
    print(f"\n{len(cases) - len(bad_cases)}/{len(cases)} pairs verified"
          f"{'' if not bad_cases else ': ' + ', '.join(bad_cases)}")
    return 1 if bad_cases else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "score", "cases"):
        p = sub.add_parser(name)
        p.add_argument("--show", type=int, default=25)
        if name == "cases":
            p.add_argument("--verify", action="store_true",
                           help="load every pair and check the centerline is inside")
        else:
            p.add_argument("--out", default=str(DEFAULT_OUT))
        if name == "run":
            p.add_argument("--baseline", default="main",
                           help="git revision to compare against")
            p.add_argument("--limit", type=int, default=10_000)
            p.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    return {"run": cmd_run, "score": cmd_score, "cases": cmd_cases}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
