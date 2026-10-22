"""Run the nub-in preset over the thoracic cohort, and score what came out.

    python tools/nub_cohort.py plan  --k 0.5        # plan only, no cutting (fast)
    python tools/nub_cohort.py sweep               # plan at several k, pick one
    python tools/nub_cohort.py run   --k 0.5       # cut, write meshes + sidecars
    python tools/nub_cohort.py report              # read the last run

PLAN is cheap and answers the question that decides k: how much runway does each
vessel actually have? A nub has to sit past the aortic wall and before whatever ends
the vessel -- its tip, or, for a trunk like the brachiocephalic, its split. k_max is
that runway in aortic radii, and the cohort's WORST k_max is the largest k every
vessel can take.

RUN cuts and scores each opening the way tools/cut_bench.py does, because a nub that
tore is worse than no nub: the whole point of the preset is a set of geometries that
can be measured without anyone checking them by hand first.

    flat     the nub's rim lies on the plane that made it (<= 0.5 mm)
    foreign  surface taken from a vessel this cut was not cutting (<= 1 mm2)
    proud    how far the nub actually stands past the aortic wall

Only the thoracic cohort. Meshes and centerlines only; results are written here,
never into the cohort.
"""
import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import branch_cut as BC                                            # noqa: E402
import nub_cut as NC                                               # noqa: E402
from cut_bench import CL_PREF, MESH_PREF, _first, build, cohort_roots   # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_JSONL = HERE / "nub_last.jsonl"
DEFAULT_OUT = ROOT / "output" / "nubs"

# The thoracic cohort. cohort_roots() also carries the abdominal cases, which have no
# arch vessels and are not what this preset is for.
THORACIC = "GT"

FLAT_MM = 0.5
FOREIGN_OK_MM2 = 1.0


def thoracic_cases(limit=None, only=None):
    """Cases under the thoracic root only.

    Deliberately not cut_bench.find_cases(), which walks every root in
    cohorts.local.json. This preset has no business in the abdominal cohort and no
    reason to be listing it, so the abdominal root is never opened at all.
    """
    roots = [r for r in cohort_roots() if r.name == THORACIC]
    if not roots:
        raise SystemExit(f"no {THORACIC} root in tools/cohorts.local.json")
    cases = []
    for root in roots:
        if not root.exists():
            raise SystemExit(f"cohort root missing: {root}")
        for case in sorted(p for p in root.iterdir() if p.is_dir()):
            mat, csv = _first(case, MESH_PREF), _first(case, CL_PREF)
            if mat and csv:
                cases.append((f"{root.name}/{case.name}", csv, mat))
                continue
            for sub in sorted(p for p in case.iterdir() if p.is_dir()):
                mat, csv = _first(sub, MESH_PREF), _first(sub, CL_PREF)
                if mat and csv:
                    cases.append((f"{root.name}/{sub.name}", csv, mat))
    if only:
        want = {o.strip() for o in only.split(",")}
        cases = [c for c in cases if c[0].split("/")[-1] in want or c[0] in want]
    return cases[:limit] if limit else cases


# What the desktop tool does to a jagged centerline, so the batch does it too. The
# window probes the line as supplied, and if its turning exceeds JAGGED_DEG_PER_MM it
# TICKS ITS OWN SMOOTHING BOX and reloads at this many millimetres. Six of the thirty
# thoracic cases are jagged, so on those six the window was cutting a different
# centerline from the one the batch scored -- same offsets, different planes, and on
# three of them the window then refused to export a cut the batch called clean.
GUI_SMOOTH_MM = 1.0

# A thoracic aorta has TWO trunk ends, a root and a distal end; an abdominal one has a
# single inlet. The tool reads this from the relative size of the largest openings and
# gets it wrong on KK24_1 and KK40_3, calling them abdominal. That is not only a naming
# error: the count decides which terminals are trunks, so on those two ONE AORTIC END
# WAS BEING CLASSIFIED AS A BRANCH -- and therefore given a nub instead of the 8 mm
# horizontal cut every other aortic end gets. It also renamed their arch vessels after
# celiac anatomy, which is the visible half of the same mistake.
#
# DETECTED, not declared. Forcing two was a mistake: it was papering over centerlines
# that have since been re-exported, and on a case that genuinely has one aortic end in
# the traced volume it PROMOTES A BRANCH TO A TRUNK -- KK40_3's second-largest opening
# is a 10.8 mm-radius vessel in the middle of the model against a 24.9 mm aorta, and
# forcing two ends made the tool try to square that branch off with a horizontal plane
# instead of nubbing it, so it would not cut at all.
#
# The tool's own test is the size ratio of the two largest openings, and on the current
# centerlines it separates them cleanly: 28 cases sit at 1.02-1.64 (two comparable
# aortic ends) and the two one-ended cases at 2.11 and 2.30. Detected per case, then
# RECORDED, so the window classifies each case exactly as the batch did rather than
# re-deciding.
TRUNK_ENDS = None

# Where the tool's own test is known to be wrong, say so per case rather than let it
# guess. EMPTY, and that is the point: it held KK40_3 until the test itself was fixed.
#
# KK40_3 is plainly a thoracic arch -- ascending, four arch branches, descending -- but
# its descending is dilated to r = 24.9 mm against a normal ascending at r = 10.8. The
# old test compared those two openings' calibre, read the pair as one dominant inlet,
# and named the arch vessels after celiac anatomy. Declaring the answer here was a
# patch over a rule that could not see an aneurysm.
#
# classify_terminals now asks a structural question instead -- does the trunk END by
# dividing, as an abdominal aorta does at the iliacs (see PAIR_AREA_FRACTION) -- which
# an aneurysm does not disturb. KK40_3 reads 0.17 against a 0.25 cut-off and comes out
# thoracic unaided, so the declaration is gone rather than left in agreeing with it.
#
# KP35_1 is the case that still needs it, and it is the tool's one remaining family
# error over 74 cases (0.31 against a 0.25 cut-off -- a near call the window now
# reports as one). Its arch carries FOUR first-order vessels, so the outflow divides
# more evenly there than an arch usually does and it reads as one-ended. Cut that
# way, its second aortic end was nubbed as a branch: the exported mesh carries a
# 40.4 mm-diameter opening labelled "Left gastric artery", which is not a thing.
#
# Confirmed by eye on the rendered centerline before being written down. Add a case
# here only on the same standard.
ANATOMY = {
    "KP35_1": 2,          # thoracic arch, four first-order vessels; reads one-ended
}


def detected_trunk_ends(tree):
    return sum(1 for t in tree.terminals if getattr(t, "is_trunk", False))


def centerline_smoothing(csv):
    """The smoothing the desktop tool would apply to this centerline, in mm."""
    from branch_tree import load_tree
    probe = load_tree(str(csv), smooth_mm=0.0)
    return GUI_SMOOTH_MM if probe.jagged else 0.0


def case_id(csv):
    """The case folder name, which is what ANATOMY is keyed by."""
    parts = Path(str(csv)).parts
    return parts[-3] if len(parts) >= 3 else ""


def prepare(csv, mat, verbose=False):
    smooth = centerline_smoothing(csv)
    declared = ANATOMY.get(case_id(csv))
    mesh, tree = build(BC, csv, mat, smooth_mm=smooth,
                       n_trunk_ends=declared if declared else TRUNK_ENDS)
    mesh = NC.drop_debris(mesh, verbose=verbose)
    if declared and verbose:
        print(f"  anatomy declared: {declared} trunk end(s)")
    if verbose and smooth:
        print(f"  centerline is jagged, smoothed {smooth:.1f} mm as the window would")
    # Before anything reads a radius. Seven of the thirty cases arrive with their
    # centerline broken at a non-junction, and the short piece left at the end of the
    # model measures its radius against the mesh's end cap instead of its wall.
    NC.heal_split_segments(tree, mesh, verbose=verbose)
    tnames, vnames, _ = BC.name_anatomy(tree, mesh)
    tree._nub_smooth_mm = float(smooth)
    tree._nub_trunk_ends = detected_trunk_ends(tree)
    return mesh, tree, tnames, vnames


def plan_case(name, csv, mat, ks):
    """One row per vessel per k. No cutting, so this is load + ownership only."""
    mesh, tree, tnames, vnames = prepare(csv, mat)
    jg = BC._junction_graph(tree)
    route = BC.aorta_route(tree, jg)
    ref = NC.aorta_reference(mesh, tree, route=route, jg=jg)
    vessels = NC.first_order_vessels(tree, jg=jg, route=route, mesh=mesh)

    rows = []
    for k in ks:
        plan, _used = NC.nub_plan(mesh, tree, k=k, ref=ref, vessels=vessels,
                                  names=vnames)
        for r in plan:
            rows.append({"case": name, "k": float(k), **_clean(r)})
    return rows, len(vessels)


def _clean(row):
    """Plan rows carry numpy scalars; JSON does not."""
    out = {}
    for key, val in row.items():
        if isinstance(val, (np.floating, np.integer)):
            val = val.item()
        if isinstance(val, float) and not np.isfinite(val):
            val = None
        if isinstance(val, set):
            val = sorted(val)
        out[key] = val
    return out


def score_openings(mesh_out, tree, offsets, terms, plan):
    """Measure every nub the cut actually produced: rim flatness and stray area."""
    loops = BC._loops_by_terminal(mesh_out, tree, offsets, terminals=terms)
    verts = np.asarray(mesh_out.vertices)
    by_seg = {r["segment"]: r for r in plan}

    for t in terms:
        row = by_seg.get(t.segment.name)
        if row is None or t.label not in offsets:
            continue
        off = float(offsets[t.label])
        pt, nrm = t.cut_plane(off)
        loop = loops.get(t.label)
        if loop is None or len(loop) == 0:
            row["flat_mm"] = None
            row["cut_ok"] = False
            continue
        ring = verts[np.asarray(loop)]
        dev = float(np.abs((ring - pt) @ nrm).max())
        row["flat_mm"] = dev
        row["rim_points"] = int(len(ring))
        row["cut_ok"] = bool(dev <= FLAT_MM)
        # Where the nub's mouth ended up and how wide it is. The reviewer's window
        # reads these straight out of the sidecar rather than re-deriving them from
        # the STL, so what is labelled on screen is what was scored here.
        row["rim_centroid"] = [float(x) for x in ring.mean(axis=0)]
        row["plane_point"] = [float(x) for x in pt]
        row["plane_normal"] = [float(x) for x in nrm]
        row["mouth_ead_mm"], row["own_plane_mm"] = _mouth(ring, nrm)
    return plan


def _mouth(ring, nrm):
    """The nub's mouth measured from the RIM, and how planar that rim is.

    Not by sectioning the cut mesh: a plane laid on an open mesh's own rim slices it
    tangentially and returns whatever else it crosses -- it reported 83 mm mouths on
    6 mm branches, because the plane also cut the aorta and the section came back as
    an open polyline whose shoelace area is meaningless. The rim IS the mouth, and it
    is already in hand.

    The second number is the rim's spread about its OWN best-fit plane. It separates
    a torn opening from a badly ORIENTED one: a rim that is 0.6 mm off the cut plane
    but 0.02 mm off its own is a flat disc that the plane normal does not quite
    match, which is a different complaint from a rim that is genuinely ragged.
    """
    c = ring.mean(axis=0)
    q = ring - c
    a = np.array([0.0, 0.0, 1.0])
    if abs(float(nrm @ a)) > 0.9:
        a = np.array([1.0, 0.0, 0.0])
    u = np.cross(nrm, a)
    u = u / max(np.linalg.norm(u), 1e-12)
    v = np.cross(nrm, u)
    x, y = q @ u, q @ v
    # boundary_loops returns a loop's vertices sorted by INDEX, not walked around the
    # ring, so they have to be put in order before any area is taken from them. A rim
    # encircles its vessel axis, so polar angle about the centroid is that order.
    o = np.argsort(np.arctan2(y, x))
    x, y = x[o], y[o]
    area = 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))
    ead = float(2.0 * np.sqrt(area / np.pi)) if area > 1e-9 else None

    # best-fit plane: the smallest singular direction of the centred rim
    try:
        n_fit = np.linalg.svd(q, full_matrices=False)[2][-1]
        own = float(np.abs(q @ n_fit).max())
    except Exception:                                              # noqa: BLE001
        own = None
    return ead, own


def run_case(name, csv, mat, k, outdir, write_mesh=True):
    mesh, tree, tnames, vnames = prepare(csv, mat)
    t0 = time.time()
    out, plan, notes, info = NC.cut_nubs(mesh, tree, k=k, anat=vnames, verbose=False)
    secs = time.time() - t0

    plan = score_openings(out, tree, info["offsets"], info["terminals"], plan)

    loops = BC.boundary_loops(out)
    case_id = name.split("/")[-1]
    # Exactly what was cut, terminal by terminal, so the desktop tool can be opened on
    # this case showing THIS cut rather than a re-derivation of it. A rescue rewrites
    # an offset and an obliquity fallback rewrites a mode, so neither is recoverable
    # from the rule and k alone.
    cuts = []
    for t in info["terminals"]:
        if t.label not in info["offsets"]:
            continue
        cuts.append({
            "label": t.label, "segment": t.segment.name, "name": t.name,
            "offset_mm": float(info["offsets"][t.label]), "mode": t.mode,
            "interior": bool(getattr(t, "is_interior", False)),
            "trunk": bool(getattr(t, "is_trunk", False)),
            "band_mm": (None if getattr(t, "band_mm", None) is None
                        else float(t.band_mm)),
            # Where this opening IS. Names are not a safe key: the loader's weld step
            # splits a host segment into <name>_A/<name>_B when it attaches a loose
            # piece mid-segment, and it did that in the window and not in the batch on
            # KP35_1 -- Line14 against Line14_A. The match then missed and the
            # ascending aorta silently kept a 6 mm default that falls outside the mesh.
            # A tip does not get renamed.
            "tip": [float(x) for x in t.tip],
            "tilt_a_deg": float(getattr(t, "tilt_a_deg", 0.0) or 0.0),
            "tilt_b_deg": float(getattr(t, "tilt_b_deg", 0.0) or 0.0),
            "subtree": sorted(getattr(t, "subtree_segments", None) or []),
        })

    rec = {"case": name, "k": float(k), "secs": round(secs, 1),
           "csv": str(csv), "mat": str(mat),
           "smooth_mm": float(getattr(tree, "_nub_smooth_mm", 0.0)),
           "n_trunk_ends": int(getattr(tree, "_nub_trunk_ends", 0)) or None,
           "faces_in": int(len(mesh.faces)), "faces_out": int(len(out.faces)),
           "openings": len(loops), "euler": int(out.euler_number),
           "validate": info["validate"], "notes": notes, "cuts": cuts,
           "vessels": [_clean(r) for r in plan]}

    if write_mesh:
        outdir.mkdir(parents=True, exist_ok=True)
        stem = f"{case_id}_nub_k{k:g}"
        out.export(str(outdir / f"{stem}.stl"))
        (outdir / f"{stem}.json").write_text(json.dumps(rec, indent=2))
        rec["mesh"] = f"{stem}.stl"
    return rec


def cmd_plan(args):
    ks = [float(x) for x in args.k.split(",")]
    cases = thoracic_cases(args.limit, args.only)
    print(f"{len(cases)} thoracic cases, k = {ks}")
    rows = []
    for i, (name, csv, mat) in enumerate(cases, 1):
        try:
            got, nv = plan_case(name, csv, mat, ks)
        except Exception as exc:                                   # noqa: BLE001
            print(f"  [{i}/{len(cases)}] {name}: FAILED {exc}")
            if args.traceback:
                traceback.print_exc()
            continue
        rows.extend(got)
        worst = min((r["k_max"] for r in got
                     if r.get("k_max") is not None), default=float("nan"))
        print(f"  [{i}/{len(cases)}] {name}: {nv} vessels, worst k_max {worst:.2f}")
    Path(args.out).write_text("\n".join(json.dumps(r) for r in rows))
    print(f"\nwrote {len(rows)} rows to {args.out}")
    summarise_plan(rows, ks)
    return 0


def summarise_plan(rows, ks):
    print(f"\n{'k':>6} {'vessels':>8} {'fits':>7} {'%':>6} {'median proud':>13} "
          f"{'worst k_max':>12}")
    print("-" * 58)
    for k in ks:
        sub = [r for r in rows if abs(r["k"] - k) < 1e-9]
        if not sub:
            continue
        fits = [r for r in sub if r["ok"]]
        proud = [r["proud_mm"] for r in fits if r.get("proud_mm") is not None]
        kmax = [r["k_max"] for r in sub if r.get("k_max") is not None]
        print(f"{k:6g} {len(sub):8d} {len(fits):7d} {100 * len(fits) / len(sub):5.0f}% "
              f"{np.median(proud) if proud else float('nan'):12.1f}  "
              f"{min(kmax) if kmax else float('nan'):12.2f}")

    tight = sorted(((r["k_max"], r["case"], r["name"]) for r in rows
                    if r.get("k_max") is not None and abs(r["k"] - ks[0]) < 1e-9))[:12]
    if tight:
        print("\ntightest vessels in the cohort (these set the ceiling on k):")
        for kmax, case, nm in tight:
            print(f"  {kmax:5.2f}  {case:16s} {nm}")


def cmd_run(args):
    cases = thoracic_cases(args.limit, args.only)
    outdir = Path(args.outdir)
    print(f"{len(cases)} thoracic cases at k = {args.k}  ->  {outdir}")
    recs = []
    for i, (name, csv, mat) in enumerate(cases, 1):
        try:
            rec = run_case(name, csv, mat, float(args.k), outdir,
                           write_mesh=not args.no_mesh)
        except Exception as exc:                                   # noqa: BLE001
            print(f"  [{i}/{len(cases)}] {name}: FAILED {exc}")
            if args.traceback:
                traceback.print_exc()
            recs.append({"case": name, "k": float(args.k), "error": str(exc),
                         "vessels": []})
            continue
        nubs = [v for v in rec["vessels"] if v.get("cut_ok")]
        bad = [v for v in rec["vessels"] if v["ok"] and not v.get("cut_ok")]
        grown = [v for v in rec["vessels"] if v.get("rescued_mm")]
        extra = ""
        if bad:
            extra += f", {len(bad)} torn"
        if grown:
            extra += f", {len(grown)} grown"
        dropped = [v for v in rec["vessels"] if not v["ok"] and v.get("s_wall_mm")]
        if dropped:
            extra += f", {len(dropped)} uncut"
        print(f"  [{i}/{len(cases)}] {name}: {len(nubs)} clean nubs{extra}"
              f"  {rec['secs']}s")
        recs.append(rec)
    Path(args.out).write_text("\n".join(json.dumps(r) for r in recs))
    print(f"\nwrote {len(recs)} cases to {args.out}")
    report(recs)
    return 0


def report(recs):
    vessels = [v for r in recs for v in r.get("vessels", [])]
    planned = [v for v in vessels if v["ok"]]
    clean = [v for v in planned if v.get("cut_ok")]
    rescued = [v for v in planned if v.get("rescued_mm")]
    skipped = [v for v in vessels if not v["ok"]]
    proud = [v["proud_mm"] for v in clean if v.get("proud_mm") is not None]

    print(f"\ncases {len(recs)}   vessels {len(vessels)}")
    print(f"  planned          {len(planned)}")
    print(f"  cut clean        {len(clean)}"
          f"  ({100 * len(clean) / max(len(planned), 1):.0f}% of planned)")
    print(f"  grown to cut     {len(rescued)}")
    clamped = [v for v in planned if v.get("clamped")]
    marginal = [v for v in planned if v.get("marginal")]
    if clamped:
        print(f"  shortened        {len(clamped)}  (vessel too short for the full k)")
    if marginal:
        print(f"  barely a nub     {len(marginal)}  (under "
              f"{NC.MIN_USEFUL_PROUD_MM:g} mm proud -- look at these)")
    print(f"  not placed       {len(skipped)}")
    if proud:
        print(f"  proud length     median {np.median(proud):.1f} mm, "
              f"range {min(proud):.1f}-{max(proud):.1f} mm")
    why = {}
    for v in skipped:
        key = (v.get("why") or "").split("--")[0].split("(")[0].strip()[:60]
        why[key] = why.get(key, 0) + 1
    for key, n in sorted(why.items(), key=lambda kv: -kv[1]):
        print(f"    {n:3d}  {key}")


def cmd_kfit(args):
    """What k the cohort will actually accept, read off a run rather than guessed.

    A nub is cut at the length asked for only if a plane there produces a rim that
    lies on it. Where it does not, the run grows the nub until one does, and THAT
    length is a measurement: the shortest stub this vessel can have and still be cut
    cleanly. It is not arbitrary -- a branch leaves its aorta as a funnel and only
    becomes a tube some way out, and a plane through a funnel cannot come back flat.

    So the growth is the answer to "how long should a nub be". Ask for less than a
    vessel can give and the run silently gives you more, and the uniform stub length
    the preset is for is lost. Pick k at a high percentile of what the cohort needed
    and most vessels are cut at the length asked for.
    """
    recs = [json.loads(l) for l in Path(args.out).read_text().splitlines() if l.strip()]
    vs = [v for r in recs for v in r.get("vessels", [])
          if v.get("ok") and v.get("cut_ok")]
    if not vs:
        raise SystemExit(f"no clean nubs in {args.out}")

    # k the vessel actually ended up with, whether it was asked for or grown into
    need = []
    for v in vs:
        r_ao, proud = v.get("r_aorta_mm"), v.get("proud_mm")
        if r_ao and proud is not None and r_ao > 0:
            need.append((proud / r_ao, proud, v))
    need.sort()
    ks = np.array([n[0] for n in need])
    proud = np.array([n[1] for n in need])
    grown = [v for v in vs if v.get("rescued_mm")]

    print(f"{len(vs)} clean nubs over {len(recs)} cases   "
          f"asked for k={recs[0]['k']:g}")
    print(f"  cut at the length asked for   {len(vs) - len(grown)}")
    print(f"  had to be grown               {len(grown)}"
          f"  ({100 * len(grown) / len(vs):.0f}%)")
    print(f"\nk each vessel actually needed (proud length / aortic radius):")
    for p in (50, 75, 90, 95, 100):
        print(f"  {p:3d}th percentile   k = {np.percentile(ks, p):.2f}"
              f"   ({np.percentile(proud, p):5.1f} mm proud)")
    print(f"\nchoosing k at the 90th percentile means ~90% of vessels are cut at the "
          f"length\nasked for and the rest are grown a little past it.")

    worst = need[-8:][::-1]
    print("\nvessels that needed the most:")
    for kk, pr, v in worst:
        print(f"  k={kk:5.2f}  {pr:5.1f} mm proud   {v['name'][:44]:46s}"
              f"{'  grown' if v.get('rescued_mm') else ''}")
    return 0


def cmd_report(args):
    recs = [json.loads(l) for l in Path(args.out).read_text().splitlines() if l.strip()]
    if recs and "vessels" in recs[0] and "case" in recs[0] and "secs" in recs[0]:
        report(recs)
    else:
        ks = sorted({r["k"] for r in recs})
        summarise_plan(recs, ks)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("plan")
    p.add_argument("--k", default="0.1,0.2,0.3,0.5,0.75,1.0")
    p.add_argument("--out", default=str(HERE / "nub_plan.jsonl"))
    p.add_argument("--limit", type=int)
    p.add_argument("--only")
    p.add_argument("--traceback", action="store_true")

    p = sub.add_parser("sweep")
    p.add_argument("--k", default="0.1,0.2,0.3,0.5,0.75,1.0")
    p.add_argument("--out", default=str(HERE / "nub_plan.jsonl"))
    p.add_argument("--limit", type=int)
    p.add_argument("--only")
    p.add_argument("--traceback", action="store_true")

    p = sub.add_parser("run")
    p.add_argument("--k", default="0.5")
    p.add_argument("--out", default=str(DEFAULT_JSONL))
    p.add_argument("--outdir", default=str(DEFAULT_OUT))
    p.add_argument("--limit", type=int)
    p.add_argument("--only")
    p.add_argument("--no-mesh", action="store_true")
    p.add_argument("--traceback", action="store_true")

    p = sub.add_parser("report")
    p.add_argument("--out", default=str(DEFAULT_JSONL))

    p = sub.add_parser("kfit")
    p.add_argument("--out", default=str(DEFAULT_JSONL))

    args = ap.parse_args()
    if args.cmd == "sweep":
        args.cmd = "plan"
    return {"plan": cmd_plan, "run": cmd_run, "report": cmd_report,
            "kfit": cmd_kfit}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
