"""Run a folder of cut meshes through CAMP and record delta K against area.

    python tools/camp_cohort.py run                  # every mesh, resumable
    python tools/camp_cohort.py run --workers 6
    python tools/camp_cohort.py plot                 # read the results, draw the curve

WHAT IT DOES, and the settings that matter

Each mesh goes weld -> ACVD uniform remesh -> CAMP partition -> area-weighted
variance of integrated Gaussian curvature. That is shape_metrics.delta_k, which is
the same path the website uses, so these numbers and an upload's are comparable.

Two settings are pinned here rather than defaulted, because delta K is EXTENSIVE
and partition-dependent -- cut extent, face count and trim all move it, and a
number computed at a different operating point cannot be plotted on the same axes:

    target_faces = 100,000     matches the reference cohort, which was fixed-FACES
                               and not fixed-density (CAMP_INTEGRATION_GUIDE Q4)
    trim_frac    = 0.0         no percentile trimming of per-patch integrated K

Everything else is CAMP's default (R2>=0.9, eta<=0.075, min_patch_size=18,
seed_patches=100). The partition itself does not depend on the trim, so a run at
trim 0 can be re-trimmed afterwards from the stored per-patch table without
re-running -- but the stored k_fluct is the trim-0 one.

Runtime is ~200-220 s per mesh single-core, so this is a queued job. Workers are
processes, each pinned to one thread: CAMP's split loop is single-core anyway and
letting BLAS spawn 12 threads inside 6 workers makes it slower, not faster.
"""
import argparse
import json
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

HERE = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "output" / "nubs_final"
DEFAULT_OUT = HERE / "camp_nubs_100k_trim0.jsonl"

TARGET_FACES = 100_000
TRIM_FRAC = 0.0


FOLD_DEG = 90.0          # a dihedral past this is the surface doubling back
REPAIR_RINGS = 2         # how far the repair blends outward from a folded edge
REPAIR_ITERS = 12


def repair_folds(mesh, fold_deg=FOLD_DEG, rings=REPAIR_RINGS, iters=REPAIR_ITERS):
    """Smooth ONLY the vertices caught in folded edges, and report what moved.

    A fold is a segmentation defect -- the surface doubles back on itself -- and a
    quadratic fitted across one returns a curvature that is not a curvature. On
    KP35_1 that produced a patch at K = -41.7 carrying 54% of the case's delta K,
    against a next-worst of -5.4; the folds are present in the source .mat and are
    not made by cutting.

    The repair is deliberately tiny and local: seed on the vertices of folded
    edges, grow a couple of rings so it blends, hold the outermost ring fixed as a
    seam, and Laplacian-smooth the interior. On KP35_1 that moves 0.73% of the
    vertices by 0.53 mm on average and changes total surface area by 0.06%, while
    delta K falls 1.442 -> 0.024 with no trimming at all. Fixing the mesh where it
    is broken beats discarding patches everywhere it is not.

    Returns (repaired_mesh, info).
    """
    import numpy as np
    import trimesh

    V0 = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces)
    fa = np.asarray(mesh.face_adjacency)
    ang = np.degrees(np.asarray(mesh.face_adjacency_angles))
    folded = fa[ang > fold_deg]
    info = {"folded_edges": int(len(folded)), "max_dihedral_before": float(
        ang.max()) if len(ang) else 0.0, "moved_vertices": 0,
        "max_move_mm": 0.0, "mean_move_mm": 0.0, "area_change_pct": 0.0,
        "folded_after": int(len(folded)), "max_dihedral_after": float(
            ang.max()) if len(ang) else 0.0}
    if not len(folded):
        return mesh, info

    adj = mesh.vertex_neighbors
    sel = set(np.unique(F[folded.ravel()]).tolist())
    for _ in range(rings):
        for v in list(sel):
            sel.update(adj[v])
    free = np.zeros(len(V0), bool)
    free[np.array(sorted(sel))] = True
    # The outer ring stays put, so the repaired area is stitched to what surrounds
    # it. Decided against a SNAPSHOT: testing free[] while clearing free[] in the
    # same pass lets each removal disqualify its neighbours in turn, and the region
    # collapses -- it left 10 of 182 vertices free on KP35_1 and made KY160_1 and
    # KR52_1 measurably worse rather than better.
    interior = free.copy()
    for v in np.flatnonzero(interior):
        if not all(interior[n] for n in adj[v]):
            free[v] = False

    V = V0.copy()
    for _ in range(iters):
        nxt = V.copy()
        for v in np.flatnonzero(free):
            nb = adj[v]
            if len(nb):
                nxt[v] = V[nb].mean(axis=0)
        V = nxt

    out = trimesh.Trimesh(vertices=V, faces=F, process=False)
    moved = np.linalg.norm(V - V0, axis=1)
    ang2 = np.degrees(np.asarray(out.face_adjacency_angles))
    info.update(moved_vertices=int(free.sum()),
                max_move_mm=float(moved.max()),
                mean_move_mm=float(moved[moved > 0].mean()) if (moved > 0).any()
                else 0.0,
                area_change_pct=float(100.0 * (out.area / mesh.area - 1.0)),
                folded_after=int((ang2 > fold_deg).sum()),
                max_dihedral_after=float(ang2.max()) if len(ang2) else 0.0)
    return out, info


def _one(path_str, target_faces, trim_frac, repair=False):
    """One mesh, in its own process with BLAS pinned to a single thread."""
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = "1"
    sys.path.insert(0, str(ROOT))
    import time

    import numpy as np
    import trimesh

    import shape_metrics as SM

    path = Path(path_str)
    t0 = time.time()
    mesh = trimesh.load(str(path), process=False)
    # Topology BEFORE anything touches it, so a remesh that opens or closes a hole
    # is an error rather than a quietly different number. delta_k asserts on it.
    welded = SM.weld(mesh)
    chi = int(welded.euler_number)
    fix = {}
    if repair:
        welded, fix = repair_folds(welded)
        if int(welded.euler_number) != chi:
            raise RuntimeError(f"fold repair changed euler {chi} -> "
                               f"{welded.euler_number}")
        mesh = welded
    res = SM.delta_k(mesh, target_faces=target_faces, trim_frac=trim_frac,
                     make_figure=False, expect_euler=chi)
    # Keep the per-patch table. The PARTITION does not depend on the trim, so with
    # these two arrays every other trim is a few microseconds of arithmetic instead
    # of another 45-second CAMP run -- which is what makes "is this scatter real or
    # is it the trim?" a question worth asking.
    stem = path.name.replace("_nub.stl", "").replace(".stl", "")
    patches = Path(path_str).parent.parent / "camp_patches"
    patches.mkdir(parents=True, exist_ok=True)
    # Face count in the name: the partition depends on it, so two densities are two
    # different tables and overwriting one with the other loses the comparison.
    suffix = "" if target_faces == TARGET_FACES else f"_{target_faces // 1000}k"

    # How far each patch sits from a CUT RIM. A patch straddling an opening has a
    # one-sided neighbourhood, so its local quadratic is fitted against half a
    # surface and its <kg> is not a curvature measurement -- on KP35_1 one such
    # patch carries 54% of delta K with R2 = 0.35. The rims are boundaries THIS TOOL
    # created, so excluding a collar around them costs no anatomy; storing the
    # distance lets the collar width be chosen later without re-running CAMP.
    from scipy.spatial import cKDTree
    rm = res["mesh"]
    edges = np.asarray(rm.edges_sorted).reshape(-1, 2)
    uq, cnt = np.unique(edges, axis=0, return_counts=True)
    rim_pts = np.asarray(rm.vertices)[np.unique(uq[cnt == 1])]
    labels = np.asarray(res["face_labels"])
    centres = np.asarray(rm.triangles_center, float)
    n_patch = len(np.asarray(res["integrated_K"]))
    rim_dist = np.full(n_patch, np.inf)
    if len(rim_pts):
        kd = cKDTree(rim_pts)
        d_face = kd.query(centres)[0]
        for p in range(n_patch):
            sel = labels == p
            if sel.any():
                rim_dist[p] = float(d_face[sel].min())

    # The worst dihedral angle inside each patch. A fold is a surface doubling back
    # on itself, and a quadratic fitted across one is not measuring curvature --
    # KP35_1 carries 28 of its 46 folded edges within 8 mm of the patch that holds
    # 54% of its delta K, and those folds are present in the SOURCE segmentation,
    # unchanged by cutting. This is a defect flag, not a shape measurement.
    fa = np.asarray(rm.face_adjacency)
    fang = np.degrees(np.asarray(rm.face_adjacency_angles))
    max_fold = np.zeros(n_patch)
    lab_a, lab_b = labels[fa[:, 0]], labels[fa[:, 1]]
    for p in range(n_patch):
        inside = (lab_a == p) & (lab_b == p)
        if inside.any():
            max_fold[p] = float(fang[inside].max())

    np.savez_compressed(patches / f"{stem}{suffix}.npz",
                        integrated_K=np.asarray(res["integrated_K"], float),
                        patch_area=np.asarray(res["patch_area"], float),
                        avg_kg=np.asarray(res["avg_kg"], float),
                        eta_patch=np.asarray(res["eta_patch"], float),
                        R2_patch=np.asarray(res["R2_patch"], float),
                        rim_dist=rim_dist, max_fold=max_fold)
    return {
        "case": path.name.replace("_nub.stl", "").replace(".stl", ""),
        "file": path.name,
        "k_fluct": res["k_fluct"],
        "A_mm2": res["A_mm2"],
        "n_patches": res["n_patches"],
        "n_kept": res["n_kept"],
        "pct_at_min_size": res["pct_at_min_size"],
        "pct_failed_quality": res["pct_failed_quality"],
        "euler_number": res["euler_number"],
        "euler_input": chi,
        "faces": res["faces"],
        "target_faces": res["target_faces"],
        "trim_frac": res["trim_frac"],
        "density_verts_per_mm2": res["density_verts_per_mm2"],
        "in_reference_range": res["in_reference_range"],
        "seconds": round(time.time() - t0, 1),
        **{f"repair_{k}": v for k, v in fix.items()},
    }


def cmd_run(args):
    meshes = sorted(Path(args.indir).glob("*.stl"))
    if not meshes:
        raise SystemExit(f"no .stl under {args.indir}")
    faces = int(getattr(args, "faces", None) or TARGET_FACES)
    out_path = Path(args.out)
    done = {}
    if out_path.exists() and not args.fresh:
        for ln in out_path.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done[r["file"]] = r
    todo = [m for m in meshes if m.name not in done]
    print(f"{len(meshes)} mesh(es), {len(done)} already done, {len(todo)} to run")
    print(f"  target_faces={faces:,}  trim_frac={TRIM_FRAC}  "
          f"workers={args.workers}")
    if not todo:
        return cmd_plot(args)

    fails = []
    with out_path.open("a", encoding="utf-8") as fh, \
            ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_one, str(m), faces, TRIM_FRAC,
                            bool(getattr(args, "repair", False))): m for m in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            m = futs[fut]
            try:
                rec = fut.result()
            except Exception as e:                              # noqa: BLE001
                fails.append((m.name, f"{type(e).__name__}: {e}"))
                print(f"  [{i}/{len(todo)}] {m.name}: FAILED {type(e).__name__}: "
                      f"{str(e)[:110]}", flush=True)
                if args.verbose:
                    traceback.print_exc()
                continue
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print(f"  [{i}/{len(todo)}] {rec['case']:11s} "
                  f"A={rec['A_mm2']:9.0f} mm2  dK={rec['k_fluct']:.4g}  "
                  f"{rec['n_patches']:5d} patches  chi={rec['euler_number']:3d}  "
                  f"{rec['seconds']:.0f}s", flush=True)
    if fails:
        print(f"\n{len(fails)} failed:")
        for name, why in fails:
            print(f"  {name}: {why[:150]}")
    return cmd_plot(args)


def _load(out_path):
    """Results, one per mesh, latest wins.

    The file is append-only so a run can be interrupted and resumed without losing
    what it already computed. That means a --fresh re-run leaves two records for
    every mesh, and averaging or double-plotting them would be silent -- so the last
    record for each file is the one that counts.
    """
    latest = {}
    for ln in Path(out_path).read_text(encoding="utf-8").splitlines():
        if ln.strip():
            r = json.loads(ln)
            latest[r["file"]] = r
    return sorted(latest.values(), key=lambda r: r["case"])


def cmd_plot(args):
    import numpy as np

    rows = _load(args.out)
    if not rows:
        raise SystemExit(f"nothing in {args.out}")
    A = np.array([r["A_mm2"] for r in rows], float)
    K = np.array([r["k_fluct"] for r in rows], float)

    # --trim re-reads the stored per-patch tables rather than re-running CAMP: the
    # partition does not depend on the trim, only the variance taken over it does.
    trim = args.trim if getattr(args, "trim", None) is not None else None
    if trim is not None and abs(trim - float(rows[0]["trim_frac"])) > 1e-12:
        patches = Path(args.indir).parent / "camp_patches"
        sys.path.insert(0, str(ROOT))
        import shape_metrics as SM
        cu = SM.load_campsrc()
        vals = []
        for r in rows:
            d = np.load(patches / f"{r['case']}.npz")
            iK, area = d["integrated_K"], d["patch_area"]
            keep = cu._trim_mask_integrated(iK, trim)
            q, a = iK[keep], area[keep]
            w = a / a.sum()
            vals.append(float(np.sum(w * (q - np.sum(w * q)) ** 2)))
        K = np.array(vals, float)
        for r, v in zip(rows, K):
            r["k_fluct"], r["trim_frac"] = v, trim

    # Khabaz-style normalisation: report delta K as a multiple of the NORMAL
    # cohort's mean, so a case reads "1.8x normal" rather than "0.0096". Only the
    # scale changes -- the ordering, the spread and the fit are untouched.
    norm_group = getattr(args, "norm_group", None)
    norm_mean = None
    if norm_group:
        sel = [i for i, r in enumerate(rows) if r["case"].startswith(norm_group)]
        if not sel:
            raise SystemExit(f"no case starts with {norm_group!r}")
        norm_mean = float(np.mean(K[sel]))
        print(f"  normalising by the mean of {len(sel)} {norm_group}* case(s): "
              f"{norm_mean:.6g}")
        K = K / norm_mean
        for r, v in zip(rows, K):
            # Onto k_fluct itself, not a second key: everything below -- the point
            # positions AND the labels beside them -- reads this one field, and
            # updating only half of it puts every label off the axes.
            r["k_fluct_raw"] = r["k_fluct"]
            r["k_fluct"] = float(v)

    # Named for what it CONTAINS. Writing a re-trimmed or normalised table back over
    # camp_nubs_100k_trim0.csv leaves a file whose name says trim 0 and whose numbers
    # are not, which is the kind of thing that ends up in a paper.
    _t = float(rows[0]["trim_frac"])
    _stem = f"camp_nubs_100k_trim{_t:g}".replace(".", "p")
    if norm_mean is not None:
        _stem += f"_norm{norm_group}"
    csv = (Path(args.out).with_suffix(".csv")
           if (_t == TRIM_FRAC and norm_mean is None) else
           Path(args.out).with_name(_stem + ".csv"))
    cols = ["case", "A_mm2", "k_fluct", "n_patches", "pct_at_min_size",
            "euler_number", "faces", "density_verts_per_mm2", "trim_frac"]
    if norm_mean is not None:
        cols.insert(3, "k_fluct_raw")
    lines = [",".join(cols)]
    for r in rows:
        lines.append(",".join(str(r.get(c, "")) for c in cols))
    csv.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # log-log, because delta K is extensive: it scales with area, and the reference
    # is a straight line of slope ~1 in these axes. A linear plot of the same points
    # is a fan that hides where a case actually sits.
    good = np.isfinite(A) & np.isfinite(K) & (A > 0) & (K > 0)
    slope = icept = r2 = float("nan")
    if good.sum() >= 2:
        x, y = np.log10(A[good]), np.log10(K[good])
        slope, icept = np.polyfit(x, y, 1)
        pred = slope * x + icept
        ss_res = float(((y - pred) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    print(f"\n{len(rows)} mesh(es)")
    print(f"  area      {A.min():,.0f} - {A.max():,.0f} mm2  (median {np.median(A):,.0f})")
    print(f"  delta K   {K.min():.4g} - {K.max():.4g}  (median {np.median(K):.4g})")
    print(f"  log-log fit: log10(dK) = {slope:.3f} log10(A) {icept:+.3f}   R2={r2:.3f}")
    print(f"  wrote {csv.name}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not installed, no figure written)")
        return 0

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8))
    for ax, logscale in zip(axes, (False, True)):
        ax.scatter(A, K, s=46, c="#1f77b4", edgecolor="white", linewidth=0.7,
                   zorder=3)
        for r in rows:
            ax.annotate(r["case"], (r["A_mm2"], r["k_fluct"]),
                        textcoords="offset points", xytext=(4, 3), fontsize=6.5,
                        color="#444")
        if logscale and good.sum() >= 2:
            xs = np.linspace(np.log10(A[good].min()), np.log10(A[good].max()), 50)
            ax.plot(10 ** xs, 10 ** (slope * xs + icept), "-", color="#d62728",
                    lw=1.4, zorder=2,
                    label=f"slope {slope:.2f}, $R^2$={r2:.2f}")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.legend(frameon=False, fontsize=9)
        if norm_mean is not None:
            # 1.0 is the normal cohort's mean, so the eye can read "above/below
            # normal" straight off the axis.
            ax.axhline(1.0, color="#888", lw=1.0, ls="--", zorder=1)
        ax.set_xlabel("surface area  $A_T$  (mm$^2$)")
        ax.set_ylabel(
            rf"$\delta K\,/\,\langle\delta K\rangle_{{{norm_group}}}$"
            if norm_mean is not None else
            r"$\delta K$  (k_fluct, area-weighted variance of $K_j$)")
        ax.set_title("log-log" if logscale else "linear", fontsize=10)
        ax.grid(alpha=0.25, zorder=0)
    shown_trim = float(rows[0]["trim_frac"])
    norm_note = (f"  |  normalised by mean {norm_group} (n={len(sel)}), "
                 f"{norm_mean:.4g}" if norm_mean is not None else "")
    fig.suptitle(f"Thoracic nub-in cohort, n={len(rows)}  |  CAMP at "
                 f"{TARGET_FACES:,} faces, trim {shown_trim:g}{norm_note}",
                 fontsize=11)
    fig.tight_layout()
    stem = f"camp_nubs_100k_trim{shown_trim:g}".replace(".", "p")
    if norm_mean is not None:
        stem += f"_norm{norm_group}"
    png = (Path(args.out).with_suffix(".png")
           if (shown_trim == TRIM_FRAC and norm_mean is None) else
           Path(args.out).with_name(stem + ".png"))
    fig.savefig(png, dpi=170)
    print(f"  wrote {png.name}")
    return 0


TRIMS = (0.0, 0.005, 0.01, 0.02, 0.05)

# Patients whose scans define the normal baseline. Everything else is placed
# against them, so this is the one list to get right.
NORMAL_COHORTS = ("KY", "SA")

# Cool hues for the normal cohorts, warm for the rest. Shape carries the same
# information as hue, so the plot survives being printed or read by someone who
# cannot separate the reds from the greens.
COHORT_STYLE = {
    "KY": ("#2A78D6", "circle",        "o"),
    "SA": ("#1BAF7A", "square",        "s"),
    "KK": ("#EB6834", "triangle-up",   "^"),
    "KR": ("#E34948", "diamond",       "D"),
    "DJ": ("#4A3AA7", "triangle-down", "v"),
    "KP": ("#008300", "cross",         "P"),
}


def cohort_of(patient_id):
    """Leading letters of the patient ID, e.g. KY153_1 -> KY."""
    import re
    return re.match(r"^[A-Za-z]+", patient_id).group(0)


def _k_at_trim(case, trim, patch_dir, method="trim"):
    """delta K over the patches a given method keeps.

    "trim"     percentile trimming on integrated K -- the quantity the reference
               cohorts were measured with, so the default and the comparable one.
    "converged" keep only patches whose local quadratic fit met CAMP's own gate
               (R2 >= 0.9, eta <= 0.075). A patch failing it was kept only because
               it hit the subdivision floor, so its <kg> is a least-squares
               accident rather than a curvature. Measured over this cohort those
               patches are 5-12% of the total and carry ~88% of delta K.

    Filters on GEOMETRY were tried and do not work: a collar around the cut rims
    (1-3 mm) discards 5% of the surface, barely moves the worst case, and makes its
    normalised value worse by stripping the normal cohort too; excluding folded
    patches (>90 deg dihedral) rescues KK15_1 but not KP35_1, whose dominant patch
    is badly fitted without being folded.
    """
    import numpy as np
    sys.path.insert(0, str(ROOT))
    import shape_metrics as SM
    d = np.load(Path(patch_dir) / f"{case}.npz")
    iK, area = d["integrated_K"], d["patch_area"]
    if method == "converged":
        keep = (d["R2_patch"] >= 0.9) & (d["eta_patch"] <= 0.075)
    else:
        keep = SM.load_campsrc()._trim_mask_integrated(iK, trim)
    q, a = iK[keep], area[keep]
    w = a / a.sum()
    return float(np.sum(w * (q - np.sum(w * q)) ** 2))


def cmd_cohort(args):
    """Area against delta K, normalised to the normal cohorts and coloured by group."""
    import numpy as np

    rows = _load(args.out)
    trim = args.trim if args.trim is not None else float(rows[0]["trim_frac"])
    patch_dir = Path(args.indir).parent / "camp_patches"
    method = getattr(args, "method", "trim") or "trim"
    for r in rows:
        r["cohort"] = cohort_of(r["case"])
        if abs(trim - float(r["trim_frac"])) > 1e-12 or method != "trim":
            r["k_fluct"] = _k_at_trim(r["case"], trim, patch_dir, method)
            r["trim_frac"] = trim
        r["method"] = method
        # Every scan here cleared remesh topology and CAMP convergence -- delta_k
        # raises otherwise, so a failed mesh never reaches this file. Kept explicit
        # so a future checkpoint has somewhere to say no.
        r["passed"] = True

    ref = [r for r in rows if r["cohort"] in NORMAL_COHORTS and r["passed"]]
    if not ref:
        raise SystemExit(f"no passing scans in the normal cohorts {NORMAL_COHORTS}")
    DK_REF = float(np.mean([r["k_fluct"] for r in ref]))
    for r in rows:
        r["k_fluct_norm"] = r["k_fluct"] / DK_REF
    print(f"Normal reference = mean delta-K of {len(ref)} passing "
          f"{'/'.join(NORMAL_COHORTS)} scans = {DK_REF:.6g}   (trim {trim:g})")
    excluded = sum(1 for r in rows
                   if r["cohort"] in NORMAL_COHORTS and not r["passed"])
    if excluded:
        print(f"({excluded} normal scan(s) excluded from the reference for "
              f"failing a checkpoint)")

    order = [c for c in COHORT_STYLE if any(r["cohort"] == c for r in rows)]
    order += sorted({r["cohort"] for r in rows} - set(order))
    print(f"\n  {'cohort':8s} {'n':>3} {'median':>8} {'mean':>8}")
    for c in order:
        v = np.array([r["k_fluct_norm"] for r in rows if r["cohort"] == c])
        print(f"  {c:8s} {len(v):3d} {np.median(v):8.3f} {v.mean():8.3f}"
              + ("   (normal)" if c in NORMAL_COHORTS else ""))

    stem = ((f"camp_nubs_100k_converged" if method == "converged" else
             f"camp_nubs_100k_trim{trim:g}".replace(".", "p"))
            + f"_norm{''.join(NORMAL_COHORTS)}_bycohort")
    csv = Path(args.out).with_name(stem + ".csv")
    cols = ["case", "cohort", "A_mm2", "k_fluct", "k_fluct_norm", "n_patches",
            "euler_number", "trim_frac", "passed"]
    csv.write_text("\n".join([",".join(cols)] + [
        ",".join(str(r.get(c, "")) for c in cols) for r in rows]) + "\n",
        encoding="utf-8")
    print(f"\n  wrote {csv.name}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not installed, no figure)")
        return 0

    fig, ax = plt.subplots(figsize=(9.6, 7.0))
    ax.axhline(1.0, ls="--", color="#898781", lw=1.2, zorder=1)
    ax.annotate(f"normal ({' + '.join(NORMAL_COHORTS)}) mean", xy=(0.012, 1.0),
                xycoords=("axes fraction", "data"), va="bottom", ha="left",
                fontsize=9, color="#6b6a66")
    for c in order:
        sub = [r for r in rows if r["cohort"] == c]
        colour, _sym, mk = COHORT_STYLE.get(c, ("#898781", "circle", "o"))
        tag = " (normal)" if c in NORMAL_COHORTS else ""
        ax.scatter([r["A_mm2"] for r in sub], [r["k_fluct_norm"] for r in sub],
                   s=95, c=colour, marker=mk, edgecolor="white", linewidth=1.0,
                   label=f"{c}{tag} (n={len(sub)})", zorder=3)
    if args.label:
        for r in rows:
            ax.annotate(r["case"], (r["A_mm2"], r["k_fluct_norm"]),
                        textcoords="offset points", xytext=(6, 4), fontsize=6.5,
                        color="#555")
    ax.set_xlabel("Surface area (mm²)")
    ax.set_ylabel("δK / normal mean")
    how = ("converged patches only (R²≥0.9, η≤0.075), no trim"
           if method == "converged" else f"trim {trim:g}")
    ax.set_title(f"GT cohort: area vs normalized δK (normals = 1.0)\n"
                 f"CAMP at {TARGET_FACES:,} faces, {how}", fontsize=11)
    ax.grid(alpha=0.25, zorder=0)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    png = Path(args.out).with_name(stem + ".png")
    fig.savefig(png, dpi=170)
    print(f"  wrote {png.name}")

    try:
        import plotly.graph_objects as go
    except ImportError:
        return 0
    f = go.Figure()
    f.add_hline(y=1.0, line_dash="dash", line_color="#898781",
                annotation_text=f"normal ({' + '.join(NORMAL_COHORTS)}) mean",
                annotation_position="top left")
    for c in order:
        sub = [r for r in rows if r["cohort"] == c]
        colour, symbol, _mk = COHORT_STYLE.get(c, ("#898781", "circle", "o"))
        tag = " (normal)" if c in NORMAL_COHORTS else ""
        f.add_trace(go.Scatter(
            x=[r["A_mm2"] for r in sub], y=[r["k_fluct_norm"] for r in sub],
            mode="markers", name=f"{c}{tag} (n={len(sub)})",
            marker=dict(size=11, color=colour, symbol=symbol,
                        line=dict(width=1, color="white")),
            hovertext=[f"{r['case']}: delta-K={r['k_fluct']:.5g} -> "
                       f"{r['k_fluct_norm']:.2f}x normal" for r in sub],
            hoverinfo="text"))
    f.update_layout(title="GT Cohort: Area vs normalized δK (normals = 1.0)",
                    xaxis_title="Surface area (mm²)",
                    yaxis_title="δK / normal mean", width=880, height=650)
    html = Path(args.out).with_name(stem + ".html")
    f.write_html(str(html))
    print(f"  wrote {html.name}")
    return 0


def _fit(A, K):
    import numpy as np
    good = np.isfinite(A) & np.isfinite(K) & (A > 0) & (K > 0)
    if good.sum() < 3:
        return float("nan"), float("nan")
    x, y = np.log10(A[good]), np.log10(K[good])
    slope, icept = np.polyfit(x, y, 1)
    pred = slope * x + icept
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(((y - pred) ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(r2)


def cmd_trims(args):
    """delta K at several trims, from the stored per-patch tables.

    delta K is a VARIANCE, so a single extreme patch can carry it. Trimming exists
    to bound that, and the reference cohort ran at 0.01 -- so a run at trim 0 is not
    wrong, but it is a different operating point and it is worth seeing by how much.
    """
    import numpy as np

    rows = _load(args.out)
    patches = Path(args.indir).parent / "camp_patches"
    cu = None
    table = {}
    for r in rows:
        npz = patches / f"{r['case']}.npz"
        if not npz.exists():
            continue
        if cu is None:
            sys.path.insert(0, str(ROOT))
            import shape_metrics as SM
            cu = SM.load_campsrc()
        d = np.load(npz)
        iK, area = d["integrated_K"], d["patch_area"]
        vals = {}
        for t in TRIMS:
            keep = cu._trim_mask_integrated(iK, t)
            q, a = iK[keep], area[keep]
            w = a / a.sum()
            vals[t] = float(np.sum(w * (q - np.sum(w * q)) ** 2))
        table[r["case"]] = (r["A_mm2"], vals, float(iK.min()), float(iK.max()))
    if not table:
        raise SystemExit(f"no per-patch tables under {patches}; re-run with --fresh")

    print(f"\ndelta K by trim, n={len(table)}")
    hdr = f"  {'case':11s} {'A_mm2':>8s}" + "".join(f"{t:>11.3f}" for t in TRIMS) \
          + f"{'min K_j':>10s}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for case, (A, vals, lo, _hi) in sorted(table.items(),
                                           key=lambda kv: -kv[1][1][0.0]):
        print(f"  {case:11s} {A:8.0f}"
              + "".join(f"{vals[t]:11.4g}" for t in TRIMS)
              + f"{lo:10.1f}")

    print(f"\n  log-log fit of delta K against area, by trim")
    print(f"  {'trim':>7s} {'slope':>8s} {'R2':>7s}")
    A = np.array([v[0] for v in table.values()], float)
    cases = list(table)
    for t in TRIMS:
        K = np.array([v[1][t] for v in table.values()], float)
        slope, r2 = _fit(A, K)
        print(f"  {t:7.3f} {slope:8.3f} {r2:7.3f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not installed, no figure written)")
        return 0

    show = [t for t in TRIMS if t in (0.0, 0.005, 0.01, 0.02)]
    fig, axes = plt.subplots(1, len(show), figsize=(5.0 * len(show), 5.0),
                             sharex=True)
    worst = max(cases, key=lambda c: table[c][1][0.0])
    for ax, t in zip(np.atleast_1d(axes), show):
        K = np.array([table[c][1][t] for c in cases], float)
        slope, r2 = _fit(A, K)
        hi = np.array([c == worst for c in cases])
        ax.scatter(A[~hi], K[~hi], s=44, c="#1f77b4", edgecolor="white",
                   linewidth=0.7, zorder=3)
        ax.scatter(A[hi], K[hi], s=95, c="#d62728", edgecolor="white",
                   linewidth=0.8, zorder=4, label=worst)
        xs = np.linspace(np.log10(A.min()), np.log10(A.max()), 50)
        ax.plot(10 ** xs, 10 ** (slope * xs + np.polyfit(
            np.log10(A), np.log10(K), 1)[1]), "-", color="#333", lw=1.3, zorder=2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"trim {t:g}   slope {slope:.2f},  $R^2$={r2:.2f}", fontsize=10)
        ax.set_xlabel("surface area  $A_T$  (mm$^2$)")
        ax.grid(alpha=0.25, zorder=0)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    np.atleast_1d(axes)[0].set_ylabel(r"$\delta K$  (area-weighted variance of $K_j$)")
    fig.suptitle(f"Thoracic nub-in cohort, n={len(cases)}, CAMP at 100,000 faces "
                 f"-- effect of trimming outlier patches", fontsize=11)
    fig.tight_layout()
    png = Path(args.out).with_name("camp_nubs_trim_sweep.png")
    fig.savefig(png, dpi=170)
    print(f"  wrote {png.name}")

    csv = Path(args.out).with_name("camp_nubs_trim_sweep.csv")
    head = ["case", "A_mm2"] + [f"k_fluct_trim{t:g}" for t in TRIMS] + ["min_K_j"]
    lines = [",".join(head)]
    for c in cases:
        A_c, vals, lo, _hi = table[c]
        lines.append(",".join([c, f"{A_c:.1f}"]
                              + [f"{vals[t]:.6g}" for t in TRIMS] + [f"{lo:.3f}"]))
    csv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  wrote {csv.name}")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "plot", "trims", "cohort"):
        p = sub.add_parser(name)
        p.add_argument("--indir", default=str(DEFAULT_IN))
        p.add_argument("--out", default=str(DEFAULT_OUT))
        p.add_argument("--workers", type=int, default=6)
        p.add_argument("--faces", type=int, default=None,
                       help="target face count for the remesh (default 100,000)")
        p.add_argument("--trim", type=float, default=None,
                       help="re-trim from the stored per-patch tables (no re-run)")
        p.add_argument("--norm-group", default=None,
                       help="report delta K as a multiple of this group's mean, "
                            "matched on the case-name prefix (e.g. KY)")
        p.add_argument("--fresh", action="store_true",
                       help="ignore existing results and recompute everything")
        p.add_argument("--repair", action="store_true",
                       help="smooth folded edges (segmentation defects) before "
                            "measuring; see repair_folds")
        p.add_argument("--method", choices=("trim", "converged"), default="trim",
                       help="'trim' percentile-trims integrated K (default); "
                            "'converged' keeps only patches whose fit met CAMP's "
                            "R2/eta gate")
        p.add_argument("--label", action="store_true",
                       help="write the case name beside each point")
        p.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    return {"run": cmd_run, "plot": cmd_plot, "trims": cmd_trims,
            "cohort": cmd_cohort}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
