# delta K (k_fluct) for an arbitrary surface, via the lab's CAMP partition.
#
# CAMP is loaded BY PATH from the project's own camp/ folder (cardiac_cut/camp), the
# canonical campsrc.py. By path, not by name: three files named campsrc.py existed in
# the source project and they are NOT the same method; whichever reached sys.path
# first would answer, and a wrong one returns a different number without complaining.
# This folder's own camp/ package is the same code (line endings aside) and is no
# longer what answers.
#
# CAMP_SRC can name another campsrc.py, but only one under the project's camp/.
#
# Nothing in here needs the cut tool. It takes any triangular surface, which is what
# lets the website offer shape metrics as its own section for people who already have
# geometries.

import importlib.util
import os
import tempfile
from pathlib import Path

import numpy as np
import trimesh
from scipy.io import savemat

# The project's camp/ folder, beside this one, and the canonical build in it.
PROJECT_CAMP = Path(__file__).resolve().parent.parent / "camp"
CANONICAL_CAMP = PROJECT_CAMP / "camp_updated_2026-08-26" / "campsrc.py"
CAMP_SRC = os.environ.get("CAMP_SRC") or str(CANONICAL_CAMP)

# The reference cohort's operating point, measured on the 30 thoracic nub geometries
# in tools/camp_cohort_final_trim0p0025.csv. delta K is extensive and partition
# dependent, so cut extent, face count, trim and topology all move it -- an upload is
# only comparable to the reference when it matches on every one of them.
#
# trim_frac is 0.0025 for a measured reason rather than a stylistic one, and it is the
# single most important thing to understand about this number. Remeshing the SAME
# surface to 98k, 100k and 102k faces changes nothing about the anatomy, but it
# rebuilds the partition from scratch -- and at trim 0 the answer then moves by a
# median factor of 2.6, and by 23.5 on the worst case. See CURVATURE.md.
REFERENCE = {
    "n_cases": 30,
    "target_faces": 100_000,
    "trim_frac": 0.0025,
    "area_mm2": (15_921.0, 87_190.0),
    "pct_at_min_size": (9.8, 29.1),
}

# How far apart two delta K values on the SAME surface land when the only thing that
# differs is where the remesher put its vertices: trim -> (median, worst) ratio over
# 30 cases x 98k/100k/102k. Produced by tools/camp_noise.py, data in
# tools/camp_noise.jsonl. A difference smaller than the median is not a finding.
REPEATABILITY = {0.0: (2.64, 23.48), 0.001: (1.88, 10.88), 0.0025: (1.59, 3.81),
                 0.005: (1.44, 2.74), 0.01: (1.28, 1.78)}

PARAMS = dict(R2_thresh=0.9, eta_thresh=0.075, min_patch_size=18, seed_patches=100)

# PARTITION NOISE, measured 2026-08-20 on 30 aortas x 8 k-means seeds x 3 densities.
#
# delta K is carried by ~22 patches out of ~850, sited where curvature changes fastest
# -- the branch ostia on a branched surface. Where k-means puts a boundary there
# decides the number, so ONE partition is one draw from a distribution, not the answer.
# Re-running the same mesh with a different seed moves delta K by ~16%.
#
# The draws are independent: over 18 case-density cells, the spread of a 4-seed mean
# matched CV/sqrt(4) with a median ratio of 1.03 (range 1.01-1.13). So averaging works
# exactly as the arithmetic says, and n_seeds buys precision at sqrt(n) cost.
#
# This SUPERSEDES the REPEATABILITY table above for anything run with the connectivity
# split and fan gate on. That table was measured before both, mixes partition noise
# with defects since fixed, and is far more pessimistic than the pipeline now is.
PARTITION_CV = 0.162          # one partition, at a fixed mesh
PARTITION_CV_RANGE = (0.07, 0.30)

# What averaging does NOT fix: delta K drifts systematically with mesh density, because
# min_patch_size is a VERTEX count, so the smallest patch covers 16 mm2 at 100k faces
# and 4 mm2 at 400k. Those are measurements at different scales, not repeat draws, and
# the mean of 8 seeds at each density is as far apart as a single seed was (CV 0.22 vs
# 0.21). Numbers are comparable only at a fixed target_faces.
DENSITY_CV = 0.21


def partition_band(n_seeds=1):
    """Expected relative uncertainty on delta K averaged over *n_seeds* partitions."""
    n = max(int(n_seeds), 1)
    return PARTITION_CV / (n ** 0.5)

# A dihedral past this is the surface doubling back through itself -- a segmentation
# defect, not anatomy. Reported rather than repaired: see surface_defects.
FOLD_DEG = 90.0
SLIVER_Q = 0.05           # 1.0 is equilateral; below this a face normal is noise


def repeatability_band(trim_frac=None):
    """(median, worst) ratio two runs of the same surface can differ by, at *trim_frac*.

    Interpolation between the measured trims is not offered: the numbers come from a
    sweep at five settings and the nearest measured one is quoted, because a smooth
    curve through five points would imply a precision the sweep does not have.
    """
    t = float(REFERENCE["trim_frac"] if trim_frac is None else trim_frac)
    key = min(REPEATABILITY, key=lambda k: abs(k - t))
    return REPEATABILITY[key], key

_campsrc = None


class ShapeError(RuntimeError):
    pass


def load_campsrc():
    """The CAMP module, loaded once, by path, from the project's camp/ folder.

    CAMP_SRC picks a different campsrc.py, but it must still be under camp/: a CAMP
    from anywhere else is refused rather than quietly used.
    """
    global _campsrc
    if _campsrc is not None:
        return _campsrc

    p = Path(CAMP_SRC).resolve()
    if not p.is_relative_to(PROJECT_CAMP.resolve()):
        raise ShapeError(f"CAMP must come from {PROJECT_CAMP}; CAMP_SRC points at {p}.")
    if not p.exists():
        raise ShapeError(f"CAMP source {p} does not exist.")
    spec = importlib.util.spec_from_file_location("campsrc_project", str(p))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _campsrc = mod
    return _campsrc


def campsrc_origin():
    """Where CAMP is being loaded from -- worth showing in a UI, since the number
    depends on which build produced it."""
    p = Path(CAMP_SRC).resolve()
    try:
        return str(p.relative_to(PROJECT_CAMP.resolve().parent))
    except ValueError:
        return str(p)


def weld(mesh, decimals=5):
    """Merge duplicated STL corners into a shared-vertex mesh, keeping face count
    and order. STL stores every facet independently, and without this CAMP's
    connectivity repair quietly becomes a no-op -- you get a number, and it is wrong.
    """
    v = np.asarray(mesh.vertices, dtype=float)
    f = np.asarray(mesh.faces, dtype=np.int64)
    _, first, inv = np.unique(np.round(v, decimals), axis=0,
                              return_index=True, return_inverse=True)
    return trimesh.Trimesh(vertices=v[first],
                           faces=np.asarray(inv).reshape(-1)[f], process=False)


def surface_defects(mesh, fold_deg=FOLD_DEG):
    """What the segmentation left behind, before any of it reaches CAMP.

    Folds -- neighbouring faces meeting at more than *fold_deg* -- are the headline.
    They are geometry the vessel does not have, they arrive in the source
    segmentation rather than being introduced by remeshing, and 14 of the 30
    reference cases carry them.

    What they are NOT is an explanation for delta K being partition sensitive. Across
    those 30 cases fold count and the spread over 98k/100k/102k correlate at -0.07,
    and the single worst case for spread (23.5x) has no folded edges at all. So this
    is reported as a property of the input, and never as a correction to the number.

    *fold_span_mm* is how far the folded faces sit from their own centroid: a
    segmentation defect is local and gives a small span, whereas noise everywhere
    gives a large one.
    """
    w = weld(mesh)
    ang = np.degrees(np.asarray(w.face_adjacency_angles))
    fa = np.asarray(w.face_adjacency)
    folded = ang > fold_deg

    tri = w.triangles
    e = np.linalg.norm(np.diff(np.concatenate([tri, tri[:, :1]], axis=1),
                               axis=1), axis=2)
    q = 4 * np.sqrt(3) * w.area_faces / np.maximum(np.sum(e ** 2, axis=1), 1e-12)

    uq, cnt = np.unique(w.edges_sorted.reshape(-1, 2), axis=0, return_counts=True)

    span = None
    if folded.any():
        c = np.asarray(w.triangles_center)[np.unique(fa[folded].ravel())]
        span = float(np.linalg.norm(c - c.mean(axis=0), axis=1).max())

    return {
        "faces": int(len(w.faces)),
        "folds": int(folded.sum()),
        "max_dihedral": float(ang.max()) if len(ang) else 0.0,
        "fold_span_mm": span,
        "slivers": int((q < SLIVER_Q).sum()),
        "min_quality": float(q.min()),
        "nonmanifold_edges": int((cnt > 2).sum()),
        "boundary_edges": int((cnt == 1).sum()),
        "euler": int(w.euler_number),
        "winding_ok": bool(w.is_winding_consistent),
        "area_mm2": float(w.area),
    }


def delta_k(mesh, target_faces=None, trim_frac=None, make_figure=True,
            post_remesh_repair=None,
            expect_euler=None):
    """Partition a surface and return its curvature fluctuation.

    *expect_euler* is the topology before remeshing. ACVD may move vertices but must
    not open or close a hole; passing the expected value turns that into an assertion
    rather than something nobody checks.
    """
    from cutclean import acvd_uniform_remesh

    cu = load_campsrc()
    target_faces = int(target_faces or REFERENCE["target_faces"])
    trim_frac = float(REFERENCE["trim_frac"] if trim_frac is None else trim_frac)

    if isinstance(mesh, (str, Path)):
        mesh = trimesh.load(str(mesh), process=False)
    welded = weld(mesh)
    chi_in = int(welded.euler_number)
    if expect_euler is not None and chi_in != int(expect_euler):
        raise ShapeError(f"the uploaded surface has euler number {chi_in}, "
                         f"but {int(expect_euler)} was expected")

    remeshed = acvd_uniform_remesh(welded, target_faces=target_faces)

    # REPAIR AFTER THE REMESH, NOT ONLY BEFORE IT.
    #
    # Callers repair folds on the surface they hand in, but the remesh happens here,
    # and ACVD reintroduces folds that the repaired input did not have: KK28_11 goes
    # from 0 folded edges to 2 at 100k and 6 at 200k, KK20_2 from 0 to 7 and 4. The
    # count swings with density, so it is not a fixed property of the geometry -- it
    # is noise that changes every time the mesh is rebuilt. A quadratic fitted across
    # a fold returns a curvature that is not a curvature, which is how a single patch
    # comes to carry half a case's delta K.
    #
    # The repair is supplied by the caller rather than imported, because it lives in
    # tools/camp_cohort.py which imports this module.
    if post_remesh_repair is not None:
        remeshed, repair_info = post_remesh_repair(remeshed)

    chi_out = int(remeshed.euler_number)
    if chi_out != chi_in:
        raise ShapeError(
            f"remeshing changed the topology (euler {chi_in} -> {chi_out}): an "
            f"opening was created or closed, so the result would not be comparable")

    tmp = Path(tempfile.mkdtemp()) / "surface.mat"
    savemat(str(tmp), {
        "mesh_surface_vertices": np.asarray(remeshed.vertices, dtype=float),
        "mesh_surface_triangles": np.asarray(remeshed.faces, dtype=np.int64) + 1,
    })
    try:
        # trim_frac=0.0 deliberately. camp() returns a FULL face_labels array
        # alongside TRIMMED per-patch arrays, so letting it trim internally leaves
        # the two misaligned and every patch colour wrong. Trim afterwards instead;
        # the partition does not depend on it.
        out = cu.camp(mat_path=str(tmp), scan_prefix=None, trim_frac=0.0,
                      color_metric="integrated_K", color_scale="Plasma",
                      make_figure=make_figure, **PARAMS)
    finally:
        tmp.unlink(missing_ok=True)
        try:
            tmp.parent.rmdir()
        except OSError:
            pass

    iK = np.asarray(out["K_integrated_patch"], dtype=float)
    area = np.asarray(out["patch_area"], dtype=float)
    keep = cu._trim_mask_integrated(iK, trim_frac)
    q, a = iK[keep], area[keep]
    w = a / a.sum()
    k_fluct = float(np.sum(w * (q - np.sum(w * q)) ** 2))

    A = float(out["mesh"].area)
    lo, hi = REFERENCE["area_mm2"]
    return {
        "k_fluct": k_fluct,
        "A_mm2": A,
        "n_patches": int(len(iK)),
        "n_kept": int(keep.sum()),
        "pct_at_min_size": float(out["pct_at_min_size"]),
        "pct_failed_quality": float(out["pct_failed_quality"]),
        "euler_number": chi_out,
        "faces": int(len(remeshed.faces)),
        "target_faces": target_faces,
        "trim_frac": trim_frac,
        "density_verts_per_mm2": len(remeshed.vertices) / A if A > 0 else None,
        "in_reference_range": bool(lo <= A <= hi),
        # per-patch table, for colouring patches back onto the mesh
        "face_labels": np.asarray(out["face_labels"]),
        "integrated_K": iK,
        "patch_area": area,
        "avg_kg": np.asarray(out["avg_kg"]),
        "eta_patch": np.asarray(out["eta_patch"]),
        "R2_patch": np.asarray(out["R2_patch"]),
        "keep": keep,
        "mesh": out["mesh"],
        "figure": out["figure"],
        # The leaf patches themselves, each carrying the VERTEX indices its Monge
        # fit was taken over. Needed to re-derive a patch's curvature from a subset
        # of its vertices -- which is how a fold can be excluded from the fit
        # without deleting the patch, trimming a percentile, or moving geometry.
        "patches": out["patches"],
    }


def delta_k_averaged(mesh, target_faces=None, trim_frac=None, seeds=(42,),
                     progress=None, **kw):
    """delta K averaged over several independent k-means partitions of one mesh.

    One partition is one draw. Re-seeding the same mesh moves delta K by ~16% because
    the number is carried by a handful of patches sited where curvature changes fastest,
    and where k-means puts a boundary there is arbitrary. The draws are independent, so
    the mean of n has an error of PARTITION_CV/sqrt(n) -- verified against disjoint
    halves, median ratio to theory 1.03.

    The ARRAYS returned are those of the FIRST seed, not an average: patch labels,
    per-patch curvature and the patch list belong to one concrete partition and cannot
    be meaningfully averaged across partitions that carve the surface differently. So
    the picture on screen is one partition; the NUMBER is the mean of them all. Anything
    drawn from the arrays describes seed[0], and k_fluct_seed is kept alongside so the
    two are never confused.

    seeds=(42,) reproduces the single-partition behaviour exactly, byte for byte.
    """
    seeds = tuple(seeds) or (42,)
    vals, first = [], None
    for i, sd in enumerate(seeds):
        if progress is not None:
            progress(i, len(seeds))
        orig = PARAMS.copy()
        PARAMS["seed"] = int(sd)
        try:
            r = delta_k(mesh, target_faces=target_faces, trim_frac=trim_frac,
                        make_figure=(kw.pop("make_figure", True) if i == 0 else False),
                        **kw)
        finally:
            PARAMS.clear()
            PARAMS.update(orig)
        vals.append(float(r["k_fluct"]))
        if first is None:
            first = r
    v = np.asarray(vals, dtype=float)
    out = dict(first)
    out["k_fluct"] = float(v.mean())
    out["k_fluct_seed"] = float(v[0])
    out["k_fluct_all"] = [float(x) for x in v]
    out["k_fluct_sd"] = float(v.std(ddof=1)) if len(v) > 1 else 0.0
    out["k_fluct_se"] = float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0
    out["k_fluct_cv"] = float(v.std(ddof=1) / v.mean()) if len(v) > 1 else 0.0
    out["n_seeds"] = len(seeds)
    out["seeds"] = [int(x) for x in seeds]
    return out


def patch_table(result):
    """The per-patch numbers as rows, for download."""
    import pandas as pd
    return pd.DataFrame({
        "patch_id": np.arange(len(result["integrated_K"])),
        "area_mm2": result["patch_area"],
        "avg_kg_per_mm2": result["avg_kg"],
        "integrated_K": result["integrated_K"],
        "eta": result["eta_patch"],
        "R2": result["R2_patch"],
        "kept_after_trim": result["keep"],
    })
