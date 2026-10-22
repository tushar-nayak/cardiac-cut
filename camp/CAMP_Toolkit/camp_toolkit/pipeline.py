"""End-to-end CAMP driver for collaborators.

    STL (from vascular-cut)  ->  weld  ->  ACVD uniform remesh  ->  CAMP  ->  outputs

The ONLY preprocessing between your STL and CAMP is one uniform remeshing step
(required — CAMP's patch statistics assume a roughly uniform triangulation).
Nothing else is done to the geometry.

Public API:
    run_camp(inputs, out_dir, target_faces=..., trim_frac=..., **camp_params)
        inputs may be a single .stl file OR a folder of .stl files.
    run_stl(stl_path, out_dir, ...)   -- one mesh, returns a summary dict.
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import trimesh
from scipy.io import savemat

from . import campsrc
from .remesh import acvd_uniform_remesh
from .visualize import save_patch_html, save_patch_png

WELD_DECIMALS = 5

# CAMP quality thresholds (the operating point used throughout the project).
DEFAULT_PARAMS = dict(
    R2_thresh=0.9,       # min fraction of patch variance the quadratic must explain
    eta_thresh=0.075,    # max curvature coefficient-of-variation within a patch
    min_patch_size=18,   # vertex floor: a patch this small is never split further
    seed_patches=100,    # number of initial KMeans seed patches
)


# ── mesh preparation ─────────────────────────────────────────────────────────

def weld_mesh(mesh, decimals=WELD_DECIMALS):
    """Merge duplicated STL corners so triangles share edges.

    STL stores each facet independently, so ``trimesh.load(process=False)`` gives a
    "triangle soup" with no shared edges and an empty ``face_adjacency`` — which CAMP
    needs.  Welding only remaps vertex indices; it preserves face COUNT and ORDER
    exactly.  This is NOT remeshing.
    """
    v = np.asarray(mesh.vertices, dtype=float)
    f = np.asarray(mesh.faces, dtype=np.int64)
    key = np.round(v, decimals)
    _, idx_first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = np.asarray(inv).reshape(-1)
    return trimesh.Trimesh(vertices=v[idx_first], faces=inv[f], process=False)


def _write_temp_mat(mesh, stem, tmp_dir):
    tmp_dir = Path(tmp_dir); tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp = tmp_dir / f"{stem}.mat"
    savemat(str(tmp), {
        "mesh_surface_vertices":  np.asarray(mesh.vertices, dtype=float),
        "mesh_surface_triangles": np.asarray(mesh.faces, dtype=np.int64) + 1,  # 1-based for CAMP
    })
    return tmp


# ── shape metric ─────────────────────────────────────────────────────────────

def _k_fluct(integrated_K, patch_area, keep):
    """Area-weighted variance of integrated Gaussian curvature over the kept patches.

    This is CAMP's headline shape descriptor (higher = more heterogeneous surface
    shape).  Matches the definition used across the project.
    """
    q = np.asarray(integrated_K, float)[keep]
    a = np.asarray(patch_area, float)[keep]
    valid = np.isfinite(q) & np.isfinite(a) & (a > 0)
    q, a = q[valid], a[valid]
    at = float(a.sum())
    if at <= 0 or q.size == 0:
        return np.nan, at
    w = a / at
    qbar = float(np.sum(w * q))
    return float(np.sum(w * (q - qbar) ** 2)), at


# ── one mesh ─────────────────────────────────────────────────────────────────

def run_stl(stl_path, out_dir, target_faces=100_000, trim_frac=0.01,
            params=None, save_html=True, save_png=False, verbose=True,
            campsrc_module=None):
    """Run the full pipeline on ONE stl and write its outputs under ``out_dir/<stem>/``.

    Returns a one-row summary dict (also collected by :func:`run_camp`).

    Outputs written per mesh:
        <stem>_facepatches.csv   face_id, patch_id           (on the remeshed mesh)
        <stem>_patch_metrics.csv patch_id, avg_kg, area_mm2, integrated_K, eta, R2, kept
        <stem>_camp.npz          V, F, face_labels + per-patch arrays (for re-plotting)
        <stem>_patches.html      interactive viewer            (if save_html)
        <stem>_integratedK.png   static render                 (if save_png)
    """
    params = {**DEFAULT_PARAMS, **(params or {})}
    camp_source = campsrc_module or campsrc
    stl_path = Path(stl_path)
    stem = stl_path.stem
    mesh_out = Path(out_dir) / stem
    mesh_out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    row = {"stem": stem, "status": "ok", "error": ""}

    try:
        # 1) load raw STL (no auto-processing) and weld duplicated corners
        raw = trimesh.load(str(stl_path), process=False)
        welded = weld_mesh(raw)
        row["faces_in"] = int(len(welded.faces))
        row["euler_in"] = int(welded.euler_number)
        if verbose:
            print(f"  [{stem}] welded: {len(welded.vertices)} v, {len(welded.faces)} f, "
                  f"chi={row['euler_in']}")

        # 2) uniform remesh to the target density  (the ONE preprocessing step)
        if target_faces and target_faces > 0:
            rm = acvd_uniform_remesh(welded, target_faces=target_faces)
        else:
            rm = welded                                        # advanced: skip remeshing
        row["faces_remesh"] = int(len(rm.faces))
        row["euler_remesh"] = int(rm.euler_number)
        if verbose:
            print(f"  [{stem}] remeshed -> {len(rm.faces)} faces (target {target_faces})")

        # 3) CAMP  (run once at trim 0 -> full partition; trim is applied to the metric below)
        tmp = _write_temp_mat(rm, stem, mesh_out / "_tmp")
        try:
            out = camp_source.camp(mat_path=str(tmp), scan_prefix=None, make_figure=False,
                                   trim_frac=0.0, **params)
        finally:
            try: tmp.unlink()
            except OSError: pass
            try: (mesh_out / "_tmp").rmdir()
            except OSError: pass

        mesh = out["mesh"]
        V = np.asarray(mesh.vertices, float)
        F = np.asarray(mesh.faces, np.int64)
        face_labels = np.asarray(out["face_labels"], np.int64)
        avg_kg   = np.asarray(out["avg_kg"], float)
        area     = np.asarray(out["patch_area"], float)
        integ_K  = np.asarray(out["K_integrated_patch"], float)
        eta      = np.asarray(out["eta_patch"], float)
        R2       = np.asarray(out["R2_patch"], float)

        # 4) apply the requested trim to the SHAPE METRIC (partition/plots keep all patches)
        keep = camp_source._trim_mask_integrated(integ_K, trim_frac)
        if not np.any(keep):
            keep = np.ones_like(integ_K, dtype=bool)
        kf, area_total = _k_fluct(integ_K, area, keep)

        n_patches = int(len(integ_K))
        row.update({
            "surface_area_mm2": float(mesh.area),
            "n_faces":          int(len(F)),
            "n_patches":        n_patches,
            "n_kept":           int(keep.sum()),
            "n_trimmed":        int(n_patches - keep.sum()),
            "k_fluct":          kf,
            "target_faces":     int(target_faces) if target_faces else 0,
            "trim_frac":        float(trim_frac),
            "pct_at_min_size":     float(out["pct_at_min_size"]),
            "pct_failed_quality":  float(out["pct_failed_quality"]),
            "runtime_s":        round(time.time() - t0, 1),
        })

        # 5) write per-mesh files
        pd.DataFrame({"face_id": np.arange(len(F)), "patch_id": face_labels}) \
          .to_csv(mesh_out / f"{stem}_facepatches.csv", index=False)

        pd.DataFrame({
            "patch_id":     np.arange(n_patches),
            "avg_kg":       avg_kg,
            "area_mm2":     area,
            "integrated_K": integ_K,
            "eta":          eta,
            "R2":           R2,
            "kept":         keep,
        }).to_csv(mesh_out / f"{stem}_patch_metrics.csv", index=False)

        np.savez_compressed(
            mesh_out / f"{stem}_camp.npz",
            V=V, F=F, face_labels=face_labels.astype(np.int32),
            avg_kg=avg_kg, patch_area=area, K_integrated_patch=integ_K,
            eta_patch=eta, R2_patch=R2, keep=keep)

        if save_html:
            p = save_patch_html(mesh_out / f"{stem}_patches.html", V, F, face_labels, integ_K,
                                title=f"{stem}  ·  {n_patches} patches  ·  k_fluct={kf:.4g}")
            row["html"] = p
        if save_png:
            row["png"] = save_patch_png(mesh_out / f"{stem}_integratedK.png",
                                        V, F, face_labels, integ_K)

        if verbose:
            print(f"  [{stem}] DONE: {n_patches} patches, k_fluct={kf:.4g}, "
                  f"{row['runtime_s']}s -> {mesh_out}")
        return row

    except Exception as e:
        import traceback
        traceback.print_exc()
        row.update({"status": "error", "error": str(e),
                    "runtime_s": round(time.time() - t0, 1)})
        return row


# ── file or folder ───────────────────────────────────────────────────────────

def _collect_stls(inputs, pattern="*.stl", recursive=False):
    p = Path(inputs)
    if p.is_file():
        return [p]
    if p.is_dir():
        globber = p.rglob if recursive else p.glob
        return sorted(q for q in globber(pattern)
                      if q.is_file() and not q.name.startswith("._"))
    raise FileNotFoundError(f"No such file or folder: {inputs}")


def run_camp(inputs, out_dir="camp_outputs", target_faces=100_000, trim_frac=0.01,
             params=None, save_html=True, save_png=False,
             pattern="*.stl", recursive=False, verbose=True):
    """Run CAMP on a single STL or every STL in a folder.

    Parameters
    ----------
    inputs : str | Path
        A ``.stl`` file, or a folder containing ``.stl`` files.
    out_dir : str | Path
        Where results are written (one subfolder per mesh + a top-level summary.csv).
    target_faces : int
        Uniform remesh density before CAMP.  Project defaults: ~50k for plain
        thoracic/abdominal aortas, ~100k when branches are present.  Set 0 to skip
        remeshing (advanced; only if your mesh is already uniform).
    trim_frac : float
        Symmetric quantile trim on integrated-K when computing the k_fluct metric
        (0.01 = production).  Does not change the saved partition or the plots.
    params : dict, optional
        Overrides for CAMP thresholds: R2_thresh, eta_thresh, min_patch_size,
        seed_patches.  See ``DEFAULT_PARAMS``.
    save_html / save_png : bool
        Which visualizations to emit per mesh.

    Returns
    -------
    pandas.DataFrame
        One row per mesh (also saved to ``out_dir/summary.csv``).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stls = _collect_stls(inputs, pattern=pattern, recursive=recursive)
    if not stls:
        raise FileNotFoundError(f"No STL files matched '{pattern}' under {inputs}")

    print(f"[CAMP] {len(stls)} mesh(es) | remesh->{target_faces} faces | "
          f"trim={trim_frac} | out={out_dir}")
    rows = []
    for n, stl in enumerate(stls, 1):
        print(f"\n[{n}/{len(stls)}] {stl.name}")
        rows.append(run_stl(stl, out_dir, target_faces=target_faces, trim_frac=trim_frac,
                            params=params, save_html=save_html, save_png=save_png,
                            verbose=verbose))

    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "summary.csv", index=False)
    ok = int((summary["status"] == "ok").sum())
    print(f"\n[CAMP] finished: {ok}/{len(summary)} ok  ->  {out_dir / 'summary.csv'}")
    return summary
