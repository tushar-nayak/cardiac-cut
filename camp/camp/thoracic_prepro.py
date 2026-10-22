# Thoracic aorta preprocessing pipeline for centerline generation and mesh cleaning

from pathlib import Path
import numpy as np
import pandas as pd
import SimpleITK as sitk
import trimesh
from scipy.io import savemat, loadmat

from cl import *
from cutclean import *

def load_surface_from_mat(mat_path, scan_prefix=None):
    """Load trimesh surface from MAT containing *_surface_vertices and *_surface_triangles."""
    data = loadmat(mat_path)
    keys = [k for k in data.keys() if not k.startswith("__")]

    if scan_prefix is not None:
        v_key = f"{scan_prefix}_surface_vertices"
        f_key = f"{scan_prefix}_surface_triangles"
        if v_key not in data or f_key not in data:
            raise KeyError(f"Could not find '{v_key}' and '{f_key}' in {mat_path}")
    else:
        v_candidates = [k for k in keys if k.endswith("_surface_vertices")]
        f_candidates = [k for k in keys if k.endswith("_surface_triangles")]
        if len(v_candidates) != 1 or len(f_candidates) != 1:
            raise KeyError(
                "MAT file must contain exactly one vertices/triangles pair, "
                "or provide scan_prefix explicitly."
            )
        v_key, f_key = v_candidates[0], f_candidates[0]

    verts = np.asarray(data[v_key], dtype=float)
    faces = np.asarray(data[f_key], dtype=int)

    if faces.min() == 1:
        faces = faces - 1

    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    mesh.remove_unreferenced_vertices()
    return mesh, (v_key, f_key)

def image_volume_mm3(mhd_img):
    vol = sitk.GetArrayFromImage(mhd_img)
    sx, sy, sz = mhd_img.GetSpacing()
    voxel_vol_mm3 = sx * sy * sz
    n_in = np.count_nonzero(vol > 0)
    return float(n_in * voxel_vol_mm3)

def preprocess(mhd_img, iso_spacing_mm=3.0):
    orig_spacing = np.array(mhd_img.GetSpacing(), float)
    orig_size = np.array(mhd_img.GetSize(), int)
    new_spacing = (float(iso_spacing_mm), float(iso_spacing_mm), float(iso_spacing_mm))
    new_size = [int(np.round(sz * sp / iso_spacing_mm)) for sz, sp in zip(orig_size, orig_spacing)]
    resampler = sitk.ResampleImageFilter()
    resampler.SetOutputSpacing(new_spacing)
    resampler.SetSize(new_size)
    resampler.SetOutputOrigin(mhd_img.GetOrigin())
    resampler.SetOutputDirection(mhd_img.GetDirection())
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    mhd_iso = resampler.Execute(mhd_img)
    vol_mm3 = image_volume_mm3(mhd_img)
    vol_mm3_iso = image_volume_mm3(mhd_iso)
    
    info = {
        "iso_spacing_mm": float(iso_spacing_mm),
        "mhd_volume_mm3": float(vol_mm3),
        "mhd_resampled_volume_mm3": float(vol_mm3_iso),
    }

    return mhd_iso, info

def _foreground_native_spread_xyz(mhd_img):
    """Return foreground spread along native x/y/z axes in mm.

    This intentionally ignores direction cosines and uses index-space extents
    with spacing, because downstream centerline/cap logic assumes a dominant
    coordinate axis frame where z is the long axis.
    """
    vol = sitk.GetArrayFromImage(mhd_img)
    fg = np.argwhere(vol > 0)  # (k, j, i)
    if fg.size == 0:
        return np.zeros(3, dtype=float)

    k_min, j_min, i_min = fg.min(axis=0)
    k_max, j_max, i_max = fg.max(axis=0)

    sx, sy, sz = [float(v) for v in mhd_img.GetSpacing()]
    dx = (int(i_max) - int(i_min) + 1) * sx
    dy = (int(j_max) - int(j_min) + 1) * sy
    dz = (int(k_max) - int(k_min) + 1) * sz
    return np.array([dx, dy, dz], dtype=float)

def orient_mhd_long_axis_to_z(mhd_img):
    """Permute image axes so the largest foreground spread maps to native z."""
    spread_xyz = _foreground_native_spread_xyz(mhd_img)
    largest_axis = int(np.argmax(spread_xyz))

    if largest_axis == 2:
        return mhd_img, {
            "changed": False,
            "order": (0, 1, 2),
            "inverse_order": (0, 1, 2),
            "spread_xyz_mm": spread_xyz,
            "largest_axis_before": largest_axis,
        }

    order = tuple(ax for ax in (0, 1, 2) if ax != largest_axis) + (largest_axis,)
    inv_order = tuple(int(x) for x in np.argsort(order))
    reoriented = sitk.PermuteAxes(mhd_img, list(order))

    return reoriented, {
        "changed": True,
        "order": order,
        "inverse_order": inv_order,
        "spread_xyz_mm": spread_xyz,
        "largest_axis_before": largest_axis,
    }

def permute_mesh_axes(mesh, order):
    """Return a mesh copy with vertex coordinate axes permuted by order."""
    out = mesh.copy()
    out.vertices = np.asarray(out.vertices)[:, list(order)]
    return out

def permute_centerline_axes(centerline_pts, centerline_tangents, order):
    """Permute centerline point/tangent coordinate axes by order."""
    pts = np.asarray(centerline_pts, dtype=float)[:, list(order)]
    tang = np.asarray(centerline_tangents, dtype=float)[:, list(order)]
    return pts, tang

def permute_endpoint_hint_axes(endpoint_hint, order):
    """Permute stored crude endpoint coordinates by the same mesh/image axis order."""
    if endpoint_hint is None:
        return None
    out = dict(endpoint_hint)
    if "crude_start_world" in out:
        out["crude_start_world"] = np.asarray(out["crude_start_world"], dtype=float)[list(order)]
        out["crude_start_z"] = float(out["crude_start_world"][2])
    if "crude_end_world" in out:
        out["crude_end_world"] = np.asarray(out["crude_end_world"], dtype=float)[list(order)]
        out["crude_end_z"] = float(out["crude_end_world"][2])
    return out

def flip_image_native_z(mhd_img):
    """Flip a binary image along native z (index k) while keeping metadata."""
    arr = sitk.GetArrayFromImage(mhd_img)
    arr_flipped = arr[::-1, :, :]
    out = sitk.GetImageFromArray(arr_flipped)
    out.SetSpacing(mhd_img.GetSpacing())
    out.SetOrigin(mhd_img.GetOrigin())
    out.SetDirection(mhd_img.GetDirection())
    return out

def flip_mesh_native_z(mesh, mhd_img):
    """Mirror mesh coordinates along native z about the image z-center."""
    out = mesh.copy()
    verts = np.asarray(out.vertices, dtype=float).copy()

    _, _, sz = [float(v) for v in mhd_img.GetSpacing()]
    _, _, oz = [float(v) for v in mhd_img.GetOrigin()]
    _, _, nz = [int(v) for v in mhd_img.GetSize()]
    z_center = oz + 0.5 * float(nz - 1) * sz

    verts[:, 2] = 2.0 * z_center - verts[:, 2]
    out.vertices = verts
    return out

def thoracic_preprocessing_batch(
    cohorts=["GE"],
    base_root=r"Z:\aorta\thoracic",
    iso_spacing_mm=3.0,
    buffer_mm=5.0,
    target_faces=50_000,
    individ_scans=None,
    include_centerline_in_log=False,
    continue_on_error=True,
):
    """
    For each scan that has both an MHD file and an existing smoothed and meshed MAT mesh:
      1. preprocess MHD  -> mhd_iso_img + info
      2. load_surface_from_mat   -> existing mesh
      3. centerline_for_mesh -> centerline in the mesh frame
      4. cut_clean (using existing MAT mesh + mesh-frame centerline to open mesh and remesh with ACVD)
      5. savemat (writes a *_o<k> MAT alongside the input)

    Expected layout:
      base_root/cohort/patient/scan/mhd/<scan_id>.mhd
      base_root/cohort/patient/scan/mesh/<scan_id>.mat

    Returns a dataframe with per-scan status.

        individ_scans:
            optional iterable of scan_id strings. If provided, only these
            scans are queued/processed (across the provided cohort folders).

        include_centerline_in_log:
            if True, adds object columns with final centerline points/tangents
            used for cutting as nested python lists per scan.
    """

    base_root = Path(base_root)
    jobs = []
    scan_filter = None
    if individ_scans is not None:
        scan_filter = {str(s).strip() for s in individ_scans if str(s).strip()}
        if len(scan_filter) == 0:
            print("[warn] individ_scans was provided but empty after cleaning; nothing will be queued.")

    for cohort in cohorts:
        cohort_dir = base_root / cohort
        if not cohort_dir.exists():
            print(f"[warn] Cohort folder not found: {cohort_dir}")
            continue
        for patient_dir in sorted([p for p in cohort_dir.iterdir() if p.is_dir()]):
            for scan_dir in sorted([s for s in patient_dir.iterdir() if s.is_dir()]):
                scan_id = scan_dir.name
                if scan_filter is not None and scan_id not in scan_filter:
                    continue
                mhd_dir = scan_dir / "mhd"
                mesh_dir = scan_dir / "mesh"

                # MHD file
                preferred_mhd = mhd_dir / f"{scan_id}.mhd"
                if preferred_mhd.exists():
                    mhd_file = preferred_mhd
                elif mhd_dir.exists():
                    candidates = sorted(mhd_dir.glob("*.mhd"))
                    if not candidates:
                        continue
                    mhd_file = candidates[0]
                else:
                    continue

                # MAT file
                preferred_mat = mesh_dir / f"{scan_id}.mat"
                if preferred_mat.exists():
                    mat_file = preferred_mat
                elif mesh_dir.exists():
                    mats = sorted(mesh_dir.glob("*.mat"))
                    if not mats:
                        continue
                    mat_file = mats[0]
                else:
                    continue

                jobs.append({
                    "cohort":     cohort,
                    "patient_id": patient_dir.name,
                    "scan_id":    scan_id,
                    "scan_dir":   scan_dir,
                    "mhd_file":   mhd_file,
                    "mat_file":   mat_file,
                })

    if scan_filter is not None:
        found_ids = {j["scan_id"] for j in jobs}
        missing = sorted(scan_filter - found_ids)
        if missing:
            print(f"[warn] individ_scans not found with required MHD+MAT under cohorts {cohorts}: {missing}")

    print(f"Found {len(jobs)} scans with MHD + MAT across cohorts: {cohorts}")

    rows = []

    for i, job in enumerate(jobs, start=1):
        scan_id = job["scan_id"]
        mhd_file = Path(job["mhd_file"])
        mat_file = Path(job["mat_file"])

        print(f"[{i}/{len(jobs)}] Preprocessing {scan_id}")
        out_mat_file = mat_file
        resid = np.nan

        try:
            # Step 1: Preprocess MHD -> isotropic image + info
            mhd_iso_img, _ = preprocess(mhd_img=sitk.ReadImage(str(mhd_file)),
                                        iso_spacing_mm=iso_spacing_mm,
                                        )

            # Keep centerline + cap removal stable by forcing the dominant extent to z.
            mhd_proc_img, orient_info = orient_mhd_long_axis_to_z(mhd_iso_img)
            if orient_info["changed"]:
                axis_name = {0: "x", 1: "y", 2: "z"}[orient_info["largest_axis_before"]]
                spread = np.round(orient_info["spread_xyz_mm"], 2).tolist()
                print(
                    f"[orientation] Reoriented MHD before centerline generation and cap removal "
                    f"(largest spread was along {axis_name}; spread xyz mm={spread}; "
                    f"axis order={orient_info['order']})."
                )

            # Step 2: Load existing MAT mesh
            mesh_in, _ = load_surface_from_mat(str(mat_file), scan_prefix=None)
            mesh_proc = permute_mesh_axes(mesh_in, orient_info["order"])

            # Step 3: Generate centerline
            pos, tang, T, resid, endpoint_hint = generate_cl(
                mhd_proc_img,
                n_points=1200,
                endpoint_axis_mode="auto",
                return_endpoint_hint=True,
            )
            if orient_info["changed"]:
                pos, tang = permute_centerline_axes(pos, tang, orient_info["order"])
                endpoint_hint = permute_endpoint_hint_axes(endpoint_hint, orient_info["order"])
            final_pos = pos
            final_tang = tang
                                                                 
            # Step 4: Cut and clean the mesh using the centerline and target remeshing density
            remeshed = cut_clean(mesh_proc, 
                                 pos,
                                 buffer_mm=buffer_mm, 
                                 mhd_iso_img=mhd_proc_img, 
                                 target_faces=target_faces,
                                 endpoint_hint=endpoint_hint,
                                 )

            if orient_info["changed"] and int(remeshed.euler_number) != 0:
                print(
                    "[orientation-retry] Reoriented case still NOT OPEN; "
                    "flipping native z and retrying centerline + cap removal."
                )
                mhd_flip = flip_image_native_z(mhd_proc_img)
                mesh_flip = flip_mesh_native_z(mesh_proc, mhd_proc_img)

                pos2, tang2, T2, resid2, endpoint_hint2 = generate_cl(
                    mhd_flip,
                    n_points=1200,
                    endpoint_axis_mode="auto",
                    return_endpoint_hint=True,
                )
                if orient_info["changed"]:
                    pos2, tang2 = permute_centerline_axes(pos2, tang2, orient_info["order"])
                    endpoint_hint2 = permute_endpoint_hint_axes(endpoint_hint2, orient_info["order"])
                remeshed2 = cut_clean(mesh_flip,
                                      pos2,
                                      buffer_mm=buffer_mm,
                                      mhd_iso_img=mhd_flip,
                                      target_faces=target_faces,
                                      endpoint_hint=endpoint_hint2,
                                      )

                # Per request, use second-pass result regardless of final openness.
                remeshed = remeshed2
                resid = float(resid2)
                final_pos = pos2
                final_tang = tang2

            # Keep the saved output in the reoriented frame so downstream
            # centerline/cap assumptions (long axis == z) continue to hold.

            # Step 5: Save MAT with suffix derived from target_faces
            mat_suffix = f"_o{target_faces // 1000}K"
            out_mat_file = mat_file.parent / f"{scan_id}{mat_suffix}.mat"
            out_mat_file.parent.mkdir(parents=True, exist_ok=True)
            mdict = {
                f"{scan_id}_surface_vertices":  remeshed.vertices,
                f"{scan_id}_surface_triangles": remeshed.faces + 1,
            }
            savemat(str(out_mat_file), mdict)

            rows.append({
                "cohort":          job["cohort"],
                "patient_id":      job["patient_id"],
                "scan_id":         scan_id,
                "mat_file":        str(out_mat_file),
                "n_faces_in":      int(mesh_in.faces.shape[0]),
                "n_faces_out":     int(remeshed.faces.shape[0]),
                "euler_number":    int(remeshed.euler_number),
                "reg_residual_mm": round(float(resid), 4),
                "centerline_points_xyz": final_pos.tolist() if include_centerline_in_log else None,
                "centerline_tangents_xyz": final_tang.tolist() if include_centerline_in_log else None,
                "status":          "ok",
                "error":           "",
            })

        except Exception as e:
            print(f"  ERROR: {e}")
            rows.append({
                "cohort":          job["cohort"],
                "patient_id":      job["patient_id"],
                "scan_id":         scan_id,
                "mat_file":        str(out_mat_file),
                "n_faces_in":      np.nan,
                "n_faces_out":     np.nan,
                "euler_number":    np.nan,
                "reg_residual_mm": float(resid) if resid == resid else np.nan,
                "status":          "error",
                "error":           str(e),
            })
            if not continue_on_error:
                raise

    results_df = pd.DataFrame(rows)
    if not results_df.empty:
        print("\nNew-preprocessing batch summary")
        print(results_df["status"].value_counts(dropna=False))

    return results_df