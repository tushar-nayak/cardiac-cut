# CAMP Code

import numpy as np
import pandas as pd 
import trimesh
from scipy.io import savemat, loadmat
import plotly.graph_objects as go
from pathlib import Path

def PerformRotation(xyz, M):
    x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
    xr = M[0,0]*x + M[0,1]*y + M[0,2]*z
    yr = M[1,0]*x + M[1,1]*y + M[1,2]*z
    zr = M[2,0]*x + M[2,1]*y + M[2,2]*z
    return xr, yr, zr

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

def _collect_mat_jobs_for_cohorts(base_root, cohorts):
    """Collect MAT files from cohort/patient/scan/mesh layout.

    Expected layout:
      base_root/cohort/patient/scan/mesh/<scan_id>.mat
    """
    base_root = Path(base_root)
    jobs = []

    for cohort in cohorts:
        cohort_dir = base_root / cohort
        if not cohort_dir.exists():
            print(f"[warn] Cohort folder not found: {cohort_dir}")
            continue

        for patient_dir in sorted([p for p in cohort_dir.iterdir() if p.is_dir()]):
            for scan_dir in sorted([s for s in patient_dir.iterdir() if s.is_dir()]):
                scan_id = scan_dir.name
                mesh_dir = scan_dir / "mesh"
                if not mesh_dir.exists() or not mesh_dir.is_dir():
                    continue

                preferred_mat = mesh_dir / f"{scan_id}.mat"
                if preferred_mat.exists():
                    mat_file = preferred_mat
                else:
                    mats = sorted(mesh_dir.glob("*.mat"))
                    if not mats:
                        continue
                    mat_file = mats[0]

                jobs.append(
                    {
                        "cohort": cohort,
                        "patient_id": patient_dir.name,
                        "scan_id": scan_id,
                        "mat_file": mat_file,
                    }
                )

    return jobs

def PlaneFinder(x,y,z,tri_normal=None):
    ones=np.ones_like(x)
    n=np.dot(ones,ones)
    xav=np.dot(ones,x)/n
    yav=np.dot(ones,y)/n
    zav=np.dot(ones,z)/n
    svd_input = np.array([x-xav,y-yav,z-zav])
    normal = np.linalg.svd(svd_input)[0].T[2]
    if tri_normal is not None and np.dot(normal, tri_normal) < 0:
        normal *= -1
    axis = np.argmax(np.abs(normal))
    return svd_input, normal, axis

def RotateVec(x1, y1, z1, x2, y2, z2):
    v1 = np.array([x1, y1, z1])
    v2 = np.array([x2, y2, z2])
    v1_unit = v1 / np.linalg.norm(v1)
    v2_unit = v2 / np.linalg.norm(v2)

    if np.allclose(v1_unit, v2_unit):
        return np.eye(3)
    if np.allclose(v1_unit, -v2_unit):
        ortho = np.array([1, 0, 0]) if not np.allclose(v1_unit, [1, 0, 0]) else np.array([0, 1, 0])
        axis = np.cross(v1_unit, ortho)
    else:
        axis = np.cross(v1_unit, v2_unit)

    axis /= np.linalg.norm(axis)
    angle = np.arccos(np.clip(np.dot(v1_unit, v2_unit), -1.0, 1.0))

    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    x, y, z = axis
    return np.array([
        [t*x*x + c,    t*x*y - s*z,  t*x*z + s*y],
        [t*x*y + s*z,  t*y*y + c,    t*y*z - s*x],
        [t*x*z - s*y,  t*y*z + s*x,  t*z*z + c]])

def Polyfit2dZ(xr, yr, zr):
    U = np.column_stack([np.ones_like(xr), xr, yr, xr**2, xr*yr, yr**2])
    V = np.dot(U.T, zr)
    M = np.dot(U.T, U)
    coefs = np.linalg.solve(M, V)
    c1, c2 = coefs[1], coefs[2]
    normal_mag = np.sqrt(1 + c1**2 + c2**2)
    nx, ny, nz = -c1 / normal_mag, -c2 / normal_mag, 1 / normal_mag
    correction_factor = abs(nz)
    err = np.sum(((zr - np.dot(U, coefs)) / correction_factor)**2) / len(xr)
    return coefs, err, nx, ny, nz

def evaluate_z(x, y, coefs):
    return coefs[0] + coefs[1]*x + coefs[2]*y + coefs[3]*x**2 + coefs[4]*x*y + coefs[5]*y**2

def MongeGaussian(coefs):
    """
    The coefs satisfy z=coefs[0]+coefs[1]*x+coefs[2]*y+coefs[3]*x**2+coefs[4]*x*y+coefs[5]*y**2
    the implicit assumption is that the origin of the x,y axis is the point of interest.
    """
    num=4*coefs[3]*coefs[5]-coefs[4]**2
    denom=(1+coefs[1]**2+coefs[2]**2)**2
    K=num/denom
    return K

def MongeMean(coefs):
    """
    The coefs satisfy z=coefs[0]+coefs[1]*x+coefs[2]*y+coefs[3]*x**2+coefs[4]*x*y+coefs[5]*y**2
    the implicit assumption is that the origin of the x,y axis is the point of interest.
    """
    num=(1+coefs[1]**2)*coefs[5]+(1+coefs[2]**2)*coefs[3]-coefs[1]*coefs[2]*coefs[4]
    denom=(1+coefs[1]**2+coefs[2]**2)**(1.5)
    H=num/denom
    return H

def face_patch_labels(faces, vertex_labels, n_patches):
    """Assign each face to the majority patch label of its vertices."""
    out = np.zeros(len(faces), dtype=int)
    for i, tri in enumerate(faces):
        labels = vertex_labels[tri]
        labels = labels[labels >= 0]
        if labels.size == 0:
            out[i] = 0
            continue
        cnt = np.bincount(labels.astype(int), minlength=n_patches)
        out[i] = int(np.argmax(cnt))
    return out

def make_plotly_patch_figure(
    mesh,
    face_values,
    boundary_segments,
    title,
    color_bounds=None,
    color_scale="Turbo",
    face_mask=None,
):
    """Plot patch-colored surface with optional clipping bounds.

    face_mask: optional boolean array of shape (n_faces,). When provided, only
               faces where face_mask is True are rendered (excluded faces are
               hidden from the plot entirely).
    """
    v = mesh.vertices
    f = mesh.faces
    face_values = np.asarray(face_values, dtype=float)

    if face_mask is not None:
        face_mask = np.asarray(face_mask, dtype=bool)
        f = f[face_mask]
        face_values = face_values[face_mask]

    cmin = None
    cmax = None
    values_for_plot = face_values.copy()

    if color_bounds is not None:
        lo, hi = float(color_bounds[0]), float(color_bounds[1])
        if lo > hi:
            lo, hi = hi, lo
        values_for_plot = np.clip(values_for_plot, lo, hi)
        cmin, cmax = lo, hi

    mesh_trace = go.Mesh3d(
        x=v[:, 0],
        y=v[:, 1],
        z=v[:, 2],
        i=f[:, 0],
        j=f[:, 1],
        k=f[:, 2],
        intensity=values_for_plot,
        intensitymode="cell",
        colorscale=color_scale,
        cmin=cmin,
        cmax=cmax,
        colorbar=dict(title="Patch metric"),
        opacity=1.0,
        name="Patch metric",
    )

    bx, by, bz = [], [], []
    for p0, p1 in boundary_segments:
        bx.extend([p0[0], p1[0], None])
        by.extend([p0[1], p1[1], None])
        bz.extend([p0[2], p1[2], None])

    border_trace = go.Scatter3d(
        x=bx, y=by, z=bz,
        mode="lines",
        line=dict(color="black", width=5),
        name="Patch borders",
        showlegend=False,
    )

    fig = go.Figure(data=[mesh_trace, border_trace])
    fig.update_layout(
        title=title,
        scene=dict(aspectmode="data", xaxis_title="X", yaxis_title="Y", zaxis_title="Z"),
        margin=dict(l=0, r=0, b=0, t=40),
        width=500, 
        height=500,
    )
    return fig

def patch_pass(
    verts, faces, fit_monge_patch, patch_needs_split, patch_priority,
    face_patch_labels_fn, edge_connected_face_labels_fn,
    seed_patches, min_patch_size,
    enforce_edge_connected_faces, min_face_component_size,
    max_children=8, child_rule="fixed", eta_thresh=None,
    split_disconnected=True, fan_thresh=None, orphan_adjacency=False, seed=42,
):
    """
    Convergent adaptive patch algorithm.

    Splits until ALL patches either pass both quality thresholds
    (R2 >= R2_thresh AND eta <= eta_thresh) or are too small to split further
    (< 2 * min_patch_size vertices).  No target_patches cap.

    Returns (leaf_patches, f_labels, patch_list, n_splits_done).
    """
    from sklearn.cluster import KMeans
    import heapq

    n_verts = len(verts)

    # A PATCH MUST BE ONE PIECE OF SURFACE.
    #
    # KMeans clusters vertex COORDINATES. It has no notion of the surface, so two
    # sheets that are close in space but far apart across the surface land in the same
    # cluster -- at the base of a branch, the branch wall and the aortic wall beside it
    # can be 1.3 mm apart through the air and half a circumference apart along the
    # surface. Measured on the arch cohort: 2.3% of patches came out in two or more
    # pieces, one of them 32 faces and 16 faces sharing no vertex at all.
    #
    # A Monge patch fitted across that gap is not describing a neighbourhood of
    # anything. So split every cluster into its connected components BEFORE fitting,
    # which also keeps the per-patch arrays aligned with the labels -- doing it
    # afterwards would add labels the curvature arrays have no entry for.
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components as _cc

    _F = np.asarray(faces, dtype=int)
    _e = np.vstack([_F[:, [0, 1]], _F[:, [1, 2]], _F[:, [2, 0]]])
    _adj = coo_matrix(
        (np.ones(2 * len(_e), dtype=np.int8),
         (np.concatenate([_e[:, 0], _e[:, 1]]),
          np.concatenate([_e[:, 1], _e[:, 0]]))),
        shape=(n_verts, n_verts)).tocsr()

    def connected_parts(idx):
        """The pieces of a vertex set that are actually joined across the surface."""
        if not split_disconnected or len(idx) < 2:
            return [idx]
        n_comp, lab = _cc(_adj[idx][:, idx], directed=False)
        if n_comp <= 1:
            return [idx]
        return [idx[lab == c] for c in range(n_comp)]

    # A PATCH MUST ALSO BE SINGLE-VALUED.
    #
    # A Monge patch is a height field over its tangent plane, so every point of it has
    # to project onto that plane uniquely. Measure how far the surface normals fan out
    # from the patch's own mean normal: on a tube the normals rotate with the surface,
    # so the fan is half the arc the patch spans. Past 90 degrees some of the patch is
    # tilted beyond vertical relative to its own tangent plane, the projection folds
    # over, and no height field exists -- the fit then returns whatever compromises
    # between the two sheets, which is where the huge fake curvatures come from.
    #
    # Measured on the nub cohort: patches over 90 degrees are 0.6% of the total and
    # carry 35% of delta K; on SA21_1, a paediatric aorta that should be the smoothest
    # case in the cohort, two such patches carry 91% of it.
    _fv = np.zeros((n_verts, 3), dtype=float)
    _tri = np.asarray(verts)[_F]
    _fnrm = np.cross(_tri[:, 1] - _tri[:, 0], _tri[:, 2] - _tri[:, 0])
    for _k in range(3):
        np.add.at(_fv, _F[:, _k], _fnrm)
    _fv /= np.maximum(np.linalg.norm(_fv, axis=1, keepdims=True), 1e-12)

    def patch_fan(idx):
        """Largest angle, in degrees, between a patch's normals and their mean."""
        if len(idx) < 3:
            return 0.0
        n = _fv[idx]
        m = n.mean(axis=0)
        nm = np.linalg.norm(m)
        if nm < 1e-9:
            return 180.0
        return float(np.degrees(np.arccos(np.clip(n @ (m / nm), -1.0, 1.0))).max())

    # Offset by 997 per step away from the default so the per-depth seeds of one run
    # cannot collide with another's, and so seed=42 reproduces the original numbers
    # exactly -- the split seed was depth+1 before this was a parameter at all.
    #
    # The seed is a parameter so a case can be partitioned several independent ways.
    # delta K is carried by ~22 patches out of ~850, and which patches those are
    # depends on where k-means lands, so one partition is one sample of a noisy
    # process rather than the answer.
    km_global = KMeans(n_clusters=seed_patches, n_init=5, random_state=seed)
    labels = km_global.fit_predict(verts)

    patches = []
    for p in range(seed_patches):
        for idx in connected_parts(np.where(labels == p)[0]):
            if len(idx) < min_patch_size:
                continue
            coefs, R2, eta, M, center, K_center, H_center = fit_monge_patch(idx)
            patches.append({
                'indices': idx, 'coefs': coefs, 'R2': R2, 'eta': eta,
                'M': M, 'center': center, 'K_center': K_center, 'H_center': H_center,
                'fan': patch_fan(idx),
                'depth': 0, 'children': [],
            })

    patch_list = patches.copy()
    split_heap = []
    heap_counter = 0
    for patch in patches:
        if patch_needs_split(patch):
            heapq.heappush(split_heap, (-patch_priority(patch), heap_counter, patch))
            heap_counter += 1

    n_splits_done = 0
    while split_heap:
        _, _, patch = heapq.heappop(split_heap)
        if not patch_needs_split(patch):
            continue

        idx = patch['indices']
        local = verts[idx] - patch['center']
        xr, yr, zr = PerformRotation(local, patch['M'])
        # SPLIT WHERE THE SHEETS DIFFER, not where they overlap.
        #
        # Clustering in (x, y) alone discards the height -- and height is exactly what
        # separates the two sheets of a wrapped patch, which lie on top of each other
        # in that projection. KMeans then cuts the overlapping disc into pie slices and
        # every child inherits both sheets, so it splits forever without ever resolving
        # the wrap, until it hits the vertex floor.
        #
        # Adding the height and the surface normal gives it the coordinates that DO
        # differ: on a wrap the two sheets face opposite ways. The normal is scaled by
        # the patch's own extent so the two halves of the feature are commensurate,
        # which leaves a flat patch clustering essentially as before.
        pts2d = np.stack([xr, yr], axis=1)
        if fan_thresh is not None:
            nx, ny, nz = PerformRotation(_fv[idx], patch['M'])
            scale = float(np.sqrt(np.mean(xr ** 2 + yr ** 2))) or 1.0
            pts2d = np.column_stack([xr, yr, zr,
                                     scale * nx, scale * ny, scale * nz])

        # How many children does this patch NEED?
        #   'fixed'    : the original rule -- as many children as the patch can AFFORD.
        #                A large patch always splits max_children ways, so it overshoots
        #                its eta-permitted size and lands on the min_patch_size floor.
        #   'adaptive' : eta ~ r^1.78 (sphere calibration), so patch AREA must shrink by
        #                (eta/eta_thresh)^(2/1.78) = ^1.12.  That ratio IS the number of
        #                children needed to reach the gate in a single split.
        cap = max(2, min(len(idx) // min_patch_size, int(max_children)))
        if child_rule == "adaptive":
            if eta_thresh and patch['eta'] is not None and patch['eta'] > eta_thresh:
                need = (patch['eta'] / eta_thresh) ** 1.12
            else:
                need = 2.0          # R2-only failure: no size target, step down gently
            n_children = int(np.clip(round(need), 2, cap))
        else:
            n_children = cap
        if n_children < 2:
            continue

        km = KMeans(n_clusters=n_children, n_init=3, random_state=int(patch['depth'] + 1
                                        + 997 * abs(int(seed) - 42)))
        sub_labels = km.fit_predict(pts2d)

        valid_children = 0
        for c in range(n_children):
            # Children are clustered in the parent's 2D tangent projection, which can
            # separate a cluster into pieces that do not touch just as easily as the
            # seeding can. Same rule applies to them.
            for sub_idx in connected_parts(idx[sub_labels == c]):
                if len(sub_idx) < min_patch_size:
                    continue
                coefs, R2, eta, M2, center2, K_center2, H_center2 = fit_monge_patch(sub_idx)
                child = {
                    'indices': sub_idx, 'coefs': coefs, 'R2': R2, 'eta': eta,
                    'M': M2, 'center': center2, 'K_center': K_center2,
                    'H_center': H_center2, 'fan': patch_fan(sub_idx),
                    'depth': patch['depth'] + 1, 'children': [],
                }
                patch['children'].append(child)
                patch_list.append(child)
                valid_children += 1
                if patch_needs_split(child):
                    heapq.heappush(split_heap,
                                   (-patch_priority(child), heap_counter, child))
                    heap_counter += 1
        if valid_children > 0:
            n_splits_done += 1

    leaf_patches = [p for p in patch_list if not p['children']]
    n_leaf = len(leaf_patches)
    if n_leaf == 0:
        raise ValueError("No valid leaf patches were produced.")

    final_labels = np.full(n_verts, -1, dtype=int)
    for i, patch in enumerate(leaf_patches):
        final_labels[patch['indices']] = i

    missing = np.where(final_labels < 0)[0]
    if len(missing) > 0:
        pending = missing
        if orphan_adjacency:
            # A split whose children fall under min_patch_size discards them, leaving
            # their vertices unlabelled. Sending those to the leaf with the nearest
            # CENTRE is surface-blind -- next to a branch ostium the nearest centre is
            # routinely across the gap, so a patch on the aortic wall can adopt branch
            # wall facing the other way. Filling from the labelled boundary inward
            # instead keeps every orphan on the surface it came from. Sweeps are
            # batched rather than applied as they are found, so the outcome does not
            # depend on vertex order.
            for _ in range(64):
                if len(pending) == 0:
                    break
                chosen, still = {}, []
                for v in pending:
                    nb = _adj.indices[_adj.indptr[v]:_adj.indptr[v + 1]]
                    lb = final_labels[nb]
                    lb = lb[lb >= 0]
                    if len(lb) == 0:
                        still.append(v)
                        continue
                    vals, cnts = np.unique(lb, return_counts=True)
                    chosen[int(v)] = int(vals[np.argmax(cnts)])
                if not chosen:
                    break
                for v, L in chosen.items():
                    final_labels[v] = L
                pending = np.asarray(still, dtype=int)
        # Nearest centre for whatever is left: everything, when the flag is off, and
        # only vertices touching nothing labelled when it is on.
        if len(pending) > 0:
            leaf_centers = np.vstack([p['center'] for p in leaf_patches])
            delta = verts[pending, None, :] - leaf_centers[None, :, :]
            nearest = np.argmin(np.sum(delta * delta, axis=2), axis=1)
            final_labels[pending] = nearest.astype(int)

    f_labels = face_patch_labels_fn(faces, final_labels, n_leaf)
    if enforce_edge_connected_faces:
        f_labels, _ = edge_connected_face_labels_fn(
            trimesh.Trimesh(vertices=verts, faces=faces, process=False),
            f_labels, n_patches=n_leaf,
            min_component_faces=min_face_component_size,
        )

    return leaf_patches, f_labels, patch_list, n_splits_done

def camp(
    mat_path,
    scan_prefix=None,
    color_metric="integrated_K",
    R2_thresh=0.9,
    eta_thresh=0.075,
    min_patch_size=18,
    seed_patches=100,
    trim_frac=0.0,
    color_bounds=None,
    color_scale="Turbo",
    enforce_edge_connected_faces=True,
    min_face_component_size=2,
    split_disconnected=True,
    fan_thresh=90.0,
    orphan_adjacency=False,
    seed=42,
    max_children=8,
    child_rule="fixed",
    make_figure=True,
):
    """
    Convergent Monge-based adaptive patch partitioning (no target_patches).

    Runs KMeans-seed + priority-queue splitting until all patches satisfy
    both R2 and eta thresholds, or are too small to split.  The final patch
    count is determined entirely by the data and quality thresholds.
    """
    mesh, keys = load_surface_from_mat(mat_path, scan_prefix=scan_prefix)
    verts = mesh.vertices.copy()
    faces = mesh.faces.copy()
    print(f"[input] Loaded mesh: {len(verts)} verts, {len(faces)} faces")

    # ── inner helpers ─────────────────────────────────────────────────────────

    def fit_monge_patch(idx):
        verts_patch = verts[idx]
        if len(verts_patch) < 6:
            return None, None, None, None, None, None, None
        center = verts_patch.mean(axis=0)
        local = verts_patch - center
        try:
            _, normal, _ = PlaneFinder(local[:, 0], local[:, 1], local[:, 2])
        except Exception:
            return None, None, None, None, None, None, None
        M = RotateVec(*normal, 0, 0, 1)
        xr, yr, zr = PerformRotation(local, M)
        try:
            coefs, err, nx, ny, nz = Polyfit2dZ(xr, yr, zr)
        except Exception:
            return None, None, None, None, None, None, None
        z_pred = evaluate_z(xr, yr, coefs)
        SS_res = np.sum((zr - z_pred) ** 2)
        SS_tot = np.sum((zr - np.mean(zr)) ** 2)
        R2 = float(1.0 - SS_res / (SS_tot + 1e-12))

        def _shifted(c, x0, y0):
            c0, c1, c2, c3, c4, c5 = c
            return np.array([
                c0+c1*x0+c2*y0+c3*x0**2+c4*x0*y0+c5*y0**2,
                c1+2*c3*x0+c4*y0, c2+c4*x0+2*c5*y0,
                c3, c4, c5,
            ], dtype=float)

        K_vals = np.array([
            float(MongeGaussian(_shifted(coefs, float(xi), float(yi))))
            for xi, yi in zip(xr, yr)
        ])
        K_center = float(MongeGaussian(coefs))
        H_center = float(MongeMean(coefs))
        eta = float(np.std(K_vals) / (abs(K_center) + 1e-12))
        return coefs, R2, eta, M, center, K_center, H_center

    def patch_needs_split(patch):
        if len(patch['indices']) < 2 * min_patch_size:
            return False
        if patch['M'] is None:
            return False
        r2_ok  = (patch['R2']  is not None) and (patch['R2']  >= R2_thresh)
        eta_ok = (patch['eta'] is not None) and (patch['eta'] <= eta_thresh)
        # A wrapped patch is not a poor fit, it is a patch no height field can
        # describe. R2 and eta can both pass on one, so it needs its own gate.
        fan_ok = (fan_thresh is None) or (patch.get('fan') is None)                  or (patch['fan'] <= fan_thresh)
        return not (r2_ok and eta_ok and fan_ok)

    def patch_priority(patch):
        r2_t  = 1.0 if patch['R2']  is None else min(1.0, max(0.0, (R2_thresh  - patch['R2'])  / max(R2_thresh,  1e-12)))
        eta_t = 1.0 if patch['eta'] is None else min(1.0, max(0.0, (patch['eta'] - eta_thresh) / max(eta_thresh, 1e-12)))
        return float(r2_t + eta_t)

    def edge_connected(mesh_obj, face_labels, n_patches, min_component_faces=2):
        face_labels = np.asarray(face_labels, dtype=int).copy()
        n_faces = len(face_labels)
        neigh = [[] for _ in range(n_faces)]
        for f0, f1 in mesh_obj.face_adjacency:
            neigh[int(f0)].append(int(f1))
            neigh[int(f1)].append(int(f0))
        changed = 0
        for lbl in range(int(n_patches)):
            members = np.where(face_labels == lbl)[0]
            if members.size == 0:
                continue
            member_set = set(members.tolist())
            seen = set()
            components = []
            for seed in members:
                if int(seed) in seen:
                    continue
                comp = []
                stack = [int(seed)]
                seen.add(int(seed))
                while stack:
                    f = stack.pop()
                    comp.append(f)
                    for nb in neigh[f]:
                        if nb in member_set and nb not in seen:
                            seen.add(nb)
                            stack.append(nb)
                components.append(comp)
            if len(components) <= 1:
                continue
            components.sort(key=len, reverse=True)
            for comp in components[1:]:
                if len(comp) >= min_component_faces:
                    continue
                for f in comp:
                    votes = {}
                    for nb in neigh[f]:
                        nb_lbl = int(face_labels[nb])
                        if nb_lbl == lbl:
                            continue
                        votes[nb_lbl] = votes.get(nb_lbl, 0) + 1
                    if votes:
                        new_lbl = max(votes.items(), key=lambda kv: kv[1])[0]
                    else:
                        tri_c = np.mean(mesh_obj.vertices[mesh_obj.faces[f]], axis=0)
                        best_lbl, best_d2 = lbl, np.inf
                        for cl in range(int(n_patches)):
                            cf = np.where(face_labels == cl)[0]
                            if cf.size == 0:
                                continue
                            cc = np.mean(mesh_obj.triangles_center[cf], axis=0)
                            d2 = float(np.sum((tri_c - cc) ** 2))
                            if d2 < best_d2:
                                best_d2, best_lbl = d2, cl
                        new_lbl = int(best_lbl)
                    if new_lbl != lbl:
                        face_labels[f] = int(new_lbl)
                        changed += 1
        return face_labels, changed

    # ── convergent patch pass ─────────────────────────────────────────────────
    leaf_patches, f_labels, _, n_splits = patch_pass(
        verts, faces,
        fit_monge_patch=fit_monge_patch,
        patch_needs_split=patch_needs_split,
        patch_priority=patch_priority,
        face_patch_labels_fn=face_patch_labels,
        edge_connected_face_labels_fn=edge_connected,
        seed_patches=seed_patches,
        min_patch_size=min_patch_size,
        enforce_edge_connected_faces=enforce_edge_connected_faces,
        min_face_component_size=min_face_component_size,
        max_children=max_children,
        child_rule=child_rule,
        eta_thresh=eta_thresh,
        split_disconnected=split_disconnected,
        fan_thresh=fan_thresh,
        orphan_adjacency=orphan_adjacency,
        seed=seed,
    )
    n_leaf = len(leaf_patches)

    # ── every patch must be ONE piece of surface, measured on FACES ───────────
    #
    # Splitting the vertex clusters before fitting is not enough. face_patch_labels
    # assigns each face from its three vertices, so a face straddling two patches is
    # resolved by a rule that can leave a label's FACE set in pieces even when its
    # vertex set was connected -- and two triangle fans meeting at a single shared
    # vertex are vertex-connected but face-disconnected by definition.
    #
    # So do it again here, on the labels that actually decide the areas. A stray
    # component big enough to fit becomes its own patch, with its own Monge fit, since
    # its curvature is a different measurement from the piece it was lumped with. One
    # too small to fit is handed to the neighbour it shares the most edges with, which
    # is always adjacent to it and therefore cannot leave that neighbour disconnected.
    if split_disconnected:
        F_arr = np.asarray(faces, dtype=int)
        neigh = [[] for _ in range(len(F_arr))]
        for f0, f1 in mesh.face_adjacency:
            neigh[int(f0)].append(int(f1))
            neigh[int(f1)].append(int(f0))

        def _face_components(members):
            member_set = set(members.tolist())
            seen, comps = set(), []
            for seed in members:
                s0 = int(seed)
                if s0 in seen:
                    continue
                comp, stack = [], [s0]
                seen.add(s0)
                while stack:
                    f = stack.pop()
                    comp.append(f)
                    for nb in neigh[f]:
                        if nb in member_set and nb not in seen:
                            seen.add(nb)
                            stack.append(nb)
                comps.append(comp)
            return comps

        leaf_patches = list(leaf_patches)
        made, merged = 0, 0
        for _sweep in range(6):          # a merge can disconnect nothing, but re-check
            moved = 0
            for lbl in range(len(leaf_patches)):
                members = np.where(f_labels == lbl)[0]
                if members.size == 0:
                    continue
                comps = _face_components(members)
                if len(comps) <= 1:
                    continue
                comps.sort(key=len, reverse=True)
                for comp in comps[1:]:
                    vids = np.unique(F_arr[comp].ravel())
                    if len(vids) >= min_patch_size:
                        try:
                            c2, R22, eta2, M2, cen2, Kc2, Hc2 = fit_monge_patch(vids)
                        except Exception:
                            c2 = None
                        if c2 is not None:
                            leaf_patches.append({
                                'indices': vids, 'coefs': c2, 'R2': R22, 'eta': eta2,
                                'M': M2, 'center': cen2, 'K_center': Kc2,
                                'H_center': Hc2, 'depth': 0, 'children': [],
                            })
                            f_labels[comp] = len(leaf_patches) - 1
                            made += 1
                            moved += 1
                            continue
                    votes = {}
                    for f in comp:
                        for nb in neigh[f]:
                            nl = int(f_labels[nb])
                            if nl != lbl:
                                votes[nl] = votes.get(nl, 0) + 1
                    if votes:
                        f_labels[comp] = max(votes.items(), key=lambda kv: kv[1])[0]
                        merged += 1
                        moved += 1
            if moved == 0:
                break
        n_leaf = len(leaf_patches)
        if made or merged:
            print(f"[connected] {made} stray component(s) became their own patch, "
                  f"{merged} too small and merged into a neighbour")

    # ── per-patch statistics ──────────────────────────────────────────────────
    face_areas = mesh.area_faces
    avg_kg     = np.zeros(n_leaf, dtype=float)
    patch_area = np.zeros(n_leaf, dtype=float)
    eta_patch  = np.full(n_leaf, np.nan, dtype=float)
    R2_patch   = np.full(n_leaf, np.nan, dtype=float)
    for i, patch in enumerate(leaf_patches):
        K_c = patch['K_center']
        avg_kg[i]     = float(K_c) if K_c is not None else 0.0
        patch_area[i] = float(np.sum(face_areas[f_labels == i]))
        eta_patch[i]  = float(patch['eta']) if patch['eta'] is not None else np.nan
        R2_patch[i]   = float(patch['R2'])  if patch['R2']  is not None else np.nan
    integrated_K = avg_kg * patch_area

    # ── trimming on integrated K (patch-wise), not discrete K ────────────────
    keep_patch = _trim_mask_integrated(integrated_K, trim_frac=trim_frac)
    if not np.any(keep_patch):
        print("[trim] All patches were trimmed; reverting to keep-all.")
        keep_patch = np.ones_like(integrated_K, dtype=bool)

    n_patches_trimmed = int(np.count_nonzero(~keep_patch))
    area_pre = float(np.sum(patch_area))
    area_kept = float(np.sum(patch_area[keep_patch]))
    area_frac_trimmed = float(1.0 - (area_kept / area_pre)) if area_pre > 0 else 0.0

    avg_kg_kept = avg_kg[keep_patch]
    patch_area_kept = patch_area[keep_patch]
    integrated_K_kept = integrated_K[keep_patch]
    eta_patch_kept = eta_patch[keep_patch]
    R2_patch_kept = R2_patch[keep_patch]
    leaf_patches_kept = [p for i, p in enumerate(leaf_patches) if keep_patch[i]]
    
    # ── convergence diagnostics ───────────────────────────────────────────────
    n_failed_quality = sum(
        1 for p in leaf_patches
        if not ((p['R2'] is not None and p['R2'] >= R2_thresh) and
                (p['eta'] is not None and p['eta'] <= eta_thresh))
    )
    n_at_min_size = sum(
        1 for p in leaf_patches
        if len(p['indices']) < 2 * min_patch_size
    )
    pct_failed_quality = 100.0 * n_failed_quality / n_leaf if n_leaf > 0 else np.nan
    pct_at_min_size    = 100.0 * n_at_min_size    / n_leaf if n_leaf > 0 else np.nan

    print(f"[converge] failed quality (R2+eta): {n_failed_quality}/{n_leaf} ({pct_failed_quality:.1f}%)")
    print(f"[converge] at min size (unsplittable): {n_at_min_size}/{n_leaf} ({pct_at_min_size:.1f}%)")
    print(f"[done] {n_leaf} patches after splits={n_splits}")

    # ── boundary segments + figure (skipped when make_figure=False) ───────────
    if make_figure:
        face_keep = keep_patch[f_labels]
        boundary_segments = []
        for f0, f1 in mesh.face_adjacency:
            if not (face_keep[f0] or face_keep[f1]):
                continue
            if f_labels[f0] == f_labels[f1]:
                continue
            shared = np.intersect1d(mesh.faces[f0], mesh.faces[f1])
            if len(shared) == 2:
                a, b = shared
                boundary_segments.append((verts[a], verts[b]))

        metric_map = {
            "integrated_K": integrated_K,
            "avg_kg":        avg_kg,
            "eta":           eta_patch,
            "R2":            R2_patch,
            "r2":            R2_patch,
        }
        if color_metric not in metric_map:
            raise ValueError(f"Unknown color_metric='{color_metric}'. "
                             f"Use one of {list(metric_map.keys())}.")
        face_values = metric_map[color_metric][f_labels]

        fig = make_plotly_patch_figure(
            mesh,
            face_values,
            boundary_segments,
            title=f"CAMP — {color_metric}",
            color_bounds=color_bounds,
            color_scale=color_scale,
            face_mask=face_keep,
        )
    else:
        fig = None

    print(f"[trim] trimmed patches: {n_patches_trimmed}/{n_leaf} (area frac={area_frac_trimmed:.4f})")
    print(f"[kg] avg_kg      : [{np.nanmin(avg_kg_kept):.4g}, {np.nanmax(avg_kg_kept):.4g}]")
    print(f"[K] integrated_K: [{np.nanmin(integrated_K_kept):.4g}, {np.nanmax(integrated_K_kept):.4g}]")
    print(f"[eta] range     : [{np.nanmin(eta_patch_kept):.4g}, {np.nanmax(eta_patch_kept):.4g}]")
    print(f"[R2] range      : [{np.nanmin(R2_patch_kept):.4g}, {np.nanmax(R2_patch_kept):.4g}]")

    return {
        "mesh":               mesh,
        "vertex_labels":      np.full(len(verts), -1, dtype=int),
        "face_labels":        f_labels,
        "avg_kg":             avg_kg_kept,
        "patch_area":         patch_area_kept,
        "K_integrated_patch": integrated_K_kept,
        "kg_avg_patch":       avg_kg_kept,
        "eta_patch":          eta_patch_kept,
        "R2_patch":           R2_patch_kept,
        "pct_failed_quality": pct_failed_quality,
        "pct_at_min_size":    pct_at_min_size,
        "n_patches_trimmed":  n_patches_trimmed,
        "area_frac_trimmed":  area_frac_trimmed,
        "trim_keep_mask":     keep_patch,
        "figure":             fig,
        "patches":            leaf_patches_kept,
    }

# ── BATCH and FEATURIZATION CODE ──────────────────────────────

def _trim_mask_integrated(integrated_K, trim_frac=0.0):
    """Keep-mask by percentile trimming on integrated K values.

    trim_frac=f removes both tails outside [f, 1-f] quantiles.
    """
    q = np.asarray(integrated_K, dtype=float)
    n = len(q)
    keep = np.ones(n, dtype=bool)

    if trim_frac is None:
        return keep

    f = float(trim_frac)
    if f <= 0.0:
        return keep
    if f >= 0.5:
        raise ValueError("trim_frac must be in [0, 0.5).")

    finite = np.isfinite(q)
    if not np.any(finite):
        return keep

    lo, hi = np.quantile(q[finite], [f, 1.0 - f])
    keep = finite & (q >= lo) & (q <= hi)
    return keep


def camp_batch(
    cohorts,
    base_root=r"Z:\aorta\thoracic",
    r2_thresh=0.9,
    eta_thresh=0.075,
    min_patch_size=18,
    seed_patches=100,
    trim_frac=0.0,
    color_metric="integrated_K",
    color_bounds=None,
    color_scale="Turbo",
    continue_on_error=True,
    mat_suffix=None,
    make_figure=False,
):
    """Run convergent adaptive patching on saved MAT meshes and aggregate stats.

    Parameters
    ----------
    mat_suffix : str or None
        If provided (e.g. '_o50'), only MAT files whose stem contains this
        substring are considered.  If None, _collect_mat_jobs_for_cohorts is used.

    Returns a dataframe with:
      total area, integrated-K fluctuation (area-weighted variance),
      and trimming diagnostics.
    """
    from pathlib import Path

    if mat_suffix is not None:
        base_root_path = Path(base_root)
        jobs = []
        for cohort in cohorts:
            cohort_dir = base_root_path / cohort
            if not cohort_dir.exists():
                print(f"[warn] Cohort folder not found: {cohort_dir}")
                continue
            for patient_dir in sorted([p for p in cohort_dir.iterdir() if p.is_dir()]):
                for scan_dir in sorted([s for s in patient_dir.iterdir() if s.is_dir()]):
                    scan_id = scan_dir.name
                    mesh_dir = scan_dir / "mesh"
                    if not mesh_dir.exists():
                        continue
                    matches = sorted(f for f in mesh_dir.glob("*.mat") if mat_suffix in f.stem)
                    if not matches:
                        continue
                    jobs.append({
                        "cohort":     cohort,
                        "patient_id": patient_dir.name,
                        "scan_id":    scan_id,
                        "mat_file":   matches[0],
                    })
        print(f"Found {len(jobs)} MAT scans (suffix='{mat_suffix}') across cohorts: {cohorts}")
    else:
        jobs = _collect_mat_jobs_for_cohorts(base_root=base_root, cohorts=cohorts)
        print(f"Found {len(jobs)} MAT scans across cohorts: {cohorts}")

    rows = []
    errors = []

    for i, job in enumerate(jobs, start=1):
        scan_id = job["scan_id"]
        mat_file = Path(job["mat_file"])
        print(f"[{i}/{len(jobs)}] Convergent patching {scan_id}  ({mat_file.name})")

        try:
            result_i = camp(
                mat_path=str(mat_file),
                scan_prefix=None,
                color_metric=color_metric,
                R2_thresh=r2_thresh,
                eta_thresh=eta_thresh,
                min_patch_size=min_patch_size,
                seed_patches=seed_patches,
                trim_frac=trim_frac,
                color_bounds=color_bounds,
                color_scale=color_scale,
                make_figure=make_figure,
            )

            mesh_area    = float(result_i["mesh"].area)
            k_integrated = np.asarray(result_i["K_integrated_patch"], dtype=float)
            patch_areas  = np.asarray(result_i["patch_area"], dtype=float)
            n_patches    = int(len(k_integrated))

            valid = np.isfinite(k_integrated) & np.isfinite(patch_areas) & (patch_areas > 0)
            if not np.any(valid):
                raise ValueError("No valid patch values/areas for fluctuation calculation")

            q = k_integrated[valid]
            a = patch_areas[valid]
            area_total = float(np.sum(a))
            w = a / area_total
            q_bar = float(np.sum(w * q))
            k_fluct = float(np.sum(w * (q - q_bar) ** 2))

            euler_number = int(result_i["mesh"].euler_number)
            is_straw = (euler_number == 0)
            if not is_straw:
                print(f"  [topology] NOT a straw: euler_number={euler_number}")

            row = {
                "scan_id":                  scan_id,
                "surface_area_mm2":         mesh_area,
                "area_total":               area_total,
                "k_fluct":                  k_fluct,
                "n_patches_trimmed":        int(result_i["n_patches_trimmed"]),
                "area_frac_trimmed":        float(result_i["area_frac_trimmed"]),
                "n_patches":                n_patches,
                "euler_number":             euler_number,
                "is_straw":                 is_straw,
                "pct_failed_quality":       float(result_i["pct_failed_quality"]),
                "pct_at_min_size":          float(result_i["pct_at_min_size"]),
            }
            rows.append(row)

        except Exception as e:
            msg = f"{scan_id}: {e}"
            print(f"  ERROR: {msg}")
            errors.append(msg)
            if not continue_on_error:
                raise

    results_df = pd.DataFrame(rows)
    results_df.attrs["errors"] = errors

    if not results_df.empty:
        print("\nConvergent patch batch complete")
        print(f"Rows: {len(results_df)}")
        if "is_straw" in results_df.columns:
            n_bad = int((~results_df["is_straw"].astype(bool)).sum())
            print(f"Non-straw meshes (euler_number != 0): {n_bad}/{len(results_df)}")
    if errors:
        print(f"Errors: {len(errors)}")

    return results_df