# Centerline Functions

# Skeleton-based centerline extraction using a DT-weighted geodesic shortest path. 
# Includes BFS graph helpers, endpoint detection, and lumen-interior recentering.

import numpy as np
import SimpleITK as sitk
from scipy import ndimage
from skimage.morphology import skeletonize
import trimesh
from scipy.spatial import cKDTree
from skimage.graph import route_through_array
from collections import deque
from pathlib import Path
import pandas as pd

def _voxel_to_world(ijk, vox_spacing, offset, direction=None):
    """Convert voxel (k,j,i) coords to world (x,y,z).

    Formula: world[x,y,z] = direction @ ([i,j,k] * spacing) + offset
    ijk from argwhere is (k,j,i); reversed to (i,j,k) for the dot product.
    """
    if direction is None:
        direction = np.eye(3, dtype=float)
    direction = np.asarray(direction, dtype=float).reshape(3, 3)

    xyz = np.atleast_2d(ijk)[:, ::-1].astype(float)
    phys = xyz * np.asarray(vox_spacing, dtype=float)[None, :]
    world = phys @ direction.T + np.asarray(offset, dtype=float)[None, :]
    return world if np.ndim(ijk) > 1 else world[0]

def _world_to_voxel(world, vox_spacing, offset, direction=None):
    """Inverse of _voxel_to_world: world (x,y,z) -> voxel (k,j,i).

    For direction vectors (not points), pass offset=np.zeros(3).
    """
    if direction is None:
        direction = np.eye(3, dtype=float)
    direction = np.asarray(direction, dtype=float).reshape(3, 3)
    inv_direction = np.linalg.inv(direction)

    rel = np.atleast_2d(world).astype(float) - np.asarray(offset, dtype=float)[None, :]
    phys = rel @ inv_direction.T
    xyz = phys / np.asarray(vox_spacing, dtype=float)[None, :]
    ijk = xyz[:, ::-1]  # flip (i,j,k) -> (k,j,i)
    return ijk if np.ndim(world) > 1 else ijk[0]

def _skeleton_to_graph(skel_coords):
    """Build adjacency graph from skeleton voxel coordinates using 26-connectivity."""
    tree = cKDTree(skel_coords)
    pairs = tree.query_pairs(r=1.75)
    n = len(skel_coords)
    adj = [[] for _ in range(n)]
    for i, j in pairs:
        adj[i].append(j)
        adj[j].append(i)
    return adj

def _bfs_dist(adj, start):
    """BFS from a single node, return distance array (-1 = unreached)."""
    n = len(adj)
    dist = np.full(n, -1, dtype=int)
    dist[start] = 0
    q = deque([start])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist

def _connected_components(adj):
    """Return list of node-index lists, one per connected component."""
    n = len(adj)
    seen = np.zeros(n, dtype=bool)
    comps = []

    for s in range(n):
        if seen[s]:
            continue
        q = deque([s])
        seen[s] = True
        comp = []
        while q:
            u = q.popleft()
            comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    q.append(v)
        comps.append(comp)

    return comps

def _select_endpoint_pair(adj, skel_world):
    """Pick endpoints on the largest skeleton component via exact endpoint diameter."""
    if len(adj) == 0:
        raise ValueError("Empty skeleton graph.")

    comps = _connected_components(adj)
    comps.sort(key=len, reverse=True)
    largest = np.asarray(comps[0], dtype=int)

    deg = np.array([len(a) for a in adj], dtype=int)
    endpoints = largest[deg[largest] <= 1]
    if len(endpoints) < 2:
        endpoints = largest

    best_u = int(endpoints[0])
    best_v = int(endpoints[0])
    best_d = -1

    for s in endpoints:
        dist = _bfs_dist(adj, int(s))
        reachable = endpoints[dist[endpoints] >= 0]
        if len(reachable) == 0:
            continue
        t = int(reachable[np.argmax(dist[reachable])])
        d = int(dist[t])
        if d > best_d:
            best_d = d
            best_u = int(s)
            best_v = t

    if best_u == best_v and len(endpoints) >= 2:
        ep_world = skel_world[endpoints]
        best_u = int(endpoints[np.argmin(ep_world[:, 2])])
        best_v = int(endpoints[np.argmax(ep_world[:, 2])])

    return best_u, best_v, largest, endpoints

def _bfs_path(adj, start, end):
    """BFS from *start* to *end*, return path as list of node indices."""
    n = len(adj)
    dist = np.full(n, -1, dtype=int)
    parent = np.full(n, -1, dtype=int)
    dist[start] = 0
    parent[start] = start
    q = deque([start])
    while q:
        u = q.popleft()
        if u == end:
            break
        for v in adj[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                parent[v] = u
                q.append(v)
    path = []
    node = end
    while parent[node] != node:
        path.append(node)
        node = parent[node]
    path.append(node)
    path.reverse()
    return path

def _safe_normalize(v, axis=-1, eps=1e-12):
    """Normalize vectors along *axis*, guarding against zero-length."""
    norm = np.linalg.norm(v, axis=axis, keepdims=True)
    return v / np.maximum(norm, eps)

def _push_point_to_mask_boundary(
    mask,
    point_world,
    outward_world,
    vox_spacing,
    offset,
    direction=None,
    max_dist_mm=30.0,
    step_mm=None,
    outside_margin_mm=0.0,
    return_outside=False,
):
    """Move a point along outward direction to boundary (or just outside with margin)."""
    p = np.asarray(point_world, dtype=float)
    d = np.asarray(outward_world, dtype=float)
    n = np.linalg.norm(d)
    if n < 1e-9:
        return p
    d = d / n

    if step_mm is None:
        step_mm = 0.5 * float(np.min(np.asarray(vox_spacing, dtype=float)))
    step_mm = max(float(step_mm), 0.2)
    max_steps = max(int(float(max_dist_mm) / step_mm), 1)

    def _inside(pt_world):
        kji = _world_to_voxel(pt_world, vox_spacing, offset, direction=direction)
        k, j, i = (int(round(v)) for v in kji)
        if 0 <= k < mask.shape[0] and 0 <= j < mask.shape[1] and 0 <= i < mask.shape[2]:
            return bool(mask[k, j, i])
        return False

    cur = p.copy()
    last_inside = cur.copy() if _inside(cur) else None
    first_outside = None

    for _ in range(max_steps):
        nxt = cur + d * step_mm
        if _inside(nxt):
            last_inside = nxt.copy()
            cur = nxt
        else:
            first_outside = nxt.copy()
            break

    if return_outside:
        if first_outside is None:
            first_outside = cur.copy() if not _inside(cur) else (cur + d * step_mm)
        m = max(float(outside_margin_mm), 0.0)
        return first_outside + d * m

    return last_inside if last_inside is not None else cur

def _recenter_point_in_terminal_plane(
    mask,
    anchor_world,
    tangent_world,
    vox_spacing,
    offset,
    direction=None,
    radius_mm=28.0,
    slab_mm=8.0,
    min_voxels=24,
):
    """Recenter a point using a tangent-orthogonal slab and medial-axis weighting."""
    if radius_mm is None:
        radius_mm = 28.0
    if slab_mm is None:
        slab_mm = 8.0

    t = np.asarray(tangent_world, dtype=float)
    t_norm = np.linalg.norm(t)
    if t_norm < 1e-9:
        return np.asarray(anchor_world, dtype=float)
    t = t / t_norm

    p = np.asarray(anchor_world, dtype=float)

    spacing_kji = np.array([vox_spacing[2], vox_spacing[1], vox_spacing[0]], dtype=float)
    p_kji = _world_to_voxel(p, vox_spacing, offset, direction=direction)
    p_kji = np.round(p_kji).astype(int)

    r_vox = np.maximum(np.ceil(float(radius_mm) / np.maximum(spacing_kji, 1e-6)).astype(int), 2)
    k0, j0, i0 = p_kji
    kmin, kmax = max(0, k0 - r_vox[0]), min(mask.shape[0], k0 + r_vox[0] + 1)
    jmin, jmax = max(0, j0 - r_vox[1]), min(mask.shape[1], j0 + r_vox[1] + 1)
    imin, imax = max(0, i0 - r_vox[2]), min(mask.shape[2], i0 + r_vox[2] + 1)

    local_mask = mask[kmin:kmax, jmin:jmax, imin:imax]
    if not np.any(local_mask):
        return p

    local_vox = np.argwhere(local_mask)
    full_vox = local_vox + np.array([kmin, jmin, imin], dtype=int)
    world_pts = _voxel_to_world(full_vox, vox_spacing, offset, direction=direction)

    rel = world_pts - p
    axial = rel @ t
    radial_vec = rel - np.outer(axial, t)
    radial = np.linalg.norm(radial_vec, axis=1)

    slab_half = 0.5 * float(slab_mm)
    keep = (np.abs(axial) <= slab_half) & (radial <= float(radius_mm))
    if np.count_nonzero(keep) < min_voxels:
        return p

    sel_local_vox = local_vox[keep]
    sel_world = world_pts[keep]
    sel_axial = axial[keep]

    select_mask = np.zeros_like(local_mask, dtype=bool)
    select_mask[sel_local_vox[:, 0], sel_local_vox[:, 1], sel_local_vox[:, 2]] = True

    labeled, nlab = ndimage.label(select_mask)
    if nlab < 1:
        return p

    p_local = p_kji - np.array([kmin, jmin, imin], dtype=int)
    nearest_idx = int(np.argmin(np.sum((sel_local_vox - p_local[None, :]) ** 2, axis=1)))
    seed = sel_local_vox[nearest_idx]
    comp_id = labeled[seed[0], seed[1], seed[2]]
    if comp_id <= 0:
        return p

    comp_keep = labeled[sel_local_vox[:, 0], sel_local_vox[:, 1], sel_local_vox[:, 2]] == comp_id
    cap_world = sel_world[comp_keep]
    cap_axial = sel_axial[comp_keep]
    cap_local = sel_local_vox[comp_keep]

    if len(cap_world) < min_voxels:
        return p

    dt_local = ndimage.distance_transform_edt(local_mask, sampling=spacing_kji)
    dt_w = dt_local[cap_local[:, 0], cap_local[:, 1], cap_local[:, 2]]
    dt_w = np.maximum(dt_w, 1e-3) ** 2
    axial_w = np.exp(-(cap_axial / max(slab_half, 1e-6)) ** 2)
    w = dt_w * axial_w

    return np.average(cap_world, axis=0, weights=w)

def generate_cl(
    mhd_iso_img,
    n_points=1200,
    spacing=None,
    downsample=1,
    cap_radius_mm=None,
    cap_slab_mm=None,
    endpoint_recenter=True,
    endpoint_recenter_radius_mm=30.0,
    endpoint_recenter_slab_mm=8.0,
    endpoint_axis_mode="auto",
    start_outside_margin_mm=2.0,  # 1-5 mm
    return_endpoint_hint=False,
):
    """Centerline via DT-based geodesic shortest path (scipy + skimage only)."""
    
    vol = sitk.GetArrayFromImage(mhd_iso_img)
    vox_spacing = np.array(mhd_iso_img.GetSpacing(), dtype=float)
    offset = np.array(mhd_iso_img.GetOrigin(), dtype=float)
    direction = np.array(mhd_iso_img.GetDirection(), dtype=float).reshape(3, 3)
    mask = vol > 0

    if cap_radius_mm is None or cap_slab_mm is None:
        vox = np.argwhere(mask)
        if len(vox) > 0:
            spacing_kji = np.array([vox_spacing[2], vox_spacing[1], vox_spacing[0]], dtype=float)
            ext_kji_mm = (vox.max(axis=0) - vox.min(axis=0) + 1) * spacing_kji
            cross_section_scale_mm = float(np.mean(np.sort(ext_kji_mm)[:2]))
        else:
            cross_section_scale_mm = 30.0

        auto_cap_radius_mm = float(np.clip(0.90 * cross_section_scale_mm, 10.0, 60.0))
        auto_cap_slab_mm = float(np.clip(0.30 * cross_section_scale_mm, 4.0, 30.0))

        if cap_radius_mm is None:
            cap_radius_mm = auto_cap_radius_mm
        if cap_slab_mm is None:
            cap_slab_mm = auto_cap_slab_mm

    cap_radius_mm = float(cap_radius_mm)
    cap_slab_mm = float(cap_slab_mm)

    ds = max(int(downsample), 1)
    mask_ds = mask[::ds, ::ds, ::ds]

    struct = ndimage.generate_binary_structure(3, 1)
    mask_ds = ndimage.binary_closing(mask_ds, structure=struct, iterations=2)
    mask_ds = ndimage.binary_fill_holes(mask_ds)

    skel = skeletonize(mask_ds).astype(np.uint8)
    skel_ijk = np.argwhere(skel > 0)
    if len(skel_ijk) < 10:
        raise ValueError(f"Skeleton too small ({len(skel_ijk)} voxels).")

    spacing_ds = vox_spacing * ds
    skel_world = _voxel_to_world(skel_ijk, spacing_ds, offset, direction=direction)
    adj = _skeleton_to_graph(skel_ijk)

    u, v, largest_nodes, largest_endpoints = _select_endpoint_pair(adj, skel_world)

    if endpoint_axis_mode in {"x", "y", "z"} and len(largest_endpoints) >= 2:
        axis_lut = {"x": 0, "y": 1, "z": 2}
        dominant_axis = axis_lut[endpoint_axis_mode]
        endpoint_world = skel_world[largest_endpoints]
        u_axis = int(largest_endpoints[np.argmin(endpoint_world[:, dominant_axis])])
        v_axis = int(largest_endpoints[np.argmax(endpoint_world[:, dominant_axis])])
        dist_u_axis = _bfs_dist(adj, u_axis)
        if dist_u_axis[v_axis] >= 0:
            u, v = u_axis, v_axis

    print(f"[centerline-endpoints] aligned and identified"
        # f"[centerline-endpoints] u={u} xyz={np.round(skel_world[u], 2)}  "
        # f"v={v} xyz={np.round(skel_world[v], 2)}"
    )

    path_indices = _bfs_path(adj, u, v)
    raw_skel_path = skel_world[path_indices]
    if raw_skel_path.shape[0] < 3:
        raise ValueError("Crude skeleton path too short.")

    endpoint_hint = {
        "crude_start_world": np.asarray(raw_skel_path[0], dtype=float),
        "crude_end_world": np.asarray(raw_skel_path[-1], dtype=float),
        "crude_start_z": float(raw_skel_path[0, 2]),
        "crude_end_z": float(raw_skel_path[-1, 2]),
    }

    n_dir = min(5, raw_skel_path.shape[0] - 1)
    start_inward = _safe_normalize((raw_skel_path[n_dir] - raw_skel_path[0])[None, :])[0]
    end_inward = _safe_normalize((raw_skel_path[-1 - n_dir] - raw_skel_path[-1])[None, :])[0]

    start_seed = np.asarray(skel_world[u], dtype=float)
    end_seed = np.asarray(skel_world[v], dtype=float)
    start_edge = _push_point_to_mask_boundary(mask, start_seed, -start_inward, vox_spacing, offset, direction=direction)
    end_edge = _push_point_to_mask_boundary(mask, end_seed, -end_inward, vox_spacing, offset, direction=direction)

    if endpoint_recenter:
        start_center = _recenter_point_in_terminal_plane(
            mask=mask,
            anchor_world=start_edge,
            tangent_world=-start_inward,
            vox_spacing=vox_spacing,
            offset=offset,
            direction=direction,
            radius_mm=endpoint_recenter_radius_mm,
            slab_mm=endpoint_recenter_slab_mm,
        )
        end_center = _recenter_point_in_terminal_plane(
            mask=mask,
            anchor_world=end_edge,
            tangent_world=-end_inward,
            vox_spacing=vox_spacing,
            offset=offset,
            direction=direction,
            radius_mm=endpoint_recenter_radius_mm,
            slab_mm=endpoint_recenter_slab_mm,
        )
    else:
        start_center = start_edge
        end_center = end_edge

    # Force start point to be outside mask by a small margin.
    start_center = _push_point_to_mask_boundary(
        mask,
        start_center,
        -start_inward,
        vox_spacing,
        offset,
        direction=direction,
        outside_margin_mm=start_outside_margin_mm,
        return_outside=True,
    )

    spacing_kji = np.array([vox_spacing[2], vox_spacing[1], vox_spacing[0]], dtype=float)
    D = ndimage.distance_transform_edt(mask, sampling=spacing_kji)

    seedA = np.round(_world_to_voxel(start_center, vox_spacing, offset, direction=direction)).astype(int)
    seedB = np.round(_world_to_voxel(end_center, vox_spacing, offset, direction=direction)).astype(int)
    seedA = np.array([
        np.clip(seedA[0], 0, mask.shape[0] - 1),
        np.clip(seedA[1], 0, mask.shape[1] - 1),
        np.clip(seedA[2], 0, mask.shape[2] - 1),
    ], dtype=int)
    seedB = np.array([
        np.clip(seedB[0], 0, mask.shape[0] - 1),
        np.clip(seedB[1], 0, mask.shape[1] - 1),
        np.clip(seedB[2], 0, mask.shape[2] - 1),
    ], dtype=int)

    def _snap_seed_to_mask_local(seed_kji, radius_vox=6):
        if mask[seed_kji[0], seed_kji[1], seed_kji[2]]:
            return seed_kji.astype(int)
        k0, j0, i0 = seed_kji.astype(int)
        k1 = max(0, k0 - radius_vox)
        k2 = min(mask.shape[0], k0 + radius_vox + 1)
        j1 = max(0, j0 - radius_vox)
        j2 = min(mask.shape[1], j0 + radius_vox + 1)
        i1 = max(0, i0 - radius_vox)
        i2 = min(mask.shape[2], i0 + radius_vox + 1)
        blk = mask[k1:k2, j1:j2, i1:i2]
        pts = np.argwhere(blk)
        if len(pts) == 0:
            pts_all = np.argwhere(mask)
            if len(pts_all) == 0:
                return seed_kji.astype(int)
            d2 = np.sum((pts_all - seed_kji[None, :]) ** 2, axis=1)
            return pts_all[int(np.argmin(d2))].astype(int)
        pts_abs = pts + np.array([k1, j1, i1], dtype=int)
        d2 = np.sum((pts_abs - seed_kji[None, :]) ** 2, axis=1)
        return pts_abs[int(np.argmin(d2))].astype(int)

    seedA = _snap_seed_to_mask_local(seedA)
    seedB = _snap_seed_to_mask_local(seedB)

    inside_cost = 1.0 / np.maximum(D, 1e-3) ** 2
    cost = np.full(mask.shape, 1e6, dtype=float)
    cost[mask] = inside_cost[mask]

    path_idx, _ = route_through_array(
        cost,
        start=tuple(seedA.tolist()),
        end=tuple(seedB.tolist()),
        fully_connected=True,
        geometric=True,
    )
    path_kji = np.asarray(path_idx, dtype=float)

    raw_path = _voxel_to_world(path_kji, vox_spacing, offset, direction=direction)
    raw_path[0] = np.asarray(start_center, dtype=float)
    raw_path[-1] = np.asarray(end_center, dtype=float)

    seg = np.linalg.norm(np.diff(raw_path, axis=0), axis=1)
    keep = np.concatenate(([True], seg > 1e-9))
    raw_path = raw_path[keep]

    if spacing is not None:
        total_len = float(np.sum(np.linalg.norm(np.diff(raw_path, axis=0), axis=1)))
        n_points = max(int(total_len / float(spacing)) + 1, 4)
    n_points = max(int(n_points), 4)

    d = np.linalg.norm(np.diff(raw_path, axis=0), axis=1)
    s = np.concatenate(([0.0], np.cumsum(d)))
    if s[-1] < 1e-9:
        raise ValueError("Centerline length is zero after endpoint centering.")
    s_new = np.linspace(0.0, float(s[-1]), n_points)

    cline_pos = np.column_stack([
        np.interp(s_new, s, raw_path[:, 0]),
        np.interp(s_new, s, raw_path[:, 1]),
        np.interp(s_new, s, raw_path[:, 2]),
    ])
    cline_div = _safe_normalize(np.gradient(cline_pos, axis=0))

    if return_endpoint_hint:
        return cline_pos, cline_div, np.eye(4), 0.0, endpoint_hint
    return cline_pos, cline_div, np.eye(4), 0.0

def save_centerline_csv(centerline_pts, centerline_tangents, csv_path):
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(
        {
            "Position x (global)": centerline_pts[:, 0],
            "Position y": centerline_pts[:, 1],
            "Position z": centerline_pts[:, 2],
            "Tangent x": centerline_tangents[:, 0],
            "Tangent y": centerline_tangents[:, 1],
            "Tangent z": centerline_tangents[:, 2],
        }
    )
    df.to_csv(csv_path, index=False)
    return csv_path