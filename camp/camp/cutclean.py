# Cutting and cleaning functions 

import numpy as np
import trimesh
from collections import deque
import pyacvd
from scipy.spatial import cKDTree

def _full_adjacency(mesh):
    n_faces = len(mesh.faces)
    full = [[] for _ in range(n_faces)]
    for a, b in mesh.face_adjacency:
        full[int(a)].append(int(b))
        full[int(b)].append(int(a))
    return full

def _vertex_face_adjacency(mesh):
    v2f = [[] for _ in range(len(mesh.vertices))]
    for fi, tri in enumerate(np.asarray(mesh.faces, dtype=int)):
        v2f[int(tri[0])].append(int(fi))
        v2f[int(tri[1])].append(int(fi))
        v2f[int(tri[2])].append(int(fi))
    return v2f

def _nearest_vertex_index(mesh, point_world, allowed_vertices=None):
    p = np.asarray(point_world, dtype=float)
    v = np.asarray(mesh.vertices, dtype=float)
    if allowed_vertices is not None:
        allowed_vertices = np.asarray(allowed_vertices, dtype=bool)
        if allowed_vertices.shape[0] != v.shape[0] or not np.any(allowed_vertices):
            return None
        idx = np.flatnonzero(allowed_vertices)
        d2 = np.sum((v[idx] - p[None, :]) ** 2, axis=1)
        return int(idx[int(np.argmin(d2))])
    d2 = np.sum((v - p[None, :]) ** 2, axis=1)
    return int(np.argmin(d2))

def _faces_all_vertices_in_z_band(mesh, zlo, zhi):
    verts_z = np.asarray(mesh.vertices, dtype=float)[:, 2]
    faces = np.asarray(mesh.faces, dtype=int)
    zf = verts_z[faces]
    return np.all((zf >= float(zlo)) & (zf <= float(zhi)), axis=1)

def _pick_seed_face_from_vertex(mesh, vertex_idx, vertex_to_faces, eligible_faces, point_world):
    if vertex_idx is None:
        return None
    candidates = [fi for fi in vertex_to_faces[int(vertex_idx)] if bool(eligible_faces[int(fi)])]
    if not candidates:
        return None
    if len(candidates) == 1:
        return int(candidates[0])

    cent = np.asarray(mesh.triangles_center, dtype=float)
    p = np.asarray(point_world, dtype=float)
    d2 = [float(np.sum((cent[int(fi)] - p) ** 2)) for fi in candidates]
    return int(candidates[int(np.argmin(d2))])

def _geodesic_patch_from_seed(seed_face, eligible_faces, full_neigh):
    if seed_face is None or not bool(eligible_faces[int(seed_face)]):
        return np.zeros(len(eligible_faces), dtype=bool)

    keep = np.zeros(len(eligible_faces), dtype=bool)
    seen = {int(seed_face)}
    q = deque([int(seed_face)])
    keep[int(seed_face)] = True
    while q:
        f = q.popleft()
        for nb in full_neigh[int(f)]:
            nbi = int(nb)
            if nbi in seen:
                continue
            if not bool(eligible_faces[nbi]):
                continue
            seen.add(nbi)
            keep[nbi] = True
            q.append(nbi)
    return keep

def _nearest_eligible_face(mesh, point_world, eligible_faces):
    idx = np.flatnonzero(np.asarray(eligible_faces, dtype=bool))
    if len(idx) == 0:
        return None
    cent = np.asarray(mesh.triangles_center, dtype=float)[idx]
    p = np.asarray(point_world, dtype=float)
    d2 = np.sum((cent - p[None, :]) ** 2, axis=1)
    return int(idx[int(np.argmin(d2))])

def remove_terminal_caps_with_buffer(mesh, centerline_pts, mhd_iso_img, buffer_mm=5.0, endpoint_hint=None):
    """Remove terminal caps using geodesic patches constrained by z-bands.

    For each endpoint, build [z - buffer_mm, z + buffer_mm], find the closest
    mesh vertex to that endpoint, choose an eligible incident face, and flood
    geodesically over adjacent faces whose three vertices all stay inside the
    same z-band.
    """
    mesh = mesh.copy()
    if len(centerline_pts) < 2:
        print("Centerline has <2 points; skipping terminal cap removal.")
        return mesh
    buffer_mm = float(max(buffer_mm, 0.0))

    start_pt = np.asarray(centerline_pts[0], dtype=float)
    end_pt = np.asarray(centerline_pts[-1], dtype=float)
    if isinstance(endpoint_hint, dict):
        if "crude_start_world" in endpoint_hint:
            start_pt = np.asarray(endpoint_hint["crude_start_world"], dtype=float)
        if "crude_end_world" in endpoint_hint:
            end_pt = np.asarray(endpoint_hint["crude_end_world"], dtype=float)

    full_neigh = _full_adjacency(mesh)
    vertex_to_faces = _vertex_face_adjacency(mesh)
    to_remove = np.zeros(len(mesh.faces), dtype=bool)

    for point, label in [(start_pt, "start"), (end_pt, "end")]:
        z_center = float(point[2])
        zlo, zhi = z_center - buffer_mm, z_center + buffer_mm
        eligible = _faces_all_vertices_in_z_band(mesh, zlo, zhi)
        n_eligible = int(np.count_nonzero(eligible))
        if n_eligible == 0:
            print(f"[{label}] no faces fully within z-band [{zlo:.1f}, {zhi:.1f}]; cap NOT removed.")
            continue

        eligible_vertices = np.zeros(len(mesh.vertices), dtype=bool)
        eligible_vertices[np.asarray(mesh.faces, dtype=int)[eligible].ravel()] = True

        seed_vertex = _nearest_vertex_index(mesh, point, allowed_vertices=eligible_vertices)
        seed_face = _pick_seed_face_from_vertex(
            mesh,
            vertex_idx=seed_vertex,
            vertex_to_faces=vertex_to_faces,
            eligible_faces=eligible,
            point_world=point,
        )
        if seed_face is None:
            seed_face = _nearest_eligible_face(mesh, point, eligible)
            if seed_face is None:
                print(f"[{label}] could not find in-band seed face in z-band "
                      f"[{zlo:.1f}, {zhi:.1f}]; cap NOT removed.")
                continue
            print(f"[{label}] fallback seed face selected by nearest in-band centroid: {seed_face}.")

        patch = _geodesic_patch_from_seed(seed_face, eligible, full_neigh)
        n_patch = int(np.count_nonzero(patch))
        if n_patch == 0:
            print(f"[{label}] seed face found but produced empty patch; cap NOT removed.")
            continue

        to_remove |= patch
        print(f"[{label}] z={z_center:.1f}mm, band=[{zlo:.1f}, {zhi:.1f}], "
              f"seed_vertex={seed_vertex}, seed_face={seed_face}, removed={n_patch} faces.")

    n_remove = int(np.count_nonzero(to_remove))
    if n_remove == 0:
        print("No cap faces selected; mesh unchanged.")
        return mesh

    mesh.update_faces(~to_remove)
    mesh.remove_unreferenced_vertices()
    print(f"Removed {n_remove} cap faces total (buffer={buffer_mm:.1f} mm in z).")
    return mesh

def keep_largest_components(mesh, n_parts=1):
    components = mesh.split(only_watertight=False)
    if len(components) <= n_parts:
        return trimesh.util.concatenate(components) if len(components) > 1 else components[0]
    largest = sorted(components, key=lambda m: len(m.faces), reverse=True)[:n_parts]
    return trimesh.util.concatenate(largest)

def keep_component_guided_by_centerline(mesh, centerline_pts):
    """Keep the component best described by the centerline trajectory.

    Selection is based on:
      1) how much of the centerline z-extent lies within the component z-range,
      2) whether both centerline endpoints lie in that z-range,
      3) mean nearest-vertex distance from centerline to component,
      4) face count as a final tie-breaker.
    """
    components = list(mesh.split(only_watertight=False))
    if len(components) == 0:
        return mesh
    if len(components) == 1:
        return components[0]

    cl = np.asarray(centerline_pts, dtype=float)
    if cl.ndim != 2 or cl.shape[0] < 2 or cl.shape[1] != 3:
        return keep_largest_components(mesh, n_parts=1)

    cl_z = cl[:, 2]
    z0 = float(cl[0, 2])
    z1 = float(cl[-1, 2])

    scored = []
    for i, comp in enumerate(components):
        v = np.asarray(comp.vertices, dtype=float)
        if v.shape[0] == 0:
            continue

        zmin = float(np.min(v[:, 2]))
        zmax = float(np.max(v[:, 2]))
        z_hits = int(np.count_nonzero((cl_z >= zmin) & (cl_z <= zmax)))
        endpoint_hits = int(zmin <= z0 <= zmax) + int(zmin <= z1 <= zmax)

        tree = cKDTree(v)
        dists, _ = tree.query(cl, k=1)
        mean_dist = float(np.mean(dists))

        score = (z_hits, endpoint_hits, -mean_dist, int(len(comp.faces)))
        scored.append((score, i, comp, zmin, zmax, mean_dist))

    if not scored:
        return keep_largest_components(mesh, n_parts=1)

    scored.sort(key=lambda x: x[0], reverse=True)
    best = scored[0]
    best_comp = best[2]

    print(
        "[cut-select] keeping component index="
        f"{best[1]} faces={len(best_comp.faces)} z=[{best[3]:.1f}, {best[4]:.1f}] "
        f"mean_cl_to_surface={best[5]:.2f} mm"
    )
    return best_comp

def trimesh_to_pyvista(tm_mesh):
    import pyvista as pv

    faces = np.hstack([
        np.full((len(tm_mesh.faces), 1), 3, dtype=np.int64),
        np.asarray(tm_mesh.faces, dtype=np.int64),
    ]).ravel()
    return pv.PolyData(np.asarray(tm_mesh.vertices), faces).triangulate()

def pyvista_to_trimesh(pd):
    f = pd.faces.reshape(-1, 4)[:, 1:].astype(np.int64)
    v = np.asarray(pd.points)
    out = trimesh.Trimesh(vertices=v, faces=f, process=False)
    out.remove_unreferenced_vertices()
    out.merge_vertices(digits_vertex=6)
    out.remove_unreferenced_vertices()
    return trimesh.Trimesh(vertices=out.vertices, faces=out.faces, process=False)

def acvd_uniform_remesh(mesh, target_faces=50_000):
    pmesh = trimesh_to_pyvista(mesh)
    target_vertices = int(target_faces // 2)
    clus = pyacvd.Clustering(pmesh)
    clus.subdivide(int(2))
    clus.cluster(target_vertices)
    remeshed = clus.create_mesh().triangulate()
    out = pyvista_to_trimesh(remeshed)
    return out

def cut_clean(precut_mesh, centerline_pts, buffer_mm, mhd_iso_img, target_faces, endpoint_hint=None):
    open_mesh = remove_terminal_caps_with_buffer(
        precut_mesh,
        centerline_pts,
        mhd_iso_img=mhd_iso_img,
        buffer_mm=buffer_mm,
        endpoint_hint=endpoint_hint,
    )
    cut = keep_component_guided_by_centerline(open_mesh, centerline_pts)
    euler = int(cut.euler_number)
    print(f"[cut] final euler_number={euler} ({'OPEN' if euler == 0 else 'NOT OPEN'})")

    cleaned = acvd_uniform_remesh(cut, target_faces=target_faces)

    return cleaned