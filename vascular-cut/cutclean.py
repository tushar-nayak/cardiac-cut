# Cutting and cleaning functions 

import numpy as np
import trimesh
from collections import deque
import pyacvd

def _extend_centerline_endpoints(centerline_pts, mhd_iso_img, overshoot_mm=20.0, tangent_span=8):
    """Extend both terminal points along their terminal tangent until they are
    OUTSIDE the MHD volume (so a ray cast back inward must cross the cap).

    Returns (start_outside, start_dir_inward, end_outside, end_dir_inward) where
    each *_dir_inward points from the outside point back toward the mesh.
    """
    cl = np.asarray(centerline_pts, float)

    def _inside(p):
        # True if world point p maps inside the MHD foreground mask.
        try:
            idx = mhd_iso_img.TransformPhysicalPointToIndex(tuple(float(v) for v in p))
        except Exception:
            return False
        sz = mhd_iso_img.GetSize()
        if not all(0 <= idx[k] < sz[k] for k in range(3)):
            return False
        return mhd_iso_img.GetPixel(idx) > 0

    spacing = float(min(mhd_iso_img.GetSpacing()))
    step = max(0.5 * spacing, 0.25)
    max_steps = int(np.ceil(overshoot_mm / step)) + 4

    def _push(p_end, tangent_outward):
        d = tangent_outward / (np.linalg.norm(tangent_outward) + 1e-12)
        p = np.array(p_end, float)
        # March outward until we are clearly outside the volume.
        for _ in range(max_steps):
            if not _inside(p):
                p = p + d * overshoot_mm  # clear the boundary with margin
                break
            p = p + d * step
        return p, -d  # inward direction = back toward the mesh

    n = min(tangent_span, len(cl) - 1)
    start_tan_out = cl[0] - cl[n]           # points outward at the start
    end_tan_out = cl[-1] - cl[-1 - n]       # points outward at the end

    s_out, s_in = _push(cl[0], start_tan_out)
    e_out, e_in = _push(cl[-1], end_tan_out)
    return s_out, s_in, e_out, e_in

def _geodesic_faces_within_z(mesh, seed_face, z_center, z_half, full_neigh):
    """Faces geodesically connected to seed_face whose centroid z is within
    [z_center - z_half, z_center + z_half]. Connectivity is on the surface;
    the z-band only constrains which faces are eligible."""
    centroids = np.asarray(mesh.triangles_center)
    zlo, zhi = z_center - z_half, z_center + z_half
    in_band = (centroids[:, 2] >= zlo) & (centroids[:, 2] <= zhi)

    if not in_band[seed_face]:
        # Walk from the seed to the nearest in-band face along the surface.
        seen = {seed_face}
        q = deque([seed_face])
        seed_face2 = None
        while q:
            f = q.popleft()
            if in_band[f]:
                seed_face2 = f
                break
            for nb in full_neigh[f]:
                if nb not in seen:
                    seen.add(nb)
                    q.append(nb)
        if seed_face2 is None:
            return np.zeros(len(mesh.faces), dtype=bool)
        seed_face = seed_face2

    keep = np.zeros(len(mesh.faces), dtype=bool)
    seen = {seed_face}
    q = deque([seed_face])
    keep[seed_face] = True
    while q:
        f = q.popleft()
        for nb in full_neigh[f]:
            if nb not in seen and in_band[nb]:
                seen.add(nb)
                keep[nb] = True
                q.append(nb)
    return keep

def _ray_first_hit(mesh, origin, direction):
    """First face index where a ray from origin along direction hits the mesh."""
    locs, ray_idx, tri_idx = mesh.ray.intersects_location(
        ray_origins=np.asarray([origin], float),
        ray_directions=np.asarray([direction], float),
        multiple_hits=False,
    )
    if len(locs) == 0:
        return None, None
    j = int(np.argmin(np.linalg.norm(locs - np.asarray(origin, float), axis=1)))
    return int(tri_idx[j]), locs[j]

def _full_adjacency(mesh):
    n_faces = len(mesh.faces)
    full = [[] for _ in range(n_faces)]
    for a, b in mesh.face_adjacency:
        full[int(a)].append(int(b))
        full[int(b)].append(int(a))
    return full

def remove_terminal_caps_with_buffer(mesh, centerline_pts, mhd_iso_img, buffer_mm=5.0):
    """Remove the two synthetic flat caps by ray-casting the (overshot) centerline
    terminals back toward the mesh, then deleting faces within +/- buffer_mm in
    GLOBAL z of each hit that are geodesically connected to it."""
    mesh = mesh.copy()
    if len(centerline_pts) < 2:
        print("Centerline has <2 points; skipping terminal cap removal.")
        return mesh
    buffer_mm = float(max(buffer_mm, 0.0))

    s_out, s_in, e_out, e_in = _extend_centerline_endpoints(centerline_pts, mhd_iso_img)
    full_neigh = _full_adjacency(mesh)
    to_remove = np.zeros(len(mesh.faces), dtype=bool)

    for origin, direction, label in [(s_out, s_in, "start"), (e_out, e_in, "end")]:
        hit_face, hit_pt = _ray_first_hit(mesh, origin, direction)
        if hit_face is None:
            print(f"[{label}] ray never re-entered the mesh; cap NOT removed "
                  f"(will surface in the euler check).")
            continue
        band = _geodesic_faces_within_z(
            mesh, hit_face, z_center=float(hit_pt[2]), z_half=buffer_mm, full_neigh=full_neigh,
        )
        to_remove |= band
        print(f"[{label}] cap hit at z={hit_pt[2]:.1f}mm; "
              f"removed {int(band.sum())} geodesically-connected faces "
              f"(z in [{hit_pt[2]-buffer_mm:.1f}, {hit_pt[2]+buffer_mm:.1f}]).")

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

"""Subdivisions applied before clustering. ACVD clusters vertices, so it can
never make more clusters than the mesh has vertices; subdividing first supplies
them. Two levels is ~16x and is what every published run used, so it stays the
starting point, raising it for a mesh that does not need it would repartition
that mesh and move any curvature measured on it."""
SUBDIVIDE_BASE = 2
SUBDIVIDE_MAX = 6
# Vertices per cluster needed before clustering is worth attempting. Measured:
# 39,284 vertices clustered to 25,000 fine (1.6x) but silently capped at 50,000.
SUBDIVIDE_HEADROOM = 1.5


def acvd_uniform_remesh(mesh, target_faces=50_000, verify=True):
    """Uniform ACVD remesh to approximately target_faces.

    A small mesh cannot always reach a high target: a 2,552-vertex region
    subdivided 16x offers 39,284 vertices, and asking for 50,000 clusters
    returns 77,536 faces rather than 100,000, with no error. So deepen the
    subdivision when the ratio demands it, and check the result afterwards
    instead of assuming it.

    Meshes that already had enough vertices keep SUBDIVIDE_BASE and are
    unchanged, so this does not perturb anything that worked before.
    """
    pmesh = trimesh_to_pyvista(mesh)
    target_vertices = int(target_faces // 2)

    depth = SUBDIVIDE_BASE
    available = pmesh.n_points * (4 ** depth)
    while available < target_vertices * SUBDIVIDE_HEADROOM and depth < SUBDIVIDE_MAX:
        depth += 1
        available *= 4

    clus = pyacvd.Clustering(pmesh)
    clus.subdivide(int(depth))
    clus.cluster(target_vertices)
    remeshed = clus.create_mesh().triangulate()
    out = pyvista_to_trimesh(remeshed)

    if verify and len(out.faces) < 0.9 * target_faces:
        raise ValueError(
            f"ACVD reached {len(out.faces):,} faces, not the {target_faces:,} "
            f"requested: {pmesh.n_points:,} input vertices subdivided {depth}x "
            f"give {available:,}, short of the {target_vertices:,} clusters "
            f"needed. Lower target_faces, or raise SUBDIVIDE_MAX.")
    return out

def cut_clean(precut_mesh, centerline_pts, buffer_mm, mhd_iso_img, target_faces):
    open_mesh = remove_terminal_caps_with_buffer(
        precut_mesh, centerline_pts, mhd_iso_img=mhd_iso_img, buffer_mm=buffer_mm,
    )
    cut = keep_largest_components(open_mesh, n_parts=1) 
    euler = int(cut.euler_number)
    print(f"[cut] final euler_number={euler} ({'OPEN' if euler == 0 else 'NOT OPEN'})")

    cleaned = acvd_uniform_remesh(cut, target_faces=target_faces)

    return cleaned