"""Uniform surface remeshing (ACVD) — the single preprocessing step before CAMP.

CAMP's patch statistics assume a roughly UNIFORM triangulation, so every mesh is
remeshed to the same target face count before partitioning.  This is a clean,
self-contained extract of the project's ``acvd_uniform_remesh`` (no SimpleITK,
no centerline code) so collaborators only need ``pyvista`` + ``pyacvd``.
"""
import numpy as np
import trimesh


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


def acvd_uniform_remesh(mesh, target_faces=100_000):
    """Uniformly remesh ``mesh`` to approximately ``target_faces`` triangles.

    Uses pyacvd (Approximated Centroidal Voronoi Diagrams) for an isotropic,
    near-equilateral retriangulation.  Topology (number of holes / branches) is
    preserved; only the triangulation is resampled.
    """
    import pyacvd
    pmesh = trimesh_to_pyvista(mesh)
    target_vertices = int(target_faces // 2)          # closed-ish surface: F ~ 2V
    clus = pyacvd.Clustering(pmesh)
    clus.subdivide(int(2))                            # densify first so clusters are well-formed
    clus.cluster(target_vertices)
    remeshed = clus.create_mesh().triangulate()
    return pyvista_to_trimesh(remeshed)
