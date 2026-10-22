# Cardiac geometry: split a closed left-ventricle myocardium into its inner
# (endocardial) and outer (epicardial) walls, and cut the base off with a plane.
#
# The cardiac counterpart of branch_cut.py. branch_cut cuts tubes across a
# centerline; a ventricle has no centerline, so it is handled here instead:
#
#   1. The largest connected component is the shell; detached islands are counted.
#   2. Inner and outer seeds come from the face normals around an estimated
#      interior point, and every face is assigned to the nearer seed along the
#      surface (a geodesic Voronoi seam). No labels, no hand-placed plane.
#   3. A plane is fitted to the basal rim -- the faces where the outer wall folds
#      over into the cavity -- so a tilted base is followed, and the whole shell
#      is cut there. With the rim gone the walls only touch through holes in the
#      wall, so two pieces means a clean split and the pieces ARE the walls.
#   4. Where holes still join them, the seam from step 2 cuts each hole's tube.
#      The holes are left OPEN by default: CAMP does not need a closed surface,
#      and a synthetic patch adds curvature the anatomy does not have
#      (CAMP_INTEGRATION_GUIDE.md, Q1). Patching is available, and recorded.
#
# Nothing here imports Qt or trame. HeartSession is the headless layer the web app
# drives, as vascular_api.CutSession is for vessels.

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import trimesh
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra


class HeartError(RuntimeError):
    pass


# ── seeds and the geodesic seam ──────────────────────────────────────────────
def _dominant_region(mask: np.ndarray, adjacency: np.ndarray, n_faces: int,
                     face_area: np.ndarray) -> tuple[np.ndarray, float]:
    """Return the connected region with the most surface area and its share."""
    pairs = adjacency[mask[adjacency[:, 0]] & mask[adjacency[:, 1]]]
    graph = coo_matrix(
        (np.ones(2 * len(pairs)),
         (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
        shape=(n_faces, n_faces),
    ).tocsr()
    n_components, labels = connected_components(graph, directed=False)
    ids = np.flatnonzero(mask)
    if not len(ids):
        return np.zeros(n_faces, dtype=bool), 0.0
    areas = np.bincount(labels[ids], weights=face_area[ids], minlength=n_components)
    largest = int(np.argmax(areas))
    return mask & (labels == largest), float(areas[largest] / areas.sum())


def _candidate_centers(mesh: trimesh.Trimesh):
    center = np.asarray(mesh.center_mass, dtype=float)
    span = np.ptp(mesh.vertices, axis=0)
    for offset in itertools.product(np.linspace(-0.5, 0.5, 5), repeat=3):
        yield center + np.asarray(offset) * span


def _choose_center(mesh: trimesh.Trimesh, seed_threshold: float = 0.08):
    centers = mesh.triangles_center
    normals = mesh.face_normals
    areas = mesh.area_faces
    adjacency = mesh.face_adjacency
    edge_vertices = mesh.face_adjacency_edges
    edge_lengths = np.linalg.norm(
        mesh.vertices[edge_vertices[:, 0]] - mesh.vertices[edge_vertices[:, 1]], axis=1
    )

    # Cheap pass over the candidate grid; rank by a short separator and reasonable
    # area balance. Only the best candidates need the more expensive connectivity test.
    ranked = []
    for point in _candidate_centers(mesh):
        score = np.einsum("ij,ij->i", normals, centers - point)
        inner = score < 0
        inner_area = float(areas[inner].sum())
        outer_area = float(areas[~inner].sum())
        frac = inner_area / (inner_area + outer_area)
        boundary = float(edge_lengths[inner[adjacency[:, 0]] != inner[adjacency[:, 1]]].sum())
        cut_cost = boundary / np.sqrt(max(inner_area * outer_area, 1e-12))
        balance_penalty = max(0.0, 0.25 - frac) + max(0.0, frac - 0.50)
        ranked.append((cut_cost + 2.0 * balance_penalty, frac, cut_cost, point, score))
    ranked.sort(key=lambda x: x[0])

    best = None
    for _, frac, cut_cost, point, score in ranked[:20]:
        inner_seed, inner_coherence = _dominant_region(
            score < -seed_threshold, adjacency, len(areas), areas
        )
        outer_seed, outer_coherence = _dominant_region(
            score > seed_threshold, adjacency, len(areas), areas
        )
        coherence = min(inner_coherence, outer_coherence)
        quality_score = (
            cut_cost
            + 2.0 * (max(0.0, 0.25 - frac) + max(0.0, frac - 0.50))
            + 0.2 * (1.0 - coherence)
        )
        candidate = (quality_score, frac, cut_cost, point, score,
                     inner_seed, outer_seed, inner_coherence, outer_coherence)
        if best is None or candidate[0] < best[0]:
            best = candidate
    if best is None or not np.any(best[5]) or not np.any(best[6]):
        raise HeartError("could not identify both inner and outer seed regions")
    return best


def _geodesic_partition(mesh: trimesh.Trimesh, inner_seed: np.ndarray,
                        outer_seed: np.ndarray) -> np.ndarray:
    """Extend two seed surfaces to a shared geodesic cut on the original mesh."""
    n_faces = len(mesh.faces)
    adjacency = mesh.face_adjacency
    centers = mesh.triangles_center
    edge_length = np.linalg.norm(
        centers[adjacency[:, 0]] - centers[adjacency[:, 1]], axis=1
    )
    inner_ids = np.flatnonzero(inner_seed)
    outer_ids = np.flatnonzero(outer_seed)
    # Add two zero-cost super-source nodes, each connected to every face in one seed.
    rows = np.r_[adjacency[:, 0], adjacency[:, 1],
                 np.full(len(inner_ids), n_faces), inner_ids,
                 np.full(len(outer_ids), n_faces + 1), outer_ids]
    cols = np.r_[adjacency[:, 1], adjacency[:, 0], inner_ids,
                 np.full(len(inner_ids), n_faces), outer_ids,
                 np.full(len(outer_ids), n_faces + 1)]
    weights = np.r_[edge_length, edge_length,
                    np.zeros(2 * (len(inner_ids) + len(outer_ids)))]
    graph = csr_matrix((weights, (rows, cols)), shape=(n_faces + 2, n_faces + 2))
    distances = dijkstra(graph, directed=False, indices=[n_faces, n_faces + 1])
    return distances[0, :n_faces] <= distances[1, :n_faces]


def _face_components(mask: np.ndarray, adjacency: np.ndarray, n_faces: int) -> np.ndarray:
    pairs = adjacency[mask[adjacency[:, 0]] & mask[adjacency[:, 1]]]
    graph = coo_matrix(
        (np.ones(2 * len(pairs)),
         (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
        shape=(n_faces, n_faces),
    )
    return connected_components(graph, directed=False)[1]


# ── the basal rim and the plane ──────────────────────────────────────────────
def _outer_lip_mask(main: trimesh.Trimesh, outer: trimesh.Trimesh,
                    axial_threshold: float = 0.5,
                    radial_threshold: float = -0.2) -> np.ndarray:
    """Flag outer faces on the basal rim where the wall folds over into the cavity.

    The long axis is the closed shell's first principal axis, oriented toward the
    outer surface's open boundary (the base). Rim faces have normals pointing toward
    the base or back toward the axis. Only rim regions connected to the seam are
    removed, so isolated bumps on the epicardial wall stay. Anything the removal
    detaches from the main wall is removed too.
    """
    center = main.vertices.mean(axis=0)
    axis = np.linalg.eigh(np.cov((main.vertices - center).T))[1][:, -1]
    counts = np.bincount(outer.edges_unique_inverse, minlength=len(outer.edges_unique))
    seam_vertices = np.unique(outer.edges_unique[counts == 1])
    axis *= np.sign(axis @ (outer.vertices[seam_vertices].mean(axis=0) - center))

    offset = outer.triangles_center - center
    radial = offset - np.outer(offset @ axis, axis)
    radial /= np.linalg.norm(radial, axis=1)[:, None] + 1e-12
    normals = outer.face_normals
    candidate = ((normals @ axis > axial_threshold)
                 | (np.einsum("ij,ij->i", normals, radial) < radial_threshold))

    n_faces = len(outer.faces)
    adjacency = outer.face_adjacency
    on_seam = np.isin(outer.faces, seam_vertices).any(axis=1)
    labels = _face_components(candidate, adjacency, n_faces)
    lip = candidate & np.isin(labels, np.unique(labels[candidate & on_seam]))

    keep = ~lip
    wall, _ = _dominant_region(keep, adjacency, n_faces, outer.area_faces)
    return ~wall


def _basal_plane(main: trimesh.Trimesh, outer: trimesh.Trimesh,
                 margin_mm: float) -> tuple[np.ndarray, np.ndarray]:
    """Fit a plane to the basal rim and move it toward the apex, below the rim.

    The rim faces come from `_outer_lip_mask`. The plane normal is the rim's
    least-variance direction (oriented to the base), so a tilted base is followed.
    The plane passes through the rim's 1st-percentile depth, then `margin_mm`
    further toward the apex, so no rim remains on either surface.
    """
    lip = _outer_lip_mask(main, outer)
    points, weights = outer.triangles_center[lip], outer.area_faces[lip]
    if len(points) < 3:
        raise HeartError("could not locate the basal rim for the plane cut")
    origin = np.average(points, axis=0, weights=weights)
    normal = np.linalg.eigh(np.cov((points - origin).T, aweights=weights))[1][:, 0]
    normal *= np.sign(normal @ (origin - main.vertices.mean(axis=0)))
    depth = np.percentile((points - origin) @ normal, 1)
    return origin + normal * (depth - margin_mm), normal


def _clip_below(mesh: trimesh.Trimesh, origin: np.ndarray, normal: np.ndarray) -> trimesh.Trimesh:
    """Keep the part of `mesh` on the apex side of the plane, splitting crossing faces."""
    vertices, faces = trimesh.intersections.slice_faces_plane(
        mesh.vertices, mesh.faces, -np.asarray(normal, float), np.asarray(origin, float))[:2]
    return trimesh.Trimesh(vertices, faces, process=True)


# ── holes ────────────────────────────────────────────────────────────────────
def _fill_small_holes(surface: trimesh.Trimesh, max_perimeter_mm: float):
    """Patch every boundary loop except the longest (the basal cut) up to a size.

    Holes through the myocardial wall join the inner and outer walls by a small
    tube; the split cuts that tube, leaving a small loop on each surface. Each
    loop up to `max_perimeter_mm` is closed by a fan of triangles to its
    centroid, wound to match the surface. Returns the surface and the perimeters
    of filled and unfilled holes. With `max_perimeter_mm` 0 nothing is filled,
    and every hole is reported as unfilled.
    """
    # Directed boundary edges: an oriented face edge whose reverse is not used.
    directed = surface.faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2)
    present = {tuple(edge) for edge in directed}
    boundary = np.array([edge for edge in directed if (edge[1], edge[0]) not in present])
    if len(boundary) == 0:
        return surface, [], []
    used, local = np.unique(boundary, return_inverse=True)
    local = local.reshape(-1, 2)
    graph = coo_matrix((np.ones(len(local)), (local[:, 0], local[:, 1])),
                       shape=(len(used), len(used)))
    loop_of_edge = connected_components(graph, directed=False)[1][local[:, 0]]
    lengths = np.linalg.norm(surface.vertices[boundary[:, 0]]
                             - surface.vertices[boundary[:, 1]], axis=1)
    perimeters = np.bincount(loop_of_edge, weights=lengths)
    basal = int(np.argmax(perimeters))

    vertices = [surface.vertices]
    faces = [surface.faces]
    filled, unfilled = [], []
    next_vertex = len(surface.vertices)
    for loop, perimeter in enumerate(perimeters):
        if loop == basal:
            continue
        if perimeter > max_perimeter_mm:
            unfilled.append(float(perimeter))
            continue
        edges = boundary[loop_of_edge == loop]
        vertices.append(surface.vertices[np.unique(edges)].mean(axis=0, keepdims=True))
        faces.append(np.column_stack([edges[:, 1], edges[:, 0],
                                      np.full(len(edges), next_vertex)]))
        next_vertex += 1
        filled.append(float(perimeter))
    if not filled:
        return surface, filled, unfilled
    patched = trimesh.Trimesh(np.vstack(vertices), np.vstack(faces), process=True)
    return patched, filled, unfilled


# ── checks ───────────────────────────────────────────────────────────────────
def _extract_surface(mesh: trimesh.Trimesh, keep_faces: np.ndarray) -> trimesh.Trimesh:
    return mesh.submesh([np.flatnonzero(keep_faces)], append=True, repair=False)


def _surface_qc(surface: trimesh.Trimesh) -> dict:
    components = surface.split(only_watertight=False)
    # Boundary edges are edges used by exactly one triangle.
    counts = np.bincount(surface.edges_unique_inverse, minlength=len(surface.edges_unique))
    boundary = surface.edges_unique[counts == 1]
    if len(boundary):
        used_vertices = np.unique(boundary)
        remap = np.full(len(surface.vertices), -1, dtype=np.int64)
        remap[used_vertices] = np.arange(len(used_vertices))
        edge = remap[boundary]
        degree = np.bincount(edge.ravel(), minlength=len(used_vertices))
        graph = coo_matrix(
            (np.ones(2 * len(edge)),
             (np.r_[edge[:, 0], edge[:, 1]], np.r_[edge[:, 1], edge[:, 0]])),
            shape=(len(used_vertices), len(used_vertices)),
        ).tocsr()
        loops = int(connected_components(graph, directed=False)[0])
        simple_loops = bool(np.all(degree == 2))
    else:
        loops = 0
        simple_loops = False
    return {
        "faces": len(surface.faces),
        "vertices": len(surface.vertices),
        "surface_area_mm2": float(surface.area),
        "components": len(components),
        "boundary_edges": len(boundary),
        "boundary_loops": loops,
        "boundary_is_simple_loops": simple_loops,
        "watertight": bool(surface.is_watertight),
        # Handles (tunnels/bridges) per surface: chi = 2C - 2g - b.
        "handles": int((2 * len(components) - surface.euler_number - loops) // 2),
        "euler": int(surface.euler_number),
    }


def _boundary_signature(surface: trimesh.Trimesh) -> np.ndarray:
    """Sorted endpoint coordinates for each open boundary edge."""
    counts = np.bincount(surface.edges_unique_inverse, minlength=len(surface.edges_unique))
    edges = surface.edges_unique[counts == 1]
    endpoints = np.round(surface.vertices[edges], decimals=5)
    a, b = endpoints[:, 0], endpoints[:, 1]
    swap = ((a[:, 0] > b[:, 0])
            | ((a[:, 0] == b[:, 0]) & (a[:, 1] > b[:, 1]))
            | ((a[:, 0] == b[:, 0]) & (a[:, 1] == b[:, 1]) & (a[:, 2] > b[:, 2])))
    endpoints[swap] = endpoints[swap, ::-1]
    flat = endpoints.reshape(len(endpoints), 6)
    return flat[np.lexsort(tuple(flat[:, i] for i in reversed(range(6))))]


def _rotate(vec, axis, deg):
    """Rodrigues rotation of `vec` about unit `axis` by `deg` degrees."""
    a = np.radians(float(deg))
    k = np.asarray(axis, float)
    v = np.asarray(vec, float)
    return v * np.cos(a) + np.cross(k, v) * np.sin(a) + k * (k @ v) * (1 - np.cos(a))


# ── the session the web app drives ───────────────────────────────────────────
# How much of the outer wall a cut may take before it is called too deep. Across the
# 34 study meshes a rim-fitted cut 3 mm below the rim removed 19-33%; the meshes that
# went past 40% were ones where the rim fit itself had failed.
MAX_REMOVED_FRACTION = 0.40
# Below this share of the cut shell's area a piece is debris from the slice, not a wall.
PIECE_FRACTION = 0.02


class HeartSession:
    """One ventricle, split once and cut as often as the plane moves.

    The expensive part -- seeds and the geodesic seam -- runs once, on load. A recut
    only slices: about a second on a 300k-face shell, so it runs when a slider is
    released, and the disc alone follows the drag.

    The plane is the rim-fitted plane, moved `depth_mm` toward the apex along its
    own normal and tilted by `tilt_a` / `tilt_b` degrees about two axes across it.
    """

    def __init__(self, path, depth_mm=3.0, fill_holes_mm=0.0):
        self.path = Path(path)
        if not self.path.is_file():
            raise HeartError(f"No such file: {self.path}")
        mesh = trimesh.load(str(self.path), force="mesh", process=True)
        if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
            raise HeartError(f"{self.path.name} has no triangles")
        parts = sorted(mesh.split(only_watertight=False),
                       key=lambda c: len(c.faces), reverse=True)
        self.main = parts[0]
        self.n_islands = len(parts) - 1
        self.input_faces = int(len(mesh.faces))
        self.source_sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()
        self.closed = bool(self.main.is_watertight)
        self.notes = []
        if self.n_islands:
            self.notes.append(f"{self.n_islands} detached island(s) ignored; the largest "
                              f"shell is split.")
        if not self.closed:
            self.notes.append("The shell is not closed. The split expects a closed "
                              "myocardium, so check the result closely.")

        (_, self.seed_inner_fraction, _, _, _, inner_seed, outer_seed,
         self.inner_coherence, self.outer_coherence) = _choose_center(self.main)
        self.inner_mask = _geodesic_partition(self.main, inner_seed, outer_seed)
        self.inner_full = _extract_surface(self.main, self.inner_mask)
        self.outer_full = _extract_surface(self.main, ~self.inner_mask)

        self.rim_fitted = True
        try:
            self.rim_origin, self.rim_normal = _basal_plane(self.main, self.outer_full, 0.0)
        except HeartError:
            # No rim found: square to the long axis, a tenth of the length below the
            # base end, and say so -- this start is a guess, not a fit.
            self.rim_fitted = False
            c = self.main.vertices.mean(axis=0)
            axis = np.linalg.eigh(np.cov((self.main.vertices - c).T))[1][:, -1]
            counts = np.bincount(self.outer_full.edges_unique_inverse,
                                 minlength=len(self.outer_full.edges_unique))
            seam = self.outer_full.vertices[np.unique(
                self.outer_full.edges_unique[counts == 1])]
            if len(seam):
                axis *= np.sign(axis @ (seam.mean(axis=0) - c))
            h = (self.main.vertices - c) @ axis
            self.rim_normal = axis
            self.rim_origin = c + axis * (h.max() - 0.1 * np.ptp(h))
            self.notes.append("No basal rim was found, so the plane starts square to "
                              "the long axis. Place it by eye.")
        n = self.rim_normal / np.linalg.norm(self.rim_normal)
        ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        self.axis_a = np.cross(n, ref)
        self.axis_a /= np.linalg.norm(self.axis_a)
        self.axis_b = np.cross(n, self.axis_a)
        self.rim_normal = n
        self.radius = 0.6 * float(np.ptp(self.main.vertices, axis=0).max())

        self.depth_mm = float(depth_mm)
        self.tilt_a = 0.0
        self.tilt_b = 0.0
        self.fill_holes_mm = float(fill_holes_mm)
        self.inner = self.outer = self.removed = None
        self.recut()

    # ── the plane ───────────────────────────────────────────────────────────
    def plane(self, depth_mm=None, tilt_a=None, tilt_b=None):
        """(origin, normal). The normal points toward the base: what is on its
        positive side is cut off."""
        depth = self.depth_mm if depth_mm is None else float(depth_mm)
        ta = self.tilt_a if tilt_a is None else float(tilt_a)
        tb = self.tilt_b if tilt_b is None else float(tilt_b)
        normal = _rotate(_rotate(self.rim_normal, self.axis_a, ta), self.axis_b, tb)
        normal /= np.linalg.norm(normal)
        origin = self.rim_origin - self.rim_normal * depth
        return origin, normal

    def disc(self, **kw):
        origin, normal = self.plane(**kw)
        return {"origin": origin, "normal": normal, "radius": self.radius}

    def set_plane(self, depth_mm=None, tilt_a=None, tilt_b=None):
        if depth_mm is not None:
            self.depth_mm = float(depth_mm)
        if tilt_a is not None:
            self.tilt_a = float(tilt_a)
        if tilt_b is not None:
            self.tilt_b = float(tilt_b)

    def reset_plane(self, depth_mm=3.0):
        self.depth_mm, self.tilt_a, self.tilt_b = float(depth_mm), 0.0, 0.0

    # ── the cut ─────────────────────────────────────────────────────────────
    def recut(self):
        """Cut, split, handle holes, and judge the result. Sets inner, outer,
        removed, info, ok and problem."""
        origin, normal = self.plane()
        cut = _clip_below(self.main, origin, normal)
        self.removed = _clip_below(self.main, origin, -normal)
        problems = []
        notes = list(self.notes)

        if len(cut.faces) == 0:
            self.inner = self.outer = None
            self.ok, self.problem = False, "The plane removed the whole ventricle."
            self.info = {"notes": notes}
            return self

        pieces = sorted(cut.split(only_watertight=False), key=lambda p: p.area,
                        reverse=True)
        total = sum(p.area for p in pieces)
        big = [p for p in pieces if p.area >= PIECE_FRACTION * total]
        self.separated = len(big) == 2
        if self.separated:
            # The walls only touch through holes, so two pieces is a clean split.
            # The outer wall's normals point away from the axis, the inner's toward.
            def outward(p):
                off = p.triangles_center - origin
                radial = off - np.outer(off @ normal, normal)
                return np.average(np.einsum("ij,ij->i", p.face_normals, radial),
                                  weights=p.area_faces)
            inner, outer = sorted(big, key=outward)
        else:
            # Still joined below the plane: the geodesic seam cuts each hole's tube.
            inner = max(_clip_below(self.inner_full, origin, normal)
                        .split(only_watertight=False), key=lambda p: p.area)
            outer = max(_clip_below(self.outer_full, origin, normal)
                        .split(only_watertight=False), key=lambda p: p.area)
            if len(big) == 1:
                notes.append("The walls are still joined below the plane, through "
                             "holes in the wall; the automatic seam cuts them.")

        inner, inner_filled, inner_open = _fill_small_holes(inner, self.fill_holes_mm)
        outer, outer_filled, outer_open = _fill_small_holes(outer, self.fill_holes_mm)
        qi, qo = _surface_qc(inner), _surface_qc(outer)
        removed_fraction = 1.0 - outer.area / max(self.outer_full.area, 1e-12)
        holes_open = len(inner_open) + len(outer_open)
        holes_filled = len(inner_filled) + len(outer_filled)

        if cut.area > 0.995 * self.main.area:
            problems.append("The plane sits above the ventricle: nothing was cut. "
                            "Move it toward the apex.")
        if qi["components"] != 1 or qo["components"] != 1:
            problems.append("A wall came out in more than one piece. Move the plane "
                            "or tilt it.")
        if removed_fraction > MAX_REMOVED_FRACTION:
            problems.append(f"The cut takes {100 * removed_fraction:.0f}% of the outer "
                            f"wall, more than {100 * MAX_REMOVED_FRACTION:.0f}%. Move "
                            f"the plane toward the base.")
        if holes_open:
            notes.append(f"{holes_open} hole(s) left open: the reference method leaves "
                         f"them open rather than adding surface.")
        if holes_filled:
            notes.append(f"{holes_filled} hole(s) patched, largest "
                         f"{max(inner_filled + outer_filled):.0f} mm around.")

        self.inner, self.outer = inner, outer
        self.ok = not problems
        self.problem = problems[0] if problems else None
        self.info = {
            "separated": bool(self.separated),
            "pieces": len(big),
            "removed_pct": round(float(100 * removed_fraction), 1),
            "inner_area_mm2": round(float(inner.area), 1),
            "outer_area_mm2": round(float(outer.area), 1),
            "inner_faces": qi["faces"], "outer_faces": qo["faces"],
            "inner_euler": qi["euler"], "outer_euler": qo["euler"],
            "inner_loops": qi["boundary_loops"], "outer_loops": qo["boundary_loops"],
            "inner_handles": qi["handles"], "outer_handles": qo["handles"],
            "holes_open": holes_open, "holes_filled": holes_filled,
            "hole_perimeters_mm": [round(float(p), 1) for p in
                                   inner_open + outer_open + inner_filled + outer_filled],
            "rim_fitted": self.rim_fitted,
            "problems": problems,
            "notes": notes,
        }
        return self

    # ── output ──────────────────────────────────────────────────────────────
    def sidecar(self):
        origin, normal = self.plane()
        return {
            "source": self.path.name,
            "source_sha256": self.source_sha256,
            "method": ("normal seeds + geodesic seam; rim-fitted basal plane cut; "
                       "holes " + ("patched up to %g mm" % self.fill_holes_mm
                                   if self.fill_holes_mm > 0 else "left open")),
            "plane": {"origin": origin.tolist(), "normal": normal.tolist(),
                      "depth_mm": self.depth_mm, "tilt_a_deg": self.tilt_a,
                      "tilt_b_deg": self.tilt_b, "rim_fitted": self.rim_fitted,
                      "rim_origin": self.rim_origin.tolist(),
                      "rim_normal": self.rim_normal.tolist()},
            "fill_holes_mm": self.fill_holes_mm,
            "seed_coherence": {"inner": float(self.inner_coherence),
                               "outer": float(self.outer_coherence)},
            "input_faces": self.input_faces,
            "islands_ignored": self.n_islands,
            "valid": bool(self.ok),
            "problem": self.problem,
            **{k: v for k, v in self.info.items() if k not in ("notes", "problems")},
            "notes": self.info.get("notes", []),
        }

    def export(self, out_dir, stem):
        """Write <stem>_inner.stl, <stem>_outer.stl and <stem>_cut.json."""
        if not self.ok:
            raise HeartError(f"The cut is not valid: {self.problem}")
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        paths = {}
        for side, surface in (("inner", self.inner), ("outer", self.outer)):
            paths[side] = out_dir / f"{stem}_inner.stl" if side == "inner" \
                else out_dir / f"{stem}_outer.stl"
            surface.export(paths[side])
        paths["sidecar"] = out_dir / f"{stem}_cut.json"
        paths["sidecar"].write_text(json.dumps(self.sidecar(), indent=2) + "\n")
        return paths
