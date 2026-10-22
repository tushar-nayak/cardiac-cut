"""boundary_loops must give the same loops as it always has, and not scale with vertex count.

THE BUG. boundary_loops finds the connected components of a graph whose nodes are EVERY
vertex of the mesh but whose edges are only the boundary edges. Every interior vertex is
therefore a component of its own, so the loop that gathered each component's members --
one Python iteration per component, each scanning the whole boundary -- ran about once
per VERTEX, not once per opening. Cost: vertices x boundary vertices.

WHAT IT COST. It is called on the whole mesh for every trial cut: 56 times in the first
cut of a 159k-face aorta, where it was 6.3 s of the 17 s. At 637k faces the first cut
took 77 s and the load 22 s; grouping the boundary vertices by component label in one
sort took them to 48 s and 18 s. Measured with the web app's own sequence
(vascular_api.CutSession, then nub_preset + recut) on KY160 and JK5 scaled 1x, 4x and
16x: 290 real calls compared array for array, all identical, and the cut mesh, presets
and openings identical byte for byte on all six.

WHAT THIS PINS.
  1. The answer is unchanged: every loop, in the same order, with the same members in
     the same order and the same dtype, against the old implementation kept verbatim
     below -- on synthetic edge cases, on random holes, and on both sample aortas.
  2. The cost no longer scales with vertices x boundary. Timed against the topology
     pass the function itself starts with (group_rows over the edges, on a fresh mesh),
     so the bound is a ratio on the same machine rather than a wall-clock guess: the new
     code runs at ~1.1x that reference on the sheet below, the old one at ~55x.
"""
import sys
import time
from pathlib import Path

import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import branch_cut as bc                                            # noqa: E402

fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def reference(mesh):
    """The implementation this replaced, verbatim. Correct, and O(V x B)."""
    lone = trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)
    if len(lone) == 0:
        return []
    edges = mesh.edges_sorted[lone]
    n = len(mesh.vertices)
    g = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n, n))
    ncomp, labels = connected_components(g, directed=False)
    used = np.unique(edges)
    loops = []
    for c in range(ncomp):
        members = used[labels[used] == c]
        if len(members) >= 3:
            loops.append(members)
    return loops


def fresh(m):
    """A new Trimesh over the same arrays, so neither side reads the other's cache."""
    return trimesh.Trimesh(np.array(m.vertices), np.array(m.faces), process=False)


def same(a, b):
    return (len(a) == len(b)
            and all(x.dtype == y.dtype and np.array_equal(x, y) for x, y in zip(a, b)))


def agrees(name, m):
    a, b = bc.boundary_loops(fresh(m)), reference(fresh(m))
    check(name, same(a, b), f"{len(b)} loops expected, got {len(a)}"
          if len(a) != len(b) else f"{len(b)} loops")
    return a


def punch(m, rng, holes, radius):
    """Delete the faces within `radius` of `holes` random vertices."""
    c = m.vertices[rng.choice(len(m.vertices), holes, replace=False)]
    d = np.linalg.norm(m.triangles_center[:, None, :] - c[None], axis=2).min(axis=1)
    return trimesh.Trimesh(m.vertices, m.faces[d > radius], process=False)


def perforated_sheet(n, every):
    """An n x n vertex grid of triangles with a one-cell hole every `every` cells."""
    i, j = np.meshgrid(np.arange(n - 1), np.arange(n - 1), indexing="ij")
    keep = ~((i % every == every // 2) & (j % every == every // 2))
    i, j = i[keep], j[keep]
    v0 = i * n + j
    v1, v2, v3 = v0 + 1, v0 + n, v0 + n + 1
    faces = np.concatenate([np.column_stack([v0, v1, v3]),
                            np.column_stack([v0, v3, v2])])
    x, y = np.meshgrid(np.arange(n, dtype=float), np.arange(n, dtype=float),
                       indexing="ij")
    return trimesh.Trimesh(np.column_stack([x.ravel(), y.ravel(), np.zeros(n * n)]),
                           faces, process=False)


print("=" * 78)
print("1. same loops as before: synthetic edge cases")
print("=" * 78)
sphere = trimesh.creation.icosphere(subdivisions=3)
check("a closed surface has no openings", bc.boundary_loops(fresh(sphere)) == [])
tri = trimesh.Trimesh([[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 1, 2]], process=False)
agrees("a single triangle is one 3-vertex loop", tri)
tube = trimesh.creation.cylinder(radius=5, height=20, sections=32)
tube = trimesh.Trimesh(tube.vertices, tube.faces[np.abs(tube.face_normals[:, 2]) < 0.5],
                       process=False)
loops = agrees("an open tube", tube)
check("  ...has exactly its two ends", len(loops) == 2, f"{len(loops)}")
agrees("a sheet with 25 holes", perforated_sheet(60, 12))
# two triangles meeting at one vertex: both rims pass through it, so they are one
# component and one loop, as they always were
bow = trimesh.Trimesh([[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
                      [[0, 1, 2], [0, 3, 4]], process=False)
agrees("a bow-tie pinched at one vertex", bow)
# a mesh whose vertex array carries vertices no face uses
spare = trimesh.Trimesh(np.vstack([tube.vertices, np.random.default_rng(0).random((50, 3))]),
                        tube.faces, process=False)
agrees("unreferenced vertices in the array", spare)
# a fan with an edge shared by three faces (non-manifold): the lone edges still close
fan = trimesh.Trimesh([[0, 0, 0], [0, 0, 1], [1, 0, 0], [0, 1, 0], [-1, 0, 0]],
                      [[0, 1, 2], [0, 1, 3], [0, 1, 4]], process=False)
agrees("an edge shared by three faces", fan)
agrees("int32 faces", trimesh.Trimesh(tube.vertices, tube.faces.astype(np.int32),
                                      process=False))

print("\n" + "=" * 78)
print("2. same loops as before: random holes, and both sample aortas")
print("=" * 78)
rng = np.random.default_rng(20260925)
big = trimesh.creation.icosphere(subdivisions=4)
bad = 0
for k in range(120):
    m = punch(big, rng, holes=int(rng.integers(1, 12)), radius=float(rng.uniform(0.05, 0.4)))
    if not same(bc.boundary_loops(fresh(m)), reference(fresh(m))):
        bad += 1
check("120 spheres with random holes", bad == 0, f"{bad} differ")

for name, mat in (("KY160_1", ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"),
                  ("JK5_1", ROOT / "JK5_1_M10.mat")):
    if not mat.exists():
        check(f"{name} sample mesh present", False, str(mat))
        continue
    m = bc.load_mat_mesh(mat)
    agrees(f"{name} as loaded", m)
    for k in range(3):
        agrees(f"{name} with {3 + 4 * k} holes punched",
               punch(m, rng, holes=3 + 4 * k, radius=4.0))

print("\n" + "=" * 78)
print("3. cost does not scale with vertices x boundary")
print("=" * 78)
N, EVERY, BOUND = 600, 6, 8.0


def best_of(fn, reps=3):
    t = []
    for _ in range(reps):
        m = perforated_sheet(N, EVERY)            # fresh: no cached edges either side
        t0 = time.perf_counter()
        fn(m)
        t.append(time.perf_counter() - t0)
    return min(t)


t_ref = best_of(lambda m: trimesh.grouping.group_rows(m.edges_sorted, require_count=1))
t_new = best_of(bc.boundary_loops)
m = perforated_sheet(N, EVERY)
print(f"  sheet: {len(m.vertices):,} vertices, {len(m.faces):,} faces, "
      f"{len(bc.boundary_loops(m)):,} openings")
print(f"  reference topology pass {t_ref:.3f} s, boundary_loops {t_new:.3f} s "
      f"({t_new / t_ref:.1f}x; the per-vertex loop was ~55x)")
check(f"boundary_loops within {BOUND:.0f}x the edge pass it starts with",
      t_new <= BOUND * t_ref, f"{t_new / t_ref:.1f}x")

print()
if fails:
    print(f"{len(fails)} FAILURE(S):")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
