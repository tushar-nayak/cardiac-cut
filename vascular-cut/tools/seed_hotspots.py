"""WHERE does delta K change when only the partition seed changes?

The remesh is deterministic and the seed only touches CAMP's k-means, so every seed
sees the identical mesh and face indices line up one to one. That makes it possible to
ask the question spatially rather than in aggregate.

Each patch's contribution to delta K is spread over its own faces in proportion to
area, giving a per-face contribution map per seed. Across seeds:

    mean   where delta K comes from at all
    SD     where it MOVES -- the instability, localised

Then the high-variance faces are located against the mesh's openings, since the
suspicion all along has been that the ostia are responsible.

    python tools/seed_hotspots.py CASE [faces]
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tools"))
import trimesh                                                     # noqa: E402
import shape_metrics as SM                                         # noqa: E402
from camp_cohort import repair_folds                               # noqa: E402

CID = sys.argv[1] if len(sys.argv) > 1 else "KK20_2"
FACES = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000
SEEDS = [42, 43, 44, 45, 46, 47, 48, 49]

base = SM.weld(trimesh.load(str(ROOT / "output" / "nubs_final" / f"{CID}_nub.stl"),
                            process=False))
chi = int(base.euler_number)
base, _ = repair_folds(base)

maps, mesh = [], None
for sd in SEEDS:
    orig = SM.PARAMS.copy(); SM.PARAMS["seed"] = sd
    try:
        r = SM.delta_k(base, target_faces=FACES, trim_frac=0.0, make_figure=False,
                       expect_euler=chi)
    finally:
        SM.PARAMS.clear(); SM.PARAMS.update(orig)
    K = np.asarray(r["integrated_K"], float); A = np.asarray(r["patch_area"], float)
    lab = np.asarray(r["face_labels"], int)
    m = r["mesh"]; fa = np.asarray(m.area_faces, float)
    w = A / A.sum(); contrib = w * (K - float((w * K).sum())) ** 2
    share = contrib / max(contrib.sum(), 1e-30)
    # spread each patch's share over its faces, by area
    per_face = share[lab] * fa / np.maximum(A[lab], 1e-30)
    maps.append(per_face)
    if mesh is None:
        mesh = m
    print(f"  seed {sd}: {len(K)} patches, dK {contrib.sum():.5f}", flush=True)

M = np.vstack(maps)
mu, sd = M.mean(axis=0), M.std(axis=0, ddof=1)
fa = np.asarray(mesh.area_faces, float)
order = np.argsort(sd)[::-1]
cum_var = np.cumsum(sd[order] ** 2) / np.sum(sd[order] ** 2)
for frac in (0.001, 0.005, 0.01, 0.05, 0.10):
    n = int(frac * len(sd))
    print(f"  top {100*frac:5.1f}% of faces ({n:6,}) hold "
          f"{100*cum_var[max(n-1,0)]:5.1f}% of the across-seed variance, "
          f"{100*fa[order[:n]].sum()/fa.sum():4.1f}% of area")

# where are they, relative to the openings?
V = np.asarray(mesh.vertices); F = np.asarray(mesh.faces)
e = np.sort(np.vstack([F[:, [0,1]], F[:, [1,2]], F[:, [2,0]]]), axis=1)
uq, cnt = np.unique(e, axis=0, return_counts=True)
be = uq[cnt == 1]
vb = np.unique(be); ix = {int(x): i for i, x in enumerate(vb)}
par = np.arange(len(vb))
def find(x):
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x
for a, b in be:
    ra, rb = find(ix[int(a)]), find(ix[int(b)])
    if ra != rb: par[ra] = rb
loops = {}
for v in vb: loops.setdefault(find(ix[int(v)]), []).append(int(v))
rims = []
for _, vs in loops.items():
    pts = V[vs]; c = pts.mean(axis=0)
    rims.append((c, float(np.linalg.norm(pts - c, axis=1).mean()), np.array(vs)))
rims.sort(key=lambda t: t[1])
branch = rims[:-2]          # the two widest are inlet and outlet
tc = np.asarray(mesh.triangles_center, float)
from scipy.spatial import cKDTree
d_branch = cKDTree(np.vstack([V[b[2]] for b in branch])).query(tc)[0]
d_rim = cKDTree(np.vstack([V[r[2]] for r in rims])).query(tc)[0]
top = order[:max(int(0.005 * len(sd)), 1)]
print(f"\n  {len(rims)} openings; {len(branch)} treated as branch ostia")
print(f"  distance to nearest BRANCH ostium: top-variance faces "
      f"{np.median(d_branch[top]):.1f} mm vs {np.median(d_branch):.1f} mm for all faces")
print(f"  distance to nearest opening rim  : top-variance faces "
      f"{np.median(d_rim[top]):.1f} mm vs {np.median(d_rim):.1f} mm for all faces")
print(f"  share of top-variance faces within 10 mm of a branch: "
      f"{100*np.mean(d_branch[top] < 10):.0f}%  (all faces: {100*np.mean(d_branch < 10):.0f}%)")
np.savez_compressed(ROOT / "output" / "camp_diag" / f"{CID}_seed_hotspots.npz",
                    verts=V.astype(np.float32), faces=F.astype(np.int32),
                    mean=mu.astype(np.float32), sd=sd.astype(np.float32),
                    d_branch=d_branch.astype(np.float32))
print(f"  -> output/camp_diag/{CID}_seed_hotspots.npz")
