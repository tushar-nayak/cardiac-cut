# CAMP → Vascular-Cut-Tool integration guide

Written for the chat building the vascular-cut web tool. It answers your 8 questions, but
**read this first section**, most of the questions assume the *Khabaz* method, and our CAMP
works differently. Getting that difference straight resolves about half your confusion.

> **This is the design document, written before the cohort was run.** Three things in
> it were later settled by measurement and are now superseded — see
> [`CURVATURE.md`](CURVATURE.md) for the data:
>
> - **`trim_frac=0.01` in the code samples below is not the operating point.** It is
>   **0.0025**, chosen because at trim 0 the same surface remeshed to 98k/100k/102k
>   faces gives δK values differing by 2.6× typically and 23.5× at worst.
> - **"We did not pre-normalise" still holds for what is stored** (raw `k_fluct` +
>   `A_mm2` per mesh), but the cohort figure divides by the normal groups' mean
>   (KY+SA, n=12) so the groups are readable on a linear axis.
> - **δK is not reproducible to better than a factor**, which is the single most
>   important caveat for anyone reading a value off the site. The guide below does not
>   mention it because nobody had measured it yet.

---

## 0. CAMP is NOT the Khabaz partition. This is the key thing.

The Khabaz et al. 2023 δK paper and *our* CAMP compute the **same descriptor** (fluctuation of
per-patch integrated Gaussian curvature) but partition the surface **completely differently**:

| | Khabaz δK (the paper) | **CAMP (our code)** |
|---|---|---|
| partition | k-means Voronoi, **fixed k** = A_T/ℓ², ℓ = aortic radius | **convergent adaptive subdivision**, no target k |
| patch count | you choose it (via ℓ) | **an OUTPUT**, decided by the geometry + quality gates |
| patch size | ~uniform (equal-area) | **varies with local shape** (small at ostia, large on smooth wall) |
| stopping rule | k reached | every patch passes R²≥0.9 AND η≤0.075, or hits a vertex floor |
| variance | unweighted over equal-area cells | **area-weighted** (because our patches are unequal area) |

So when you ask "what sets k" or "k = A_T/ℓ²?", **CAMP has no k to set.** You give it two quality
thresholds and it grows/splits patches until each one is a good local quadratic fit. The patch
count falls out (typically ~800 at 50k faces, ~950–1300 at 100k for a thoracic aorta).

Everything else below follows from this.

---

## 1. Where the code is (copy these)

All under `c:\Users\cmdavis2\Projects\camp\`.

**The CAMP method, most up-to-date version:**
`UpdatedCamp\campsrc.py`  ← THIS one, not the root `campsrc.py` (that's older, seed_patches=250, different trimming)
- `camp(mat_path, ...)`, entry point (line 358). Runs the whole partition + curvature + figure.
- `patch_pass(...)` (line 252), the convergent subdivision loop.
- `make_plotly_patch_figure(...)` (line 176), the patch figure builder.
- `camp_batch(...)` (line 677), reference for the **k_fluct formula** (lines 766–771).
- `_trim_mask_integrated(...)` (line 629), percentile trimming on integrated K.

**Mesh preparation (run BEFORE camp):**
`cutclean.py`
- `acvd_uniform_remesh(mesh, target_faces=50_000)` (line 173), uniform ACVD remesh via `pyacvd`.
- `trimesh_to_pyvista` / `pyvista_to_trimesh` (lines 155–171), the converters it uses.

`run_camp_abdominal.py`
- `weld_mesh(mesh, decimals=5)` (line 47), merges STL "triangle soup" into a shared-vertex mesh.
  **Required** or CAMP silently misbehaves (see Q8).
- `_write_temp_mat(welded, stem)` (line 82), serialises a mesh to the `.mat` camp() expects.

**The end-to-end pipeline to copy the flow from:**
`run_branches_camp.py` → `_worker()`. This is exactly your use case: STL → weld → ACVD remesh →
camp → k_fluct at several trims. Read this function; it's the template.

**Patch visualisations (the integrated-K plasma views you want):**
`view_branch_patches.py`
- `figure(stem, metric="integrated_K", colorscale="Plasma", ...)`, interactive Plotly figure,
  patches coloured by integrated_K, patch borders drawn, colorscale dropdown + cmin/cmax sliders.
- `contact_sheet(stems, png, ...)`, static grid on a shared colour scale.
`camp_recolor.py`, `build_patch_figure` / `build_bounds_slider_figure` recolour a cached result
without re-running camp (instant; good for a "play with the color axis" UI).

---

## 2. The pipeline, as a drop-in wrapper

`camp()` only reads a `.mat`. Your tool has STL/trimesh. Here's a thin wrapper that takes a mesh
and returns everything the website needs (scalar δK, per-patch table, and the figure). This IS the
integration point, wire your cut-tool output straight into `camp_delta_k(mesh)`.

```python
import numpy as np, trimesh
import importlib.util
from pathlib import Path

# load the UP-TO-DATE campsrc explicitly (avoid grabbing the old root one)
_spec = importlib.util.spec_from_file_location(
    "campsrc_new", r"c:\Users\cmdavis2\Projects\camp\UpdatedCamp\campsrc.py")
cu = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(cu)

import run_camp_abdominal as R      # weld_mesh, _write_temp_mat
import cutclean as C                # acvd_uniform_remesh

PARAMS = dict(R2_thresh=0.9, eta_thresh=0.075, min_patch_size=18, seed_patches=100)

def camp_delta_k(mesh, target_faces=100_000, trim_frac=0.01, make_figure=True):
    """mesh: a trimesh (from your cut STL).  Returns dict with scalar + per-patch table + figure."""
    welded, winfo = R.weld_mesh(mesh)                       # triangle soup -> shared vertices
    rm = C.acvd_uniform_remesh(welded, target_faces=target_faces)   # uniform remesh
    tmp = R._write_temp_mat(rm, "upload")
    try:
        out = cu.camp(mat_path=str(tmp), scan_prefix=None, trim_frac=trim_frac,
                      color_metric="integrated_K", color_scale="Plasma",
                      make_figure=make_figure, **PARAMS)
    finally:
        Path(tmp).unlink(missing_ok=True)

    iK   = np.asarray(out["K_integrated_patch"], float)     # per-patch  K_j = <kg>_j * A_j
    area = np.asarray(out["patch_area"], float)
    keep = np.isfinite(iK) & (area > 0)
    q, a = iK[keep], area[keep]
    w = a / a.sum()
    k_fluct = float(np.sum(w * (q - np.sum(w*q))**2))       # area-weighted variance  ( = delta K )

    return dict(
        k_fluct     = k_fluct,                    # THE scalar (a variance; see Q6)
        A_mm2       = float(out["mesh"].area),    # total surface area  (the x-axis of the dK-vs-A curve)
        n_patches   = int(keep.sum()),
        euler       = int(winfo["euler_number"]), # topology QC (see Q8)
        # per-patch table for colouring patches onto the mesh in the viewer:
        face_labels = np.asarray(out["face_labels"]),   # face -> patch id  (on the REMESHED mesh)
        avg_kg      = np.asarray(out["avg_kg"]),         # <kg>_j : K at each patch centroid
        patch_area  = area,                              # A_j
        integrated_K= iK,                                # K_j = <kg>_j * A_j
        eta_patch   = np.asarray(out["eta_patch"]),
        R2_patch    = np.asarray(out["R2_patch"]),
        vertices    = np.asarray(out["mesh"].vertices),  # remeshed geometry for the viewer
        faces       = np.asarray(out["mesh"].faces),
        figure      = out["figure"],                     # ready-made Plotly figure (plasma patches)
    )
```

Runtime: ~90 s/mesh at 50k, ~200–220 s/mesh at 100k, single core (see Q8). **Queue it, don't inline.**

---

## Your 8 questions

### 1. Branch ostia, open holes, lid, or excluded?

**Leave them as OPEN holes. Do not lid, do not fill.** That's what our validated meshes are:
Euler characteristic −4, i.e. **six open boundary loops** (aortic inlet + aortic outlet + 4
supra-aortic cuts), single connected component.

- CAMP is topology-agnostic, it fits local Monge charts, so it does not care whether a hole is
  open. Open boundaries are fine; the fit near a rim just uses the faces that are there.
- **Lidding is actively harmful.** A synthetic cap adds fake flat/curved patches with their own
  curvature, and the ostia are exactly where the real signal lives, you'd be pouring noise onto
  your best feature. Don't do it.
- Excluding the branches entirely (cutting them off flush) is the *Khabaz* convention and is a
  legitimate alternative, but then you're not doing the branched-aorta thing that's the whole
  point. Keep them, open.

One caveat we measured: the ostia carry the highest |K|, and at uniform 100k the branch tubes
below ~6.5 mm local radius are under-resolved (they hit the vertex floor). They still contribute;
they're just noisier there. That's a known operating limit, not a reason to lid.

### 2. End cuts define the size metric, standardise the axial extent? **(the important one)**

**Yes. The cut protocol is part of the method.** Lead with this on the website.

`k_fluct` is **extensive**, it is a variance of K_j = ⟨kg⟩_j·A_j, and total surface area A_T is
literally the x-axis you plot δK against. Where you place the root and distal cuts sets A_T, and it
moves δK. We measured this directly: just lengthening the branches (same patient, more tube)
changed δK by ~2.35× at trim 0.01. So δK is **not** insensitive to how much tube you include.

Two consequences for the site:
- For any cross-patient comparison you **must** fix an axial extent (e.g. sinotubular junction →
  a defined distal landmark) and the same supra-aortic cut convention. Bake it into the cut tool.
- The intended comparison is not "same δK" but **position on the δK-vs-A master curve**, you plot
  the upload as a point in (A_T, δK) space against the reference cohort. That absorbs *some* size
  variation, but only if the cut convention is consistent; it does not license arbitrary cuts.

### 3. What sets k in practice?

**Nothing does, there is no k.** (See §0.) CAMP uses convergent adaptive subdivision:
1. seed with `seed_patches=100` k-means clusters on the vertices,
2. split any patch that fails **R² ≥ 0.9** (fit quality) **or** **η ≤ 0.075** (curvature homogeneity),
3. stop when every patch passes both or is too small to split (`min_patch_size = 18` vertices).

The final patch count is an emergent output. Typical for a thoracic aorta: **~800 patches at 50k
faces, ~950–1300 at 100k.** If you want the Khabaz-style fixed-k behaviour instead, that's a
different function (root `campsrc.py`), but the up-to-date method is adaptive and you should use it.

The one number with the flavour of ℓ: the η gate is equivalent to **patch radius ≤ 0.39 × local
radius of curvature** (we derived and verified this). So patch size *is* tied to local geometry -
just adaptively, per patch, not via a single global ℓ.

### 4. Is 50k a floor or a density? Faces or vertices?

**Faces.** `acvd_uniform_remesh(mesh, target_faces=N)` targets a **face count** (≈ N/2 vertices),
a fixed count per mesh, **not** a density.

**IMPORTANT correction (supersedes an earlier draft of this guide):** the branched thoracic
reference cohort ran at a **fixed target_faces = 100,000**, NOT a fixed density. Across that cohort
(A = 26k–98k mm²) the density therefore *varied* from **0.51 to 1.90 verts/mm² (3.7×), median
1.01**, spacing 0.73–1.40 mm. So to place an upload **on that existing curve**, match what the
cohort did: **target_faces = 100,000**, for any geometry of similar extent. Do **not** scale faces
to hold density constant, that would partition uploads more finely/coarsely than the reference and
shift δK off the curve.

The reason this "works" despite variable density: `pct_at_min_size` (floor-pinned fraction) stayed
**~29% across the whole 3.7× density range** (see Q2), so the partition sat in one convergence
regime and the master curve absorbed the residual density scatter.

The mesh-dependence you worried about is real but bounded: it bites hard only when density leaves
the reference's 0.5–1.9 verts/mm² band, which is exactly what happens for a small isolated region
at fixed 100k faces (see Q2/Q3). For a fresh cohort built from scratch, constant density would be
the cleaner design; but this reference is fixed-faces, so uploads onto it must be fixed-faces too.

### 5. Why thoracic only?

Not a code constraint, CAMP runs on any surface (we've run it on abdominal aortas too). It's two
things:
- **Reference cohort.** The normal / dissection labels and the δK–A curve we trust are thoracic.
- **Where the signal lives.** The arch is a *bent tube*, so it has genuinely non-zero Gaussian
  curvature (torus-like: positive on the outer curve, negative on the inner). A straight infrarenal
  segment is developable (K ≈ 0) except at the sac, less for a Gaussian-curvature descriptor to
  grab. So δK is *more meaningful* on the arch.

The site *can* offer δK on infrarenal/other geometries; you just have no reference cohort to place
the point against, and the interpretation differs. Gate it behind "thoracic (reference available)"
vs "other (raw δK only)".

### 6. Entry point, return, and which quantity? **(units, read carefully)**

- **Entry:** `camp()` in `UpdatedCamp\campsrc.py`. Use the `camp_delta_k()` wrapper in §2.
- **Returns the full per-patch table**, not just the scalar. Keys you need for the viewer:
  `face_labels` (face→patch on the remeshed mesh), `avg_kg` (⟨kg⟩_j, K at patch centroid),
  `patch_area` (A_j), `K_integrated_patch` (K_j = ⟨kg⟩_j·A_j), `eta_patch`, `R2_patch`, plus
  `mesh` and a ready-made `figure`. So you can colour patches onto the mesh directly.
- **The scalar is a VARIANCE, not its square root.** Exactly:
  `k_fluct = Σ_j w_j (K_j − K̄)²`,  `w_j = A_j / A_T`,  `K_j = ⟨kg⟩_j · A_j`.
  This matches Khabaz's `δK = ⟨K²⟩ − ⟨K⟩²` in *form* (both are the variance of per-patch integrated
  curvature), **with two differences that will silently break a naive port:**
  1. **Area-weighted**, because our patches are unequal area. Khabaz's cells are ~equal area so his
     unweighted mean ≈ our weighted mean; don't "fix" this to unweighted.
  2. It is the variance of **integrated K (K_j = ⟨kg⟩_j·A_j, dimensionless)**, *not* the variance of
     the curvature ⟨kg⟩_j itself. If you compute variance of ⟨kg⟩ you'll get a different quantity
     with units mm⁻⁴ and the normalisation will be wrong. Use K_j.
  Units of k_fluct: K_j is dimensionless (mm⁻²·mm²), so k_fluct is **dimensionless**.
  If a UI wants a "δK" with the same units as K_j, take √k_fluct, but the stored/normalised
  quantity is the variance.

### 7. Reference cohort format

Khabaz normalises by `⟨δK⟩_norm` = mean δK of the normal aortas. **We did not pre-normalise**, our
deliverables store **raw** k_fluct + A_mm2 per mesh, and we place points on the raw δK-vs-A curve.

For the website, the reference cohort should be, in order of usefulness:
1. **The fitted master curve** (log δK vs log A regression: slope ≈ 1) **plus the normal-cohort
   scatter**, so an upload is shown as a point against the band. This is the right UX.
2. If you want a single Khabaz-style normaliser: `⟨k_fluct⟩` over the normal (KY) cohort at your
   chosen trim, and report `k_fluct / ⟨k_fluct⟩_norm`.

The raw values are in `Branches_Metrics\branches_kfluct_long_100k.csv` (thoracic branched, per-mesh,
per-trim) and `AAA_Metrics\CAMP_kfluct_trim_sweep.xlsx` (abdominal). Ship a small CSV of
(A_mm2, k_fluct, group) at your production trim as the reference, computed by the same pipeline the
site uses, do NOT mix trims or face counts between reference and uploads.

### 8. Failure behaviour and runtime

**Failure modes to guard against, CAMP does NOT validate input and can return a quietly-wrong number:**
- **Triangle soup (un-welded STL).** STL stores each facet independently (~3× duplicate vertices,
  empty face-adjacency). If you skip `weld_mesh`, CAMP's connectivity repair becomes a silent no-op
  and patch borders vanish, you'll get a number, and it'll be wrong. **Always weld before remesh.**
  (ACVD also needs connectivity to cluster, so welding is doubly required.)
- **Non-manifold / self-intersecting uploads.** ACVD usually cleans these, but can occasionally
  change topology (we saw Euler number shift on a couple of hard cases). Do a **topology QC after
  remesh**: check `euler_number` and connected-component count, and that face count landed near
  target. Reject or flag meshes where remesh changed the opening count.
- **Units.** CAMP assumes millimetres (K reported in mm⁻²). Verify/scale on upload.

**Runtime (single core, `threadpool_limits(1)`):** ~90 s/mesh at 50k faces, **~200–220 s/mesh at
100k.** So at the 100k density you want, this is a **~3–4 minute queued job**, not an inline
request. It parallelises trivially across uploads (we run 6 at once with joblib/loky). The CAMP
split loop dominates; the ACVD remesh is ~15–25 s of that.

---

## TL;DR for the other chat

1. Use `UpdatedCamp\campsrc.py`, not the root `campsrc.py`.
2. Pipeline = **weld → ACVD remesh (density-scaled, ~100k-equivalent) → camp() → area-weighted
   variance of integrated_K**. Copy `run_branches_camp.py::_worker`, or use the `camp_delta_k()`
   wrapper above.
3. There is **no k**, CAMP is adaptive. Forget `k = A_T/ℓ²`; that's Khabaz, not us.
4. δK is a **variance of integrated K (K_j = ⟨kg⟩_j·A_j), area-weighted, dimensionless.**
5. **Ostia stay open; cut extent is standardised and part of the method; scale face count to
   surface area; queue the job.**
6. For the viewer, `camp()` hands back `face_labels` + `integrated_K` + a ready plasma `figure` -
   colour patches straight onto the remeshed mesh.
