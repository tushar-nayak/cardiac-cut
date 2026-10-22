# CAMP: what was changed here, and what was measured

Four modifications to `campsrc.py`, made 19–20 August 2026, plus a characterisation of
δK's reliability that the modifications made possible. **All four default to the
previous behaviour or are verified byte-identical against it**, so nothing here changes
a number that was not meant to change.

Written for the lab's CAMP work. Everything below is reproducible from files in this
repository; the command that produced each figure is named beside it.

---

> **STATUS: this is the canonical CAMP.** Decided by CAMP's author, 31 August 2026.
> This lineage — `campsrc.py` beside this file, with the four changes below — is the
> one to build on, cite and measure with. Where it disagrees with any other copy,
> including `Y:\aorta\python_code\camp`, this one wins. The older share copy still
> holds a parameter-sweep harness, `mesh_volume_mm3` and an `n_jobs` runner that are
> not here; if those are wanted, port them **forward onto this** rather than the other
> way round. See `README.md`.

## Summary

| # | change | default | effect |
|---|---|---|---|
| 1 | connectivity split | **on** | 0 disconnected patches, was up to 85% of a case's δK |
| 2 | fan gate, `fan_thresh=90` | **on** | wrapped patches split; SA21_1 fell 18× |
| 3 | `orphan_adjacency` | **off** | surface-aware orphan re-homing; tested, no benefit |
| 4 | `seed` parameter | **42** | makes the partition reproducible *and* averageable |

Verification that the defaults are unchanged: `seed=42` with flags 1–2 on and 3 off
reproduces `tools/camp_nub_bothfix_100k.jsonl` exactly (`DJ18_2` 0.008699, `SA21_1`
0.003754, `KK20_2` 0.022496, all to 1e-12).

---

## 1. Patches were coming out in two disconnected pieces

CAMP seeds with k-means on vertex **coordinates**. That has no notion of the surface,
so two sheets close in space but far apart *across* the surface land in the same
cluster. At the base of a branch, the branch wall and the aortic wall beside it are
~1.3 mm apart through the air and half a circumference apart along the surface.

A Monge patch fitted across that gap is not describing a neighbourhood of anything.

**Fix.** Split every cluster into its connected components before fitting, at vertex
level (`scipy.sparse.csgraph.connected_components` over the edge graph) and again at
face level after the partition, promoting large stray components to their own patch
with their own fit and merging sub-`min_patch_size` ones into whichever neighbour they
share the most edges with.

**Measured.** 2.3% of patches were affected on the arch cohort. On the 30-case nub
cohort at 100k faces the fix took disconnected patches to **0/30 cases**, and on
individual cases the affected patches had been carrying **up to 85%** of that case's
δK. `KY160_1`, a normal reading 3.85 — inside the diseased range — corrected to 1.45.

---

## 2. Patches were wrapping around the surface

A patch whose normals span more than 90° has folded around something. No height field
can describe it, so the Monge fit — which *is* a height field over a tangent plane — is
meaningless there. **R² and η can both pass on such a patch**, which is why it needed
its own gate: a quadratic fitted to a folded-over neighbourhood can look like a good
fit by residual while describing nothing.

**Fix.** `fan` = the largest angle between any vertex normal and the patch's mean
normal. `patch_needs_split` fails a patch whose fan exceeds `fan_thresh` (90°),
whatever its R² and η.

The gate alone was not enough. The splitter clusters in the patch's tangent plane
`(x, y)`, and a wrapped patch's two sheets lie *on top of each other* in that
projection — k-means cuts the overlapping disc into pie slices and every child
inherits both sheets, splitting forever without ever separating them. So the split
feature space gained height and the surface normal, `(x, y, z, s·nx, s·ny, s·nz)`,
scaled by the patch's own extent. On a wrap the two sheets face opposite ways, and the
normals are what tells them apart.

**Measured.** `SA21_1` (paediatric) 0.0671 → 0.0038, an **18×** drop, moving it to the
bottom of the cohort where a paediatric aorta belongs. `KY7_1`, a control with no
wrapped patches, barely moved (0.0057 → 0.0040). Cohort effect below.

### Known limitation

The gate judges **vertex** sets at fit time, but δK is computed from **face** labels
assigned afterwards by `face_patch_labels_fn`. A patch's face set reaches a ring beyond
its vertex set, so at a branch ostium it can pick up faces pointing the other way. The
patch that was *fitted* is a clean 27–62° disc; the patch that is *drawn and measured*
can span the junction at 94–137°.

Measured on `DJ18_2`: 0 patches wrap at fit time, 7 wrap on the final labels. On the
arch cohort, 58 such patches across 21 of 30 cases, carrying a median **0.5%** of δK
(worst case `KR6_1`, 42.7%). Removing them makes the cohort *worse* on every measure,
so they are left in place and recorded here rather than trimmed.

---

## 3. Combined effect of 1 and 2

30 nub geometries, 100,000 faces, **trim 0** — no trimming at all:

| | normals band | CV | diseased median | below normals p75 |
|---|---|---|---|---|
| before | 0.14 – 3.47 | 1.04 | 1.89 | 8 of 18 |
| connectivity only | 0.21 – 2.53 | 0.78 | 2.44 | 5 of 18 |
| **both** | **0.66 – 1.37** | **0.23** | 1.64 | 2 of 18 |

For comparison, the previous reference trim of 0.0025 with only the connectivity fix
gave CV 0.30. **Both fixes at trim 0 beat that**: δK no longer needs the trim to
suppress a structural defect, because the defect is gone.

Reproduce: `tools/camp_cohort.py`, plotted by `tools/nub_bothfix_trim0_labelled.png`.

### The fixes are specific: the unbranched cohort barely moves

Both changes were characterised on branched geometry, where they are enormous. The
obvious follow-up question — do they also rescale everything else? — is answered no,
and that is a result rather than a reassurance.

The same 30 patients cut as plain open tubes, run twice at 100k faces with the **same
seed 42**, changing only the two flags. Pairing on the seed cancels the partition
noise, so what is left is the flags alone.

| | branched, fan gate alone | unbranched, both fixes |
|---|---|---|
| median ratio | **0.234** | **0.980** |
| cases moving > 2% | 30 of 30 | 18 of 30 |
| cases moving > 10% | 28 of 30 | **1 of 30** |
| worst | 97% (`KK15_1`) | 13.4% (`KR47_1`) |
| rank ρ before vs after | — | **+0.996** |

The cohort keeps its shape exactly: the normals band goes 0.67–1.72 → 0.64–1.79, CV
0.386 → 0.394. **It does not tighten**, and that is the tell — on branched geometry
these same fixes took the normals CV from 1.04 to 0.23 because there was a defect to
remove. On an unbranched tube there is not one. Both defects require two surfaces
close through the air but far apart *across* the surface, or a patch wrapped around a
junction. That is a branch ostium, and these meshes have none.

For scale, the median 2% shift is well below the **±16%** that re-seeding the same
mesh produces (§4). The fixes move an unbranched aorta less than the partition noise
does.

Reported rather than rounded away: the normal-vs-diseased separation dips 4.02× →
3.79×, about 6% relative, in the same direction as the median shift. And this is one
partition per condition, which is right for isolating the flags but means these
absolute values should not be quoted against the 8-seed figures elsewhere in this
document.

**Consequence.** Unbranched δK measured before these fixes remains comparable with
unbranched δK measured after them. Branched δK does not — see §3 above.

Reproduce: `tools/unbranched_fix_effect.py`, plotted by
`tools/unbranched_fix_effect.png`.

---

## 4. One partition is one draw

This is the most consequential finding, and it is a property of CAMP rather than a bug
in it.

δK is carried by a **small number of patches**. Over 30 cases at 100k faces the
effective count (inverse Herfindahl of the contribution shares) has a median of
**22.5 out of ~867 patches**; the top 20 patches — 2.3% of them — carry a median
**58%** of δK. Those patches sit where curvature changes fastest, which on a branched
aorta means the branch ostia: 0.5% of faces hold ~85% of the seed-to-seed variance,
and 92–99% of it lies within 20 mm of an ostium. The inlet and outlet rims contribute
**0.0–0.1%** — the cut planes are not implicated.

Where k-means puts a boundary in those few places decides the number. So:

```
re-running the same mesh with a different seed moves δK by ~16%
```

**The draws are independent.** Over 18 case×density cells, the spread of a 4-seed mean
matched `CV/√4` with a median ratio of **1.03** (range 1.01–1.13), tested on disjoint
4-vs-4 halves. Averaging therefore delivers exactly what the arithmetic promises:

| partitions | error |
|---|---|
| 1 | ±16% |
| 4 | ±8% |
| 8 | ±6% |
| 16 | ±4% |

`seed` is now a parameter (`PARAMS["seed"]`, default 42, and `seed=42` reproduces every
historic number byte for byte). `shape_metrics.delta_k_averaged(mesh, seeds=...)` runs
one partition per seed and returns the mean with `k_fluct_se`, `k_fluct_sd`,
`k_fluct_cv` and `k_fluct_all`. **The arrays it returns are the first seed's** — patch
labels cannot be averaged across partitions that carve the surface differently, so the
picture is one partition and the number is the mean of all of them.

Reproduce: `tools/cohort_seeds.py`, `tools/seed_average.py`, `tools/seed_hotspots.py`.

### What averaging does *not* fix

δK drifts **systematically** with mesh density, and averaging cannot touch it because
it is not noise. `min_patch_size` is a **vertex** count, so the smallest allowed patch
covers ~16 mm² at 100k faces and ~4 mm² at 400k — different physical scales, not repeat
draws of the same measurement. Four of six cases drift monotonically downward as
density rises.

| | 1 seed | mean of 8 seeds |
|---|---|---|
| CV across 100k/200k/400k | 0.22 | 0.21 |

**δK values are comparable only at a fixed `target_faces`.**

---

## 5. Operating point

**100,000 faces**, chosen on evidence rather than convention, over six densities from
25k to 400k (`tools/camp_robust.py`, plotted by `tools/camp_robust_plot.py`):

| faces | normals CV | % of vertex ceiling | rank ρ vs next |
|---|---|---|---|
| 25,000 | 0.35 | 68.7% | +0.61 |
| 50,000 | 0.42 | 57.0% | +0.81 |
| 75,000 | 0.32 | 40.0% | +0.66 |
| **100,000** | **0.23** | **31.1%** | **+0.91** |
| 200,000 | 0.34 | 17.8% | +0.87 |
| 300,000 | 0.34 | 13.4% | +0.90 |
| 400,000 | 0.40 | 10.6% | — |

Below 100k the partition is limited by the `min_patch_size` floor rather than by
geometry, and agreement between neighbouring densities is erratic (+0.61, +0.81,
+0.66). From 100k up it is steady at +0.87 to +0.91 and **never converges further** —
more mesh does not buy more agreement. 100k also has the tightest normals band.

**Jitter.** Independent per-vertex Gaussian displacement of the source surface, mesh
held at 100k (`tools/camp_robust.py --mode sweep`):

| σ | rank ρ vs unperturbed |
|---|---|
| 0.05 mm | +0.77 |
| 0.10 mm | +0.73 |
| 0.20 mm | −0.20 |
| 0.40 mm | −0.33 |

It holds to 0.1 mm and collapses between 0.1 and 0.2 mm. **This is a worst case**:
white noise is the harshest possible model of segmentation disagreement, which in
reality is spatially correlated and much smoother. Treat 0.1 mm as a lower bound on
robustness, not an estimate of it.

**What δK supports.** With 8-partition averaging the median SE is ±6%. Ranking 30 cases
across a 0.68–2.5 range puts neighbours about 6% apart, so **only 1 of 29
adjacent-in-rank pairs is resolvable at 2 SE**. δK distinguishes populations and large
individual differences; it does not support fine ordering between neighbouring cases.

---

## 6. Branched vs unbranched: the branch offset

The same 30 patients, cut both with branch nubs and as plain open tubes, 100k faces,
8 partitions each (`tools/cohort_seeds.py`, plotted by
`tools/branched_vs_unbranched_compare.png`):

```
branched δK = 0.90 × unbranched δK + 0.00384        R² 0.828
```

Slope ≈ 1 with a large positive intercept: **the branches add rather than scale.** The
offset is essentially independent of the aorta (ρ with the unbranched value = −0.23),
and at 0.00384 it is **twice the median unbranched δK of a normal aorta** (0.00191).

| | area | share of branched δK |
|---|---|---|
| the three nubs | 3.9% of the surface | **40%** (median) |

A **10× over-representation** — and it falls unevenly. For a *normal* aorta the nubs
are **50–80%** of the entire measurement; for a diseased one, 13–30%. So the offset
inflates normals far more than diseased, and the discrimination suffers:

| | normals | diseased | separation | p |
|---|---|---|---|---|
| branched | 1.02 | 1.68 | 1.7× | 0.0003 |
| branched − offset | 1.04 | 2.82 | 2.7× | 0.0003 |
| unbranched | 0.78 | 2.77 | **3.6×** | 0.00004 |

**Consequences.**

- Branched and unbranched δK are **different quantities**. Do not pool them or compare
  across studies that differ in this.
- For normal-vs-diseased discrimination the unbranched form is **better**.
- If the object of study *is* branch anatomy, the offset is the subject rather than
  contamination.

Two caveats. The unbranched diseased are also more spread (≈1.2–4.7 against 1.0–2.5),
and nothing here separates real biological variation the branches were masking from
amplified noise. And two cases (`KK32_0`, `KK24_1`) have *negative* nub area — their
unbranched mesh is larger — so those two were not cut at matching stations and their
offsets are not purely "the nubs".

---

## 7. What was tried and did not work

Recorded because a write-up needs the negative results, and because each cost real time.

| attempt | result |
|---|---|
| divide δK by mean patch area² | over-corrects; normals ρ vs area +0.87 → **−0.84** |
| hold patch **area** constant across cases | reshuffles the cohort (rank ρ +0.23); rejected |
| scale patch count as `N ∝ area^α` | rejected on principle — count must follow geometry |
| `orphan_adjacency` (surface-aware orphans) | ρ(100k,200k) +0.906 → **+0.896**; left off |
| `min_patch_size` 18 → 12, 9, 6 | non-monotonic and within noise; 6 destroys the fit |
| post-remesh fold repair, 90/75/60/45° | **every threshold worse than doing nothing** |
| averaging partitions to fix density drift | no effect (CV 0.22 → 0.21) — it is not noise |

The fold-repair result is worth expanding, because the hypothesis was reasonable and
the data killed it. The remesh *does* reintroduce folds that the repaired source did
not have (`KK28_11` 0 → 2 at 100k → 6 at 200k), fold and sharp-crease counts *do*
correlate with per-case instability (ρ +0.56 to +0.61), and the mesh is otherwise
clean — 0 non-manifold edges, 0 duplicate faces, 0 unreferenced vertices across 90
scans. But repairing them made the error bars **wider** at every threshold, because the
repair Laplacian-smooths exactly where curvature is sharpest, which is exactly where δK
comes from. The defects are a **marker** of sharp anatomy, not its cause.

---

## Reproducing any of this

```
tools/camp_cohort.py          the cohort at one setting
tools/unbranched_fix_effect.py  the two fixes off vs on, unbranched, paired seed
tools/cohort_seeds.py         a cohort × n partitions          [--dir --suffix --tag]
tools/camp_robust.py          density / jitter / seed sweeps   [--mode ...]
tools/camp_robust_plot.py     the sweep figures
tools/seed_average.py         the 1/sqrt(n) verification
tools/seed_hotspots.py        where the seed variance sits, per face
tools/mesh_defects.py         defect inventory of the SOURCE surfaces
tools/remesh_defects.py       the same, on the REMESHED surface at each density
tools/arch_lock.py            freeze a geometry set + SHA-256 manifest
camp_viewer.py --patches      integrated-K surface, per case
pair_viewer.py                branched and unbranched side by side
```

`CURVATURE.md` holds the earlier characterisation. **Its `REPEATABILITY` table is
superseded** for anything run with changes 1 and 2: those numbers were measured before
both fixes and mix partition noise with defects since corrected, making them far more
pessimistic than the pipeline now is. `shape_metrics.PARTITION_CV` and `DENSITY_CV`
carry the current figures.
