# CAMP

The Pocivavsek Lab's convergent adaptive patch partition, behind the δK (`k_fluct`)
shape metric — with three additions to the partition and one to reproducibility, made
while measuring 30 thoracic aortas in August 2026.

`CHANGES.md` is the full record: every change with its measured effect, the operating
point the measurements settled on, and the things that were tried and did not work.

---

## Provenance

`campsrc.py` here descends from the lab's copy taken **22 July 2026**, git blob SHA-1
`893d0d2d028ee9cd9274cd7f4b0a9ffc75376bd2`, preserved unmodified in
`reference/upstream_campsrc.py`. `reference/campsrc.patch` is the difference — 8 hunks,
+248/−32 lines.

### THIS IS THE CANONICAL VERSION — decided 31 August 2026

By CAMP's author. This file is the one to build on, cite, and measure with. Where it
disagrees with any other copy, **this one wins**.

That is a decision about which lineage carries forward, not a claim that this file
contains everything — it does not, and the difference still matters when porting.
Checked against `Y:\aorta\python_code\camp\campsrc.py` (20 July 2026) on 26 August,
the two had diverged in **both** directions:

| only in the older share copy | only here (canonical) |
|---|---|
| `mesh_volume_mm3` reporting | `max_children`, `child_rule`, `eta_thresh` in `patch_pass` |
| `n_jobs` parallelism | the adaptive `cap` children rule |
| `camp_param_sweep`, `_run_one_trim`, `_process_one_job` | the four additions above |

Our four additions were made on the `UpdatedCamp` copy of 22 July, so **the sweep
harness, the volume reporting and the parallel job runner are not in this file**, and
applying `campsrc.patch` to the older copy fails 3 of its 8 hunks.

**The direction of travel is now one-way.** If you need the sweep harness, the volume
reporting or the parallel runner, port them **forward onto this file** — do not port
these additions back onto the older copy and use that instead. They are independent of
one another and touch `patch_pass` and `camp` only, which is not where the sweep
harness lives, so the port is tractable. `verify_upstream_equivalence.py` checks the
result: with the additions switched off it must still reproduce
`reference/upstream_campsrc.py` to 1e-12, and that check is what every number in
`CHANGES.md` rests on.

---

## What is different

| addition | default | off switch |
|---|---|---|
| connected-component split before fitting | on | `split_disconnected=False` |
| fan gate on wrapped patches | `fan_thresh=90.0` | `fan_thresh=None` |
| surface-aware orphan re-homing | **off** | — (off is upstream) |
| `seed` parameter for k-means | `42` | `seed=42` is upstream |

**With all of them off, this returns exactly what upstream returns.** Not approximately
— exactly. `verify_upstream_equivalence.py` runs both copies on a surface you supply and
compares δK to 1e-12:

```
python verify_upstream_equivalence.py <surface.mat>
```

That check is the foundation for everything in `CHANGES.md`. If it ever fails, the
additions are not additions and comparisons to prior work are broken.

---

## The three additions, briefly

**A patch must be one connected piece of surface.** k-means clusters vertex
*coordinates*, which has no notion of the surface, so at a branch base the branch wall
and the aortic wall ~1.3 mm away through the air land in the same cluster and the Monge
fit spans a gap. Splitting clusters into connected components before fitting took
disconnected patches to **zero across 30 cases**; on affected cases those patches had
been carrying **up to 85%** of that case's δK.

**A patch whose normals span past 90° has wrapped.** A Monge patch *is* a height field
over a tangent plane, so no such fit describes it — but **R² and η both pass on one**,
because a quadratic fits the residuals fine while describing nothing. Hence a separate
gate. Splitting it also required the height and the surface normal in the clustering
features: the two sheets sit on top of each other in the tangent-plane projection, so
k-means was cutting the disc into pie slices and giving every child both sheets.

Wrapping is a **branch phenomenon** — zero occurrences on unbranched aortas, 0–6 patches
out of ~850 on branched ones. Rare, but not unimportant: SA21_1 has exactly two, and
gating them moved its δK by **18×**.

Together, at trim 0 with no trimming at all, the normals band goes from **CV 1.04 to
0.23** — better than the previous reference trim achieved.

**The k-means seed is a parameter.** δK is carried by a median **22 patches out of
~850** — 0.5% of faces hold 85% of the seed-to-seed variance, sited at the branch ostia
— so one partition is a draw rather than the answer. Re-seeding the same mesh moves δK
by **~16%**. The draws are independent (a 4-seed mean matches CV/√4, median ratio 1.03
over 18 case×density cells), so averaging *n* partitions divides the error by √n:

| partitions | error |
|---|---|
| 1 | ±16% |
| 4 | ±8% |
| 8 | ±6% |

---

## Two things to know before using this

**These change δK for every branched geometry.** Anything already published from this
code will move. Reproduce the old numbers with the switches above.

**Averaging does not make densities comparable.** δK drifts ~21% across a 4× change in
face count, because `min_patch_size` is a *vertex* count — the smallest patch covers
16 mm² at 100k faces and 4 mm² at 400k. Those are measurements at different scales, not
repeat draws, and averaging cannot remove it. Compare only at a fixed target face count.

---

## Not included

The averaging wrapper itself. `seed` is what makes averaging possible, but running *n*
partitions and returning mean ± SE lives in the calling project
(`shape_metrics.delta_k_averaged`). It is a dozen lines and deliberately outside CAMP.
