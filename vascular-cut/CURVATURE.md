# delta K: what the number is, and how far you can trust it

Everything here was measured on the 30 thoracic nub geometries, and every figure in
this document can be recomputed from files in this repository:

| what | where |
|---|---|
| the 180-run noise sweep | `tools/camp_noise.jsonl`, produced by `tools/camp_noise.py` |
| the reference cohort | `tools/camp_cohort_final_trim0p0025.csv`, produced by `tools/camp_cohort.py` |
| the mesh defect inventory | `tools/mesh_defects.py` |
| the figure sent to collaborators | `tools/camp_final.py` |

The operating point those measurements settled on lives in `shape_metrics.REFERENCE`
and is what both web front ends default to.

---

## The one thing to know first

**delta K is not a property of the surface alone.** It is a property of the surface
*and* of a partition that CAMP rebuilds from scratch on every run. Remesh the same
geometry to 98,000 faces instead of 100,000 -- the anatomy is identical, not one
vertex of the original is discarded -- and the partition comes out different, so the
number does too.

How different, over 30 cases at three face counts each:

| trim | typical spread | worst case | worst |
|---|---|---|---|
| **0** | **2.64x** | **23.48x** | DJ18_2 |
| 0.001 | 1.88x | 10.88x | KK16_1 |
| **0.0025** | **1.59x** | **3.81x** | KK16_1 |
| 0.005 | 1.44x | 2.74x | KY244_1 |
| 0.01 | 1.28x | 1.78x | KY244_1 |

Read that as a noise floor. Two delta K values that differ by less than the typical
spread are the same value. The web app prints this band beside the trim control and
again beneath the result, because the value on its own invites a comparison it cannot
support.

## Why the reference trim is 0.0025

Trim 0 is the untouched number and it is the one people ask for. It is also the least
reproducible setting on the table -- a factor of 2.6 typically, and on DJ18_2 a factor
of **23.5** between runs of the same surface. Nothing about the vessel changed.

The trim is a quantile on integrated K and it removes **both tails**, so the fraction
typed is per tail: 0.0025 discards 0.5% of patches, about 1% of area. That is enough
to cut the spread from 2.64x to 1.59x, and it is small enough that the shape of the
cohort curve does not move -- the ordering of cases is preserved and the separation
between the normal and dilated groups survives.

Trims of 0.005 and 0.01 are steadier still. They were not chosen because at that
point the discard stops being an outlier rule: 1% per tail is 2% of patches and
nearly 2% of area, and there is no longer a clean answer to "what did you remove and
why".

**0.05 was dropped from the front ends entirely.** A tenth of all patches is not a
trim, it is a different metric.

## Patches that fail to converge are never dropped

CAMP marks a patch as failed when its Monge fit misses `R2 >= 0.9` or `eta <= 0.075`.
Those patches are counted in delta K, always, and there is no option to exclude them.

This is deliberate and it is not a detail. A patch that will not take a quadratic is
precisely where the surface is least like a quadratic -- which is the thing delta K
exists to measure. Excluding them does not clean the number up; it removes most of
it. The trim is a quantile on the *value*, which is a different operation: it drops
the extremes wherever they came from, rather than dropping a category.

## Folded faces: a real defect, and not the explanation

A dihedral above 90 degrees means the surface doubles back through itself. That is
geometry the vessel does not have, and it arrives in the source segmentation -- it is
not introduced by remeshing or by the cut.

**14 of the 30 reference meshes carry folds.** The worst are KP35_1 (46 folded edges,
maximum dihedral 174.6 degrees, 28 of them inside one 8 mm region), KK40_3 (35),
KY160_1 (31), KR52_1 (20), KR6_1 (15), KR27_2 (12), KR47_1 (10). On KP35_1 the folded
region produced a single patch at integrated K = -41.7 carrying 54% of that case's
entire delta K, against a next-worst patch of -5.4.

It is very tempting to conclude that folds are what makes delta K unstable. **They
are not.**

- correlation between fold count and the 98k/100k/102k spread: **-0.07**
- the worst case for spread, DJ18_2 at 23.5x, has **zero** folded edges
- among the 16 fold-free meshes, the worst spread is that same 23.5x

So folds are reported as a property of the input -- both front ends warn before the
run and again with the result -- and are never applied as a correction. Fix them
upstream if the value matters.

### Repairing folds is not a default, and should not become one

`tools/camp_cohort.py --repair` smooths only the vertices caught in folded edges, two
rings out. Judged against each case's own noise band, over the 30 cases:

| verdict | n | cases |
|---|---|---|
| real improvement | 2 | KP35_1 (0.21x), KY244_1 (0.34x) |
| real worsening | 1 | KR73_1 (2.02x) |
| inside the noise | 11 | KK16_1, KK20_2, KK28_11, KK40_3, ... |
| no folds to repair | 16 | DJ18_2, DJ8_1, KK15_1, ... |

Two clear wins, one clear loss, and eleven cases where it changed nothing you could
distinguish from remeshing twice. That is not a repair worth applying silently.

### Two approaches that were tried and did not work

Recorded so they are not tried again:

- **Rim collar** -- excluding patches within a fixed distance of the open rim. Ate 5%
  of the area and left KP35_1 at 1.26; its normalised value got *worse*, 33x to 57x.
- **Fold masking** -- excluding only patches containing folded faces. A patch can be
  badly fitted without containing a fold, and folded faces sit in patches that are
  otherwise fine, so the mask both over- and under-reaches.

## The comparison settings

A delta K is only comparable to another delta K produced the same way. All of these
must match, and both front ends warn when they do not:

- **face count** -- 100,000. delta K moves with it; that is the whole point of the
  table above.
- **trim** -- 0.0025, both tails.
- **cut extent** -- delta K is extensive. A longer piece of aorta has a larger one.
- **topology** -- the reference cases are all at euler -3. A closed surface still
  partitions, so nothing fails; it answers a different question.
- **CAMP build** -- vendored under `camp/`, hash in `camp/README.md`. Three files
  named `campsrc.py` existed in the source project and they are **not** the same
  method. Whichever reached `sys.path` first would answer, and a wrong one returns a
  different number without complaining. `CAMP_SRC` overrides it deliberately.

The reference cohort at those settings: 30 cases, area 15,921 to 87,190 mm2, 9.8% to
29.1% of patches at the vertex floor, raw delta K 0.00325 to 0.03598.

## Normalisation

Raw delta K is reported. For the cohort figure it is also divided by the mean of the
normal groups (`KY` and `SA`, n = 12, mean 0.007232), which puts 1.0 at "normal" and
makes the groups readable on a linear axis:

| cohort | n | median normalised |
|---|---|---|
| SA | 1 | 0.524 |
| KY | 11 | 0.703 |
| KR | 5 | 1.824 |
| KK | 9 | 1.885 |
| KP | 2 | 2.186 |
| DJ | 2 | 2.415 |

The normal groups separate from the rest. **Bovine arch does not**: n = 9 bovine at
median 1.344 against n = 21 non-bovine at 1.537, which is well inside the noise band.
That is a negative result and it is worth keeping as one.

## Reproducing any of it

```
python tools/camp_cohort.py run  --indir output/nubs_final --faces 100000
python tools/camp_cohort.py cohort --trim 0.0025
python tools/camp_noise.py                      # 180 runs, about 40 minutes
python tools/mesh_defects.py --indir output/nubs_final
python tools/camp_final.py --by cohort          # the collaborator figure
```

Every export carries a sidecar recording face count, trim and CAMP build. A delta K
without those cannot be reproduced or compared, which is why it is not optional.
