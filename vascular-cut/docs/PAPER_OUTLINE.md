# Paper outline — VascularCut: a browser tool for measurable vascular surfaces

**Target: *Scientific Reports* (Nature Portfolio), as a tool / resource paper.**
Section order there is **Introduction → Results → Discussion → Methods** — Methods come
*last*, so Results must read before the reader has been told how anything was made.

Tool papers sit squarely in Sci Reports' scope, and their "importance is not assessed"
policy stops being a shield: for a tool, what matters is that it **works and is
available**, both of which are demonstrable.

Working file for drafting. Every number is in this repository and the source is named,
so nothing has to be re-derived while writing.

---

## The spine, in one sentence

**A browser tool that turns a closed vascular surface into an open one with clean
planar inlets and outlets — the form every downstream geometric measurement needs — and
a demonstration that what it produces is trustworthy enough to measure curvature on,
including at branch junctions where that had not previously been possible.**

The product is the **website**. The δK work is not the headline; it is the evidence
that the tool's output means something.

> **Hard requirement: the site must be live and usable by a stranger at submission.**
> Reviewers of a tool paper will click the link. Today they can get data neither *in*
> (upload is unbuilt; the app takes a server-side filesystem path) nor *out* (export
> writes to the server's own `output/`). See "Blocking work" at the foot of this file.

---

## 1. Introduction

- Geometric analysis of vessels needs an **open** surface with planar openings. Getting
  there from a segmentation is manual, slow, and inconsistent between operators — and
  the inconsistency propagates into every number measured afterwards.
- Branch junctions are the hard part and the interesting part: you cannot answer
  questions about branch anatomy on a surface with the branches removed.
- What this provides: a tool anyone can use in a browser, the same geometry engine as a
  desktop application, and a characterisation of how far the resulting measurements can
  be trusted.

---

## 2. Results

### 2.1 The tool, and what a user does

The workflow in one figure: load a surface and a centerline → openings are found,
classified and named → adjust or accept the cut stations → isolate a region or trim a
branch group → export an open mesh with a sidecar recording every plane.

Points worth making concrete rather than listing features:
- **Anatomical names are inferred**, not typed, so region presets follow automatically
  (aortic neck, infrarenal aorta, infrarenal aorta with iliacs, arch).
- **Every cut is validated before it is handed back.** Topology alone is not enough — a
  cut that leaks into the parent vessel still leaves one boundary loop and the Euler
  number stays right. Openings are checked against their own cut plane.
- **The tool refuses rather than producing quietly bad output**, and a refusal carries a
  workable alternative offset rather than "try something else".

### 2.2 It works on real cohorts — including in someone else's hands

The strongest evidence in the paper, and it is not ours: **an independent user cut a
49-case abdominal cohort, 48 of 49 successfully.** The single failure, `JK25_1`, is a
documented anatomical edge case (an aneurysmal common iliac makes the aortic
bifurcation divide so unevenly that the family call goes the wrong way; it also names no
renal in any recorded run). Their reported problems were otherwise cosmetic and are
fixed.

Automated benchmarks, most recent runs:

| | |
|---|---|
| abdominal anatomy, 70 cases | family correct **69/70**, renals named **67**, 4+ region presets **62**, presets offered **320**, openings named **663** |
| thoracic arches | **30/30** cut with the arch preset, each with recorded cut planes |
| cutter probe | 432 cut/offset combinations scored for rim flatness, topology and foreign surface |
| test suite | **50 invocations, all green** |

*Sources: `tools/anat_bench_after.log`, `output/arch_final/`, `tools/cut_bench.py`.*

### 2.3 Anatomical variation is handled automatically

- **Family call** — thoracic vs abdominal decided on whether the trunk *ends by
  dividing*, not on calibre. The old calibre rule placed 58 of 74 correctly and failed
  in both directions; the current rule scores 77/79. Where it is close, the tool says so
  on screen.
- **Bovine arch detected structurally** — `bovine = len(kids) == 3`, from the
  brachiocephalic feeding three vessels. **8 of 30**, no human labelling, independent of
  the anatomical naming.
- **Same representation from different anatomy:** a first-order vessel feeding three or
  more hands its stubs to its daughters, so bovine and normal arches both yield three
  stubs. Before that rule, six bovine cases fell through a fallback while two kept their
  stem and produced *two* stubs — same anatomy, different geometry, decided by an
  accident of stem length.

**Reported as a null:** bovine vs non-bovine δK, median 0.553 vs 0.634, **p = 0.60**,
n = 8 vs 22. Underpowered. Report it as a null; imply a trend and the section will be
taken apart.

**Limits to state:** the test recognises only `[3,1]` topology — `KK28_11` is `[2,2]`,
`KP35_1` has a fourth first-order branch, neither bovine by the test and **neither
confirmed by eye**.

*Sources: `tools/ANATOMY.md`, `tools/NUBS.md`.*

### 2.4 What the tool guarantees about the surface, and records about itself

This is what separates it from "a script that cuts meshes", and it is what makes §2.5
possible.

- **Openings are planar or the cut is refused.** Rim deviation is measured against the
  cut's own plane. Above 3 mm it refuses; below 0.5 mm it is clean; in between the cut
  is returned *and the deviation is recorded per opening* (`rim_dev_mm`, `rim_flawed`).
  This matters downstream: δK is an area-weighted variance of curvature, and a torn rim
  is real area with curvature the anatomy does not have.
- **Provenance in the sidecar**: every cut plane's origin and normal, which names were
  inferred versus hand-declared, and any plane the tool moved and why. An inferred side
  rests on the coordinate frame; a declared one does not, and a reader must be able to
  tell them apart.
- **The tool says when it does not know.** Left/right is decided for a whole case from a
  single anatomical witness (the splenic/hepatic pair; the brachiocephalic). Where no
  witness exists it falls back to the coordinate frame — **11 of 59 abdominal cases** —
  and now says so rather than presenting a guess identically to a reading.

### 2.5 Validating the downstream measurement: δK on branched geometry

The tool exists so measurements can be made. This is the evidence that they can.

**Two defects in the curvature partition, both specific to branched geometry.**
A partition that clusters on spatial **coordinates** has no notion of the surface: a
branch wall and the aortic wall beside it are ~1.3 mm apart through the air and half a
circumference apart across the surface. And a patch straddling an ostium *wraps* — no
height field describes it, and **R² and η both pass** on it, so goodness of fit is not a
validity check.

| 30 branched aortas, 100k faces, trim 0 | normals band | CV | diseased below normals p75 |
|---|---|---|---|
| before | 0.14 – 3.47 | 1.04 | 8 of 18 |
| **both fixes** | **0.66 – 1.37** | **0.23** | **2 of 18** |

Complete overlap before; separation after. `KY160_1`, a normal reading 3.85 — inside the
diseased range — corrects to 1.45.

**The fixes are specific, not a rescaling.** Paired design, same patients cut as plain
tubes, same seed, only the flags changed: branched median ratio **0.234**, unbranched
**0.980**, rank ρ **+0.996**. Unbranched CV goes 0.386 → 0.394 — *it does not tighten*,
because there was no defect to remove.

**One partition is one draw.** δK is carried by ~22 effective patches of ~850, sited at
the ostia: 0.5% of faces hold ~85% of the seed-to-seed variance. Re-seeding moves δK
**~16%**. The draws are independent (4-seed mean matched CV/√4, median ratio 1.03), so
averaging gives **±6% at 8 partitions**. *The inlet and outlet rims contribute 0.0–0.1%
— the cut planes are not implicated*, which pre-empts the obvious objection.

**Comparable only at fixed mesh density** — ~21% drift across a 4× change, because
`min_patch_size` is a vertex count. Operating point 100k faces, chosen across seven
densities; agreement plateaus at ρ +0.87–0.91 and never converges further.

*Sources: `camp/CHANGES.md` §1–5.*

### 2.6 What the tool enables: the branch contribution

```
branched δK = 0.90 × unbranched δK + 0.00384        R² 0.828
```

Slope ≈ 1 with a large intercept: **branches add rather than scale.** The three nubs are
**3.9% of the surface and 40% of δK** — a 10× over-representation, and the same fact
that explains why branched geometry was where the defects bit.

| | normals | diseased | separation | p |
|---|---|---|---|---|
| branched | 1.02 | 1.68 | 1.7× | 0.0003 |
| branched − offset | 1.04 | 2.82 | 2.7× | 0.0003 |
| unbranched | 0.78 | 2.77 | 3.6× | 0.00004 |

> **Reviewer will ask why not strip the branches.** Three-part answer, all needed:
> the compression is *explained* (the offset is 50–80% of a normal's δK, 13–30% of a
> diseased one, so it squeezes the gap); subtracting it recovers **2.7×**, proving
> compression not degradation; and **you cannot study branch anatomy on a surface with
> the branches removed** — unbranched separates better because it discarded the subject.

### 2.7 Availability

Browser tool (URL), desktop application, and source. Same geometry engine behind both.

---

## 3. Discussion

- What the tool is for and what it is not: it produces surfaces for measurement, it does
  not diagnose.
- What δK supports — populations and large individual differences — and what it does
  not: across 30 cases spanning 0.68–2.5, only **1 of 29 adjacent-in-rank pairs is
  resolvable at 2 SE**. State plainly; it is a strength.
- Branched and unbranched are **different quantities**. Do not pool them.
- **Prior δK values measured before these fixes are affected.** `KY160_1` moved
  3.85 → 1.45, crossing a normal/diseased boundary. Handle deliberately — "we found the
  defect, quantified it, report corrected values" — rather than letting a reader find it.
- Generalisable beyond this tool: coordinate-based clustering fails at any branch
  junction; a stochastic partition makes a shape metric a draw rather than a value; and
  **goodness of fit does not establish validity** for a manifold fit.
- Limitations: n=30 for the curvature cohort, single site, single segmenter; jitter
  robustness holds to 0.1 mm and collapses by 0.2 mm (white noise is the worst case, so
  a lower bound not an estimate); the arch cohort never reached the nub cohort's density
  stability (best ρ +0.60 vs +0.87–0.91), which is partly why the bovine test is a null.

### Negative results — a subsection, load-bearing

| attempt | result |
|---|---|
| divide δK by mean patch area² | over-corrects; normals ρ vs area +0.87 → **−0.84** |
| hold patch **area** constant | reshuffles the cohort (ρ +0.23); rejected |
| surface-aware orphan re-homing | ρ(100k,200k) +0.906 → **+0.896**; left off |
| `min_patch_size` 18 → 12, 9, 6 | non-monotonic, within noise; 6 destroys the fit |
| post-remesh fold repair, 4 thresholds | **every one worse than doing nothing** |
| averaging partitions to fix density drift | no effect — it is not noise |

Expand the fold one: folds *do* correlate with instability (ρ +0.56–0.61), but repair
smooths exactly where curvature is sharpest, which is exactly where δK comes from. The
defects are a **marker** of sharp anatomy, not its cause.

---

## 4. Methods *(last, per Scientific Reports)*

**4.1 Segmentation and centerlines** — collaborator's text, edit and trim.

**4.2 The cutting engine** — plane placement, cut modes, the obliquity and tear gates,
validation, the region and nub-in presets. `n = s_wall + k·(D_aorta/2)` with `s_wall`
**measured**, not assumed (assuming is wrong by up to 4.7 mm on oblique take-offs —
enough to put the plane inside the aortic lumen).

**4.3 Web application** — architecture, the shared API layer, what runs server-side, and
the deployment.

**4.4 Remeshing — do not omit.** ACVD to a fixed 100,000 faces. This is *why* δK is
comparable only at a fixed face count.

**4.5 CAMP and δK** — the estimator, unchanged and locked. Then the four partition
changes, each defaulting to previous behaviour, and the equivalence proof: with the
additions off, identical to the July 2026 upstream (0.397519623517 both ways on KK2_5).

**4.6 Reliability methods** — seed averaging, density sweep, jitter, the paired design.

**4.7 Statistics** — Mann–Whitney, Spearman, Fisher-z CIs. δK reported as a mean over 8
partitions with SE.

---

## 5. Declarations — start now, cannot be retrofitted

- **Ethics / IRB** and consent, for patient-derived scans. **Required.**
- **Data availability.**
- **Code availability** — `https://github.com/SurgBioMech/CAMP` (canonical lineage,
  settled 31 Aug 2026) plus the tool's own repository and the live URL.

---

---

## Related work — read this before drafting the Introduction

Checked 31 Aug 2026. **The field is more crowded than it looks, and three integrated
aorta platforms have appeared in the last eighteen months.** A reviewer will know them.
A claim that nothing like this exists is not available and would be caught.

### What exists

| tool | what it is | overlap with this work |
|---|---|---|
| **AortaAnalyzer** (Computers & Graphics, Nov 2025, open access, doi 10.1016/j.cag.2025.104415) | end-to-end CTA pipeline: 3D segmentation network (Dice 0.95), interactive correction, surface extraction, centerline computation, **inlet/outlet capping**, clinical metrics. Single GUI. Seven-expert evaluation; processing time benchmarked against 3D Slicer and SimVascular | **the closest.** Caps vessel ends where we open them |
| **AortaExplorer** (Med Biol Eng Comput, 2026) | AI-driven aortic analysis from CT | same neighbourhood, imaging-led |
| **AortaPlot** (medRxiv, Jan 2026) | desktop Python; 20 standardised landmarks, diameter (sliding-axis), arc length, regional volume, radiomics. Validated on ~556 studies | landmark profiling, **no curvature**, no branch naming |
| **VMTK** | the standard library: centerlines, `vmtkbranchclipper`, surface analysis | genuinely clips branched surfaces; a command-line library, no anatomical naming, no opening validation |
| **SimVascular** | full CFD pipeline; trimming gives "a plane or box to cut the model at an arbitrary location or perpendicular to a path" | the cut, driven manually |
| **svMorph** | interactive virtual sculpting of vascular geometries | editing/morphing, not measurement preparation |

### Where this work sits

Two searches returned **no browser-based tool** for any of this — every comparator is a
desktop GUI or a library. That gap looks real, but "we found nothing" is a weak claim on
its own and should not be the load-bearing one.

The defensible distinctions, in descending order of strength:

1. **A curvature metric with a characterised reliability.** The others quantify
   diameter, volume, landmarks and radiomics. None reports a shape metric's measurement
   error. δK's ±16% single-partition draw, the 1/√n averaging law and the fixed-density
   constraint are, as far as this search goes, without parallel in this space.
2. **Opening rather than capping.** AortaAnalyzer closes vessel ends so a solver has
   boundary conditions; this opens a closed surface with planar rims so it can be
   measured. Same geometric problem, opposite output, different purpose. Say it plainly.
3. **Branch-level anatomical naming driving region presets.** The others locate
   landmarks *along* the aorta; this names individual branches and builds anatomically
   defined regions from those names.
4. **Delivered in a browser**, no install.
5. **It refuses rather than producing quietly bad output**, and a refusal carries a
   workable alternative offset.

### The framing consequence — important

AortaAnalyzer occupies "integrated aorta analysis platform" convincingly, with expert
evaluation and timing comparisons against the obvious competitors. **Pitched as a
general aorta tool, this paper arrives second into a space someone has just filled
well.** Pitched as *a browser tool that makes this curvature metric usable and
trustworthy on branched geometry*, it is not competing with them at all, and every
distinctive result here — the two partition defects, the one-draw finding, the branch
offset — sits where none of those three papers goes.

Still tool-led. But the tool is **for the metric**, not a general platform.

### Two things worth copying from them

- **AortaAnalyzer evaluated with seven domain experts.** Cheap credibility for a tool
  paper, and it partly substitutes for the inter-operator study this work lacks.
- **They benchmarked processing time against named alternatives.** Worth doing, on the
  axis that suits: not "faster than Slicer" but how long a cut and a δK partition take.

### And a precedent for the data statement

AortaAnalyzer's own data availability reads: training data public on Figshare, "the
additional dataset used in this study is not publicly available due to patient privacy
concerns and institutional data use agreements" — in an open-access paper in a good
venue. That is the model for ours.

### Still to check before submission

Nobody has read AortaAnalyzer or AortaExplorer in full — only abstracts and highlights.
**Confirm whether either does branch-level naming or any curvature measure** before the
Introduction commits to distinction 1 or 3.

## Figures

| # | figure | file |
|---|---|---|
| 1 | the tool: workflow, and a screenshot of the browser app | *to make* |
| 2 | one case cut — before/after, openings named, region isolated | *to make* |
| 3 | the two partition defects, before/after on `KY160_1` | *to make from* `camp_viewer.py --patches` |
| 4 | cohort δK vs area, before and after both fixes | `tools/nub_bothfix_trim0_labelled.png` |
| 5 | specificity: unbranched barely moves | `tools/unbranched_fix_effect.png` |
| 6 | partition noise and 1/√n | `tools/seed_stability_summary.png`, `cohort_seedavg_SE.png` |
| 7 | density sweep | `tools/camp_robust_density.png` |
| 8 | branched vs unbranched, the offset | `tools/branched_vs_unbranched_compare.png` |
| 9 | arch cuts including bovine | `output/arch_final/`, `arch_viewer.py` |

**House style for the cohort plot** — asked for repeatedly, do not substitute: δK
against surface area, linear axes, one point per case coloured by group, **case names
labelled on every point**, no fitted regression lines.

---

## Blocking work — the site must be usable by a stranger

Re-triaged 31 Aug when the product was settled as the website. These were previously
recorded as deferred; that was right for a methods paper and wrong for this one.

1. **Uploads.** The app takes a filesystem path typed into a text box and the server
   reads it — arbitrary file read on the host. It cannot be safely hosted as it stands.
   This is both the security gate and what makes the tool usable by someone who is not
   on the machine. `tools/HOSTING.md` job 1.
2. **Download.** Export writes the STL and sidecar into the server's `output/`. A remote
   user currently has no way to retrieve their own result.
3. **Hosting** — a reachable URL, stable through review and after.

Nothing else in `tools/FREEZE.md`'s deferred list blocks submission.
