# Frozen for the write-up — 28 August 2026

Development on the vascular cut tool and on CAMP is **paused** so the work can be
written up. This file records what state it was frozen in, what is deliberately left
open, and what each open item would need. It exists so that picking this up in three
months does not start with rediscovering all of it.

**What "frozen" means here.** No new features. Correct a fault that makes a published
number wrong, or that stops a collaborator working, and say so here. Everything else
waits.

## State at the freeze

| | |
|---|---|
| Test suite | **50 invocations, all green**, all returning a meaningful exit code |
| Desktop | `branch_gui.py`, the reference implementation |
| Web | `webapp/cut_app.py` (trame), `webapp/app.py` (Streamlit, superseded) |
| Release | `dist/vascular-cut-2026-08-28.zip`, verified to carry every fix below |
| Docs | `docs/Whats-New.pdf` rebuilt; `docs/User-Guide.pdf` has no source and is unchanged |
| CAMP | `camp/campsrc.py` + `camp/CHANGES.md`; also at `SurgBioMech/CAMP` and on the lab share |
| Benchmarks | abdominal 48/49 family, 146/182 end to end; thoracic 29/30, 87/99 |

Two things worth knowing about the suite. Until this week **17 of 41 scripts printed
`FAILURES` and exited 0**, so any runner read them as green — one assertion had been
silently broken for some time and one was passing vacuously. And two scripts take an
argument (`test_anatomy_gui.py` / `test_region_gui.py`, `jk5|ky160`); a runner that
calls them bare will report `IndexError` and that is the runner's fault, not theirs.

## Fixed in the last round

All from a collaborator who cut a 49-case abdominal cohort and reported four things.

- **An open-ended region closed itself when its start plane moved.** The End control is
  disabled for such a region, so its number is not the user's — and both the desktop and
  the web read it back anyway. `branch_gui._span_from_controls`, and
  `cut_app._sync_region_boxes` / `commit_region`. The iliacs went 4 → 0 on JK5_1.
- **The web could not select that preset at all** — `float(None)` raised, and the
  Streamlit app crashed on `sorted(s.region)` at the same preset.
- **Left/right is now visibly qualified.** One scalar decides side for a whole case; when
  no witness is found it falls back to the coordinate frame and used to say nothing.
  11 of 59 abdominal cases are in that position. Both front ends now warn, and both draw
  an L/P/S triad so a side can be checked against something.
- **`Whats-New.pdf` was quoting superseded δK noise figures** (2.6× typical, 23.5× worst)
  measured before the two CAMP partition fixes, and a fold correlation of −0.07 that is
  now +0.56 to +0.61. Corrected, and the document now tells readers to discard the
  figures in any copy printed earlier.
- **`HOSTING.md` said δK was dropped from the hosted build.** It is not, and the flag it
  describes as job 2 was never built.

## Open, and what each needs

### Needs the abdominal cohort on the share drive

1. ~~**Iliac sides are judged against the wrong midline.**~~ **DONE, 28 Aug** — measured
   over 70 cases before and after in the same harness. `side_of()` compared each limb
   to `mid_x`, a median taken over an aorta path that *includes one of the limbs*; it
   now judges the two limbs against each other, which needs no midline. **Openings
   named 645 → 663, everything else identical, 8 cases gained, 0 lost, every gained
   name an iliac.** JK36_2 recovered all six (the `_unique` double-withdrawal); JK17_1
   gained the right iliacs, the case the lab's QC log flags as "did not detect right
   iliac correctly". Cases naming both common iliacs: 54 → 62 of 70. Logs in
   `tools/anat_bench_baseline.log` and `tools/anat_bench_after.log`.
2. ~~**Two thin-branch cuts still owed**~~ — **INVESTIGATED 28 Aug, no defect.** JK2_1
   Branch 6 at 80% and JK24_1 Branch 4 at 50%, both very thin (r 1.8 and 1.5 mm).
   - The band they sit in is **documented policy**, not an oversight: above
     `TORN_REFUSE_MM` (3 mm) a cut is refused, below `RIM_FLAW_MM` (0.5) it is clean,
     and in between it is handed back deliberately — measured over 2,927 cuts on 65
     aortas, see the comment at `TORN_REFUSE_MM`.
   - `suggest_pullback` finds a clean cut for **every one** of them: 0.6 mm back on
     JK2_1, 0.8 mm back on JK24_1.
   - And the user-facing paths never reach the flawed cut at all. `flat_trunk_offset`
     backs a bad cut off automatically and records the move in the sidecar as
     `adjustments.vessel_end`; on JK24_1 the session pulled 10.8 → 8.1 mm by itself.
     `cut_bench` probes raw `cut_terminal` at fixed fractions, which is harsher than
     the tool as used, so those two rows were the benchmark being strict.

   One genuine gap remains and is **deferred**: `validate_cut` reports a flawed rim
   clearly but says "Move this cut" without a number, while the refusal path appends a
   workable offset. Closing that needs the base mesh passed into `validate_cut`, which
   is a signature change — after the freeze.

### Needs one fact from the person who reported it

3. **Which mechanism actually crossed the iliacs.** Side is one global scalar, so an
   inverted sign inverts *everything*. **Were the renals crossed too?** Both crossed → it
   is `left_sign` (the silent frame fallback, or a wrong splenic/hepatic pick). Renals
   right, iliacs wrong → it is not `left_sign`, and the likely cause is reading a side off
   an anterior view, where patient-left sits on the viewer's right. The warning and the
   triad shipped because they are correct either way; the root fix waits on this.

   Note the midline fix above does **not** close this. A biased `mid_x` cannot produce a
   swap: the two limbs straddle the true midline, so shifting it either leaves them
   correctly opposite or pushes both onto the same side, which is the withdrawal case.
   A swap needs `left_sign` wrong, and the pair-relative test still multiplies by it.

### Known desktop/web gaps, from the 28 Aug parity audit

Verified against source. None is a regression; all predate this round.

**Fixed on 28 Aug, each with a test** — every one was a control that reported success
and did nothing, which is the same fault class the cohort reported:

- ~~Tilt / shift / plane size / reset on a **region's two boundary cuts** wrote to
  throwaway `Terminal` objects.~~ `CutSession._bound_planes` now carries them, as the
  window has since `test_region_planes.py`.
- ~~Trims offered for branch groups **outside** the region, recorded and then silently
  dropped by `_trims_within`.~~ `trim_candidates` is scoped to the region.
- ~~Send-cut-to-mark could select a region boundary and report success having moved
  nothing.~~ The guard tested `is_region_bound`, which **nothing ever set**; `region_plan`
  now sets it where boundaries are made.
- ~~The API sidecar dropped `names` and `named_by_hand`.~~ Both recorded. This mattered
  most for provenance: an inferred side rests on the coordinate frame and a declared one
  does not, and a reader had no way to tell them apart.
- ~~Web `Reset` left trims and hand-declared arteries standing.~~ Scope now lives in
  `CutSession.reset_all`, called by both front ends so they cannot drift again.
- ~~A hand-set family carried into the next case on the web.~~ Resets on a path change,
  as `branch_gui._reset_family` does.

**BLOCKING THE PAPER — re-triaged 31 Aug 2026**

The product of this work is **the website** (`tools/HOSTING.md`: "The web tool is for a
paper: readers upload a surface and a centerline, get a cut mesh back"). The paper is
tool-led, so **reviewers will open the link and try to use it**. Two things stop them,
and both were previously filed here as deliberately deferred. That was right for a
methods paper and wrong for this one.

- ~~**No upload.**~~ **Built 24 Sep 2026.** Public mode has no path box; files come in
  over HTTP, are checked (bytes, faces, format) before anything is kept, and are stored
  under a generated id with the notice shown first. `HOSTING.md` jobs 1, 4 and 7.
- ~~**No download.**~~ **Built 24 Sep 2026.** The cut and its sidecar come back to the
  browser as a zip through a one-time link.
- **Visitors kept apart, 25 Sep 2026.** One visitor per process in public mode, wiped
  when they leave (`HOSTING.md` job 5). What still stands between this and a public
  link is in `HOSTING.md`: a container that renders headless, the notice's last two
  blanks, institutional sign-off, and enough processes for the visitors expected.

**The freeze is lifted for `webapp/` and `vascular_api.py` only.** Everything else stays
frozen: `branch_cut.py` in particular is untouched without a `cut_bench` measurement
either side. Whoever does this work should update this section as it lands.

**Still open — behaviour drift**

- No live cross-section readout while dragging an offset; region discs do not follow a
  slider until release.
- Export is STL only on the web; the desktop offers STL / PLY / OBJ.

**Cosmetic**

- Boundary-loop colouring: the web never passes the cut modes, and the mapping in
  `scene.py` is inverted. Plane discs and opening rows are not branch-coloured.
- Pinned marks have no numbered labels, a fixed rather than model-scaled radius, and no
  depth bias, so they can be hidden inside the mesh.
- Marks and the measure disc stay on screen over the CAMP partition in the shape view.
- Base mesh colour differs from the desktop and flips when branch colouring is toggled.
- The two web front ends offer different partition counts.

### Known and settled — do not reopen

- **The family rule is at its ceiling.** Thoracic reaches 0.249 and abdominal falls to
  0.225 on smoothed centerlines, so JK25_1 (0.205) and KP35_1 (0.309) are on the wrong
  side of *every* cut-off relative to each other. 0.25 scores 77/79. See `ANATOMY.md`.
- **JK25_1 cannot be rescued by fixing the family.** It names no renal in any of 22
  recorded runs, and the infrarenal presets are gated on one. It needs a hand family call
  *and* hand-named renals.
- The plane-reach ladder must not go below 0.4×; station nudging rescues nothing under
  `validate_cut`; the BOUNDS proxy has been chased three times and left end-to-end flat.
- CAMP: the δK estimator is **locked**. Only partition rules may change,
  `min_patch_size` stays 18, and patch count stays geometry-driven. See
  `camp/CHANGES.md`, including its table of approaches that were tried and failed.

## The cut_bench delta against main — RESOLVED, not a regression

`cut_bench run --limit 12` against `main` reports **ALL clean 296 → 276 (69% → 65%)**
and lists 20 cuts clean on main and not now. Every one is a **trunk** cut on a GT case.
Chased down 28 Aug with `tools/session_vs_raw.py`; all 20 are accounted for, and none
is a loss in what a person gets.

**12 of 20 — the correction moved up a layer.** This branch added `flat_trunk_offset`,
which main has no equivalent of: a session detects a trunk cut that comes out badly and
backs it off to where the vessel is flat, recording the move as
`adjustments.vessel_end`. All 12 come out VALID through a `CutSession`, and all 12 were
moved to get there (e.g. `KK20_2` root 65.1 → 48.8 mm). `cut_bench` probes raw
`cut_terminal` and never builds a session, so it cannot see any of this.

**8 of 20 — a deliberate convention change, working as intended.** These are refused by
the obliquity gate: *"the cut plane meets the vessel at 56°, it is slicing ALONG the
vessel rather than across it"*. The gate is **identical on both** revisions
(`MAX_OBLIQUITY_DEG = 50.0`, and `obliquity_deg` is byte-for-byte the same function).
What changed is the **cut mode on aortic trunk ends**:

| | main | this branch |
|---|---|---|
| `Aorta distal`, `Aorta root` | `normal` | **`z`** |

Centerline-normal is perpendicular by construction, so main's obliquity is always 0°
and it accepts anything. This branch cuts aortic ends **horizontally** on purpose (see
"keep the aortic ends horizontal"), which is the reporting convention — and at 80–90%
along a descending aorta that runs horizontal, a horizontal plane really is slicing
lengthwise. Refusing is correct, and the message names the remedy: *"Switch it to
Centerline-normal, which is perpendicular by construction."*

**So the benchmark delta is an artefact of what `cut_bench` measures**, not a defect.
It scores `dev`/`euler`/`foreign` on raw `cut_terminal` and never calls `validate_cut`,
so it sees neither the session's recovery nor the obliquity rule. Worth knowing when
reading any future `--baseline main` run.

## CAMP: which lineage is canonical — SETTLED 31 Aug 2026

**Ours is canonical.** Decided by CAMP's author. The copy in `camp/` — `campsrc.py`
with the connectivity split, the fan gate, the orphan flag and the seed parameter — is
the version to build on, cite and measure with. Where it disagrees with any other copy,
including `Y:\aorta\python_code\camp`, ours wins. A paper's code-availability statement should
point at `https://github.com/SurgBioMech/CAMP`.

Re-verified on the day the decision was made, because the claim rests on it:
`verify_upstream_equivalence.py` on KK2_5 gives **0.397519623517** from the July 2026
upstream and **0.397519623517** from ours with the additions switched off — identical.
With them on, 0.016950636521.

The lineages **had** diverged in both directions and that still matters for porting,
just no longer for choosing. The older share copy holds a parameter-sweep harness,
`mesh_volume_mm3` and an `n_jobs` runner that were never in the `UpdatedCamp` copy our
additions were made on, and `campsrc.patch` fails 3 of its 8 hunks against it. **The
direction is now one-way:** port those forward onto ours if they are wanted, never the
reverse. They touch `patch_pass` and `camp` only, which is not where the sweep harness
lives.

All three copies carry the decision and are in sync as of 31 Aug: `camp/` here,
`SurgBioMech/CAMP` on GitHub, and `camp_updated_2026-08-26/` on the share. Nothing in
`Y:\aorta\python_code\camp` has been modified at any point.

## branch_cut.boundary_loops made fast, output unchanged — 25 Sep 2026

A speed fix inside the freeze, because on the web it is what makes a large mesh
unusable. `boundary_loops` grouped the boundary vertices with one Python pass per
connected component of a graph whose nodes are every vertex, so it ran once per
**vertex**; it is called on the whole mesh for every trial cut, ~60 times in a first
cut. It now groups them in one sort. Measured with the web app's own sequence
(`CutSession`, then `nub_preset` + `recut`):

| mesh | faces | load | first cut |
|---|---|---|---|
| KY160 x4 | 159k | 4.8 → 4.4 s | 17.2 → 11.9 s |
| JK5 x16 | 589k | 40.3 → 30.1 s | 46.6 → 30.7 s |
| KY160 x16 | 637k | 22.1 → 18.0 s | 77.3 → 47.7 s |

**The measurement either side, all identical:**

- 290 real calls during load and first cut on KY160 and JK5 at 1x, 4x and 16x, each
  answered by both implementations and compared array for array: 0 differ. The cut
  mesh, presets and openings byte-identical on all six.
- `cut_bench run --baseline HEAD` on the two local sample cases (KY160_1, JK5_1): 96
  cuts, ALL clean 59 → 59, every record identical field for field before and after.
  **Not yet run on the two cohorts**, which needs the cohort paths for that request.
- The suite, one process per script, before and after: identical verdicts on every
  script (still red either side: `test_region_station` 10, `test_nub_session` 3). The
  four scripts that read the thoracic cohort were not run.
- `tests/test_boundary_loops.py` holds the old implementation verbatim as the
  reference, checks equality on synthetic edge cases, random holes and both sample
  aortas, and times the call against the edge pass it starts with: 56.6x on the old
  code (fails), 1.2x on the new.

What is left of the time is the trial-and-error itself: `suggest_pullback` probes up to
16 offsets with two full-mesh cuts each, and the nub search grows its cuts about 12
times, each rebuilding face adjacency from scratch. Cutting on the region around the
plane rather than the whole mesh is the next lever, and would not be byte-identical.
