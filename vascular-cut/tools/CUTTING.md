# Working on the cutter

Standing rules for unattended work on `branch_cut.py`. The point of writing them down
is that a change can be accepted or thrown away without anyone being asked.

## Where data may come from

Read nothing outside this project folder except the two cohort roots recorded in
`tools/cohorts.local.json`, and inside those, only meshes and centerlines. Imaging
directories (`dicom`, `dicom_sip`, `sip`, `mhd`) are off limits. If a case looks
incomplete, say what is missing and stop -- do not go looking for it somewhere else,
and do not widen the search to a parent or a sibling folder. The cohorts are
read-only: load from them, write results here.

## Getting oriented

The cutter work described here is **merged into `main`**: normalised ownership, the
connected-region rule, the mismatch preflight, and the offset suggestion. Start new
work on a branch off `main`, never on `main` itself.

```powershell
python tools/cut_bench.py cases --verify   # DO THIS FIRST -- see below
python tools/cut_bench.py score            # last run, no waiting
python tools/cut_bench.py run --limit 8    # a quick pass, ~4 min
python tools/cut_bench.py run              # all 60 cases, ~35 min
```

`cases --verify` loads every pair and checks the centerline lies inside the mesh. Run
it before trusting any number, because the way this goes wrong is silent: a trunk-only
surface paired with a branch centerline still loads, still cuts, and still scores. The
file name is only a hint -- `mesh/<id>.mat` is the trunk in GT and the branched surface
in JK -- so the geometric check is the one that settles it. All 60 pairs verify today.

`tools/bench_last.jsonl` holds the most recent run and is gitignored, so it survives
a branch switch but not a clean checkout. If it is missing, `run` rebuilds it.

The desktop suites are standalone scripts, **one per process** -- two `QtInteractor`
windows in one interpreter segfault VTK, which is why there is no pytest entry point:

```powershell
foreach ($t in Get-ChildItem tests\test_*.py) { .venv\Scripts\python $t.FullName }
```

Some take a model argument (`jk5` or `ky160`); see `tests/README.md`. They print
`ALL PASS` or a list of failures, except `test_empty_region`, which prints `NO BUG`.

## Never lose the working cutter

The shipped cutter is preserved three ways, because one of them is not enough:

| where | what |
|---|---|
| `git tag cutter-shipped-2026-08-04` | a name that does not move, even if `main` does |
| `origin/main` and the pushed tag | off this machine entirely |
| `_backup/cutter/branch_cut.SHIPPED.py` | a plain file, recoverable with no git at all |

`branch_cut.py` at that tag is blob `8c25713`. `_backup/cutter/REVISION.txt` records
the commit. To get back to it:

```powershell
git checkout cutter-shipped-2026-08-04 -- branch_cut.py branch_tree.py
```

**Hard rules.** Never force-push. Never `reset --hard` on `main`. Never delete a
branch or move that tag. Never commit to `main` directly -- merge into it, and only
when the accept rules below are met. `_backup/` lives inside the project folder and
shares its fate, so the pushed copy is the real backup.

## The objective

Maximise **CLEAN** cuts on the cohort. A cut is clean when all three hold:

| | |
|---|---|
| **flat** | the opening's rim sits within 0.5 mm of the plane that made it |
| **euler** | the result has exactly one more opening than the input |
| **foreign** | at most 1 mm² of surface removed from a vessel this cut was not cutting |

Nothing else counts as success. In particular a cut that produces a flat rim by
eating part of the parent is a failure, and no topology check can see it: it still
removes one connected disk, still leaves one boundary loop, and keeps the euler
number exactly right. That is how 727,389 mm² of collateral damage went unnoticed.

Refusing is not failure. A refusal is a cut declining to produce something it cannot
make flat, which was chosen deliberately over reaching onto a neighbour to close a
rim and reporting success. But refusals in the **usable range** (≤80% of vessel
length) are a cost to be driven down, because that is where a cut has to work.

## Running it

```powershell
python tools/cut_bench.py run                 # every case, ~35 min
python tools/cut_bench.py run --limit 8       # a quick pass while iterating
python tools/cut_bench.py score               # re-read the last run
```

Cohort roots live in `tools/cohorts.local.json` (gitignored). The benchmark judges
**output only** — never which faces an implementation picked — so it stays honest
across a rewrite.

## Accept or reject

`main` moves as work lands, and the default baseline IS `main`, so a change that has
already been committed there gets compared against itself -- silently, and the output
looks like a normal run. That has happened once already. Tags exist so a comparison
can be pinned:

```powershell
python tools/cut_bench.py run --baseline cutter-shipped-2026-08-04  # vs the old cutter
python tools/cut_bench.py run --baseline cutter-merged-2026-08-05   # vs the first merge
python tools/cut_bench.py run --baseline main                       # vs what is shipped now
```

Judge a NEW change against `main` -- from a branch, with the change NOT yet on `main`.
If a comparison comes back with the two columns identical, suspect this before
believing the change did nothing.

Accept a change when, against the revision it is meant to improve on:

1. **No cut goes from clean to not-clean** without being understood and written down.
   The scorecard prints these individually for exactly this reason.
2. CLEAN does not fall in any of the four buckets (branch/trunk × usable/junction).
3. Foreign surface removed does not rise.
4. Refusals in the usable range do not rise without a matching gain in CLEAN.
5. The 23 desktop suites still pass, one per process.

Reject otherwise, and record what was measured. A rule that fixes forty cases and
breaks five is not obviously good; look at the five.

**`--baseline` defaults to `main`.** On a branch that is not the comparison you want:
the "base" column becomes main, so the scorecard attributes the whole branch's
accumulated work to whatever you are testing right now. That misread a good change as
a flat rejection -- 67 broken thoracic cuts, none of them caused by the change under
test. Pass `--baseline HEAD` to ask "what did the edit in my working tree do".

Two tells that you are comparing against the wrong thing, both of which were visible
and both of which were ignored:

  * two runs with a real code change between them come back byte-identical;
  * a named regression refuses with the SAME deviation to two decimals whether the
    change is in or out. Probe one directly before believing the list.

**Rule 3 needs a decomposition, not a total.** A refused cut removes nothing, so it
contributes 0 mm2 of foreign surface. A change that turns refusals into cuts therefore
raises the total while taking nothing from anybody, and rule 3 fires on it. Split the
runs into cuts that succeeded BOTH ways and cuts that changed existence before
concluding anything:

| | cuts | foreign mm2 |
|---|---|---|
| succeeded both ways | 2,906 | 12,044 -> 12,521 (+476) |
| exist only with the change | 41 | +2,986 |
| exist only without it | 7 | -108 |

86% of that "+3,355 mm2 regression" was coverage. The rule still stands -- the +476 is
a real rise and was accepted deliberately, not waved through -- but judge it on the
like-for-like figure.

## What is already known

Do not re-derive these. Each cost a day.

- **Ownership by nearest centerline is wrong beside a fat vessel.** A trunk's axis
  sits one radius deep in the lumen, so trunk WALL can be nearer a thin branch's
  centerline than its own. 304 mm² of aorta claimed by a 3.7 mm branch on KY160_1,
  deleted as part of it, opening 17 mm off plane on a cut 32 mm clear of any
  junction. `vessel_coords` now normalises by each vessel's own radius.
- **A tear sits at the band's edge but is not caused by the band's width.** Off-plane
  rim points sit at a median 0.94 of `band_w` from the cut point, which reads like the
  bound binding -- and it is not. Growing the width on a ladder to 4x, keeping the
  grown width only where it actually closed the tear, fixed **one** cut in 432 and
  raised foreign surface 10%. The walk stops there for OWNERSHIP, not distance:
  `mine & near_cut` is false past the vessel's own wall however far `near_cut` reaches.
  The correlation with `band_w` is real and the causation is not. Do not re-tune the
  width; it has been measured and it is not the lever.

  A tear is detectable without doing the slice, which is worth keeping: a face that is
  DELETED while touching a face that STAYS leaves a raw edge, whereas a CLIPPED face
  leaves its edge on the plane. So `delete` adjacent to `untouched` is exactly the
  question "is this opening torn".
- **Cuts leak onto the sibling.** Branch → the strip of aortic wall between two ostia
  → branch next door, all below the plane and inside the band. On JK5_1, branch 9's
  rim came back with 23 of 98 points on branch 10, and branch 10's returned 22 on
  branch 9.
- **The refusals cost no CLEAN cuts.** Of the 564 cuts the branch refuses and `main`
  produced, **zero** were clean on `main`. Every one was already torn or contaminated.
  So a refusal is not a lost good cut; it is an honest answer replacing a bad one.
  Worth reasserting whenever the refusal count looks alarming.
- **Why the band comes up empty**, over all 237 usable-range refusals:
  172 (73%) the ring at the plane is owned by the PARENT, near the ostium, so
  ownership excludes it; 36 (15%) the plane meets no surface within `band_w`;
  29 (12%) the ring is this vessel's and still unreachable by the walk. Bucket B is
  the one worth solving, and note that a high `frac` means close to the junction.
- **Admitting the parent's ring cannot work, and the reason is conserved.** Three
  variants, each measured against the branch rather than `main`, 12 cases:

  | variant | flat | euler | CLEAN | foreign |
  |---|---|---|---|---|
  | (branch as it stands) | 290 | 382 | 263 | 8,165 |
  | ball of 2 radii around the cut point | 318 | 388 | **263** | 17,926 |
  | same ball, straddling faces only | 300 | 380 | **263** | 9,756 |
  | no ownership on the band, full `band_w` | 331 | 364 | **263** | 17,600 |

  CLEAN does not move. Flatness rises and the euler count falls by about as much,
  every time: the ring does close, and clipping the neighbour's surface to close it
  opens a second boundary loop, so the cut trades one failure for another. This is not
  a tuning problem and more surface will not fix it. Admitting a neighbour's surface
  always adds a loop, because the neighbour is still attached on both sides of the cut.

  So a cut at an ostium has no clean answer available by widening. Either it moves to
  where the vessel still has a cross-section of its own, or it refuses. The refusal is
  the honest half of that; the missing half is telling the user where the cut WOULD
  work, which is what the validator already says in words ("Pull it back") and what
  the tool could compute.

  **The same trade governs the usable range**, which is the more important finding.
  357 cuts there produce a single-loop opening that is not flat. Relaxing ownership
  makes them flat and costs the same amount of euler; tightening keeps euler and
  leaves them torn.

  **And the reason is not the budget, it is the size of the bite.** Every one of the
  41 cuts that became flat under the ungated band failed on foreign surface, and the
  median amount taken from the neighbour was **94 mm²** -- smallest 4, largest 1,034.
  Closing those rims is not clipping a sliver off the parent, it is removing a piece
  of it. The 1 mm² bar is not too strict; it is measuring something real. Capping the
  holes left behind does not help either: it recovered 3 of 18 lost euler counts, and
  a nick is not what those failures are.

  So this direction is closed. Seven rules were tried -- a neck of 2 radii, the same
  restricted to straddling faces, an ungated band, band growth on a ladder to 4x, the
  same keeping the width only where it healed, nick capping alone, and capping
  combined with the ungated band. CLEAN came back as exactly 263 every single time.
  A cut at an ostium cannot be made clean by choosing which surface to take, because
  the surface it would have to take belongs to the parent and there is a lot of it.
- **Two cases are not enough.** Every threshold tried during the investigation won on
  KY160_1 and lost on JK5_1, or the reverse. The trunk-versus-branch split only became
  visible at ~60 cases. Do not tune on the samples.
- **Mesh and centerline must be the same anatomy.** GT keeps the branched mesh in
  `branch_model/` (25 of 27 cases) or `branch_mesh/` (2); `mesh/<id>.mat` is the
  trunk-only surface. Pairing wrongly puts every branch centerline outside the mesh
  and is silent. `preflight` now catches it.

## Where it stood when this was written

3,143 cuts over 60 cases, `main` versus the `cut-leak-fix` branch, with the
trunk/branch split:

| bucket | CLEAN main → branch | foreign mm² main → branch | refused |
|---|---|---|---|
| branch usable | 57% → **66%** | 206,540 → 11,288 | 197 |
| branch junction | 5% → **10%** | 311,408 → 9,019 | 369 |
| trunk usable | 58% → **61%** | 129,087 → 49,742 | 40 |
| trunk junction | 22% → **24%** | 83,795 → 25,733 | 33 |
| **all** | 41% → **48%** | 730,830 → **95,782** | 639 |

CLEAN is up in every bucket and foreign surface is down 87%. One cut went from clean
to not-clean: `JK/JK4_1 Branch 1` at 70%, rim 0.80 mm.

## The two region "failures" were the fix working

`tests/test_region_gui.py ky160` flagged two cuts that pass where it expected a
refusal. That suite asserts a refusal for "Ishimaru zone 3" and "Aortic arch" because
on the shipped cutter their rims tore 7.6 mm off plane, and its own comment says a
pass means either the geometry improved or a check stopped firing -- go and look.

Looked. Measured off the geometry rather than asked of the validator, every rim of
both regions is **0.000 mm** out of its own best-fit plane, over 109-146 rim points
each. The branch fixed them, and normalised ownership is why: that boundary sits
beside the arch, where a trunk's axis lies a radius deep in the lumen, so aortic wall
used to read as nearer a thin branch's centerline than its own. The expectation was
updated, with the measurement recorded beside it.

Two things worth keeping from this. **A benchmark is not the test suite** -- the
cohort scores 3,143 terminal cuts and never exercises region planning, so it could
not have seen this either way. And **a suite going red is not the same as a
regression**; this one was the tool getting better, and the only way to tell was to
measure the thing the assertion stood for.

## The ownership guard, and the two cuts it costs

Against `cutter-merged-2026-08-05`, over 3,215 cuts on 62 cases: CLEAN 1,563 to 1,599,
refusals 664 to 606, foreign surface **91,543 to 40,565 mm²**. Most of that is the
trunk buckets, where trunk-junction foreign falls 30,400 to 1,666. All 23 suites pass.

Two cuts go from clean to not-clean, both understood and both accepted:

| cut | before | after |
|---|---|---|
| `GT/KK15_1 Branch 3` @85% | flat, euler 1, no foreign | rim 0.62 mm, euler 0, 3.2 mm² of Line9 |
| `JK/JK6_2 Branch 8` @80% | rim 0.48 mm, no foreign | rim **0.11 mm**, 8.4 mm² of Line22/Line20 |

Same mechanism in both: the guard makes a face HARDER to take from the centerline
nearest it, so more faces stay where they started, and in these two the branch keeps a
sliver its neighbour used to be given. Note the second is *flatter* than before and
fails only on foreign surface. Two cuts against 36 gained and 51,000 mm² recovered.

## Still one radius per vessel, and that is the next thing

The guard arbitrates better; it still arbitrates using a single median radius per
segment, which is wrong exactly where cuts are hard:

- KY160's **left subclavian** is 6.5 mm at its flared origin and its median says 3.69,
  so 97% of its own wall there is flagged doubtful.
- KY160's **brachiocephalic trunk** measures 12.64 mm because radius is taken as the
  distance to the nearest wall, and near a confluence that distance crosses the arch.
  `length >= 4r` then hides the trim-a-group menu entry for it -- the arithmetic
  excludes it, not the anatomy.
- The aorta is three segments, each with its own radius, so a face near a boundary is
  judged against whichever claims it and the two judgements disagree.

The fix is a radius PROFILE along each centerline, estimated from the faces rather
than by ray distance from the centerline, and used everywhere the scalar is now: the
doubtful test, the reassignment, `band_w`, and the trim filter. The trim filter should
stop guessing entirely and probe offsets with the tear test, so the menu offers what
actually cuts.

## Open, in rough priority order

1. ~~Tell the user where the cut would work.~~ Done. A refusal now names the largest
   offset that has a clean cross-section of its own, checked with the tear test and
   held to the same 1 mm² foreign bar a cut is. 42 of 42 suggested offsets cut CLEAN.
2. **The 29 "own ring unreachable" refusals are not a bug.** Closed. The walk really
   is fenced off from a ring that belongs to it -- on `GT/KK5_1 Branch 2`, 110 of the
   114 faces between them are admitted and two are not, being ostium wall that reads
   as the parent's. A goal-directed crossing of up to three foreign faces, along one
   shortest path to that ring, does reach it: refusals fell 48 to 44 on twelve cases.
   Every one of the four came out TORN, because a four-face ring is a fragment of a
   cross-section and not a cross-section. Gated on the tear test the crossing never
   fires at all, which is the correct answer: **turning a refusal into a torn cut is
   worse than refusing**, since a torn opening measures as a real one.

   This is the ninth rule to leave CLEAN at 263 on those twelve cases, and it is the
   one that looked most like a genuine defect rather than anatomy. Treat the junction
   band as closed unless something changes upstream of the cutter -- a better
   centerline, or attribution that is not per-face nearest-neighbour.

3. **`JK/JK4_1 Branch 1` at 70%** — understood, and left alone deliberately. 30 of its
   32 rim points sit exactly on the plane and two stray 0.80 mm, on a 1.73 mm branch;
   `main` returns 29 points all exact. A torn rim strays everywhere, so this is two
   unclipped vertices, not a tear. One cut in 3,143 and not worth a rule.
2. **The 36 refusals where the plane meets no surface within `band_w`**, and the
   **29 where the ring is this vessel's own and the walk still cannot reach it**.
   Both are small, both are probably real bugs rather than anatomy.
3. **`JK/JK4_1 Branch 1` at 70%** — the single regression, r=1.7 mm, rim 0.80 mm.
4. **Whether to merge.** CLEAN up in all four buckets and foreign down 87% is a strong
   case; the cost is 639 refusals, none of which was a clean cut on `main`.
