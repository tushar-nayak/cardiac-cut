# Standing brief: work the cohorts, find what breaks, fix it

For autonomous or looped sessions. The instruction is not "do the next thing on a
list" -- it is "go and find where this tool fails on real geometry, and move a
number". Read this, pick your own target, and leave the repository better measured
than you found it.

## The two benchmarks are the objective

| | what it scores | where |
|---|---|---|
| `tools/anat_bench.py` | family call, aortic route, naming, region presets | `--indir <cohort>` |
| `tools/region_bench.py` | region presets the TOOL accepts, end to end | `[cohort]` |
| `tools/cut_bench.py` | the cutter: clean openings, foreign surface, refusals | see `tools/CUTTING.md` |

**`anat_bench --cut` is a proxy and it flatters.** It cuts only a region's two
boundary planes; a real region cut also cuts every branch inside it and validates the
result. On the abdominal cohort the proxy says 177 of 182 and the tool accepts 146.
Use it to find out WHY a region fails; use `region_bench.py` for whether it does.

Both append every run to a JSONL and print the delta against the previous run. A
change is worth keeping only if a benchmark moves, and only if nothing else moves the
wrong way. Run the relevant one **before** you change anything, so you have a
baseline in the same harness rather than a number from a different script.

```bash
python tools/anat_bench.py --indir <cohort> --labels <labels.csv>
python tools/anat_bench.py --compare tools/anat_bench.jsonl tools/anat_bench.jsonl
python tools/anat_bench.py --indir <cohort> --case JK4_1 --verbose   # one case
```

To baseline a change you have already made, stash just that file:

```bash
git stash push branch_cut.py -q && python tools/anat_bench.py --indir <cohort>
git stash pop -q               && python tools/anat_bench.py --indir <cohort>
```

## Rules of evidence

These are not style preferences. Every one of them is here because ignoring it
produced a wrong conclusion in this repository, and the wrong conclusion got written
down as fact before anyone caught it.

1. **Never generalise from a handful.** Five hand-picked scans produced a root-cause
   analysis for half the abdominal cohort that was simply wrong, in a document that
   already warned against doing exactly that. Measure over every case you have.
2. **A holdout with no distractors measures nothing.** Six mutually contradictory
   renal rules all scored 27/27 against scans that had one candidate each. Widening
   the window to manufacture competitors dropped five of them to 3-5 out of 9.
3. **Baseline in the same harness.** Comparing a new script's number against an old
   script's number is not a delta.
4. **Check what the harness is actually reading -- and PROCESSING.** The first
   `anat_bench` run scored 49 of 49 cases as thoracic-with-no-branches because it
   globbed the single-line centerline instead of the branched one. Later it read the
   right file and still measured something else: the tool SMOOTHS a jagged
   centerline before cutting and the benchmark did not, so every number it reported
   was computed on a line the tool never uses. Correcting it moved the abdominal
   family call from 49/49 to 48/49 and "bounds severing a bystander" from 0 to 4 --
   none of which was ever true.
5. **Report the negative results.** Bovine arch does not separate on delta K; folds
   do not predict delta K's instability; the route was not why the renals went
   unnamed. Those cost real time to establish and are worth as much as the fixes.
6. **If a rule cannot be validated, do not ship it.** Offer the person a way to say
   it instead -- that is what manual artery naming is for.
7. **Judge a fix by what the TOOL accepts, not by an easier proxy.** A boundary-plane
   fitter scored on opening count "rescued" three region cuts that the very next
   recut still refused, because validate_cut checks rim flatness and the probe did
   not. Every search must be judged by the same test the product applies.

## Where to look for work

In rough order of how much is known to be wrong:

- **`tools/ABDOMINAL.md`** -- the live work list, with per-scan counts. Items 4-6 are
  open: the proximal cut control, the renal-caught-in-the-cut seven, the centerline
  pre-flight lint.
- **`tools/abdominal_log.csv`** -- the lab's own per-scan complaints, one row per
  scan with the six complaint families as columns. This is ground truth for what a
  user hit. Cross-reference any benchmark result against it: a scan the log calls
  broken and the benchmark calls fine means one of the two is wrong, and finding out
  which is real work.
- **The flags column of `anat_bench`** -- short routes, close family calls, renals
  that were a choice, cases with too few presets.
- **Cases that RAISE.** An exception is a bug with no argument attached to it.
- **`tools/mesh_defects.py`** over any cohort -- folds, slivers, non-manifold edges.
  14 of 30 thoracic nubs carry folded geometry from the segmentation.

## What good work looks like here

- A commit whose message states the measurement, both directions, over the whole
  cohort. Not "improves naming" but "30 -> 44 renals named, 28 -> 33 with 4+ presets,
  0 worse, 31 suites pass".
- A test in `tests/` that would fail if the fix regressed, driving the real pipeline
  rather than a mock.
- The reasoning in a comment where the constant is, not in the commit message. The
  next person reads the code.
- A negative result written into the relevant `.md` when something does not work, so
  nobody spends the afternoon again.

## Hard constraints

- **Files outside the working folder are off limits** unless the user pointed at that
  exact path in that request. Cohort folders reachable today were granted explicitly;
  do not wander to siblings or parents, and never open `dicom`, `sip` or `mhd`.
- `tools/CUTTING.md` governs merges: never commit straight to `main`, never
  force-push, never move the shipped tags.
- Run the suites one per process (`tests/README.md`); two Qt windows in one
  interpreter segfault VTK.
- Do not edit source while a suite is running against it.
