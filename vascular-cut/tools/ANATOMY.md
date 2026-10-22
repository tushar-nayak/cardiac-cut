# Thoracic or abdominal: how the tool decides, and where it is wrong

Every anatomical name and every region preset hangs off one call: does this aorta
run **end to end** through the model, or does it enter once and **terminate by
dividing**? A thoracic aorta does the first — in at the root, out at the distal
descending end, two trunk ends. An abdominal aorta does the second — in at the
inlet, out through the iliac bifurcation, one trunk end.

Get it wrong and nothing downstream can be right. `region_presets` offers Ishimaru
zones to one family and the aortic neck, infrarenal aorta and iliacs to the other
(`branch_cut.py`, `region_presets`), and `name_anatomy` names a branch after celiac
anatomy or arch anatomy on the same evidence. A misread abdominal case shows up to
the user as *"I can't find the option for cutting the iliacs"* and a list of
openings called `Branch 4 … Branch 10`.

## What it used to do, and why that failed

The second-largest opening was called a trunk end if its calibre was comparable to
the largest — `TRUNK_RATIO = 0.6` on the ratio of mid-segment radii.

That is the wrong measurement, and it fails in **both** directions:

- A common iliac is legitimately comparable to the aorta feeding it. **JK1_1**
  measures 0.657 against a 0.60 cut-off and was read as thoracic — an abdominal
  aorta with celiac, SMA, both renals and four iliac tips, offered an *Aortic arch*
  region.
- An **aneurysm** makes a thoracic aorta's two ends wildly unequal. **KK40_3** is
  plainly an arch — ascending, four arch branches, descending — but its descending
  is dilated to r = 24.9 mm against a normal ascending at r = 10.8, a ratio of
  0.435, and it was read as abdominal.

Measured over 44 abdominal (JK) and 30 thoracic (GT) cases, the calibre rule placed
**58 of 74** correctly. **14 of the 44 abdominals were treated as thoracic.** That
is the fault the lab reported, and it was systematic, not a one-off.

## What it does now

Walk down the trunk from the fattest terminal, always continuing into whichever
side carries the most outflow — which is the trunk by construction, so the walk
runs the length of the aorta whichever end it started from. At every junction, ask
how evenly the outflow **divides**:

> the lesser side's summed terminal cross-section, over the greater side's

Flow divides with **area**, not with radius, which is why the sides are compared by
summed r². A true bifurcation sends a quarter or more of the area down its lesser
side. A branch leaving the aorta never does. The threshold is
`PAIR_AREA_FRACTION = 0.25`.

JK1_1 divides 0.97 at its aortic bifurcation (two common iliacs). KY160's arch
divides 0.14 at its most even junction (aorta against a great vessel).

| | median | worst | next worst |
|---|---|---|---|
| thoracic (n=30) | 0.108 | 0.309 | 0.214 |
| abdominal (n=44) | 0.887 | 0.269 | 0.526 |

Excluding the single crossover below, the classes separate cleanly at 0.214/0.269
and the cut-off sits between them. **73 of 74** land in the right family.

Starting the walk from *every* terminal was tried and is much worse — 59 of 74.
Begin at a renal and the first junction offers "up the aorta" against "down the
aorta", two large comparable sides, which reads as a bifurcation on a thoracic case.

## Radii on an open surface

`segment_radius` measured with `signed_distance`, which needs a **closed** surface
to tell inside from outside. On an open one the sign is meaningless and the radius
comes back negative, sorting that vessel to the bottom of every ranking.

Two ordinary things produce an open mesh:

- a segmentation that arrived with holes — 3 of the 44 abdominal cases here
  (JK5_1, JK6_1, JK6_2), with centerline stretches reading to −89 mm;
- **a mesh this tool cut and the user loaded back in to carry on**, which is open
  by construction.

Reloading JK1_1's own export turned its inlet from +7.98 mm to −14.61 mm, dropped
the aorta out of the trunk ranking, and took every anatomical name with it — the
reported *"I cut the renal, exported, uploaded it again, and it does not recognize
the vessels"*.

The distance is still right when the sign is not, so an open mesh is now measured
with unsigned distance. Closed meshes take the original path and read exactly as
before, bit for bit. This alone recovered JK6_1.

A reloaded cut mesh still warns: its centerlines now run past the new rims, so
radii read long there (JK1_1's inlet comes back as 14.95 mm against a true 7.98).
The **ranking** is restored, which is what the family call and the naming need, but
the numbers on those stubs are not trustworthy and the tool says so.

## Failure modes that remain

**KP35_1 (thoracic, read as abdominal, 0.309).** Its first-order vessel feeds three
great vessels, so the arch genuinely divides more evenly than an arch usually does.
It was misread by the old rule too. The window now says the call was close — *"and
it was a CLOSE call (0.31 against 0.25) — worth confirming"* — and the family can be
set by hand under Centerline.

**JK25_1 (abdominal, 0.205, MISREAD as thoracic).** One common iliac is aneurysmal
at r = 17.7 mm against an 11.0 mm inlet, so the bifurcation is lopsided and divides
far less evenly than a healthy one. Aneurysms push cases toward the threshold from
both sides; this is the same weakness the calibre rule had, reduced rather than
removed. **This is the 1 of 49 in the abdominal family baseline.**

> This entry used to read *"abdominal, 0.269, correct but nearest the line"*. That
> was measured on the UNSMOOTHED centerline, which is not what the tool cuts: a
> jagged line is reloaded smoothed first, and on the smoothed line the split is
> lower. The benchmark records the flip exactly — `anat_bench_abd.jsonl` runs 1–9
> give family 1, split 0.269; from run 10 (13 Aug 06:31), where `jagged` first
> appears, it is family 2, split 0.225; a later `segment_radius` change took it to
> 0.205 (`anat_bench.jsonl`, 15 Aug). Corrected 28 Aug 2026 after a user reported
> this as the one case of 49 they could not cut.

Do **not** move `PAIR_AREA_FRACTION` for it. Thoracic reaches 0.249 and abdominal
falls to 0.225 on smoothed lines, so JK25_1 (0.205) and KP35_1 (0.309) sit on the
wrong side of *every* cut-off relative to each other. 0.25 scores 77/79 and that is
the rule's ceiling.

**And fixing the family would not rescue the case**, which is the part worth knowing
before anyone tries. JK25_1 names **no renal in any run on record** — `n_renals` is 0
across all 22 benchmark runs, in both eras. The three iliac and infrarenal presets
are gated on a named renal, so even in the runs where the family was read correctly
the only region offered was *Whole abdominal aorta*. To cut this case a user must
declare the family by hand **and** name both renals by hand. The lab's own note on
the scan gives the reason: resolution too poor to trace the coeliac branches.

**A single traced iliac.** If a centerline stops at the aortic bifurcation with only
one iliac followed, there is no pair to find and the case will read thoracic. The
geometry genuinely does not contain the evidence; declare the family by hand.

**Anything with no junction to walk** — an unbranched tube — falls back to the old
calibre ratio, which is correct there: two ends and no branches is two trunk ends.

## Checking it

    python tests/test_family_gui.py          # the declaration, through the window
    python tests/test_anatomy_gui.py jk5     # naming on an abdominal case
    python tests/test_region_planes.py jk5   # a region's own boundary planes

The 74-case sweep is not committed — it reads two lab cohort roots named in
`tools/cohorts.local.json` and is reproduced by classifying every pair and counting
trunk ends against the folder each case came from.
