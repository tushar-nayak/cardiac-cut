# The nub-in preset

Every aortic branch cut back to a short stub standing proud of the aortic wall, the
same rule applied to every branch of every case. Built for the branches paper, where
what is wanted is not a useful length of each vessel but a consistent, minimal
representation of where it leaves and how big it is.

    nub_cut.py            the rule and the cut
    tools/nub_cohort.py   run it over the thoracic cohort and score what came out
    nub_viewer.py         click through the results and record a verdict on each

`branch_cut.py` is **not touched by any of this** -- the preset is a caller. So the
cohort benchmark does not apply to it, and the 23 desktop suites cannot regress from
it. If a change here starts to need one in the cutter, that is a separate piece of
work with `tools/CUTTING.md`'s rules attached.

## The rule

    n = s_wall + k * (D_aorta / 2)

`n` is measured **along the branch centerline, starting at the junction node**, which
sits on the aortic centerline. The walk begins on the aortic axis and goes outward up
the branch; `n` is where the cut plane lands.

| | |
|---|---|
| `s_wall` | how far up the branch before it is clear of the aorta -- **measured** |
| `D_aorta` | the aorta's local diameter where that branch leaves it |
| `k` | how much stub past the wall, in aortic radii |

The cutter takes offsets measured back from a vessel's far end, so the preset hands it
`offset = length - n`. That is the same expression for an ordinary branch and for a
trunk cut before it splits, because `n` is measured from the aortic end and the cutter
measures from the other one.

### Why `s_wall` is measured and not assumed

The rule was first written `n = D/2 + k*(D/2)`: travel one aortic radius to reach the
wall, then `k` of another for the nub. That is exactly right when a branch leaves
perpendicular, and wrong by the obliquity when it does not -- and the arch vessels do
not. On KY160_1, walking each branch out from its junction node:

| vessel | `D_aorta/2` predicts | actually clear at |
|---|---|---|
| Branch 1 | 13.3 mm | 18.0 mm |
| Branch 3 | 13.3 mm | 10.0 mm |
| Branch 2 | 18.3 mm | 21.0 mm |
| Line1 | 22.0 mm | 21.0 mm |

Out by up to 4.7 mm. At `k = 0.1` the intended nub is 1.3-2.2 mm, shorter than that
error, so on two of the four the plane landed **inside the aortic lumen** and there was
no nub to have. Measuring `s_wall` makes the proud length mean what it says at every
`k`, and where a take-off really is perpendicular the two forms agree exactly.

Over the whole thoracic cohort `s_wall` runs 10.3 to 33.9 mm, median 18.1.

### How `s_wall` is measured

Walk the branch centerline out from the junction node in 0.25 mm steps. At each
station, the distance to the aortic centerline against the aorta's own measured radius
at the nearest station. Clear when the distance exceeds the radius and **stays** clear
for 2 mm -- a run test, because the aortic radius is a percentile over binned wall
faces and wobbles, so a lone station over the line is noise rather than an ostium.

Two other signals were tried and rejected, both on KY160_1:

| signal | why not |
|---|---|
| where the branch's own faces begin | Fires 3-5 mm early. A branch OWNS its flare, the funnel where its mouth meets the aorta, and that flare lies inside the aortic envelope. Right answer to "where does this vessel's wall begin", wrong one for "where am I clear of the aorta". |
| where the perpendicular section collapses | Confounded by obliquity. Near an ostium the plane slices the aorta at a shallow angle, so the section shrinks as the branch turns away rather than as it leaves. Put Branch 1's wall at 9 mm, inside a 13.3 mm-radius aorta. |

## What the cohort says about `k`

30 thoracic cases, 84 first-order vessels. Planning is cheap -- no cutting -- so the
runway question can be answered before anything is cut:

```powershell
.venv\Scripts\python tools\nub_cohort.py plan          # sweeps k, no cutting
.venv\Scripts\python tools\nub_cohort.py run --k 0.5   # cuts, writes meshes
.venv\Scripts\python tools\nub_cohort.py kfit          # what k the cohort wanted
.venv\Scripts\python nub_viewer.py --dir output\nubs_k0.5
```

**0.1 is far too short**, as suspected, and for a reason worth stating: it is smaller
than the error in the term it was added to. It is not a nub length, it is a rounding.

**A nub cannot be shorter than the vessel's funnel.** A branch leaves its aorta as a
flare and only becomes a tube some way out, and a plane through a funnel does not come
back flat. Where the requested length is too short the run grows the nub until the rim
lies on its plane, and that grown length is a measurement: the shortest stub this
vessel can have and still be cut cleanly. At `k = 0.5`, 20 of 82 vessels had to grow.
The `k` each vessel actually needed:

| percentile | k | proud |
|---|---|---|
| 50th | 0.50 | 9.0 mm |
| 75th | 0.50 | 11.3 mm |
| 90th | 0.74 | 13.6 mm |
| 95th | 0.85 | 15.9 mm |
| max | 1.66 | 26.5 mm |

So `k` around **0.75** is where the cohort stops being grown behind your back and the
stub lengths become as uniform as the rule promises. Ask for less and the run silently
gives you more on a quarter of the vessels, which is the one thing the preset exists
to avoid. That is why `NUB_K` defaults to 0.75 -- it is the cohort's answer, not a
preference, and it is one argument to change.

At 0.5, 20 of 82 vessels are grown. At 0.75, 8. Both cut 82 of 82 planned nubs clean,
so the choice is about uniformity, not about success.

## The two anatomical cases the rule does not cover on its own

**A trunk that splits before it leaves the aorta.** Six of the cohort's first-order
vessels -- four bovine arches among them -- are 7 to 25 mm long inside an aorta 14 to
28 mm in radius, so no station on them is outside it. `d - r` is still negative at
their last station: -2.6, -5.8, -2.6, -10.2, -8.3, -2.9 mm. There is no nub to cut.
The preset hands the job **down to their daughters**, which are 34 to 185 mm long and
clear easily, so the great vessels are still represented, one generation lower. The
trunk's own row is kept in the plan so the reason is visible rather than inferred.

**A trunk with some runway but not enough.** KY160's brachiocephalic clears the aorta
at 21.2 mm and splits at 29.3, so it can hold 6.6 mm of nub against the 10.1 mm that
`k = 0.5` asks for. It is **shortened rather than handed down**: a 6.6 mm stub of
brachiocephalic is still a brachiocephalic, and descending would drop the vessel from
the set and replace it with the two it feeds. Shortened nubs are flagged, and one
under 2 mm is flagged harder -- KY244's clears 1.0 mm before it splits, which is
placeable and still not worth keeping without someone looking at it.

## Results

30 thoracic cases, 95 vessels considered, 85 nubs planned.

| | k = 0.75 | k = 0.5 |
|---|---|---|
| **cut clean** | **85 (100% of planned)** | **85 (100% of planned)** |
| cut at the length asked for | 75 | 64 |
| grown to the shortest that cuts | 8 | 20 |
| shortened, vessel too short | 2 | 1 |
| proud length | median 13.1, range 7.7-26.5 mm | median 9.0, range 5.1-26.5 mm |

The 10 vessels not placed at either k are 6 trunks that split inside the aorta and
handed their nubs to their daughters -- not losses, the great vessels are still there
-- and **4 that cannot be cut cleanly at any offset**: KK40_3's Branch 5, KY238_1's
Branch 3 and Branch 2, and KY244_1's bovine brachiocephalic. The last has 2.4 mm of
vessel outside the aorta before it splits. The other three have room, so the obstacle
is the surface rather than the length; they want doing by hand.

16 of the 30 cases carry something worth a look -- a grown nub, a shortened one, an
uncuttable vessel, an aorta end that had to be cut perpendicular instead of horizontal,
or one that would not open at all. `nub_viewer.py --flagged` shows exactly those.

**Watch the vessel count, not the percentage.** When the 8 mm horizontal rule first
went in it read as a perfect run -- 100% of planned nubs cut clean -- while three
cases had been lost whole. The percentage was perfect *because* those cases had
stopped contributing to the denominator. Vessels 95 -> 83 and planned 85 -> 75 was
the tell, and a change that only touches the aorta's two ends cannot honestly move
either.

Clean means the same three things `tools/cut_bench.py` means: the rim lies on the
plane that made it (<= 0.5 mm), the result has exactly one more opening than it
started with, and the cut took no more than a square millimetre of any vessel it was
not cutting.

## The aorta's own two ends

They are **not** nubbed. The ascending and descending ends are cut **8 mm in from
each end on a horizontal plane** (`TRUNK_OFFSET_MM`, `TRUNK_MODE = "z"`), which is
how the Simpleware exports square these off and what keeps the two ends looking the
same from case to case. The mode is set rather than assumed: `classify_terminals`
already picks `z` for the trunk ends, but only when asked, and a case where it did not
would be cut perpendicular to a centerline that is 17-25 degrees off horizontal there
-- one case in thirty with a differently shaped end, which is exactly what nobody spots
until the figures are drawn.

**The horizontal plane cannot be unconditional.** Where the aorta itself is steeply
oblique at its end, a z plane stops going across the vessel and starts running along
it -- 52, 56 and 58 degrees against the 50 the cutter allows, on KK24_1, KK40_3 and
KP35_1. `preflight` refuses before anything is cut, so each of those took its **whole
case** down, branches and all, the moment the 8 mm z rule went in. Those three ends now
fall back to a centerline-perpendicular cut at the same 8 mm, and say so.

Preflight's complaint is now recorded rather than raised, for the same reason. It
catches configurations worth knowing about, but every individual cut is defended and
every opening is measured afterwards, so a complaint is something to show the reviewer
-- not a reason to return no geometry at all for a case whose branches are all
perfectly cuttable.

A refusal at an aortic end must not be fatal either, and used to be:
DJ18_2's distal aorta comes out 15 mm off plane at every offset from 2 to 10 mm, and
its root only survives at 2 mm, both properties of that case's mesh and nothing to do
with the preset. Re-raising cost the **whole case**, all three of its nubs with it. It
now tries to place the end somewhere workable, and failing that leaves it closed and
says so. Three nubs and an unopened aorta is a case for the reviewer; no case at all is
a case nobody sees.

## Four things that were silent, and are not any more

**The probe must be run against the UNCUT mesh.** `suggest_pullback` accepts a
candidate partly on `euler == 1`, which says the surface carries exactly one opening
-- true of a single cut on the closed mesh, false the moment a second exists. Run
against the running mesh it rejected every nub after the first, four of seven on the
first three cases, for the crime of having neighbours.

**And the assembled model must be checked again afterwards.** The probe judges one cut
at a time, and that is not the model that gets built: a neighbouring cut can already
have taken faces this one's band needed, so a nub that probes clean can still assemble
ragged. Two nubs in 80, both on trunks whose cut also removes a subtree. A second pass
that grows only the torn ones by 2 mm took the cohort from 80 clean to 82.

**And the growth must be searched over nub LENGTH, not over offset fractions.**
`suggest_pullback` is the right tool for "pull this back until it works" from a
hand-driven slider and the wrong one here. It tries 0.85, 0.7, 0.55... of the offset
and then bisects, so where it lands depends on where it started -- and the acceptance
test is not monotone in the offset, a cut can pass at 176 mm, fail at 179 and pass
again at 152. KY160's left subclavian grew **6.9 mm when asked for k = 0.5 and 43.0 mm
when asked for k = 0.75**: a longer nub from a longer request. `grow_to_fit` walks
outward in 2 mm steps and refines at 0.5 mm, so the answer is the shortest nub that
cuts cleanly whatever was asked for. The same subclavian now lands at 16.2 mm from
k = 0.75 and 16.0 mm from k = 0.5; KR73's brachiocephalic at 26.5 and 26.8. Two routes
to the same vessel-determined length is the property the preset needs, and it is what
makes the grown numbers worth reading as measurements.

Nub cuts sit at 85-95% of the way along each branch, which is the regime the cohort
benchmark scores worst -- 15% CLEAN on branch junctions against 73% in the usable
range. That 100% is not the cutter being better than the benchmark says; it is the
probe-and-retry refusing to accept anything else, at the cost of a longer nub where it
had to.

## Reviewing the output

```powershell
.venv\Scripts\python nub_viewer.py --dir output\nubs_k0.5     # click through
.venv\Scripts\python nub_viewer.py --dir output\nubs_k0.5 --flagged   # only the ones that moved
.venv\Scripts\python nub_viewer.py --dir output\nubs_k0.5 --render    # a PNG per case
```

One window, one case at a time. Rims are ringed and labelled with what was measured,
coloured by it and not by eye: **green** cut at the length asked for, **orange** grown
or shortened, **red** rim off its plane or barely a nub. `g` keeps, `b` flags, `f`
jumps to the next unjudged case, `r` resets the view; verdicts are written to
`review.json` on every keystroke so a review can be stopped and picked up.

**Every case arrives in the same place, at the same size, in whatever orientation the
last one was left in.** Only the camera's DIRECTION carries over between cases, never
its position: these geometries sit at their own scanner coordinates, tens of
centimetres and in one pair over 1.6 m apart in z, so carrying the whole camera over
aimed it at empty space and the next aorta came up blank. It was on screen the whole
time -- just not this screen.

## Arch variants, and what one stub stands for

A stub carries a first-order aortic vessel. On a normal arch three vessels leave the
aorta and the brachiocephalic feeds two, so **every stub stands for at most two great
vessels**. On a BOVINE arch two leave the aorta and the shared stem feeds three, so a
stub on the stem would stand for three and the left common carotid would vanish inside
it.

So a first-order vessel feeding `SHARED_ORIGIN_VESSELS` (3) or more hands the stubs to
its daughters. **Decided on what the vessel feeds, not on whether it happens to be
cuttable** -- which is what made this inconsistent before it was written down. Eight of
the thirty thoracic cases are bovine by the structural test (two first-order vessels,
one feeding three great vessels: KK24_1, KK5_1, KR27_2, KY160_1, KY215_1, KY238_1,
KY244_1, KY68_1) and they were being handled two different ways. Six had a stem too
short to cut and were handed down by the older "cannot be placed" fallback; KY160_1 and
KY215_1 had a stem long enough, kept it, and came out with **two** branch stubs where
every other case has three. Same anatomy, different geometry, decided by an accident of
length.

|  | normal (21) | bovine (8) |
|---|---|---|
| stub 1 | brachiocephalic, feeding RSA + RCCA | brachiocephalic proper, feeding RSA + RCCA |
| stub 2 | left common carotid | left common carotid |
| stub 3 | left subclavian | left subclavian |

The one difference that remains is that a bovine case keeps its shared stem, because
the stem exists. That is a fact about the anatomy, not about the method.

**Failure modes worth stating in a write-up.** The bovine test is structural and does
not depend on the anatomical naming, but it only recognises the `[3,1]` topology. Two
cohort cases match neither pattern: KK28_11 is `[2,2]`, two first-order vessels each
feeding two, and KP35_1 carries a fourth first-order branch and so ends with four stubs
rather than three. Neither is bovine by this test, and neither has been confirmed by
eye. Bovine stubs also come out slightly shorter -- median 11.7 mm against 13.5 --
because a daughter is measured from where IT clears the aorta, and it leaves the stem
already outside.

## Working through the cohort by hand

Cutting all 30 cleanly is not the same as cutting all 30 *right*: a nub that had to
grow is a nub nobody chose. The workbench walks the cohort in the ordinary desktop
tool, each case opening with its batch cut ALREADY MADE -- trims added, sliders where
the preset put them, aortic ends squared off.

```powershell
.venv\Scripts\python tools\nub_workbench.py --list     # what is done, and why each is flagged
.venv\Scripts\python tools\nub_workbench.py            # work through the flagged ones
.venv\Scripts\python tools\nub_workbench.py --all      # every case
.venv\Scripts\python tools\nub_workbench.py --case KY160_1
```

Adjust, Export, close the window, and the next case opens. Exports land in
`output/nubs_final/` as `<case>_nub.stl`; a case counts as done when that file exists,
so the pass survives being interrupted -- close without exporting and it comes round
again. Only the flagged cases are offered by default, since the rest were cut at the
length asked for with a flat rim.

Nothing about the window changes. `branch_gui.py` gains three optional arguments
(`--nubs`, `--out`, `--name`) that are inert when absent, and `nub_session.py` sets
the same things a person would set. Three properties of the window it has to respect,
each of which fails silently:

- **`_build_rows` writes the current row values back into `_offset_memory` before it
  rebuilds**, so seeding the memory and rebuilding hands the seeds straight back to
  whatever the sliders already showed. Set the rows, not the memory.
- **Adding a trim rebuilds every row**, so all trims go in before any slider is touched.
- **Offsets are matched by SEGMENT, not by label**: the batch calls an interior nub
  `<seg>.nub` and the window calls the same cut `<seg>.trim`.

### Check the cohort opens before sitting down to a pass

```powershell
.venv\Scripts\python tools
ub_precheck.py
```

Opens all 30 the way the workbench will -- a real window, the real files, the batch
cut restored through `nub_session` -- and reports which arrive exportable. The window
refuses to export while its status bar is red, so this is the only check that says
whether a case can actually be signed off. The batch scoring its own mesh cannot see
it. **26 of 30 open ready.** One case per process, because two QtInteractor windows in
one interpreter segfault VTK.

### Both routes must produce the same mesh

`tests/test_nub_session.py` cuts a case headlessly, restores that cut into the window,
and requires an **identical face set**. It caught two ways they differed, neither of
which shows up as an error:

- The window's sliders move in **tenths of a millimetre**, so an offset off that grid
  cannot be represented there. KY160's left subclavian at 178.76 mm became 178.80 the
  moment it was opened -- 0.04 mm, 66 faces. Offsets are now placed on that grid
  (`OFFSET_GRID_MM`), which costs nothing these meshes can resolve.
- The window **snaps near-coincident rim vertices** unless its "clean" box is cleared,
  and the batch did not. Same offsets, different mesh. The batch now follows the
  window.

And a third, found only by opening the cases rather than scoring them:

- **The ownership cache is keyed by segment names and lengths, not radii.** Ownership
  is decided by distance normalised by each segment's radius, so the split-line repair
  changes its answer -- but it changes no name and no length, so the cache did not
  notice. The batch repairs first and computes ownership after, and got DJ18_2's
  aortic root at 0.478 mm, inside the standard. The window cuts once as it LOADS,
  before the repair, caching ownership from the 6.28 mm radius; the repair then
  changed the radius and not the cache, and the same case came up at 0.7 mm and
  refused to export. Same case, same offsets, two verdicts. The repair now clears
  that cache.

A mesh that changes the moment it is opened is a mesh nobody can vouch for, and the
paper's set has to be one or the other, not both.

### Two cases are named as abdominal anatomy

KK24_1 and KK40_3 come back with **only** abdominal artery names -- supraceliac aorta,
left gastric, splenic, common hepatic, celiac trunk -- and no thoracic ones. The other
28 are named thoracic. The names are internally contradictory, so they are wrong
rather than merely surprising: `name_anatomy` calls KK40_3's `Line15` the celiac trunk
while that segment measures 13.4 mm in radius -- 27 mm across, aorta-sized, against
about 8 mm for a real celiac -- and calls a 7.4 mm-radius vessel the left gastric.

**No cut moves because of this.** Nothing in the preset reads a name: the aortic route,
which vessels are first-order, where each wall crossing is and where each plane goes
are all geometric, and `name_anatomy` runs afterwards and feeds back into nothing. But
these are the names the viewer and the workbench print, and they would be the names in
a figure, so they cannot be left as they are.

Both cases are already in the problem set -- both needed the obliquity fallback on an
aortic end, and KK40_3 is one of the four that will not export. Whatever makes the tool
read them as abdominal is worth understanding before trusting either. `name_anatomy`
lives in `branch_cut.py`, so fixing it is a change with the cohort benchmark attached.

### The window and the batch must read the same centerline

**The window smooths a jagged centerline and the batch did not.** It probes the line as
supplied and, if its turning exceeds the jagged threshold, TICKS ITS OWN SMOOTHING BOX
and reloads at whatever is in its settings. Six of the thirty are jagged -- KK20_2,
KR27_2, KY183_1, KY240_1, KY7_1, SA21_1, up to 63 deg/mm against 5-10 for the rest --
so on those six the two were cutting different centerlines. Same offsets, different
planes.

Three of them refused to export a cut the batch had called clean, reporting that "the
centerline for this opening runs OUTSIDE the mesh" -- a true statement about a curve
the batch had never seen. The other three opened READY while being different geometry
from what was scored, which is the worse half: a case that fails loudly gets fixed, a
case that quietly disagrees goes into the paper.

The batch now smooths where the window would, records what it used, and the session
sets the window explicitly from that record rather than letting it decide or letting a
per-user setting decide.

### One vessel could not be spanned by its own plane

KK40_3's blocking vessel measures 7.43 mm in radius against a perpendicular section
33-42 mm across, so the automatic band -- 2.5 radii, 18.6 mm -- could not close a ring
around it. The rim came out 1.0-2.9 mm off plane at EVERY offset, which reads as "this
vessel cannot be cut" when the truth is that the plane was too small. At 30 mm it cuts
at 0.000 mm. `BAND_TRIES` widens the plane before giving up, exactly as a person would
with the Plane size control, and records the width so the window is given the same one.

### Mesh topology is not always a plain closed surface

The acceptance test asked for `euler == 1` after a cut. That is the statement "one
opening was made" only for a closed surface of genus 0, and two cohort meshes are not
one: **KY238_1 is watertight with euler 0** -- a tunnel through it -- and **KK40_3 is
watertight with euler 4**, two separate closed pieces. Cutting either gives an euler
correct for the shape and wrong for the test, so every candidate was rejected.
KY238_1's two branches were being refused while producing a rim **0.000 mm** off plane.

The test now counts openings: one more than the mesh already had. True of any topology.
KK40_3's second piece turned out to be a 4-face tetrahedron 0.2 mm across floating
clear of the aorta, which `drop_debris` now removes before anything measures against
it.

### Horizontal ends, wherever a horizontal cut is possible

Obliquity is a property of the STATION, not of the vessel: KK40_3's aortic end meets a
horizontal plane at 49 degrees 2 mm in and at 17 degrees 40 mm in. Conceding a
perpendicular plane on the first number throws away a horizontal cut that exists a
little further along. The preset now searches the vessel for a station where a
horizontal cut is both within the obliquity limit and actually clean, and only falls
back to perpendicular when there is nowhere. Two ends were saved that way -- KK24_1's
moved to 27.1 mm and KK40_3's to 12.6 mm -- and no case now needs a perpendicular
aortic end.

### The anatomy is detected, not declared

A thoracic aorta has TWO trunk ends, a root and a distal end; an abdominal one has a
single inlet. The tool reads which from the relative size of the largest openings, and
on KK24_1 and KK40_3 it reads abdominal. Their arch vessels then came back named after
celiac anatomy -- supraceliac aorta, left gastric, splenic, common hepatic -- which is
the visible half of the mistake.

**The invisible half is geometric.** That count decides which terminals are trunks, so
on those two ONE AORTIC END WAS CLASSIFIED AS A BRANCH and given a nub instead of the
8 mm horizontal cut every other aortic end gets. Both cases cut clean, validated clean,
and would have gone into the set with an aortic end shaped like a branch stub. The wrong
vessel name was the only thing on screen that gave it away.

This is the thoracic cohort, so `TRUNK_ENDS = 2` is declared rather than inferred, and
the session sets the window's Anatomy dropdown to match. Checked against all thirty: it
changes nothing on the twenty-eight already read correctly, and both offending cases now
name their arch vessels correctly and square off both aortic ends.

### The same root cause, three times

`segment_radius` returns a badly wrong number on certain vessels, and everything scaled
from it inherits the error: band width (a cut that cannot span its vessel), ownership
(a cut scored as damaging a vessel that is really itself), and the "another vessel"
test. Every fix above is a workaround in this preset's loading path. Fixing it at
source is a `branch_cut.py` change with the cohort benchmark attached, and on this
evidence it is worth doing properly.

## Where data may come from

Same rule as everything else here. The thoracic cohort root in
`tools/cohorts.local.json`, meshes and centerlines only, nothing else and no imaging
directories. `nub_cohort.py` deliberately does not call `cut_bench.find_cases()`,
which walks every root -- the abdominal cohort has no arch vessels, is not what this
preset is for, and is never opened.
