# The abdominal cohort: what breaks, and why

Read from the lab's own per-scan log (`overview.xlsx`, sheet `JK`, 54 scans) plus
reproduction against the geometries. This is the work list for making the tool
usable on abdominal aortas.

**9 of 49 processed scans exported with no complaint recorded.**

## What the log says

| complaint | scans |
|---|---|
| "does not detect branches anatomy (therefore isolate region functions not working)" | **24** |
| ...of which also note "left iliac without bifurcation" | 18 |
| a renal caught in the cut ("too big cutting plane as the renal goes straight down") | 7 |
| renals not detected — polar / accessory variants | 5 |
| proximal cut placed too low, or cannot be moved without losing the iliacs | 5 |
| unusable input (5 mm slices, no contrast) | 4 |

The first row is half the cohort and it is one bug.

## Root cause of "does not detect branches anatomy"

**A first version of this document blamed the aorta route, on the evidence of five
hand-picked scans. That was wrong.** Measured across all 47 processed scans, with
the route swapped for four alternative rules: *renal named but not on the route,
0 scans, under every rule.* The route is never why a renal fails to be placed. The
number of region presets was identical under all four rules. Recorded here because
the mistake is instructive -- it is precisely the "do not generalise from a handful"
failure this document already warned about, made by the person who wrote the warning.

The actual gate is in `region_presets()`. An abdominal aorta gets three of its five
presets -- Aortic neck, Infrarenal aorta, Infrarenal aorta with iliacs -- only if
`find("renal")` matches a landmark. So the question is entirely whether
`name_anatomy()` names a renal, and on **18 of 47 scans it does not**:

| why no renal is named | scans |
|---|---|
| several lateral pairs fit, so the rule refuses to choose | **14** |
| no lateral pair survives at all | 4 |

The renal rule pairs two unnamed branches leaving at the same level (within 1.5
aortic radii in z) on opposite sides, and names them **only if exactly one such pair
exists** -- because an accessory renal (~30% of people) or a lumbar pair otherwise
steals the name. That refusal is deliberate and correct: a confident wrong label is
worse than no label. JK4_1 has four competing pairs, JK24_2 has five.

The four with no pair at all: JK23_3 (one unnamed branch, so no pair is possible),
JK14_1 (the two candidates are 1.5 radii apart in z), JK25_1 and JK25_3 (seven
unnamed branches each, all rejected on side or level).

### Why this cannot be fixed by a better rule

There is no ground truth in the cohort for *which* pair is the renals. Any automatic
disambiguation written against these scans could not be validated, only asserted.
That is what makes manual naming the fix rather than a fallback -- and it is what the
lab asked for unprompted: "the problem of branch type would be solved by manually
choosing which branch is which".

### The route, separately

Still worth correcting, just not worth 24 scans. The structural rule takes the FIRST
`pre_bifurcation_candidates()` entry with two openings past it, from a list sorted by
fewest-openings then LONGEST segment -- so where the iliacs never bifurcate it lands
on whatever else splits in two. JK20_1 got a 47 mm ruler on a 206 mm aorta.

Fixed by trying both that rule and an outflow walk (`trunk_walk`, the same test that
tells abdominal from thoracic) and taking the **more distal** answer, since the aorta
ends at its bifurcation and there is nothing below it. Neither rule alone suffices:
the outflow walk taken alone stops at 72 mm on JK5, where the SMA carries 0.324 of
the remaining outflow and reads as a division.

Measured over all 47: **1 better (JK20_1, 47 mm -> 150 mm, 1 preset -> 2), 0 worse.**
JK30_4 stays at 50 mm; its CT stops above the bifurcation and there is no ruler to
recover. Real reasons the iliacs do not bifurcate in this cohort: a stent (JK6_2), a
femoro-femoral bypass (JK10_2), low contrast (JK20_1), the CT stopping short (JK30_4).

## The rest of the list

**Manual vessel naming.** Requested by the lab, and the only thing that closes the
tail: no geometric rule will reliably identify an accessory renal, and polar renals
alone account for 5 scans. Naming drives the region presets, so a hand-named renal
immediately restores "Aortic neck" and "Infrarenal aorta". For a public website
this matters more than any detection accuracy figure -- it turns a dead end into a
dropdown.

**Proximal cut control.** JK17_3 "too low proximal cut"; JK26_1/2/3 "cut very low
below the renals"; JK1_2 "when you want to move the proximal cut on the aorta you
are not able to include the iliacs anymore". The proximal bound needs to be settable
relative to the lowest renal without giving up the distal end.

**A renal caught in the cut.** JK2_1, JK3_1, JK3_3, JK12_2, JK14_2, JK17_2. Measured,
and it took three wrong measurements to find the right one:

1. `cut_bench`'s foreign surface does not reproduce it at all -- the reported six have
   LESS foreign area than the rest (median 7 vs 12 mm²). It scores per-branch cuts,
   and this is a REGION BOUNDARY plane.
2. Foreign area on the region cut is ~10,000 mm² on every case, reported or not,
   because these presets begin AT the lowest renal: removing renal surface above that
   boundary is the region working, not clipping it.
3. Renal area removed is ~1,000 mm² on every one of 62 region cuts, for the same
   reason. Neither number can see the defect.

What does see it: **two planes through a closed surface must leave exactly two
openings.** More means a plane reached sideways through a vessel it was not aimed at
-- a renal running parallel to the aorta rather than leaving laterally, which is
exactly the logged "too big cutting plane as the renal goes straight down".

Over 182 region presets on 49 scans: **167 cut cleanly, 6 sever a bystander**, and the
failure is enriched five-fold in the reported scans (2 of 10 against 2 of 52 in the
first 64-cut sample). `anat_bench --cut` now counts this.

**A smaller plane fixes it, and the window is narrow.** Sweeping the boundary plane's
reach on the two worst:

| case | automatic | severs a bystander | clean | refuses, rim off plane |
|---|---|---|---|---|
| JK17_2 | 35.3 mm | ≥ 21.2 | **17.7** | ≤ 14.1 |
| JK4_2 | 63.9 mm | ≥ 51.1 | **38.4 / 32.0** | ≤ 25.6 |

**The tool already refuses these rather than exporting them** -- `ok=False`, "3
openings, expected 2". So this is not a silent bad export. The gap is that a branch
cut offers pull-back when it refuses, and a region boundary offers nothing, even
though a smaller plane demonstrably works.

The obvious fix -- search the reach the way `suggest_pullback` searches the offset --
was attempted and withdrawn unshipped: `band_width(t)` alone returns 35.3 mm while
`cut_terminal` sizes its own band radius-aware and actually uses 12.5, so a helper
built on the former makes planes BIGGER, not smaller. Understand that interaction
before trying again.

JK2_2 is the mirror image (a renal NOT included that should have been), so the reach
has to move both ways.

**Typed input path.** The output folder is a text box; the input is browse-only.
54 scans is 54 dialogs.

**Pre-flight lint.** A dozen scans carry notes like "centrelines not ok?", "no
contrast", "too bad original resolution (5 mm)". `mesh_defects.py` covers the mesh
side; the centerline half is designed and unbuilt. Saying "this centerline stops
before the iliac bifurcation" up front beats a silent 47 mm route.

## Order

1. ~~the route fix~~ **done** -- 1 scan better, 0 worse. Not the 24 it was billed as.
2. ~~manual vessel naming~~ **done** -- the real fix, in both the desktop tool and the
   web app. Each opening's name box is now a dropdown of the artery names this
   anatomy can carry; picking one declares the vessel and the region presets follow.
   Declaring is not renaming: `rename()` sets a label, `declare_artery()` sets what
   `region_presets` matches on. A declaration displaces the automatic name of the
   same artery, so the aortic neck cannot land on the vessel being corrected.
3. ~~typed input path~~ **done** -- both boxes accept a pasted path, quotes and
   trailing spaces trimmed, and say so immediately if it is not there.
4. **proximal cut control** -- next. JK17_3, JK26_1/2/3, JK1_2. Worth re-testing
   against manual naming first: several of these are "cut very low below the renals",
   and if the low cut came from the wrong vessel being named the renal, declaring the
   right one already moves it.
5. re-test the renal-in-the-cut seven against the plane-size control
6. the lint

Also done alongside, though not from this list: the **nub-in preset** is now a button
in both front ends. It was validated over 82 vessels in 30 thoracic cases and was
reachable only from `tools/nub_cohort.py`.

Measure anything here by re-running all 47 scans and counting presets, not by
inspecting a few. Do not generalise from a handful -- four separate conclusions have
now been wrong for exactly that reason, including the original root cause above.

## The plain path, which is what a website visitor does

`tools/smoke_bench.py` scores load -> cut -> validate -> export with no region at
all. It is the only number that decides whether the tool is safe in front of
strangers, and nothing measured it until now.

**33 of 49 abdominal scans cut and exported. 0 raised an exception.** Nothing crashes,
which is the part that matters most for a hosted app; every refusal comes with a
reason. Exports round-trip: the file reloads with the topology the tool claimed.

### 12 of 49 cannot cut the aorta's own inlet, and no control fixes it

**Mostly solved -- 6 remain.** The cause is below, under "The inlet, solved": a KD-tree
clamping past the end of a centerline, and the clamped distance being read as a radius.
The six still refusing are the ones where `segment_radius` ALSO reads the model's end
cap instead of the lumen; that is a separate fault, measured and not yet fixed. The
survey below is left as written because the sweep it records is what ruled out every
control on the row, and that conclusion still holds.

The single concentrated defect. `cut_terminal` refuses with "cannot be cut here, the
opening comes out N mm out of plane" on the inlet -- 9 cases directly, 3 more as a
knock-on. Every control the tool offers was swept over the failing set:

| remedy | rescued |
|---|---|
| `suggest_pullback` (shallower offsets, the existing button) | 2 of 16 |
| deeper offsets, 8 to 40 mm | 1 of 16 |
| switching the plane mode, z <-> centerline-normal | 1 of 16 |
| plane reach, 0.4x to 3.0x | **0 of 16** |

So it is not the offset, not the plane orientation, and not the plane size. Obliquity
is not the cause either: most of these are already perpendicular at under 1 degree,
well inside the 50 degree limit. Summary geometry does not separate them from cases
that cut cleanly -- tip-to-wall distance, segment length, radius and first-branch
station all overlap the working set.

A promising reading -- that the cut produces no opening near its plane and the rim
check then blames a distant one -- was WRONG, and is recorded so nobody re-derives
it: `cut_terminal` raises, so the cut never completes, and the mesh being measured is
the one from before it. The openings 21-35 mm away were the previous cuts.

What has not been tried: the mesh itself at the inlet. `tools/mesh_defects.py` over
these 12, and a look at whether the refusal tracks a local defect rather than a
cutting decision.

### JK1_2's inlet: the root cause, and a fix that was worse than the fault

Reported from the tool as "it looks like it should work there". It should, and the
reason it does not is upstream of every control on the row -- offset either way,
plane mode and plane size were all swept and none of them helps.

The chain:

1. the inlet segment is 17.8 mm long, and the aorta there flares from 22.5 mm to
   31.6 mm equivalent diameter over 10 mm;
2. `segment_radius` takes the MEDIAN over the segment and returns 8.59 mm, a 17 mm
   diameter for a vessel that is 31 mm across;
3. `vessel_attribution` normalises each face by its owning vessel's radius, so the
   fatter next-segment-down wins the contest for the inlet's own wall;
4. **31% of the surface -- 240 mm2 -- in the very 3 mm the cut was asked to remove is
   owned by that neighbour**, so the cut leaves it, sitting above the plane;
5. the rim is then 25-35% ragged, with a deviation that grows exactly in step with the
   offset (1.24 mm at 1.5, 2.73 at 3.0) because those faces never move.

Median rim deviation is 0.000: most of the rim is perfectly on the plane, and several
contiguous arcs are stuck above it.

**The obvious fix is wrong.** Letting the cut also take faces beyond its own plane and
inside its band, whoever owns them, gains one case on the plain path and destroys the
cutter: foreign surface 14,275 -> 121,942 mm2, CLEAN 2162 -> 2092, and 71 cuts that
were clean no longer are. "Beyond the plane and inside the band" is not a tight enough
guard, because on a branch cut the band routinely spans the parent.

The fault is in the ATTRIBUTION and has to be fixed there -- most likely by sizing a
segment's radius locally rather than as a median, so a flaring inlet is not handed to
its neighbour. That is a change to the core of the cutter and needs `cut_bench` run
against it, not a benchmark of region presets.

### The inlet fault generalises, and two fixes for it do not

Measured across all 49 abdominal scans, the twelve that refuse their inlet are
exactly the ones where the segment holding it has a median radius far below its local
calibre:

| | local radius / median | surface near the tip owned elsewhere |
|---|---|---|
| refuses the inlet (12) | **1.88** | **27%** |
| cuts fine (37) | 1.10 | 0% |

One fault, twelve scans. `segment_radius` takes a median over a flaring vessel,
`vessel_attribution` normalises by it, and the inlet's own wall goes to the next
aortic segment -- which the cut then leaves sitting above its plane.

Two fixes were built and both failed `cut_bench`:

| attempt | foreign mm2 | CLEAN | cuts broken |
|---|---|---|---|
| any terminal, beyond the plane and inside the band | 14,275 -> 121,942 | 2162 -> 2092 | 71 |
| trunk ends only, aortic segments only | 14,275 -> 54,122 | 2162 -> 2099 | 64 |

Both gained exactly one case on the plain path. A half-space is the wrong guard on a
curved vessel: an aorta bends back over its own cut plane, so "beyond the plane and
inside the band" reaches around the arch and takes wall other cuts need.

**Connectivity was not needed.** The fault was never in the cutter. It is a
measurement error in `vessel_coords`, and the geometry was as simple as the lab said
it looked.

### The inlet, solved: a KD-tree clamps, and the clamped distance was read as a radius

Look at what is above the plane before theorising about the cutter. On JK1_2 at the
default 6 mm it is **220 faces, 496 mm2, in ONE connected piece**, and slicing it with
a plain plane gives a rim **0.000 mm out of plane** -- at every offset, in both plane
modes. There is nothing up there but the aorta's own end.

The refusal was self-inflicted:

1. the aorta is split at z=517.6, `Line09` below and `Line06` above;
2. every one of the 86 disputed faces sits ~16 mm ABOVE `Line09`'s last centerline
   point, and `_segment_kd` CLAMPS -- it answers with the distance to that last
   point, 18.74 mm, which is mostly a distance ALONG the vessel, not across it;
3. the wall contest compares that to `Line09`'s radius there, 15.67 mm, for a miss of
   3.06 -- beating the 5.11 mm miss against `Line06`, whose stub they actually are;
4. so `Line09` won a quarter of the surface the cut had to remove, the cutter will
   not remove another vessel's wall, and the tongue left behind held the rim 25-35%
   off plane at every setting.

Past its endpoint a segment has no wall to be near. `_wall_distance` measures the
component across the vessel and refuses a claim more than one centerline point-spacing
past the end. The tolerance is the data's own resolution, not a fitted number.

**Only at a JUNCTION.** An endpoint means two things: at a junction the vessel stops
and another carries on; at a FREE end it does not stop at all -- the centerline ran
out before the mesh did and the cap past it is that vessel's own to claim. Gating both
alike is byte-identical over 4,145 cuts here, because past a free end there is only
the cap and nearest-centerline has already awarded it; the distinction is kept as the
true statement, not as a measured gain.

It is not a JK1_2 quirk. **All 49 abdominal scans** have surface claimed past some
endpoint, a median of 5.9% and up to 23.8%; on JK5_1 alone ten segments won wall up to
13.3 mm past their own ends.

#### Two ways this investigation nearly went wrong, both worth remembering

**`cut_bench --baseline` defaults to `main`, not to HEAD.** Run on a branch, the
"base" column is main, so the scorecard reports the whole branch's accumulated work as
if the change under test had caused it. That read as a flat rejection -- 67 broken
thoracic cuts -- and was wrong. Two tells were visible and ignored: two runs with a
real code change between them came back byte-identical, and the named regressions
refused with the same deviation to 2 decimal places with the change on and off. Use
`--baseline HEAD`.

**A refused cut contributes 0 mm2 of foreign surface.** So turning refusals into cuts
raises the total while leaking nothing, and rule 3 fires on a change that took nothing
from anybody. Decompose it before believing it:

| | cuts | foreign mm2 |
|---|---|---|
| succeeded BOTH ways (like for like) | 2,906 | 12,044 -> 12,521 (**+476**, 0.16 per cut) |
| exist only with the fix | 41 | +2,986 |
| exist only without it | 7 | -108 |

### Fixed: `segment_radius` measured the model's end cap, not the lumen

The second fault behind the inlet. `segment_radius` takes
the median distance from the centerline to the surface over the middle 20-80% of a
segment, which its own docstring says is to dodge the cap at a terminal. On a SHORT end
segment the middle is still beside the cap, so the dodge fails and the number it
returns is the distance to the model's end, not the vessel's radius.

The cleanest possible demonstration, same centerline, same 25 stations, same segment,
the only difference being whether the cap is present:

| JK1_2 `Line06` (17.8 mm long) | radius |
|---|---|
| on the input mesh, cap present | **8.586 mm** |
| on the same mesh after the inlet is cut away | **14.673 mm** |

8.586 is the value everything downstream uses. Across the abdominal cohort **39 of 891
segments (4.4%)** are measured against a cap -- always exactly one per case, always the
short (13-25 mm) segment carrying the inlet, understated 1.24x where any wall sample
survives at all. On ten of them EVERY sample lands on the cap, so there is no wall
reading to fall back on.

The fix, cheap in the common case: classify each sample by whether the vector to its
closest surface point is more along the centerline than across it, and take the median
of the WALL samples. Only when too few survive -- about one segment per case -- fall
back to rays cast perpendicular to the tangent, which a cap is broadside to and cannot
answer. A NEGATIVE median is returned untouched, because that is the mesh/centerline
mismatch signal `mismatched()` and `_no_radius_reason()` are built on, and rays fired
from outside a surface would come back with a confident positive number.

It reaches much further than the inlet, because a radius feeds trunk-versus-branch
ranking, naming and `band_width` as well as attribution:

| | before | after |
|---|---|---|
| plain path | 40/49 | **43/49** |
| region presets offered | 182 | **195** |
| regions end to end | 156/182 | **170/195** |
| abdominals with renals named | 43 | **46** |
| cases with 4+ presets | 34 | **37** |
| openings named | 432 | **434** |
| route under half the model | 1 | **0** |
| family called correctly | 48/49 | 48/49 |
| cutter CLEAN | 2139 | 2139 |
| cuts clean -> not clean | -- | **0** |

Accept rule 3 breaks again and again it is coverage plus already-broken cuts, not
damage to anything the tool would accept. Split like for like by whether the cut was
ALREADY over the 1 mm2 foreign tolerance:

| | cuts | foreign delta |
|---|---|---|
| base already clean on foreign | 2,277 | **+17.6 mm2** (0.008 per cut) |
| base already leaking | 662 | +1,896.6 mm2 |

Two cuts carry 1,470 mm2 of that: JK25_3's inlet at 80% and 85% depth, where a
`band_width` of 2.5x a radius that just doubled reaches into everything around a short
segment. They were leaking 107 and 112 mm2 before, so no clean cut was lost -- but a
ceiling on the band for short segments is the obvious next thing to look at.

### Regions: 90% -> 98%, and what the last few are

A region preset refusing was the second-largest defect after the inlet. Grouped over
both cohorts, 30 refusals were two faults wearing different messages.

**A boundary undercutting a branch's ostium (fixed).** "N openings, expected N-1" --
one opening too many. The extra loop is not a torn rim and not a severed neighbour: on
JK17_2's infrarenal bound the aortic cross-section is 163 points, 15.46 mm across and
0.00 mm off its plane, and beside it sits a 19-point hole 1.79 mm across, 35 mm away.
That is a branch's ostium. The bound sits a millimetre or two from where the branch
leaves, the removed side takes wall from around its base, and the hole is what is left.
Every face taken is the cut's own -- zero foreign surface -- which is why nothing in
the cutter touched it. One to three millimetres of clearance fixes it.

**A boundary with no clean cross-section where it was put (mostly fixed).** Refuses
outright, or comes back not flat. Same remedy, but the clear side depends on which
side the branch is on, so the search has to run both ways: JK25_3's aortic neck needs
-6 mm and refuses outright at +1; JK6_1's visceral bound needs -8.

| | abdominal | thoracic |
|---|---|---|
| before | 190/210 (90%) | 106/116 (91%) |
| after | **210/215 (98%)** | **108/116 (93%)** |

Four explanations were built and measured away first. Do not rebuild them:

| tried | result |
|---|---|
| the boundary plane's reach | swept 55 mm down to 6 mm on six presets, no effect at all |
| ownership of the faces taken | zero foreign surface in every case |
| the subtree being removed in part | taking it whole changes nothing |
| the subtree wrongly containing an inside branch | none does |

**Still refusing, and why.** Abdominal 5: JK14_2 x2 and JK30_4 are the celiac hole
below; JK25_3's visceral bound and JK6_1's aortic neck have no flat offset anywhere in
-8..+12 mm. Thoracic 8: KK4_2's bounds refuse at every offset swept, KY238_1 is a
folded centerline, KR6_1's bound cuts flat at its own offset and the preset still
fails -- which means that one is not the bound at all and wants its own look.

### Open: a cut takes a neighbour's ostium, and the hole is an extra opening

JK14_2's inlet cut removes 109 faces, 170 mm2, of the CELIAC TRUNK. That is the "12
openings, expected 11" it refuses with on the plain path, and one opening too many in
two of its region presets. JK30_4 has the same shape.

The mechanism is `ambiguous` in `cut_terminal`: a face more than RADIAL_SUSPECT radii
from the axis of whatever owns it may be walked through by any cut, as the ostium fuzz
a rim needs to close across. At an ostium a branch FLARES, so its own wall stands at
1.8 of its median radii and qualifies -- all 109 faces did.

Two fixes were built and measured, and neither is the answer:

* **"does the face fit ITS OWNER's wall"** -- withdraws the fuzz a cut legitimately
  needs, breaking JK14_3 and costing JK1_2 its default offset. Reverted.
* **"does the face fit THIS vessel's wall"** -- no regressions anywhere, and no gain
  either, because a branch's ostium sits ON the parent's wall surface and fits both.
  Reverted; it cost a KD query per cut for nothing.

Distance cannot separate them, and attribution already knows the answer. The thing
left to try is limiting how FAR from its own plane a cut may walk on ambiguous faces:
closing a rim needs fuzz within a few millimetres of the cut, and these holes are 15
to 20 mm away, well inside a band that is 2.5 radii and so 37 mm on an aorta.

### Folded centerlines

`centerline_defects()` reports a centerline that doubles back on itself: DJ18_2's
descending aorta reverses four times and touches itself at fourteen places, which is
why it draws as a thick tangle, and it is the one thoracic aortic end that cannot take
a horizontal cut plane. KK24_1 has twelve reversals. 2 of 79 cases, both on pieces the
loader had to split, no abdominal case at all. Reported, never repaired.
