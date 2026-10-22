# Isolating the Ishimaru zones

How the tool cuts an aorta into zones 0–4, what was tried and measured to death first,
and what is still open. Numbers are on the 30-case thoracic cohort unless stated.

## What a zone has to be

Five properties, all of them checkable, none of them sufficient alone:

1. **They tile.** Zone 1 ends where zone 2 begins. No wall claimed twice, none left over.
2. **Each is one piece**, and **no zone sits as an island inside another**.
3. **Each branch belongs to exactly one zone** — the zone its ostium is in.
4. **Boundaries are as planar as the geometry allows.**
5. **Zone 3 and zone 4 own no branch**, because their bounds are distances, not vessels.

Every one of these passed on versions that were then rejected on sight, which is why
`tools/zone_map_page.py` exists: the last judge is a person looking at 30 cases.

## The coordinate: harmonic, not projected

The obvious way to place a boundary is to project each vertex onto the aortic
centerline and cut at an arc length. It does not work on an arch, and the reason is
exact rather than incidental. Nearest-point projection onto a curve is single-valued
only within `1/κ` of it — the medial-axis condition. Measured on these cases, `κ·d` at
the arch has a median of 1.41 (KK20_2) and a peak of 13.6 (KR73_1): the station there is
genuinely multivalued, and no filtering fixes a question with two answers. On KY160_1
that put Ishimaru zone 3 at 129.7–158.9 mm for a span asked as 131.3–151.3.

So the wall carries a **harmonic coordinate** instead: solve Laplace on the vessel wall
with the inlet ring at 0 and the outlet ring at 1. It is smooth, monotone, its level
sets are closed rings, and it has no reach condition. The branches are stripped and
their mouths capped with triangle fans before the solve (`harmonic_coordinate_capped`),
so the holes do not bend the level sets.

Boundaries are then pulled onto their own best-fit plane by `flatten_boundary`, fitting
the ring itself rather than a plane through the landmark — a landmark plane cuts the
wall in two or three places wherever `κ·d > 1`. That takes arch boundaries from 7.8 mm
off plane to 1.6, and the cohort sits at a median of 2.46 mm over 107 boundaries.

## Branches: the ostium is a saddle

A branch belongs to the zone of its **ostium**, half-open at the distal end
(`branch_inside`), so zone 0 owns the brachiocephalic, zone 1 the left carotid, zone 2
the left subclavian — which is what the zone names mean.

Applying that needs to know which vertices are the branch, and that is the hard part. A
great vessel often **curls back over the arch and hovers above the next zone**. The
part of it the branch mask missed then took the wall's labelling, and the vessel came
out one colour for most of its length with a streak of the next zone down one flank.

Three families of discriminator were tried. **The first two are dead and the numbers are
here so they stay dead:**

| family | why it fails | measured |
|---|---|---|
| distance from the aortic centerline | the populations overlap | curled stub `d/r` 1.05–1.45, p10 0.91; wall 1.00, p90 1.06 |
| frontier width of a flood from the cut loop | a pure trade-off, no sweet spot | pullback 0→5 mm: 2934→269 wall vertices dragged, 4→58 stubs still split. Jump 1.15→2.20: 418→2934 dragged, 56→4 still split |

Both fail for the same reason: **once a stub curls back, nothing about where it is or
how wide it is says it is still a branch.**

What does is the shape of the surface at its base. A vessel wall is convex; the join
where a branch leaves the aorta is the union of two convex tubes, and such a union turns
inward all the way round. **The mouth is a closed ring where the surface is at its least
convex, and it is there whether the vessel stands up off the arch or lies back along
it.** That is local to the surface — no centerline, no radius, no distance.

*Least convex, not concave.* On a smoothed segmentation the fillet rarely crosses zero:
KR47_1's wall sits at −0.050 mm⁻¹ with a p90 of −0.013, and the ring its brachiocephalic
must be crossed at is about −0.045. Only the relative maximum matters, which is what the
bottleneck flood reads, so nothing depends on the crease being concave in absolute
terms.

`region_station.stub_regions` uses it:

- **Concavity** (`surface_concavity`) is the discrete mean-curvature normal from the
  cotangent Laplacian, read against the outward normal and negated. Checked against a
  sphere of radius 10, where it returns −0.1995 against a true −2H of −0.2. Smoothed
  twice: a single stray concave vertex inside a stub is enough to stop a flood, and
  smoothing *increases* coverage rather than blurring the crease away — KR6_1's second
  stub goes from 50 vertices at smooth 0 to 235 at smooth 2.
- **A bottleneck flood** from each stub's own cut loop orders the surface by the worst
  crease that must be climbed to reach it. Ordinary shortest path would not do: it adds
  cost up, so it grows with length, which is distance again under another name.
- **Where to stop is a statement about shape, not a threshold on concavity** — that
  number varies with calibre and mesh resolution. While the region is the stub its
  boundary is a ring the size of the stub; the instant it spills, the boundary is an
  aortic ring several times longer. So the level rises until the region's boundary first
  exceeds 1.5× the stub's own cut loop. Both lengths belong to the stub, so the rule
  carries no length scale of its own.
- **First failure, not last success.** Raise the level far enough and the region is the
  whole mesh, whose boundary is nothing at all, so the bound comes good again on the far
  side of the leak.
- **Levels are spaced by rank, not by value.** Concavity is nearly flat through the tube
  and climbs steeply at the mouth, so an evenly-spaced ladder spends itself on the flat
  part: on KR47_1's brachiocephalic, 25 of 40 rungs left the region at 33 vertices and
  the last two took it from 154 to 2704. The mouth was between them and could not be
  seen.
- **And the step must add vessel.** The length bound asks whether the rim is short; it
  cannot ask whether the rim is in the right place. On a vessel with a wide opening —
  KR47_1's brachiocephalic cut loop is 80 mm round — half again is generous enough to let
  a splash of aortic wall through underneath it, which showed up as an arrow-shaped
  tongue on the ascending aorta. So each step is also judged by the median distance of
  its **new** vertices from the aortic centerline, in local aortic radii, and growth
  stops below 1.10.

  **This is the radius test, which is dead per-vertex** — a curled stub sits at 1.05–1.45
  and the wall at 1.00 with p90 1.06, so no threshold splits them. A median over a step
  is not the same statistic: each step adds dozens of vertices at once, and on KR47_1
  every legitimate step adds at a median of 1.18–1.51 while the step laying down the
  splash adds at 1.02. It is a question about a set, not about a vertex. No single value
  of the length bound could do this — at 1.5 that case leaked, at 1.2 a different vessel
  on the same case dropped to nothing but its rim.

`own_stubs` then spreads each stub's ostium zone over the whole stub and repairs
connectivity **with the stubs frozen** — without that guard the repair hands the stub
straight back to whichever zone surrounds its mouth.

## What it costs, and how that is known

Applied after the boundaries are flattened, so it can only repaint a stub:

| | before | after |
|---|---|---|
| stubs carrying two zones | 25 | **0** of 89 |
| zone islands inside another zone | — | **0** on 0 of 30 cases |
| zone pieces | 137 | 137 |
| vertices repainted | — | 1,327 over 24 cases |
| boundary flatness (107 boundaries) | 2.46 mm median | unchanged |
| overlap faces / unclaimed faces | 0 / 0 | unchanged |

Trunk damage is checked against a vertex being nearer the aortic centerline than 0.95×
the local aortic radius, the radius read off the vessel's own cross sections so branches
cannot inflate it. It reports 28 vertices, all on KK20_2.

**Read that number knowing it is no longer arms-length.** The step rule above draws on
the same `aortic_radius_profile`, so finder and audit now share an estimator. It was
independent when it mattered most — an earlier attempt was reported as changing no wall
and was not, and this test found 168 of 193 changed vertices on KR6_1 to be aortic — but
it is not independent now. The checks that remain independent of the finder are the
topological ones: bounded, one piece, no enclosed pockets, no zone islanded inside
another, zone-piece count unchanged. And a person looking at 30 cases, which is what
caught the arrow.

`tests/test_stub_regions.py --prove` reinstates three defects — stopping at the last
success, skipping the pocket fill, skipping the connectivity repair — and fails if the
suite does not catch each. The third needs an arch whose brachiocephalic mouth straddles
the 0/1 boundary; neither repo-local case has one, so KR52_1 is pulled from the cohort
when it is mounted and the run says so when it is not.

## Where zone 3 begins: a plane, not a level set

The clinical figure draws every zone boundary as a transverse cut, with the branch at
the distal end of its zone and the next zone starting immediately past it. The harmonic
coordinate cannot draw that. **A level set follows the wall, so on a bend it runs
oblique**, and no amount of sliding it along the vessel makes it transverse. Three
placements were tried and each was rejected on sight. What was wrong with all three was
not where they put the line but what shape the line was.

So the 2/3 boundary is a centerline-normal plane taken off the mouth itself
(`mouth_distal_plane`, `split_by_plane`). One rim vertex is projected onto the
centerline — a single point projects cleanly where a whole ring would not — and the
plane through it, normal to the tangent, separates zone 2 from zone 3. Zone 3's far side
is that plane carried 20 mm on, which is Ishimaru's own definition.

### Neither lip: the plane that clears the mouth

Both anchors below were a lip of the mouth plus a constant, and both were fitted to the
five hand-marked cases. Shown across thirty maps, each failed **in both directions at
once**: overshooting the branch where the mouth is short, cutting through it where the
mouth is long. The mouth's own span runs from 1 to 20 mm of centerline in this cohort,
so no constant can track it — the 8 mm below was never a discovery, it was the span of
the compact mouths it was fitted to.

The quantity that does track it was found by measuring eight cases judged by eye against
every candidate. It is the plane's **signed clearance** over the ostium:

| clearance | verdict |
|---|---|
| −17.6 mm (whole mouth distal) | "encroaching on the LSA" |
| −7.3, −5.2 | zone-3 faces showing on a branch |
| **+0.05** | **"very close, nice job"** |
| +9.2 | "way too far" |
| +13.7, +18.6 | "comes up short, should be closer to that last branch" |

One number, monotone, with the verdicts not overlapping. So `mouth_clearing_plane` asks
for clearance zero directly: sweep the centerline and keep the FIRST station whose
normal plane leaves every rim vertex behind it. No percentile, no offset, no fitted
constant, and **no harmonic coordinate** — which is how KK5_1 and KY183_1 came to run
the rule at all. `harmonic_coordinate_capped` returns None on those two (a singular
Laplacian from stray disconnected trunk pieces) and `plane_boundaries` returned the
labelling untouched, so they kept the old level set and no change to the rule could
reach them. One of them was among the eight complaints.

Because the plane clears the ostium by construction it cannot cut a branch, which also
removes the faces that rendered as zone 3 sitting on a vessel: 677 faces over 16 cases
down to 29 on one. Those were never repainted stub vertices — zone 3 was a single
connected piece throughout and no crease-region vertex ever changed zone — but a face
with one corner inside a mouth and one corner past a plane that had sliced between them
draws as branch in the wrong colour. Freezing the crease regions' 1-ring for the plane
cut, while leaving the repair pass on the narrow mask, removes the remainder; the two
masks must differ, because a repair pass given the wider freeze strands vertices it can
no longer absorb and leaves new islands.

**What it costs.** Against the marks it is worse: pooled mean 5.2 mm against 3.0, median
2.8 against 1.5. The three compact-mouthed marked cases do not move; the two elongated
ones do — DJ18_2 2.5 → 6.6, KK16_1 4.9 → 9.1. That is structural. On those two the marks
sit INSIDE the mouth's span, so the marked line cuts the subclavian: on KK16_1 the median
station of the mark ring equals the ostium rim's median station exactly. What was drawn
is not a plane. It hugs the ostium's distal margin on one flank, and any plane doing that
must slice the other. Satisfying the marks and the eight judgements together needs a
boundary that follows the mouth and closes transversely elsewhere — a change of shape,
not of anchor, and not built.

### Which lip, and why it is the proximal one — SUPERSEDED, see above

The obvious choice is the mouth's distal lip: the vessel's last hold on the aorta, with
the next zone beginning just past it. Scored against hand-placed marks it is exact where
the takeoff is square and late where it is not, in direct proportion:

| case | mouth spans | distal lip, against the marks |
|---|---|---|
| DJ8_1 | 7.0 mm | +0.5 mm |
| KK32_0 | 7.9 mm | +1.2 |
| KK24_1 | 8.8 mm | +1.0 |
| DJ18_2 | 15.0 mm | +5.5 |
| KK16_1 | 19.7 mm | +10.6 |

A shallow takeoff trails downstream along the wall and the distal lip follows it. That is
the whole of the "zone 3 starts too late, especially on the diseased arches" complaint —
and because the lag is a property of each ostium's obliquity rather than a fixed
distance, **no constant nudge could ever have fixed it**. Several were tried; all
regressed.

The proximal lip does not move with obliquity: it is where the vessel first parts from
the aorta. Measured from there the same marks sit at +6.5, +6.8, +7.8, +9.2 and +9.4 mm —
one number, sd 1.3 against the distal lip's 4.2. So the cut goes a mouth's width past
where the vessel leaves. Over 167 marks: mean 3.0 mm / trimmed 1.7 / median 1.5, against
5.1 / 3.5 / 2.8 for the distal lip, with no marked case worse and the two the plane
already fitted unchanged to the decimal.

### A region the size of its own cut loop is not a stub

The length bound weighs a region's rim against the loop the flood started from. On a
vessel that flares just below where the nub cut crosses it, the first rung already
carries a rim 1.68 times that loop, so the ladder breaks before it has reported anything
and the region comes out as the seed loop plus one vertex. Every guard in
`tests/test_stub_regions.py` passed it — bounded, whole, unfreckled, spilling no wall,
because it is nothing. But its rim is then the cut opening partway up the vessel instead
of the ostium, so the mouth every later rule reads is in the wrong place: on KK16_1,
11 mm downstream, carrying the zone boundary with it.

That failure announces itself, so it is retried rather than met by loosening the bound
for every stub. Over the cohort 31 of 90 stubs stop on the length test and **not one** of
them stops on a step whose new vertices look like aortic wall; the ratio test catches all
59 leaks unaided. The retry drops the length bound and leaves the ratio test to stop the
flood. It fires on four stubs in thirty cases — KK16_1, KK32_0, KP35_1, KR6_1 — and the
cohort separates cleanly: those four add exactly one vertex beyond their cut loop, and
the smallest genuine stub adds 17.

## Measured and rejected

Kept here so they are not re-tried:

- **Fitting the full ring including branch vertices** when flattening: flatness median
  2.35 → 3.20 mm, max 4.02 → 23.78. A mouth's vertices sit off in the branch's direction
  and drag the plane.
- **Majority-label smoothing** of the boundary: raggedness 1.63 → 1.59, flatness
  1.89 → 11.12 mm.
- **Cutting at one aortic radius** to find stubs: detection fell from 81 of 90 stubs to
  3, because the radius estimate was itself inflated by branch vertices.
- **Freezing the stubs while flattening**, now that a reliable mask exists: a wash.
  Raggedness median 2.326 → 2.294, flatness 2.74 → 2.82 mm, islands 0 either way.
- **Placing boundaries on the landmark's own cross-section** (`CROSS_SECTION_PLACEMENT`).
  This one is worth reading twice. Measured by projecting every finished boundary ring
  back onto the centerline and comparing with the station it names, it is unambiguously
  better: median displacement 0.15 mm → **0.00**, worst 15.0 → 10.8, boundary flatness
  2.52 → 2.44 mm. It also fixes KK40_3, whose 1/2 boundary otherwise sits 15 mm past the
  left carotid with a 48 mm spread so zone 1 never closes around the trunk.

  It was turned on, and rejected on sight: it grows zone 3 into a lobe over the arch
  apex and squeezes the bands either side of it, on every case. **Every number improved
  and the zone map got worse.** Landing a ring on its landmark's plane says nothing about
  whether the band it bounds is the right shape, and nothing being measured was asking
  that question. KK40_3's boundary is left open rather than trading every case's
  appearance for it.

## Measuring what ships

Three separate times, a checker in this directory was found to be measuring a program
the tool does not run: `zone_stub_bench` scored the branch-mask proxy after `zone_map`
had moved to the ostium, `label_case` placed boundaries on cross sections without
flattening onto them, and `zone_islands` reported an island on KY160_1 that does not
exist in the map. Each was a divergence between a reference implementation and the tool,
and each produced a confident number about nothing.

They all now route through the same functions in `region_station`, and where a flag
exists in two places (`CROSS_SECTION_PLACEMENT`) it is commented in both to say so. A
checker with its own copy of the rule can only ever agree with itself.

## Does this transfer to aortas nobody has looked at?

Audited deliberately rather than assumed, because "it looks right on thirty" is not
evidence about the thirty-first. The maps in this cohort are sound; the RULES behind
them are not yet portable, and the reasons are specific.

**The neck bound's denominator is not a perimeter.** `boundary_loops` returns each loop
as `used[labels[used] == c]` — its vertices in INDEX order, not in the order they sit
round the ring. `stub_regions` then walks that order with `np.roll` and calls the result
the cut loop's length. It is a zig-zag across the opening: measured against the true
perimeter over 27 loops it over-states it by **3.8x to 18.2x, median 7.6x**, and the
error grows with the loop's vertex count. So `NECK_ALPHA = 1.5` does not mean "half
again the stub's own opening" as its comment claims — it means about **11x** the
opening, varying five-fold between loops and drifting with mesh density. The one
genuinely scale-free quantity in the mouth finder is not scale-free.

This is why the sensitivity sweep finds no resolvable plateau for `NECK_ALPHA` at
+/-1.3%, and why `STEP_RATIO_MIN`'s cohort plateau is 0.4% wide with `PERP_COS_MAX`
sitting exactly on the upper edge of a 0.9%-wide one. The constants documented here as
dimensionless and safe are the fitted ones; the absolute-millimetre constants that
looked suspect (`PLANE_REACH_MM`, `BOUND_TOL_MM`, `MIN_STUB_GROWTH`) have plateaus
100-460% wide with the shipped value near the middle, and are fine.

**The mouth depends on where the branch was truncated.** Nub cut position is a
preprocessing artefact, but the flood seeds from it. The 2/3 boundary is exactly
invariant on 18 of 29 cases across +/-5 mm of cut, and on the other 11 it jumps 2-10 mm.
The response is a step function, not a drift: on a fine grid, moving the cut 0.15-0.88 mm
can jump the boundary 2.0-6.75 mm, and KY244_1 flip-flops 5 mm out and back over 2.9 mm
of cut movement. Because `mouth_clearing_plane` sits exactly at the rim's distal extreme,
ostium error passes into the boundary at gain 1, with no margin to absorb it.

**Mesh resolution moves the boundary up to 6 mm** on geometrically identical surfaces,
and `nub_preset` stalls outright at any scale other than 1.0 — so the cut stage is not
scale-portable even though the labelling stage very nearly is (median boundary movement
well under 1 mm at 0.7x and 1.4x, and rigid motion changes nothing).

**What does hold.** The topological half is solid on all 30 with no tuning: tiling, one
zone per branch, no islands, connectivity, ordering. The clearing plane is stated as a
definition rather than a fit and has no constant to be wrong. Rotation and translation
are exact.

**What a checker would say.** Run cohort-wide, the existing guards would flag about 4 of
30; a checker covering the geometric half would flag 19 of 30 rather than pass them
silently. The floor on silently-wrong is ~10%: KK24_1, KY160_1 and KY238_1 each carry a
stub whose mouth is 0-17 vertices and every check that runs today passes them.

Order of work if this is ever taken to a new cohort: fix the perimeter first, since
`NECK_ALPHA` cannot be re-derived until its denominator means something — and expect
re-tuning it to move every region, so the thirty maps have to be re-checked after.

## Open

- **Bovine arches.** On KK24_1, KK5_1 and KY238_1 two neighbouring vessels sit close
  enough that one crease region swallows both.
- **A crowded opening the retry cannot rescue.** JK5_1's sixth stub is a failed flood —
  it holds its cut loop and one vertex — and the retry declines it, because the very
  first rung reaches a neighbouring branch's opening. It is the same crowding as the
  bovine cases, on the abdominal side. The `GREW` guard in `tests/test_stub_regions.py`
  fails on it, and is meant to: the region is wrong, and it was wrong before the retry
  existed. Nothing downstream reads that mouth today, since the zone rules are thoracic.
- **Residue.** Scored on the older frontier flood's wider notion of a stub, 884 vertices
  still disagree with their stub's zone. That set runs past the mouth onto the aorta,
  where carrying the next zone is correct, so the number is an upper bound and not yet
  attributed.
- **KK28_11 has no zone 2 at all.** `aorta_landmark_stations` gives the left carotid and
  the left subclavian the same station (120.058), so Ishimaru zone 2 is zero-length, no
  stub can be inside it, and the plane rule has nothing to anchor on. The case renders
  0, 1, 3, 4 and is the one marked case the mark score cannot use.
- **SA21_1 keeps 29 zone-3 faces on a mouth** where every other case is now clean.
- **Carotid naming.** On DJ8_1, KK20_2 and KY7_1 "Left common carotid" shares the right
  subclavian's station, so the two are swapped and those cases lose zone 1.
- **`region_by_station`** still cuts on the projection path, so it inherits the reach
  problem the harmonic coordinate was built to escape. `tests/test_region_station.py` is
  red on 10 assertions because of it.
