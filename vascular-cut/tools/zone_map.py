"""Paint the whole aorta by which zone each piece belongs to, and show the seams.

    .venv\\Scripts\\python tools\\zone_map.py                    # thoracic
    .venv\\Scripts\\python tools\\zone_map.py Y:/aorta/abdominal/MK/JK --tag JK

One picture per case instead of one per zone, every zone a different colour. That is
the view that can be checked at a glance, because the two properties a zone set has to
have are both visible in it:

  NO OVERLAP   the Ishimaru zones tile the aorta -- zone 1 ends where zone 2 begins --
               so no piece of surface may be claimed twice. Anything claimed by more
               than one zone is painted CRIMSON, loudly, rather than being resolved
               silently by whichever zone was tested first.
  NO GAP       between the first zone's start and the last zone's end every piece
               should belong to something. What is left over stays wall-grey.

Only the Ishimaru zones are mapped. "Aortic arch" deliberately spans zones 0 to 3 and
would collide with all of them; it is a region to isolate, not a tile of a partition.
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
import branch_cut as BC                                      # noqa: E402
from anat_bench import find_cases                            # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
# Fraction of the case's median |grad u| below which wall is read as a
# branch stub rather than trunk.
BRANCH_GRAD_FRAC = 0.25
# How wide a band of the harmonic coordinate the boundary plane is
# allowed to redecide. Wide enough to cover the slant it is replacing,
# narrow enough that the plane never reaches the other limb of a fold.
# How far out from the cross-section contour, in local radii, the plane
# is allowed to decide. Beyond that the harmonic labelling stands.
COLLAR_RADII = 2.5
# How far a contour-derived boundary may sit from the station asked for
# before it is rejected as not a cross-section at all.
CONTOUR_TOL_MM = 8.0
# Place each boundary on its landmark's own CROSS-SECTION rather than by interpolating
# the station->u curve. ON: measured over all 107 boundaries by projecting each ring
# back onto the centerline and comparing with the station asked for, it takes the
# median displacement from 0.15 mm to 0.00, the worst from 15.0 mm to 10.8, and
# boundary flatness from 2.52 mm to 2.44. That worst case is KK40_3, where the curve
# put the 1/2 boundary 15 mm distal of the left carotid with a 48 mm spread, so zone 1
# never closed around the trunk.
#
# OFF, and the displacement numbers above are exactly why a metric is not a review.
# It was reviewed and rejected on sight: it visibly moves the zones, growing zone 3
# into a lobe over the arch apex and squeezing the bands either side of it. Boundary
# displacement improved and the zone map got worse, so the metric was measuring the
# wrong thing -- landing a ring on its landmark's plane says nothing about whether the
# band it bounds is the right shape.
#
# The cost of leaving it off is real and stays open: KK40_3's 1/2 boundary sits 15 mm
# distal of the left carotid with a 48 mm spread, so zone 1 does not close cleanly
# around that trunk. That is one case's boundary, against every case's appearance.
CROSS_SECTION_PLACEMENT = False
# A BOUNDARY IS A VESSEL, SO PUT IT ON THAT VESSEL'S MOUTH. By the half-open rule a zone
# owns the vessel at its distal end and the next begins immediately after it, so the
# boundary closing zone k belongs at the far lip of the mouth of the stub zone k owns.
# The crease regions locate that mouth on the wall directly, with neither the landmark
# station nor the station->u curve in between -- and those are what put it wrong.
#
# Measured as millimetres of wall between the left subclavian's mouth and the start of
# zone 3, the cohort was carrying a median 22 mm of zone 2 past the vessel, which is
# most of zone 3's 20 mm length. The landmark station is unreliable in BOTH directions:
# KK32_0's subclavian mouth ends 8.7 mm PROXIMAL of its station, KK20_2's 4.3 mm distal.
# So this is not a correction in one direction, it is dropping an estimate for a
# measurement.
#
# OFF: it over-corrects, and the reason is geometric. A zone boundary is a level set of
# the harmonic coordinate, so moving it by a fixed amount OF COORDINATE moves it by very
# different amounts of wall -- little on the outer curvature where the level sets are
# far apart, a great deal on the inner curvature where they crowd. On KK32_0 the mouth's
# lip is 8.7 mm proximal of the landmark station in arc length, and shifting the level
# set to it swept zone 3 across the whole underside of the arch. What was wanted was the
# band moved a little way over, uniformly; what a coordinate shift gives is a wedge.
RIM_PLACEMENT = False
# CHECK THE SUBCLAVIAN'S LANDMARK AGAINST ITS OWN MOUTH. A zone boundary is asked for at
# the station where the branch's centerline meets the aorta's, and on most arches that is
# the ostium -- the median lag over this cohort is 0.1 mm, so nothing moves. On a few it
# is not: the branch centerline runs alongside the aorta before diverging and the
# junction lands well past where the vessel leaves the wall, 14 mm on KK32_0. Zone 2 then
# runs on past the subclavian and zone 3 starts late.
#
# The mouth is not an estimate -- the crease regions find it on the wall -- so where the
# two disagree by more than a couple of millimetres the mouth wins and zones 2, 3 and 4
# all slide back together. Zone 3 keeps its 20 mm.
#
# OFF -- REJECTED ON SIGHT, and it made things worse rather than merely not better:
# KK16_1 ran on far enough to swallow the branch, KK32_0 broke outright, DJ18_2 got
# worse, and KK20_2 and KK24_1 were not improved at all. The lag is real and measured,
# so the landmark IS past the mouth on those arches -- but taking it out at the bounds
# does not put the zones where they belong, which means the lag is not the whole story
# and probably not the operative part of it.
#
# THREE ATTEMPTS AT THIS HAVE NOW FAILED, each for a different reason, each defensible
# on paper and wrong on the screen: boundary-on-the-mouth swept the inner curvature, a
# geodesic collar ate the wedge's narrow side, and this. Do not try a fourth without a
# clearer statement of where the boundary is supposed to land.
LANDMARK_LAG_FIX = False
LANDMARK_LAG_MIN_MM = 2.0
# GIVE EACH STUB TO ITS OWN OSTIUM. A great vessel that curls back over the arch is
# read as trunk by the branch mask along the flank that hovers over the next zone, so
# it takes the wall labelling there and comes out with a streak of the next colour.
# The stub is found from the crease at its mouth (region_station.stub_regions), which
# is local to the surface and so unaffected by where the stub subsequently leans, and
# every vertex of it is forced to the zone its ostium already chose. Applied AFTER the
# boundaries are flattened, so it can only repaint a stub -- measured on 30 cases it
# changes no vertex that an independent radius test calls aortic wall.
STUB_OWNERSHIP = True

# Five zones that tile, plus the two states that mean something is wrong. Picked to
# stay separable on a light or a dark page and to keep the two fault colours away
# from the five good ones.
# The palette lives in region_station so the browser, the desktop viewer and these
# renders cannot drift apart. OVERLAP is local: only this checker draws it.
from region_station import ZONE_COLOURS, ZONE_UNCLAIMED as UNCLAIMED  # noqa: E402
OVERLAP = "#e11d48"
VIEWS = [np.array([0.4, -1.0, 0.15]), np.array([-0.4, 1.0, 0.15])]
UP = np.array([0.0, 0.0, 1.0])
TILING = ("Ishimaru zone 0", "Ishimaru zone 1", "Ishimaru zone 2",
          "Ishimaru zone 3", "Ishimaru zone 4")


def zone_number(name):
    """The Ishimaru number in the preset's name.

    NOT ITS POSITION IN THE LIST. A bovine arch has no zone 1 -- the left carotid
    shares the brachiocephalic's origin, so that zone has zero length and is never
    offered -- and colouring by position then paints zone 2 in zone 1's colour and
    calls it zone 1. Measured: KK5_1 and KY244_1 are both bovine, both offer zones
    0, 2, 3, 4, and both were drawn as 0, 1, 2, 3.
    """
    m = re.search(r"(\d+)\s*$", name or "")
    return int(m.group(1)) if m else 0


def hex_rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--tag", default="GT")
    ap.add_argument("--cell", type=int, default=430)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    outdir = ROOT / "tools" / f"zone_map_{args.tag}"
    outdir.mkdir(exist_ok=True)
    cases = find_cases(args.cohort)
    if args.only:
        cases = [c for c in cases if c[0] in set(args.only.split(","))]
    if args.limit:
        cases = cases[:args.limit]
    print(f"{len(cases)} case(s) -> {outdir}")

    index, t0 = [], time.time()
    for i, (cid, cl, mp) in enumerate(cases, 1):
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            st = s.recut()
            mesh = s.mesh
            route = BC.aorta_route(s.tree)
            station, _own, on_route = RS.station_field(mesh, s.tree, route)
            zones = [p for p in s.presets() if p["name"] in TILING]
            zones.sort(key=lambda p: p["s0"])
            znum = [zone_number(p["name"]) for p in zones]

            # THE HARMONIC COORDINATE, where it can be built. Nearest-point projection
            # is multivalued wherever the vessel bends tighter than it is wide, which
            # on the arch is exactly where zones 2 and 3 meet, and no amount of
            # filtering the projection fixes a question that has no answer. Laplace on
            # the wall has no reach condition: smooth, monotone, level sets that close
            # around the vessel. Falls back to the station field if the two aortic ends
            # cannot be told apart, so a case is never silently mislabelled.
            # Branch mask first: the capped solve needs to know what to strip.
            is_branch = ~on_route
            # CAP ON THE CREASE REGIONS. The projection's branch call is wrong at the
            # arch, so the capped solve was cutting wall out rather than branches --
            # 2377 vertices on KK40_3 where the stubs are 511 -- and the field paid for
            # it in hot spots: |grad u| p99 0.01824 against 0.00516. That is where
            # level sets bunch and a zone boundary runs lengthwise down the vessel,
            # which on KY68_1 was a tongue of zone 2 down the whole ascending aorta.
            g_adj = RS.surface_graph(mesh)
            reg = RS.stub_regions(mesh, route, adj=g_adj)
            cap = (reg >= 0) if reg.max() >= 0 else is_branch
            u = RS.harmonic_coordinate_capped(mesh, route, cap)
            oz = (RS.stub_ostium_zone(s.tree, route, mesh, reg, zones)
                  if reg.max() >= 0 else {})
            lag = 0.0
            if LANDMARK_LAG_FIX and reg.max() >= 0 and u is not None:
                cal = RS.calibrate_u_to_station(mesh, route, u)
                mouth = RS.landmark_lag(mesh, u, cal, reg, oz, zone=2)
                lsa = [p["s1"] for p in zones if zone_number(p["name"]) == 2]
                if mouth is not None and lsa and lsa[0] is not None:
                    d = float(lsa[0]) - mouth
                    if d > LANDMARK_LAG_MIN_MM:
                        lag = d
                        zones = RS.shift_zone_bounds(zones, lag)
                        znum = [zone_number(p["name"]) for p in zones]
            if u is None:                       # capping failed -- plain solve
                u = RS.harmonic_coordinate(mesh, route)
            ss_cal, uu_cal = (RS.calibrate_u_to_station(mesh, route, u)
                              if u is not None else (None, None))
            use_u = u is not None and ss_cal is not None
            # The far lip of each stub's mouth, per zone. max() because a zone may own
            # more than one vessel -- a bovine trunk's daughters -- and the zone ends
            # after the last of them.
            by_zone = {}
            if use_u and RIM_PLACEMENT and reg.max() >= 0:
                rim = RS.stub_rim_coord(mesh, u, reg)
                for j, z in oz.items():
                    if j in rim:
                        by_zone.setdefault(z, []).append(rim[j])
            # WHICH VERTICES GET THE BRANCH RULE. Not station_field's on_route: that
            # is a nearest-centerline call and it fails where the station does, so
            # arch wall lying under a great vessel was read as branch and took that
            # vessel's zone -- whole non-planar patches, visibly worse than before.
            # A stub carries no flux, so |grad u| collapses on it and holds up on the
            # trunk; the threshold is relative to the case's own median.


            # Claim count per vertex, so an overlap is COUNTED rather than hidden by
            # whichever zone happened to be painted last.
            nv = len(mesh.vertices)
            claims = np.zeros(nv, int)
            which = -np.ones(nv, int)
            per = []
            for zi, p in enumerate(zones):
                if use_u:
                    # The first zone starts at the vessel's own start and the last
                    # ends at its end, so those bounds are opened out rather than
                    # taken from the calibration -- which lands a hair inside both
                    # rings and would leave the inlet and outlet unclaimed.
                    # THE BOUNDARY'S u COMES FROM THE LANDMARK'S OWN CROSS-SECTION.
                    # Interpolating the station->u curve is biased on a bend: measured
                    # against the true cross-sections, boundaries placed that way sat a
                    # median +3.0 mm distal of their landmark and up to 12.6 mm, so
                    # every zone grew into the next -- while the one boundary in the
                    # straight descending aorta landed at +0.1 mm.
                    def _ub(st, zone=None):
                        if zone is not None and zone in by_zone:
                            return float(max(by_zone[zone]))
                        # The cross-section is the better answer WHERE IT EXISTS. On
                        # the tightest part of an arch a normal plane slices along the
                        # vessel instead of across it, and the ring it cuts runs tens
                        # of millimetres down the wall -- on KK20_2 that returned
                        # u = 0.376 at a station the curve puts at 0.228, so u ran
                        # backwards between two landmarks and the zone between them
                        # collapsed. Trust the contour only when it agrees with the
                        # curve to within a few millimetres of station; past that the
                        # plane is not cutting a cross-section and the curve, biased
                        # but monotone, is the safer answer.
                        base = float(np.interp(st, ss_cal, uu_cal))
                        if not CROSS_SECTION_PLACEMENT:
                            return base
                        seg = s.tree.segments[route.locate(float(st))[0]]
                        v = RS.boundary_u_at(mesh, route, u, st, seg, prior=base)
                        if v is None:
                            return base
                        back = float(np.interp(v, uu_cal, ss_cal))
                        return v if abs(back - float(st)) <= CONTOUR_TOL_MM else base
                    lo = -1.0 if zi == 0 else _ub(p["s0"], znum[zi - 1])
                    hi = (2.0 if (p["s1"] is None or zi == len(zones) - 1)
                          else _ub(p["s1"], znum[zi]))
                    # EACH COORDINATE FOR WHAT IT IS GOOD AT. The harmonic field runs
                    # smoothly along the WALL, which is what made the boundaries land
                    # as clean rings -- but it also varies across a branch mouth, so a
                    # boundary falling there splits the vessel between two zones. A
                    # branch belongs to one zone entire, decided by its OSTIUM, which
                    # is what the projection coordinate was always right about. So the
                    # wall is cut by u and every branch is assigned whole.
                    m = (~is_branch) & (u >= lo) & (u < hi)
                    inb = RS.branch_inside(np.where(is_branch, station, -1e9),
                                           p["s0"], p["s1"],
                                           float(route.length))
                    m = m | (is_branch & inb)
                    claims += m.astype(int)
                    which[m] = znum[zi]
                    per.append({"name": p["name"], "s0": p["s0"], "s1": p["s1"],
                                "verts": int(m.sum())})
                    continue
                else:
                    f = RS.zone_field(station, on_route, p["s0"], p["s1"],
                                      route_end=float(route.length))
                    lo = p["s0"] if p["s0"] is not None else -np.inf
                    hi = p["s1"] if p["s1"] is not None else np.inf
                # HALF-OPEN. Testing the distal bound inclusively paints a
                # one-vertex crimson seam along every boundary two zones share --
                # which looks exactly like the overlap this panel detects, and is not.
                m = (f >= lo) & (f < hi)
                claims += m.astype(int)
                which[m] = znum[zi]
                per.append({"name": p["name"], "s0": p["s0"], "s1": p["s1"],
                            "verts": int(m.sum())})

            # FLATTEN EACH BOUNDARY AGAINST ITS OWN RING. The boundary already
            # comes out as one closed loop everywhere -- that was measured -- it is
            # just not flat: 5 to 11 mm off its best-fit plane on the arch against
            # 2 mm on the descending aorta, which is the harmonic level set meeting
            # the wall obliquely on a bend.
            #
            # Two earlier attempts used a plane built from the LANDMARK and both
            # shattered the boundary, because where k*d > 1 such a plane cuts the wall
            # in two or three places. A plane fitted to the ring itself cannot: it
            # sits on the wall where the boundary already is. With a geodesic collar
            # around the ring it stays local, and three passes take every boundary
            # measured to 1.5-2.3 mm.
            which, moved = RS.enforce_connected_labels(mesh, which)
            if use_u:
                adj = RS.surface_graph(mesh)
                # The aorta's own two ends are not up for re-decision by a boundary a
                # zone away: a flared root is geodesically close to a boundary 30 mm up
                # the vessel, and on KK40_3 that handed 182 root vertices to zone 1 and
                # then to zone 2.
                ends = RS.vessel_end_collar(mesh, route, adj)
                flat = []
                for zi in range(len(zones) - 1):
                    # ...and it is flattened onto THAT cross-section's plane, not
                    # onto the plane of wherever the boundary currently sits. Fitting
                    # the current ring makes it flat while locking in its displacement:
                    # flattening alone moved boundaries from +1.1 mm to +2.7 mm.
                    tgt = None
                    bb = zones[zi]["s1"]
                    if bb is not None and CROSS_SECTION_PLACEMENT:
                        base = float(np.interp(bb, ss_cal, uu_cal))
                        seg = s.tree.segments[route.locate(float(bb))[0]]
                        v = RS.boundary_u_at(mesh, route, u, bb, seg, prior=base)
                        ok_c = (v is not None and
                                abs(float(np.interp(v, uu_cal, ss_cal)) - float(bb))
                                <= CONTOUR_TOL_MM)
                        if ok_c:
                            tgt = RS.cross_section_plane(mesh, route, s.tree, bb)
                    # Confining the collar ALONG the vessel as well was measured and
                    # is not needed once the ends are protected: it fixed the same root
                    # and cost KR27_2's 1/2 boundary 2.74 -> 11.29 mm off plane and
                    # KK32_0's 0/1 3.17 -> 9.58.
                    which, dev = RS.flatten_boundary(
                        mesh, which, znum[zi], znum[zi + 1], adj, target=tgt,
                        frozen=ends)
                    if dev is not None:
                        flat.append(round(dev, 2))



            which, moved = RS.enforce_connected_labels(mesh, which)
            n_stub = 0
            if STUB_OWNERSHIP:
                g = g_adj
                if reg.max() >= 0:
                    before = which.copy()
                    # THE OSTIUM DECIDES, read straight off the centerline. Reading it
                    # off the region's own branch-called vertices is a proxy, and where
                    # there are none it fell back to the wall around the mouth -- which
                    # gave a bovine trunk's daughter the arch's zone on KK5_1 and
                    # KY244_1. Byte-identical on every non-bovine case measured.
                    which = RS.own_stubs(mesh, which, reg, is_branch, adj=g,
                                         region_zone=oz)
                    # ZONE 3 IS CUT BY TWO PLANES, not by a level set. Its proximal edge
                    # is a centerline-normal plane at the distal lip of the subclavian's
                    # mouth and its far edge is that plane 20 mm on -- the band the
                    # definition describes, and the transverse cut the clinical figure
                    # draws. A level set follows the wall, so on a bend it is oblique
                    # however it is placed, which is why moving it never fixed this.
                    import zone_stub_bench as _ZB
                    if _ZB.PLANE_BOUNDARY:
                        which = _ZB.plane_boundaries(s, mesh, route, which, reg, oz, g)
                    n_stub = int((which != before).sum())
            for zi in range(len(zones)):
                per[zi]["verts"] = int((which == znum[zi]).sum())
                per[zi]["zone"] = znum[zi]

            faces = np.asarray(mesh.faces)
            fclaim = claims[faces]
            n_over = int((fclaim > 1).all(axis=1).sum())
            n_gap = int((claims[faces] == 0).all(axis=1).sum())

            faces = np.asarray(mesh.faces)
            fclaim = claims[faces]
            n_over = int((fclaim > 1).all(axis=1).sum())
            n_gap = int((claims[faces] == 0).all(axis=1).sum())

            cols = np.tile(hex_rgb(UNCLAIMED), (nv, 1))
            for z in znum:
                cols[which == z] = hex_rgb(ZONE_COLOURS[z % len(ZONE_COLOURS)])
            cols[claims > 1] = hex_rgb(OVERLAP)

            poly = pv.PolyData(np.asarray(mesh.vertices),
                               np.hstack([np.full((len(faces), 1), 3),
                                          faces]).ravel())
            poly.point_data["rgb"] = cols.astype(np.uint8)
            pl = pv.Plotter(shape=(1, len(VIEWS)), off_screen=True, border=False,
                            window_size=(args.cell * len(VIEWS),
                                         int(args.cell * 1.3)))
            for c, v in enumerate(VIEWS):
                pl.subplot(0, c)
                pl.add_mesh(poly, scalars="rgb", rgb=True, smooth_shading=True,
                            specular=0.12, diffuse=0.92, ambient=0.3)
                pl.view_vector(v, viewup=UP)
                pl.reset_camera()
                pl.camera.zoom(1.4)
            img = pl.screenshot(None, transparent_background=True, return_img=True)
            pl.close()
            rows = np.nonzero((np.asarray(img)[..., 3] > 8).any(axis=1))[0]
            if len(rows):
                img = np.asarray(img)[max(int(rows.min()) - 6, 0):
                                      min(int(rows.max()) + 7,
                                          np.asarray(img).shape[0])]
            import imageio.v3 as iio
            png = outdir / f"{cid}.png"
            iio.imwrite(str(png), np.asarray(img))

            index.append({"case": cid, "png": png.name,
                          "base_ok": bool(st["ok"]),
                          "base_msg": "" if st["ok"]
                          else str(st.get("message", ""))[:160],
                          "route_mm": round(float(route.length), 1),
                          "zones": per,
                          "overlap_faces": n_over,
                          "gap_faces": n_gap,
                          "faces": int(len(faces)),
                          "repaired": int(moved),
                          "landmark_lag_mm": round(float(lag), 2),
                          "stub_forced": int(n_stub),
                          "harmonic": bool(use_u),
                          "flatness_mm": flat if use_u else [],
})
            print(f"  [{i:3d}/{len(cases)}] {cid:12s} {len(zones)} zone(s)  "
                  f"overlap {n_over:5d}  unclaimed {n_gap:6d}  "
                  f"{time.time() - t0:6.0f}s", flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  [{i:3d}/{len(cases)}] {cid:12s} ERROR {type(e).__name__}: {e}")
            index.append({"case": cid, "error": f"{type(e).__name__}: {e}"[:200]})

    (outdir / "index.json").write_text(
        json.dumps({"colours": ZONE_COLOURS, "unclaimed": UNCLAIMED,
                    "overlap": OVERLAP, "cases": index}, indent=1), encoding="utf-8")
    print(f"\nwrote {outdir}/index.json  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
