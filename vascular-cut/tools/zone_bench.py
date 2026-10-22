"""Every zone preset, cut BOTH ways, over a cohort. Does the station cut hold up?

    .venv\\Scripts\\python tools\\zone_bench.py                      # thoracic (GT)
    .venv\\Scripts\\python tools\\zone_bench.py Y:/aorta/abdominal/MK/JK
    .venv\\Scripts\\python tools\\zone_bench.py <dir> --station-only

For each preset a case offers, two answers are recorded side by side:

  PLANAR   what the tool does today -- two planes, nudged clear of every ostium, then
           walked along a ladder of alternative stations until one cuts flat. Recorded
           as accepted/refused, and as the DRIFT: how far the boundary ended up from
           the station the anatomy asked for.
  STATION  region_station.region_by_station -- the zone as a level set of arc length.
           Recorded as accepted/refused and checked for the things that make a surface
           usable afterwards: winding, zero-area faces, boundary loops, and whether any
           original wall from outside the span survived.

The comparison that matters is not "does it cut" -- the station cut always cuts, which
is most of the point -- but whether what it produces is a sane zone. So the geometry is
measured, not just counted.

Results stream to tools/zone_bench_<cohort>.jsonl as they are produced, because a
30-case run takes tens of minutes and a crash on case 27 should not cost the first 26.
"""
import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import branch_cut as BC                                      # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
STATION_TOL_MM = 1.0


def bound_station(route, t):
    """Where a planar boundary actually sits, as a station on the route."""
    seg = t.segment
    local = (seg.length - t._region_offset) if t.at_end else t._region_offset
    return float(route.station_of(seg.name, float(local)))


def planar_drift(tree, route, jg, mesh, s0, s1):
    """How far region_plan's boundaries end up from the stations asked for."""
    try:
        bounds, _in, _rm, notes = BC.region_plan(tree, s0, s1, route, jg, mesh)
    except Exception as e:                                   # noqa: BLE001
        return None, None, [f"region_plan: {type(e).__name__}: {e}"]
    got = sorted(bound_station(route, t) for t in bounds)
    ask = [x for x in (s0, s1) if x is not None]
    d = [got[i] - ask[i] for i in range(min(len(got), len(ask)))]
    return (d[0] if d else None), (d[1] if len(d) > 1 else None), list(notes)


def station_report(base_mesh, tree, route, s0, s1):
    """Cut the zone by station and measure what came out."""
    m, info = RS.region_by_station(base_mesh, tree, s0, s1, route=route)
    st, _own, on_route = RS.station_field(m, tree, route)
    d_orig, _ = cKDTree(np.asarray(base_mesh.vertices, float)).query(
        np.asarray(m.vertices, float))
    orig = on_route & (d_orig < 1e-9)          # interpolated rim vertices excluded
    wall = st[orig]
    a = m.area_faces
    info.update({
        "winding_ok": bool(m.is_winding_consistent),
        "degenerate": int((a < 1e-9).sum()),
        "wall_min": None if not len(wall) else float(wall.min()),
        "wall_max": None if not len(wall) else float(wall.max()),
        "overrun_prox": (None if not len(wall) or s0 is None
                         else float(max(0.0, s0 - wall.min()))),
        "overrun_dist": (None if not len(wall) or s1 is None
                         else float(max(0.0, wall.max() - s1))),
    })
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--station-only", action="store_true",
                    help="skip the planar path (much faster, no comparison)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    cases = find_cases(args.cohort)
    if args.limit:
        cases = cases[:args.limit]
    tag = args.tag or Path(args.cohort).name
    out_path = ROOT / "tools" / f"zone_bench_{tag}.jsonl"
    fh = out_path.open("w", encoding="utf-8")
    print(f"{len(cases)} case(s) under {args.cohort}")
    print(f"-> {out_path}\n")

    t0 = time.time()
    for i, (cid, cl, mp) in enumerate(cases, 1):
        rec = {"case": cid}
        try:
            s = CutSession(str(cl), str(mp))
            st = s.recut()
            rec["base_ok"] = bool(st["ok"])
            rec["base_msg"] = "" if st["ok"] else str(st.get("message", ""))[:120]
            tree, route = s.tree, BC.aorta_route(s.tree)
            jg = BC._junction_graph(tree)
            # SNAPSHOT THE WHOLE-VESSEL MESH. s.mesh is whatever the last recut
            # produced, and clear_region() does not put the full surface back until
            # something re-cuts -- so reading s.mesh inside the loop hands zone N+1
            # the mesh of zone N, already clipped. That scored 54 of 125 zones EMPTY
            # and read as a failure of the method rather than of this loop.
            whole = s.mesh
            rec["route_mm"] = None if route is None else float(route.length)
            rec["n_trunk_ends"] = sum(1 for t in tree.terminals
                                      if getattr(t, "is_trunk", False))
            presets = s.presets()
            rec["presets"] = [p["name"] for p in presets]
            zones = []
            for p in presets:
                z = {"name": p["name"], "s0": p["s0"], "s1": p["s1"]}
                try:
                    z["station"] = station_report(whole, tree, route,
                                                  p["s0"], p["s1"])
                    # Not raising is not succeeding. An empty zone is the failure
                    # mode a level-set cut actually has, and it has to be counted.
                    z["station_ok"] = bool(z["station"]["faces"] > 0)
                except Exception as e:                       # noqa: BLE001
                    z["station_ok"] = False
                    z["station_err"] = f"{type(e).__name__}: {e}"[:160]
                if not args.station_only:
                    try:
                        s.apply_preset(p["name"])
                        r = s.recut()
                        z["planar_ok"] = bool(r["ok"])
                        z["planar_msg"] = ("" if r["ok"]
                                           else str(r.get("message", ""))[:120])
                    except Exception as e:                   # noqa: BLE001
                        z["planar_ok"] = False
                        z["planar_msg"] = f"{type(e).__name__}: {e}"[:120]
                    finally:
                        s.clear_region()
                    dp, dd, notes = planar_drift(tree, route, jg, s.base,
                                                 p["s0"], p["s1"])
                    z["planar_drift_prox"] = dp
                    z["planar_drift_dist"] = dd
                    z["planar_notes"] = notes[:6]
                zones.append(z)
            rec["zones"] = zones
        except Exception as e:                               # noqa: BLE001
            rec["error"] = f"{type(e).__name__}: {e}"[:200]
            rec["trace"] = traceback.format_exc()[-400:]
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        nz = len(rec.get("zones", []))
        sok = sum(1 for z in rec.get("zones", []) if z.get("station_ok"))
        pok = sum(1 for z in rec.get("zones", []) if z.get("planar_ok"))
        print(f"  [{i:3d}/{len(cases)}] {cid:12s} {nz} zone(s)  "
              f"station {sok}/{nz}  planar {pok}/{nz}"
              f"{'  ERROR ' + rec['error'] if 'error' in rec else ''}"
              f"   {time.time() - t0:6.0f}s", flush=True)
    fh.close()
    print(f"\nwrote {out_path}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
