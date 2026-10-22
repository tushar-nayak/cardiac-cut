"""Click through every arch-and-nub geometry and record a verdict on each.

    .venv\\Scripts\\python arch_viewer.py                  # output/arch_final
    .venv\\Scripts\\python arch_viewer.py --dir output/arch_nubs   # the batch's cuts
    .venv\\Scripts\\python arch_viewer.py --unjudged        # only ones not yet judged

Same window and the same keys as nub_viewer, which this reuses outright -- one case
at a time, rotate freely, every opening ringed and labelled:

    right / space / n     next case          g   good, keep it
    left / p              previous case      b   bad, needs work
    f   next unjudged                        u   undecided (clears a verdict)
    r   back to the starting orientation     s   save a screenshot
    q   quit

WHAT THIS SHOWS THAT nub_viewer DOES NOT

The two REGION BOUNDARIES. An arch export lists its branch openings in the sidecar
and its boundaries only as stations, so a viewer reading the sidecar alone rings
three nubs on a mesh with five openings -- and the two it leaves out are the ones
under review. Where the arch starts and ends is the open question; the nubs are
settled.

So the boundaries are found in the MESH rather than read from the file: every
boundary loop that no listed opening claims is a boundary of the region, ringed in
blue and labelled with the station and calibre the sidecar recorded for it. Nothing
has to be re-cut to see them.

Verdicts go to review.json in the same folder, written on every keystroke.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import nub_viewer as NV                                            # noqa: E402
from branch_cut import boundary_loops                              # noqa: E402

BOUND = "#2f6fd0"          # a region boundary: the thing being judged here
MATCH_MM = 6.0             # how near a loop's centre must be to claim a listed opening


def _loops_of(stl):
    """Every boundary loop of the mesh, as (centroid, mean radius, count).

    MERGED FIRST, and that is not optional. An STL stores every triangle with its own
    three vertices, so loaded raw the mesh is a soup with no shared edges at all and
    every edge reads as a boundary: DJ18_2 came back with 9,494 "openings" instead of
    5. Welding coincident vertices is what turns it back into a surface.
    """
    m = trimesh.load(str(stl), process=False)
    m.merge_vertices()
    v = np.asarray(m.vertices)
    out = []
    for L in boundary_loops(m):
        P = v[np.asarray(L)]
        c = P.mean(axis=0)
        out.append((c, float(np.linalg.norm(P - c, axis=1).mean()), len(L)))
    return out


def _bounds_from_mesh(stl, rec, region):
    """The rims no listed opening accounts for: the region's own two boundaries.

    Read off the geometry, not the sidecar, because the sidecar does not carry them.
    A listed opening claims the loop nearest its recorded centroid; whatever is left
    is a boundary. Named start or end by matching the calibre the region block
    recorded, which is the only thing in the file that tells them apart.
    """
    listed = [v.get("rim_centroid") for v in rec["vessels"]
              if v.get("rim_centroid") is not None]
    free = []
    for c, rad, n in _loops_of(stl):
        if listed and min(float(np.linalg.norm(np.asarray(x, float) - c))
                          for x in listed) <= MATCH_MM:
            continue
        free.append((c, rad, n))
    if not free:
        return []

    want = {"Region start": region.get("start_diameter_mm"),
            "Region end": region.get("end_diameter_mm")}
    out, used = [], set()
    for label, dia in want.items():
        if dia is None or not free:
            continue
        j = min((k for k in range(len(free)) if k not in used),
                key=lambda k: abs(free[k][1] * 2.0 - float(dia)), default=None)
        if j is None:
            continue
        used.add(j)
        c, rad, n = free[j]
        out.append({"name": label, "ok": True, "cut_ok": True, "bound": True,
                    "rim_centroid": [float(x) for x in c],
                    "mouth_ead_mm": round(rad * 2.0, 2),
                    "offset_mm": region.get(
                        "start_mm_along_aorta" if label == "Region start"
                        else "end_mm_along_aorta"),
                    "cut_mode": None, "hand_adjusted": False})
    for k, (c, rad, n) in enumerate(free):
        if k in used:
            continue
        out.append({"name": "Region boundary", "ok": True, "cut_ok": True,
                    "bound": True, "rim_centroid": [float(x) for x in c],
                    "mouth_ead_mm": round(rad * 2.0, 2),
                    "cut_mode": None, "hand_adjusted": False})
    return out


LSA_WORDS = ("left subclavian",)


def _jstations(tree, route, jg):
    """Every junction that sits ON the route, and its station."""
    out = {}
    for i, name in enumerate(route.chain):
        ns, ne = jg[name]
        seg = tree.segments[name]
        a, b = (ns, ne) if route.forward[i] else (ne, ns)
        if a is not None:
            out.setdefault(a, route.starts[i])
        if b is not None:
            out.setdefault(b, route.starts[i] + seg.length)
    return out


def lsa_distances(csv, mat, sidecar, batch=None):
    """Distance along the centerline from the LEFT SUBCLAVIAN's node to each cut.

    The node is where the left subclavian meets the aorta -- its ostium on the route,
    walked out through whatever it leaves through, so a subclavian arising from a
    trunk is measured from where that TRUNK meets the aorta and not from a junction
    floating off the route.

    Two components, because a branch cut is not on the aorta at all:

        along_route   from the LSA node to where this vessel leaves the aorta.
                      Signed: + is distal, toward the descending aorta.
        along_branch  from that ostium out to the cut plane. Zero for a cut that is
                      itself on the aorta, which is what the two region bounds are.

    The distance quoted is the PATH: |along_route| + along_branch, the length a walk
    down the centerline would cover. The signed route component is carried too,
    because "18 mm from the subclavian" means something different proximally and
    distally.

    Offsets come from the EXPORT, not from the batch, so a plane moved by hand in the
    workbench is measured where it actually ended up.
    """
    from vascular_api import CutSession
    import branch_cut as BC

    d = json.loads(Path(sidecar).read_text(encoding="utf-8"))
    by_name, interior = {}, []
    for c in ((batch or {}).get("cuts") or []):
        by_name.setdefault(str(c.get("name")), c.get("segment"))
        if c.get("interior") and not c.get("bound") and c.get("segment"):
            interior.append(c["segment"])
    s_sess = CutSession(str(csv), str(mat), family="thoracic")
    route, tree = s_sess.route, s_sess.tree
    if route is None:
        return None, []
    jg = BC._junction_graph(tree)
    jst = _jstations(tree, route, jg)

    def ostium(seg_name):
        ns, ne = jg.get(seg_name, (None, None))
        for j in (ns, ne):
            st = BC._station_on_route(jg, jst, j)
            if st is not None:
                return float(st)
        return None

    lsa = next((t for t in tree.terminals
                if any(w in str(t.name).lower() for w in LSA_WORDS)), None)
    if lsa is None:
        return None, []
    lsa_s = ostium(lsa.segment.name)
    if lsa_s is None:
        return None, []

    rows = []
    for o in d.get("openings", []):
        nm = str(o.get("name") or "")
        off = o.get("offset_mm")
        # The two boundaries are reported below, from the stations the export
        # recorded. Here they are openings with no terminal behind them -- their
        # offset is measured along whichever segment the route passes through, not
        # along the route -- so looking them up produced a "not found" line directly
        # above the real answer, on every case.
        if nm in ("Region start", "Region end"):
            continue
        t = next((x for x in tree.terminals if str(x.name) == nm), None)
        if t is None and off is not None:
            # A TRIM is an interior cut, not a terminal, and the export does not
            # record which segment it sits on -- only the name of the group it stands
            # for. The batch record does, so the segment comes from there and the
            # OFFSET still comes from the export, which is where the reviewer left it.
            seg_name = by_name.get(nm)
            if seg_name is None and interior:
                # The window renames a trim after the group it swallows -- "Trim
                # before Right common carotid + Right subclavian" -- while the batch
                # kept the vessel's own name, so the two do not match on text. An arch
                # carries at most one trim, so an unmatched opening that the batch
                # recorded as interior is that one.
                seg_name = interior.pop(0)
            if seg_name:
                try:
                    lbl = s_sess.add_trim(seg_name)
                except Exception:                                   # noqa: BLE001
                    lbl = None
                if lbl:
                    t = next((x for x in s_sess._cut_list()[0]
                              if x.label == lbl), None)
        if t is None or off is None:
            rows.append({"name": nm, "along_route": None, "along_branch": None,
                         "path": None,
                         "note": "not found in the tree or the batch record"})
            continue
        s_cut = float(t.s_cut(float(off)))
        if t.segment.name in route.chain:
            st = float(route.station_of(t.segment.name, s_cut))
            rows.append({"name": nm, "along_route": st - lsa_s,
                         "along_branch": 0.0, "path": abs(st - lsa_s)})
            continue
        ost = ostium(t.segment.name)
        if ost is None:
            rows.append({"name": nm, "along_route": None, "along_branch": None,
                         "path": None, "note": "does not reach the aorta"})
            continue
        ns, ne = jg[t.segment.name]
        L = float(t.segment.length)
        # Which end of this segment is the junction decides which way to measure.
        at_start = ns is not None and BC._station_on_route(jg, jst, ns) is not None
        along = s_cut if at_start else (L - s_cut)
        rows.append({"name": nm, "along_route": ost - lsa_s,
                     "along_branch": float(along),
                     "path": abs(ost - lsa_s) + float(along)})

    # The two region bounds, from the stations the export recorded for them.
    reg = d.get("region") or {}
    for label, key in (("Region start", "start_mm_along_aorta"),
                       ("Region end", "end_mm_along_aorta")):
        if reg.get(key) is not None:
            st = float(reg[key])
            rows.append({"name": label, "along_route": st - lsa_s,
                         "along_branch": 0.0, "path": abs(st - lsa_s)})
    return lsa_s, rows


def load_cases(folder, unjudged_only=False, review=None):
    out = []
    for js in sorted(folder.glob("*_cut_planes.json")):
        stl = js.with_name(js.name.replace("_cut_planes.json", ".stl"))
        if not stl.exists():
            continue
        rec = NV._from_export(js, stl)
        rec["case"] = stl.stem.replace("_arch", "")
        d = json.loads(js.read_text(encoding="utf-8"))
        region = d.get("region") or {}
        rec["region"] = region
        rec["vessels"] = list(rec["vessels"]) + _bounds_from_mesh(stl, rec, region)
        for n in (region.get("boundary_notes") or []):
            rec.setdefault("notes", []).append(n)
        if unjudged_only and review is not None and review.get(rec["case"]):
            continue
        out.append((stl, rec))
    # The batch's own sidecars, for a folder that has not been through the workbench.
    for js in sorted(folder.glob("*_arch_k*.json")):
        stl = js.with_suffix(".stl")
        if not stl.exists() or js.name == "review.json":
            continue
        d = json.loads(js.read_text(encoding="utf-8"))
        rec = {"case": stl.stem.split("_arch")[0], "openings": d.get("openings"),
               "euler": d.get("euler"), "area_mm2": d.get("area_mm2"),
               "export": False, "notes": list(d.get("notes") or []),
               "region": d.get("region") or {}, "vessels": []}
        if unjudged_only and review is not None and review.get(rec["case"]):
            continue
        out.append((stl, rec))
    return out


class ArchViewer(NV.Viewer):
    """nub_viewer's window, with the arch's own numbers in the corners."""

    def _header(self, rec):
        r = rec.get("region") or {}
        span = ""
        if r.get("start_mm_along_aorta") is not None:
            span = (f"   arch {r['start_mm_along_aorta']:.0f} to "
                    f"{r.get('end_mm_along_aorta', 0):.0f} mm "
                    f"({r.get('length_mm', 0):.0f} mm of a "
                    f"{r.get('aorta_route_length_mm', 0):.0f} mm aorta)")
        v = self.review.get(rec["case"])
        mark = {"good": "  GOOD", "bad": "  BAD"}.get(v, "")
        return (f"[{self.i + 1}/{len(self.cases)}]  {rec['case']}{mark}\n"
                f"{rec.get('openings', '?')} openings, "
                f"{rec.get('area_mm2', 0):,.0f} mm2{span}")

    def _footer(self, rec):
        bits = []
        r = rec.get("region") or {}
        if r.get("max_diameter_mm"):
            bits.append(f"widest {r['max_diameter_mm']:.0f} mm")
        if r.get("branches_inside"):
            bits.append(f"{len(r['branches_inside'])} branch(es) kept")
        base = "; ".join(bits)
        notes = rec.get("notes") or []
        if notes:
            base += ("\n" if base else "") + "\n".join(str(n)[:110] for n in notes[:3])
        good, bad = self.review.counts()
        return (base + f"\n\ng good   b bad   u undecided   f next unjudged   "
                       f"r reset view   s shot   q quit      "
                       f"[{good} good, {bad} bad]")


def cmd_distances(folder, run, csv_out):
    """Print, for every case, how far each cut sits from the left subclavian's node."""
    batches = {}
    if run.exists():
        for line in run.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                batches[r["case"].split("/")[-1]] = r
    sides = sorted(folder.glob("*_cut_planes.json"))
    if not sides:
        print(f"no exports in {folder}")
        return 1
    rows_out = []
    print("\n  Distance along the centerline from the LEFT SUBCLAVIAN's node.")
    print("  route  = along the aorta to where that vessel leaves it (+ is distal)")
    print("  branch = out along the vessel to the cut; 0 for a cut on the aorta")
    print("  path   = |route| + branch, the walk down the centerline\n")
    for side in sides:
        cid = side.name.replace("_arch_cut_planes.json", "")
        b = batches.get(cid)
        if b is None or not b.get("csv"):
            print(f"  {cid:10s} no batch record, cannot locate the centerline")
            continue
        try:
            lsa, rows = lsa_distances(b["csv"], b["mat"], side, batch=b)
        except Exception as e:                                      # noqa: BLE001
            print(f"  {cid:10s} {type(e).__name__}: {str(e)[:60]}")
            continue
        if lsa is None:
            print(f"  {cid:10s} no left subclavian found")
            continue
        print(f"  {cid}   LSA node at {lsa:.1f} mm along the aorta")
        for x in rows:
            if x["path"] is None:
                print(f"      {x['name'][:34]:34s}  --  {x.get('note')}")
                continue
            print(f"      {x['name'][:34]:34s}  route {x['along_route']:+7.1f}   "
                  f"branch {x['along_branch']:6.1f}   path {x['path']:7.1f} mm")
            rows_out.append((cid, x["name"], x["along_route"],
                             x["along_branch"], x["path"], lsa))
    if rows_out:
        head = "case,cut,along_route_mm,along_branch_mm,path_mm,lsa_node_mm"
        body = "\n".join(f"{c},{n},{r:.2f},{br:.2f},{p:.2f},{l:.2f}"
                         for c, n, r, br, p, l in rows_out)
        csv_out.write_text(head + "\n" + body + "\n", encoding="utf-8")
        print(f"\n  {len(rows_out)} cut(s) over {len(sides)} case(s)  ->  {csv_out}")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="output/arch_final")
    ap.add_argument("--unjudged", action="store_true")
    ap.add_argument("--distances", action="store_true",
                    help="print the distance along the centerline from the left "
                         "subclavian's node to every cut, and write a CSV")
    ap.add_argument("--run", default=str(ROOT / "tools" / "arch_run.jsonl"),
                    help="the batch record, read only to resolve a trim's segment")
    ap.add_argument("--csv", default=str(ROOT / "tools" / "arch_lsa_distances.csv"))
    args = ap.parse_args()

    folder = Path(args.dir)
    if not folder.exists():
        print(f"no such folder: {folder}")
        return 1
    if args.distances:
        return cmd_distances(folder, Path(args.run), Path(args.csv))
    review = NV.Review(folder / "review.json")
    cases = load_cases(folder, unjudged_only=args.unjudged, review=review)
    if not cases:
        print(f"nothing to show in {folder}")
        return 1
    # Blue for a region boundary, so it reads apart from the nubs at a glance.
    NV._rim_colour = (lambda v, _f=NV._rim_colour:
                      BOUND if v.get("bound") else _f(v))
    print(f"{len(cases)} case(s) from {folder}")
    ArchViewer(cases, review, folder).run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
