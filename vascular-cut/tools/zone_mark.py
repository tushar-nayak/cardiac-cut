"""Mark, by hand, where each zone boundary actually belongs.

    .venv\\Scripts\\python tools\\zone_mark.py
    .venv\\Scripts\\python tools\\zone_mark.py --case KK32_0
    .venv\\Scripts\\python tools\\zone_mark.py --boundary 1        # the 1/2 boundary

    right / space / n   next case          left / p   previous case
    click               mark where the boundary belongs
    u                   undo the last mark on this case
    x                   clear this case's marks
    g                   ghost the zone colours, so the mesh is easier to read
    r                   back to the starting orientation
    s                   screenshot
    q                   quit

WHY THIS EXISTS. Zone 3 starting too late has now survived three fixes, and each was
built on a measurement of mine that turned out to be measuring nothing: a coordinate
difference run through a calibration flagged six cases when the reviewer saw twenty; a
wall-distance flagged the right cases but returned 36 mm on one of them, which erased a
zone; a landmark-versus-mouth lag was contaminated by the same calibration that places
the boundaries, and correcting it made five of its eight cases worse.

The common fault is not arithmetic. Every one of those numbers was checked against
itself. The only time the loop closed was when the reviewer listed the bad cases by eye
and the metric disagreed with the list.

So this collects that judgement directly. Click on the wall where the boundary should
be, on as many cases as are worth marking, and the marks are written to
tools/zone_marks.json. tools/zone_score.py then scores any candidate rule by how far
its boundary lands from the marks -- so a rule is judged against where a person put the
line, and a change that improves the numbers while moving the line away is caught.

WHAT IS RECORDED. The picked point, the vertex nearest it, that vertex's harmonic
coordinate, and its arc length projected onto the centerline in 3D. Both coordinates,
because which of them a boundary should be specified in is exactly what is in dispute.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
from zone_map import ZONE_COLOURS, UNCLAIMED, UP, hex_rgb    # noqa: E402
from zone_stub_bench import label_case, ostium_zones         # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
MARKS = ROOT / "tools" / "zone_marks.json"
GHOST = "#b9c2cc"


def load(cid, cl, mp):
    """The labelling that ships, plus what is needed to record a mark against it."""
    s = CutSession(str(cl), str(mp))
    s.nub_preset()
    s.recut()
    mesh, route, zones, lab, isb, adj, _flat = label_case(s)
    reg = RS.stub_regions(mesh, route, adj=adj)
    oz = ostium_zones(s, mesh, route, reg) if reg.max() >= 0 else {}
    if reg.max() >= 0:
        lab = RS.own_stubs(mesh, lab, reg, isb, adj=adj, region_zone=oz)
    u = RS.harmonic_coordinate_capped(mesh, route, reg >= 0) if reg.max() >= 0 else None
    if u is None:
        u = RS.harmonic_coordinate(mesh, route)
    ss = np.arange(0.0, float(route.length) + 1e-9, 0.5)
    P = np.array([route.point(float(x)) for x in ss])
    znum = [int(p["name"].split()[-1]) for p in zones]
    return dict(mesh=mesh, lab=lab, reg=reg, oz=oz, u=np.asarray(u, float),
                ss=ss, P=P, zones=zones, znum=znum, route=route)


class Marker:
    def __init__(self, cases, boundary):
        self.cases = cases
        self.boundary = boundary          # mark the boundary between b and b+1
        self.i = 0
        self.ghost = False
        self.cache = {}
        self.marks = json.loads(MARKS.read_text()) if MARKS.exists() else {}
        self.p = pv.Plotter(window_size=(1100, 900))
        self.p.set_background("white")
        self.p.add_key_event("Right", lambda: self.step(+1))
        self.p.add_key_event("space", lambda: self.step(+1))
        self.p.add_key_event("n", lambda: self.step(+1))
        self.p.add_key_event("Left", lambda: self.step(-1))
        self.p.add_key_event("p", lambda: self.step(-1))
        self.p.add_key_event("u", self.undo)
        self.p.add_key_event("x", self.clear)
        self.p.add_key_event("g", self.toggle_ghost)
        self.p.add_key_event("r", self.home)
        self.p.add_key_event("s", self.shot)
        self.p.enable_surface_point_picking(
            callback=self.pick, show_message=False, show_point=True,
            color="black", point_size=14, left_clicking=True)

    # ---- data -----------------------------------------------------------
    def case(self):
        cid, cl, mp = self.cases[self.i]
        if cid not in self.cache:
            print(f"  loading {cid} ...", flush=True)
            self.cache[cid] = load(cid, cl, mp)
        return cid, self.cache[cid]

    def save(self):
        MARKS.write_text(json.dumps(self.marks, indent=1), encoding="utf-8")

    # ---- interaction ----------------------------------------------------
    def pick(self, point):
        cid, d = self.case()
        V = np.asarray(d["mesh"].vertices, float)
        j = int(np.argmin(np.linalg.norm(V - np.asarray(point, float), axis=1)))
        st = float(d["ss"][int(np.argmin(np.linalg.norm(d["P"] - V[j], axis=1)))])
        rec = {"boundary": f"{self.boundary}/{self.boundary + 1}",
               "point": [float(x) for x in V[j]], "vertex": j,
               "u": float(d["u"][j]), "station_mm": st,
               "zone_here": int(d["lab"][j])}
        self.marks.setdefault(cid, []).append(rec)
        self.save()
        print(f"    {cid}: marked {self.boundary}/{self.boundary + 1} at "
              f"u {rec['u']:.4f}, station {st:.1f} mm "
              f"(currently zone {rec['zone_here']})", flush=True)
        self.draw()

    def undo(self):
        cid, _ = self.case()
        if self.marks.get(cid):
            self.marks[cid].pop()
            self.save()
            print(f"    {cid}: undid the last mark", flush=True)
            self.draw()

    def clear(self):
        cid, _ = self.case()
        if self.marks.pop(cid, None) is not None:
            self.save()
            print(f"    {cid}: cleared", flush=True)
            self.draw()

    def step(self, d):
        self.i = (self.i + d) % len(self.cases)
        self.draw()

    def toggle_ghost(self):
        self.ghost = not self.ghost
        self.draw()

    def home(self):
        self.p.view_vector(np.array([0.4, -1.0, 0.15]), viewup=UP)
        self.p.reset_camera()

    def shot(self):
        cid, _ = self.case()
        out = ROOT / "tools" / f"mark_{cid}.png"
        self.p.screenshot(str(out))
        print(f"    wrote {out}", flush=True)

    # ---- drawing --------------------------------------------------------
    def draw(self):
        cid, d = self.case()
        cam = self.p.camera_position
        self.p.clear()
        mesh, lab, reg = d["mesh"], d["lab"], d["reg"]
        V = np.asarray(mesh.vertices, float)
        F = np.asarray(mesh.faces, int)
        e = np.asarray(mesh.edges_unique, int)
        cols = np.tile(hex_rgb(GHOST if self.ghost else UNCLAIMED), (len(V), 1))
        if not self.ghost:
            for z in range(int(lab.max()) + 1):
                cols[lab == z] = hex_rgb(ZONE_COLOURS[z % len(ZONE_COLOURS)])
        b = self.boundary
        on = (((lab[e[:, 0]] == b) & (lab[e[:, 1]] == b + 1))
              | ((lab[e[:, 0]] == b + 1) & (lab[e[:, 1]] == b)))
        if on.any():
            cols[np.unique(e[on])] = np.array([255, 255, 255], float)
        owner = [j for j, z in d["oz"].items() if z == b]
        if owner:
            sel = np.isin(reg, owner)
            cut = sel[e[:, 0]] != sel[e[:, 1]]
            cols[np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))] = np.zeros(3)
        poly = pv.PolyData(V, np.hstack([np.full((len(F), 1), 3), F]).ravel())
        poly.point_data["rgb"] = cols.astype(np.uint8)
        self.p.add_mesh(poly, scalars="rgb", rgb=True, smooth_shading=True,
                        specular=0.12, diffuse=0.92, ambient=0.3, pickable=True)
        got = self.marks.get(cid, [])
        if got:
            pts = np.array([m["point"] for m in got], float)
            self.p.add_points(pts, color="#e11d48", point_size=16,
                              render_points_as_spheres=True, pickable=False)
        self.p.add_text(
            f"{cid}   [{self.i + 1}/{len(self.cases)}]   boundary {b}/{b + 1}   "
            f"{len(got)} mark(s)\n"
            f"white = boundary now,  black = mouth of the vessel zone {b} owns\n"
            f"click where it belongs   n/p step   u undo   x clear   g ghost   q quit",
            font_size=9, color="#222222")
        if cam is not None and self.p.camera_position is not None:
            self.p.camera_position = cam
        self.p.render()

    def run(self):
        self.draw()
        self.home()
        self.p.show()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--case", default="")
    ap.add_argument("--boundary", type=int, default=2,
                    help="mark the boundary between this zone and the next")
    args = ap.parse_args()
    cases = find_cases(args.cohort)
    if args.case:
        want = set(args.case.split(","))
        cases = [c for c in cases if c[0] in want]
    if not cases:
        raise SystemExit("no cases")
    print(f"{len(cases)} case(s); marks go to {MARKS}")
    Marker(cases, args.boundary).run()


if __name__ == "__main__":
    main()
