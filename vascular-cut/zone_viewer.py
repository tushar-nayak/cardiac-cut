"""Click through every case's zones and rotate them, one window at a time.

    .venv\\Scripts\\python zone_viewer.py                       # thoracic cohort
    .venv\\Scripts\\python zone_viewer.py --cohort Y:/aorta/abdominal/MK/JK
    .venv\\Scripts\\python zone_viewer.py --case KK24_1         # straight to one

    right / space / n     next case
    left / p              previous case
    a                     all zones at once, the map
    1 2 3 4 5             show one zone, the rest ghosted
    i                     islands on/off -- the pieces a zone would have to discard
    g                     ghost of the whole vessel on/off
    r                     back to the starting orientation
    s                     save a screenshot
    q                     quit

Same shape as nub_viewer.py deliberately: one window, one case, free rotation, and the
next case arrives in the same place at the same size in whatever orientation the last
one was left in -- so a view worth having is set once and each case is a keystroke.

WHY ISLANDS ARE A MODE OF THEIR OWN. A zone has to be one continuous piece; a zone
delivered as two patches on opposite walls describes neither, and any area or delta K
taken over it averages across the gap. region_by_station enforces that by keeping the
largest piece, which means the failure is invisible in the result -- the discarded
patch is simply gone. Here it is drawn, in crimson, on the case it came from, because
whether that patch is a defect in the station field or a genuinely separate scrap of
wall is a question that can only be answered by looking at it.

The zone colouring is computed on load and cached to output/zone_cache/, so stepping
back through cases is instant. Delete that folder to force a recompute.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

THORACIC = r"Z:/Martino/branches_CAMP/GT"
CACHE = ROOT / "output" / "zone_cache"

# The five that tile, in order. Matches tools/zone_map.py so a picture here and a
# picture there are the same picture.
ZONE_COLOURS = ["#1d4ed8", "#0e9488", "#b8860b", "#c2410c", "#7d3ac1"]
UNCLAIMED = "#d9d4cf"
ISLAND = "#e11d48"
GHOST = "#e8e4df"
TILING = ("Ishimaru zone 0", "Ishimaru zone 1", "Ishimaru zone 2",
          "Ishimaru zone 3", "Ishimaru zone 4")
DEFAULT_VIEW = np.array([0.4, -1.0, 0.15])
DEFAULT_UP = np.array([0.0, 0.0, 1.0])


def hex_rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], np.uint8)


def face_components(faces, keep):
    """Component label per KEPT face, using shared-edge adjacency."""
    idx = np.nonzero(keep)[0]
    if len(idx) == 0:
        return np.array([], int), 0
    sub = faces[idx]
    e = np.sort(np.vstack([sub[:, [0, 1]], sub[:, [1, 2]], sub[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(sub)), 3)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, owner = e[order], owner[order]
    same = np.all(e[1:] == e[:-1], axis=1)
    a, b = owner[:-1][same], owner[1:][same]
    if len(a) == 0:
        return np.arange(len(sub)), len(sub)
    g = sp.coo_matrix((np.ones(len(a)), (a, b)), shape=(len(sub), len(sub)))
    k, lab = connected_components(g, directed=False)
    return lab, k


def build(cid, cl, mp):
    """Zone label and island mask for one case. Cached to disk."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{cid}.npz"
    if f.exists():
        d = np.load(f, allow_pickle=True)
        return (d["V"], d["F"], d["label"], d["island"],
                list(d["names"]), d["spans"], bool(d["base_ok"]))

    from vascular_api import CutSession
    import region_station as RS
    import branch_cut as BC

    s = CutSession(str(cl), str(mp))
    s.nub_preset()
    st = s.recut()
    mesh = s.mesh
    route = BC.aorta_route(s.tree)
    station, _own, on_route = RS.station_field(mesh, s.tree, route)
    zones = sorted((p for p in s.presets() if p["name"] in TILING),
                   key=lambda p: p["s0"])

    V = np.asarray(mesh.vertices, float)
    F = np.asarray(mesh.faces, int)
    label = -np.ones(len(V), int)
    island = np.zeros(len(F), bool)
    names, spans = [], []
    for zi, p in enumerate(zones):
        fld = RS.zone_field(station, on_route, p["s0"], p["s1"],
                            route_end=float(route.length))
        lo = p["s0"] if p["s0"] is not None else -np.inf
        hi = p["s1"] if p["s1"] is not None else np.inf
        m = (fld >= lo) & (fld < hi)              # half-open, as region_station
        label[m] = zi
        # Islands are judged on FACES, because that is the surface that would be
        # kept or discarded, and a lone vertex on the far side of a boundary is not
        # a piece of anything.
        names.append(p["name"])
        spans.append([p["s0"], -1.0 if p["s1"] is None else p["s1"]])
    # Repair first, THEN look for islands: what is left after reassignment is a piece
    # that has nowhere to go, which is the thing worth drawing in crimson.
    label, moved = RS.enforce_connected_labels(mesh, label)
    island[:] = False
    for zi in range(len(names)):
        keep = (label[F] == zi).all(axis=1)
        lab, k = face_components(F, keep)
        if k > 1:
            idx = np.nonzero(keep)[0]
            sizes = np.bincount(lab)
            island[idx[lab != int(np.argmax(sizes))]] = True
    spans = np.asarray(spans, float) if spans else np.zeros((0, 2))
    np.savez_compressed(f, V=V, F=F, label=label, island=island,
                        names=np.array(names, dtype=object), spans=spans,
                        base_ok=bool(st["ok"]))
    return V, F, label, island, names, spans, bool(st["ok"])


class Viewer:
    def __init__(self, cases, off_screen=False):
        self.cases = cases
        self.i = 0
        self.mode = "all"            # "all" or a zone index
        self.show_islands = True
        self.show_ghost = True
        self.cam = None
        self.p = pv.Plotter(window_size=(1400, 900), off_screen=off_screen)
        self.p.set_background("white")
        self.p.add_key_event("Right", lambda: self.step(+1))
        self.p.add_key_event("space", lambda: self.step(+1))
        self.p.add_key_event("n", lambda: self.step(+1))
        self.p.add_key_event("Left", lambda: self.step(-1))
        self.p.add_key_event("p", lambda: self.step(-1))
        self.p.add_key_event("a", lambda: self.set_mode("all"))
        for k in range(1, 6):
            self.p.add_key_event(str(k), lambda k=k: self.set_mode(k - 1))
        self.p.add_key_event("i", self.toggle_islands)
        self.p.add_key_event("g", self.toggle_ghost)
        self.p.add_key_event("r", self.home)
        self.p.add_key_event("s", self.shot)
        self.draw()

    # -- state ------------------------------------------------------------
    def step(self, d):
        self.cam = self.p.camera_position
        self.i = (self.i + d) % len(self.cases)
        self.draw()

    def set_mode(self, m):
        self.cam = self.p.camera_position
        self.mode = m
        self.draw()

    def toggle_islands(self):
        self.cam = self.p.camera_position
        self.show_islands = not self.show_islands
        self.draw()

    def toggle_ghost(self):
        self.cam = self.p.camera_position
        self.show_ghost = not self.show_ghost
        self.draw()

    def home(self):
        self.cam = None
        self.draw()

    def shot(self):
        cid = self.cases[self.i][0]
        out = ROOT / "output" / f"zone_view_{cid}.png"
        self.p.screenshot(str(out))
        print(f"wrote {out}")

    # -- drawing ----------------------------------------------------------
    def draw(self):
        cid, cl, mp = self.cases[self.i]
        self.p.clear()
        try:
            V, F, label, island, names, spans, base_ok = build(cid, cl, mp)
        except Exception as e:                                # noqa: BLE001
            self.p.add_text(f"{cid}\n{type(e).__name__}: {e}", font_size=11,
                            color="#b42318")
            self.p.render()
            return

        poly = pv.PolyData(V, np.hstack([np.full((len(F), 1), 3), F]).ravel())
        # Face colours, not vertex: a face belongs to a zone only if all three of its
        # corners do, so the seam between zones is a clean line rather than a blended
        # band that makes a boundary look uncertain when it is not.
        fl = np.where((label[F] == label[F][:, [0]]).all(axis=1), label[F][:, 0], -1)
        cols = np.tile(hex_rgb(UNCLAIMED), (len(F), 1))
        for zi in range(len(names)):
            if self.mode == "all" or self.mode == zi:
                cols[fl == zi] = hex_rgb(ZONE_COLOURS[zi % len(ZONE_COLOURS)])
            else:
                cols[fl == zi] = hex_rgb(GHOST)
        if self.show_islands:
            live = island if self.mode == "all" else (island & (fl == self.mode))
            cols[live] = hex_rgb(ISLAND)
        poly.cell_data["rgb"] = cols
        self.p.add_mesh(poly, scalars="rgb", rgb=True, smooth_shading=True,
                        specular=0.12, diffuse=0.92, ambient=0.3)
        if self.show_ghost:
            self.p.add_mesh(poly.copy(), color=GHOST, opacity=0.12,
                            smooth_shading=True, ambient=0.6)

        n_isl = int(island.sum())
        head = (f"{cid}   [{self.i + 1}/{len(self.cases)}]"
                f"{'' if base_ok else '   BRANCH CUT REFUSED'}")
        lines = [head, ""]
        for zi, nm in enumerate(names):
            s0, s1 = spans[zi]
            end = "end" if s1 < 0 else f"{s1:.0f}"
            mark = ">" if self.mode == zi else " "
            n_f = int((fl == zi).sum())
            isl = int((island & (fl == zi)).sum())
            lines.append(f"{mark} {nm.replace('Ishimaru ', ''):8s} "
                         f"{s0:6.0f}-{end:>5s} mm  {n_f:6d} faces"
                         + (f"   island {isl}" if isl else ""))
        lines += ["", f"islands {'on' if self.show_islands else 'off'}"
                      f"   total {n_isl} face(s)",
                  "a all  1-5 one zone  i islands  g ghost  n/p case  r reset  s shot"]
        self.p.add_text("\n".join(lines), position="upper_left", font_size=9,
                        color="#1b2733", font="courier")
        for zi, nm in enumerate(names):
            self.p.add_text("  ", position=(14, 830 - 22 * zi), font_size=9,
                            color=ZONE_COLOURS[zi % len(ZONE_COLOURS)])
        if self.cam is not None:
            self.p.camera_position = self.cam
        else:
            self.p.view_vector(DEFAULT_VIEW, viewup=DEFAULT_UP)
            self.p.reset_camera()
        self.p.render()

    def run(self):
        self.p.show()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default=THORACIC)
    ap.add_argument("--case", default=None)
    ap.add_argument("--rebuild", action="store_true",
                    help="ignore the cache and recompute every case")
    args = ap.parse_args()

    from anat_bench import find_cases
    cases = find_cases(args.cohort)
    if args.case:
        cases = [c for c in cases if c[0] == args.case] or cases
    if args.rebuild and CACHE.exists():
        for f in CACHE.glob("*.npz"):
            f.unlink()
    print(f"{len(cases)} case(s). First load of each takes a few seconds; "
          f"after that they are cached.")
    Viewer(cases).run()


if __name__ == "__main__":
    main()
