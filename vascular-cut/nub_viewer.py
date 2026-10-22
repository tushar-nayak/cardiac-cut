"""Click through every nub-in geometry and record a verdict on each.

    .venv\\Scripts\\python nub_viewer.py                    # output/nubs
    .venv\\Scripts\\python nub_viewer.py --dir output/nubs2
    .venv\\Scripts\\python nub_viewer.py --flagged          # only ones needing a look

One window, one case at a time, rotate freely. Every nub is ringed and labelled with
what it measured, so a verdict can be given from the picture rather than from a
spreadsheet next to it.

    right / space / n     next case
    left / p              previous case
    g                     good, keep it
    b                     bad, needs work
    u                     undecided (clears a verdict)
    f                     jump to the next case with no verdict
    r                     back to the starting orientation
    s                     save a screenshot next to the mesh
    q                     quit

Every case arrives in the same place at the same size, in whatever orientation the
last one was left in -- so a view worth having is set once and then each geometry is
a keystroke and a glance.

Verdicts are written to review.json in the same folder as they are given, so closing
the window never loses them and a session can be picked up where it stopped. A case
opens showing whatever verdict it already has.

Rims are coloured by what the cohort run measured, not by eye:

    green     the rim lies on its plane and the nub is the length that was asked for
    orange    the nub had to grow to cut cleanly -- longer than k asked for
    red       the rim strays off its plane: torn, not sliced
    grey      an aortic end, not a nub
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import trimesh

MESH_RGB = (0.86, 0.86, 0.88)
GOOD = "#2e9e5b"
GROWN = "#e08a1e"
TORN = "#cc3311"
PLAIN = "#9aa0a6"
VERDICT_COLOUR = {"good": GOOD, "bad": TORN, None: "#cccccc"}

# Where a case first appears. Left anterior oblique-ish, patient upright, so the arch
# and its nubs face the viewer without anything having to be rotated to be judged.
DEFAULT_VIEW = np.array([0.4, -1.0, 0.15])
DEFAULT_UP = np.array([0.0, 0.0, 1.0])


def load_cases(folder, flagged_only=False):
    """Every cut mesh in *folder* that has a sidecar, in case order.

    Reads both shapes of sidecar: the batch's (`<case>_nub_k0.75.json`, written by
    tools/nub_cohort.py) and the desktop tool's own export (`<case>_nub.stl` beside
    `<case>_nub_cut_planes.json`). The second is what a finished cohort actually
    consists of, so a last sweep over the exports has to read it.
    """
    out = []
    for js in sorted(folder.glob("*_nub_k*.json")):
        if js.name == "review.json":
            continue
        stl = js.with_suffix(".stl")
        if not stl.exists():
            continue
        rec = json.loads(js.read_text())
        if flagged_only and not _needs_a_look(rec):
            continue
        out.append((stl, rec))
    for js in sorted(folder.glob("*_cut_planes.json")):
        stl = js.with_name(js.name.replace("_cut_planes.json", ".stl"))
        if not stl.exists():
            continue
        rec = _from_export(js, stl)
        if flagged_only and not _needs_a_look(rec):
            continue
        out.append((stl, rec))
    return out


def _from_export(js, stl):
    """A desktop-tool export, in the shape the viewer draws from.

    The export records what was cut, not why: there is no proud length in it, because
    that is the preset's idea and not the tool's. What it does carry is every plane's
    position, orientation and the diameter of the opening it made, which is what a
    confirmation sweep is looking at.
    """
    d = json.loads(js.read_text())
    vessels = []
    for o in d.get("openings", []):
        moved = (abs(float(o.get("tilt_a_deg") or 0)) > 1e-9
                 or abs(float(o.get("tilt_b_deg") or 0)) > 1e-9
                 or abs(float(o.get("shift_a_mm") or 0)) > 1e-9
                 or abs(float(o.get("shift_b_mm") or 0)) > 1e-9)
        vessels.append({
            "name": o.get("name") or o.get("spline"),
            "ok": True, "cut_ok": True,
            "rim_centroid": o.get("centroid"),
            "plane_normal": o.get("plane_normal"),
            "mouth_ead_mm": o.get("diameter_mm"),
            "offset_mm": o.get("offset_mm"),
            "cut_mode": o.get("cut_mode"),
            "trimmed_group": o.get("trimmed_group"),
            "hand_adjusted": moved,
            "tilt": (float(o.get("tilt_a_deg") or 0), float(o.get("tilt_b_deg") or 0)),
        })
    notes = []
    if d.get("rim_slivers"):
        notes.append(f"{d['rim_slivers']} sliver triangle(s) on the rims")
    q = d.get("min_triangle_quality")
    if q is not None and q < 0.01:
        notes.append(f"worst triangle quality {q:.4f}")
    return {"case": stl.stem.replace("_nub", ""), "k": None,
            "openings": d.get("n_openings"), "euler": d.get("euler_number"),
            "faces_out": d.get("n_faces"), "area_mm2": d.get("area_mm2"),
            "export": True, "notes": notes, "vessels": vessels}


def _needs_a_look(rec):
    vs = rec.get("vessels", [])
    if rec.get("export"):
        return any(v.get("hand_adjusted") for v in vs) or bool(rec.get("notes"))
    return (any(v.get("rescued_mm") or v.get("retried_mm") or v.get("clamped")
                or v.get("marginal") for v in vs)
            or any(v["ok"] and not v.get("cut_ok") for v in vs)
            or any(not v["ok"] for v in vs)
            or bool(rec.get("notes")) or bool(rec.get("validate")))


def _rim_colour(v):
    if v.get("hand_adjusted"):
        return GROWN            # a plane that was tilted or shifted off default
    if v.get("ok") and not v.get("cut_ok"):
        return TORN
    if v.get("marginal"):
        return TORN            # placeable, and too short to be worth keeping
    if v.get("rescued_mm") or v.get("retried_mm") or v.get("clamped"):
        return GROWN
    return GOOD


class Review:
    """The verdict file. Written on every change, because a crash three hours into a
    review of sixty geometries must not cost the review."""

    def __init__(self, path):
        self.path = path
        self.data = {}
        if path.exists():
            try:
                self.data = json.loads(path.read_text())
            except Exception:                                      # noqa: BLE001
                self.data = {}

    def get(self, key):
        return (self.data.get(key) or {}).get("verdict")

    def set(self, key, verdict, note=None):
        if verdict is None:
            self.data.pop(key, None)
        else:
            self.data[key] = {"verdict": verdict, "note": note}
        self.path.write_text(json.dumps(self.data, indent=2, sort_keys=True))

    def counts(self):
        good = sum(1 for v in self.data.values() if v.get("verdict") == "good")
        bad = sum(1 for v in self.data.values() if v.get("verdict") == "bad")
        return good, bad


class Viewer:
    def __init__(self, cases, review, folder, off_screen=False):
        self.cases = cases
        self.review = review
        self.folder = folder
        self.i = 0
        self.off_screen = off_screen
        self.p = pv.Plotter(window_size=(1400, 900), off_screen=off_screen)
        self.p.set_background("white")
        if off_screen:
            return
        self.p.add_key_event("Right", lambda: self.step(+1))
        self.p.add_key_event("space", lambda: self.step(+1))
        self.p.add_key_event("n", lambda: self.step(+1))
        self.p.add_key_event("Left", lambda: self.step(-1))
        self.p.add_key_event("p", lambda: self.step(-1))
        self.p.add_key_event("g", lambda: self.verdict("good"))
        self.p.add_key_event("b", lambda: self.verdict("bad"))
        self.p.add_key_event("u", lambda: self.verdict(None))
        self.p.add_key_event("f", self.next_unjudged)
        self.p.add_key_event("r", self.home)
        self.p.add_key_event("s", self.shot)

    def home(self):
        """Back to the starting orientation, for when a rotation has got away."""
        self.draw(reset_camera=True, home=True)

    def next_unjudged(self):
        """Skip to the next case with no verdict, so a resumed review does not have
        to be walked back through from the start."""
        for d in range(1, len(self.cases) + 1):
            j = (self.i + d) % len(self.cases)
            if self.review.get(self.cases[j][1]["case"]) is None:
                self.i = j
                self.draw(reset_camera=True)
                return
        print("  every case has a verdict")

    # -- navigation ---------------------------------------------------------
    def key(self):
        return self.cases[self.i][1]["case"]

    def step(self, d):
        self.i = (self.i + d) % len(self.cases)
        self.draw(reset_camera=True)

    def verdict(self, v):
        self.review.set(self.key(), v)
        # Straight on to the next one: the whole point is to get through sixty of
        # these, and a keystroke that leaves you looking at what you just judged
        # costs a second every time.
        if v is not None:
            self.step(+1)
        else:
            self.draw(reset_camera=False)

    def shot(self):
        stl, _ = self.cases[self.i]
        out = stl.with_name(stl.stem + "_view.png")
        self.p.screenshot(str(out))
        print(f"  saved {out}")

    # -- drawing ------------------------------------------------------------
    def look_direction(self):
        """The direction the camera is currently looking from, and its up vector.

        Carried from case to case so a rotation you settled on survives the move.
        """
        try:
            pos, focal, up = self.p.camera_position
        except Exception:                                          # noqa: BLE001
            return DEFAULT_VIEW, DEFAULT_UP
        d = np.asarray(pos, float) - np.asarray(focal, float)
        n = float(np.linalg.norm(d))
        if not np.isfinite(n) or n < 1e-9:
            return DEFAULT_VIEW, DEFAULT_UP
        return d / n, np.asarray(up, float)

    def draw(self, reset_camera=True, home=False):
        stl, rec = self.cases[self.i]
        # Where we are looking FROM is kept; where we are looking AT is not. Every
        # case sits at its own scanner coordinates, tens of centimetres apart, so
        # carrying the whole camera over pointed it at empty space and the next
        # aorta came up blank -- it was there, just off screen. Keep the direction,
        # re-frame the distance, and each case arrives in the same place at the same
        # size, in whatever orientation the last one was left in.
        look, up = (DEFAULT_VIEW, DEFAULT_UP) if home else self.look_direction()
        cam = None if reset_camera else self.p.camera_position
        # clear_actors, NOT clear: clear() takes the lights out with the actors, and
        # every case after the first then renders as a flat silhouette with no
        # shading at all -- which for a review whose whole job is judging surfaces by
        # eye is worse than useless, and looks like a broken mesh rather than a
        # broken light.
        self.p.clear_actors()

        mesh = trimesh.load(str(stl), process=False)
        poly = pv.PolyData(np.asarray(mesh.vertices),
                           np.hstack([np.full((len(mesh.faces), 1), 3),
                                      np.asarray(mesh.faces)]).ravel())
        self.p.add_mesh(poly, color=MESH_RGB, smooth_shading=True,
                        specular=0.2, opacity=1.0)

        # rims, coloured by what was measured, plus a label at each nub's mouth
        for v in rec.get("vessels", []):
            c = v.get("rim_centroid")
            if not c or not v.get("ok"):
                continue
            colour = _rim_colour(v)
            self.p.add_mesh(pv.Sphere(radius=1.2, center=c), color=colour)
            # Lifted off the rim along the cut normal. Left on the rim itself the
            # labels of two arch vessels a centimetre apart print on top of each
            # other and neither can be read.
            nrm = np.asarray(v.get("plane_normal") or (0.0, 0.0, 1.0), float)
            anchor = np.asarray(c, float) + nrm * 7.0
            bits = [v["name"]]
            if v.get("proud_mm") is not None:
                bits.append(f"{v['proud_mm']:.1f} mm proud")
            if v.get("mouth_ead_mm"):
                bits.append(f"D {v['mouth_ead_mm']:.1f}")
            if v.get("offset_mm") is not None and v.get("proud_mm") is None:
                bits.append(f"cut {v['offset_mm']:.1f} mm ({v.get('cut_mode')})")
            if v.get("trimmed_group"):
                bits.append(f"trims {len(v['trimmed_group'])}")
            if v.get("hand_adjusted"):
                ta, tb = v.get("tilt", (0, 0))
                bits.append(f"tilted {ta:+.0f}/{tb:+.0f}" if (ta or tb) else "shifted")
            if v.get("rescued_mm"):
                bits.append(f"grown +{v['rescued_mm']:.1f}")
            if v.get("retried_mm"):
                bits.append(f"retried +{v['retried_mm']:.1f}")
            if v.get("clamped"):
                bits.append("SHORTENED, vessel too short")
            if v.get("marginal"):
                bits.append("BARELY A NUB")
            if v.get("flat_mm") is not None and not v.get("cut_ok"):
                bits.append(f"TORN {v['flat_mm']:.1f} mm")
            self.p.add_mesh(pv.Line(c, anchor), color=colour, line_width=1)
            self.p.add_point_labels(
                [anchor], ["  " + "  |  ".join(bits)], font_size=13,
                text_color=colour, shape_opacity=0.0, show_points=False,
                always_visible=True)

        self.p.add_text(self._header(rec), position="upper_left", font_size=11,
                        color="black")
        self.p.add_text(self._footer(rec), position="lower_left", font_size=9,
                        color="#555555")
        if cam is not None:
            self.p.camera_position = cam
        else:
            self.p.view_vector(look, viewup=up)
            self.p.reset_camera()
        self.p.render()

    def _header(self, rec):
        good, bad = self.review.counts()
        v = self.review.get(rec["case"])
        mark = {"good": "KEEP", "bad": "NEEDS WORK", None: "-- not judged --"}[v]
        # A finished export carries no k -- that is the preset's idea, not the tool's.
        # Say what the geometry IS instead, which is what a confirmation sweep wants.
        what = (f"k={rec['k']:g}" if rec.get("k") is not None
                else f"{rec.get('openings')} openings, euler {rec.get('euler')}")
        return (f"[{self.i + 1}/{len(self.cases)}]  {rec['case']}   {what}\n"
                f"{mark}          reviewed: {good} keep, {bad} needs work")

    def _footer(self, rec):
        lines = []
        for v in rec.get("vessels", []):
            if not v.get("ok"):
                lines.append(f"NOT CUT  {v['name']}: {v.get('why', '')}")
        for n in rec.get("notes", []):
            lines.append(n)
        lines.append("g keep   b needs work   u clear   f next unjudged   "
                     "<- ->  move   r reset view   s screenshot   q quit")
        return "\n".join(lines[-6:])

    def run(self):
        self.draw()
        self.p.show()

    def render_all(self):
        """Every case to a PNG, drawn by the same code the window uses.

        A contact sheet to flip through, and the only way to check what the window
        draws without opening one -- rims, labels and all, from the same call.
        """
        for i in range(len(self.cases)):
            self.i = i
            self.draw(reset_camera=True, home=True)
            # No zoom past 1: the nubs are the smallest things on screen and they sit
            # at the top of the arch, which is exactly what a zoom crops off.
            self.p.camera.zoom(0.9)
            out = self.cases[i][0].with_name(self.cases[i][0].stem + "_view.png")
            self.p.screenshot(str(out))
            print(f"  {out.name}")
        self.p.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="output/nubs")
    ap.add_argument("--flagged", action="store_true",
                    help="only cases with a grown, torn or unplaced nub")
    ap.add_argument("--render", action="store_true",
                    help="write a PNG per case instead of opening the window")
    args = ap.parse_args()

    folder = Path(args.dir)
    if not folder.exists():
        raise SystemExit(f"no such folder: {folder}  "
                         f"(run tools/nub_cohort.py run first)")
    cases = load_cases(folder, args.flagged)
    if not cases:
        raise SystemExit(f"nothing to review in {folder}")

    review = Review(folder / "review.json")
    print(f"{len(cases)} geometries in {folder}")
    if args.render:
        Viewer(cases, review, folder, off_screen=True).render_all()
        return 0
    Viewer(cases, review, folder).run()

    good, bad = review.counts()
    print(f"\nreviewed: {good} keep, {bad} needs work, "
          f"{len(cases) - good - bad} not judged")
    print(f"verdicts in {folder / 'review.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
