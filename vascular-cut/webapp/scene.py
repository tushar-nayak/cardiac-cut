# The web app's 3-D scene, built to match the desktop window's.
#
# Both render the same thing with the same library: a PyVista plotter holding the
# cut mesh, the centerline and one disc per cut plane. The desktop window embeds
# that plotter in Qt; here the same kind of render window goes to trame, which
# renders on the server and streams pixels to the browser. The desktop tool does not
# import this module -- it has to stand alone in the distribution -- so the colours
# and the update strategy are copied deliberately rather than shared.
#
# Actors are named and updated in place: rebuilding the whole scene on every slider
# tick is what made the desktop version stutter before it drew planes separately
# from the mesh.
#
# Three contents, one window: the cut view (mesh, centerline, cut discs), the
# cardiac view (inner and outer walls, the removed cap, the plane) and the shape
# view (the CAMP partition, coloured by patch). Switching sections takes the others'
# actors down; see clear_cut, clear_heart and clear_patches.

import numpy as np
import pyvista as pv


def to_polydata(mesh):
    """trimesh -> PolyData, without a copy of the face array per call."""
    tri = np.asarray(mesh.faces, dtype=np.int64)
    faces = np.hstack([np.full((len(tri), 1), 3, dtype=np.int64), tri]).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


class Scene:
    """Owns a plotter and the named actors on it.

    Two update paths on purpose, mirroring the desktop tool:
      set_mesh()   - expensive, after a re-cut
      set_planes() - cheap, every drag frame
    """

    def __init__(self, plotter=None, off_screen=True):
        self.plotter = plotter or pv.Plotter(off_screen=off_screen)
        # The same anatomical triad the desktop window carries, for the same reason:
        # every Left/Right name follows from the coordinate frame, and an anterior
        # view puts the patient's LEFT on the viewer's RIGHT. L/P/S are the +x/+y/+z
        # of branch_cut.ANATOMY_FRAME. Parity matters here -- a person moving between
        # the two front ends must not lose the one cue that lets them check a side.
        try:
            self.plotter.add_axes(xlabel="L", ylabel="P", zlabel="S",
                                  color="#333333")
        except Exception:                                       # noqa: BLE001
            pass                 # never worth failing a scene over
        self._actors = {}
        # Taken from the desktop window rather than chosen here. Someone moving
        # between the two should not have to re-learn what a colour means.
        self.style = dict(
            background="white",
            mesh="#d8b7a0",
            mesh_opacity=0.55,       # translucent, so the centerline reads through
            centerline="#1f2d5a",    # dark navy
            plane="#e69f00",         # amber, translucent
            plane_active="#d55e00",  # the selected row's disc: darker orange
            measure_plane="#8e24aa",  # purple: a measurement, not a cut
            mark="#0a7d2c",          # pinned points, as on the desktop
            plane_opacity=0.45,
            loop_z="#0070c0",        # boundary loops, coloured by their cut's mode
            loop_normal="#00b050",
            patch_trimmed="#bdbdbd",  # patches the trim excluded, still drawn
            patch_edge="#3a3a3a",     # the partition, when outlines are on
        )
        # A soft vertical wash, lavender-white at the top to lavender-grey at the
        # bottom, so the model sits in a lit space rather than on a blank page. Every
        # channel stays at or above 234: the privacy checks count anything darker than
        # 231 as model ink, and an empty scene must read as empty.
        self.style["background"] = "#EAECF8"
        self.style["background_top"] = "#FBFBFF"
        self.plotter.set_background(self.style["background"],
                                    top=self.style["background_top"])

    # ── actor bookkeeping ────────────────────────────────────────────────────
    def _replace(self, key, adder):
        old = self._actors.pop(key, None)
        if old is not None:
            try:
                self.plotter.remove_actor(old, render=False)
            except Exception:
                pass
        actor = adder()
        if actor is not None:
            self._actors[key] = actor
        return actor

    def clear(self):
        for key in list(self._actors):
            self._replace(key, lambda: None)
        self._drop_scalar_bar()

    def clear_cut(self):
        """Everything the cut view puts on screen, and nothing the shape view does.

        The two views share one render window, so switching between them has to take
        the other one's actors down. The camera is deliberately left alone: it is the
        same geometry either way, and resetting it would make a change of colouring
        look like a change of subject.
        """
        # The measure plane and the pinned marks are the cut view's too. Left out, they
        # were drawn over the delta K patches.
        for key in [k for k in self._actors
                    if k in ("mesh", "centerline", "measure")
                    or k.startswith(("plane:", "loop:", "mark:"))]:
            self._replace(key, lambda: None)

    # ── content ──────────────────────────────────────────────────────────────
    def set_mesh(self, mesh, valid=True, face_rgb=None):
        """The cut surface. Translucent, so the centerline and the cut discs inside
        it stay visible -- the desktop window does the same, and it is the only way
        to see that a plane sits where you think it does."""
        if mesh is None:
            self._replace("mesh", lambda: None)
            return
        poly = to_polydata(mesh)
        colour = self.style["mesh"] if valid else "#8a8a8a"
        # Per-face colour when the caller supplies it: the surface is tinted by which
        # branch owns each face, so a branch's wall, its disc and its row all read as
        # the same colour. An invalid cut stays grey -- the point of grey is that
        # something is wrong, and colouring over it would hide that.
        if face_rgb is not None and valid and len(face_rgb) == len(mesh.faces):
            self._replace("mesh", lambda: self.plotter.add_mesh(
                poly, scalars=np.asarray(face_rgb, np.uint8), rgb=True,
                preference="cell", smooth_shading=True,
                opacity=self.style["mesh_opacity"], specular=0.25,
                specular_power=15, name="mesh", render=False, show_scalar_bar=False))
            return
        self._replace("mesh", lambda: self.plotter.add_mesh(
            poly, color=colour, smooth_shading=True, opacity=self.style["mesh_opacity"],
            specular=0.25, specular_power=15, name="mesh", render=False))

    def set_loops(self, mesh, modes=None):
        """A sphere cloud on each open boundary, coloured by its cut's mode.

        This is what makes a bad cut obvious: an opening that should not be there
        shows up as a ring of spheres somewhere unexpected, long before you think
        to read the euler number.
        """
        for key in [k for k in self._actors if k.startswith("loop:")]:
            self._replace(key, lambda: None)
        if mesh is None:
            return
        try:
            from branch_cut import boundary_loops
            loops = boundary_loops(mesh)
        except Exception:
            return
        verts = np.asarray(mesh.vertices, dtype=float)
        for i, loop in enumerate(loops):
            pts = verts[np.asarray(loop, dtype=np.int64)]
            if len(pts) == 0:
                continue
            mode = (modes or {}).get(i, "z")
            colour = (self.style["loop_z"] if mode == "z"
                      else self.style["loop_normal"])
            cloud = pv.PolyData(pts)
            self._replace(
                f"loop:{i}",
                lambda c=cloud, col=colour: self.plotter.add_mesh(
                    c, color=col, point_size=7, render_points_as_spheres=True,
                    pickable=False, render=False))

    def set_centerline(self, tree, show=True):
        if tree is None or not show:
            self._replace("centerline", lambda: None)
            return
        lines = []
        for seg in tree.segments.values():
            pts = np.asarray(seg.pts, dtype=float)
            if len(pts) > 1:
                lines.append(pv.lines_from_points(pts))
        if not lines:
            self._replace("centerline", lambda: None)
            return
        merged = lines[0] if len(lines) == 1 else lines[0].merge(lines[1:])
        # Not pickable, like every actor except the mesh -- otherwise a hover picker
        # latches onto the centerline instead of the wall behind it.
        self._replace("centerline", lambda: self.plotter.add_mesh(
            merged, color=self.style["centerline"], line_width=3,
            render_lines_as_tubes=True, pickable=False, render=False))

    def set_planes(self, planes):
        """One disc per cut. `planes` is a list of dicts with keys
        origin, normal, radius, active."""
        for key in [k for k in self._actors if k.startswith("plane:")]:
            self._replace(key, lambda: None)
        for i, p in enumerate(planes or []):
            disc = pv.Disc(center=p["origin"], normal=p["normal"],
                           inner=0.0, outer=float(p["radius"]), c_res=48)
            colour = (self.style["plane_active"] if p.get("active")
                      else self.style["plane"])
            self._replace(
                f"plane:{i}",
                lambda d=disc, c=colour: self.plotter.add_mesh(
                    d, color=c, opacity=self.style["plane_opacity"],
                    smooth_shading=True, pickable=False, render=False))

    def set_marks(self, points, radius=1.6):
        """Pinned points on the surface. They are data, not part of the cut, so they
        survive a re-cut and are drawn from this list every time."""
        for key in [k for k in self._actors if k.startswith("mark:")]:
            self._replace(key, lambda: None)
        for i, p in enumerate(points or []):
            sph = pv.Sphere(radius=radius, center=[float(v) for v in p])
            self._replace(f"mark:{i}", lambda s=sph: self.plotter.add_mesh(
                s, color=self.style["mark"], pickable=False, render=False))

    def set_measure_plane(self, plane):
        """The station a measurement was taken at, or None to take it away.

        Its own colour, and sized from the section rather than from a cut's reach,
        so it is never read as a cut plane: this one removes nothing.
        """
        if not plane:
            self._replace("measure", lambda: None)
            return
        disc = pv.Disc(center=plane["origin"], normal=plane["normal"],
                       inner=0.0, outer=float(plane["radius"]), c_res=48)
        self._replace("measure", lambda: self.plotter.add_mesh(
            disc, color=self.style["measure_plane"], opacity=0.55,
            smooth_shading=True, pickable=False, render=False))

    # ── the cardiac view ─────────────────────────────────────────────────────
    # Inner wall blue and outer wall the cut view's tan, so the two read apart at a
    # glance; what the plane removes is grey and faint, never gone, so you can see
    # what you cut. Holes are drawn as rings in the hole colour.
    HEART = dict(inner="#4c78a8", outer="#d8b7a0", removed="#9e9e9e",
                 hole="#d62728", rim="#0070c0")

    def clear_heart(self):
        for key in [k for k in self._actors if k.startswith("heart:")]:
            self._replace(key, lambda: None)

    def set_heart(self, inner=None, outer=None, removed=None, show=None,
                  outer_opacity=0.55):
        """The two walls and the removed cap. `show` maps inner/outer/removed to
        bool; a wall not shown is taken down rather than hidden, so its loops go
        with it."""
        show = {"inner": True, "outer": True, "removed": True, **(show or {})}
        for key in [k for k in self._actors if k.startswith("heart:")
                    and k != "heart:plane"]:
            self._replace(key, lambda: None)
        for name, mesh, colour, opacity in (
                ("inner", inner, self.HEART["inner"], 1.0),
                ("outer", outer, self.HEART["outer"], outer_opacity),
                ("removed", removed, self.HEART["removed"], 0.18)):
            if mesh is None or not show.get(name) or len(mesh.faces) == 0:
                continue
            poly = to_polydata(mesh)
            self._replace(f"heart:{name}", lambda p=poly, c=colour, o=opacity:
                          self.plotter.add_mesh(p, color=c, opacity=o,
                                                smooth_shading=True, specular=0.2,
                                                specular_power=12, pickable=False,
                                                render=False))
            if name == "removed":
                continue
            # Open boundaries: the longest is the basal cut, the rest are holes.
            try:
                from branch_cut import boundary_loops
                loops = boundary_loops(mesh)
            except Exception:
                loops = []
            verts = np.asarray(mesh.vertices, dtype=float)
            sized = []
            for loop in loops:
                pts = verts[np.asarray(loop, dtype=np.int64)]
                if len(pts):
                    per = float(np.linalg.norm(np.diff(np.vstack([pts, pts[:1]]), axis=0),
                                               axis=1).sum())
                    sized.append((per, pts))
            sized.sort(key=lambda t: -t[0])
            for i, (_, pts) in enumerate(sized):
                colour_l = self.HEART["rim"] if i == 0 else self.HEART["hole"]
                cloud = pv.PolyData(pts)
                self._replace(f"heart:loop:{name}:{i}",
                              lambda c=cloud, col=colour_l: self.plotter.add_mesh(
                                  c, color=col, point_size=6,
                                  render_points_as_spheres=True, pickable=False,
                                  render=False))

    def set_heart_plane(self, plane):
        """The cut plane as a disc, or None to take it down. Cheap: this is what
        follows a slider drag, while the walls re-cut only on release."""
        if not plane:
            self._replace("heart:plane", lambda: None)
            return
        disc = pv.Disc(center=plane["origin"], normal=plane["normal"],
                       inner=0.0, outer=float(plane["radius"]), c_res=64)
        self._replace("heart:plane", lambda: self.plotter.add_mesh(
            disc, color=self.style["plane_active"], opacity=self.style["plane_opacity"],
            smooth_shading=True, pickable=False, render=False))

    # ── the CAMP partition ───────────────────────────────────────────────────
    def _drop_scalar_bar(self):
        try:
            for title in list(self.plotter.scalar_bars.keys()):
                self.plotter.remove_scalar_bar(title, render=False)
        except Exception:
            pass

    def clear_patches(self):
        for key in [k for k in self._actors if k.startswith("patch")]:
            self._replace(key, lambda: None)
        self._drop_scalar_bar()

    def set_patches(self, verts, faces, scalars=None, clim=None, cmap="plasma",
                    n_colors=256, title="", show_bar=True, kept_mask=None,
                    edges=None):
        """The partitioned surface, one colour per patch.

        Opaque, unlike the cut view. There is nothing behind it worth seeing here, and
        a translucent surface blends the colour of every patch with the ones behind
        it, which is how you end up reading a value off the picture that no patch has.

        Trimmed patches are drawn separately in grey rather than being dropped. They
        are still part of the surface, and a hole where one used to be would read as a
        defect in the mesh instead of as an excluded outlier.
        """
        self.clear_patches()
        if verts is None or faces is None:
            return

        verts = np.asarray(verts, dtype=float)
        faces = np.asarray(faces, dtype=np.int64)

        def poly_of(sub):
            return pv.PolyData(
                verts, np.hstack([np.full((len(sub), 1), 3, dtype=np.int64),
                                  sub]).ravel())

        keep = None if kept_mask is None else np.asarray(kept_mask, dtype=bool)
        if keep is not None and not keep.all():
            trimmed = poly_of(faces[~keep])
            self._replace("patch:trim", lambda: self.plotter.add_mesh(
                trimmed, color=self.style["patch_trimmed"], smooth_shading=True,
                pickable=False, render=False))
            faces_main = faces[keep]
            scal = None if scalars is None else np.asarray(scalars)[keep]
        else:
            faces_main, scal = faces, (None if scalars is None
                                       else np.asarray(scalars))

        main = poly_of(faces_main)
        if scal is not None:
            main.cell_data["v"] = np.asarray(scal, dtype=float)

        bar_args = dict(title=title or " ", color="black", n_labels=5,
                        vertical=True, title_font_size=15, label_font_size=12,
                        width=0.045, height=0.46, position_x=0.925,
                        position_y=0.28, fmt="%.3g")
        self._replace("patch:main", lambda: self.plotter.add_mesh(
            main, scalars=("v" if scal is not None else None),
            # Before a run there is nothing to colour by, and the surface still has to
            # be on screen -- an empty viewport reads as a failure to load.
            color=(None if scal is not None else self.style["mesh"]),
            cmap=cmap, clim=clim, n_colors=int(n_colors),
            show_scalar_bar=bool(show_bar and scal is not None),
            scalar_bar_args=bar_args, smooth_shading=True,
            specular=0.15, specular_power=12, pickable=False, render=False))

        if edges is not None and len(edges):
            segs = np.asarray(edges, dtype=np.int64)
            lines = np.hstack([np.full((len(segs), 1), 2, dtype=np.int64),
                               segs]).ravel()
            outline = pv.PolyData(verts, lines=lines)
            self._replace("patch:edges", lambda: self.plotter.add_mesh(
                outline, color=self.style["patch_edge"], line_width=1.2,
                pickable=False, render=False))

    # ── camera ───────────────────────────────────────────────────────────────
    def reset_camera(self, mesh=None):
        self.plotter.reset_camera(render=False)
        if mesh is not None:
            self.plotter.view_isometric(render=False)

    def render(self):
        try:
            self.plotter.render()
        except Exception:
            pass


def planes_for(session, active_label=None):
    """Cut-plane discs for every opening in a session.

    Two details that are easy to get wrong and both matter:

    The disc sits at cut_ANCHOR, not cut_plane. The two differ by the in-plane
    shift, and the anchor is what the band selection actually uses -- drawing at
    cut_plane means the Shift sliders move the cut without the disc following,
    which is worse than not drawing it at all.

    The radius is band_width, because the disc IS the cut's reach. Sizing it from
    the measured opening diameter would draw a plausible-looking disc that does not
    show the thing you need to see: how far the cut grabs, which is exactly what you
    are adjusting when one cut punches two holes.
    """
    from branch_cut import band_width

    out = []
    active, offs = session._cut_list()
    for t in active:
        try:
            origin, normal = t.cut_anchor(float(offs[t.label]))
        except Exception:
            continue
        out.append({
            "label": t.label,
            "origin": np.asarray(origin, dtype=float),
            "normal": np.asarray(normal, dtype=float),
            "radius": max(float(band_width(t)), 2.0),
            "active": (t.label == active_label),
        })
    return out
