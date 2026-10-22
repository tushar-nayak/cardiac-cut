"""Rendering of CAMP results, with an adjustable integrated-K color axis.

Outputs (all optional, individually guarded so a rendering failure never loses your
CSV/NPZ data):

* ``build_patch_figure`` / ``save_patch_html`` — a self-contained interactive Plotly
  view.  Baked-in controls (no server needed):
      - color-axis sliders  (drag the min / max of the integrated-K scale)
      - colorscale dropdown (Plasma / Viridis / Turbo / ...)
      - a button to toggle between the curvature map and random per-patch colors
* ``save_patch_png`` — a static PyVista render (needs a GL context; auto-skipped).
* ``view`` — convenience: load a result ``.npz`` (or an output folder) and show it,
  colored by any per-patch metric with an optional fixed color axis.
"""
import numpy as np

# per-patch metric -> (npz key, colorbar label)
METRIC = {
    "integrated_K": ("K_integrated_patch", "∫K dA"),
    "avg_kg":       ("avg_kg",             "K (patch center)"),
    "eta":          ("eta_patch",          "η"),
    "R2":           ("R2_patch",           "R²"),
}

COLORSCALES = ["Plasma", "Viridis", "Turbo", "Inferno", "Magma", "Cividis", "RdBu", "Jet"]


def _percentile_clip(values, lo=2, hi=98):
    v = np.asarray(values, float)
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        return -1.0, 1.0
    a, b = np.percentile(finite, [lo, hi])
    if a == b:
        a, b = float(finite.min()), float(finite.max())
    if a == b:
        a, b = a - 1.0, b + 1.0
    return float(a), float(b)


def _random_palette(n, seed=11):
    import matplotlib.colors as mcolors
    g = np.random.default_rng(seed)
    hsv = np.column_stack([g.uniform(0, 1, n), g.uniform(0.5, 0.85, n), g.uniform(0.78, 1, n)])
    rgb = (mcolors.hsv_to_rgb(hsv) * 255).astype(int)
    return ["rgb(%d,%d,%d)" % tuple(c) for c in rgb]


def _axis_controls(face_val, cmin, cmax, colorscale):
    """Build Plotly sliders (color min / max) + a colorscale dropdown for trace 0.

    The sliders span the full data range; their initial (active) positions are the
    nearest grid ticks to the requested cmin/cmax.  The exact cmin/cmax the caller
    asked for stay on the trace itself (set in ``build_patch_figure``), so an
    explicit axis is honored precisely and only snaps to the grid once dragged.
    """
    finite = np.asarray(face_val, float)
    finite = finite[np.isfinite(finite)]
    lo_full = float(finite.min()) if finite.size else -1.0
    hi_full = float(finite.max()) if finite.size else 1.0
    # widen the grid so an explicit axis set outside the data range is still reachable
    lo_full = min(lo_full, float(cmin)); hi_full = max(hi_full, float(cmax))
    if lo_full == hi_full:
        lo_full, hi_full = lo_full - 1.0, hi_full + 1.0
    grid = np.linspace(lo_full, hi_full, 25)

    def nearest(v):
        return int(np.argmin(np.abs(grid - v)))

    cmin0, cmax0 = nearest(cmin), nearest(cmax)
    cmin_steps = [dict(method="restyle", label=f"{g:.3g}", args=[{"cmin": [float(g)]}, [0]])
                  for g in grid]
    cmax_steps = [dict(method="restyle", label=f"{g:.3g}", args=[{"cmax": [float(g)]}, [0]])
                  for g in grid]
    sliders = [
        dict(active=cmin0, x=0.0, y=0, len=0.46, pad={"t": 40, "b": 10},
             currentvalue={"prefix": "color min = ", "font": {"size": 12}}, steps=cmin_steps),
        dict(active=cmax0, x=0.54, y=0, len=0.46, pad={"t": 40, "b": 10},
             currentvalue={"prefix": "color max = ", "font": {"size": 12}}, steps=cmax_steps),
    ]
    updatemenus = [
        dict(type="buttons", direction="right", x=0.0, xanchor="left", y=1.06, yanchor="bottom",
             showactive=True, buttons=[
                 dict(label="integrated K", method="restyle", args=[{"visible": [True, False]}]),
                 dict(label="patches", method="restyle", args=[{"visible": [False, True]}]),
             ]),
        dict(type="dropdown", direction="down", x=1.0, xanchor="right", y=1.06, yanchor="bottom",
             showactive=True, buttons=[
                 dict(label=cs, method="restyle", args=[{"colorscale": cs}, [0]])
                 for cs in ([colorscale] + [c for c in COLORSCALES if c != colorscale])
             ]),
    ]
    return sliders, updatemenus


def build_patch_figure(V, F, face_labels, values, value_label="∫K dA",
                       cmin=None, cmax=None, colorscale="Plasma", clip=(2, 98),
                       title="", controls=True):
    """Interactive Plotly figure of a CAMP partition, colored per face by ``values``.

    ``values`` is one value PER PATCH (indexed by ``face_labels``); pass any per-patch
    metric.  ``cmin``/``cmax`` fix the color axis (default: 2–98th percentile).  When
    ``controls`` is True the figure carries color-axis sliders, a colorscale dropdown,
    and an integrated-K/patches toggle.  Call ``fig.show()`` to display inline.
    """
    import plotly.graph_objects as go
    V = np.asarray(V, float); F = np.asarray(F, int)
    fl = np.asarray(face_labels, int)
    face_val = np.asarray(values, float)[fl]                 # per-face value
    auto_lo, auto_hi = _percentile_clip(face_val, *clip)
    lo = auto_lo if cmin is None else float(cmin)
    hi = auto_hi if cmax is None else float(cmax)

    sliders, updatemenus = [], []
    if controls:
        # sliders start at the nearest tick; the exact lo/hi stay on the trace below
        sliders, updatemenus = _axis_controls(face_val, lo, hi, colorscale)

    i, j, k = F[:, 0], F[:, 1], F[:, 2]
    curv = go.Mesh3d(
        x=V[:, 0], y=V[:, 1], z=V[:, 2], i=i, j=j, k=k,
        intensity=face_val, intensitymode="cell",
        colorscale=colorscale, cmin=lo, cmax=hi,
        colorbar=dict(title=value_label), showscale=True, name=value_label,
        flatshading=True, visible=True,
    )
    pal = np.array(_random_palette(int(fl.max()) + 1))
    patches = go.Mesh3d(
        x=V[:, 0], y=V[:, 1], z=V[:, 2], i=i, j=j, k=k,
        facecolor=pal[fl], showscale=False, name="patches",
        flatshading=True, visible=False,
    )
    fig = go.Figure([curv, patches])
    fig.update_layout(
        title=title,
        scene=dict(aspectmode="data",
                   xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
        margin=dict(l=0, r=0, t=60, b=80 if controls else 10),
        sliders=sliders, updatemenus=updatemenus,
    )
    return fig


def save_patch_html(path, V, F, face_labels, values, value_label="∫K dA",
                    cmin=None, cmax=None, colorscale="Plasma", clip=(2, 98), title=""):
    """Write a standalone interactive HTML viewer.  Returns path, or None on failure."""
    try:
        fig = build_patch_figure(V, F, face_labels, values, value_label=value_label,
                                 cmin=cmin, cmax=cmax, colorscale=colorscale, clip=clip,
                                 title=title, controls=True)
        # write UTF-8 explicitly: fig.write_html uses the platform codec (cp1252 on
        # Windows), which chokes on non-ASCII glyphs like the integral sign.
        html = fig.to_html(include_plotlyjs="inline", full_html=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(html)
        return str(path)
    except Exception as e:                                    # pragma: no cover
        print(f"      [html] skipped ({e})")
        return None


def save_patch_png(path, V, F, face_labels, values, cmin=None, cmax=None,
                   colorscale="plasma", clip=(2, 98)):
    """Static PyVista render of the color map.  Returns path or None."""
    try:
        import pyvista as pv
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
    except Exception as e:                                    # pragma: no cover
        print(f"      [png] skipped (pyvista/matplotlib unavailable: {e})")
        return None
    try:
        pv.OFF_SCREEN = True
        V = np.asarray(V, float); F = np.asarray(F, int)
        fl = np.asarray(face_labels, int)
        face_val = np.asarray(values, float)[fl]
        auto_lo, auto_hi = _percentile_clip(face_val, *clip)
        lo = auto_lo if cmin is None else float(cmin)
        hi = auto_hi if cmax is None else float(cmax)
        rgb = (cm.get_cmap(colorscale)(mcolors.Normalize(lo, hi)(face_val))[:, :3] * 255).astype(np.uint8)

        faces_pv = np.hstack([np.full((len(F), 1), 3, np.int64), F.astype(np.int64)]).ravel()
        m = pv.PolyData(V, faces_pv)
        m.cell_data["rgb"] = rgb
        pl = pv.Plotter(off_screen=True, window_size=(1100, 1100), border=False)
        pl.add_mesh(m, scalars="rgb", rgb=True)
        pl.view_vector([1, -0.2, 0.3]); pl.camera.up = (0, 0, 1); pl.camera.zoom(1.4)
        pl.screenshot(str(path))
        pl.close()
        return str(path)
    except Exception as e:                                    # pragma: no cover
        print(f"      [png] skipped ({e})")
        return None


def _resolve_npz(source):
    from pathlib import Path
    p = Path(source)
    if p.is_dir():
        hits = sorted(p.glob("*_camp.npz")) or sorted(p.glob("**/*_camp.npz"))
        if not hits:
            raise FileNotFoundError(f"No *_camp.npz found under {p}")
        return hits[0]
    if p.suffix != ".npz":
        raise ValueError(f"Expected a .npz file or a folder, got {p}")
    return p


def view(source, metric="integrated_K", cmin=None, cmax=None, colorscale="Plasma",
         clip=(2, 98), title=None, controls=True):
    """Load a CAMP result and return an interactive figure.

    ``source`` is a ``*_camp.npz`` file or an output folder containing one.
    ``metric`` is one of ``integrated_K`` (default), ``avg_kg``, ``eta``, ``R2``.

    Example
    -------
    >>> from camp_toolkit import view
    >>> view("results/KY160_1_open_1").show()                 # adjustable color axis
    >>> view("results/KY160_1_open_1", cmin=-5, cmax=40).show()  # fixed axis
    """
    from pathlib import Path
    if metric not in METRIC:
        raise ValueError(f"metric must be one of {list(METRIC)}")
    key, label = METRIC[metric]
    npz = _resolve_npz(source)
    z = np.load(npz)
    if title is None:
        title = f"{Path(npz).stem.replace('_camp','')}   ·  {label}"
    return build_patch_figure(z["V"], z["F"], z["face_labels"], z[key],
                              value_label=label, cmin=cmin, cmax=cmax,
                              colorscale=colorscale, clip=clip, title=title, controls=controls)
