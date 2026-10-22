# Bare-bones web front end for the vascular cut tool and the shape metrics.
#
# SUPERSEDED by webapp/cut_app.py, which has both of these sections with the live cut
# planes and a 3-D view of the CAMP partition. Kept because it needs no VTK render
# window and no server-side rendering, so it still runs where that one will not, and
# because it is a second opinion: if the two ever disagree on a number, one of them is
# wrong and it is worth knowing which.
#
# Two sections that do not depend on each other: cut a geometry, or measure one you
# already have. Deliberately plain -- this exists so the pipeline is reachable from a
# browser and so UI work has something to sit on, not because these are the visuals.
#
#   streamlit run webapp/app.py
#
# Uploads are written to a per-session temporary directory and deleted when the
# session ends. Nothing is retained.

import atexit
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import streamlit as st
import trimesh

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vascular_api import CutSession, CutError            # noqa: E402
import shape_metrics as SM                               # noqa: E402

st.set_page_config(page_title="Vascular geometry tools", layout="wide")

PREVIEW_FACE_BUDGET = 60_000   # above this a browser starts to struggle with Mesh3d


# ── per-session scratch, cleaned up on exit ──────────────────────────────────
def session_dir():
    if "tmpdir" not in st.session_state:
        d = tempfile.mkdtemp(prefix="vascweb_")
        st.session_state.tmpdir = d
        atexit.register(shutil.rmtree, d, True)
    return Path(st.session_state.tmpdir)


def save_upload(upload):
    """Uploads live only for this session."""
    p = session_dir() / upload.name
    p.write_bytes(upload.getbuffer())
    return p


def mesh_figure(mesh, title=""):
    v = np.asarray(mesh.vertices, float)
    f = np.asarray(mesh.faces, int)
    fig = go.Figure(go.Mesh3d(
        x=v[:, 0], y=v[:, 1], z=v[:, 2],
        i=f[:, 0], j=f[:, 1], k=f[:, 2],
        color="#c44e52", opacity=1.0, flatshading=False,
        lighting=dict(ambient=0.45, diffuse=0.85, specular=0.25, roughness=0.6),
        hoverinfo="skip",
    ))
    fig.update_layout(
        title=title, margin=dict(l=0, r=0, t=30, b=0), height=560,
        scene=dict(aspectmode="data", xaxis_visible=False,
                   yaxis_visible=False, zaxis_visible=False),
    )
    return fig


st.title("Vascular geometry tools")
tab_cut, tab_shape = st.tabs(["Geometry: cut and export", "Shape metrics: delta K"])

# ═════════════════════════════════════════════════════════════════════════════
# SECTION 1: cut a closed mesh into an open one
# ═════════════════════════════════════════════════════════════════════════════
with tab_cut:
    st.caption("A closed surface plus its centerline, opened at every vessel end.")
    c1, c2, c3 = st.columns([2, 2, 1])
    up_cl = c1.file_uploader("Centerline (.csv)", type=["csv"], key="cl")
    up_mesh = c2.file_uploader("Mesh (.mat or .stl)", type=["mat", "stl"], key="mesh")
    family = c3.selectbox("Anatomy", ["auto", "thoracic", "abdominal"],
                          help="Auto reads it from the geometry. Override when "
                               "the openings come out named as the wrong vessels.")

    if up_cl and up_mesh:
        sig = (up_cl.name, up_cl.size, up_mesh.name, up_mesh.size, family)
        if st.session_state.get("cut_sig") != sig:
            with st.spinner("Loading and cutting…"):
                try:
                    st.session_state.session = CutSession(
                        save_upload(up_cl), save_upload(up_mesh), family=family)
                    st.session_state.cut_sig = sig
                except CutError as e:
                    st.session_state.session = None
                    st.error(str(e))

        s = st.session_state.get("session")
        if s is not None:
            anat = s.anatomy_summary()
            st.info(f"Treated as **{anat['kind']}** ({anat['how']}). "
                    f"{anat['named_openings']} opening(s) named anatomically, "
                    f"assuming an LPS frame. Confirm left/right on a new case.")

            left, right = st.columns([1, 2])

            with left:
                st.subheader("Region")
                names = ["Whole vessel"] + [p["name"]
                                            for p in s.presets(planar=True)]
                pick = st.selectbox("Isolate", names, key="region_pick")
                if pick == "Whole vessel":
                    s.clear_region()
                else:
                    s.apply_preset(pick)
                    # s1 None is an OPEN distal end -- one plane, everything below it
                    # kept, which is how "Infrarenal aorta with iliacs" keeps the
                    # iliacs (their ostia sit at the bifurcation, outside any
                    # two-plane region). sorted() cannot order None against a float,
                    # so this used to raise TypeError and take the page down on the
                    # one preset an abdominal user most wants. Unpack, do not sort.
                    lo, hi = s.region
                    if hi is None:
                        a = st.slider("Start (mm along the vessel) — everything "
                                      "below this is kept",
                                      0.0, float(s.route.length), float(lo),
                                      step=0.5)
                        if a != lo:
                            s.set_region(a, None, pick, manual=True)
                    else:
                        lo, hi = min(lo, hi), max(lo, hi)
                        a, b = st.slider("Stations (mm along the vessel)",
                                         0.0, float(s.route.length),
                                         (float(lo), float(hi)), step=0.5)
                        if (a, b) != (lo, hi):
                            s.set_region(a, b, pick, manual=True)

                st.subheader("Openings")
                st.caption("How far back from each vessel end to cut.")
                for o in s.openings():
                    v = st.number_input(
                        o["name"], min_value=1.0, max_value=200.0,
                        value=float(o["offset_mm"]), step=0.5,
                        key=f"off_{o['label']}",
                        disabled=o["is_region_bound"],
                        help="Set by the region" if o["is_region_bound"] else None)
                    if not o["is_region_bound"]:
                        s.set_offset(o["label"], v)

                s.clean_rims = st.checkbox(
                    "Sliver-free rim", value=s.clean_rims,
                    help="Snaps the rim onto the cut plane. Slower, but removes the "
                         "near-degenerate triangles that upset curvature fitting.")

            with st.spinner("Cutting…"):
                status = s.recut()

            with right:
                if status["ok"]:
                    st.success(status["message"])
                    m1, m2, m3 = st.columns(3)
                    m1.metric("Surface area", f"{status['area_mm2']:,.0f} mm²")
                    m2.metric("Openings", status["n_openings"])
                    m3.metric("Faces", f"{status['n_faces']:,}")
                else:
                    st.error(status["message"])

                if s.mesh is not None and len(s.mesh.faces) <= PREVIEW_FACE_BUDGET:
                    st.plotly_chart(mesh_figure(s.mesh), use_container_width=True)
                elif s.mesh is not None:
                    st.warning(f"{len(s.mesh.faces):,} faces is too many to draw "
                               f"interactively; the download is unaffected.")

            if status["ok"]:
                st.divider()
                stem = Path(up_mesh.name).stem
                d1, d2, d3 = st.columns(3)
                d1.download_button(
                    "Download open mesh (.stl)",
                    data=s.mesh.export(file_type="stl"),
                    file_name=f"{stem}_open.stl", mime="model/stl")
                import json
                d2.download_button(
                    "Download sidecar (.json)",
                    data=json.dumps(s.sidecar(), indent=2),
                    file_name=f"{stem}_cut_planes.json", mime="application/json")
                if d3.button("Send to shape metrics →"):
                    p = session_dir() / f"{stem}_open.stl"
                    s.mesh.export(str(p))
                    st.session_state.handoff = str(p)
                    st.session_state.handoff_euler = int(s.mesh.euler_number)
                    st.success("Ready in the shape-metrics tab.")
    else:
        st.info("Upload a centerline and a mesh to begin.")

# ═════════════════════════════════════════════════════════════════════════════
# SECTION 2: curvature fluctuation on any surface
# ═════════════════════════════════════════════════════════════════════════════
with tab_shape:
    st.caption("Any open triangular surface. It does not have to come from the "
               "cut tool.")
    # Which build produced the number matters -- several versions of CAMP exist and
    # they are different methods, so say which one is answering.
    st.caption(f"CAMP source: `{SM.campsrc_origin()}`")

    handoff = st.session_state.get("handoff")
    if handoff:
        st.success(f"Using **{Path(handoff).name}** from the geometry tab. "
                   f"Upload another file to replace it.")
    up_stl = st.file_uploader("Surface (.stl)", type=["stl"], key="shape_stl")

    c1, c2, c3 = st.columns(3)
    target = c1.select_slider(
        "Target faces", options=[25_000, 50_000, 100_000, 200_000],
        value=SM.REFERENCE["target_faces"],
        help="The reference cohort ran at 100,000. δK moves with face count, so "
             "matching it is what makes a comparison meaningful.")
    trim = c2.select_slider(
        "Trim fraction", options=[0.0, 0.001, 0.0025, 0.005, 0.01],
        value=SM.REFERENCE["trim_frac"],
        help="Drops the highest and lowest integrated-K patches — both tails, so the "
             "number here is per tail. The reference is 0.0025, chosen because at "
             "trim 0 the same surface remeshed to 98k/100k/102k gives δK values that "
             "differ by 2.6× typically and 23× at worst. Reference and upload must "
             "match.")
    n_seeds = c3.select_slider(
        "Partitions", options=[1, 2, 4, 8, 16], value=1,
        help="How many independent k-means partitions to average. One partition is "
             "one draw: δK is carried by ~22 patches out of ~850, sited where "
             "curvature changes fastest, and where k-means puts those boundaries is "
             "arbitrary. Re-seeding moves δK by ~16%. The error falls as 1/√n.")

    st.caption(
        f"**Partition uncertainty**: one partition carries about "
        f"±{100 * SM.PARTITION_CV:.0f}%; averaging {n_seeds} brings that to "
        f"±{100 * SM.partition_band(n_seeds):.0f}%. Measured over 30 aortas × 8 seeds "
        f"× 3 densities, and the draws are independent — a 4-seed mean matched "
        f"CV/√4 with a median ratio of 1.03. Two δK values closer together than "
        f"2× the error bar are not distinguishable. "
        f"**Comparable only at a fixed face count**: δK drifts ~"
        f"{100 * SM.DENSITY_CV:.0f}% across a 4× density change because "
        f"min_patch_size is a vertex count, and averaging does not remove that. "
        f"See CURVATURE.md.")

    source = save_upload(up_stl) if up_stl else handoff
    if source and st.button("Compute δK", type="primary"):
        est = 10 * n_seeds
        with st.spinner(f"Partitioning at {target:,} faces × {n_seeds}, "
                        f"roughly {est} to {est * 4} seconds..."):
            try:
                res = SM.delta_k_averaged(source, target_faces=target, trim_frac=trim,
                                          seeds=range(42, 42 + n_seeds),
                                          make_figure=True)
                st.session_state.shape_result = res
            except Exception as e:
                st.session_state.shape_result = None
                st.error(f"{type(e).__name__}: {e}")

    res = st.session_state.get("shape_result")
    if res:
        m = st.columns(4)
        se = res.get("k_fluct_se", 0.0)
        # "No error bar" was the wrong way to say this. A single partition has an
        # error bar; it just cannot be estimated FROM ONE DRAW, so the cohort figure
        # stands in -- PARTITION_CV, over 30 aortas x 8 seeds. Saying "no error bar"
        # reads as "exact", which is the opposite of the truth.
        m[0].metric("δK (k_fluct)", f"{res['k_fluct']:.4g}",
                    delta=(f"± {se:.2g} SE ({100 * se / res['k_fluct']:.0f}%)"
                           if se else
                           f"± {100 * SM.PARTITION_CV:.0f}% — one partition, "
                           f"cohort figure"),
                    delta_color="off")
        m[1].metric("Surface area", f"{res['A_mm2']:,.0f} mm²")
        m[2].metric("Patches", f"{res['n_patches']:,}")
        m[3].metric("At vertex floor", f"{res['pct_at_min_size']:.1f}%")

        if res.get("n_seeds", 1) > 1:
            vals = res.get("k_fluct_all", [])
            st.caption(
                f"Mean of {res['n_seeds']} independent partitions: "
                + ", ".join(f"{v:.4g}" for v in vals)
                + f".  The patches drawn below are the **first** partition "
                  f"({res['k_fluct_seed']:.4g}) — labels cannot be averaged across "
                  f"partitions that carve the surface differently, so the picture is "
                  f"one draw and the number is the mean.")

        lo, hi = SM.REFERENCE["area_mm2"]
        if not res["in_reference_range"]:
            st.warning(
                f"This surface is {res['A_mm2']:,.0f} mm², outside the reference "
                f"cohort's {lo:,.0f}–{hi:,.0f} mm². δK is still computed, but there "
                f"is no reference curve underneath it to place the value against.")
        f_lo, f_hi = SM.REFERENCE["pct_at_min_size"]
        if not (f_lo <= res["pct_at_min_size"] <= f_hi):
            st.warning(
                f"{res['pct_at_min_size']:.1f}% of patches sit at the vertex floor, "
                f"outside the reference's {f_lo}–{f_hi}%. The partition is coarser "
                f"or finer than the cohort's. Try a different target face count.")

        # Known before the run and worth saying after it too: a fold is geometry the
        # vessel does not have, and no amount of trimming turns it into anatomy.
        try:
            d = SM.surface_defects(trimesh.load(source, process=False))
        except Exception:                                   # noqa: BLE001
            d = None
        if d and d["folds"]:
            st.warning(
                f"{d['folds']} folded edge(s) in the source surface, worst dihedral "
                f"{d['max_dihedral']:.0f}°. The surface doubles back through itself "
                f"there — a segmentation defect, not anatomy. δK is not corrected for "
                f"it: across the reference cohort, fold count and how much δK moves "
                f"between runs correlate at −0.07, so folds flag a bad input rather "
                f"than explain a high value.")

        if res["figure"] is not None:
            st.plotly_chart(res["figure"], use_container_width=True)

        st.download_button(
            "Download per-patch table (.csv)",
            data=SM.patch_table(res).to_csv(index=False),
            file_name="patches.csv", mime="text/csv")
    elif not source:
        st.info("Upload a surface, or cut one in the geometry tab first.")
