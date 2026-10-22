"""delta K against surface area, with each case's geometry above its own point.

A scatter says there is a trend; a montage says the tool ran on everything; neither
shows WHICH anatomy sits where on the trend. This puts both in one frame.

Defaults to the MESH SWEEP: each case's delta K averaged over 100k / 200k / 300k /
400k faces, with the SD across those densities as the error bar. That is a stronger
statement than any single density -- it shows the cohort keeps its order across a
fourfold change in mesh resolution, which is the property a reader should care about.

    .venv\\Scripts\\python tools\\dk_area_thumbs.py                     # mesh sweep
    .venv\\Scripts\\python tools\\dk_area_thumbs.py --preset seeds      # 8 partitions

PLACING THE THUMBNAILS -- two failed approaches, recorded because the reasoning
generalises to any figure of this shape:

  1. Evenly spaced in a strip with leader lines down to the points. The leaders cross
     badly. The strip is EVENLY spaced while the axis is by ACTUAL area, and half this
     cohort sits under 40,000 mm2, so the bunched left-hand cases fan out.
  2. Even spacing and true position cannot both be had. True position is the one that
     carries meaning, so thumbnails sit at their real x and the bunching is handled by
     staggering over rows -- consecutive cases in area order alternate, so neighbours
     are never side by side.

With the thumbnails at true x a STRICTLY VERTICAL drop line makes the mapping explicit
and cannot cross anything. That is what is drawn.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import trimesh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.offsetbox import AnnotationBbox, OffsetImage        # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# Serif throughout, including the maths, so the delta and the axis numerals belong to
# the same family as the labels. STIX is Times-metric, which is what mathtext needs.
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.edgecolor": "#333333",
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

VIEW = np.array([-1.0, 0.0, 0.0])        # left lateral: arch apex top-left, descending
                                         # sweeping down-right. Matches the reference
                                         # the lab uses for compartment figures.
UP = np.array([0.0, 0.0, 1.0])
# Coloured and marked by PATIENT PREFIX, not by disease group -- the same six-way
# scheme as tools/camp_robust_plot.py and the mean-of-four-densities figure, so the
# paper's plots agree with each other and a reader can follow one case between them.
# Shared verbatim rather than re-picked: two palettes for the same cohort is a way to
# make the same figure look like two different results.
STYLE = {
    "KY": ("#1f77b4", "o", "KY (normal)"), "SA": ("#2ca02c", "s", "SA (peds/normal)"),
    "KK": ("#ff7f0e", "^", "KK"),          "KR": ("#d62728", "D", "KR"),
    "DJ": ("#7d3ac1", "v", "DJ"),          "KP": ("#2a9d54", "P", "KP"),
}


def prefix(case):
    """KY7_1 -> KY. Two letters, which is what every case id in this cohort uses."""
    return case[:2].upper()


def style_of(case):
    return STYLE.get(prefix(case), ("#666666", "o", prefix(case)))
MESHES = ROOT / "output" / "nubs_final"

# Ten cases rather than thirty, so each can be shown large enough to read INSIDE the
# axes. Chosen to span both axes and to carry the cases that make an argument, not by
# rank alone:
#   KY7_1, KY68_1   the small healthy end
#   KY228_1         a healthy aorta of middling size
#   DJ18_2, KK15_1  diseased of near-identical area, delta K 1.60 against 1.00 --
#                   area alone does not decide the number
#   KY160_1         a NORMAL as large as the diseased around it, sitting at 1.22 while
#                   they sit at 1.6-2.1. This read 3.85 before the partition fixes,
#                   inside the diseased range; it is the correction made visible
#   KK28_11, KK20_2 the diseased top
#   KR27_2          large and diseased but low -- the counterexample to "bigger is more"
#   KK32_0          the largest surface in the cohort
SELECTION = ["KY7_1", "KY228_1", "DJ18_2", "KK15_1", "KY160_1",
             "KK28_11", "KR27_2", "KK32_0"]

# Where each thumbnail sits relative to its point, in axes fractions. Placed by hand:
# ten points is few enough to position deliberately and too few for an automatic
# rule to beat an eye. dx, dy are offsets from the point.
OFFSETS = {
    # Close to their points now, and no connector lines: at this distance the pairing
    # reads by proximity alone, and the lines were the busiest thing in the frame.
    # Direction chosen per case to sit in whichever side is clear.
    "KY7_1":   (0.026, 0.135), "KY228_1": (-0.026, 0.140),
    "DJ18_2":  (-0.030, 0.130), "KK15_1": (0.028, -0.125),
    "KY160_1": (-0.030, -0.130), "KK28_11": (0.030, 0.115),
    "KR27_2":  (0.028, -0.130), "KK32_0":  (-0.046, -0.020),
}

PRESETS = {
    "sweep": dict(csv="tools/nub_meanSD_100kplus.csv", value="mean_dk_norm", err="sd",
                  ylab="$\\delta K$, normalised",
                  note="Mean over four mesh densities, ± SD.",
                  # Wider y range than the seeds preset, so the fractions that place a
                  # thumbnail into the empty lower band are smaller here.
                  offsets={"KK15_1": (0.028, -0.085), "KY160_1": (-0.030, -0.095),
                           "KR27_2": (0.028, -0.095), "KK32_0": (-0.046, -0.020)}),
    "seeds": dict(csv="tools/cohort_seedavg.csv", value="norm_mean8", err=None,
                  ylab="$\\delta K$, normalised",
                  note="Mean of 8 CAMP partitions at 100k faces, ± SE."),
}


def thumbnail(cid, colour, px=240):
    """One case rendered small on a transparent ground, in its point's colour."""
    stl = MESHES / f"{cid}_nub.stl"
    if not stl.exists():
        return None
    m = trimesh.load(str(stl), process=False)
    m.merge_vertices()
    poly = pv.PolyData(np.asarray(m.vertices),
                       np.hstack([np.full((len(m.faces), 1), 3),
                                  np.asarray(m.faces)]).ravel())
    pl = pv.Plotter(off_screen=True, window_size=(px, int(px * 1.5)))
    pl.add_mesh(poly, color=colour, smooth_shading=True, specular=0.18,
                diffuse=0.92, ambient=0.30)
    pl.view_vector(VIEW, viewup=UP)
    pl.reset_camera()
    pl.camera.zoom(1.15)
    img = pl.screenshot(None, transparent_background=True, return_img=True)
    pl.close()
    return np.asarray(img)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", choices=sorted(PRESETS), default="sweep")
    ap.add_argument("--thumb-px", type=int, default=240)
    ap.add_argument("--zoom", type=float, default=0.26)
    ap.add_argument("--all", action="store_true",
                    help="every case; default is the SELECTION below")
    ap.add_argument("--title", default="$\delta K$ against surface area")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    spec = PRESETS[args.preset]
    out = args.out or str(ROOT / "docs" / "figures" / f"dk_area_thumbs_{args.preset}.png")

    rows = list(csv.DictReader(open(ROOT / spec["csv"], encoding="utf-8")))
    rows.sort(key=lambda r: float(r["area_mm2"]))
    n = len(rows)
    area = np.array([float(r["area_mm2"]) for r in rows])
    val = np.array([float(r[spec["value"]]) for r in rows])
    if spec["err"]:
        err = np.array([float(r[spec["err"]]) for r in rows])
    else:                                    # seeds preset: SE from the 8-seed SD
        sd = np.array([float(r["sd8"]) for r in rows])
        m8 = np.array([float(r["dk_mean8"]) for r in rows])
        err = (sd / np.sqrt(8.0)) * np.where(m8 > 0, val / np.where(m8 > 0, m8, 1), 0)
    groups = [r["group"] for r in rows]
    rows_by_case = {r["case"]: r for r in rows}
    print(f"{n} cases, {args.preset}: {spec['value']} +/- {spec['err'] or 'SE'}")

    fig, ax = plt.subplots(figsize=(17, 11))
    fig.subplots_adjust(left=0.075, right=0.985, top=0.905, bottom=0.085)

    cases = [r["case"] for r in rows]
    for pref in STYLE:
        m = [i for i, c in enumerate(cases) if prefix(c) == pref]
        if not m:
            continue
        col, mk, lab = STYLE[pref]
        ax.errorbar(area[m], val[m], yerr=err[m], fmt=mk, ms=9, color=col,
                    ecolor=col, elinewidth=1.4, capsize=3.5, zorder=4,
                    markeredgecolor="white", markeredgewidth=0.9,
                    label=f"{lab} (n = {len(m)})")
    # The normal mean, which every value is expressed against. Without it "normalised"
    # is a word in the axis label rather than a line on the page.
    ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1.0, zorder=1)

    ax.set_xlabel("Surface area (mm$^2$)", fontsize=17, labelpad=9)
    ax.set_ylabel(spec["ylab"], fontsize=17, labelpad=9)
    ax.tick_params(labelsize=14)
    ax.grid(axis="y", alpha=0.16, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(loc="upper left", frameon=True, framealpha=0.95, fontsize=12.5,
              edgecolor="#cccccc", handletextpad=0.6)

    # Headroom so a thumbnail placed above a high point is not clipped.
    ylo, yhi = ax.get_ylim()
    # Headroom below as well as above, for the thumbnails placed into the empty
    # lower-right -- but NEVER below zero. delta K is an area-weighted variance, so a
    # negative value does not exist and an axis that shows one is simply wrong. The
    # sweep preset has a much wider range and drove the floor to -1 before this.
    # Modest margins now that nothing is parked far from its point. Never below
    # zero: delta K is an area-weighted variance and cannot be negative.
    ax.set_ylim(max(0.0, ylo - (yhi - ylo) * 0.05), yhi + (yhi - ylo) * 0.10)
    xlo, xhi = ax.get_xlim()
    ax.set_xlim(xlo, xhi + (xhi - xlo) * 0.02)

    # Positions hand-placed with tools/place_thumbs.py, if any have been saved.
    _json = ROOT / "tools" / "thumb_offsets.json"
    placed = (json.loads(_json.read_text()).get(args.preset, {})
              if _json.exists() else {})
    if placed:
        print(f"  using {len(placed)} hand-placed position(s)")

    show = set(rows_by_case) if args.all else set(SELECTION)
    missing = show - {r["case"] for r in rows}
    if missing:
        print(f"  not in the data: {sorted(missing)}")
    print(f"rendering {len(show & {r['case'] for r in rows})} thumbnails...")

    for r, x, y in zip(rows, area, val):
        cid = r["case"]
        if cid not in show:
            continue
        img = thumbnail(cid, style_of(cid)[0], args.thumb_px)
        if img is None:
            print(f"  missing mesh: {cid}")
            continue
        # Hand-placed positions from tools/place_thumbs.py win over the built-in
        # ones. Delete tools/thumb_offsets.json to go back to these defaults.
        dx, dy = placed.get(cid) or spec.get("offsets", {}).get(
            cid, OFFSETS.get(cid, (0.0, 0.16)))
        # Offsets are axes fractions, so a hand-placed figure survives a change of
        # data range or figure size without every position needing to be redone.
        tx = x + dx * (ax.get_xlim()[1] - ax.get_xlim()[0])
        ty = y + dy * (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.add_artist(AnnotationBbox(
            OffsetImage(img, zoom=args.zoom), (tx, ty), frameon=False,
            box_alignment=(0.5, 0.5), zorder=5, annotation_clip=False))
        ax.annotate(cid, (tx, ty), fontsize=12, color=style_of(cid)[0],
                    xytext=(0, -40), textcoords="offset points",
                    ha="center", zorder=6)

    fig.text(0.075, 0.962, args.title, fontsize=20, color="#1b2733")
    fig.text(0.075, 0.936, spec["note"], fontsize=11.5, color="#66707a")

    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
