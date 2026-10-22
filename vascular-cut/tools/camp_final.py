"""The cohort figure: repair only where it is verified, trim only far-out patches.

    python tools/camp_final.py

Two decisions, each made on measurement rather than preference:

REPAIR is applied only to cases where fixing the folded geometry moved delta K by
more than that case's own partition noise (camp_noise.py measures the noise by
re-running the SAME surface at 98k/100k/102k). Over 30 cases that is two: KP35_1
(0.02x against a 4.06x noise band, 46 folds, 174.6 deg dihedral) and KY244_1
(0.29x against 1.39x). One case got verifiably WORSE (KK28_11) and is left alone;
for eleven more the change sat inside the noise and is not evidence either way.

TRIM is 0.0025 -- a quarter of a percent per tail, about 6 patches and 0.5% of the
surface. Not chosen to make anything look better: it is where the measurement
becomes reproducible. Worst-case noise across the cohort falls from 23.5x at trim
0 to 3.8x here, and the folds do not explain that instability (correlation between
fold count and noise: r = -0.07; the worst offender, DJ18_2 at 23.5x, has none).
So the trim is doing a different job from the repair, and both are needed.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from camp_cohort import COHORT_STYLE, NORMAL_COHORTS, cohort_of   # noqa: E402

REPAIRED = ("KP35_1", "KY244_1")     # verified outside their own noise band
TRIM = "0.0025"
FACES = 100_000

# Arch type, read from the lab's own classification rather than inferred here --
# tools/bovine.csv is the sheet, one row per scan. It covers scans that are not in
# the cut cohort (the KT series), so the join is by scan id and anything missing on
# either side is reported rather than assumed.
BOVINE_CSV = HERE / "bovine.csv"
BOVINE_STYLE = {                       # label -> (colour, matplotlib marker)
    "bovine": ("#EB6834", "^"),
    "non-bovine": ("#2A78D6", "o"),
}


def load_bovine(path=BOVINE_CSV):
    out = {}
    for i, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines()):
        if not line.strip() or i == 0:
            continue
        scan, flag = [x.strip() for x in line.split(",")[:2]]
        out[scan] = int(flag)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--noise", default=str(HERE / "camp_noise.jsonl"))
    ap.add_argument("--trim", default=TRIM)
    ap.add_argument("--label", action="store_true")
    ap.add_argument("--by", choices=("cohort", "bovine"), default="cohort",
                    help="colour by patient cohort (default) or by arch type")
    ap.add_argument("--mark-repaired", action="store_true",
                    help="ring the cases whose folded geometry was repaired "
                         "(internal review; off for anything shared)")
    args = ap.parse_args()

    rows = [json.loads(l) for l in
            Path(args.noise).read_text(encoding="utf-8").splitlines() if l.strip()]
    by = {(r["case"], r["faces_target"], r["repaired"]): r for r in rows}
    cases = sorted({r["case"] for r in rows})

    data = []
    for c in cases:
        rep = c in REPAIRED
        r = by.get((c, FACES, rep)) or by.get((c, FACES, False))
        data.append({"case": c, "cohort": cohort_of(c), "repaired": rep,
                     "A_mm2": r["A_mm2"], "k_fluct": r["dk"][args.trim],
                     "folds": by[(c, FACES, True)].get("repair_folded_edges", 0)})

    ref_rows = [d for d in data if d["cohort"] in NORMAL_COHORTS]
    DK_REF = float(np.mean([d["k_fluct"] for d in ref_rows]))
    for d in data:
        d["k_fluct_norm"] = d["k_fluct"] / DK_REF
    print(f"Normal reference = mean of {len(ref_rows)} "
          f"{'/'.join(NORMAL_COHORTS)} scans = {DK_REF:.6g}  "
          f"(trim {args.trim}, {FACES:,} faces)")
    print(f"Repaired: {', '.join(REPAIRED)}  (verified outside their noise band)\n")

    if args.by == "bovine":
        bov = load_bovine()
        missing = [d["case"] for d in data if d["case"] not in bov]
        unused = [s for s in bov if s not in {d["case"] for d in data}]
        for d in data:
            d["bovine"] = bov.get(d["case"])
            d["group"] = ("bovine" if bov.get(d["case"]) == 1 else
                          "non-bovine" if bov.get(d["case"]) == 0 else "unlabelled")
        if missing:
            print(f"  NOT IN {BOVINE_CSV.name}: {', '.join(missing)}")
        if unused:
            print(f"  labelled but not in the cut cohort: {', '.join(sorted(unused))}")
        order = [g for g in BOVINE_STYLE if any(d["group"] == g for d in data)]
        order += sorted({d["group"] for d in data} - set(order))
        style = BOVINE_STYLE
        header = "arch"
    else:
        for d in data:
            d["group"] = d["cohort"]
        order = [c for c in COHORT_STYLE if any(d["group"] == c for d in data)]
        order += sorted({d["group"] for d in data} - set(order))
        style = {k: (v[0], v[2]) for k, v in COHORT_STYLE.items()}
        header = "cohort"

    print(f"  {header:12s} {'n':>3} {'median':>8} {'mean':>8}")
    for g in order:
        v = np.array([d["k_fluct_norm"] for d in data if d["group"] == g])
        note = ("   (normal)" if args.by == "cohort" and g in NORMAL_COHORTS else "")
        print(f"  {g:12s} {len(v):3d} {np.median(v):8.3f} {v.mean():8.3f}{note}")

    stem = (f"camp_cohort_final_trim{args.trim.replace('.', 'p')}"
            + ("_bovine" if args.by == "bovine" else ""))
    cols = ["case", "cohort", "A_mm2", "k_fluct", "k_fluct_norm", "repaired",
            "folds"] + (["bovine"] if args.by == "bovine" else [])
    (HERE / f"{stem}.csv").write_text(
        "\n".join([",".join(cols)]
                  + [",".join(str(d[k]) for k in cols) for d in data]) + "\n",
        encoding="utf-8")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.6, 7.0))
    ax.axhline(1.0, ls="--", color="#898781", lw=1.2, zorder=1)
    ax.annotate(f"normal ({' + '.join(NORMAL_COHORTS)}) mean", xy=(0.012, 1.0),
                xycoords=("axes fraction", "data"), va="bottom", ha="left",
                fontsize=9, color="#6b6a66")
    for g in order:
        sub = [d for d in data if d["group"] == g]
        colour, mk = style.get(g, ("#898781", "o"))
        tag = (" (normal)" if args.by == "cohort" and g in NORMAL_COHORTS else "")
        ax.scatter([d["A_mm2"] for d in sub], [d["k_fluct_norm"] for d in sub],
                   s=95, c=colour, marker=mk, edgecolor="white", linewidth=1.0,
                   label=f"{g}{tag} (n={len(sub)})", zorder=3)
    for d in data:
        # Ringing the two repaired cases is a DIAGNOSTIC, not part of the result --
        # it invites the reader to treat them as a separate class when they are
        # simply two meshes whose segmentation defects were corrected before
        # measurement, like any other preprocessing. Off by default; --mark-repaired
        # brings it back for internal review.
        if d["repaired"] and args.mark_repaired:
            ax.scatter([d["A_mm2"]], [d["k_fluct_norm"]], s=260, facecolors="none",
                       edgecolors="#333", linewidth=1.2, zorder=4)
        if args.label:
            ax.annotate(d["case"], (d["A_mm2"], d["k_fluct_norm"]),
                        textcoords="offset points", xytext=(6, 4), fontsize=6.5,
                        color="#555")
    ax.set_xlabel("Surface area (mm²)")
    ax.set_ylabel("δK / normal mean")
    ax.set_title(f"GT cohort: surface area vs normalized δK   (normals = 1.0)\n"
                 f"CAMP at {FACES:,} faces, trim {args.trim}, n={len(data)}",
                 fontsize=11)
    ax.grid(alpha=0.25, zorder=0)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / f"{stem}.png", dpi=170)
    print(f"\n  wrote {stem}.png and {stem}.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
