"""Delta K of the arch, bovine against non-bovine.

    python tools/arch_bovine_plot.py                       # read the 50k run, plot
    python tools/arch_bovine_plot.py --run tools/camp_arch_50k.jsonl
    python tools/arch_bovine_plot.py --exports output/arch_final

WHAT IS BEING COMPARED, and what a reader has to keep in mind

Delta K is EXTENSIVE. It is an area-weighted variance of integrated Gaussian
curvature over a partition, so it grows with how much surface there is and with how
that surface was partitioned. Two numbers are only comparable when cut extent, face
count and trim all match -- which is why the run this reads is pinned to one operating
point (50,000 faces, trim 0) and why the plot plots delta K AGAINST AREA rather than
on its own. A group difference that is really an area difference would show up as the
two clouds sitting at different places along x, not as a vertical separation.

WHERE THE GROUPS COME FROM

Not from a spreadsheet. name_anatomy already decides it while naming the arch -- a
brachiocephalic with THREE daughters instead of two is a bovine arch, the left common
carotid sharing the innominate origin -- and records it as a note in the export
sidecar. So the label travels with the geometry it describes, and a case relabelled by
a better naming rule is relabelled here too.

The test is a difference in DISTRIBUTION, not in means: eight against twenty-two, and
delta K is skewed by construction. Mann-Whitney U on the ranks, reported with the
group medians and the overlap, and no claim made beyond what that supports.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent


def is_bovine(sidecar):
    """A bovine arch, as name_anatomy decided it when the case was named."""
    d = json.loads(Path(sidecar).read_text(encoding="utf-8"))
    notes = (d.get("anatomy") or {}).get("notes") or []
    names = [o.get("name", "") for o in d.get("openings", [])]
    return any("bovine" in str(x).lower() for x in list(notes) + names)


def load(run, exports):
    rows = [json.loads(l) for l in Path(run).read_text(encoding="utf-8").splitlines()
            if l.strip()]
    out = []
    for r in rows:
        case = str(r.get("case") or Path(r["file"]).stem).replace("_arch", "")
        side = Path(exports) / f"{case}_arch_cut_planes.json"
        if not side.exists():
            print(f"  {case}: no sidecar, skipped (cannot say bovine or not)")
            continue
        out.append({"case": case, "dk": float(r["k_fluct"]),
                    "area": float(r["A_mm2"]), "faces": int(r["target_faces"]),
                    "patches": int(r["n_patches"]),
                    "folded": int(r.get("repair_folded_edges") or 0),
                    "bovine": is_bovine(side)})
    return out


def mann_whitney(a, b):
    """U and a normal-approximation p, two-sided, with ties handled by mid-ranks."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan")
    allv = np.concatenate([a, b])
    order = allv.argsort()
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # mid-ranks for ties
    for v in np.unique(allv):
        m = allv == v
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u = min(u1, n1 * n2 - u1)
    mu = n1 * n2 / 2.0
    _, counts = np.unique(allv, return_counts=True)
    tie = (counts ** 3 - counts).sum()
    n = n1 + n2
    sd = np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1))))
    if sd == 0:
        return u, float("nan")
    z = (abs(u - mu) - 0.5) / sd
    from math import erfc, sqrt
    return u, float(erfc(z / sqrt(2.0)))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=str(HERE / "camp_arch_50k.jsonl"))
    ap.add_argument("--exports", default=str(ROOT / "output" / "arch_final"))
    ap.add_argument("--out", default=str(HERE / "arch_bovine_50k.png"))
    ap.add_argument("--csv", default=str(HERE / "arch_bovine_50k.csv"))
    args = ap.parse_args()

    rows = load(args.run, args.exports)
    if not rows:
        print("nothing to plot")
        return 1
    bov = [r for r in rows if r["bovine"]]
    non = [r for r in rows if not r["bovine"]]
    faces = {r["faces"] for r in rows}
    print(f"\n  {len(rows)} arches at {'/'.join(f'{f:,}' for f in sorted(faces))} "
          f"faces: {len(bov)} bovine, {len(non)} non-bovine")

    def stat(g, key):
        v = np.array([r[key] for r in g], float)
        return (np.median(v), np.percentile(v, 25), np.percentile(v, 75),
                v.min(), v.max())

    for label, g in (("bovine", bov), ("non-bovine", non)):
        if not g:
            continue
        m, q1, q3, lo, hi = stat(g, "dk")
        am = stat(g, "area")[0]
        print(f"    {label:11s} n={len(g):2d}  delta K median {m:8.4f} "
              f"[{q1:.4f}, {q3:.4f}]  range {lo:.4f}-{hi:.4f}   area median "
              f"{am:,.0f} mm2")

    u, p = mann_whitney([r["dk"] for r in bov], [r["dk"] for r in non])
    print(f"\n    Mann-Whitney U = {u:.1f}, two-sided p = {p:.3f}")
    if p == p and p > 0.05:
        print("    Not separable at this n. Reported as measured, not as a finding.")

    # Area is the confounder to rule out before reading anything into delta K.
    ua, pa = mann_whitney([r["area"] for r in bov], [r["area"] for r in non])
    print(f"    area, same test:      p = {pa:.3f}"
          + ("   (the groups differ in SIZE, so any delta K gap may be that)"
             if pa == pa and pa <= 0.05 else "   (areas are comparable)"))

    Path(args.csv).write_text(
        "case,bovine,delta_k,area_mm2,patches,folded_edges\n"
        + "\n".join(f"{r['case']},{int(r['bovine'])},{r['dk']:.6f},{r['area']:.1f},"
                    f"{r['patches']},{r['folded']}" for r in sorted(
                        rows, key=lambda r: (not r["bovine"], r["case"]))) + "\n",
        encoding="utf-8")
    print(f"\n  table -> {args.csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                      # noqa: BLE001
        print(f"  (no matplotlib: {e}; the table above is the result)")
        return 0

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.6),
                                 gridspec_kw={"width_ratios": [2, 1]})
    for g, c, mk, lab in ((non, "#4c72b0", "o", "non-bovine"),
                          (bov, "#c44e52", "^", "bovine")):
        if not g:
            continue
        ax.scatter([r["area"] for r in g], [r["dk"] for r in g], c=c, marker=mk,
                   s=46, alpha=0.85, edgecolor="white", linewidth=0.6, label=lab)
    ax.set_xlabel("surface area (mm$^2$)")
    ax.set_ylabel("delta K  (area-weighted variance of integrated K)")
    ax.set_yscale("log")
    ax.set_title(f"Aortic arch with nubs, {sorted(faces)[0]:,} faces, trim 0")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, linewidth=0.6)

    data = [[r["dk"] for r in non], [r["dk"] for r in bov]]
    bx.boxplot(data, tick_labels=[f"non-bovine\nn={len(non)}", f"bovine\nn={len(bov)}"],
               widths=0.55, showfliers=False)
    for i, g in enumerate((non, bov), start=1):
        x = np.random.default_rng(0).normal(i, 0.045, len(g))
        bx.scatter(x, [r["dk"] for r in g], s=26, alpha=0.8,
                   c="#4c72b0" if i == 1 else "#c44e52", zorder=3)
    bx.set_yscale("log")
    bx.set_ylabel("delta K")
    bx.set_title(f"p = {p:.3f} (Mann-Whitney)")
    bx.grid(alpha=0.25, axis="y", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"  plot  -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
