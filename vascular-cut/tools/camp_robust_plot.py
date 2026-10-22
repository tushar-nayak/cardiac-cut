"""Robustness, as the cohort plot repeated across conditions.

    python tools/camp_robust_plot.py density
    python tools/camp_robust_plot.py sweep

Each panel is normalised to its OWN KY+SA mean, which is the point: delta K is expected
to shift with mesh density and with surface noise, and a shift that moves every case
together is not a failure. Removing it per panel leaves only the thing that has to
hold -- whether the cohort keeps its shape.

The trajectory panel carries one line per case across the conditions. Lines that stay
flat and parallel mean the ordering survived; lines that cross mean it did not, and
crossings are legible in a way a single correlation is not.
"""
import glob
import json
import sys

import numpy as np
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402

MODE = sys.argv[1] if len(sys.argv) > 1 else "density"
TAG = f"_{sys.argv[2]}" if len(sys.argv) > 2 else ""
STYLE = {'KY': ("#1f77b4", "o", "KY (normal)"), 'SA': ("#2ca02c", "s", "SA (peds/normal)"),
         'KK': ("#ff7f0e", "^", "KK"), 'KR': ("#d62728", "D", "KR"),
         'DJ': ("#7d3ac1", "v", "DJ"), 'KP': ("#2a9d54", "P", "KP")}
TRIM = "0.0"

# An untagged cohort must not swallow a tagged one: 'camp_robust_density_*' also
# matches 'camp_robust_density_arch_0'. Shard files always end in a digit, so an
# untagged run is pinned to that.
if TAG:
    pats = [f'tools/camp_robust_{MODE}{TAG}_*.jsonl']
else:
    pats = [f'tools/camp_robust_{MODE}_[0-9]*.jsonl']
    if MODE == "density":
        pats.append('tools/camp_robust_hidensity_[0-9]*.jsonl')
rows = [json.loads(l) for pat in pats for f in glob.glob(pat)
        for l in open(f, encoding='utf-8') if l.strip()]
ok = [r for r in rows if 'error' not in r]
print(f"{len(rows)} runs, {len(rows)-len(ok)} errors")

if MODE == "density":
    key = lambda r: r['faces']                                     # noqa: E731
    label = lambda v: f"{v:,} faces"                               # noqa: E731
    ref_cond = 100_000
else:
    key = lambda r: r['sigma_mm']                                  # noqa: E731
    label = lambda v: f"sigma {v:g} mm"                            # noqa: E731
    ref_cond = None

groups = {}
for r in ok:
    groups.setdefault(key(r), {}).setdefault(r['case'], []).append(r['dk'][TRIM])
conds = sorted(groups)
# Only conditions with a usable set of cases, and only cases present in ALL of them --
# a partially finished run must degrade to a smaller plot rather than crash, and a
# panel normalised on a different set of cases than its neighbour is not comparable.
cases = sorted(set.intersection(*[set(groups[c]) for c in conds]))
if not cases or not any(c.startswith(('KY', 'SA')) for c in cases):
    counts = ", ".join(f"{label(c)}: {len(groups[c])}" for c in conds)
    raise SystemExit(f"  not enough overlap yet ({counts}); "
                     f"{len(cases)} cases common to all conditions, "
                     f"{sum(1 for c in cases if c.startswith(('KY','SA')))} of them normals")
area = {}
for r in ok:
    area[r['case']] = r['area']
A = np.array([area[c] for c in cases])
g = np.array([c.startswith(('KY', 'SA')) for c in cases])

# The reference column is the 100k run of whichever cohort this is. For the arches
# there is no prior 100k file, so fall back to this run's own 100k condition.
try:
    ref = {r['case']: r['dk'][TRIM] for r in (json.loads(l) for l in
           open('tools/camp_nub_bothfix_100k.jsonl', encoding='utf-8') if l.strip())
           if 'error' not in r}
    REF = np.array([ref[c] for c in cases])
except (FileNotFoundError, KeyError):
    ref_cond = 100_000 if 100_000 in groups else conds[-1]
    REF = np.array([np.mean(groups[ref_cond][c]) for c in cases])

n = len(conds)
fig, axes = plt.subplots(2, n, figsize=(6.6 * n, 13),
                         gridspec_kw={"height_ratios": [3, 2]})
if n == 1:
    axes = axes.reshape(2, 1)

norm = {}
for j, cond in enumerate(conds):
    ax = axes[0, j]
    Y = np.array([np.mean(groups[cond][c]) for c in cases])
    Yn = Y / Y[g].mean()
    norm[cond] = Yn
    for pre, (col, mk, lab) in STYLE.items():
        m = np.array([c.startswith(pre) for c in cases])
        if not m.any():
            continue
        ax.scatter(A[m], Yn[m], c=col, marker=mk, s=80, zorder=3,
                   label=f"{lab} (n={m.sum()})", edgecolors="white", linewidths=0.5)
    for i, (c, x, y) in enumerate(zip(cases, A, Yn)):
        ax.annotate(c, (x, y), textcoords="offset points",
                    xytext=(0, 9 if i % 2 == 0 else -15), ha="center",
                    fontsize=6.5, color=STYLE[c[:2]][0])
    ax.axhline(1.0, color="grey", ls="--", lw=1.1)
    rho = spearmanr(REF, Y)[0]
    ax.set_title(f"{label(cond)}\nnormals {Yn[g].min():.2f}-{Yn[g].max():.2f}, "
                 f"CV {Yn[g].std(ddof=1)/Yn[g].mean():.2f}   |   rank vs 100k ref "
                 f"rho {rho:+.2f}", fontsize=11)
    ax.set_xlabel("Surface area (mm$^2$)")
    ax.grid(alpha=.3)
    ax.set_ylim(0, 4.4)
    if j == 0:
        ax.set_ylabel(r"$\delta K$ / normal mean (within panel)", fontsize=11)
        ax.legend(frameon=True, fontsize=8, loc="upper left")

gs = axes[1, 0].get_gridspec()
for a in axes[1, :]:
    a.remove()
axt = fig.add_subplot(gs[1, :])
X = np.arange(len(conds))
for i, c in enumerate(cases):
    col = STYLE[c[:2]][0]
    y = [norm[cond][i] for cond in conds]
    axt.plot(X, y, "-o", color=col, lw=1.3, ms=4.5, alpha=0.8)
    axt.annotate(c, (X[-1], y[-1]), textcoords="offset points", xytext=(6, 0),
                 fontsize=7, color=col, va="center")
axt.axhline(1.0, color="grey", ls="--", lw=1.1)
axt.set_xticks(X)
axt.set_xticklabels([label(c) for c in conds], fontsize=10)
axt.set_ylabel(r"$\delta K$ / normal mean", fontsize=11)
axt.set_xlim(-0.2, len(conds) - 0.45)
axt.grid(alpha=.3)
pairs = [f"{label(conds[i])[:9]}->{label(conds[i+1])[:9]} {spearmanr([np.mean(groups[conds[i]][c]) for c in cases], [np.mean(groups[conds[i+1]][c]) for c in cases])[0]:+.2f}"
         for i in range(len(conds) - 1)]
span = spearmanr([np.mean(groups[conds[0]][c]) for c in cases],
                 [np.mean(groups[conds[-1]][c]) for c in cases])[0]
axt.set_title("one line per case, normalised within each condition   |   "
              "adjacent rank rho: " + ",  ".join(pairs) +
              f"   |   full span {span:+.2f}", fontsize=10.5)

fig.suptitle(f"CAMP robustness — {MODE}, trim 0, connectivity + fan gate", fontsize=14)
fig.tight_layout()
out = f"tools/camp_robust_{MODE}{TAG}.png"
fig.savefig(out, dpi=135)
print(f"  -> {out}")
for i in range(len(conds) - 1):
    a = [np.mean(groups[conds[i]][c]) for c in cases]
    b = [np.mean(groups[conds[i + 1]][c]) for c in cases]
    print(f"  {label(conds[i]):>16s} -> {label(conds[i+1]):<16s} rank rho "
          f"{spearmanr(a, b)[0]:+.2f}")
print(f"  full span {label(conds[0])} -> {label(conds[-1])}: rho {span:+.2f}")
