# CAMP Toolkit

Run **CAMP** (Convergent Adaptive Monge Patching) on STL surface meshes — e.g. the
outputs of the vascular-cut tool — and get a per-face patch partition plus the
integrated–Gaussian–curvature shape metric for each mesh.

Everything is driven from one notebook: point it at a file or folder of STLs, set a
few parameters, and Run All.

```
STL  →  weld duplicated corners  →  uniform ACVD remesh  →  CAMP  →  outputs
```

The uniform remesh is the **only** preprocessing step. Nothing else is done to your
geometry.

---

## 1. Install

Pick one:

```bash
# pip
pip install -r requirements.txt

# or conda
conda env create -f environment.yml
conda activate camp
```

## 2. Run

```bash
jupyter notebook run_camp.ipynb      # (or open it in VS Code / JupyterLab)
```

In the notebook, edit the single **① Settings** cell:

- `INPUT` — a `.stl` file **or** a folder of `.stl` files
- `OUT_DIR` — where results are written
- `TARGET_FACES` — remesh density (~50k plain aortas, ~100k branched)
- `TRIM_FRAC` — quantile trim for the `k_fluct` metric (0.01 = production)
- `CAMP_PARAMS` — CAMP quality thresholds (defaults are the project operating point)

Then **Kernel ▸ Restart & Run All**.

### Or call it from Python directly

```python
from camp_toolkit import run_camp

summary = run_camp(
    inputs="my_stls/",           # file or folder
    out_dir="results",
    target_faces=100_000,
    trim_frac=0.01,
    params=dict(R2_thresh=0.9, eta_thresh=0.075, min_patch_size=18, seed_patches=100),
)
print(summary)                    # one row per mesh; also written to results/summary.csv
```

## 3. Outputs

Per mesh, in `OUT_DIR/<mesh_name>/`:

| file | contents |
|---|---|
| `<name>_facepatches.csv`   | `face_id, patch_id` for every face of the remeshed mesh |
| `<name>_patch_metrics.csv` | per patch: `avg_kg, area_mm2, integrated_K, eta, R2, kept` |
| `<name>_camp.npz`          | `V, F, face_labels` + per-patch arrays (for re-plotting) |
| `<name>_patches.html`      | interactive 3-D viewer with an **adjustable integrated-K color axis** |
| `<name>_integratedK.png`   | static curvature image |

Plus `OUT_DIR/summary.csv` — one row per mesh with `n_patches`, `k_fluct`,
`surface_area_mm2`, `euler_in`/`euler_remesh`, `pct_at_min_size`, `runtime_s`, and the
run parameters.

**`k_fluct`** is the headline shape descriptor: the area-weighted variance of
integrated Gaussian curvature across patches (higher = more shape heterogeneity).

## 4. Visualize / re-color

Each `*_patches.html` opens in any browser with built-in controls — **color min/max
sliders** (drag to set the integrated-K axis), a **colorscale dropdown**, and an
**integrated-K ↔ patches** toggle. No server or install needed; the file is shareable
as-is.

Re-open or re-color any result from Python (or in the notebook):

```python
from camp_toolkit import view

view("results/KY160_1_open_1").show()                    # integrated K, adjustable axis
view("results/KY160_1_open_1", cmin=-5, cmax=40).show()  # fixed axis (compare meshes)
view("results/KY160_1_open_1", metric="eta").show()      # or avg_kg / R2 / integrated_K
```

## Parameters, briefly

| knob | default | effect |
|---|---|---|
| `TARGET_FACES` | 100k | remesh density; more faces = finer subdivision near ostia/arch, slower. **Keep fixed across a cohort.** |
| `TRIM_FRAC` | 0.01 | trims extreme integrated-K patches before `k_fluct`; metric only, not the partition |
| `R2_thresh` | 0.9 | min variance a patch's quadratic fit must explain; lower → coarser |
| `eta_thresh` | 0.075 | max curvature variation within a patch (≈ radius ≤ 0.39× local R); lower → finer |
| `min_patch_size` | 18 | vertex floor; patches never split below this |
| `seed_patches` | 100 | initial KMeans seeds; final count is data-driven |

## What's in here

```
CAMP_Toolkit/
├── run_camp.ipynb            # ← start here
├── requirements.txt
├── environment.yml
├── README.md
└── camp_toolkit/
    ├── __init__.py           # run_camp, run_stl, weld_mesh, ...
    ├── pipeline.py           # weld → remesh → CAMP → write outputs
    ├── remesh.py             # uniform ACVD remeshing
    ├── visualize.py          # interactive HTML + static PNG
    └── campsrc.py            # the CAMP algorithm (authoritative source)
```

## Notes

- **Input meshes** should be single-component open surfaces (as produced by
  vascular-cut). STL "triangle soup" is handled automatically (welding).
- The static PNG needs a working OpenGL context; if one isn't available it's skipped
  and the interactive HTML is still written.
- Runtime is roughly 1–3 minutes per mesh at 100k faces.
