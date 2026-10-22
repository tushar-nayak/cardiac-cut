# Cardiac surface review and CAMP

This repository includes a browser based tool for inspecting a heart surface,
choosing an inner/outer split, reviewing the basal cut, and optionally running
CAMP. Meshes are local inputs: the app does not require a dataset to be stored
in or committed to the repository.

## Standalone app: clone, install, run

### 1. Clone the repository

In a terminal with Git installed:

```bash
git clone https://github.com/tushar-nayak/cardiac-cut.git
cd cardiac-cut
```

### 2. Create the Conda environment

Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
or Mamba if it is not already available, then create and activate the environment
defined by the repository:

```bash
conda env create -f camp/CAMP_Toolkit/environment.yml
conda activate camp
```

If the `camp` environment already exists, activate it instead. The environment
includes the mesh, visualization, and CAMP dependencies used by the app.

### 3. Add a mesh

The repository contains an empty `input/` folder. Copy one or more STL, OBJ, or
PLY surface meshes into it. For example:

```bash
cp /path/to/heart.stl input/
```

Meshes should use millimeter coordinates. The `input/` folder is gitignored, so
your meshes remain local and are not uploaded to GitHub. You can also skip the
copy and use **Upload mesh** in the app; uploads are saved in `input/` and remain
available after restarting it. OBJ and PLY uploads are normalized to STL.

### 4. Start the app

From the repository root, run:

```bash
python code/review_cardiac_split_web.py --port 8767 --no-browser
```

Open [http://localhost:8767](http://localhost:8767) in your browser. The mesh
dropdown scans `input/` at startup. If you add a file while the app is open,
click **Refresh inputs**. In WSL, run the command in the VS Code WSL terminal and
open the localhost link in your Windows browser.

Choose **Center + normals** for the automatic split or **Manual plane** to
separate the surfaces with the plane controls. Inspect both views, then choose
**Save this mesh** to save the reviewed surfaces and provenance record under
`results/manual_surface_split/`. Choose **Save cut + Run CAMP** to also run CAMP
on both surfaces; CAMP results are stored under `results/cardiac_delta_k_gui/`.
All files under `results/` are local outputs and are gitignored.

Stop the server with **Ctrl+C** in its terminal. If port 8767 is occupied, use
another port, for example `--port 8768`, and open that port in the browser.

### Run with a different input folder

By default the app scans the repository's `input/` directory. To scan another
folder instead, pass `--input`:

```bash
python code/review_cardiac_split_web.py --input /path/to/meshes --port 8767 --no-browser
```

Uploads go to the selected input folder too. `--dataset` remains available as a
backwards-compatible alias for `--input`.

## Optional batch analysis

The repository also has command line scripts for mesh inventory, automatic
splitting, and batch CAMP analysis. They are separate from the browser app and
read STL files from a local `dataset/` folder. No such dataset is included.
The geometry based split is a heuristic: topology checks do not prove that the
surfaces match the anatomical endocardial and epicardial walls. Review results
visually and treat them as exploratory until the anatomy and CAMP fit quality
are validated.

For batch runs, the expected input layout is:

```text
dataset/*.stl
    ├── inventory and mesh checks
    ├── automatic geometric split ──> *_inner.stl + *_outer.stl
    ├── CAMP on each surface ───────> patches + per-surface metrics
    └── filename-group summaries ───> phase spread and between-group spread
```

Batch scripts ignore macOS `._*` metadata files. Names matching
`CG<patient>_<scan>_<phase>_MY.stl` are used as metadata in summaries; names do
not determine whether a mesh is processed.

## 1. Inventory the STL files

Create an inventory with mesh format, triangle count, area, bounds, and volume
estimates:

```bash
python code/inventory_dataset.py
```

Add welded topology checks (watertightness, winding, connected components,
Euler number, and a density-based mass estimate):

```bash
python code/inventory_dataset.py --topology
```

The mesh-derived volume and mass are geometry checks. They are not replacements
for voxel-based mass measurements.

## 2. Automatically split each mesh

Run just the split and its quality checks:

```bash
python code/run_cardiac_delta_k.py --split-only
```

For each input, `code/split_cardiac_surfaces.py`:

1. Loads the mesh, finds connected components, and uses the largest component
   as the primary shell. Smaller detached components are retained as separate
   satellite STL files.
2. Builds a fixed grid of candidate interior points around the mesh center of
   mass. For each point it scores triangles by the direction of their normals
   relative to the point.
3. Ranks candidate cuts by separator length and inner/outer area balance. It
   tests the best candidates for coherent connected inner and outer seed
   regions.
4. Extends the two seed regions over the face-adjacency graph using shortest
   geodesic distances between triangle centers. Each face goes to the nearer
   seed; ties go to the inner side. This makes a curved cut on the mesh rather
   than using a manually placed plane.
5. Exports `_inner.stl` and `_outer.stl` and checks connectedness, simple
   boundary loops, matching shared boundary edges, area balance, and seed
   coherence.
6. Cuts off the basal top with a plane (`--top-cut plane`, the default). Rim
   faces are outer faces connected to the seam whose normals point toward the
   base (within 60° of the long axis) or back toward the axis. A plane is
   fitted to them, so a tilted base is followed, and placed at the rim's
   1st-percentile depth plus `--plane-margin` (default 3 mm) toward the apex.
   The closed shell is cut there. With the rim gone, the walls only touch
   through holes in the wall, so if the cut leaves exactly two pieces they
   become the inner and outer surfaces directly (outer = normals pointing away
   from the axis). Otherwise the geodesic split from step 4 is clipped by the
   same plane, which cuts each hole's tube. `outer_lip_area_fraction` records
   the outer area removed (typically 19–33%); above 0.40 goes to review.
   `--top-cut angle` instead drops only the rim faces from the outer surface,
   and `--top-cut none` / `--keep-outer-lip` exports untrimmed surfaces.
   `inner_handles` / `outer_handles` count tunnels or bridges per surface.
7. Leaves holes through the wall **open** (`--fill-holes 0`, the default), as
   CAMP's reference method does: CAMP does not need a closed surface, and a
   patch adds curvature the wall does not have. `holes_left_open` counts them.
   `--fill-holes <mm>` patches holes up to that perimeter instead, for
   comparison; a split then goes to review if any hole was too big to patch.

The geometry lives in `vascular-cut/heart_split.py`, shared with the
interactive app's **Cardiac** section (see `code/README.md`), so a batch split
and an interactive cut are the same method.

The procedure is deterministic for the same triangle mesh and software
versions; it uses no random sampling or per-file manual choices. The exact STL
bytes can still vary across dependency versions or mesh-processing changes.

Split outputs are written to `results/automatic_surface_split/`; the per-file
checks are in `results/automatic_surface_split.csv`. A `auto_split_qc_pass`
means the topology and heuristic thresholds passed, not that anatomy was
validated.

## 3. Run CAMP on both surfaces

Run the full workflow. By default this reruns the split, then runs CAMP on each
inner and outer surface:

```bash
python code/run_cardiac_delta_k.py
```

The default CAMP inputs are remeshed to a fixed target of 50,000 faces per
surface. The integrated-curvature trim fraction is 0.0025 (trim the lowest and
highest 0.25% of patch integrated-curvature values when computing `k_fluct`),
the operating point vascular-cut measured: at it, two remeshes of one surface
differ by 1.6× typically (`vascular-cut/CURVATURE.md`). Results from before this
change used 0.01 and are not comparable. To reuse the saved split files rather
than recomputing them:

```bash
python code/run_cardiac_delta_k.py --reuse-split --target-faces 50000 --trim-frac 0.0025
```

δK is computed by `vascular-cut/shape_metrics.delta_k`, the same function the
web app's Shape δK section runs, so a batch value and an app value for the same
surface and settings are identical. It welds the STL, refuses a remesh that
changes the topology, and loads CAMP only from `camp/`.

Add `--save-html` to the run command to ask CAMP to save its per-surface HTML
visualizations. The number of processed surfaces depends on the STL files in
your local `dataset/` directory (two surfaces per mesh).

### What CAMP writes

Each surface has a folder under `results/cardiac_delta_k/`. CAMP writes:

- `*_camp.npz`: remeshed vertices and faces, face-to-patch labels, and patch
  measurements.
- `*_facepatches.csv`: patch ID assigned to each face.
- `*_patch_metrics.csv`: patch area, mean Gaussian curvature, integrated
  Gaussian curvature, fit `eta`, fit `R2`, and whether the patch was retained
  after trimming.

Recount patch failures by CAMP's default limits and split them by failure type:

```bash
python code/audit_cardiac_camp_quality.py
```

This writes `quality_patch_audit.csv` under the CAMP results folder. The local
3D viewer also offers a green-pass/red-fail quality field using those defaults.

The main `summary.csv` has one row per surface. Its `k_fluct` field is the
workflow's δK measure: the area-weighted variance of patch integrated Gaussian
curvature over retained patches,

```text
δK = sum_p (A_p / A_total) * (K_p - weighted_mean(K))^2
```

where `A_p` is patch area and `K_p` is integrated Gaussian curvature for patch
`p`. CAMP also records `pct_failed_quality` and `pct_at_min_size`. High values
mean many patch fits failed quality criteria or reached CAMP's minimum patch
size; a numeric `k_fluct` is not by itself evidence of a reliable fit.

For each patch, CAMP fits a local quadratic surface (a Monge patch) to the
surface vertices. Its two fit diagnostics are:

- **R²** measures how well that quadratic surface fits the patch point heights
  in CAMP's local coordinate frame: `1 - residual sum of squares / total sum of
  squares`. Higher is better. CAMP's default split target is `R² >= 0.90`.
- **η (eta)** measures relative variation in the Gaussian curvature estimated
  across the patch: `std(K across patch) / abs(K at patch center)`. Lower means
  curvature is more uniform within the patch. CAMP's default target is
  `η <= 0.075`. When center curvature is near zero, η can become very large
  because it divides by a small number.

CAMP subdivides patches while either target is missed and the patch is large
enough to split. A patch can remain in the output while failing one or both
targets if it reaches the minimum size. These are geometric fit and
curvature-uniformity diagnostics, not anatomical confidence scores.

## 4. Summarize phase variation

When all CAMP surfaces finish successfully, the runner writes:

- `within_scan_delta_k.csv`: per filename-group phase count, mean δK, standard
  deviation, range, and coefficient of variation.
- `within_between_summary.csv`: a descriptive pooled within-group spread and
  between-group spread, separately for inner and outer surfaces.

Regenerate these files from an existing `summary.csv` with:

```bash
python code/summarize_cardiac_delta_k.py
```

Grouping uses the patient and scan labels parsed from filenames. It includes
incomplete groups, and reports how many phases are present. These summaries
describe the files in this folder; they do not establish a conserved biological
effect.

## 5. View patches in 3D

An example paired viewer is already available:

[CG10 phase 10: inner and outer CAMP patches side by side](results/cardiac_delta_k/viewer/CG10_1_10_MY_inner_outer.html)

Open it in a browser. Drag each surface to rotate, scroll to zoom, use the
dropdown to color by a CAMP field, and hover for patch ID and value. Each side
uses its own color scale.

Export a paired viewer for another filename stem, without `_inner` or `_outer`:

```bash
python code/view_cardiac_camp.py \
  --case CG10_1_10_MY \
  --side-by-side \
  --export results/cardiac_delta_k/viewer/CG10_1_10_MY_inner_outer.html
```

The script also has a local selector app for all surfaces:

```bash
python code/view_cardiac_camp.py
```

It serves the viewer at `http://localhost:8765`. In environments that prohibit
local network ports, use the standalone HTML export command instead. Viewer
HTML embeds the mesh data and can be large.

For a paired green/red CAMP quality map, export with:

```bash
python code/view_cardiac_camp.py \
  --case CG2_1_2_MY \
  --side-by-side \
  --field quality_status \
  --export results/split_visual_qc/CG2_1_2_MY_quality.html
```

## Small remesh-resolution check

The focused sensitivity run covers the three CG2 phases, holding the split,
CAMP fit thresholds, and trim fixed while changing remesh target:

```bash
python code/run_cardiac_camp_sensitivity.py
```

It compares 25k and 100k remesh targets against the existing 50k baseline. The
current outputs are under `results/cardiac_delta_k_sensitivity/`. In this
sample, 100k lowered median failed-fit rates by about 13 percentage points on
inner surfaces and 10 points on outer surfaces, but increased patch counts by
roughly 2× and shifted median `k_fluct` by about -58% and -39%, respectively.
Thus the fit diagnostics improved while δK remained strongly
resolution-sensitive; this does not justify choosing 100k as a final setting
or rerunning the full dataset.

The per-surface comparison is `comparison_to_50000.csv`; `run_summary.csv`
contains the raw CAMP rows for each alternative target.

## Current diagnostic review

Representative meshes were selected across the observed inner-surface area
fraction range and inspected as 3D partitions and as three short-axis slices
along each mesh's PCA long axis. In the middle slices of these three examples,
the inner candidate traces appear inside the outer candidate contour. Some
basal/apical slices are incomplete or irregular. This is encouraging geometry
evidence for the split, but it is not anatomical validation and does not cover
all meshes. The figures are in `results/split_visual_qc/`:

- `representative_splits.png` — 3D partition and isolated candidate surfaces.
- `split_cross_sections.png` — inner/outer contours at 25%, 50%, and 75% of
  the PCA long axis for representative meshes.
- `camp_fit_quality_CG2.png` — red faces fail at least one default fit target;
  green faces pass both, for the three CG2 phases at the 50k baseline.
- `CG2_1_2_MY_quality.html` — interactive paired 3D fit-quality map.

The quality audit's median passing-patch fraction is about 12% on inner
surfaces and 34% on outer surfaces. Low R² is the most frequent failure mode;
many patches also reach CAMP's minimum-size stop. The 100k target improves
these rates in the CG2 sample, but its δK values differ substantially from
50k, so neither the split nor the metric is ready for a biological conclusion.

## Output map

```text
results/
├── dataset_inventory.csv
├── automatic_surface_split.csv
├── automatic_surface_split/
│   ├── <stem>_inner.stl
│   ├── <stem>_outer.stl
│   └── <stem>_satellite_*.stl        # only if detached components exist
└── cardiac_delta_k/
    ├── summary.csv
    ├── within_scan_delta_k.csv
    ├── within_between_summary.csv
    ├── <stem>_inner/<stem>_inner/    # CAMP outputs
    ├── <stem>_outer/<stem>_outer/
    └── viewer/*.html
```

## Current interpretation limits

- The automatic split uses mesh geometry and face normals; it has not been
  validated against labeled endocardial/epicardial anatomy.
- Matching cut boundaries and connected surfaces are topology checks, not
  anatomical validation.
- Current CAMP outputs have high patch fit-failure and minimum-size rates,
  particularly for inner surfaces. Review those fields in `summary.csv` before
  interpreting δK.
- The current summaries are exploratory and do not yet show that cardiac-cycle
  shape variation is conserved or smaller than between-subject variation.
