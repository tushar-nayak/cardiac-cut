# Cardiac δK workflow

This workflow takes every real STL in `dataset/`, automatically splits each
closed myocardial mesh into inner and outer open surfaces, checks the cut, and
runs CAMP on both surfaces. Filenames are used as metadata only; incomplete and
additional scan groups remain in the run.

## Run

For the automatic split and geometry checks only:

```bash
python code/run_cardiac_delta_k.py --split-only
```

For the full split plus CAMP analysis, create and activate the existing CAMP
environment, then run:

```bash
conda env create -f camp/CAMP_Toolkit/environment.yml
conda activate camp
python code/run_cardiac_delta_k.py
```

The default CAMP remesh target is 50,000 faces per surface, held fixed across
the dataset. To choose a different fixed target or reuse existing split files:

```bash
python code/run_cardiac_delta_k.py --reuse-split --target-faces 50000
```

Add `--save-html` to write interactive CAMP patch viewers. A full run processes
68 surfaces (inner and outer for each of the 34 STLs) and can take a while. CAMP
quality diagnostics (`pct_failed_quality`, `pct_at_min_size`) are retained in
the summary because a returned δK does not by itself mean the patches converged
well on these cardiac surfaces.

## Interactive CAMP viewer

After CAMP has completed, launch the local browser viewer, which lets you switch
among all inner and outer surfaces and CAMP fields:

```bash
.conda-env/bin/python code/view_cardiac_camp.py
```

It opens a browser at `http://localhost:8765` with selectors for every inner and
outer surface and for CAMP patch ID, mean Gaussian curvature, integrated
curvature, eta, R², patch area, and pass/fail status under CAMP's default R²/eta
limits. Drag to rotate, scroll to zoom, hover for patch measurements, and toggle
patch boundary lines. Only the selected mesh is loaded into the browser. The
viewer is read-only and serves files on localhost.

To export one standalone offline HTML view instead, choose a surface name from
the output folders and run:

```bash
.conda-env/bin/python code/view_cardiac_camp.py \
  --case CG10_1_10_MY \
  --side-by-side \
  --export results/cardiac_delta_k/viewer/CG10_1_10_MY_inner_outer.html
```

The paired viewer puts inner and outer surfaces next to each other and uses one
field dropdown for both. Each surface has its own color scale so differences in
scale do not flatten one side. Open the HTML in a browser; patch boundaries
remain visible. The viewer HTML embeds mesh data, so it is large and should be
shared only through approved channels.

For a red/green quality pass overlay, add `--field quality_status` to the export
command. Recount per-patch pass and failure categories across all CAMP outputs:

```bash
python code/audit_cardiac_camp_quality.py
```

For the small resolution sensitivity check on the three CG2 phases:

```bash
python code/run_cardiac_camp_sensitivity.py
```

This keeps CAMP's quality thresholds fixed and compares 25k/100k remesh targets
against the existing 50k results. See the top-level `README.md` for the current
diagnostic findings and limits.

## Automatic split

For each mesh, the script estimates an interior point from its geometry and
scores faces by their normal direction relative to that point. It chooses
coherent inner and outer seed regions, then assigns the remaining faces by
geodesic distance to those regions. The shared cut follows the mesh, rather
than a manually placed plane. No file-specific labels are required.

The split report checks that each result is connected, has simple boundary
loops, and shares the same boundary edges with its paired surface. It also
records the estimated center, surface-area split, and seed coherence. These are
geometry and topology checks; they do not prove anatomical correctness. The
method is an exploratory automatic splitter and its first full-dataset outputs
should be visually inspected before interpreting δK scientifically.

## Outputs

- `results/dataset_inventory.csv` — face count, area, bounds, volume, mass
  estimate, and source topology for every input STL.
- `results/automatic_surface_split/` — paired `_inner.stl` and `_outer.stl`
  surfaces, plus any disconnected source components saved as satellite STLs.
- `results/automatic_surface_split.csv` — automatic split QC for every file.
- `results/cardiac_delta_k/summary.csv` — one CAMP row per source surface when
  the full run completes.
- `results/cardiac_delta_k/within_scan_delta_k.csv` and
  `within_between_summary.csv` — phase spread by filename group and the
  pooled within-scan versus between-scan comparison.

## Interactive cardiac mesh review

Start the dual-view browser app from the WSL terminal:

```bash
.conda-env/bin/python code/review_cardiac_split_web.py \
  --dataset dataset \
  --split-dir results/gui_normals_split \
  --output results/manual_surface_split \
  --camp-results results/cardiac_delta_k_gui \
  --port 8767 --no-browser
```

Open `http://localhost:8767/`. The app lists the dataset meshes and also accepts
STL, OBJ, and PLY uploads. Uploaded meshes are interpreted in millimeters and
stored separately in `results/cardiac_uploads/`; the source dataset is left alone.
Choose automatic center/normal splitting or a manual plane, adjust the basal cut,
and review the inner and outer walls side by side. The review panel reports
connected components, boundary loops, Euler number, and removed area. Its
cross-section tool measures area, perimeter, and approximate diameters on both
surfaces. Saved outputs include the two STLs and a `review.json` provenance and
quality record. `Save cut + Run CAMP` runs CAMP on both saved surfaces.

Both viewers use the same fixed mesh bounds and share rotations, pan, and zoom.
Moving or tilting the cut plane preserves the live camera during both the
immediate plane update and the completed mesh preview. Camera synchronization
also follows the short settling motion after turntable rotation.

The app binds to localhost and needs no desktop window. Stop it with Ctrl+C.

The STL-derived mass estimate uses the paper's density of 1.055 g/cm³ and is a
mesh-based check, not a replacement for the paper's voxel-based measurements.
