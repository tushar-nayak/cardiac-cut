# Tests

Regression suites for the vascular cut tool. They drive the real GUI and the real
cutting pipeline against the geometries in this repository, there are no mocks,
because the failures worth catching are geometric.

## Running

Each suite is a standalone script, run one at a time:

    python tests/test_region_gui.py ky160

**Run them one per process.** Two `QtInteractor` windows in one interpreter
segfault VTK, which is why there is no single pytest entry point.

Suites taking a model argument accept `jk5` or `ky160`. The rest take none.

## What they cover

- `test_ui_fixes`, `test_layout`, `test_timer_resurrect`, panel layout, the
  collapsible sections, and the debounced commit timers
- `test_slider_commit`, `test_slider_gradual`, that a slider changed by keyboard,
  wheel or groove click actually re-cuts, and that region bounds track smoothly
- `test_region_escape`, `test_empty_region`, that a region can always be left, and
  that a region which collapses to nothing is refused rather than exported
- `test_region`, `test_region_gui`, `test_arch_preset`, region planning, presets,
  and the aortic arch's distal-bound search
- `test_region_planes`, that an isolated region's own two boundary cuts can be
  tilted and resized, and that those settings survive a nudge of a region slider
  (which rebuilds the bounds from scratch)
- `test_anatomy_gui`, `test_family_gui`, anatomical naming and the abdominal /
  thoracic declaration
- `test_trim_gui`, `test_trim_fixes`, cutting a vessel before it splits
- `test_newgeo`, `test_features_ab`, `test_new_features`, the cutting core
- `test_angle`, `test_gui_measure`, diameter along the vessel and the arc angle of
  the seal zone, both from the geometry and through the window
- `test_open_region`, `test_multitrim`, one-sided regions, and trimming more than
  one branch inside an isolated one
- `test_shape_job`, delta K over the CAMP partition: the worker process, the arrays
  that colour the patches, the comparability warnings, and both failure paths
- `test_region_station`, that a zone cut by station lands where the anatomy says and
  arrives as a usable surface: winding, no zero-area faces, one piece, half-open
  branch membership
- `test_stub_regions`, that a branch stub is found from the crease at its own mouth --
  bounded, whole, no enclosed pockets, off the aortic wall by a radius test the finder
  never sees -- and that giving each stub its ostium's zone leaves no zone islanded
  inside another. Its `--prove` run needs one cohort case (`KR52_1`) to exercise the
  island guard, and says so rather than passing quietly when the share is not mounted

## Not included here

Four suites validate against a larger private cohort held on lab network shares
(`test_order_dependence`, `test_repair`, `test_repaired_cut`, `test_neck_jk6_2`),
and two read the output of a 40-case sweep over that cohort. They are omitted
because they embed paths to data that is not part of this repository.
- `test_manual_naming`, `test_api_naming_nub`, that declaring an artery by hand
  drives the region presets and displaces the automatic name of the same vessel,
  through the tree directly and through the session the web app uses
- `test_nub_preset`, that the nub-in cut can be applied to a model already open,
  and matches what the batch produces for it
- `test_typed_paths`, that the mesh and centerline boxes take a pasted path
