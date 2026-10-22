# Outputs on the share

Heavy run roots from this project live on the lab share, under
`$PROJECTS_SHARE/branches/`, not here. This file lists them by date. Add to it by
appending a dated line at the bottom; never rewrite it. The share also holds a one-time snapshot of
this folder from 2026-09-14, which is not kept in sync.

| date | run root under `branches/` | what wrote it |
|---|---|---|
| 2026-09-25 to 28 | `_cache/check_privacy/runs/<stamp>/`, `_cache/check_privacy/fixtures/` | `tests/web/check_privacy.py`: server log, store, screenshots per run; the 637k-face race mesh `wfPRIVE_KY160_sub2.stl` (30 MB) |
| 2026-09-24 to 28 | `_cache/check_upload/runs/<stamp>/`, `_cache/check_upload/fixtures/` | `tests/web/check_upload.py`: server log, store, downloads, screenshots per run |
| 2026-09-24 | `_cache/face-scaling/meshes_20260924/` | KY160 and JK5 subdivided to x1, x4, x16, x64 (up to 127 MB), for load and cut timings |
| 2026-09-24 | `_cache/freeze_probe/`, `_cache/limits_live/`, `_cache/http-route/` | one-off probes of the upload freeze fix and the upload limits |
| 2026-09-28 | `_cache/preview/<stamp>/` | a public-mode instance started for looking at the site: server log, and the store for anything uploaded to it |
| 2026-09-28 | `_cache/ui_redesign/before_<stamp>/`, `_cache/ui_redesign/after_<stamp>/` | screenshots of the web page in each state, before and after the UI redesign; `compare.html` in the latest after folder sets them side by side |
| 2026-10-02 | `ZF_Example/zones_20261002_165211/`, `ZF_Example/zones_20261002_165404/` | KK16_6 Ishimaru zones 0-3 isolated from the app's zone labelling: one STL per zone, the labelled mesh (.vtp), a preview and summary.json. 165404 is the clean _2 mesh with the original centerline (all zones); 165211 the _2 centerline, which has no zone 2 |
| 2026-10-02 | `ZF_Example/zones_20261002_170822/` | the same KK16_6 zones 0-3 (clean _2 mesh, original centerline) with each zone one piece: a 6-face stranded piece of zone 2 joined zone 1; supersedes 165404 |
