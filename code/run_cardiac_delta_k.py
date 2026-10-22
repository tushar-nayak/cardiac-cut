#!/usr/bin/env python3
"""Run the automatic myocardial split and CAMP δK analysis on every dataset STL."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


NAME_RE = re.compile(r"^(?P<patient>CG\d+)_(?P<scan>\d+)_(?P<phase>\d+)_MY$")
ROOT = Path(__file__).resolve().parents[1]
CAMP_DIR = ROOT / "camp"
CANONICAL_CAMP = CAMP_DIR / "camp_updated_2026-08-26" / "campsrc.py"
# delta K is computed by vascular-cut's shape_metrics.delta_k -- the same function the
# web app's Shape section runs -- so a batch number and an interactive number for the
# same surface and settings are the same number.
sys.path.insert(0, str(ROOT / "vascular-cut"))
# The δK operating point vascular-cut measured: at trim 0.0025 two remeshes of one
# surface differ by 1.6x typically (CURVATURE.md). The app defaults to it too.
DEFAULT_TRIM = 0.0025


def use_camp(camp_source: Path = CANONICAL_CAMP):
    """Point shape_metrics at a campsrc.py under camp/ (it refuses anything else)."""
    camp_source = Path(camp_source).resolve()
    if not camp_source.is_relative_to(CAMP_DIR.resolve()):
        raise ValueError(f"CAMP source must be under {CAMP_DIR}: {camp_source}")
    os.environ["CAMP_SRC"] = str(camp_source)
    import shape_metrics as SM
    SM.CAMP_SRC = str(camp_source)
    SM._campsrc = None
    SM.load_campsrc()
    return SM


def run_surface(stl: Path, out_dir: Path, target_faces: int, trim_frac: float,
                save_html: bool = False) -> dict:
    """δK for one surface via shape_metrics.delta_k, written in the layout the viewer
    and the quality audit read: <stem>/<stem>_camp.npz, _patch_metrics.csv and
    _facepatches.csv. Returns one summary row; a failure is a row, not an exception.
    """
    import pandas as pd
    import trimesh
    import shape_metrics as SM
    stl = Path(stl)
    stem = stl.stem
    mesh_out = Path(out_dir) / stem
    mesh_out.mkdir(parents=True, exist_ok=True)
    row = {"stem": stem, "status": "ok", "error": ""}
    t0 = time.time()
    try:
        welded = SM.weld(trimesh.load(str(stl), process=False))
        row["faces_in"] = int(len(welded.faces))
        row["euler_in"] = int(welded.euler_number)
        defects = SM.surface_defects(welded)
        # delta_k checks remeshing did not change the topology, and that the remesh
        # reached the target face count.
        r = SM.delta_k(welded, target_faces=target_faces, trim_frac=trim_frac,
                       make_figure=False, expect_euler=row["euler_in"])
        mesh = r["mesh"]
        V, F = np.asarray(mesh.vertices, float), np.asarray(mesh.faces, np.int64)
        labels = np.asarray(r["face_labels"], np.int64)
        n = len(r["integrated_K"])
        row.update({
            "faces_remesh": r["faces"], "euler_remesh": r["euler_number"],
            "surface_area_mm2": r["A_mm2"], "n_faces": int(len(F)),
            "n_patches": n, "n_kept": r["n_kept"], "n_trimmed": int(n - r["n_kept"]),
            "k_fluct": r["k_fluct"], "target_faces": int(target_faces),
            "trim_frac": float(trim_frac),
            "pct_at_min_size": r["pct_at_min_size"],
            "pct_failed_quality": r["pct_failed_quality"],
            "folds": defects["folds"], "campsrc": SM.campsrc_origin(),
        })
        pd.DataFrame({"face_id": np.arange(len(F)), "patch_id": labels}) \
            .to_csv(mesh_out / f"{stem}_facepatches.csv", index=False)
        pd.DataFrame({"patch_id": np.arange(n), "avg_kg": r["avg_kg"],
                      "area_mm2": r["patch_area"], "integrated_K": r["integrated_K"],
                      "eta": r["eta_patch"], "R2": r["R2_patch"], "kept": r["keep"]}) \
            .to_csv(mesh_out / f"{stem}_patch_metrics.csv", index=False)
        np.savez_compressed(mesh_out / f"{stem}_camp.npz", V=V, F=F,
                            face_labels=labels.astype(np.int32), avg_kg=r["avg_kg"],
                            patch_area=r["patch_area"],
                            K_integrated_patch=r["integrated_K"],
                            eta_patch=r["eta_patch"], R2_patch=r["R2_patch"],
                            keep=r["keep"])
        if save_html:
            sys.path.insert(0, str(CAMP_DIR / "CAMP_Toolkit"))
            from camp_toolkit.visualize import save_patch_html
            row["html"] = save_patch_html(
                mesh_out / f"{stem}_patches.html", V, F, labels, r["integrated_K"],
                title=f"{stem}  ·  {n} patches  ·  k_fluct={r['k_fluct']:.4g}")
    except Exception as e:                                          # noqa: BLE001
        row.update({"status": "error", "error": f"{type(e).__name__}: {e}"})
    row["runtime_s"] = round(time.time() - t0, 1)
    print(f"  [{stem}] {row['status']}: k_fluct={row.get('k_fluct', float('nan')):.4g}, "
          f"{row['runtime_s']}s {row['error']}", flush=True)
    return row


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=root / "dataset")
    parser.add_argument("--results", type=Path, default=root / "results" / "cardiac_delta_k")
    parser.add_argument("--target-faces", type=int, default=50_000,
                        help="fixed CAMP remesh size per surface (default: 50000)")
    parser.add_argument("--trim-frac", type=float, default=DEFAULT_TRIM,
                        help=f"CAMP integrated-K trim fraction (default: {DEFAULT_TRIM:g}, "
                             f"the web app's and vascular-cut's operating point)")
    parser.add_argument("--save-html", action="store_true",
                        help="save an interactive CAMP viewer for every surface")
    parser.add_argument("--split-only", action="store_true",
                        help="run the automatic split and QC, without CAMP")
    parser.add_argument("--reuse-split", action="store_true",
                        help="use the existing automatic split outputs")
    parser.add_argument("--reviewed-splits", type=Path,
                        help="run CAMP on wall pairs exported by the app's Cardiac section "
                             "(results/cardiac_walls)")
    parser.add_argument("--camp-source", type=Path,
                        default=CANONICAL_CAMP,
                        help="CAMP source file to use (default: canonical Aug 26 version)")
    args = parser.parse_args()

    split_rows = []
    jobs = []
    if args.reviewed_splits:
        for folder in sorted(p for p in args.reviewed_splits.iterdir() if p.is_dir()):
            inner = folder / f"{folder.name}_inner.stl"
            outer = folder / f"{folder.name}_outer.stl"
            if inner.exists() and outer.exists():
                jobs.extend((folder.name, side, path) for side, path in
                            (("inner", inner), ("outer", outer)))
        if not jobs:
            raise FileNotFoundError(f"No reviewed inner/outer surface pairs found under {args.reviewed_splits}")
    else:
        split_dir = root / "results" / "automatic_surface_split"
        split_report = root / "results" / "automatic_surface_split.csv"
        if not args.reuse_split:
            subprocess.run([
                sys.executable, str(Path(__file__).with_name("split_cardiac_surfaces.py")),
                "--dataset", str(args.dataset), "--output", str(split_dir),
                "--report", str(split_report),
            ], check=True)
        if not split_report.exists():
            raise FileNotFoundError(f"Missing split report: {split_report}")
        with split_report.open(newline="") as f:
            split_rows = list(csv.DictReader(f))
        failures = [r["file"] for r in split_rows if r.get("status") != "auto_split_qc_pass"]
        if failures:
            raise RuntimeError(
                f"Automatic split QC did not pass for {len(failures)} files: "
                + ", ".join(failures)
            )
        if args.split_only:
            print(f"Automatic split and QC completed for {len(split_rows)} files.")
            return
        for row in split_rows:
            stem = Path(row["file"]).stem
            for side, column in (("inner", "inner_stl"), ("outer", "outer_stl")):
                jobs.append((stem, side, split_dir / row[column]))

    use_camp(args.camp_source)

    args.results.mkdir(parents=True, exist_ok=True)
    rows = []
    n_jobs = len(jobs)
    for stem, surface_type, stl in jobs:
        match = NAME_RE.match(stem)
        output_dir = args.results / f"{stem}_{surface_type}"
        print(f"[CAMP {len(rows) + 1}/{n_jobs}] {stl.name}", flush=True)
        result = run_surface(stl, output_dir, args.target_faces, args.trim_frac,
                             save_html=args.save_html)
        result.update({
            "source_file": f"{stem}.stl",
            "patient_label": match.group("patient") if match else "",
            "scan_label": match.group("scan") if match else "",
            "phase_label": match.group("phase") if match else "",
            "surface_type": surface_type,
            "split_inner_area_fraction": "",
        })
        rows.append(result)

    report = args.results / "summary.csv"
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with report.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    ok = sum(row.get("status") == "ok" for row in rows)
    print(f"CAMP complete: {ok}/{len(rows)} succeeded; results in {report}")
    if ok == len(rows):
        subprocess.run([
            sys.executable, str(Path(__file__).with_name("summarize_cardiac_delta_k.py")),
            "--input", str(report), "--output-dir", str(args.results),
        ], check=True)


if __name__ == "__main__":
    main()
