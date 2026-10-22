#!/usr/bin/env python3
"""Summarize within-scan phase spread and between-scan δK variation."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def _number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "results" / "cardiac_delta_k" / "summary.csv")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "cardiac_delta_k")
    args = parser.parse_args()

    with args.input.open(newline="") as f:
        rows = list(csv.DictReader(f))
    groups = defaultdict(list)
    for row in rows:
        value = _number(row.get("k_fluct"))
        if row.get("status") != "ok" or value is None:
            continue
        key = (row.get("surface_type", ""), row.get("patient_label", ""), row.get("scan_label", ""))
        groups[key].append((row.get("phase_label", ""), value))

    group_rows = []
    for (surface, patient, scan), values in sorted(groups.items()):
        delta = [v for _, v in values]
        mean = statistics.mean(delta)
        sd = statistics.stdev(delta) if len(delta) > 1 else None
        group_rows.append({
            "surface_type": surface,
            "patient_label": patient,
            "scan_label": scan,
            "n_phases": len(delta),
            "phases": ";".join(sorted(phase for phase, _ in values)),
            "mean_delta_k": mean,
            "sd_delta_k": sd if sd is not None else "",
            "range_delta_k": max(delta) - min(delta),
            "cv_percent": 100.0 * sd / abs(mean) if sd is not None and mean else "",
        })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    group_path = args.output_dir / "within_scan_delta_k.csv"
    with group_path.open("w", newline="") as f:
        columns = ["surface_type", "patient_label", "scan_label", "n_phases", "phases",
                   "mean_delta_k", "sd_delta_k", "range_delta_k", "cv_percent"]
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(group_rows)

    summary_rows = []
    for surface in ("inner", "outer"):
        complete = [r for r in group_rows if r["surface_type"] == surface and int(r["n_phases"]) > 1]
        means = [float(r["mean_delta_k"]) for r in complete]
        ss_within = sum((int(r["n_phases"]) - 1) * float(r["sd_delta_k"]) ** 2 for r in complete)
        df_within = sum(int(r["n_phases"]) - 1 for r in complete)
        pooled_within = math.sqrt(ss_within / df_within) if df_within else ""
        between_sd = statistics.stdev(means) if len(means) > 1 else ""
        summary_rows.append({
            "surface_type": surface,
            "n_source_meshes": sum(r["surface_type"] == surface for r in rows),
            "n_scan_groups": sum(r["surface_type"] == surface for r in group_rows),
            "n_groups_with_multiple_phases": len(complete),
            "pooled_within_scan_sd": pooled_within,
            "mean_within_scan_sd": statistics.mean(float(r["sd_delta_k"]) for r in complete) if complete else "",
            "mean_within_scan_range": statistics.mean(float(r["range_delta_k"]) for r in complete) if complete else "",
            "between_scan_mean_sd": between_sd,
            "within_to_between_sd_ratio": pooled_within / between_sd if pooled_within != "" and between_sd else "",
        })
    summary_path = args.output_dir / "within_between_summary.csv"
    with summary_path.open("w", newline="") as f:
        columns = list(summary_rows[0]) if summary_rows else []
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote {group_path} and {summary_path}")


if __name__ == "__main__":
    main()
