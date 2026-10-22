#!/usr/bin/env python3
"""Summarize CAMP patch pass/fail counts using the toolkit's default limits."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=root / "results" / "cardiac_delta_k")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--r2-threshold", type=float, default=0.90)
    parser.add_argument("--eta-threshold", type=float, default=0.075)
    args = parser.parse_args()
    summary_path = args.results / "summary.csv"
    with summary_path.open(newline="") as stream:
        surfaces = list(csv.DictReader(stream))

    rows = []
    for surface in surfaces:
        stem = surface["stem"]
        archives = list((args.results / stem).rglob(f"{stem}_camp.npz"))
        if surface.get("status") != "ok" or len(archives) != 1:
            rows.append({"stem": stem, "status": "missing_or_failed_output"})
            continue
        with np.load(archives[0], allow_pickle=False) as data:
            r2 = np.asarray(data["R2_patch"], dtype=float)
            eta = np.asarray(data["eta_patch"], dtype=float)
            keep = np.asarray(data["keep"], dtype=bool)
        valid = np.isfinite(r2) & np.isfinite(eta)
        r2_ok = valid & (r2 >= args.r2_threshold)
        eta_ok = valid & (eta <= args.eta_threshold)
        passed = r2_ok & eta_ok
        r2_only = valid & ~r2_ok & eta_ok
        eta_only = valid & r2_ok & ~eta_ok
        both = valid & ~r2_ok & ~eta_ok
        n = len(r2)
        rows.append({
            "stem": stem,
            "source_file": surface.get("source_file", ""),
            "surface_type": surface.get("surface_type", ""),
            "status": "ok",
            "n_patches": n,
            "n_pass": int(passed.sum()),
            "n_fail_r2_only": int(r2_only.sum()),
            "n_fail_eta_only": int(eta_only.sum()),
            "n_fail_both": int(both.sum()),
            "n_unavailable": int((~valid).sum()),
            "pct_pass": 100 * float(passed.sum()) / n if n else "",
            "pct_fail_r2_only": 100 * float(r2_only.sum()) / n if n else "",
            "pct_fail_eta_only": 100 * float(eta_only.sum()) / n if n else "",
            "pct_fail_both": 100 * float(both.sum()) / n if n else "",
            "pct_kept_for_k_fluct": 100 * float(keep.sum()) / len(keep) if len(keep) else "",
            "reported_pct_at_min_size": surface.get("pct_at_min_size", ""),
            "R2_threshold": args.r2_threshold,
            "eta_threshold": args.eta_threshold,
        })

    output = args.output or args.results / "quality_patch_audit.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output} ({sum(row.get('status') == 'ok' for row in rows)}/{len(rows)} surfaces audited)")


if __name__ == "__main__":
    main()
