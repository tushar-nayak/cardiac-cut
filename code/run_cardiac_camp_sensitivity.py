#!/usr/bin/env python3
"""Compare CAMP remesh resolutions on a small repeated-phase scan group.

Only target face count changes. CAMP fit thresholds, trim fraction, and the
saved automatic split surfaces remain fixed. The comparison is exploratory.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


DEFAULT_STEMS = ("CG2_1_2_MY", "CG2_1_7_MY", "CG2_1_10_MY")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stems", nargs="+", default=DEFAULT_STEMS,
                        help="filename stems without .stl (default: three CG2 phases)")
    parser.add_argument("--targets", nargs="+", type=int, default=[25_000, 100_000],
                        help="alternative remesh targets; 50000 baseline already exists")
    parser.add_argument("--split-dir", type=Path, default=root / "results" / "automatic_surface_split")
    parser.add_argument("--split-report", type=Path, default=root / "results" / "automatic_surface_split.csv")
    parser.add_argument("--baseline", type=Path, default=root / "results" / "cardiac_delta_k" / "summary.csv")
    parser.add_argument("--output", type=Path, default=root / "results" / "cardiac_delta_k_sensitivity")
    parser.add_argument("--trim-frac", type=float, default=0.0025)
    parser.add_argument("--camp-source", type=Path,
                        default=root / "camp" / "camp_updated_2026-08-26" / "campsrc.py",
                        help="CAMP source file under camp/ (default: canonical Aug 26 version)")
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_cardiac_delta_k import run_surface, use_camp
    use_camp(args.camp_source)

    split_rows = {Path(row["file"]).stem: row for row in read_csv(args.split_report)}
    baseline_rows = read_csv(args.baseline)
    baseline = {row["stem"]: row for row in baseline_rows if row.get("status") == "ok"}
    run_rows = []
    compare_rows = []
    for source_stem in args.stems:
        if source_stem not in split_rows:
            raise ValueError(f"No split record for {source_stem}")
        split = split_rows[source_stem]
        for side, column in (("inner", "inner_stl"), ("outer", "outer_stl")):
            split_path = args.split_dir / split[column]
            surface_stem = f"{source_stem}_{side}"
            base = baseline.get(surface_stem)
            if base is None:
                raise ValueError(f"No successful 50000-face baseline for {surface_stem}")
            for target in args.targets:
                out_dir = args.output / f"target_{target}" / surface_stem
                print(f"[{source_stem} {side} @ {target} faces]", flush=True)
                result = run_surface(split_path, out_dir, target, args.trim_frac)
                result.update({
                    "source_file": split["file"],
                    "surface_type": side,
                    "comparison_group": "50k baseline",
                    "baseline_k_fluct": base.get("k_fluct", ""),
                    "baseline_n_patches": base.get("n_patches", ""),
                    "baseline_pct_failed_quality": base.get("pct_failed_quality", ""),
                    "baseline_pct_at_min_size": base.get("pct_at_min_size", ""),
                })
                run_rows.append(result)
                if result.get("status") == "ok":
                    b_k = float(base["k_fluct"])
                    a_k = float(result["k_fluct"])
                    compare_rows.append({
                        "source_file": split["file"],
                        "surface_type": side,
                        "target_faces": target,
                        "status": result["status"],
                        "baseline_50k_k_fluct": b_k,
                        "candidate_k_fluct": a_k,
                        "k_fluct_change_percent": 100 * (a_k - b_k) / abs(b_k) if b_k else "",
                        "baseline_n_patches": int(base["n_patches"]),
                        "candidate_n_patches": int(result["n_patches"]),
                        "n_patches_change": int(result["n_patches"]) - int(base["n_patches"]),
                        "baseline_pct_failed_quality": float(base["pct_failed_quality"]),
                        "candidate_pct_failed_quality": float(result["pct_failed_quality"]),
                        "failed_quality_change_pp": float(result["pct_failed_quality"]) - float(base["pct_failed_quality"]),
                        "baseline_pct_at_min_size": float(base["pct_at_min_size"]),
                        "candidate_pct_at_min_size": float(result["pct_at_min_size"]),
                        "at_min_size_change_pp": float(result["pct_at_min_size"]) - float(base["pct_at_min_size"]),
                        "baseline_target_faces": 50_000,
                        "trim_frac": args.trim_frac,
                        "R2_threshold": 0.90,
                        "eta_threshold": 0.075,
                    })

    write_csv(args.output / "run_summary.csv", run_rows)
    write_csv(args.output / "comparison_to_50000.csv", compare_rows)
    print(f"Wrote {args.output / 'run_summary.csv'} and {args.output / 'comparison_to_50000.csv'}")


if __name__ == "__main__":
    main()
