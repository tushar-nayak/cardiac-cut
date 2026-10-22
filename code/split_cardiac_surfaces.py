#!/usr/bin/env python3
"""Automatically split closed LV myocardial STLs into endocardial/epicardial walls.

The splitter uses the closed mesh's normal field to find coherent inner and outer
regions. It estimates an interior point per mesh, then extends the two regions to a
shared, geodesic Voronoi seam. No per-file labels or hand-placed planes are used.

This is an exploratory geometry-driven method. Each output is accompanied by QC
measurements; the CSV should be reviewed before interpreting CAMP results.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import trimesh

# The geometry lives in vascular-cut/heart_split.py, shared with the web app's
# Cardiac section, so the batch run and the interactive cut are the same method.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "vascular-cut"))
from heart_split import (  # noqa: E402
    _basal_plane, _boundary_signature, _choose_center, _clip_below, _dominant_region,
    _extract_surface, _face_components, _fill_small_holes, _geodesic_partition,
    _outer_lip_mask, _surface_qc,
)


def split_one(path: Path, output_dir: Path, top_cut: str = "plane",
              plane_margin_mm: float = 3.0, fill_holes_mm: float = 0.0) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    mesh = trimesh.load(str(path), force="mesh", process=True)
    source_components = sorted(mesh.split(only_watertight=False),
                               key=lambda component: len(component.faces), reverse=True)
    if not source_components:
        raise ValueError("mesh has no connected components")
    main = source_components[0]
    input_faces = len(mesh.faces)
    # The primary closed shell is split. Detached source islands are preserved as
    # separate files and counted in the QC record instead of being silently dropped.
    for i, component in enumerate(source_components[1:], 1):
        component.export(output_dir / f"{path.stem}_satellite_{i:02d}.stl")

    (score, seed_fraction, cut_cost, center, face_score,
     inner_seed, outer_seed, inner_coherence, outer_coherence) = _choose_center(main)
    inner_mask = _geodesic_partition(main, inner_seed, outer_seed)
    inner = _extract_surface(main, inner_mask)
    outer = _extract_surface(main, ~inner_mask)
    inner_qc = _surface_qc(inner)
    outer_qc = _surface_qc(outer)
    inner_boundary = _boundary_signature(inner)
    outer_boundary = _boundary_signature(outer)
    shared_boundary = (
        inner_boundary.shape == outer_boundary.shape
        and np.array_equal(inner_boundary, outer_boundary)
    )

    ratio = inner.area / (inner.area + outer.area)
    # Topology checks above use the untrimmed partition so the shared seam is
    # verified; the exported surfaces then have the basal top removed.
    outer_untrimmed = outer
    plane_pieces = 0
    plane_separates = False
    if top_cut == "plane":
        origin, normal = _basal_plane(main, outer, plane_margin_mm)
        pieces = sorted(_clip_below(main, origin, normal).split(only_watertight=False),
                        key=lambda piece: piece.area, reverse=True)
        total = sum(piece.area for piece in pieces)
        plane_pieces = sum(piece.area >= 0.02 * total for piece in pieces)
        plane_separates = plane_pieces == 2
        if plane_separates:
            # With the rim gone, the walls only touch through tunnels; two pieces
            # means a clean split. The outer wall's normals point away from the
            # long axis, the inner wall's toward it.
            def outward(piece):
                offset = piece.triangles_center - origin
                radial = offset - np.outer(offset @ normal, normal)
                return np.average(np.einsum("ij,ij->i", piece.face_normals, radial),
                                  weights=piece.area_faces)
            inner, outer = sorted(pieces[:2], key=outward)
        else:
            # Walls still joined below the plane: clip the geodesic split instead.
            inner = max(_clip_below(inner, origin, normal).split(only_watertight=False),
                        key=lambda piece: piece.area)
            outer = max(_clip_below(outer, origin, normal).split(only_watertight=False),
                        key=lambda piece: piece.area)
    elif top_cut == "angle":
        outer = _extract_surface(outer, ~_outer_lip_mask(main, outer))
    inner_holes = outer_holes = ([], [])
    if top_cut != "none":
        # With fill_holes_mm 0 nothing is patched and every hole is reported open.
        inner, *inner_holes = _fill_small_holes(inner, fill_holes_mm)
        outer, *outer_holes = _fill_small_holes(outer, fill_holes_mm)
    trimmed_qc = _surface_qc(outer)
    inner_trimmed_qc = _surface_qc(inner)
    lip_fraction = 1.0 - outer.area / outer_untrimmed.area

    output_dir.mkdir(parents=True, exist_ok=True)
    inner_path = output_dir / f"{path.stem}_inner.stl"
    outer_path = output_dir / f"{path.stem}_outer.stl"
    outer_untrimmed_path = output_dir / f"{path.stem}_outer_untrimmed.stl"
    inner.export(inner_path)
    # Keep the full geometric outer wall for the interactive basal-rim review.
    # The historical heuristic-trimmed output remains available for comparison.
    outer_untrimmed.export(outer_untrimmed_path)
    outer.export(outer_path)

    accepted = (
        inner_qc["components"] == 1 and outer_qc["components"] == 1
        and inner_qc["boundary_is_simple_loops"] and outer_qc["boundary_is_simple_loops"]
        and inner_qc["boundary_loops"] == outer_qc["boundary_loops"]
        and inner_qc["boundary_edges"] == outer_qc["boundary_edges"]
        and shared_boundary
        and inner_qc["boundary_loops"] >= 1 and outer_qc["boundary_loops"] >= 1
        and 0.12 <= ratio <= 0.55
        and min(inner_coherence, outer_coherence) >= 0.70
        and trimmed_qc["components"] == 1
        and inner_trimmed_qc["components"] == 1
        and (top_cut != "plane" or lip_fraction <= 0.40)
        # Holes left open are accepted (the default, as CAMP's reference method
        # does); when patching is asked for, every hole must have been patched.
        and (top_cut == "none" or fill_holes_mm <= 0
             or not (inner_holes[1] or outer_holes[1]))
        and (top_cut != "angle" or lip_fraction <= 0.25)
    )
    return {
        "file": path.name,
        "status": "auto_split_qc_pass" if accepted else "auto_split_qc_review",
        "input_faces": input_faces,
        "main_component_faces": len(main.faces),
        "source_components": len(source_components),
        "satellite_faces": input_faces - len(main.faces),
        "estimated_inner_center_x": float(center[0]),
        "estimated_inner_center_y": float(center[1]),
        "estimated_inner_center_z": float(center[2]),
        "seed_inner_area_fraction": float(seed_fraction),
        "seed_cut_cost": float(cut_cost),
        "inner_seed_coherence": float(inner_coherence),
        "outer_seed_coherence": float(outer_coherence),
        "inner_area_fraction": float(ratio),
        "inner_faces": inner_qc["faces"],
        "outer_faces": outer_qc["faces"],
        "top_cut": top_cut,
        "plane_margin_mm": plane_margin_mm if top_cut == "plane" else "",
        "plane_pieces": plane_pieces,
        "plane_separates_walls": plane_separates,
        "outer_lip_area_fraction": float(lip_fraction),
        "inner_trimmed_faces": inner_trimmed_qc["faces"],
        "inner_trimmed_boundary_loops": inner_trimmed_qc["boundary_loops"],
        "inner_handles": inner_trimmed_qc["handles"],
        "inner_holes_filled": len(inner_holes[0]),
        "outer_holes_filled": len(outer_holes[0]),
        "max_filled_hole_perimeter_mm": max(inner_holes[0] + outer_holes[0], default=0.0),
        "holes_left_open": len(inner_holes[1]) + len(outer_holes[1]),
        "outer_handles": trimmed_qc["handles"],
        "outer_trimmed_faces": trimmed_qc["faces"],
        "outer_trimmed_components": trimmed_qc["components"],
        "outer_trimmed_boundary_loops": trimmed_qc["boundary_loops"],
        "inner_boundary_edges": inner_qc["boundary_edges"],
        "outer_boundary_edges": outer_qc["boundary_edges"],
        "inner_boundary_loops": inner_qc["boundary_loops"],
        "outer_boundary_loops": outer_qc["boundary_loops"],
        "boundary_is_simple_loops": inner_qc["boundary_is_simple_loops"] and outer_qc["boundary_is_simple_loops"],
        "shared_boundary_edges_match": shared_boundary,
        "inner_components": inner_qc["components"],
        "outer_components": outer_qc["components"],
        "inner_stl": inner_path.name,
        "outer_stl": outer_path.name,
        "outer_untrimmed_stl": outer_untrimmed_path.name,
        "error": "",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--dataset", type=Path, default=root / "dataset")
    parser.add_argument("--output", type=Path, default=root / "results" / "automatic_surface_split")
    parser.add_argument("--report", type=Path, default=root / "results" / "automatic_surface_split.csv")
    parser.add_argument("--top-cut", choices=("plane", "angle", "none"), default="plane",
                        help="basal top removal: rim-fitted plane (default), outer-rim "
                             "normal-angle trim, or none")
    parser.add_argument("--plane-margin", type=float, default=3.0,
                        help="mm to move the plane toward the apex beyond the rim (default: 3)")
    parser.add_argument("--fill-holes", type=float, default=0.0,
                        help="patch holes up to this perimeter in mm after the top cut; "
                             "0 (default) leaves them open, as CAMP's reference method does")
    parser.add_argument("--keep-outer-lip", action="store_true",
                        help="same as --top-cut none (untrimmed surfaces for manual review)")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    paths = sorted(p for p in args.dataset.glob("*.stl") if not p.name.startswith("._"))
    rows = []
    for index, path in enumerate(paths, 1):
        print(f"[{index}/{len(paths)}] {path.name}", flush=True)
        try:
            row = split_one(path, args.output,
                            top_cut="none" if args.keep_outer_lip else args.top_cut,
                            plane_margin_mm=args.plane_margin,
                            fill_holes_mm=args.fill_holes)
        except Exception as exc:
            row = {"file": path.name, "status": "auto_split_error", "error": str(exc)}
        rows.append(row)
        print(f"  {row['status']}: {row.get('inner_area_fraction', float('nan')):.3f} inner area; "
              f"loops {row.get('inner_boundary_loops', '?')}/{row.get('outer_boundary_loops', '?')}; "
              f"outer area removed {row.get('outer_lip_area_fraction', float('nan')):.1%}; "
              f"plane pieces {row.get('plane_pieces', '?')}; "
              f"handles inner/outer {row.get('inner_handles', '?')}/{row.get('outer_handles', '?')}; "
              f"holes filled {row.get('inner_holes_filled', '?')}/{row.get('outer_holes_filled', '?')}, "
              f"open {row.get('holes_left_open', '?')}", flush=True)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    with args.report.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nProcessed {len(rows)} files. Report: {args.report}")
    print(f"QC pass: {sum(r['status'] == 'auto_split_qc_pass' for r in rows)}; "
          f"review: {sum(r['status'] != 'auto_split_qc_pass' for r in rows)}")


if __name__ == "__main__":
    main()
