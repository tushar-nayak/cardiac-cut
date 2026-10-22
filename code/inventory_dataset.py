#!/usr/bin/env python3
"""Inventory every STL in the cardiac dataset without loading meshes into RAM.

Binary STL files are scanned in chunks. Per-file surface area and signed volume
are geometry checks only; signed volume is meaningful only for a closed,
consistently oriented mesh.
"""

from __future__ import annotations

import argparse
import csv
import re
import struct
from collections import defaultdict
from pathlib import Path

import numpy as np


NAME_RE = re.compile(r"^(?P<patient>CG\d+)_(?P<scan>\d+)_(?P<phase>\d+)_MY$")
RECORD_DTYPE = np.dtype([
    ("normal", "<f4", (3,)),
    ("vertices", "<f4", (3, 3)),
    ("attribute", "<u2"),
])


def inspect_binary_stl(path: Path, chunk_faces: int = 1_000_000) -> dict:
    size = path.stat().st_size
    with path.open("rb") as f:
        header = f.read(84)
    if len(header) < 84:
        raise ValueError("shorter than binary STL header")
    n_faces = struct.unpack_from("<I", header, 80)[0]
    expected_size = 84 + 50 * n_faces
    if expected_size != size:
        raise ValueError(f"not binary STL size (header implies {expected_size} bytes)")

    bounds_min = np.full(3, np.inf)
    bounds_max = np.full(3, -np.inf)
    area = 0.0
    signed_volume = 0.0
    with path.open("rb") as f:
        f.seek(84)
        remaining = n_faces
        while remaining:
            n = min(chunk_faces, remaining)
            records = np.fromfile(f, dtype=RECORD_DTYPE, count=n)
            if len(records) != n:
                raise ValueError("truncated face data")
            v = records["vertices"].astype(np.float64)
            bounds_min = np.minimum(bounds_min, v.min(axis=(0, 1)))
            bounds_max = np.maximum(bounds_max, v.max(axis=(0, 1)))
            cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
            area += 0.5 * np.linalg.norm(cross, axis=1).sum(dtype=np.float64)
            signed_volume += np.einsum("ij,ij->i", v[:, 0], cross).sum(dtype=np.float64) / 6.0
            remaining -= n

    return {
        "format": "binary_stl",
        "faces": n_faces,
        "surface_area_mm2": area,
        "signed_volume_mm3_unverified": signed_volume,
        "bounds_min_x": bounds_min[0], "bounds_min_y": bounds_min[1], "bounds_min_z": bounds_min[2],
        "bounds_max_x": bounds_max[0], "bounds_max_y": bounds_max[1], "bounds_max_z": bounds_max[2],
    }


def inspect_ascii_stl(path: Path, chunk_vertices: int = 300_000) -> dict:
    bounds_min = np.full(3, np.inf)
    bounds_max = np.full(3, -np.inf)
    area = 0.0
    signed_volume = 0.0
    vertices = []
    n_faces = 0

    def consume(batch):
        nonlocal bounds_min, bounds_max, area, signed_volume, n_faces
        v = np.asarray(batch, dtype=np.float64).reshape(-1, 3, 3)
        bounds_min = np.minimum(bounds_min, v.min(axis=(0, 1)))
        bounds_max = np.maximum(bounds_max, v.max(axis=(0, 1)))
        cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
        area += 0.5 * np.linalg.norm(cross, axis=1).sum(dtype=np.float64)
        signed_volume += np.einsum("ij,ij->i", v[:, 0], cross).sum(dtype=np.float64) / 6.0
        n_faces += len(v)

    with path.open("rt", errors="replace") as f:
        for line in f:
            if not line.lstrip().startswith("vertex "):
                continue
            parts = line.split()
            if len(parts) != 4:
                raise ValueError("malformed vertex line")
            vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            if len(vertices) >= chunk_vertices:
                if len(vertices) % 3:
                    raise ValueError("vertex count is not a multiple of three")
                consume(vertices)
                vertices.clear()
    if vertices:
        if len(vertices) % 3:
            raise ValueError("vertex count is not a multiple of three")
        consume(vertices)
    if n_faces == 0:
        raise ValueError("no facets found")
    return {
        "format": "ascii_stl",
        "faces": n_faces,
        "surface_area_mm2": area,
        "signed_volume_mm3_unverified": signed_volume,
        "bounds_min_x": bounds_min[0], "bounds_min_y": bounds_min[1], "bounds_min_z": bounds_min[2],
        "bounds_max_x": bounds_max[0], "bounds_max_y": bounds_max[1], "bounds_max_z": bounds_max[2],
    }


def inspect_topology(path: Path) -> dict:
    """Weld STL triangle soup and check whether it encloses a solid."""
    import trimesh

    mesh = trimesh.load(str(path), force="mesh", process=True)
    watertight = bool(mesh.is_watertight)
    winding = bool(mesh.is_winding_consistent)
    volume = float(mesh.volume) if watertight and winding else float("nan")
    parts = mesh.split(only_watertight=False)
    # Density used by the paper: 1.055 g/cm^3. Convert mm^3 to cm^3.
    mass_g = volume * 1.055 / 1000 if np.isfinite(volume) and volume > 0 else float("nan")
    return {
        "topology_vertices": len(mesh.vertices),
        "topology_faces": len(mesh.faces),
        "watertight": watertight,
        "winding_consistent": winding,
        "connected_components": len(parts),
        "euler_number": int(mesh.euler_number),
        "mesh_volume_mm3": volume,
        "estimated_mass_g_density_1_055": mass_g,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=Path(__file__).resolve().parents[1] / "dataset")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "results" / "dataset_inventory.csv")
    parser.add_argument("--topology", action="store_true", help="weld each STL and check watertightness/components")
    parser.add_argument("--topology-only", action="store_true", help="add topology checks to an existing inventory CSV")
    args = parser.parse_args()

    paths = sorted(p for p in args.dataset.glob("*.stl") if not p.name.startswith("._"))
    if args.topology_only:
        with args.output.open(newline="") as f:
            rows = list(csv.DictReader(f))
    else:
        rows = []
        for i, path in enumerate(paths, 1):
            match = NAME_RE.match(path.stem)
            row = {
                "file": path.name,
                "patient_label": match.group("patient") if match else "",
                "scan_label": match.group("scan") if match else "",
                "phase_label": match.group("phase") if match else "",
                "size_bytes": path.stat().st_size,
                "status": "ok",
                "error": "",
            }
            try:
                row.update(inspect_binary_stl(path))
            except Exception as exc:  # Keep every file in the inventory even if malformed.
                try:
                    row.update(inspect_ascii_stl(path))
                except Exception as ascii_exc:
                    row["status"] = "needs_review"
                    row["error"] = f"binary: {exc}; ASCII: {ascii_exc}"
            rows.append(row)

    if args.topology or args.topology_only:
        path_by_name = {p.name: p for p in paths}
        for i, row in enumerate(rows, 1):
            path = path_by_name.get(row["file"])
            if path is None or row["status"] != "ok":
                continue
            try:
                row.update(inspect_topology(path))
            except Exception as exc:
                row["topology_error"] = str(exc)
            print(f"[topology {i}/{len(rows)}] {path.name}: {row.get('watertight', 'error')}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "file", "patient_label", "scan_label", "phase_label", "size_bytes",
        "format", "faces", "surface_area_mm2", "signed_volume_mm3_unverified",
        "bounds_min_x", "bounds_min_y", "bounds_min_z",
        "bounds_max_x", "bounds_max_y", "bounds_max_z", "status", "error",
        "topology_vertices", "topology_faces", "watertight", "winding_consistent",
        "connected_components", "euler_number", "mesh_volume_mm3",
        "estimated_mass_g_density_1_055", "topology_error",
    ]
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    groups = defaultdict(list)
    for row in rows:
        groups[(row["patient_label"], row["scan_label"])].append(row["phase_label"])
    print(f"\nInventoried {len(rows)} STL files; wrote {args.output}")
    print("Groups by filename (informational; no files are excluded):")
    for (patient, scan), phases in sorted(groups.items()):
        print(f"  {patient or '?'} scan {scan or '?'}: phases {', '.join(sorted(phases))}")


if __name__ == "__main__":
    main()
