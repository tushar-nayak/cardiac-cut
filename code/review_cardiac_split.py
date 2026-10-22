#!/usr/bin/env python3
"""Interactive review and manual basal trim of automatically split cardiac meshes.

Run the automatic splitter first with --keep-outer-lip so this viewer starts from
the untrimmed outer wall. In the viewer, move/rotate the plane on the outer wall,
inspect the kept region, then save the plane and surfaces as a JSON/STL bundle.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def _qt_app():
    from PySide6 import QtWidgets
    from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
    from vtkmodules.vtkCommonColor import vtkNamedColors
    from vtkmodules.vtkCommonCore import vtkVersion
    from vtkmodules.vtkCommonDataModel import vtkPlane
    from vtkmodules.vtkFiltersCore import vtkClipPolyData
    from vtkmodules.vtkIOGeometry import vtkSTLReader, vtkSTLWriter
    from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer
    from vtkmodules.vtkInteractionWidgets import (
        vtkImplicitPlaneRepresentation, vtkImplicitPlaneWidget2,
    )

    class SplitReview(QtWidgets.QMainWindow):
        def __init__(self, split_dir: Path, output_dir: Path):
            super().__init__()
            self.split_dir = split_dir
            self.output_dir = output_dir
            self.setWindowTitle("Cardiac surface split review")
            self.resize(1500, 850)
            self.cases = sorted(p.name.removesuffix("_inner.stl")
                                for p in split_dir.glob("*_inner.stl"))
            if not self.cases:
                raise RuntimeError(f"No *_inner.stl files found in {split_dir}")
            self.readers = {}
            self.clipper = vtkClipPolyData()
            self.plane = vtkPlane()
            self.clipper.SetClipFunction(self.plane)
            self.clipper.InsideOutOff()
            self.clipper.GenerateClipScalarsOff()
            self.clipper.Update()
            self._build_ui()
            self.case_box.setCurrentIndex(0)
            self.load_case(self.cases[0])

        def _reader(self, path: Path):
            key = str(path)
            if key not in self.readers:
                reader = vtkSTLReader()
                reader.SetFileName(key)
                reader.Update()
                if reader.GetOutput().GetNumberOfPoints() == 0:
                    raise RuntimeError(f"Could not read mesh: {path}")
                self.readers[key] = reader
            return self.readers[key]

        def _paths(self, case):
            inner = self.split_dir / f"{case}_inner.stl"
            raw_outer = self.split_dir / f"{case}_outer_untrimmed.stl"
            outer = raw_outer if raw_outer.exists() else self.split_dir / f"{case}_outer.stl"
            return inner, outer, raw_outer.exists()

        def _outer_reader(self):
            if not hasattr(self, "case_box"):
                return self._reader(self.split_dir / f"{self.cases[0]}_outer_untrimmed.stl")
            _, path, _ = self._paths(self.case_box.currentText() or self.cases[0])
            return self._reader(path)

        def _build_ui(self):
            root = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(root)
            controls = QtWidgets.QHBoxLayout()
            controls.addWidget(QtWidgets.QLabel("Mesh:"))
            self.case_box = QtWidgets.QComboBox()
            self.case_box.addItems(self.cases)
            controls.addWidget(self.case_box, 1)
            self.side_box = QtWidgets.QComboBox()
            self.side_box.addItems(["Keep positive side", "Keep negative side"])
            # The initial plane is near the PCA-axis maximum, so default to
            # retaining the larger, negative side of that cut.
            self.side_box.setCurrentIndex(1)
            controls.addWidget(self.side_box)
            self.save_button = QtWidgets.QPushButton("Save reviewed surfaces")
            controls.addWidget(self.save_button)
            layout.addLayout(controls)
            views = QtWidgets.QHBoxLayout()
            self.inner_view = QVTKRenderWindowInteractor(root)
            self.outer_view = QVTKRenderWindowInteractor(root)
            views.addWidget(self.inner_view, 1)
            views.addWidget(self.outer_view, 1)
            layout.addLayout(views, 1)
            self.status = QtWidgets.QLabel(
                "Left: inner surface. Right: outer surface; drag the plane handles to trim the basal rim."
            )
            layout.addWidget(self.status)
            self.setCentralWidget(root)

            self.inner_renderer = vtkRenderer()
            self.outer_renderer = vtkRenderer()
            colors = vtkNamedColors()
            for renderer, color in ((self.inner_renderer, "SlateGray"),
                                    (self.outer_renderer, "LightSteelBlue")):
                renderer.SetBackground(colors.GetColor3d("WhiteSmoke"))
            self.inner_view.GetRenderWindow().AddRenderer(self.inner_renderer)
            self.outer_view.GetRenderWindow().AddRenderer(self.outer_renderer)
            self.inner_view.Initialize()
            self.outer_view.Initialize()

            self.case_box.currentTextChanged.connect(self.load_case)
            self.side_box.currentIndexChanged.connect(self._side_changed)
            self.save_button.clicked.connect(self.save_case)
            self.plane_rep = vtkImplicitPlaneRepresentation()
            self.plane_widget = vtkImplicitPlaneWidget2()
            self.plane_widget.SetInteractor(self.outer_view.GetRenderWindow().GetInteractor())
            self.plane_widget.SetRepresentation(self.plane_rep)
            self.plane_widget.AddObserver("InteractionEvent", self._plane_changed)

        @staticmethod
        def _actor(mapper, color):
            actor = vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(color)
            actor.GetProperty().SetInterpolationToPhong()
            actor.GetProperty().SetSpecular(0.18)
            actor.GetProperty().SetSpecularPower(18)
            return actor

        def load_case(self, case):
            if not case:
                return
            inner_path, outer_path, has_untrimmed = self._paths(case)
            inner_mapper = vtkPolyDataMapper()
            inner_mapper.SetInputConnection(self._reader(inner_path).GetOutputPort())
            self.inner_renderer.RemoveAllViewProps()
            self.inner_renderer.AddActor(self._actor(inner_mapper, (0.72, 0.35, 0.34)))
            outer_reader = self._reader(outer_path)
            self.clipper.SetInputConnection(outer_reader.GetOutputPort())
            outer_mapper = vtkPolyDataMapper()
            outer_mapper.SetInputConnection(self.clipper.GetOutputPort())
            self.outer_renderer.RemoveAllViewProps()
            self.outer_renderer.AddActor(self._actor(outer_mapper, (0.30, 0.55, 0.78)))
            self.outer_reader = outer_reader
            self._set_initial_plane()
            self.inner_renderer.ResetCamera()
            self.outer_renderer.ResetCamera()
            self.inner_view.GetRenderWindow().Render()
            self.outer_view.GetRenderWindow().Render()
            self.plane_widget.SetEnabled(1)
            self.status.setText(
                f"{case} — left: inner. Right: kept outer preview. "
                + ("Using untrimmed outer surface." if has_untrimmed else
                   "WARNING: untrimmed outer is missing; regenerate with --keep-outer-lip.")
            )

        def _set_initial_plane(self):
            poly = self.outer_reader.GetOutput()
            bounds = poly.GetBounds()
            # Initially orient perpendicular to the longest PCA axis. The clinician
            # can rotate and translate it with the VTK plane handles.
            pts = np.array([poly.GetPoint(i) for i in range(poly.GetNumberOfPoints())])
            mean = pts.mean(axis=0)
            _, vecs = np.linalg.eigh(np.cov((pts - mean).T))
            normal = vecs[:, -1]
            dominant_axis = int(np.argmax(np.abs(normal)))
            if normal[dominant_axis] < 0:
                normal *= -1
            projection = (pts - mean) @ normal
            center = mean + normal * (projection.max() - 0.12 * (projection.max() - projection.min()))
            self.plane.SetOrigin(center)
            self.plane.SetNormal(normal)
            self.plane_rep.SetPlaceFactor(1.15)
            self.plane_rep.PlaceWidget(bounds)
            self.plane_rep.SetOrigin(center)
            self.plane_rep.SetNormal(normal)
            self._sync_plane()

        def _sync_plane(self):
            self.plane_rep.GetPlane(self.plane)
            self.clipper.SetValue(0.0)
            self.clipper.SetInsideOut(self.side_box.currentIndex() == 1)
            self.clipper.Update()
            self.outer_view.GetRenderWindow().Render()

        def _plane_changed(self, *_):
            self._sync_plane()

        def _side_changed(self, _):
            self._sync_plane()

        def save_case(self):
            case = self.case_box.currentText()
            self._sync_plane()
            self.clipper.Update()
            output = self.output_dir / case
            output.mkdir(parents=True, exist_ok=True)
            inner_in, outer_in, has_untrimmed = self._paths(case)
            if not has_untrimmed:
                QtWidgets.QMessageBox.warning(
                    self, "Untrimmed outer surface missing",
                    "Regenerate the split files with --keep-outer-lip before reviewing."
                )
                return
            # Copy the inner surface byte-for-byte; write the manually clipped outer.
            inner_out = output / f"{case}_inner.stl"
            inner_out.write_bytes(inner_in.read_bytes())
            outer_out = output / f"{case}_outer.stl"
            writer = vtkSTLWriter()
            writer.SetFileName(str(outer_out))
            writer.SetInputConnection(self.clipper.GetOutputPort())
            if writer.Write() != 1:
                QtWidgets.QMessageBox.critical(self, "Save failed", "Could not write clipped STL.")
                return
            origin = self.plane.GetOrigin()
            normal = self.plane.GetNormal()
            digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
            metadata = {
                "source_inner_stl": str(inner_in.resolve()),
                "source_outer_stl": str(outer_in.resolve()),
                "source_inner_sha256": digest(inner_in),
                "source_outer_sha256": digest(outer_in),
                "inner_stl": inner_out.name,
                "outer_stl": outer_out.name,
                "plane_origin": origin,
                "plane_normal": normal,
                "kept_side": "negative" if self.side_box.currentIndex() == 1 else "positive",
                "saved_at_utc": datetime.now(timezone.utc).isoformat(),
                "split_method": "automatic geodesic split; clinician-adjusted planar outer trim",
                "camp_source": "camp/camp_updated_2026-08-26/campsrc.py",
                "vtk_version": vtkVersion.GetVTKVersion(),
            }
            (output / "review.json").write_text(json.dumps(metadata, indent=2) + "\n")
            self.status.setText(f"Saved reviewed inner/outer surfaces and plane: {output}")

    return QtWidgets, SplitReview


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--split-dir", type=Path, default=root / "results" / "manual_split_inputs")
    parser.add_argument("--output", type=Path, default=root / "results" / "manual_surface_split")
    args = parser.parse_args()
    QtWidgets, Window = _qt_app()
    app = QtWidgets.QApplication(sys.argv)
    window = Window(args.split_dir, args.output)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
