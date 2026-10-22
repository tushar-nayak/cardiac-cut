"""CAMP Toolkit — run Convergent Adaptive Monge Patching on vascular-cut STL outputs.

Typical use (see run_camp.ipynb):

    from camp_toolkit import run_camp
    summary = run_camp("my_stls/", out_dir="results", target_faces=100_000, trim_frac=0.01)

The pipeline is:  STL -> weld -> uniform ACVD remesh -> CAMP -> outputs.
The remesh is the only preprocessing; everything else is CAMP with your parameters.
"""
from .pipeline import (
    run_camp,
    run_stl,
    weld_mesh,
    DEFAULT_PARAMS,
)
from .remesh import acvd_uniform_remesh
from .visualize import view, build_patch_figure, save_patch_html, save_patch_png
from . import campsrc

__all__ = [
    "run_camp",
    "run_stl",
    "weld_mesh",
    "acvd_uniform_remesh",
    "view",
    "build_patch_figure",
    "save_patch_html",
    "save_patch_png",
    "DEFAULT_PARAMS",
    "campsrc",
]
