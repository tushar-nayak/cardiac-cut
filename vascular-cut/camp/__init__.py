"""CAMP, the lab's convergent adaptive patch partition, vendored.

Import from here rather than reaching for a file path:

    from camp import camp, trim_mask_integrated

`campsrc.py` beside this file is a verbatim copy of the lab's source; see
README.md for where it came from and how to refresh it. Nothing in this package
is modified, so a diff against upstream stays meaningful.
"""

from .campsrc import (  # noqa: F401
    camp,
    load_surface_from_mat,
    make_plotly_patch_figure,
    _trim_mask_integrated,
)

# The upstream name is private by convention but is part of how CAMP is actually
# used -- the partition does not depend on the trim, so callers run camp() once at
# trim_frac=0 and apply the mask themselves. Give it a public alias rather than
# have every caller reach for an underscore.
trim_mask_integrated = _trim_mask_integrated

__all__ = [
    "camp",
    "trim_mask_integrated",
    "load_surface_from_mat",
    "make_plotly_patch_figure",
]
