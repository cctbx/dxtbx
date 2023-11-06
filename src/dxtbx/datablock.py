from __future__ import annotations

import warnings

from dxtbx.model.experiment_list import (
    BeamComparison,
    DetectorComparison,
    FormatChecker,
    GoniometerComparison,
)

__all__ = [
    "BeamComparison",
    "DetectorComparison",
    "FormatChecker",
    "GoniometerComparison",
]

# Remove this module after DIALS 3.17
warnings.warn(
    "dxtbx.datablock is deprecated; please import from dxtbx.experiment_list instead",
    UserWarning,
    stacklevel=2,
)
