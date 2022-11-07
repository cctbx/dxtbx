from __future__ import annotations

import json
import logging

import libtbx

from dxtbx.format.Registry import get_format_class_for_file

logger = logging.getLogger(__name__)


class FormatChecker:
    """A helper class to speed up identifying the correct image format by first
    trying the last format that was used."""

    def __init__(self):
        """Set the format class to none."""
        self._format_class = None

    def find_format(self, filename):
        """Search the registry for the image format class.
        Where possible use the last seen format class as a prioritisation hint.
        """
        if self._format_class:
            self._format_class = get_format_class_for_file(
                filename, format_hint=self._format_class.__name__
            )
        else:
            self._format_class = get_format_class_for_file(filename)
        if self._format_class:
            logger.debug("Using %s for %s", self._format_class.__name__, filename)
        else:
            logger.debug("No format class found for %s", filename)
        return self._format_class

    def iter_groups(self, filenames):
        group_format = None
        group_fnames = []
        for filename in filenames:
            fmt = self.find_format(filename)
            if fmt == group_format:
                group_fnames.append(filename)
            else:
                if group_fnames:
                    yield group_format, group_fnames
                group_fnames = [filename]
                group_format = fmt
            if fmt is not None:
                logger.debug("Using %s for %s", fmt.__name__, filename)
        if group_fnames:
            yield group_format, group_fnames


class AutoEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, libtbx.AutoType):
            return "Auto"
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


class BeamComparison:
    """A class to provide simple beam comparison"""

    def __init__(
        self,
        wavelength_tolerance=1e-6,
        direction_tolerance=1e-6,
        polarization_normal_tolerance=1e-6,
        polarization_fraction_tolerance=1e-6,
    ):
        self.wavelength_tolerance = wavelength_tolerance
        self.direction_tolerance = direction_tolerance
        self.polarization_normal_tolerance = polarization_normal_tolerance
        self.polarization_fraction_tolerance = polarization_fraction_tolerance

    def __call__(self, a, b):
        if a is None and b is None:
            return True
        return a.is_similar_to(
            b,
            wavelength_tolerance=self.wavelength_tolerance,
            direction_tolerance=self.direction_tolerance,
            polarization_normal_tolerance=self.polarization_normal_tolerance,
            polarization_fraction_tolerance=self.polarization_fraction_tolerance,
        )


class DetectorComparison:
    """A class to provide simple detector comparison"""

    def __init__(
        self, fast_axis_tolerance=1e-6, slow_axis_tolerance=1e-6, origin_tolerance=1e-6
    ):
        self.fast_axis_tolerance = fast_axis_tolerance
        self.slow_axis_tolerance = slow_axis_tolerance
        self.origin_tolerance = origin_tolerance

    def __call__(self, a, b):
        if a is None and b is None:
            return True
        return a.is_similar_to(
            b,
            fast_axis_tolerance=self.fast_axis_tolerance,
            slow_axis_tolerance=self.slow_axis_tolerance,
            origin_tolerance=self.origin_tolerance,
        )


class GoniometerComparison:
    """A class to provide simple goniometer comparison"""

    def __init__(
        self,
        rotation_axis_tolerance=1e-6,
        fixed_rotation_tolerance=1e-6,
        setting_rotation_tolerance=1e-6,
    ):
        self.rotation_axis_tolerance = rotation_axis_tolerance
        self.fixed_rotation_tolerance = fixed_rotation_tolerance
        self.setting_rotation_tolerance = setting_rotation_tolerance

    def __call__(self, a, b):
        if a is None and b is None:
            return True
        elif a is None or b is None:
            return False
        return a.is_similar_to(
            b,
            rotation_axis_tolerance=self.rotation_axis_tolerance,
            fixed_rotation_tolerance=self.fixed_rotation_tolerance,
            setting_rotation_tolerance=self.setting_rotation_tolerance,
        )
