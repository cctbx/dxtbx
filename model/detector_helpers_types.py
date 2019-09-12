"""
Helpers for the detector class... this time enumerating all of the common
detector types, hashed by the sensor type, image dimensions and pixel
dimensions.
"""

from __future__ import absolute_import, division, print_function

import io
import os
import sys
from builtins import object

import dxtbx
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.detector_helpers import detector_helper_sensors


class detector_helpers_types(object):
    """A singleton class to help with identifying specific detectors used for
    macromolecular crystallography."""

    def __init__(self):
        detector_lib = os.path.join(
            os.path.split(dxtbx.__file__)[0], "data", "detectors.lib"
        )

        if not os.path.exists(detector_lib):
            raise RuntimeError("detector library not found")

        self._detectors = {}

        with io.open(detector_lib, "r", encoding="ascii") as fh:
            for record in fh:
                if record.startswith(("Sensor", "-----")):
                    continue

                text = record.split("#")[0].strip()

                if not text:
                    continue

                tokens = text.split()

                assert len(tokens) == 6

                sensor = DetectorFactory.sensor(tokens[0])
                fast, slow, df, ds = map(int, tokens[1:5])

                self._detectors[(sensor, fast, slow, df, ds)] = tokens[5]

    def get(self, sensor, fast, slow, df, ds):
        """Look up a name for a detector with this sensor type (listed in
        detector_helpers) these image dimensions in the fast and slow
        directions, these corresponding pixel sizes, in microns (integers).
        If the sensor is unknown, all sensor types will be tested - be warned
        if there are duplicates."""

        sensor = DetectorFactory.sensor(sensor)

        if sensor == detector_helper_sensors.SENSOR_UNKNOWN:
            for s in detector_helper_sensors.all():
                try:
                    return self.get(s, fast, slow, df, ds)
                except ValueError:
                    pass

            raise ValueError(
                "detector %s %d %d %d %d unknown" % (sensor, fast, slow, df, ds)
            )

        if (sensor, fast, slow, df, ds) in self._detectors:
            return self._detectors[(sensor, fast, slow, df, ds)]

        # OK allow for small variations in the recorded pixel dimensions

        for ddf in -2, -1, 1, 2:
            for dds in -2, -1, 1, 2:
                if (sensor, fast, slow, df + ddf, ds + dds) in self._detectors:
                    return self._detectors[(sensor, fast, slow, df + ddf, ds + dds)]

        raise ValueError(
            "detector %s %d %d %d %d unknown" % (sensor, fast, slow, df, ds)
        )


detector_helpers_types = detector_helpers_types()

if __name__ == "__main__":
    sensor = sys.argv[1]
    fast, slow, df, ds = map(int, sys.argv[2:6])

    print(detector_helpers_types.get(sensor, fast, slow, df, ds))
