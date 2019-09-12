from __future__ import absolute_import, division, print_function

import iotbx.cif.model

import dxtbx.data.beamline_defs


class get_definition(dxtbx.data.beamline_defs.template):
    def __init__(self, timestamp=None, **kwargs):
        self._timestamp = timestamp

    def CIF_block(self):
        """Interface function to generate a CIF block for this detector."""
        return self._identify_time()(mmcif=False)

    def mmCIF_block(self):
        """Interface function to generate an mmCIF block for this detector."""
        return self._identify_time()(mmcif=True)

    def _identify_time(self):
        """Determine detector environment based on timestamp."""
        if not self._timestamp:  # default to I19
            return self._at_I19

        if self._timestamp >= self._date_to_epoch(2015, 11, 1):
            return self._at_I19  # moved to I19 on 01.11.2015

        return self._at_I04_1  # before was on I04-1

    def _base(self, mmcif=False):
        """Generates
        1. a CIF/mmCIF block that contains information that
        is always true about the detector.
        2. a lookup function for CIF/mmCIF strings."""
        # prepare string lookup table
        lookup = self._lookup(mmcif)

        b = iotbx.cif.model.block()
        b[lookup("df.detector")] = "Photon counting pixel array"
        b[lookup("df.rad.type")] = "Synchrotron"

        return b, lookup

    def _at_I19(self, mmcif=False):
        b, lookup = self._base(mmcif)

        b[lookup("df.m.dev")] = "Fixed \\c 3-circle diffractometer"
        b[lookup("df.m.dev_type")] = "Fluid Film Devices"
        b[lookup("df.m.method")] = "shutterless scans"
        b[lookup("df.m.spec_supp")] = "MiTeGen MicroMount"
        b[lookup("df.rad.source")] = "Diamond Light Source Beamline I19-1"
        b[lookup("df.rad.mono")] = "Silicon 111"

        return b

    def _at_I04_1(self, mmcif=False):
        b, lookup = self._base(mmcif)

        b[lookup("df.rad.source")] = "Diamond Light Source Beamline I04-1"

        return b
