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
        """Detector has always been on I19."""
        return self._at_I19

    def _base(self, mmcif=False):
        """Generates
        1. a CIF/mmCIF block that contains information that
        is always true about the detector.
        2. a lookup function for CIF/mmCIF strings."""
        # prepare string lookup table
        lookup = self._lookup(mmcif)

        b = iotbx.cif.model.block()
        b[lookup("df.detector")] = "Photon counting pixel array"
        b[lookup("df.detector.type")] = "Dectris PILATUS 300K"
        b[lookup("df.src")] = "Synchrotron"

        return b, lookup

    def _at_I19(self, mmcif=False):
        b, lookup = self._base(mmcif)

        b[lookup("df.m.dev")] = "4-circle \\k-geometry diffractometer"
        b[lookup("df.m.dev.type")] = "Newport IS4CCD"
        b[lookup("df.m.method")] = "shutterless scans"
        b[lookup("df.rad.mono")] = "Silicon 111"
        b[lookup("df.src.details")] = "Nowell et al. (2012)"
        b[lookup("df.src.type")] = "Diamond Light Source Beamline I19-2"
        b[
            lookup("references")
        ] = "Nowell, H. et al. (2012) J. Synchrotron Rad. 19, 435-441."
        b[lookup("sw.collection")] = "GDA - generic data acquisition software"

        return b
