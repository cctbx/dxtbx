#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
"""
An implementation of the SMV image reader for ADSC images. Inherits from
FormatSMVADSC, customised for example on Australian Synchrotron SN 457
which has reversed phi.
"""

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSNSN457(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 457."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 457."""

        # check this is detector serial number 457

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        return int(header["DETECTOR_SN"]) == 457

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header."""

        return self._goniometer_factory.single_axis_reverse()
