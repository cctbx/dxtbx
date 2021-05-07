"""An implementation of the CBF image reader for Eiger images"""


import os
import sys

from iotbx.detectors.pilatus_minicbf import PilatusImage

from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.format.FormatPilatusHelpers import get_vendortype_eiger as gv


class FormatCBFMiniEigerPhotonFactory(FormatCBFMini):
    """A class for reading mini CBF format Eiger images, and correctly
    constructing a model for the experiment from this.

    Specific for 2-Eiger setup at Photon Factory based on images sent by
    Yusuke Yamada. Geometry currently hard-coded as no geometry information in
    image headers. Assumes image filenames contain the string '_upper_' or
    '_lower_' to distinguish the two detectors.

    Only works if environment variable ENABLE_PHOTON_FACTORY_TWO_EIGER is set.
    """

    @staticmethod
    def understand(image_file):
        # no way to uniquely identify given example images

        if "ENABLE_PHOTON_FACTORY_TWO_EIGER" in os.environ:
            return True
        return False

    def _goniometer(self):
        return self._goniometer_factory.make_goniometer(
            (1, 0, 0), (1, 0, 0, 0, 1, 0, 0, 0, 1)
        )

    def _detector(self):

        if "_lower_" in self._image_file:
            # Detector:
            # Panel:
            # pixel_size:{0.075,0.075}
            # image_size: {2070,2167}
            # trusted_range: {0,1e+07}
            # thickness: 0
            # material:
            # mu: 0
            # fast_axis: {0.999997,-0.00228845,0.00127629}
            # slow_axis: {-0.0026197,-0.862841,0.505469}
            # origin: {-79.3767,2.47119,-169.509}

            detector = self._detector_factory.complex(
                "PAD",
                origin=(-79.3767, 2.47119, -169.509),
                fast=(0.999997, -0.00228845, 0.00127629),
                slow=(-0.0026197, -0.862841, 0.505469),
                pixel=(0.075, 0.075),
                size=(2070, 2167),
                trusted_range=(-1, 1e7),  # XXX FIXME
            )

        elif "_upper_" in self._image_file:
            # Detector:
            # Panel:
            # pixel_size:{0.075,0.075}
            # image_size: {2070,2167}
            # trusted_range: {0,1e+07}
            # thickness: 0
            # material:
            # mu: 0
            # fast_axis: {-0.999993,0.00338011,-0.00156153}
            # slow_axis: {0.00213163,0.863573,0.504219}
            # origin: {75.1912,14.1117,-124.34}

            detector = self._detector_factory.complex(
                "PAD",
                origin=(75.1912, 14.1117, -124.34),
                fast=(-0.999993, 0.00338011, -0.00156153),
                slow=(0.00213163, 0.863573, 0.504219),
                pixel=(0.075, 0.075),
                size=(2070, 2167),
                trusted_range=(-1, 1e7),  # XXX FIXME
            )

        else:
            raise RuntimeError("Don't understand image")

        return detector

    def _beam(self):
        # Beam:
        # wavelength: 1.1
        # sample to source direction : {-0.00573979,-0,0.999984}
        # divergence: 0
        # sigma divergence: 0
        # polarization normal: {0,1,0}
        # polarization fraction: 0.999

        return self._beam_factory.complex(
            sample_to_source=(-0.00573979, -0, 0.999984),
            polarization_fraction=0.999,
            polarization_plane_normal=(0, 1, 0),
            wavelength=1.1,
        )

    def _scan(self):
        exposure_time = 1  # XXX
        osc_start = float(self._cif_header_dictionary["Start_angle"].split()[0])
        osc_range = 0.1
        timestamp = 1  # XXX

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, timestamp
        )

    def detectorbase_start(self):

        self.detectorbase = PilatusImage(self._image_file)
        self.detectorbase.readHeader()

        self.detectorbase.parameters["SIZE1"] = 2167
        self.detectorbase.parameters["SIZE2"] = 2070

    def get_vendortype(self):
        return gv(self.get_detector())


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniEigerPhotonFactory.understand(arg))
