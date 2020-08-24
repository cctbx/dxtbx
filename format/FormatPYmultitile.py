from __future__ import absolute_import, division, print_function

import sys
from calendar import timegm
from time import strptime

import six
import six.moves.cPickle as pickle

from scitbx.matrix import col
from xfel.cftbx.detector.cspad_detector import CSPadDetector

from dxtbx.format.FormatPY import FormatPY
from dxtbx.model import Detector

from iotbx.detectors.npy import image_dict_to_unicode


class FormatPYmultitile(FormatPY):
    @staticmethod
    def understand(image_file):
        try:
            with FormatPYmultitile.open_file(image_file, "rb") as fh:
                if six.PY3:
                    data = pickle.load(fh, encoding="bytes")
                    data = image_dict_to_unicode(data)
                else:
                    data = pickle.load(fh)
        except IOError:
            return False

        wanted_header_items = ["TILES", "METROLOGY"]
        if any(header_item not in data for header_item in wanted_header_items):
            return False

        unwanted_header_items = ["SIZE1", "SIZE2"]
        if any(header_item in data for header_item in unwanted_header_items):
            return False

        return True

    def _start(self):
        # this Format class depends on stuff from detectorbase
        self.detectorbase_start()

    def detectorbase_start(self):
        self.detectorbase = CSPadDetector(self._image_file)
        self.detectorbase.readHeader()
        self._metrology_params = self.detectorbase._metrology_params
        self._tiles = self.detectorbase._tiles

    def _goniometer(self):
        return self._goniometer_factory.single_axis()

    def _detector(self):
        """The _detector() function returns a model for a CSPAD detector as
        used at LCLS's CXI and XPP endstations.  It converts the
        metrology information in the pure Python object extracted from
        the image pickle to DXTBX-style transformation vectors.  Only
        ASIC:s are considered, since DXTBX metrology is not concerned
        with hierarchies.

        Merged from xfel.cftbx.detector.cspad_detector.readHeader() and
        xfel.cftbx.detector.metrology.metrology_as_dxtbx_vectors().
        """

        # XXX Introduces dependency on cctbx.xfel!  Should probably be
        # merged into the code here!
        from xfel.cftbx.detector.metrology import _transform, get_projection_matrix

        # Apply the detector distance to the translation of the root
        # detector object.
        d = self._metrology_params.detector
        Tb_d = _transform(
            col(d.orientation).normalize(),
            col(d.translation) + col((0, 0, -self._metrology_params.distance * 1e-3)),
        )[1]

        self._raw_data = []
        detector = Detector()

        for p in d.panel:
            Tb_p = (
                Tb_d * _transform(col(p.orientation).normalize(), col(p.translation))[1]
            )

            for s in p.sensor:
                Tb_s = (
                    Tb_p
                    * _transform(col(s.orientation).normalize(), col(s.translation))[1]
                )

                for a in s.asic:
                    Tb_a = (
                        Tb_s
                        * _transform(
                            col(a.orientation).normalize(), col(a.translation)
                        )[1]
                    )

                    Pb = get_projection_matrix(a.pixel_size, a.dimension)[1]

                    # The DXTBX-style metrology description consists of three
                    # vectors for each ASIC.  The origin vector locates the
                    # (0, 0)-pixel in the laboratory frame in units of mm.
                    # The second and third vectors give the directions to the
                    # pixels immediately next to (0, 0) in the fast and slow
                    # directions, respectively, in arbitrary units.
                    origin = Tb_a * Pb * col((0, 0, 1))
                    fast = Tb_a * Pb * col((0, a.dimension[0], 1)) - origin
                    slow = Tb_a * Pb * col((a.dimension[1], 0, 1)) - origin

                    # Convert vector units from meter to millimeter.  The
                    # default, SimplePxMmStrategy applies here.  XXX Due to
                    # dark subtraction, a valid pixel intensity may be
                    # negative, and this is currently not reflected by
                    # trusted_range.
                    key = (d.serial, p.serial, s.serial, a.serial)

                    panel = detector.add_panel()
                    panel.set_type("PAD")
                    panel.set_name("%d:%d:%d:%d" % key)
                    panel.set_local_frame(
                        [t * 1e3 for t in fast.elems[0:3]],
                        [t * 1e3 for t in slow.elems[0:3]],
                        [t * 1e3 for t in origin.elems[0:3]],
                    )

                    panel.set_pixel_size([t * 1e3 for t in a.pixel_size])
                    panel.set_image_size(a.dimension)
                    panel.set_trusted_range((0, a.saturation))

                    self._raw_data.append(self._tiles[key])

        return detector

    def get_raw_data(self):
        """Return a tuple of flex arrays, one for each asic"""

        return tuple(self._raw_data)

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(self._metrology_params.wavelength)

    def _scan(self):
        """Return the scan information for this image."""

        # Convert textual ISO 8601 timestamp in UTC to
        # millisecond-precision Unix epoch.
        str_min = self._metrology_params.timestamp[:16] + "UTC"
        str_sec = self._metrology_params.timestamp[17:]
        epoch = timegm(strptime(str_min, "%Y-%m-%dT%H:%M%Z")) + float(str_sec)

        return self._scan_factory.make_scan(
            image_range=(1, 1),
            exposure_times=[1e-15 * self._metrology_params.pulse_length],
            oscillation=(0, 0),
            epochs={1: epoch},
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPYmultitile.understand(arg))
