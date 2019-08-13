from __future__ import absolute_import, division, print_function

from builtins import range
import copy
import os

from dxtbx.format.FormatPY import FormatPY
import six
import six.moves.cPickle as pickle
from past.builtins import basestring
from iotbx.detectors.cspad_detector_formats import reverse_timestamp
from iotbx.detectors.cspad_detector_formats import detector_format_version
from dxtbx import IncorrectFormatError
from spotfinder.applications.xfel import cxi_phil
from scitbx.array_family import flex
from iotbx.detectors.npy import NpyImage
import sys


class FormatPYunspecified(FormatPY):
    @staticmethod
    def understand(image_file):
        """Seems like the static method wastes a lot of effort here; it's not possible to
        just read the first few bytes; instead understand() reads the entire first data
        item in the file; an entire binary image.  This data is then read again in the
        _start() method and again in the detectorbase constructor."""
        try:
            with FormatPYunspecified.open_file(image_file, "rb") as fh:
                if six.PY3:
                    data = pickle.load(fh, encoding="bytes")
                    headers = {key.decode("ascii") for key in data}
                else:
                    data = pickle.load(fh)
                    headers = set(data)
        except IOError:
            return False

        wanted_header_items = {"SIZE1", "SIZE2", "TIMESTAMP"}

        return wanted_header_items.issubset(headers)

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatPY.__init__(self, image_file, **kwargs)

    def detectorbase_start(self):
        pass

    def _start(self):
        if isinstance(self._image_file, basestring) and os.path.isfile(
            self._image_file
        ):
            with FormatPYunspecified.open_file(self._image_file, "rb") as fh:
                if six.PY3:
                    data = pickle.load(fh, encoding="bytes")
                    data = {
                        key.decode("ascii"): value.decode("latin-1")
                        if isinstance(value, bytes)
                        else value
                        for key, value in data.items()
                    }
                else:
                    data = pickle.load(fh)
        else:
            data = self._image_file

        if "DETECTOR_ADDRESS" not in data:
            # legacy format; try to guess the address
            self.LCLS_detector_address = "CxiDs1-0|Cspad-0"
            if "DISTANCE" in data and data["DISTANCE"] > 1000:
                # downstream CS-PAD detector station of CXI instrument
                self.LCLS_detector_address = "CxiDsd-0|Cspad-0"
        else:
            self.LCLS_detector_address = data["DETECTOR_ADDRESS"]

        self._timesec = reverse_timestamp(data["TIMESTAMP"])[0]

        version_lookup = detector_format_version(
            self.LCLS_detector_address, self._timesec
        )
        self.start_helper(
            version_token="distl.detector_format_version=%s" % version_lookup
        )

    def start_helper(self, version_token):

        is_file = isinstance(self._image_file, basestring) and os.path.isfile(
            self._image_file
        )

        if is_file:
            file_name = self._image_file
        else:
            file_name = "inmem"

        args = [
            file_name,
            version_token,
            "viewer.powder_arcs.show=False",
            "viewer.powder_arcs.code=3n9c",
        ]

        params = cxi_phil.cxi_versioned_extract(args)
        horizons_phil = params.persist.commands

        if is_file:
            image = NpyImage(file_name)
        else:
            print(
                "This is not a file; assume the data are in the defined dictionary format"
            )
            image = NpyImage(file_name, source_data=self._image_file)
        image.readHeader(horizons_phil)
        image.translate_tiles(horizons_phil)
        # necessary to keep the phil parameters for subsequent calls to get_tile_manager()
        image.horizons_phil_cache = copy.deepcopy(horizons_phil)
        self.detectorbase = image

    def _goniometer(self):

        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector"""
        trusted_range = (
            self.detectorbase.parameters.get("MIN_TRUSTED_VALUE", 0),
            self.detectorbase.saturation,
        )

        return self._detector_factory.simple(
            sensor="PAD",
            distance=self.detectorbase.distance,
            beam_centre=(self.detectorbase.beamx, self.detectorbase.beamy),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(self.detectorbase.pixel_size, self.detectorbase.pixel_size),
            image_size=(self.detectorbase.size2, self.detectorbase.size1),
            trusted_range=trusted_range,
            mask=[],
        )

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(self.detectorbase.wavelength)

    def _scan(self):
        """Return the scan information for this image."""

        if (
            self.detectorbase.osc_start is not None
            and self.detectorbase.deltaphi is not None
            and self.detectorbase.deltaphi > 0
        ):
            format = ""

            exposure_time = self.detectorbase.parameters.get("TIME", 1)
            osc_start = self.detectorbase.osc_start
            osc_range = self.detectorbase.deltaphi
            timestamp = self._timesec

            return self._scan_factory.single(
                self._image_file, format, exposure_time, osc_start, osc_range, timestamp
            )

        else:
            return self._scan_factory.make_scan(
                image_range=(1, 1),
                # femtosecond X-ray pulse, set this to a dummy argument
                exposure_times=[1.0],
                oscillation=(0.0, 0.0),
                epochs={1: self._timesec},
            )

    def get_static_mask(self):
        """Creates a mask merging untrusted pixels with active areas."""
        detector_base = self.detectorbase
        # get effective active area coordinates
        tile_manager = detector_base.get_tile_manager(detector_base.horizons_phil_cache)
        tiling = tile_manager.effective_tiling_as_flex_int(
            reapply_peripheral_margin=True
        )
        if tiling is None or len(tiling) == 0:
            return None

        # get the raw data to get the size of the mask
        data = self.get_raw_data()

        # set the mask to the same dimensions as the data
        mask = flex.bool(flex.grid(data.focus()))

        # set active areas to True so they are not masked
        for i in range(len(tiling) // 4):
            x1, y1, x2, y2 = tiling[4 * i : (4 * i) + 4]
            sub_array = flex.bool(flex.grid(x2 - x1, y2 - y1), True)
            mask.matrix_paste_block_in_place(sub_array, x1, y1)
        return mask


class FormatPYunspecifiedInMemory(FormatPYunspecified):
    """Overrides the Format object's init method to accept an image dictionary
    instead of a file name. Used with XFELs when it is desirable to never write
    a file to disk, but to process it only in memory.
    """

    @staticmethod
    def understand(image_file):
        """ If it's an image dictionary, we understand this """
        wanted_header_items = ["SIZE1", "SIZE2", "TIMESTAMP"]
        try:
            if any(
                header_item not in image_file for header_item in wanted_header_items
            ):
                return False

            return True

        except AttributeError:
            return False

    def __init__(self, data, **kwargs):
        """ @param data In memory image dictionary, alredy initialized """
        FormatPYunspecified.__init__(self, data, **kwargs)

        self._image_file = copy.deepcopy(data)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPYunspecified.understand(arg))
