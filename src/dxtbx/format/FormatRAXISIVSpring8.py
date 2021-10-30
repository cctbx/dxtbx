import calendar
import datetime
import struct

from dxtbx.format.Format import Format
from dxtbx.format.FormatRAXIS import RAXISHelper


class FormatRAXISIVSPring8(RAXISHelper, Format):
    """Format class for R-AXIS4 images. Currently the only example we have is
    from SPring-8, which requires a reverse axis goniometer. It is not clear how
    to distinguish this detector from others that produce 'R-AXIS4' images that we
    may come across in future. This will be a problem if other detectors do not
    also require reverse axis.
    """

    @staticmethod
    def understand(image_file):
        try:
            with Format.open_file(image_file, "rb") as fh:
                header = fh.read(1024)
        except OSError:
            return False

        # A few items expected to be the same from image to image that we can use
        # as a fingerprint for this instrument
        if header[0:7] != b"R-AXIS4":
            return False

        # We expect an invalid goniometer section, indicated by wrong magic number
        if struct.unpack(">i", header[852:856]) == 1:
            return False

        if header[812:822].strip() != b"IRIS":
            return False

        return True

    def _goniometer(self):
        return self._goniometer_factory.single_axis_reverse()

    def _detector(self):
        """Return a model for the detector as defined in the image header,
        with the additional knowledge about how things are arranged i.e. that
        the principle rotation axis vector points from the sample downwards."""

        values = self._detector_helper()

        return self._detector_factory.simple(
            "IMAGE_PLATE",
            values["distance"],
            values["beam"],
            "+x",
            "+y",
            (values["dx"], values["dy"]),
            (values["nx"], values["ny"]),
            (0, 1000000),
            [],
        )

    def _beam(self):
        """Return a simple model for the beam."""

        wavelength = struct.unpack(self._f, self._header_bytes[292:296])[0]

        return self._beam_factory.simple(wavelength)

    def _scan(self):
        """Return the scan information for this image."""
        f = self._f
        header = self._header_bytes

        exposure_time = struct.unpack(f, header[536:540])[0]

        y, m, d = map(int, header[256:268].strip().split(b"-"))

        epoch = calendar.timegm(datetime.datetime(y, m, d, 0, 0, 0).timetuple())

        # For this instrument, the header goniometer section is invalid. 3 floats
        # starting at byte 520 specify phi0, phis and phie.
        # (http://www.rigaku.com/downloads/software/readimage.html)

        osc_dat = struct.unpack(f, header[520:524])[0]
        osc_start = osc_dat + struct.unpack(f, header[524:528])[0]
        osc_end = osc_dat + struct.unpack(f, header[528:532])[0]

        osc_range = osc_end - osc_start

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )
