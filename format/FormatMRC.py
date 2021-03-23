"""Classes for Formats that recognise data stored in the MRC format, an open
standard used in electron microscopy
(http://www.ccpem.ac.uk/mrc_format/mrc2014.php)"""

from __future__ import absolute_import, division, print_function

import logging
import os
import re
from math import sqrt

import mrcfile
from mrcfile.constants import MAP_ID, MAP_ID_OFFSET_BYTES

from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.model import ScanFactory

logger = logging.getLogger("dials")


class FormatMRC(Format):
    @staticmethod
    def understand(image_file):
        with FormatMRC.open_file(image_file, "rb") as fh:
            start = fh.read(MAP_ID_OFFSET_BYTES + len(MAP_ID))
            if start[-len(MAP_ID) :] == MAP_ID:
                return True
        return False

    def _start(self):
        """Open the MRC file, read the metadata into an internal dictionary
        self._header_dictionary, add FEI extended metadata if available"""

        with mrcfile.open(self._image_file, header_only=True) as mrc:
            h = mrc.header
            xh = mrc.extended_header

        self._header_dictionary = self._unpack_header(h)
        if xh:
            self._extend_header(xh)

    @staticmethod
    def _unpack_header(header):
        hd = {}
        # What do we need from the header?
        fields = ("nx", "ny", "nz", "mx", "my", "mz")
        for key in fields:
            hd[key] = int(header[key])
        hd["exttyp"] = header["exttyp"].item()

        # For image stacks, NX==MX etc. should always be true. Assert this
        # to ensure we fail on an MRC file of the wrong type.
        assert hd["nx"] == hd["mx"]
        assert hd["ny"] == hd["my"]
        assert hd["nz"] == hd["mz"]

        return hd

    def _extend_header(self, xh):

        # Extensions from FEI headers only
        if not self._header_dictionary["exttyp"].startswith(b"FEI"):
            return

        self._header_dictionary["alphaTilt"] = xh["Alpha tilt"][0]
        self._header_dictionary["integrationTime"] = xh["Integration time"][0]
        self._header_dictionary["tilt_axis"] = xh["Tilt axis angle"][0]
        self._header_dictionary["pixelSpacing"] = xh["Pixel size X"][0]
        assert self._header_dictionary["pixelSpacing"] == xh["Pixel size Y"][0]
        self._header_dictionary["acceleratingVoltage"] = xh["HT"][0]
        self._header_dictionary["camera"] = xh["Camera name"][0]
        self._header_dictionary["binning"] = xh["Binning Width"][0]
        self._header_dictionary["noiseReduction"] = xh["Ceta noise reduction"][0]
        if b"Ceta" in self._header_dictionary["camera"]:
            # Does this ever differ from the Binning Width from the header?
            assert (
                self._header_dictionary["binning"]
                == 4096 / self._header_dictionary["nx"]
            )

        # Make an assumption for Thermo Fisher (Ceta-D and Falcon) detectors.
        # Is this stored anywhere else?
        self._header_dictionary["physicalPixel"] = 14e-6

        def cal_wavelength(V0):
            h = 6.626e-34  # Js, Planck's constant
            m = 9.109e-31  # kg, electron mass
            e = 1.6021766208e-19  # C, electron charge
            c = 3e8  # m/s^2, speed

            # Default to e- wavelength at 200 keV if voltage set to zero
            if V0 == 0:
                V0 = 200000
            return (
                h / sqrt(2 * m * e * V0 * (1 + e * V0 / (2 * m * c * c))) * 1e10
            )  # return wavelength in Angstrom

        self._header_dictionary["wavelength"] = cal_wavelength(
            self._header_dictionary["acceleratingVoltage"]
        )
        self._header_dictionary["cameraLength"] = (
            self._header_dictionary["physicalPixel"]
            * self._header_dictionary["binning"]
        ) / (
            self._header_dictionary["pixelSpacing"]
            * self._header_dictionary["wavelength"]
            * 1e-10
        )

        # Include FEI2 items. This should be detectable by xh["Metadata version"] >= 2
        # but in practice this is not used in all files I have seen.
        if self._header_dictionary["exttyp"] != b"FEI2":
            return
        self._header_dictionary["scanRotation"] = xh["Scan rotation"]
        self._header_dictionary["diffractionPatternRotation"] = xh[
            "Diffraction pattern rotation"
        ]
        self._header_dictionary["imageRotation"] = xh["Image rotation"]
        self._header_dictionary["scanModeEnum"] = xh["Scan mode enumeration"]
        self._header_dictionary["acquisitionTimeStamp"] = xh["Acquisition time stamp"]
        self._header_dictionary["detectorCommercialName"] = xh[
            "Detector commercial name"
        ]
        self._header_dictionary["startTiltAngle"] = xh["Start tilt angle"]
        self._header_dictionary["endTiltAngle"] = xh["End tilt angle"]
        self._header_dictionary["tiltPerImage"] = xh["Tilt per image"]
        self._header_dictionary["tiltSpeed"] = (xh["Tilt speed"],)
        self._header_dictionary["beamCentreXpx"] = xh["Beam center X pixel"]
        self._header_dictionary["beamCentreYpx"] = xh["Beam center Y pixel"]
        self._header_dictionary["cfegFlashTimestamp"] = xh["CFEG flash timestamp"]
        self._header_dictionary["phasePlatePositionIndex"] = xh[
            "Phase plate position index"
        ]
        self._header_dictionary["objectiveApertureName"] = xh["Objective aperture name"]

        return

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer."""

        direction = matrix.col((0.0, -1.0, 0.0))
        rot_by_deg = self._header_dictionary.get("tilt_axis", 0.0)
        direction.rotate(axis=matrix.col((0.0, 0.0, 1.0)), angle=rot_by_deg, deg=True)

        return self._goniometer_factory.known_axis(direction)

    def _detector(self):
        """Dummy detector"""

        image_size = (self._header_dictionary["nx"], self._header_dictionary["ny"])

        # Get pixel size, defaulting to 14 um for the Ceta if unknown
        physical_pixel = self._header_dictionary.get("physicalPixel", 1.4e-5)
        binning = self._header_dictionary.get("binning")
        if binning is None:
            if image_size == (2048, 2048):
                binning = 2.0
            else:
                binning = 1.0
        pixel_size = physical_pixel * 1000.0 * binning

        # The best distance measure is calculated from the calibrated pixel
        # size. If this is not available then default to the nominal camera
        # length, or finally to 2.0m
        distance = self._header_dictionary.get("cameraLength", 2.0) * 1000
        try:
            calibrated_pixel_size = self._header_dictionary["pixelSpacing"]  # 1/m
            wavelength = self._header_dictionary["wavelength"] * 1e-10  # m
            distance = pixel_size / (calibrated_pixel_size * wavelength)  # mm
        except KeyError:
            pass

        # Get detector-specific details for TF detectors as discussed with
        # Lingbo Yu. Ceta has gain of > 26 and Ceta and Falcon III both saturate
        # at about 8000.0 for binning=1
        camera = self._header_dictionary.get("camera", b"").lower()
        if b"ceta" in camera:
            gain = 26.0
            saturation = 8000 * binning ** 2
        elif b"falcon" in camera:
            gain = 1.0
            saturation = 8000 * binning ** 2
        else:
            gain = 1.0
            saturation = 1e6
        trusted_range = (-1000, saturation)

        # Beam centre not in the header - set to the image centre
        beam_centre = [(pixel_size * i) / 2 for i in image_size]
        detector = self._detector_factory.simple(
            "PAD",
            distance,
            beam_centre,
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
            trusted_range,
        )

        for panel in detector:
            panel.set_gain(gain)
        return detector

    def _beam(self):
        """Unpolarized beam model"""

        # Default to 200 keV
        wavelength = self._header_dictionary.get("wavelength", 0.02508)
        return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=wavelength,
            polarization=(0, 1, 0),
            polarization_fraction=0.5,
        )


class FormatMRCimages(FormatMRC):
    @staticmethod
    def understand(image_file):
        with mrcfile.open(image_file, header_only=True) as mrc:
            return int(mrc.header["nz"]) == 1

    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        Format.__init__(self, image_file, **kwargs)

    def _scan(self):
        """Scan model for this image, filling out any unavailable items with
        dummy values"""

        alpha = self._header_dictionary.get("alphaTilt", 0.0)
        dalpha = self._header_dictionary.get("tiltPerImage", 1.0)
        exposure = self._header_dictionary.get("integrationTime", 0.0)
        oscillation = (alpha, dalpha)
        fname = os.path.split(self._image_file)[-1]
        # assume that the final number before the extension is the image number
        s = fname.split("_")[-1].split(".")[0]
        try:
            index = int(re.match(".*?([0-9]+)$", s).group(1))
        except AttributeError:
            index = 1
        return ScanFactory.make_scan((index, index), exposure, oscillation, {index: 0})

    def get_raw_data(self):

        with mrcfile.open(self._image_file) as mrc:
            image = flex.double(mrc.data.astype("double"))

        return image


class FormatMRCstack(FormatMultiImage, FormatMRC):
    @staticmethod
    def understand(image_file):
        with mrcfile.open(image_file, header_only=True) as mrc:
            return int(mrc.header["nz"]) > 1

    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatMultiImage.__init__(self, **kwargs)
        Format.__init__(self, image_file, **kwargs)

    def get_num_images(self):
        return self._header_dictionary["nz"]

    def get_goniometer(self, index=None):
        return Format.get_goniometer(self)

    def get_detector(self, index=None):
        return Format.get_detector(self)

    def get_beam(self, index=None):
        return Format.get_beam(self)

    def get_scan(self, index=None):
        if index is None:
            return Format.get_scan(self)
        else:
            scan = Format.get_scan(self)
            return scan[index]

    def get_image_file(self, index=None):
        return Format.get_image_file(self)

    def get_raw_data(self, index):

        # Note MRC files use z, y, x ordering
        with mrcfile.mmap(self._image_file) as mrc:
            image = flex.double(mrc.data[index, ...].astype("double"))

        return image

    def _scan(self):
        """Scan model for this stack, filling out any unavailable items with
        dummy values"""

        alpha = self._header_dictionary.get("alphaTilt", 0.0)
        dalpha = self._header_dictionary.get("tiltPerImage", 1.0)
        exposure = self._header_dictionary.get("integrationTime", 0.0)

        oscillation = (alpha, dalpha)
        nframes = self.get_num_images()
        image_range = (1, nframes)
        epochs = [0] * nframes

        return self._scan_factory.make_scan(
            image_range, exposure, oscillation, epochs, deg=True
        )
