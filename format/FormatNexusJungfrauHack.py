from __future__ import absolute_import, division, print_function

import string

import h5py
import numpy

from cctbx.eltbx import attenuation_coefficient
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.format.nexus import (
    BeamFactory,
    DataFactory,
    MaskFactory,
    NXmxReader,
    convert_units,
)
from dxtbx.model import Detector, Panel, ParallaxCorrectedPxMmStrategy, Scan


def clean_string(input):
    if hasattr(string, "letters"):
        letters = string.letters
    else:
        letters = string.ascii_letters
    return "".join(i for i in input if i in letters)


class FormatNexusJungfrauHack(FormatNexus):
    @staticmethod
    def understand(image_file):
        try:
            with h5py.File(image_file, "r") as handle:
                return "/entry/instrument/JF1M" in handle
        except IOError:
            return False

    def _start(self):

        # Read the file structure
        self._reader = reader = NXmxReader(self._image_file)

        # Only support 1 set of models at the moment
        assert len(reader.entries) == 1, "Currently only supports 1 NXmx entry"
        assert len(reader.entries[0].data) == 1, "Currently only supports 1 NXdata"
        assert (
            len(reader.entries[0].instruments) == 1
        ), "Currently only supports 1 NXinstrument"
        assert len(reader.entries[0].samples) == 1, "Currently only supports 1 NXsample"
        assert (
            len(reader.entries[0].samples[0].beams) == 1
        ), "Currently only supports 1 NXbeam"

        # Get the NXmx model objects
        entry = reader.entries[0]
        self.instrument = instrument = entry.instruments[0]
        detector = instrument.detectors[0]
        sample = entry.samples[0]
        beam = sample.beams[0]
        data = entry.data[0]

        # Construct the models
        self._beam_factory = BeamFactory(beam).model
        self._beam_factory.load_model(0)

        self._setup_detector(detector, self._beam_factory.model)
        self._setup_gonio_and_scan(sample, detector)

        if self._scan_model:
            array_range = self._scan_model.get_array_range()
            num_images = array_range[1] - array_range[0]
        else:
            num_images = 0

        self._raw_data = DataFactory(data, max_size=num_images)

    def _setup_detector(self, detector, beam):
        nx_detector = detector.handle
        nx_module = detector.modules[0].handle

        # Get the detector name and type
        if "type" in nx_detector:
            detector_type = str(nx_detector["type"][()])
        else:
            detector_type = "unknown"
        detector_name = str(nx_detector.name)

        # Get the trusted range of pixel values
        if "saturation_value" in nx_detector:
            trusted_range = (-1, float(nx_detector["saturation_value"][()]))
        else:
            trusted_range = (-1, 99999999)

        # Get the detector thickness
        thickness = nx_detector["sensor_thickness"]
        thickness_value = float(thickness[()])
        thickness_units = thickness.attrs["units"]
        thickness_value = float(convert_units(thickness_value, thickness_units, "mm"))

        # Get the detector material
        material = nx_detector["sensor_material"][()]
        if hasattr(material, "shape"):
            material = "".join(m.decode() for m in material)
        detector_material = clean_string(str(material))
        material = {
            "Si": "Si",
            numpy.string_("Si"): "Si",
            numpy.string_("Silicon"): "Si",
            numpy.string_("Sillicon"): "Si",
            numpy.string_("CdTe"): "CdTe",
            numpy.string_("GaAs"): "GaAs",
        }.get(detector_material)
        if not material:
            raise RuntimeError("Unknown material: %s" % detector_material)

        try:
            x_pixel = nx_detector["x_pixel_size"][()] * 1000.0
            y_pixel = nx_detector["y_pixel_size"][()] * 1000.0

            legacy_beam_x = float(x_pixel * nx_detector["beam_center_x"][()])
            legacy_beam_y = float(y_pixel * nx_detector["beam_center_y"][()])
            legacy_distance = float(1000 * nx_detector["detector_distance"][()])
        except KeyError:
            legacy_beam_x = 0
            legacy_beam_y = 0

        # Get the fast pixel size and vector
        fast_pixel_direction = nx_module["fast_pixel_direction"]
        fast_pixel_direction_value = float(fast_pixel_direction[()])
        fast_pixel_direction_units = fast_pixel_direction.attrs["units"]
        fast_pixel_direction_vector = fast_pixel_direction.attrs["vector"]
        fast_pixel_direction_value = convert_units(
            fast_pixel_direction_value, fast_pixel_direction_units, "mm"
        )
        fast_axis = matrix.col(fast_pixel_direction_vector).normalize()

        # Get the slow pixel size and vector
        slow_pixel_direction = nx_module["slow_pixel_direction"]
        slow_pixel_direction_value = float(slow_pixel_direction[()])
        slow_pixel_direction_units = slow_pixel_direction.attrs["units"]
        slow_pixel_direction_vector = slow_pixel_direction.attrs["vector"]
        slow_pixel_direction_value = convert_units(
            slow_pixel_direction_value, slow_pixel_direction_units, "mm"
        )
        slow_axis = matrix.col(slow_pixel_direction_vector).normalize()

        origin = matrix.col((0.0, 0.0, -legacy_distance))
        if origin.elems[0] == 0.0 and origin.elems[1] == 0:
            origin = -(origin + legacy_beam_x * fast_axis + legacy_beam_y * slow_axis)

        # Change of basis to convert from NeXus to IUCr/ImageCIF convention
        cob = matrix.sqr((-1, 0, 0, 0, 1, 0, 0, 0, -1))
        origin = cob * matrix.col(origin)
        fast_axis = cob * fast_axis
        slow_axis = cob * slow_axis

        # Ensure that fast and slow axis are orthogonal
        normal = fast_axis.cross(slow_axis)
        slow_axis = -fast_axis.cross(normal)

        # Compute the attenuation coefficient.
        # This will fail for undefined composite materials
        # mu_at_angstrom returns cm^-1, but need mu in mm^-1
        table = attenuation_coefficient.get_table(material)
        wavelength = beam.get_wavelength()
        mu = table.mu_at_angstrom(wavelength) / 10.0

        # Construct the detector model
        pixel_size = (fast_pixel_direction_value, slow_pixel_direction_value)
        # image size stored slow to fast but dxtbx needs fast to slow
        image_size = tuple(int(x) for x in reversed(nx_module["data_size"][-2:]))
        image_size = (1030, 1064)

        self.model = Detector()
        self.model.add_panel(
            Panel(
                detector_type,
                detector_name,
                tuple(fast_axis),
                tuple(slow_axis),
                tuple(origin),
                pixel_size,
                image_size,
                trusted_range,
                thickness_value,
                material,
                mu,
            )
        )

        # Set the parallax correction
        for panel in self.model:
            panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, thickness_value))
            panel.set_type("SENSOR_PAD")

        self._detector_model = self.model

    def _setup_gonio_and_scan(self, sample, detector):
        with h5py.File(self._image_file, "r") as handle:
            phi = handle["/entry/sample/goniometer/omega"][()]
        image_range = (1, len(phi))
        oscillation = (float(phi[0]), float(phi[1] - phi[0]))

        # Get the exposure time
        num_images = len(phi)
        exposure_time = flex.double(num_images, 0)
        epochs = flex.double(num_images, 0.0)
        for i in range(1, len(epochs)):
            epochs[i] = epochs[i - 1] + exposure_time[i - 1]

        self._scan_model = Scan(image_range, oscillation, exposure_time, epochs)
        self._goniometer_model = self._goniometer_factory.known_axis((-1, 0, 0))

    def _end(self):
        return

    def _goniometer(self):
        return self._goniometer_model

    def _detector(self):
        return self._detector_model

    def _beam(self, index=None):
        self._beam_factory.load_model(index)
        self._beam_model = self._beam_factory.model
        return self._beam_model

    def _scan(self):
        return self._scan_model

    def get_goniometer(self, index=None):
        return self._goniometer()

    def get_detector(self, index=None):
        return self._detector()

    def get_beam(self, index=None):
        return self._beam(index)

    def get_scan(self, index=None):
        if index is None:
            return self._scan()
        scan = self._scan()
        if scan is not None:
            return scan[index]
        return scan

    def get_raw_data(self, index):
        return self._raw_data[index]

    def get_mask(self, index=None, goniometer=None):
        return MaskFactory(self.instrument.detectors, index).mask

    def get_num_images(self):
        if self._scan() is not None:
            return self._scan().get_num_images()
        return len(self._raw_data)

    def get_image_file(self, index=None):
        return self._image_file

    def get_detectorbase(self, index=None):
        raise NotImplementedError

    @staticmethod
    def get_instrument_name(handle):
        if "short_name" in handle["/entry/instrument"].attrs:
            name = handle["/entry/instrument"].attrs["short_name"]
        elif "/entry/instrument/name" in handle:
            if "short_name" in handle["/entry/instrument/name"].attrs:
                name = handle["/entry/instrument/name"].attrs["short_name"]
            else:
                name = handle["/entry/instrument/name"][()]
        else:
            name = None
        return name
