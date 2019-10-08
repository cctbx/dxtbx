from __future__ import absolute_import, division, print_function

import h5py
from dials.array_family import flex
from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.format.nexus import (
    BeamFactory,
    DataFactory,
    DetectorFactory,
    MaskFactory,
    NXmxReader,
)
from dxtbx.model import Scan


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
        self._beam_model = BeamFactory(beam).model

        self._setup_gonio_and_scan(sample, detector)

        if self._scan_model:
            array_range = self._scan_model.get_array_range()
            num_images = array_range[1] - array_range[0]
        else:
            num_images = 0

        self._detector_model = DetectorFactory(detector, self._beam_model).model
        self._raw_data = DataFactory(data, max_size=num_images).model

    def _setup_gonio_and_scan(self, sample, detector):
        # scan

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
        self._goniometer_model = self._goniometer_factory.single_axis()

    def _end(self):
        return

    def _goniometer(self):
        return self._goniometer_model

    def _detector(self):
        return self._detector_model

    def _beam(self, index=None):
        if index is None:
            index = 0

        entry = self._reader.entries[0]
        sample = entry.samples[0]
        beam = sample.beams[0]

        self._beam_model = BeamFactory(beam, index).model
        return self._beam_model

    def _scan(self):
        return self._scan_model

    def get_goniometer(self, index=None):
        return self._goniometer()

    def get_detector(self, index=None):
        return self._detector()

    def get_beam(self, index=None):
        return self._beam()

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
