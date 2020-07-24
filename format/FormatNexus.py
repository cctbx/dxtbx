from __future__ import absolute_import, division, print_function

import sys

import h5py

from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.format.FormatMultiImageLazy import FormatMultiImageLazy
from dxtbx.format.FormatStill import FormatStill
from dxtbx.format.nexus import (
    BeamFactory,
    DataFactory,
    DetectorFactory,
    DetectorFactoryFromGroup,
    GoniometerFactory,
    MaskFactory,
    NXmxReader,
    detectorgroupdatafactory,
    generate_scan_model,
    is_nexus_file,
)


class FormatNexus(FormatHDF5):
    @staticmethod
    def understand(image_file):
        try:
            return is_nexus_file(image_file)
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
            or len(reader.entries[0].instruments[0].beams) == 1
        ), "Currently only supports 1 NXbeam"

        # Get the NXmx model objects
        entry = reader.entries[0]
        self.instrument = instrument = entry.instruments[0]
        detector = instrument.detectors[0]
        sample = entry.samples[0]
        beam = sample.beams[0] if sample.beams else instrument.beams[0]
        data = entry.data[0]

        # Construct the models
        self._beam_factory = BeamFactory(beam)
        self._beam_factory.load_model(0)

        self._setup_gonio_and_scan(sample, detector)

        if self._scan_model:
            array_range = self._scan_model.get_array_range()
            num_images = array_range[1] - array_range[0]
        else:
            num_images = 0

        if len(instrument.detector_groups) == 0:
            assert (
                len(reader.entries[0].instruments[0].detectors) == 1
            ), "Currently only supports 1 NXdetector unless in a detector group"
            assert (
                len(reader.entries[0].instruments[0].detectors[0].modules) == 1
            ), "Currently only supports 1 NXdetector_module unless in a detector group"

            self._raw_data = DataFactory(data, max_size=num_images)
            self._detector_model = DetectorFactory(
                detector, self._beam_factory.model, shape=self._raw_data.shape()
            ).model
        else:
            self._raw_data = detectorgroupdatafactory(data, instrument)
            self._detector_model = DetectorFactoryFromGroup(
                instrument, self._beam_factory.model
            ).model

    def _setup_gonio_and_scan(self, sample, detector):
        """Set up rotation-specific models"""
        self._goniometer_model = GoniometerFactory(sample).model
        self._scan_model = generate_scan_model(sample, detector)

    def _end(self):
        return

    def _goniometer(self):
        return self._goniometer_model

    def _detector(self):
        return self._detector_model

    def _beam(self, index=None):
        self._beam_model = self._beam_factory.load_model(index)
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


class FormatNexusStill(FormatMultiImageLazy, FormatNexus, FormatStill):
    @staticmethod
    def understand(image_file):
        is_nexus_still = False
        try:
            from dxtbx.format.nexus import find_entries, find_class

            # Get the file handle
            with h5py.File(image_file, "r") as handle:
                if "/entry/sample/goniometer/omega_increment" in handle:
                    return False

                for entry in find_entries(handle, "/"):
                    for sample in find_class(entry, "NXsample"):
                        if "depends_on" not in sample:
                            is_nexus_still = True
        except IOError:
            return False
        return is_nexus_still

    def _setup_gonio_and_scan(self, sample, detector):
        """No rotation-specific models for stills"""
        self._goniometer_model = None
        self._scan_model = None

    def get_num_images(self):
        return len(self._raw_data)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        if FormatNexus.understand(arg):

            format_instance = FormatNexus(arg)

            beam = format_instance.get_beam()
            detector = format_instance.get_detector()
            goniometer = format_instance.get_goniometer()
            scan = format_instance.get_scan()

            iset = FormatNexus.get_imageset(arg)
            print(beam)
            print(detector)
            print(goniometer)
            print(scan)

            print(len(iset))
