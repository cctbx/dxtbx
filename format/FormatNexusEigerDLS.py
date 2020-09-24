from __future__ import absolute_import, division, print_function

import ast
import os

import h5py
import numpy

from dxtbx.format.FormatNexusEiger import FormatNexusEiger
from dxtbx.format.nexus import (
    BeamFactory,
    DataFactory,
    DetectorFactory,
    DetectorFactoryFromGroup,
    NXmxReader,
    detectorgroupdatafactory,
)

try:
    from dxtbx_format_nexus_ext import (
        image_as_flex_double,
        image_as_flex_float,
        image_as_flex_int,
    )
except ImportError:
    # Workaround for psana build, which doesn't link HDF5 properly
    if "SIT_ROOT" not in os.environ:
        raise


def image_as_flex(dataset, index):
    if numpy.issubdtype(dataset.dtype, numpy.integer):
        return image_as_flex_int(dataset.id.id, index)
    else:
        assert numpy.issubdtype(dataset.dtype, numpy.floating)
        if dataset.dtype in [
            numpy.half,
            numpy.single,
            numpy.float_,
            numpy.float16,
            numpy.float32,
        ]:
            return image_as_flex_float(dataset.id.id, index)
        elif dataset.dtype in [
            numpy.double,
            numpy.longfloat,
            numpy.float64,
            numpy.float96,
            numpy.float128,
        ]:
            return image_as_flex_double(dataset.id.id, index)
        else:
            assert False, "unknown floating data type (%s)" % str(dataset.dtype)


class DLSEigerDataFactory(DataFactory):
    def __getitem__(self, index):
        d = self._lookup[index]
        i = index - self._offset[d]

        # a lock-free most-recently-used file handle cache based on
        # immutability of python tuples
        cached_handle = self._cache
        if self._datasets[d].file == cached_handle[0]:
            data = cached_handle[1]
        else:
            data = self._datasets[d].accessor()
            self._cache = (self._datasets[d].file, data)

        return image_as_flex(dataset=data, index=i)


def get_count_limit_from_meta(meta_file_name):
    with h5py.File(meta_file_name, "r") as f:

        config = f["/config"][()]
        config_data = ast.literal_eval(config.decode("utf-8"))

    return config_data["countrate_correction_count_cutoff"]


def get_bit_depth_from_meta(meta_file_name):
    with h5py.File(meta_file_name, "r") as f:
        config = f["/config"][()]
        config_data = ast.literal_eval(config.decode("utf-8"))

    return config_data["bit_depth_image"]


def find_meta_filename(master_like):
    meta_filename = None
    f = h5py.File(master_like, "r")

    def _local_visit(name):
        obj = f[name]
        if not hasattr(obj, "keys"):
            return None
        for k in obj.keys():
            kclass = obj.get(k, getlink=True, getclass=True)
            if kclass is h5py._hl.group.ExternalLink:
                kfile = obj.get(k, getlink=True).filename
                if kfile.split(".")[0].endswith("meta"):
                    return kfile

    master_dir = os.path.split(master_like)[0]
    meta_filename = f.visit(_local_visit)

    return os.path.join(master_dir, meta_filename)


class FormatNexusEigerDLS(FormatNexusEiger):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            name = FormatNexusEigerDLS.get_instrument_name(handle)
            if name is None:
                return False
            if name.lower() in (b"i03", b"i04", b"i24", b"vmxi"):
                return True

        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super(FormatNexusEigerDLS, self).__init__(image_file, **kwargs)
        try:
            self._meta = find_meta_filename(image_file)
            self._bit_depth_image = get_bit_depth_from_meta(self._meta)
        except Exception:
            self._bit_depth_image = 16

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

            self._raw_data = DLSEigerDataFactory(data, max_size=num_images)
            self._detector_model = DetectorFactory(
                detector, self._beam_factory.model, shape=self._raw_data.shape()
            ).model
        else:
            self._raw_data = detectorgroupdatafactory(data, instrument)
            self._detector_model = DetectorFactoryFromGroup(
                instrument, self._beam_factory.model
            ).model

    def get_detector(self, index=None):
        # workaround for https://jira.diamond.ac.uk/browse/I03-365
        # read the count limit from the meta file - if anything goes
        # wrong, do nothing

        detector = self._detector()

        try:
            limit = get_count_limit_from_meta(self._meta)
            assert limit > 0

        except Exception:
            for panel in detector:
                trusted = panel.get_trusted_range()
                panel.set_trusted_range((-1, trusted[1]))

        else:
            for panel in detector:
                panel.set_trusted_range((-1, limit))

        return detector

    def get_raw_data(self, index):
        data = self._raw_data[index]
        if self._bit_depth_image:
            # if 32 bit then it is a signed int, I think if 8, 16 then it is
            # unsigned with the highest two values assigned as masking values
            if self._bit_depth_image == 32:
                top = 2 ** 31
            else:
                top = 2 ** self._bit_depth_image
            d1d = data.as_1d()
            d1d.set_selected(d1d == top - 1, -1)
            d1d.set_selected(d1d == top - 2, -2)

        return data
