from __future__ import absolute_import, division, print_function

import sys
import uuid
from builtins import range

import h5py
import numpy

from iotbx.detectors.eiger import EIGERImage
from scitbx import matrix

from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.format.FormatPilatusHelpers import determine_eiger_mask
from dxtbx.format.FormatPilatusHelpers import get_vendortype_eiger as gv
from dxtbx.format.nexus import (
    BeamFactory,
    DataFactory,
    DataFactoryCache,
    DetectorFactory,
    GoniometerFactory,
    MaskFactory,
    NXdata,
    NXmxReader,
    generate_scan_model,
)


def find_entries(nx_file):
    """
    Find NXmx entries
    """
    if "entry" in nx_file:
        entry = nx_file["entry"]
        if "NX_class" in entry.attrs:
            if entry.attrs["NX_class"] == numpy.string_("NXentry"):
                if "definition" not in entry:
                    return entry
    return None


def is_eiger_nearly_nexus_file(filename):
    """
    A hacky function to check if this is an EIGER-flavoured nexus file
    """
    # Get the file handle
    with h5py.File(filename, "r") as handle:
        # Find the NXmx entries
        entry = find_entries(handle)
        if entry is not None:
            try:
                return (
                    numpy.string_("dectris eiger")
                    in entry["instrument"]["detector"]["description"][()].lower()
                )
            except KeyError:
                pass
        return False


class EigerNXmxFixer(object):
    """
    A hacky class to read an NXmx file
    """

    def __init__(self, input_filename, memory_mapped_name):
        # Copy the master file to the in memory handle
        handle_orig = h5py.File(input_filename, "r")
        handle = h5py.File(
            name=memory_mapped_name, driver="core", backing_store=False, mode="w"
        )
        handle_orig.copy("entry", handle)

        # Add some simple datasets
        def create_scalar(handle, path, dtype, value):
            dataset = handle.create_dataset(path, (), dtype=dtype)
            dataset[()] = value

        # Add NXmx definition
        create_scalar(handle["entry"], "definition", "S4", numpy.string_("NXmx"))

        # Add saturation value
        try:
            create_scalar(
                handle["entry/instrument/detector"],
                "saturation_value",
                "int32",
                handle[
                    "/entry/instrument/detector/detectorSpecific/countrate_correction_count_cutoff"
                ],
            )
        except Exception:
            create_scalar(
                handle["entry/instrument/detector"],
                "saturation_value",
                "int32",
                handle[
                    "entry/instrument/detector/detectorSpecific/detectorModule_000/countrate_correction_count_cutoff"
                ],
            )

        # Add detector type
        create_scalar(
            handle["entry/instrument/detector"], "type", "S4", numpy.string_("PIXEL")
        )

        # Move the beam
        # print "Copying /entry/instrument/beam to /entry/sample/beam"
        handle.copy("/entry/instrument/beam", "/entry/sample/beam")

        # Create detector module
        module_path = "/entry/instrument/detector/module"
        # print "Creating detector module %s" % (module_path)
        group = handle.create_group(module_path)
        group.attrs["NX_class"] = numpy.string_("NXdetector_module")

        # Add a module index
        create_scalar(group, "module_index", "int64", 0)

        # Create detector data origin
        dataset = group.create_dataset("data_origin", (2,), dtype="int32")
        dataset[0] = 0
        dataset[1] = 0

        # cope with badly structured chunk information i.e. many more data
        # entries than there are in real life...
        delete = []
        handle_orig_entry_properties = {}
        self.data_factory_cache = {}
        for k in sorted(handle_orig["/entry/data"]):
            try:
                handle_orig_entry = handle_orig["/entry/data/%s" % k]
                shape = handle_orig_entry.shape
            except KeyError:
                delete.append("/entry/data/%s" % k)
                continue
            handle_orig_entry_properties[k] = {
                "shape": shape,
                "length": len(handle_orig_entry),
                "filename": handle_orig_entry.file.filename,
            }
            self.data_factory_cache[k] = DataFactoryCache(
                shape=shape,
                ndim=handle_orig_entry.ndim,
                filename=handle_orig_entry.file.filename,
            )
        for d in delete:
            del handle[d]

        # Create detector data size
        dataset = group.create_dataset("data_size", (2,), dtype="int32")
        dataset[0] = handle_orig_entry_properties["data_000001"]["shape"][2]
        dataset[1] = handle_orig_entry_properties["data_000001"]["shape"][1]

        depends_on = "/entry/instrument/detector/transformations/translation"

        # Add fast_pixel_size dataset
        # print "Using /entry/instrument/detector/geometry/orientation/value as fast/slow pixel directions"
        fast_axis = handle["/entry/instrument/detector/geometry/orientation/value"][0:3]
        fast_axis = [
            fast_axis[0],
            fast_axis[1],
            -fast_axis[2],
        ]  # swap Z axis to align with Dectis/NeXus documentation
        slow_axis = handle["/entry/instrument/detector/geometry/orientation/value"][3:6]
        slow_axis = [
            slow_axis[0],
            slow_axis[1],
            -slow_axis[2],
        ]  # swap Z axis to align with Dectis/NeXus documentation
        create_scalar(
            group,
            "fast_pixel_direction",
            "float32",
            handle["/entry/instrument/detector/x_pixel_size"][()],
        )
        group["fast_pixel_direction"].attrs["transformation_type"] = numpy.string_(
            "translation"
        )
        group["fast_pixel_direction"].attrs["vector"] = fast_axis
        group["fast_pixel_direction"].attrs["offset"] = (0, 0, 0)
        group["fast_pixel_direction"].attrs["units"] = numpy.string_("m")
        group["fast_pixel_direction"].attrs["depends_on"] = numpy.string_(depends_on)

        # Add slow_pixel_size dataset
        create_scalar(
            group,
            "slow_pixel_direction",
            "float32",
            handle["/entry/instrument/detector/y_pixel_size"][()],
        )
        group["slow_pixel_direction"].attrs["transformation_type"] = numpy.string_(
            "translation"
        )
        group["slow_pixel_direction"].attrs["vector"] = slow_axis
        group["slow_pixel_direction"].attrs["offset"] = (0, 0, 0)
        group["slow_pixel_direction"].attrs["units"] = numpy.string_("m")
        group["slow_pixel_direction"].attrs["depends_on"] = numpy.string_(depends_on)

        # Add module offset dataset
        # print "Set module offset to be zero relative to detector"
        create_scalar(group, "module_offset", "float32", 0)
        group["module_offset"].attrs["transformation_type"] = numpy.string_(
            "translation"
        )
        group["module_offset"].attrs["vector"] = (0, 0, 0)
        group["module_offset"].attrs["offset"] = (0, 0, 0)
        group["module_offset"].attrs["units"] = numpy.string_("m")
        group["module_offset"].attrs["depends_on"] = numpy.string_(depends_on)

        # Create detector depends_on
        create_scalar(
            handle["/entry/instrument/detector"],
            "depends_on",
            "S%d" % len(depends_on),
            numpy.string_(depends_on),
        )

        # Add detector position
        # print "Using /entry/instrument/detector/geometry/translation/distances as transformation"
        detector_offset_vector = handle[
            "/entry/instrument/detector/geometry/translation/distances"
        ][()]
        # swap Z axis to align with Dectis/NeXus documentation
        detector_offset_vector = matrix.col(
            (
                detector_offset_vector[0],
                detector_offset_vector[1],
                -detector_offset_vector[2],
            )
        )
        group = handle.create_group("/entry/instrument/detector/transformations")
        group.attrs["NX_class"] = numpy.string_("NXtransformations")
        create_scalar(group, "translation", "float32", detector_offset_vector.length())
        group["translation"].attrs["transformation_type"] = numpy.string_("translation")
        if detector_offset_vector.length() > 0:
            group["translation"].attrs["vector"] = detector_offset_vector.normalize()
        else:
            group["translation"].attrs["vector"] = detector_offset_vector
        group["translation"].attrs["offset"] = 0
        group["translation"].attrs["units"] = numpy.string_("m")
        group["translation"].attrs["depends_on"] = numpy.string_(".")

        # Create goniometer transformations if not found
        if "/entry/sample/transformations" not in handle:
            # print "Creating group /entry/sample/transformation"
            group = handle.create_group("/entry/sample/transformations")
            group.attrs["NX_class"] = numpy.string_("NXtransformations")
        else:
            group = handle["/entry/sample/transformations"]

        # check for incomplete omega definitions dirty hack...
        if "omega" in group:
            try:
                group["omega"][()]
            except AttributeError:
                del group["omega"]

        if "omega" not in group:
            # In here assume goniometer axis is 1,0,0 unless (i) specified somewhere
            # we can know or (ii) a known special case. For (i) for this instrument
            # listed here this is at
            #
            # /entry/sample/transformations/omega->vector
            #
            # which will join the special case list once this is properly resolved
            # as this corrently returns 0, -1, 0 so needs transforming...
            #
            # special cases:
            # E-32-0105 - Max IV, vertical axis

            try:
                key = handle["/entry/instrument/detector/detector_number"][()]
                default_axis = {"E-32-0105": (0, 1, 0)}[key]
            except KeyError:
                default_axis = (-1, 0, 0)

            num_images = 0
            for name in sorted(handle["/entry/data"]):
                num_images += handle_orig_entry_properties[name]["length"]
            dataset = group.create_dataset("omega", (num_images,), dtype="float32")
            dataset.attrs["units"] = numpy.string_("degree")
            dataset.attrs["transformation_type"] = numpy.string_("rotation")
            dataset.attrs["vector"] = default_axis
            dataset.attrs["offset"] = 0
            dataset.attrs["depends_on"] = numpy.string_(".")
            omega_range_average = handle[
                "/entry/sample/goniometer/omega_range_average"
            ][()]
            omega_range_average = int(omega_range_average * 100 + 0.5) / 100.0
            for i in range(num_images):
                angle = omega_range_average * i
                dataset[i] = angle
        else:
            dataset = group["omega"]

        if "depends_on" not in handle["/entry/sample"]:
            # Create sample depends_on
            create_scalar(
                handle["/entry/sample"],
                "depends_on",
                "S%d" % len(dataset.name),
                numpy.string_(dataset.name),
            )

        # Change relative paths to absolute paths
        for name in handle["/entry/data"]:
            del handle["entry/data"][name]
            filename = handle_orig_entry_properties[name]["filename"]
            handle["entry/data"][name] = h5py.ExternalLink(filename, "entry/data/data")
            handle["entry/data"]["_filename_" + name] = filename  # Store file names

        self.handle = handle
        self.handle_orig = handle_orig


class FormatHDF5EigerNearlyNexus(FormatHDF5):
    @staticmethod
    def understand(image_file):
        try:
            return is_eiger_nearly_nexus_file(image_file)
        except IOError:
            return False

    def _start(self):
        # Read the file structure
        temp_file = "tmp_master_%s.nxs" % uuid.uuid1().hex
        fixer = EigerNXmxFixer(self._image_file, temp_file)
        reader = NXmxReader(handle=fixer.handle)

        # Only support 1 set of models at the moment
        assert len(reader.entries) == 1, "Currently only supports 1 NXmx entry"
        assert len(reader.entries[0].data) == 1, "Currently only supports 1 NXdata"
        assert (
            len(reader.entries[0].instruments) == 1
        ), "Currently only supports 1 NXinstrument"
        assert len(reader.entries[0].samples) == 1, "Currently only supports 1 NXsample"
        assert (
            len(reader.entries[0].instruments[0].detectors) == 1
        ), "Currently only supports 1 NXdetector"
        assert (
            len(reader.entries[0].instruments[0].detectors[0].modules) == 1
        ), "Currently only supports 1 NXdetector_module"
        assert (
            len(reader.entries[0].samples[0].beams) == 1
        ), "Currently only supports 1 NXbeam"

        # Get the NXmx model objects
        entry = reader.entries[0]
        self.instrument = instrument = entry.instruments[0]
        detector = instrument.detectors[0]
        sample = entry.samples[0]
        beam = sample.beams[0]

        # Use data from original master file
        data = NXdata(fixer.handle_orig[entry.data[0].handle.name])

        # Construct the models
        self._beam_factory = BeamFactory(beam)
        self._beam_factory.load_model(0)
        self._detector_model = DetectorFactory(detector, self._beam_factory.model).model
        self._goniometer_model = GoniometerFactory(sample).model
        self._scan_model = generate_scan_model(sample, detector)
        self._raw_data = DataFactory(data, cached_information=fixer.data_factory_cache)

        # update model for masking Eiger detectors
        for f0, f1, s0, s1 in determine_eiger_mask(self._detector_model):
            self._detector_model[0].add_mask(f0 - 1, s0 - 1, f1, s1)

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
        return self._scan()[index]

    def get_raw_data(self, index):
        return self._raw_data[index]

    def get_mask(self, index=None, goniometer=None):
        return MaskFactory(self.instrument.detectors, index).mask

    def get_num_images(self):
        scan = self._scan()
        if isinstance(scan, list):
            return sum(s.get_num_images() for s in scan)
        return scan.get_num_images()

    def get_image_file(self, index=None):
        return self._image_file

    def detectorbase_start(self, index=0):
        self.detectorbase = EIGERImage(self._image_file, index=index)
        self.detectorbase.readHeader(dxtbx_instance=self)

        def model_get_raw_data(ptr, index):
            return self.get_raw_data(index)

        self.detectorbase.get_raw_data_callback = model_get_raw_data

    def get_detectorbase(self, index=0):
        self.detectorbase_start(index)
        return self.detectorbase

    def get_vendortype(self):
        return gv(self.get_detector())


if __name__ == "__main__":
    f = FormatHDF5EigerNearlyNexus(sys.argv[1])
    for i in range(10):
        print(f.get_raw_data(i))
