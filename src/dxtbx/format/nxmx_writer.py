"""
Note, scans and gonios not supported here. This writer essentially writes still images

Example to write 10 cbfs to a single NXmx file:
writer = NXmxWriter("*.cbf")
writer.write_nxmx("example.h5")
"""

from __future__ import annotations

import datetime
import os
import sys

import h5py
import numpy as np

from scitbx.array_family import flex
from xfel.cftbx.detector.cspad_cbf_tbx import basis, angle_and_axis
from scitbx import matrix

from dxtbx import flumpy


# Conversion from the imgCIF coordinate system conventionally used by dxtbx to
# the McStas coordinate system as used by NeXus:
#   https://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Caxis.html
#   https://manual.nexusformat.org/design.html#design-coordinatesystem
IMGCIF_TO_MCSTAS = matrix.diag([-1, 1, -1])


class NXmxWriter:
    """Class for writing NXmx NeXus files from any dxtbx-supported format class"""

    def __init__(self, params, experiments):
        self.experiments = experiments
        self.params = params

    def get_metrology_dict(self, index=None):
        """Build a metrology dictionary.  This dictionary maps hierarchy keys to basis
        objects. A hierarchy key looks like this (0,1,2), where the entries are
        levels in a hierarchy and the numbers refer to a panel or group within that
        level"""
        metro = {}

        def recursive_setup_dict(panelgroup, key):
            metro[key] = basis(panelgroup=panelgroup)
            if panelgroup.is_panel():
                return
            for i, child in enumerate(panelgroup):
                childkey = tuple(list(key) + [i])
                recursive_setup_dict(child, childkey)

        if index is None:
            detector = self.experiments.detectors()[0]
        else:
            detector = self.experiments.detectors()[index]

        recursive_setup_dict(detector.hierarchy(), (0,))
        return metro

    def _create_scalar(self, handle, path, dtype, value):
        dataset = handle.create_dataset(path, (), dtype=dtype)
        dataset[()] = value

    def create_vector(self, handle, name, value, **attributes):
        handle.create_dataset(name, (1,), data=[value], dtype="f")
        for key, attribute in attributes.items():
            handle[name].attrs[key] = attribute

    def add_frame_shift(self, handle, basis):
        """Add an axis representing a frame shift (a rotation axis with an offset)"""
        angle, axis = angle_and_axis(basis)

        if angle == 0:
            axis = matrix.col((0, 0, 1))

        if basis.include_translation:
            translation = basis.translation
        else:
            translation = matrix.col((0, 0, 0))

        axis = IMGCIF_TO_MCSTAS * axis
        translation = IMGCIF_TO_MCSTAS * translation

        self.create_vector(
            handle,
            os.path.basename(basis.name),
            angle,
            depends_on=basis.depends_on,
            equipment="detector",
            equipment_component=basis.equipment_component,
            transformation_type="rotation",
            units="deg",
            vector=axis,
            offset=translation,
            offset_units="mm",
        )

    def get_h5py_handle(self):
        """
        Hierarchical structure of master nexus file. Format information available here
        http://download.nexusformat.org/sphinx/classes/base_classes/NXdetector_module.html#nxdetector-module
        --> entry
          --> data
          --> definition (leaf)
          --> instrument
          --> sample
        """
        # set up the metrology dictionary to include axis names, pixel sizes, and so forth
        assert (
            len(self.experiments.detectors()) == 1
        ), "Multiple detectors not supported"
        detector = self.experiments.detectors()[0]
        metro = self.get_metrology_dict()

        def panel_group_from_key(key):
            # Find the node that a hierarchy key refers to
            if len(key) == 1:
                assert key[0] == 0
                return detector.hierarchy()

            node = detector.hierarchy()
            for i in key[1:]:
                node = node[i]
            return node

        def level_string(key):
            # Example for key (0,1,2). "L0M0_L1M1_L2M2", where L is level and M is module
            return "_".join(["L%dM%d" % (l, m) for l, m in enumerate(key)])

        # set up the metrology dictionary to include axis names, pixel sizes, and so forth
        panelkeys = []
        panelnames = []

        def recursive_setup_basis_dict(key, parent_name="", panel_id=0):
            # Set up CBF axis names, including equipment components and depends_on chains
            basis = metro[key]
            node = panel_group_from_key(key)
            nodename = level_string(key)
            if basis.name:
                dname = basis.name
            else:
                basis.name = dname = "AXIS_" + nodename

            if node.is_panel():
                panelname = "PANEL_%d" % panel_id
                panelkeys.append(key)
                panelnames.append(panelname)
                panel_id += 1

            if len(key) == 1:
                assert key == (0,)  # only one detector allowed for now
                basis.depends_on = "."
            else:
                basis.depends_on = parent_name

            basis.equipment_component = "detector_level_%d" % (len(key) - 1)
            basis.axis_name = dname
            if not node.is_panel():
                for c, child in enumerate(node):
                    panel_id = recursive_setup_basis_dict(
                        tuple(list(key) + [c]), dname, panel_id
                    )
            return panel_id

        recursive_setup_basis_dict((0,))

        output_file_name = (
            self.params.output_file
            if self.params.output_file is not None
            else "converted.h5"
        )
        f = h5py.File(output_file_name, "w")
        f.attrs["NX_class"] = "NXroot"
        f.attrs["file_name"] = os.path.basename(output_file_name)
        f.attrs["file_time"] = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        f.attrs["HDF5_Version"] = h5py.version.hdf5_version
        entry = f.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        if self.params.nexus_details.start_time:
            entry["start_time"] = self.params.nexus_details.start_time
        if self.params.nexus_details.end_time:
            entry["end_time"] = self.params.nexus_details.end_time
        if self.params.nexus_details.end_time_estimated:
            entry["end_time_estimated"] = self.params.nexus_details.end_time_estimated

        # --> definition
        self._create_scalar(entry, "definition", "S4", np.string_("NXmx"))

        # --> sample
        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        if self.params.nexus_details.sample_name:
            sample["name"] = self.params.nexus_details.sample_name
        sample["depends_on"] = "."  # This script does not support scans/gonios
        # --> source
        source = entry.create_group("source")
        source.attrs["NX_class"] = "NXsource"
        source["name"] = self.params.nexus_details.source_name
        source["name"].attrs["short_name"] = self.params.nexus_details.source_short_name
        # --> instrument
        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"
        instrument["name"] = self.params.nexus_details.instrument_name
        instrument["name"].attrs[
            "short_name"
        ] = self.params.nexus_details.instrument_short_name
        beam = instrument.create_group("beam")
        beam.attrs["NX_class"] = "NXbeam"
        if self.params.nexus_details.total_flux:
            self._create_scalar(
                beam, "total_flux", "f", self.params.nexus_details.total_flux
            )
            beam["total_flux"].attrs["units"] = "Hz"

        if self.params.wavelength:
            beam.create_dataset(
                "incident_wavelength", (1,), data=self.params.wavelength, dtype="f8"
            )
        else:
            if len(self.experiments.beams()) > 1:
                wavelengths = np.array(
                    [beam.get_wavelength() for beam in self.experiments.beams()],
                    dtype="f8",
                )
                beam.create_dataset(
                    "incident_wavelength",
                    wavelengths.shape,
                    data=wavelengths,
                    dtype=wavelengths.dtype,
                )
            else:
                wavelength = self.experiments.beams()[0].get_wavelength()
                beam.create_dataset(
                    "incident_wavelength", (1,), data=wavelength, dtype="f8"
                )
        beam["incident_wavelength"].attrs["units"] = "angstrom"

        det_group = instrument.create_group("detector")
        det_group.attrs["NX_class"] = "NXdetector_group"

        det_group.create_dataset("group_index", data=list(range(1, 3)), dtype="i")
        data = [np.string_("detector"), np.string_("ELE_D0")]
        det_group.create_dataset("group_names", (2,), data=data, dtype="S12")
        det_group.create_dataset("group_parent", (2,), data=[-1, 1], dtype="i")
        det_group.create_dataset("group_type", (2,), data=[1, 2], dtype="i")
        det = instrument.create_group("ELE_D0")
        det.attrs["NX_class"] = "NXdetector"
        det["description"] = "Detector converted from DIALS models"
        det["depends_on"] = "/entry/instrument/ELE_D0/transformations/AXIS_RAIL"
        det["gain_setting"] = "auto"
        assert len(set([p.get_material() for p in detector])) == 1
        assert len(set([p.get_thickness() for p in detector])) == 1
        det["sensor_material"] = detector[0].get_material()
        self._create_scalar(
            det, "sensor_thickness", "f", detector[0].get_thickness() * 1000
        )
        self._create_scalar(det, "bit_depth_readout", "i", 16)  # XXX
        self._create_scalar(det, "count_time", "f", 10)  # XXX
        self._create_scalar(det, "frame_time", "f", 40)  # XXX
        det["sensor_thickness"].attrs["units"] = "microns"
        det["count_time"].attrs["units"] = "us"
        det["frame_time"].attrs["units"] = "us"
        transformations = det.create_group("transformations")
        transformations.attrs["NX_class"] = "NXtransformations"

        if self.params.trusted_range is None:
            assert len(set([p.get_trusted_range() for p in detector])) == 1
            trusted_min, trusted_max = detector[0].get_trusted_range()
        else:
            trusted_min, trusted_max = self.params.trusted_range
        # DIALS definitions match up with NXmx
        det.create_dataset("underload_value", (1,), data=[trusted_min], dtype="int32")
        det.create_dataset("saturation_value", (1,), data=[trusted_max], dtype="int32")

        def find_panel_id(panel):
            for i in range(len(detector)):
                if detector[i].get_name() == panel.get_name():
                    return i

        # Create a series of axis describing frame shifts from each level of the detector to the next
        for key in sorted(metro):
            basis = metro[key]
            self.add_frame_shift(transformations, basis)

            node = panel_group_from_key(key)

            if node.is_panel():
                aname = level_string(key)
                fast = IMGCIF_TO_MCSTAS * matrix.col((1.0, 0.0, 0.0))
                slow = IMGCIF_TO_MCSTAS * matrix.col((0.0, 1.0, 0.0))

                module = det.create_group(aname)
                module.attrs["NX_class"] = "NXdetector_module"

                if len(detector) > 1:
                    panel_id = find_panel_id(node)
                    module.create_dataset(
                        "data_origin", (3,), data=[panel_id, 0, 0], dtype="i"
                    )
                    module.create_dataset(
                        "data_size",
                        (3,),
                        data=[1, *list(reversed(node.get_image_size()))],
                        dtype="i",
                    )
                else:
                    module.create_dataset("data_origin", (2,), data=[0, 0], dtype="i")
                    module.create_dataset(
                        "data_size",
                        (2,),
                        data=list(reversed(node.get_image_size())),
                        dtype="i",
                    )

                self.create_vector(
                    module,
                    "fast_pixel_direction",
                    node.get_pixel_size()[1],
                    depends_on=transformations.name
                    + "/"
                    + os.path.basename(basis.name),
                    transformation_type="translation",
                    units="mm",
                    vector=fast,
                    offset=(0.0, 0.0, 0.0),
                    offset_units="mm",
                )
                self.create_vector(
                    module,
                    "slow_pixel_direction",
                    node.get_pixel_size()[0],
                    depends_on=transformations.name
                    + "/"
                    + os.path.basename(basis.name),
                    transformation_type="translation",
                    units="mm",
                    vector=slow,
                    offset=(0.0, 0.0, 0.0),
                    offset_units="mm",
                )

        return f

    def append_frame(self, handle, index=None, data=None):
        """
        Given a h5py handle, append the data
        """
        if data is None:
            if index is None:
                data = self.experiments[0].imageset[0]
            else:
                data = self.experiments[0].imageset[index]
        if not isinstance(data, tuple):
            data = (data,)

        if len(data) > 1:
            assert len(set([d.focus() for d in data])) == 1
            shape = len(data), data[0].focus()[0], data[0].focus()[1]
        else:
            shape = data[0].focus()[0], data[0].focus()[1]

        dataisint = []
        for panel_data in data:
            assert len(panel_data.focus()) == 2
            if isinstance(panel_data, flex.int):
                dataisint.append(True)
            elif isinstance(panel_data, flex.double):
                dataisint.append(False)
            else:
                raise TypeError("Ints or doubles are required")
        assert all(dataisint) or not any(dataisint), "Mix of ints and doubles found"
        dataisint = all(dataisint)

        entry = handle["entry"]
        if "data" in entry:
            data_group = entry["data"]
            dset = data_group["data"]
            dset.resize(dset.shape[0] + 1, axis=0)
        else:
            data_group = entry.create_group("data")
            data_group.attrs["NX_class"] = "NXdata"
            dset = data_group.create_dataset(
                "data",
                (1, *shape),
                maxshape=(None, *shape),
                dtype=int if dataisint else float,
            )

        if len(data) > 1:
            data = [flumpy.to_numpy(d) for d in data]
            data = np.stack(data)
        else:
            data = flumpy.to_numpy(data[0])
        dset[-1:] = data


def run(args):
    from dials.util.options import ArgumentParser, flatten_experiments

    # phil_scope = parse("""
    # """)
    from xfel.swissfel.jf16m_cxigeom2nexus import phil_scope

    parser = ArgumentParser(
        usage=None,
        sort_options=True,
        phil=phil_scope,
        read_experiments_from_images=True,
        epilog=None,
    )
    params, options = parser.parse_args(args=args)

    experiments = flatten_experiments(params.input.experiments)
    print(len(experiments))

    writer = NXmxWriter(params, experiments)
    handle = writer.get_h5py_handle()
    for imageset in experiments.imagesets():
        for i in range(len(imageset)):
            writer.append_frame(handle, data=imageset[i])


if __name__ == "__main__":
    run(sys.argv[1:])
