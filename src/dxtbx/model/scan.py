from __future__ import annotations

import os

import pycbf

import libtbx.phil
from scitbx.array_family import flex

from dxtbx.model.scan_helpers import scan_helper_image_files

try:
    from ..dxtbx_model_ext import Scan
except ModuleNotFoundError:
    from dxtbx_model_ext import Scan  # type: ignore

scan_phil_scope = libtbx.phil.parse(
    """
  scan
    .expert_level = 1
    .short_caption = "Scan overrides"
  {

    image_range = None
      .type = ints(size=2)
      .help = "Override the image range"
      .short_caption = "Image range"

    extrapolate_scan = False
      .type = bool
      .help = "When overriding the image range, extrapolate exposure and epoch information from existing images"
      .short_caption = "Extrapolate scan"

    oscillation = None
      .type = floats(size=2)
      .help = "Override the image oscillation"
      .short_caption = "Oscillation"

    batch_offset = None
      .type = int(value_min=0)
      .help = "Override the batch offset"
      .short_caption = "Batch offset"
  }
"""
)


class ScanFactory:
    """A factory for scan instances, to help with constructing the classes
    in a set of common circumstances."""

    @staticmethod
    def from_phil(params, reference=None):
        """
        Generate a scan model from phil parameters

        """
        if reference is None:
            if params.scan.image_range is None and params.scan.oscillation is None:
                return None
            if params.scan.image_range is None:
                raise RuntimeError("No image range set")
            if params.scan.oscillation is None:
                raise RuntimeError("No oscillation set")
            scan = Scan(params.scan.image_range, params.scan.oscillation)
        else:
            scan = reference

            if params.scan.image_range is not None:
                most_recent_image_index = (
                    scan.get_image_range()[1] - scan.get_image_range()[0]
                )
                scan.set_oscillation(
                    scan.get_image_oscillation(params.scan.image_range[0])
                )
                scan.set_image_range(params.scan.image_range)
                if (
                    params.scan.extrapolate_scan
                    and (params.scan.image_range[1] - params.scan.image_range[0])
                    > most_recent_image_index
                ):
                    exposure_times = scan.get_exposure_times()
                    epochs = scan.get_epochs()
                    exposure_time = exposure_times[most_recent_image_index]
                    epoch_correction = epochs[most_recent_image_index]
                    for i in range(
                        most_recent_image_index + 1,
                        params.scan.image_range[1] - params.scan.image_range[0] + 1,
                    ):
                        exposure_times[i] = exposure_time
                        epoch_correction += exposure_time
                        epochs[i] = epoch_correction
                    scan.set_epochs(epochs)
                    scan.set_exposure_times(exposure_times)

            if params.scan.oscillation is not None:
                scan.set_oscillation(params.scan.oscillation)

        if params.scan.batch_offset is not None:
            scan.set_batch_offset(params.scan.batch_offset)

        return scan

    @staticmethod
    def from_dict(d, t=None):
        """Convert the dictionary to a scan model

        Params:
            d The dictionary of parameters
            t The template dictionary to use

        Returns:
            The scan model
        """

        def add_properties_table(scan_dict, num_images):

            """
            Handles legacy case before Scan had a properties table.
            Moves oscillation, epochs, and exposure times to a properties
            table and adds this to scan_dict.
            """

            properties = {}
            if scan_dict:
                if "oscillation" in scan_dict:
                    if num_images == 1:
                        properties["oscillation_width"] = [scan_dict["oscillation"][1]]
                        properties["oscillation"] = [scan_dict["oscillation"][0]]

                    else:
                        osc = scan_dict["oscillation"]
                        properties["oscillation"] = [
                            osc[0] + (osc[1] - osc[0]) * i for i in range(num_images)
                        ]
                    del scan_dict["oscillation"]
                if "exposure_time" in scan_dict:
                    properties["exposure_time"] = scan_dict["exposure_time"]
                    del scan_dict["exposure_time"]
                if "epochs" in scan_dict:
                    properties["epochs"] = scan_dict["epochs"]
                    del scan_dict["epochs"]
            scan_dict["properties"] = make_properties_table_consistent(
                properties, num_images
            )
            return scan_dict

        def make_properties_table_consistent(properties, num_images):

            """
            Handles legacy case before Scan had a properties table.
            Ensures oscillation, epochs, and exposure times have the same length.
            """

            if not properties:
                return properties

            property_length = len(next(iter(properties.values())))
            all_same_length = all(
                len(lst) == property_length for lst in properties.values()
            )
            if all_same_length and property_length == num_images:
                return properties

            if "oscillation" in properties:
                assert len(properties["oscillation"]) > 0

                if num_images == 1 and "oscillation_width" not in properties:
                    assert len(properties["oscillation"]) > 1
                    properties["oscillation_width"] = [properties["oscillation"][1]]
                    properties["oscillation"] = [properties["oscillation"][0]]
                elif num_images > 1:
                    osc_0 = properties["oscillation"][0]
                    if "oscillation_width" in properties:
                        osc_1 = properties["oscillation_width"][0]
                        del properties["oscillation_width"]
                    else:
                        assert len(properties["oscillation"]) > 1
                        osc_1 = (
                            properties["oscillation"][1] - properties["oscillation"][0]
                        )
                    properties["oscillation"] = [
                        osc_0 + osc_1 * i for i in range(num_images)
                    ]

            if "exposure_time" in properties:
                assert len(properties["exposure_time"]) > 0

                # Assume same exposure time for each image
                properties["exposure_time"] = [
                    properties["exposure_time"][0] for i in range(num_images)
                ]

            if "epochs" in properties:
                assert len(properties["epochs"]) > 0
                # If 1 epoch, assume increasing by epochs[0]
                # Else assume increasing as epochs[1] - epochs[0]
                if len(properties["epochs"]) == 1:
                    properties["epochs"] = [
                        properties["epochs"][0] + properties["epochs"][0] * i
                        for i in range(num_images)
                    ]
                else:
                    diff = properties["epochs"][1] - properties["epochs"][0]
                    properties["epochs"] = [
                        properties["epochs"][0] + i * diff for i in range(num_images)
                    ]
            return properties

        if d is None and t is None:
            return None
        joint = t.copy() if t else {}

        # Accounting for legacy cases where t or d does not
        # contain properties dict
        num_images = None
        if "image_range" in d:
            num_images = 1 + d["image_range"][1] - d["image_range"][0]
        elif "image_range" in joint:
            num_images = 1 + joint["image_range"][1] - joint["image_range"][0]

        if "properties" in joint and "properties" in d:
            properties = t["properties"].copy()
            properties.update(d["properties"])
            joint.update(d)
            joint["properties"] = properties
        elif "properties" in d:
            joint = add_properties_table(joint, num_images)
            d_copy = d.copy()
            joint["properties"].update(d_copy["properties"])
            joint["properties"] = make_properties_table_consistent(
                joint["properties"], num_images
            )
            del d_copy["properties"]
            joint.update(d_copy)
        elif "properties" in joint:
            d = add_properties_table(d, num_images)
            d_copy = d.copy()
            joint["properties"].update(d_copy["properties"])
            joint["properties"] = make_properties_table_consistent(
                joint["properties"], num_images
            )
            del d_copy["properties"]
            joint.update(d_copy)
        else:
            joint.update(d)

        if "properties" not in d and not isinstance(joint["exposure_time"], list):
            joint["exposure_time"] = [joint["exposure_time"]]
        joint.setdefault("batch_offset", 0)  # backwards compatibility 20180205
        joint.setdefault("valid_image_ranges", {})  # backwards compatibility 20181113

        return Scan.from_dict(joint)

    @staticmethod
    def make_scan(
        image_range, exposure_times, oscillation, epochs, batch_offset=0, deg=True
    ):
        if not isinstance(exposure_times, list):
            num_images = image_range[1] - image_range[0] + 1
            exposure_times = [exposure_times for i in range(num_images)]
        else:
            num_images = image_range[1] - image_range[0] + 1
            num_exp = len(exposure_times)
            if num_exp != num_images:
                if num_exp == 0:
                    exposure_times = [0 for i in range(num_images)]
                else:
                    exposure_times = exposure_times.extend(
                        [exposure_times[-1] for i in range(num_images - num_exp)]
                    )

        epoch_list = [epochs[j] for j in sorted(epochs)]

        return Scan(
            tuple(map(int, image_range)),
            tuple(map(float, oscillation)),
            flex.double(list(map(float, exposure_times))),
            flex.double(list(map(float, epoch_list))),
            batch_offset,
            deg,
        )

    @staticmethod
    def make_scan_from_properties(image_range, properties, batch_offset=0, deg=True):

        return Scan(tuple(map(int, image_range)), properties, batch_offset, deg)

    @staticmethod
    def single_file(filename, exposure_times, osc_start, osc_width, epoch):
        """Construct an scan instance for a single image."""
        index = scan_helper_image_files.image_to_index(os.path.split(filename)[-1])
        if epoch is None:
            epoch = 0.0

        # if the oscillation width is negative at this stage it is almost
        # certainly an artefact of the omega end being 0 when the omega start
        # angle was ~ 360 so it would be ~ -360 - see dxtbx#378
        if osc_width < -180:
            osc_width += 360

        return ScanFactory.make_scan(
            (index, index), exposure_times, (osc_start, osc_width), {index: epoch}
        )

    @staticmethod
    def imgCIF(cif_file):
        """Initialize a scan model from an imgCIF file."""

        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

        return ScanFactory.imgCIF_H(cif_file, cbf_handle)

    @staticmethod
    def imgCIF_H(cif_file, cbf_handle):
        """Initialize a scan model from an imgCIF file handle, where it is
        assumed that the file has already been read."""

        exposure = cbf_handle.get_integration_time()
        timestamp = cbf_handle.get_timestamp()[0]

        gonio = cbf_handle.construct_goniometer()
        try:
            angles = tuple(gonio.get_rotation_range())
        except Exception as e:
            if str(e).strip() == "CBFlib Error(s): CBF_NOTFOUND":
                # probably a still shot -> no scan object
                return None
            raise

        # xia2-56 handle gracefully reverse turning goniometers - this assumes the
        # rotation axis is correctly inverted in the goniometer factory
        if angles[1] < 0:
            angles = -angles[0], -angles[1]

        index = scan_helper_image_files.image_to_index(cif_file)

        gonio.__swig_destroy__(gonio)

        return ScanFactory.make_scan(
            (index, index), exposure, angles, {index: timestamp}
        )

    @staticmethod
    def add(scans):
        """Sum a list of scans wrapping the slightly clumsy idiomatic method:
        sum(scans[1:], scans[0])."""
        return sum(scans[1:], scans[0])

    @staticmethod
    def search(filename):
        """Get a list of files which appear to match the template and
        directory implied by the input filename. This could well be used
        to get a list of image headers to read and hence construct scans
        from."""

        template, directory = scan_helper_image_files.image_to_template_directory(
            filename
        )

        indices = scan_helper_image_files.template_directory_to_indices(
            template, directory
        )

        return [
            scan_helper_image_files.template_directory_index_to_image(
                template, directory, index
            )
            for index in indices
        ]
