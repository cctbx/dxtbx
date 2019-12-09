from __future__ import absolute_import, division, print_function

import os
from builtins import object, range

import libtbx.phil
from scitbx.array_family import flex

import pycbf
from dxtbx.model.scan_helpers import scan_helper_image_files, scan_helper_image_formats
from dxtbx_model_ext import Scan

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


class ScanFactory(object):
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

        # Return the model
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
        if d is None and t is None:
            return None
        joint = t.copy() if t else {}
        joint.update(d)

        if not isinstance(joint["exposure_time"], list):
            joint["exposure_time"] = [joint["exposure_time"]]
        joint.setdefault("batch_offset", 0)  # backwards compatibility 20180205
        joint.setdefault("valid_image_ranges", {})  # backwards compatibility 20181113

        # Create the model from the joint dictionary
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
    def single(filename, format, exposure_times, osc_start, osc_width, epoch):
        """Construct an scan instance for a single image."""

        index = scan_helper_image_files.image_to_index(os.path.split(filename)[-1])
        if epoch is None:
            epoch = 0.0
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
                # probaby a still shot -> no scan object
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
        """Sum a list of scans wrapping the sligtly clumsy idiomatic method:
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

    @staticmethod
    def format(name):
        """Return the correct format token for a given name, for example:

        cbf, CBF
        smv, SMV
        tiff, tif, TIFF
        raxis, RAXIS
        mar, MAR

        to the appropriate static token which will be used as a handle
        everywhere else in this."""

        if name.upper() == "CBF":
            return scan_helper_image_formats.FORMAT_CBF
        elif name.upper() == "SMV":
            return scan_helper_image_formats.FORMAT_SMV
        elif name.upper() == "TIF" or name.upper() == "TIFF":
            return scan_helper_image_formats.FORMAT_TIFF
        elif name.upper() == "RAXIS":
            return scan_helper_image_formats.FORMAT_RAXIS
        elif name.upper() == "MAR":
            return scan_helper_image_formats.FORMAT_MAR

        raise RuntimeError("name %s not known" % name)
