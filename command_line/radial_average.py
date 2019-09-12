# LIBTBX_SET_DISPATCHER_NAME dxtbx.radial_average
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import absolute_import, division, print_function

import math
import os
import sys
from builtins import range

import numpy as np

import libtbx.phil
from libtbx import easy_pickle
from libtbx.utils import Sorry, Usage
from scitbx.array_family import flex
from scitbx.matrix import col
from xfel import radial_average

import dxtbx.datablock
from dxtbx.model.experiment_list import ExperimentListFactory

master_phil = libtbx.phil.parse(
    """
  file_path = None
    .type = str
    .multiple = True
  beam_x = None
    .type = float
  beam_y = None
    .type = float
  n_bins = 0
    .type = int
  verbose = True
    .type = bool
  output_bins = True
    .type = bool
  output_file = None
    .type = str
  plot_x_max = None
    .type = int
  plot_y_max = None
    .type = int
  low_max_two_theta_limit = None
    .type = float
  normalize = False
    .type = bool
  show_plots = True
    .type = bool
  mask = None
    .type = str
  median_filter_size = None
    .type = int
  x_axis = *two_theta q resolution
    .type = choice
  image_number = None
    .type = int
  panel = None
    .type = int
  max_images = None
    .type = int
"""
)


def run(args, imageset=None):
    # Parse input
    try:
        len(args)
    except Exception:
        params = args
    else:
        user_phil = []
        for arg in args:
            if "=" in arg:
                try:
                    user_phil.append(libtbx.phil.parse(arg))
                except RuntimeError as e:
                    raise Sorry(
                        "Unrecognized argument '%s' (error: %s)" % (arg, str(e))
                    )
            else:
                try:
                    user_phil.append(libtbx.phil.parse("""file_path=%s""" % arg))
                except ValueError:
                    raise Sorry("Unrecognized argument '%s'" % arg)
        params = master_phil.fetch(sources=user_phil).extract()
    if imageset is None:
        if (
            params.file_path is None
            or len(params.file_path) == 0
            or not all(os.path.isfile(f) for f in params.file_path)
        ):
            master_phil.show()
            raise Usage(
                "file_path must be defined (either file_path=XXX, or the path alone)."
            )
    assert params.n_bins is not None
    assert params.verbose is not None
    assert params.output_bins is not None

    # Allow writing to a file instead of stdout
    if params.output_file is None:
        logger = sys.stdout
    else:
        logger = open(params.output_file, "w")
        logger.write("%s " % params.output_file)

    if params.show_plots:
        from matplotlib import pyplot as plt

        colormap = plt.cm.gist_ncar
        plt.gca().set_color_cycle(
            [colormap(i) for i in np.linspace(0, 0.9, len(params.file_path))]
        )

    if params.mask is not None:
        params.mask = easy_pickle.load(params.mask)

    if imageset is None:
        iterable = params.file_path

        def load_func(x):
            try:
                obj = dxtbx.datablock.DataBlockFactory.from_filenames([x])[
                    0
                ].extract_imagesets()[0]
            except IndexError:
                try:
                    obj = dxtbx.datablock.DataBlockFactory.from_json_file(x)[
                        0
                    ].extract_imagesets()[0]
                except dxtbx.datablock.InvalidDataBlockError:
                    obj = ExperimentListFactory.from_json_file(x)[0].imageset
            return obj

    else:
        iterable = [imageset]

        def load_func(x):
            return x

    # Iterate over each file provided
    for item in iterable:
        iset = load_func(item)
        n_images = len(iset)
        if params.image_number is None:
            if params.max_images is None:
                subiterable = range(n_images)
            else:
                subiterable = range(0, min(params.max_images, n_images))
        else:
            subiterable = [params.image_number]
        for image_number in subiterable:
            beam = iset.get_beam(image_number)
            detector = iset.get_detector(image_number)
            s0 = col(beam.get_s0())

            # Search the detector for the panel farthest from the beam. The
            # number of bins in the radial average will be equal to the
            # farthest point from the beam on the detector, in pixels, unless
            # overridden at the command line
            panel_res = [p.get_max_resolution_at_corners(s0) for p in detector]
            farthest_panel = detector[panel_res.index(min(panel_res))]
            size2, size1 = farthest_panel.get_image_size()
            corners = [(0, 0), (size1 - 1, 0), (0, size2 - 1), (size1 - 1, size2 - 1)]
            corners_lab = [col(farthest_panel.get_pixel_lab_coord(c)) for c in corners]
            corner_two_thetas = [
                farthest_panel.get_two_theta_at_pixel(s0, c) for c in corners
            ]
            extent_two_theta = max(corner_two_thetas)
            max_corner = corners_lab[corner_two_thetas.index(extent_two_theta)]
            extent = int(
                math.ceil(
                    max_corner.length()
                    * math.sin(extent_two_theta)
                    / max(farthest_panel.get_pixel_size())
                )
            )
            extent_two_theta *= 180 / math.pi

            if params.n_bins < extent:
                params.n_bins = extent

            # These arrays will store the radial average info
            sums = flex.double(params.n_bins) * 0
            sums_sq = flex.double(params.n_bins) * 0
            counts = flex.int(params.n_bins) * 0

            all_data = iset[image_number]

            if not isinstance(all_data, tuple):
                all_data = (all_data,)

            for tile, (panel, data) in enumerate(zip(detector, all_data)):
                if params.panel is not None and tile != params.panel:
                    continue

                if params.mask is None:
                    mask = flex.bool(flex.grid(data.focus()), True)
                else:
                    mask = params.mask[tile]

                if hasattr(data, "as_double"):
                    data = data.as_double()

                logger.flush()
                if params.verbose:
                    logger.write(
                        "Average intensity tile %d: %9.3f\n" % (tile, flex.mean(data))
                    )
                    logger.write("N bins: %d\n" % params.n_bins)
                    logger.flush()

                x1, y1, x2, y2 = (
                    0,
                    0,
                    panel.get_image_size()[1],
                    panel.get_image_size()[0],
                )
                bc = panel.get_beam_centre_px(beam.get_s0())
                bc = int(round(bc[1])), int(round(bc[0]))

                # compute the average
                radial_average(
                    data,
                    mask,
                    bc,
                    sums,
                    sums_sq,
                    counts,
                    panel.get_pixel_size()[0],
                    panel.get_distance(),
                    (x1, y1),
                    (x2, y2),
                )

            # average the results, avoiding division by zero
            results = sums.set_selected(counts <= 0, 0)
            results /= counts.set_selected(counts <= 0, 1).as_double()

            if params.median_filter_size is not None:
                logger.write(
                    "WARNING, the median filter is not fully propagated to the variances\n"
                )
                from scipy.ndimage.filters import median_filter

                results = flex.double(
                    median_filter(
                        results.as_numpy_array(), size=params.median_filter_size
                    )
                )

            # calculate standard devations
            stddev_sel = ((sums_sq - sums * results) >= 0) & (counts > 0)
            std_devs = flex.double(len(sums), 0)
            std_devs.set_selected(
                stddev_sel,
                (
                    sums_sq.select(stddev_sel)
                    - sums.select(stddev_sel) * results.select(stddev_sel)
                )
                / counts.select(stddev_sel).as_double(),
            )
            std_devs = flex.sqrt(std_devs)

            twotheta = (
                flex.double(range(len(results))) * extent_two_theta / params.n_bins
            )
            q_vals = (
                4 * math.pi * flex.sin(math.pi * twotheta / 360) / beam.get_wavelength()
            )
            # nlmbda = 2dsin(theta)
            resolution = flex.double(len(twotheta), 0)
            nonzero = twotheta > 0
            resolution.set_selected(
                nonzero,
                beam.get_wavelength()
                / (2 * flex.asin((math.pi / 180) * twotheta.select(nonzero) / 2)),
            )

            if params.low_max_two_theta_limit is None:
                subset = results
            else:
                subset = results.select(twotheta >= params.low_max_two_theta_limit)

            max_result = flex.max(subset)

            if params.x_axis == "two_theta":
                xvals = twotheta
                max_x = twotheta[flex.first_index(results, max_result)]
            elif params.x_axis == "q":
                xvals = q_vals
                max_x = q_vals[flex.first_index(results, max_result)]
            elif params.x_axis == "resolution":
                xvals = resolution
                max_x = resolution[flex.first_index(results, max_result)]

            for i, r in enumerate(results):
                val = xvals[i]
                if params.output_bins and "%.3f" % r != "nan":
                    # logger.write("%9.3f %9.3f\n"%     (val,r))        #.xy  format for Rex.cell.
                    logger.write(
                        "%9.3f %9.3f %9.3f\n" % (val, r, std_devs[i])
                    )  # .xye format for GSASII
                # logger.write("%.3f %.3f %.3f\n"%(val,r,ds[i]))  # include calculated d spacings
            logger.write(
                "Maximum %s: %f, value: %f\n" % (params.x_axis, max_x, max_result)
            )

            if params.show_plots:
                if params.plot_x_max is not None:
                    results = results.select(xvals <= params.plot_x_max)
                    xvals = xvals.select(xvals <= params.plot_x_max)
                if params.normalize:
                    plt.plot(
                        xvals.as_numpy_array(),
                        (results / flex.max(results)).as_numpy_array(),
                        "-",
                    )
                else:
                    plt.plot(xvals.as_numpy_array(), results.as_numpy_array(), "-")
                if params.x_axis == "two_theta":
                    plt.xlabel("2 theta")
                elif params.x_axis == "q":
                    plt.xlabel("q")
                elif params.x_axis == "resolution":
                    plt.xlabel("Resolution ($\\AA$)")
                    plt.gca().set_xscale("log")
                    plt.gca().invert_xaxis()
                    plt.xlim(0, 50)
                plt.ylabel("Avg ADUs")
                if params.plot_y_max is not None:
                    plt.ylim(0, params.plot_y_max)

    if params.show_plots:
        # plt.legend([os.path.basename(os.path.splitext(f)[0]) for f in params.file_path], ncol=2)
        plt.show()

    return xvals, results


if __name__ == "__main__":
    run(sys.argv[1:])
