# LIBTBX_SET_DISPATCHER_NAME dxtbx.radial_average
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import annotations

import math
import os
import sys

import numpy as np

import cctbx.miller
import iotbx.phil
from cctbx.crystal import symmetry
from libtbx import easy_pickle
from libtbx.utils import Sorry, Usage
from scitbx.array_family import flex
from scitbx.matrix import col
from xfel import radial_average

import dxtbx
import dxtbx.util
from dxtbx.model.experiment_list import ExperimentListFactory

master_phil = iotbx.phil.parse(
    """
  file_path = None
    .type = str
    .multiple = True
    .help = Files to read
  n_bins = 0
    .type = int
    .help = Number of bins in the radial average. Auto determined if set to 0
  verbose = True
    .type = bool
    .help = Extra logging information
  output_bins = True
    .type = bool
    .help = Whether to print values for each bin
  output_file = None
    .type = str
    .help = Output file for logging results
  plot_x_max = None
    .type = int
    .help = Max value for x axis
  plot_y_max = None
    .type = int
    .help = Max value for xyaxis
  low_max_two_theta_limit = None
    .type = float
    .help = Low two theta cutoff
  normalize = False
    .type = bool
    .help = Whether to normalize the Y values to 1
  show_plots = True
    .type = bool
    .help = Whether to show the radial average plot
  mask = None
    .type = str
    .help = DIALS style pixel mask. Average will skip these pixels
  median_filter_size = None
    .type = int
    .help = If not none, applies a scipy.ndimage median_filter to the average
  x_axis = *two_theta q resolution
    .type = choice
    .help = Units for x axis
  image_number = None
    .type = int
    .help = When supplying a composite file, which image to show. Otherwise \
            shows all images in the file.
  panel = None
    .type = int
    .help = Only use data from the specified panel
  max_images = None
    .type = int
    .help = When supplying a composite file, only show up to max_images images
  reference_geometry = None
    .type = path
    .help = Apply this geometry before creating average
  unit_cell = None
    .type = unit_cell
    .help = Show positions of miller indices from this unit_cell and space \
            group.
  space_group = None
    .type = space_group
    .help = Show positions of miller indices from this unit_cell and space \
            group.
  d_min = 4
    .type = float
    .help = Resolution limit for miller index display
"""
)


def run(args=None, imageset=None):
    dxtbx.util.encode_output_as_utf8()
    args = sys.argv[1:] if args is None else args

    # Parse input
    try:
        len(args)
    except TypeError:
        params = args
    else:
        user_phil = []
        for arg in args:
            if "=" in arg:
                try:
                    user_phil.append(iotbx.phil.parse(arg))
                except RuntimeError as e:
                    raise Sorry(
                        "Unrecognized argument '%s' (error: %s)" % (arg, str(e))
                    )
            else:
                try:
                    user_phil.append(iotbx.phil.parse("""file_path=%s""" % arg))
                except ValueError:
                    raise Sorry("Unrecognized argument '%s'" % arg)
        params = master_phil.fetch(sources=user_phil).extract()
    if imageset is None:
        if (
            params.file_path is None
            or len(params.file_path) == 0
            or not all(os.path.isfile(f) for f in params.file_path)
        ):
            master_phil.show(attributes_level=2)
            raise Usage(
                "file_path must be defined (either file_path=XXX, or the path alone)."
            )
    assert params.n_bins is not None
    assert params.verbose is not None
    assert params.output_bins is not None

    if params.reference_geometry is not None:
        ref_expts = ExperimentListFactory.from_json_file(
            params.reference_geometry, check_format=None
        )
        assert (
            len(ref_expts.detectors()) == 1
        ), "Provide only one detector in the reference geometry file"
        detector = ref_expts.detectors()[0]

    # Allow writing to a file instead of stdout
    if params.output_file is None:
        logger = sys.stdout
    else:
        logger = open(params.output_file, "w")
        logger.write("%s " % params.output_file)

    if params.show_plots:
        from matplotlib import pyplot as plt

        colormap = plt.cm.gist_ncar
        plt.gca().set_prop_cycle(
            color=[colormap(i) for i in np.linspace(0, 0.9, len(params.file_path))]
        )

    if params.mask is not None and not isinstance(params.mask, tuple):
        params.mask = easy_pickle.load(params.mask)

    if imageset is None:
        iterable = params.file_path

        def load_func(x):
            try:
                obj = dxtbx.load(x).get_imageset([x])
            except (TypeError, ValueError):
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
            try:
                beam = iset.get_beam(image_number)
            except Exception:
                continue
            if params.reference_geometry is None:
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
            logger.write(f"Maximum {params.x_axis}: {max_x:f}, value: {max_result:f}\n")

            if params.show_plots:
                if params.plot_x_max is not None:
                    if params.x_axis == "resolution":
                        results = results.select(xvals >= params.plot_x_max)
                        xvals = xvals.select(xvals >= params.plot_x_max)
                    else:
                        results = results.select(xvals <= params.plot_x_max)
                        xvals = xvals.select(xvals <= params.plot_x_max)
                if params.x_axis == "resolution":
                    xvals = 1 / (xvals**2)
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
                    from matplotlib.ticker import FuncFormatter

                    plt.xlabel("Resolution ($\\AA$)")

                    def resolution(x, pos):
                        if x <= 0:
                            return "-"
                        return "%.1f" % (1 / math.sqrt(x))

                    formatter = FuncFormatter(resolution)
                    plt.gca().xaxis.set_major_formatter(formatter)
                plt.ylabel("Avg ADUs")
                if params.plot_y_max is not None:
                    plt.ylim(0, params.plot_y_max)

    if params.show_plots:
        if params.unit_cell or params.space_group:
            assert params.unit_cell and params.space_group
            sym = symmetry(
                unit_cell=params.unit_cell, space_group=params.space_group.group()
            )
            hkl_list = cctbx.miller.build_set(sym, False, d_min=params.d_min)
            tt = sym.unit_cell().two_theta(
                hkl_list.indices(), beam.get_wavelength(), deg=True
            )
            if params.x_axis == "q":
                q_vals = (
                    4 * math.pi * flex.sin(math.pi * tt / 360) / beam.get_wavelength()
                )
                vals = q_vals
            elif params.x_axis == "resolution":
                # nlmbda = 2dsin(theta)
                resolution = flex.double(len(tt), 0)
                nonzero = tt > 0
                resolution.set_selected(
                    nonzero,
                    beam.get_wavelength()
                    / (2 * flex.asin((math.pi / 180) * tt.select(nonzero) / 2)),
                )
                vals = resolution
                vals = 1 / (vals**2)
            elif params.x_axis == "two_theta":
                vals = tt

            y_max = 1 if params.normalize else 50
            for val in vals:
                plt.plot([val, val], [0, y_max], "r")

        # plt.legend([os.path.basename(os.path.splitext(f)[0]) for f in params.file_path], ncol=2)
        plt.show()

    return xvals, results


if __name__ == "__main__":
    run()
