"""
Goniometer and scan models are written when the imageset is a sequence (or when
they are passed explicitly), so rotation data can be round-tripped as well as stills

Example to write the first 10 images from an h5 file:
writer = FullCBFWriter("data.h5")
for i in range(10):
  writer.write_cbf("example_%d.cbf"%i, index=i)
"""

from __future__ import annotations

import datetime
import os
import sys

import pycbf

from scitbx.array_family import flex
from serialtbx.detector import basis

import dxtbx.format.Registry
from dxtbx.format.FormatCBFMultiTile import cbf_wrapper


def add_frame_specific_cbf_tables(
    cbf,
    wavelength,
    timestamp,
    trusted_ranges,
    diffrn_id="DS1",
    is_xfel=True,
    gain=1.0,
    flux=None,
    goniometer_axis_names=None,
    exposure_time=0.0,
):
    """Adds tables to cbf handle that won't already exsist if the cbf file is just a header
    @ param wavelength Wavelength in angstroms
    @ param timestamp String formatted timestamp for the image
    @ param trusted_ranges Array of trusted range tuples (min, max), one for each element
    @ param goniometer_axis_names Names of goniometer axes (rotation data), if any
    @ param exposure_time Exposure time of the frame in seconds
    """

    """Data items in the DIFFRN_RADIATION category describe
   the radiation used for measuring diffraction intensities,
   its collimation and monochromatization before the sample.

   Post-sample treatment of the beam is described by data
   items in the DIFFRN_DETECTOR category."""
    if flux:
        cbf.add_category(
            "diffrn_radiation", ["diffrn_id", "wavelength_id", "probe", "beam_flux"]
        )
        cbf.add_row([diffrn_id, "WAVELENGTH1", "x-ray", "%f" % flux])
    else:
        cbf.add_category("diffrn_radiation", ["diffrn_id", "wavelength_id", "probe"])
        cbf.add_row([diffrn_id, "WAVELENGTH1", "x-ray"])

    """ Data items in the DIFFRN_RADIATION_WAVELENGTH category describe
   the wavelength of the radiation used in measuring the diffraction
   intensities. Items may be looped to identify and assign weights
   to distinct wavelength components from a polychromatic beam."""
    cbf.add_category("diffrn_radiation_wavelength", ["id", "wavelength", "wt"])
    cbf.add_row(["WAVELENGTH1", str(wavelength), "1.0"])

    """Data items in the DIFFRN_MEASUREMENT category record details
   about the device used to orient and/or position the crystal
   during data measurement and the manner in which the
   diffraction data were measured."""
    cbf.add_category(
        "diffrn_measurement", ["diffrn_id", "id", "number_of_axes", "method", "details"]
    )
    goniometer_axis_names = goniometer_axis_names or []
    if goniometer_axis_names:
        measurement_id = "GONIOMETER"
        cbf.add_row(
            [
                diffrn_id,
                measurement_id,
                str(len(goniometer_axis_names)),
                "rotation",
                ".",
            ]
        )
        """Data items in the DIFFRN_MEASUREMENT_AXIS category associate
       axes with goniometers."""
        cbf.add_category("diffrn_measurement_axis", ["measurement_id", "axis_id"])
        for name in goniometer_axis_names:
            cbf.add_row([measurement_id, name])
    else:
        cbf.add_row(
            [
                diffrn_id,
                "INJECTION" if is_xfel else "unknown",
                "0",
                (
                    "electrospray"
                    if is_xfel
                    else (
                        "unknowncrystals injected by electrospray"
                        if is_xfel
                        else "unknown"
                    )
                ),
            ]
        )

    """ Data items in the DIFFRN_SCAN category describe the parameters of one
     or more scans, relating axis positions to frames."""
    cbf.add_category("diffrn_scan", ["id", "frame_id_start", "frame_id_end", "frames"])
    cbf.add_row(["SCAN1", "FRAME1", "FRAME1", "1"])

    """Data items in the DIFFRN_SCAN_FRAME category describe
   the relationships of particular frames to scans."""
    cbf.add_category(
        "diffrn_scan_frame",
        ["frame_id", "frame_number", "integration_time", "scan_id", "date"],
    )
    cbf.add_row(["FRAME1", "1", "%g" % exposure_time, "SCAN1", timestamp])

    """ Data items in the ARRAY_INTENSITIES category record the
   information required to recover the intensity data from
   the set of data values stored in the ARRAY_DATA category."""
    # More detail here: http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Carray_intensities.html
    array_names = []
    cbf.find_category(b"diffrn_data_frame")
    while True:
        try:
            cbf.find_column(b"array_id")
            array_names.append(cbf.get_value().decode())
            cbf.next_row()
        except Exception as e:
            assert "CBF_NOTFOUND" in str(e)
            break

    if not isinstance(gain, list):
        gain = [gain] * len(array_names)

    cbf.add_category(
        "array_intensities",
        [
            "array_id",
            "binary_id",
            "linearity",
            "gain",
            "gain_esd",
            "overload",
            "underload",
            "undefined_value",
        ],
    )
    for i, array_name in enumerate(array_names):
        overload = trusted_ranges[i][1] + 1
        underload = trusted_ranges[i][0]
        undefined = underload - 1
        cbf.add_row(
            [
                array_name,
                str(i + 1),
                "linear",
                "%f" % gain[i],
                "0.0",
                str(overload),
                str(underload),
                str(undefined),
            ]
        )


class FullCBFWriter:
    """Class for writing full CBF files from any dxtbx-supported format class"""

    def __init__(self, filename=None, imageset=None):
        """Provide a file name or imageset as input"""
        assert [filename, imageset].count(None) == 1, (
            "Supply either filename or imageset"
        )

        if filename is not None:
            format_class = dxtbx.format.Registry.get_format_class_for_file(filename)
            imageset = format_class.get_imageset([filename])

        self.imageset = imageset

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
            detector = self.imageset.get_detector()
        else:
            detector = self.imageset.get_detector(index)

        recursive_setup_dict(detector.hierarchy(), (0,))
        return metro

    def get_cbf_handle(
        self,
        index=None,
        header_only=False,
        detector_only=False,
        goniometer=None,
        scan=None,
    ):
        """Build a cbf handle in memory. For rotation data, the goniometer and
        scan are taken from the imageset if it is a sequence, or may be passed
        explicitly; the frame written is the one at position index of the scan
        (or its first frame if index is None)."""
        # set up the metrology dictionary to include axis names, pixel sizes, and so forth
        if index is None:
            detector = self.imageset.get_detector()
            beam = self.imageset.get_beam()
        else:
            detector = self.imageset.get_detector(index)
            beam = self.imageset.get_beam(index)

        if goniometer is None:
            try:
                goniometer = self.imageset.get_goniometer()
            except Exception:
                goniometer = None
        if scan is None:
            try:
                scan = self.imageset.get_scan()
            except Exception:
                scan = None
        goniometer_axis_names = []
        gonio_axis_rows = []  # rows for the axis table
        gonio_scan_axis_rows = []  # rows for diffrn_scan_axis
        gonio_frame_axis_rows = []  # rows for diffrn_scan_frame_axis
        exposure_time = 0.0
        timestamp = "unknown"
        if goniometer is not None and scan is not None:
            frame = 0 if index is None else index
            osc_start, osc_width = scan.get_oscillation()
            scan_angle = osc_start + frame * osc_width
            exposure_times = scan.get_exposure_times()
            if len(exposure_times) > frame:
                exposure_time = exposure_times[frame]
            # the imgCIF scan reader needs a parseable date for each frame
            epochs = scan.get_epochs()
            epoch = epochs[frame] if len(epochs) > frame and epochs[frame] > 0 else 0.0
            timestamp = datetime.datetime.fromtimestamp(
                epoch, datetime.timezone.utc
            ).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
            if hasattr(goniometer, "get_axes"):
                # multi-axis goniometer
                axes = list(goniometer.get_axes())
                angles = list(goniometer.get_angles())
                names = list(goniometer.get_names())
                scan_axis = goniometer.get_scan_axis()
            else:
                # single axis goniometer
                axes = [goniometer.get_rotation_axis()]
                angles = [0.0]
                names = ["GONIOMETER_PHI"]
                scan_axis = 0
            if not names or any(not n for n in names):
                names = ["GONIOMETER_AXIS%d" % i for i in range(len(axes))]
            # dxtbx orders axes from the crystal outwards to the goniometer base;
            # imgCIF depends_on chains from each axis to the one it sits on
            for i, (name, axis) in enumerate(zip(names, axes)):
                depends_on = names[i + 1] if i + 1 < len(names) else "."
                gonio_axis_rows.append(
                    [name, "rotation", "goniometer", depends_on]
                    + ["%.10g" % v for v in axis]
                    + ["0", "0", "0", "."]
                )
                if i == scan_axis:
                    start, rng, inc, angle = (
                        scan_angle,
                        osc_width,
                        osc_width,
                        scan_angle,
                    )
                else:
                    start, rng, inc, angle = angles[i], 0.0, 0.0, angles[i]
                gonio_scan_axis_rows.append(
                    [name, "SCAN1", "%g" % start, "%g" % rng, "%g" % inc, "0", "0", "0"]
                )
                gonio_frame_axis_rows.append([name, "FRAME1", "%g" % angle, "0"])
            goniometer_axis_names = names

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
        dname = None
        detector_axes_names = []  # save these for later
        panelkeys = []
        panelnames = []
        panelindices = []

        def recursive_setup_basis_dict(key, parent_name="", panel_id=0):
            # Set up CBF axis names, including equipment components and depends_on chains
            basis = metro[key]
            node = panel_group_from_key(key)
            nodename = level_string(key)
            if basis.name:
                dname = basis.name
            else:
                dname = "AXIS_" + nodename

            if node.is_panel():
                panelname = "PANEL_%d" % panel_id
                panelkeys.append(key)
                panelnames.append(panelname)
                dxtbx_panel_index = list(detector).index(node)
                panelindices.append(dxtbx_panel_index)
                panel_id += 1

            if len(key) == 1:
                assert key == (0,)  # only one detector allowed for now

                for a in ["_X", "_Y", "_Z", "_R"]:
                    detector_axes_names.append(dname + a)
                basis.depends_on = dname + "_X"
                dname += "_R"
            else:
                detector_axes_names.append(dname)
                basis.depends_on = parent_name

            basis.equipment_component = "detector_level_%d" % (len(key) - 1)
            basis.axis_name = detector_axes_names[-1]
            if not node.is_panel():
                for c, child in enumerate(node):
                    panel_id = recursive_setup_basis_dict(
                        tuple(list(key) + [c]), dname, panel_id
                    )
            return panel_id

        recursive_setup_basis_dict((0,))
        # The order `recursive_setup_basis_dict` visits panels does not necessarily match
        # the order returned by dxtbx's `get_detector` and `get_raw_data`.
        # See https://github.com/cctbx/dxtbx/issues/745.
        self.panelindices = panelindices  # needed in add_data_to_cbf
        sorted_panels = [detector[i] for i in panelindices]

        if index is None:
            cbf_root = self.imageset.paths()[0]
        else:
            cbf_root = self.imageset.paths()[index]
        cbf_root = os.path.splitext(os.path.basename(cbf_root))[0]

        # the data block is the root cbf node
        cbf = cbf_wrapper()
        cbf.new_datablock(cbf_root.encode())

        # Each category listed here is preceded by the imageCIF description taken from here:
        # http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/index.html

        """Data items in the DIFFRN category record details about the
     diffraction data and their measurement."""
        diffrn_id = "dxtbx_%s" % self.imageset.reader().format_class.__name__
        cbf.add_category("diffrn", ["id"])
        cbf.add_row([diffrn_id])

        """Data items in the DIFFRN_SOURCE category record details of
    the source of radiation used in the diffraction experiment."""
        cbf.add_category("diffrn_source", ["diffrn_id", "source", "type"])
        cbf.add_row([diffrn_id, "unknown", "unknown"])

        """Data items in the DIFFRN_DETECTOR category describe the
     detector used to measure the scattered radiation, including
     any analyser and post-sample collimation."""
        cbf.add_category(
            "diffrn_detector", ["diffrn_id", "id", "type", "details", "number_of_axes"]
        )
        detector_id = "dxtbx_detector"
        cbf.add_row(
            [
                diffrn_id,
                detector_id,
                "General dxtbx detector",
                ".",
                str(len(detector_axes_names)),
            ]
        )

        """Data items in the DIFFRN_DETECTOR_AXIS category associate
       axes with detectors."""
        # Note, does not include the fast and the slow axes
        cbf.add_category("diffrn_detector_axis", ["detector_id", "axis_id"])
        for name in detector_axes_names:
            cbf.add_row([detector_id, name])

        """Data items in the DIFFRN_DETECTOR_ELEMENT category record
     the details about spatial layout and other characteristics
     of each element of a detector which may have multiple elements."""
        cbf.add_category("diffrn_detector_element", ["id", "detector_id"])

        for panelname in panelnames:
            cbf.add_row(["ELE_" + panelname, detector_id])

        """Data items in the DIFFRN_DATA_FRAME category record
     the details about each frame of data."""
        cbf.add_category(
            "diffrn_data_frame", ["id", "detector_element_id", "array_id", "binary_id"]
        )

        for i, panelname in enumerate(panelnames):
            cbf.add_row(
                ["FRAME1", "ELE_" + panelname, "ARRAY_" + panelname, "%d" % (i + 1)]
            )

        if not detector_only:
            # add_frame_specific_cbf_tables creates _array_intensities.undefined_value
            # and _array_intensities.overload. The trusted_range of the detector is
            # defined as [min-trusted-value, max-trusted-value]. The CBF definition
            # of 'overload' is in fact saturation - i.e. the max-trusted-value, while
            # the undefined_value is below the min-trusted-value.
            trusted_ranges = [panel.get_trusted_range() for panel in sorted_panels]
            try:
                add_frame_specific_cbf_tables(
                    cbf,
                    beam.get_wavelength(),
                    timestamp,
                    trusted_ranges,
                    diffrn_id,
                    False,
                    gain=[panel.get_gain() for panel in sorted_panels],
                    flux=beam.get_flux(),
                    goniometer_axis_names=goniometer_axis_names,
                    exposure_time=exposure_time,
                )
            except TypeError:
                # Needed until next cctbx release
                add_frame_specific_cbf_tables(
                    cbf,
                    beam.get_wavelength(),
                    "unknown",
                    trusted_ranges,
                    diffrn_id,
                    False,
                    gain=[panel.get_gain() for panel in sorted_panels],
                )

        """Data items in the AXIS category record the information required
       to describe the various goniometer, detector, source and other
       axes needed to specify a data collection.  The location of each
       axis is specified by two vectors: the axis itself, given as a unit
       vector, and an offset to the base of the unit vector.  These vectors
       are referenced to a right-handed laboratory coordinate system with
       its origin in the sample or specimen"""
        # More detail here: http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Caxis.html
        # Note we also use two new columns not in the latest imageCIF dictionary: rotation and rotation_axis.
        # We use them to specify an translation and a rotation in a single axis to describe a change in setting
        # when moving from one frame (say a quadrant) to another (say a sensor)
        cbf.add_category(
            "axis",
            [
                "id",
                "type",
                "equipment",
                "depends_on",
                "vector[1]",
                "vector[2]",
                "vector[3]",
                "offset[1]",
                "offset[2]",
                "offset[3]",
                "equipment_component",
            ],
        )

        # Keep a list of rows to add to the scan frame axis table
        axis_settings = []
        # keep a list of rows to add to the scan axis table
        axis_names = []

        root_key = (0,)
        dname = metro[root_key].axis_name[:-2]
        eqc = metro[root_key].equipment_component

        # Create a series of axis describing frame shifts from each level of the detector to the next
        cbf.add_row(
            "AXIS_SOURCE  general     source   .        0  0  1 . . . .".split()
        )
        axis_names.append("AXIS_SOURCE")
        cbf.add_row(
            "AXIS_GRAVITY general     gravity  .        0 -1  0 . . . .".split()
        )
        axis_names.append("AXIS_GRAVITY")
        for row in gonio_axis_rows:
            cbf.add_row(row)
        cbf.add_row(
            (
                "%s_Z         translation detector .        0  0  1 . . . %s"
                % (dname, eqc)
            ).split()
        )
        axis_names.append("%s_Z" % dname)
        cbf.add_row(
            (
                "%s_Y         translation detector %s_Z     0  1  0 . . . %s"
                % (dname, dname, eqc)
            ).split()
        )
        axis_names.append("%s_Y" % dname)
        cbf.add_row(
            (
                "%s_X         translation detector %s_Y     1  0  0 . . . %s"
                % (dname, dname, eqc)
            ).split()
        )
        axis_names.append("%s_X" % dname)

        axis_settings.append(["AXIS_SOURCE", "FRAME1", "0", "0"])
        axis_settings.append(["AXIS_GRAVITY", "FRAME1", "0", "0"])
        axis_settings.append([dname + "_X", "FRAME1", "0", "0"])
        axis_settings.append([dname + "_Y", "FRAME1", "0", "0"])
        axis_settings.append([dname + "_Z", "FRAME1", "0", "0"])

        for key in sorted(metro):
            basis = metro[key]

            cbf.add_frame_shift(basis, axis_settings)
            axis_names.append(basis.axis_name)

            node = panel_group_from_key(key)

            if node.is_panel():
                axis_settings[-1][-2] = (
                    "0"  # Drop the setting change for leaves as it's encoded below
                )

                aname = level_string(key)
                fast = [str(v) for v in node.get_local_fast_axis()]
                slow = [str(v) for v in node.get_local_slow_axis()]

                cbf.add_row(
                    [
                        "AXIS_" + aname + "_S",
                        "translation",
                        "detector",
                        basis.axis_name,
                        slow[0],
                        slow[1],
                        slow[2],
                        "0",
                        "0",
                        "0",
                        basis.equipment_component,
                    ]
                )
                cbf.add_row(
                    [
                        "AXIS_" + aname + "_F",
                        "translation",
                        "detector",
                        "AXIS_" + aname + "_S",
                        fast[0],
                        fast[1],
                        fast[2],
                        "0",
                        "0",
                        "0",
                        basis.equipment_component,
                    ]
                )
                axis_names.append("AXIS_" + aname + "_F")
                axis_names.append("AXIS_" + aname + "_S")
                axis_settings.append(["AXIS_" + aname + "_F", "FRAME1", "0", "0"])
                axis_settings.append(["AXIS_" + aname + "_S", "FRAME1", "0", "0"])

        """Data items in the DIFFRN_SCAN_AXIS category describe the settings of
       axes for particular scans.  Unspecified axes are assumed to be at
       their zero points."""
        # leave all the settings zero. Levels with settings are set below.
        cbf.add_category(
            "diffrn_scan_axis",
            [
                "axis_id",
                "scan_id",
                "angle_start",
                "angle_range",
                "angle_increment",
                "displacement_start",
                "displacement_range",
                "displacement_increment",
            ],
        )
        for name in axis_names:
            cbf.add_row([name, "SCAN1", "0", "0", "0", "0", "0", "0"])
        for row in gonio_scan_axis_rows:
            cbf.add_row(row)

        """Data items in the DIFFRN_SCAN_FRAME_AXIS category describe the
       settings of axes for particular frames.  Unspecified axes are
       assumed to be at their zero points."""
        cbf.add_category(
            "diffrn_scan_frame_axis", ["axis_id", "frame_id", "angle", "displacement"]
        )
        for row in axis_settings:
            cbf.add_row(row)
        for row in gonio_frame_axis_rows:
            cbf.add_row(row)

        """Data items in the ARRAY_STRUCTURE_LIST category record the size
       and organization of each array dimension.
       The relationship to physical axes may be given."""

        cbf.add_category(
            "array_structure_list",
            [
                "array_id",
                "index",
                "dimension",
                "precedence",
                "direction",
                "axis_set_id",
            ],
        )
        for panelkey, panelname in zip(panelkeys, panelnames):
            aname = level_string(panelkey)
            panel = panel_group_from_key(panelkey)
            fast_dim, slow_dim = panel.get_image_size()
            cbf.add_row(
                [
                    "ARRAY_" + panelname,
                    "1",
                    "%d" % fast_dim,
                    "1",
                    "increasing",
                    "AXIS_" + aname + "_F",
                ]
            )
            cbf.add_row(
                [
                    "ARRAY_" + panelname,
                    "2",
                    "%d" % slow_dim,
                    "2",
                    "increasing",
                    "AXIS_" + aname + "_S",
                ]
            )

        """Data items in the ARRAY_STRUCTURE_LIST_SECTION category identify
       the dimension-by-dimension start, end and stride of each section of an
       array that is to be referenced."""
        # no sections here

        """Data items in the ARRAY_STRUCTURE_LIST_AXIS category describe
       the physical settings of sets of axes for the centres of pixels that
       correspond to data points described in the
       ARRAY_STRUCTURE_LIST category."""
        cbf.add_category(
            "array_structure_list_axis",
            ["axis_set_id", "axis_id", "displacement", "displacement_increment"],
        )
        for panelkey in panelkeys:
            aname = level_string(panelkey)
            panel = panel_group_from_key(panelkey)
            pixel_size = panel.get_pixel_size()
            cbf.add_row(
                [
                    "AXIS_" + aname + "_F",
                    "AXIS_" + aname + "_F",
                    "0.0",
                    str(pixel_size[0]),
                ]
            )
            cbf.add_row(
                [
                    "AXIS_" + aname + "_S",
                    "AXIS_" + aname + "_S",
                    "0.0",
                    str(pixel_size[1]),
                ]
            )

        if not header_only:
            self.add_data_to_cbf(cbf, index)

        return cbf

    def add_data_to_cbf(self, cbf, index=None, data=None, verbose=False):
        """
        Given a cbf handle, add the raw data and the necessary tables to support it
        """
        if data is None:
            if index is None:
                data = self.imageset.get_raw_data(0)
            else:
                data = self.imageset.get_raw_data(index)
        if not isinstance(data, tuple):
            data = (data,)
        data = tuple([data[i] for i in self.panelindices])

        array_names = []
        cbf.find_category(b"diffrn_data_frame")
        while True:
            try:
                cbf.find_column(b"array_id")
                array_names.append(cbf.get_value())
                cbf.next_row()
            except Exception as e:
                e_message = str(e)

                assert "CBF_NOTFOUND" in e_message
                break

        dataisint = flex.bool()
        for panel_data in data:
            assert len(panel_data.focus()) == 2
            if isinstance(panel_data, flex.int):
                dataisint.append(True)
            elif isinstance(panel_data, flex.double):
                dataisint.append(False)
            else:
                raise TypeError("Ints or doubles are required")

        """ Data items in the ARRAY_STRUCTURE category record the organization and
    encoding of array data in the ARRAY_DATA category."""
        cbf.add_category(
            "array_structure", ["id", "encoding_type", "compression_type", "byte_order"]
        )
        for i, array_name in enumerate(array_names):
            if dataisint[i]:
                cbf.add_row(
                    [array_name, "signed 32-bit integer", "packed", "little_endian"]
                )
            else:
                cbf.add_row(
                    [array_name, "signed 64-bit real IEEE", "packed", "little_endian"]
                )

        """ Data items in the ARRAY_DATA category are the containers for the array data
    items described in the category ARRAY_STRUCTURE. """
        cbf.add_category("array_data", ["array_id", "binary_id", "data"])

        if verbose:
            print("Compressing tiles...")

        for i, (panel_data, array_name) in enumerate(zip(data, array_names)):
            focus = panel_data.focus()
            # panel_data += 1

            cbf.add_row([array_name, str(i + 1)])

            binary_id = i + 1
            byte_str = panel_data.copy_to_byte_str()
            elements = len(panel_data)
            byteorder = "little_endian"
            dimfast = focus[1]
            dimmid = focus[0]
            dimslow = 1
            padding = 0

            if dataisint[i]:
                elsize = 4
                elsigned = 1
                byteorder = byteorder.encode()
                cbf.set_integerarray_wdims_fs(
                    pycbf.CBF_PACKED,
                    binary_id,
                    byte_str,
                    elsize,
                    elsigned,
                    elements,
                    byteorder,
                    dimfast,
                    dimmid,
                    dimslow,
                    padding,
                )
            else:
                elsize = 8
                byteorder = byteorder.encode()
                cbf.set_realarray_wdims_fs(
                    pycbf.CBF_PACKED,
                    binary_id,
                    byte_str,
                    elsize,
                    elements,
                    byteorder,
                    dimfast,
                    dimmid,
                    dimslow,
                    padding,
                )

    def write_cbf(self, filename, index=None, cbf=None):
        """Write a CBF file. If the handle is not provided, create one"""
        assert [index, cbf].count(None) in (1, 2), "Supply either index or cbf"

        if cbf is None:
            cbf = self.get_cbf_handle(index=index, header_only=True)
            self.add_data_to_cbf(cbf, index=index)

        filename = filename.encode()
        cbf.write_widefile(
            filename, pycbf.CBF, pycbf.MIME_HEADERS | pycbf.MSG_DIGEST | pycbf.PAD_4K, 0
        )


if __name__ == "__main__":
    filename = sys.argv[1]
    if len(sys.argv) > 2:
        index: int | None = int(sys.argv[2])
    else:
        index = None

    writer = FullCBFWriter(filename=filename)
    writer.write_cbf("converted.cbf", index)
