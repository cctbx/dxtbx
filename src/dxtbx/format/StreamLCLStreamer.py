import bitshuffle
import cbor2
from cctbx import factor_kev_angstrom
import datetime
from dxtbx.model import Detector, ParallaxCorrectedPxMmStrategy
from cctbx.eltbx import attenuation_coefficient
from dxtbx.model.experiment_list import ExperimentList, Experiment
from dxtbx.format.Stream import StreamClass
import numpy as np
from PSCalib.GeometryAccess import GeometryAccess
from scitbx.matrix import col


def get_jungfrau_detector_asic(file_content, wavelength):
    try:
        from PSCalib.SegGeometryStore import sgs
    except ModuleNotFoundError:
        # psana2
        from psana.pscalib.geometry.SegGeometryStore import sgs

    from serialtbx.detector.xtc import basis_from_geo

    PIXEL_SIZE = 0.075
    TRUSTED_RANGE = (-10, 2e6)
    THICKNESS = 0.32
    MATERIAL = "Si"

    geom = GeometryAccess()
    geom.load_pars_from_str(file_content)

    d = Detector()
    pg0 = d.hierarchy()
    # first deal with D0
    det_num = 0
    root = geom.get_top_geo()
    root_basis = basis_from_geo(root)
    while len(root.get_list_of_children()) == 1:
        sub = root.get_list_of_children()[0]
        sub_basis = basis_from_geo(sub)
        root = sub
        root_basis = root_basis * sub_basis
    t = root_basis.translation
    distance = t[2]
    root_basis.translation = col((t[0], t[1], -distance))

    origin = col((root_basis * col((0, 0, 0, 1)))[0:3])
    fast = col((root_basis * col((1, 0, 0, 1)))[0:3]) - origin
    slow = col((root_basis * col((0, 1, 0, 1)))[0:3]) - origin

    ###!!! This rotation is unexplained
    normal = fast.cross(slow)
    rotation = normal.axis_and_angle_as_r3_rotation_matrix(-90, deg=True)
    fast = rotation * fast
    slow = rotation * slow
    origin = rotation * origin

    pg0.set_local_frame(fast.elems, slow.elems, origin.elems)
    pg0.set_name("D%d" % (det_num))

    # Now deal with modules
    for module_num in range(len(root.get_list_of_children())):
        module = root.get_list_of_children()[module_num]
        module_basis = basis_from_geo(module)
        origin = col((module_basis * col((0, 0, 0, 1)))[0:3])
        fast = col((module_basis * col((1, 0, 0, 1)))[0:3]) - origin
        slow = col((module_basis * col((0, 1, 0, 1)))[0:3]) - origin
        pg1 = pg0.add_group()
        pg1.set_local_frame(fast.elems, slow.elems, origin.elems)
        pg1.set_name("D%dM%d" % (det_num, module_num))

        # Read the known layout of the Jungfrau 2x4 module
        sg = sgs.Create(segname=module.oname)
        xx, yy = sg.get_seg_xy_maps_um()
        xx = xx / 1000
        yy = yy / 1000

        # Now deal with ASICs
        ###!!! Information about the individual asics is not included in the geometry
        ###!!! and is hard-coded
        for asic_num in range(8):
            val = "ARRAY_D0M%dA%d" % (module_num, asic_num)
            dim_slow = xx.shape[0]
            dim_fast = xx.shape[1]
            sensor_id = asic_num // 4  # There are 2X4 asics per module
            asic_in_sensor_id = asic_num % 4  # this number will be 0,1,2 or 3
            id_slow = sensor_id * (dim_slow // 2)
            id_fast = asic_in_sensor_id * (dim_fast // 4)
            origin = col((xx[id_slow][id_fast], yy[id_slow][id_fast], 0))
            fp = col((xx[id_slow][id_fast + 1], yy[id_slow][id_fast + 1], 0))
            sp = col((xx[id_slow + 1][id_fast], yy[id_slow + 1][id_fast], 0))
            fast = (fp - origin).normalize()
            slow = (sp - origin).normalize()
            p = pg1.add_panel()
            p.set_local_frame(fast.elems, slow.elems, origin.elems)
            p.set_pixel_size((PIXEL_SIZE, PIXEL_SIZE))
            p.set_trusted_range(TRUSTED_RANGE)
            p.set_name(val)

            p.set_thickness(THICKNESS)  # mm
            p.set_material(MATERIAL)
            p.set_type("jungfrau")
            # Compute the attenuation coefficient.
            # This will fail for undefined composite materials
            # mu_at_angstrom returns cm^-1, but need mu in mm^-1
            table = attenuation_coefficient.get_table(MATERIAL)
            mu = table.mu_at_angstrom(wavelength) / 10.0
            p.set_mu(mu)
            p.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, THICKNESS))
            p.set_gain(factor_kev_angstrom / wavelength)
            p.set_image_size((dim_fast // 4, dim_slow // 2))
    return d


def get_jungfrau_detector_module(file_content, wavelength):
    """
    Copied from Fred Poitevin's work on FormatXTCJungfrau2M.py
    """
    try:
        from PSCalib.SegGeometryStore import sgs
    except ModuleNotFoundError:
        # psana2
        from psana.pscalib.geometry.SegGeometryStore import sgs

    from serialtbx.detector.xtc import basis_from_geo

    PIXEL_SIZE = 0.075
    TRUSTED_RANGE = (-10, 2e6)
    THICKNESS, MATERIAL = 0.32, "Si"
    DIM_FAST = 1030
    DIM_SLOW = 514
    # Compute the attenuation coefficient.
    # This will fail for undefined composite materials
    # mu_at_angstrom returns cm^-1, but need mu in mm^-1
    table = attenuation_coefficient.get_table(MATERIAL)
    mu = table.mu_at_angstrom(wavelength) / 10.0

    geom = GeometryAccess()
    geom.load_pars_from_str(file_content)

    d = Detector()
    pg0 = d.hierarchy()
    # first deal with D0
    det_num = 0
    root = geom.get_top_geo()
    root_basis = basis_from_geo(root)
    while len(root.get_list_of_children()) == 1:
        sub = root.get_list_of_children()[0]
        sub_basis = basis_from_geo(sub)
        root = sub
        root_basis = root_basis * sub_basis
    t = root_basis.translation
    distance = t[2]
    root_basis.translation = col((t[0], t[1], -distance))

    origin = col((root_basis * col((0, 0, 0, 1)))[0:3])
    fast = col((root_basis * col((1, 0, 0, 1)))[0:3]) - origin
    slow = col((root_basis * col((0, 1, 0, 1)))[0:3]) - origin

    ###!!! This rotation is unexplained
    normal = fast.cross(slow)
    rotation = normal.axis_and_angle_as_r3_rotation_matrix(-90, deg=True)
    fast = rotation * fast
    slow = rotation * slow
    origin = rotation * origin

    pg0.set_local_frame(fast.elems, slow.elems, origin.elems)
    pg0.set_name("D%d" % (det_num))

    # Now deal with modules
    for module_num in range(len(root.get_list_of_children())):
        module = root.get_list_of_children()[module_num]
        module_basis = basis_from_geo(module)
        origin = col((module_basis * col((0, 0, 0, 1)))[0:3])
        fast = col((module_basis * col((1, 0, 0, 1)))[0:3]) - origin
        slow = col((module_basis * col((0, 1, 0, 1)))[0:3]) - origin
        pg1 = pg0.add_panel()
        pg1.set_local_frame(fast.elems, slow.elems, origin.elems)
        pg1.set_name("D%dM%d" % (det_num, module_num))

        pg1.set_pixel_size((PIXEL_SIZE, PIXEL_SIZE))
        pg1.set_trusted_range(TRUSTED_RANGE)

        pg1.set_thickness(THICKNESS)  # mm
        pg1.set_material(MATERIAL)
        pg1.set_type("jungfrau")
        pg1.set_mu(mu)
        pg1.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, THICKNESS))
        pg1.set_gain(factor_kev_angstrom / wavelength)
        pg1.set_image_size((DIM_FAST, DIM_SLOW))
    return d


class LCLStreamer(StreamClass):
    """
    https://github.com/slac-lcls/lclstreamer/blob/features/simplon/src/lclstreamer/frontend/data_serializers/json.py
    https://confluence.slac.stanford.edu/spaces/PSDM/pages/150405475/Detector+Geometry
    """

    def __init__(
        self,
        port=None,
        ports=None,
        ip_address=None,
        socket_library=None,
        socket_type=None,
        socket_mode=None,
        zmq_context=None,
        rcvhwm=None,
        rcvbuf=None,
    ):
        super().__init__(
            port=port,
            ports=ports,
            ip_address=ip_address,
            socket_library=socket_library,
            socket_type=socket_type,
            socket_mode=socket_mode,
            zmq_context=zmq_context,
            rcvhwm=rcvhwm,
            rcvbuf=rcvbuf,
        )
        self.name = "LCLStreamer"
        self._split_modules = True

    def recv(self, copy=True):
        # The LCLStream adds a character at the beginning of the message to identify
        # The message type withouit decoding.
        #   b"c" == control message
        #   b"m" == image message
        if self.socket_library == "zmq":
            encoded_message = self.socket.recv(copy=True)
        elif self.socket_library == "nng":
            # pynng recv() returns a bytes object directly, no copy parameter
            encoded_message = self.socket.recv()
        message_type = encoded_message[:1]
        encoded_message = encoded_message[1:]
        return encoded_message

    def decode(self, encoded_message):
        message = cbor2.loads(encoded_message)
        if "run" in message.keys():
            message["series_id"] = int(message.pop("run"))
        elif "run_number" in message.keys():
            message["series_id"] = message.pop("run_number")
            message["series_id"] = int(message["series_id"])
        if "message_id" in message.keys():
            message["image_id"] = message.pop("message_id")
        if "shape" in message.keys():
            if type(message["shape"]) == str:
                message["image_shape"] = tuple(map(int, message["shape"].split("x")))
        if "datatype" in message.keys():
            message["image_dtype"] = message.pop("datatype")
        if "data_collection_rate" in message.keys():
            if type(message["data_collection_rate"]) == str:
                message["data_collection_rate"] = float(
                    message["data_collection_rate"].split("Hz")[0]
                )
        return message

    def handle_start_message(self, message, reference_experiment=None, sync_reference_geom=True, wavelength=None):
        from dials.command_line.stills_process import sync_geometry
        from dxtbx.format.nxmx_writer import phil_scope as nxmx_writer_phil_scope
        from dxtbx.model.beam import beam_phil_scope
        from dxtbx.model.beam import BeamFactory

        file_writer_params = nxmx_writer_phil_scope.extract()

        file_writer_params.dtype = message["image_dtype"]
        file_writer_params.nexus_details.instrument_name = message["beamline"]
        file_writer_params.nexus_details.instrument_short_name = message["beamline"]
        file_writer_params.nexus_details.source_name = "LCLS"
        file_writer_params.nexus_details.source_short_name = "LCLS"

        file_writer_params.nexus_details.start_time = message["start_time"]
        if message["duration"] == "N/A":
            collection_time_estimate = 5 * 60  # five minute run
        else:
            collection_time_estimate = message["duration"]

        # This is wrong. The first term should be the params.nexus_details.start_time
        end_time_estimated = datetime.datetime.now(datetime.UTC) + datetime.timedelta(
            seconds=collection_time_estimate
        )
        file_writer_params.nexus_details.end_time_estimated = (
            end_time_estimated.strftime("%Y-%m-%dT%H:%M:%SZ")
        )
        file_writer_params.nexus_details.count_time = (
            1 / message["data_collection_rate"]
        )
        file_writer_params.nexus_details.frame_time = (
            1 / message["data_collection_rate"]
        )
        file_writer_params.nexus_details.sample_name = message["experiment"]

        # Construct beam
        beam_params = beam_phil_scope.extract()
        beam_params.beam.direction = [0, 0, 1]
        beam_params.beam.divergence = None
        beam_params.beam.flux = None
        beam_params.beam.polarization_fraction = message["polarization"]["fraction"]
        beam_params.beam.polarization_normal = message["polarization"]["axis"]
        beam_params.beam.probe = message["beam_type"].lower()
        beam_params.beam.sample_to_source_distance = 0
        beam_params.beam.sigma_divergence = None
        beam_params.beam.transmission = None
        beam_params.beam.type = "monochromatic"
        if wavelength:
            beam_params.beam.wavelength = wavelength
        else:
            # Wavelength is not included in the start message ...
            beam_params.beam.wavelength = 1 # float(message["photon_wavelength"])
        beam_params.beam.wavelength_range = None
        beam = BeamFactory.from_phil(beam_params)

        if reference_experiment is None or sync_reference_geom:
            # Construct detector
            if "jungfrau" in message["detector"]["name"].lower():
                if self._split_modules:
                    detector = get_jungfrau_detector_asic(
                        message["detector"]["geometry"], beam_params.beam.wavelength
                    )
                else:
                    detector = get_jungfrau_detector_module(
                        message["detector"]["geometry"], beam_params.beam.wavelength
                    )
            else:
                assert False
            if reference_experiment and sync_reference_geom:
                sync_geometry(
                    reference_experiment[0].detector.hierarchy(),
                    detector.hierarchy(),
                )
        else:
            # Add the type of detector to the panel type field. This gets used
            # by the get_reader method to determine how to read the data.
            if "jungfrau" in message["detector"]["name"].lower():
                for panel in reference_experiment[0].detector:
                    panel.set_type("jungfrau")
            else:
                assert False
            # If the reference_experiment has an imageset, it gets removed by
            # creating a new experiment without the imageset.
            beam = reference_experiment[0].beam
            detector = reference_experiment[0].detector

        reference_experiment = ExperimentList([Experiment(beam=beam, detector=detector)])
        file_writer_params.detector.sensor_material = detector[0].get_material()
        file_writer_params.detector.sensor_thickness = detector[0].get_thickness()

        return file_writer_params, reference_experiment

    def get_data(self, message, **kwargs):
        image_data = bitshuffle.decompress_lz4(
            np.frombuffer(message["compressed_data"], dtype=np.uint8),
            shape=kwargs["image_shape"],
            dtype=np.dtype(kwargs["image_dtype"]),
            block_size=2**12,
        )

        if self._split_modules:
            image_data = self._reshape_jungfrau_asic(image_data)
        else:
            image_data = self._reshape_jungfrau_module(image_data)
        wavelength = message["photon_wavelength"]
        return image_data, wavelength

    def _reshape_jungfrau_asic(self, image_data):
        # JF16M: (32 x 512 x 1024) -> (256 x 256 x 256)
        n_modules = image_data.shape[0]
        # Reshape to (n_modules, 2, 256, 4, 256)
        reshaped = image_data.reshape(n_modules, 2, 256, 4, 256)
        # Transpose to (n_modules, 2, 4, 256, 256)
        transposed = reshaped.transpose(0, 1, 3, 2, 4)
        # Reshape to (n_modules * 8, 256, 256)
        return transposed.reshape(n_modules * 8, 256, 256)

    def _reshape_jungfrau_module(self, image_data):
        # JF16M: (32 x 512 x 1024) -> (32 x 514 x 1030)
        n_modules = image_data.shape[0]
        # Create output array with gaps
        output = np.zeros((n_modules, 514, 1030), dtype=image_data.dtype)
        for row in range(2):
            for col in range(4):
                # Calculate source position (256x256 asic)
                src_row_start = row * 256
                src_row_end = (row + 1) * 256
                src_col_start = col * 256
                src_col_end = (col + 1) * 256

                # Calculate destination position (with 2 pixel gaps)
                dst_row_start = src_row_start + 2 * row
                dst_row_end = dst_row_start + 256
                dst_col_start = src_col_start + 2 * col
                dst_col_end = dst_col_start + 256

                # Copy ASIC data
                output[:, dst_row_start:dst_row_end, dst_col_start:dst_col_end] = (
                    image_data[:, src_row_start:src_row_end, src_col_start:src_col_end]
                )

        for col in range(1, 4):
            left_index = 256 * col + 2 * (col - 1) - 1
            right_index = left_index + 3
            left = output[:, :, left_index]
            right = output[:, :, right_index]
            output[:, :, left_index + 1 : right_index] = ((left + right) / 4)[
                :, :, np.newaxis
            ]
            output[:, :, left_index] = left / 2
            output[:, :, right_index] = right / 2

        top_index = 256 - 1
        bottom_index = top_index + 3
        top = output[:, top_index, :]
        bottom = output[:, bottom_index, :]
        output[:, top_index + 1 : bottom_index] = ((top + bottom) / 4)[:, np.newaxis, :]
        output[:, top_index] = top / 2
        output[:, bottom_index] = bottom / 2

        return output

    def get_reader(self, image_data, **kwargs):
        from dials.array_family import flex
        from dxtbx.imageset import StreamReader

        image_data = tuple(
            [flex.double(image_data[i]) for i in range(image_data.shape[0])]
        )
        return StreamReader([image_data])
