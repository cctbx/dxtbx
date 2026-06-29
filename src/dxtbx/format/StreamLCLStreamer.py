from __future__ import annotations

import datetime

import bitshuffle
import cbor2
import numpy as np
from PSCalib.GeometryAccess import GeometryAccess

from cctbx import factor_ev_angstrom, factor_kev_angstrom
from cctbx.eltbx import attenuation_coefficient
from scitbx.array_family import flex
from scitbx.matrix import col

from dxtbx.format.Stream import StreamClass
from dxtbx.model import Detector, ParallaxCorrectedPxMmStrategy, Spectrum
from dxtbx.model.experiment_list import Experiment, ExperimentList


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
    # Keep the psana root translation as-is: its z is already negative (downstream),
    # which is the dxtbx/dials backend convention. (An earlier version negated z here,
    # putting the detector on the wrong (+z) side; that was only made to look correct by
    # a compensating reflect_detector_for_display, which has since been removed.)

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
    # Return the 3-level (root -> module groups -> ASIC panels) tree in the backend
    # convention (slow (0, +1, 0) at downstream z < 0); the caller uses it as-is.
    return d


def get_jungfrau_detector_module(file_content, wavelength):
    """
    Copied from Fred Poitevin's work on FormatXTCJungfrau2M.py
    """
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
    # Keep the psana root translation as-is: its z is already negative (downstream),
    # which is the dxtbx/dials backend convention. (An earlier version negated z here,
    # putting the detector on the wrong (+z) side; that was only made to look correct by
    # a compensating reflect_detector_for_display, which has since been removed.)

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

    # Now deal with modules. Each module is a group holding a single full-module
    # panel, so the constructed tree (root -> module groups -> panel) matches the
    # 3-level shape of the refined reference geometries that sync_geometry pairs
    # against. The tree is in the backend convention (slow (0, +1, 0) at downstream
    # z < 0) and is used by the caller as-is (no display reflection).
    for module_num in range(len(root.get_list_of_children())):
        module = root.get_list_of_children()[module_num]
        module_basis = basis_from_geo(module)
        origin = col((module_basis * col((0, 0, 0, 1)))[0:3])
        fast = col((module_basis * col((1, 0, 0, 1)))[0:3]) - origin
        slow = col((module_basis * col((0, 1, 0, 1)))[0:3]) - origin
        pg1 = pg0.add_group()
        pg1.set_local_frame(fast.elems, slow.elems, origin.elems)
        pg1.set_name("D%dM%d" % (det_num, module_num))

        # Identity-framed panel under the module group: the module placement lives
        # on the group, so the panel's global frame equals the group frame.
        p = pg1.add_panel()
        p.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 0))
        p.set_name("ARRAY_D%dM%d" % (det_num, module_num))
        p.set_pixel_size((PIXEL_SIZE, PIXEL_SIZE))
        p.set_trusted_range(TRUSTED_RANGE)

        p.set_thickness(THICKNESS)  # mm
        p.set_material(MATERIAL)
        p.set_type("jungfrau")
        p.set_mu(mu)
        p.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, THICKNESS))
        p.set_gain(factor_kev_angstrom / wavelength)
        p.set_image_size((DIM_FAST, DIM_SLOW))
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
        tcp_keepalive=False,
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
            tcp_keepalive=tcp_keepalive,
        )
        self.name = "LCLStreamer"
        self._split_modules = False

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
        # Strip the leading topic/type byte the wire prepends.
        encoded_message = encoded_message[1:]
        return encoded_message

    def decode(self, encoded_message):
        message = cbor2.loads(encoded_message)
        # LCLStreamer signals run-end with type "stop"; the rest of the system
        # (ControlHub dispatch, finalize protocol) keys off the canonical "end".
        # Normalize here at the reader boundary so no component carries the
        # API-specific spelling.
        if message.get("type") == "stop":
            message["type"] = "end"
        # Translate the LCLStreamer wire-format run identifier to the internal
        # "run_id" the components consume. (Start messages carry "run_number",
        # image/end messages carry "run".)
        if "run" in message.keys():
            message["run_id"] = int(message.pop("run"))
        elif "run_number" in message.keys():
            message["run_id"] = int(message.pop("run_number"))
        if "message_id" in message.keys():
            message["image_id"] = message.pop("message_id")
        if "shape" in message.keys():
            if isinstance(message["shape"], str):
                message["image_shape"] = tuple(map(int, message["shape"].split("x")))
        # The optional per-frame spectrometer array (1D trace or 2D image) carries
        # its own shape string (e.g. "2048" or "512x1024"); normalize it to a tuple
        # so get_spectrum can decompress it. Absent for streams with no spectrometer.
        if "spectrometer_shape" in message.keys():
            if isinstance(message["spectrometer_shape"], str):
                message["spectrometer_shape"] = tuple(
                    map(int, message["spectrometer_shape"].split("x"))
                )
        if "datatype" in message.keys():
            message["image_dtype"] = message.pop("datatype")
        if "data_collection_rate" in message.keys():
            if isinstance(message["data_collection_rate"], str):
                message["data_collection_rate"] = float(
                    message["data_collection_rate"].split("Hz")[0]
                )
        return message

    def handle_start_message(
        self,
        message,
        reference_experiment=None,
        sync_reference_geom=True,
        wavelength=None,
    ):
        from dials.command_line.stills_process import sync_geometry

        from dxtbx.format.nxmx_writer import phil_scope as nxmx_writer_phil_scope
        from dxtbx.model.beam import BeamFactory, beam_phil_scope

        file_writer_params = nxmx_writer_phil_scope.extract()

        file_writer_params.dtype = message["image_dtype"]
        # The LCLS start message leaves "beamline" empty, so use the experiment
        # name as the instrument identifier (nxmx_writer.validate() requires a
        # non-empty instrument_name).
        file_writer_params.nexus_details.instrument_name = message["experiment"]
        file_writer_params.nexus_details.instrument_short_name = message["experiment"]
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
            beam_params.beam.wavelength = 1  # float(message["photon_wavelength"])
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
            # The builders now return the correct backend convention directly:
            # slow axis (0, +1, 0) at the true downstream distance (z < 0). No
            # reflect_detector_for_display is applied on either path. A refined
            # reference is already in that same convention, so when present it is
            # synced onto the matching-depth constructed tree and used AS-IS.
            # (A previous version built the placeholder on the wrong (+z) side and
            # applied a compensating global reflection; that band-aid double-flipped
            # the no-reference path into looking right but left a synced reference
            # vertically mirrored in dials.image_viewer. Both have been removed.)
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

        reference_experiment = ExperimentList(
            [Experiment(beam=beam, detector=detector)]
        )
        file_writer_params.detector.sensor_material = detector[0].get_material()
        file_writer_params.detector.sensor_thickness = detector[0].get_thickness()

        return file_writer_params, reference_experiment

    def _decode_image(self, message, image_shape, image_dtype):
        # Decompress + reshape the main detector frame (bitshuffle-lz4). Factored so
        # get_data and the single-decode get_frame_data share one implementation.
        image_data = bitshuffle.decompress_lz4(
            np.frombuffer(message["compressed_data"], dtype=np.uint8),
            shape=image_shape,
            dtype=np.dtype(image_dtype),
            block_size=2**12,
        )
        if self._split_modules:
            return self._reshape_jungfrau_asic(image_data)
        return self._reshape_jungfrau_module(image_data)

    def get_data(self, message, **kwargs):
        image_data = self._decode_image(
            message, kwargs["image_shape"], kwargs["image_dtype"]
        )
        # The (image_data, wavelength) tuple feeds the read-only consumers
        # (Display, viewers); they have no spectrometer calibration, so the
        # wavelength here is the non-spectrum fallback chain. Components that need the
        # spectrum decode it once via get_frame_data instead.
        return image_data, self.get_wavelength(message)

    def get_frame_data(self, message, image_shape=None, image_dtype=None, calib=None):
        # Single-decode entry point for components: decompress the main image AND the
        # spectrometer array exactly once, then derive every per-frame model from
        # them. Returns (image_data, wavelength, spectrum, spectrometer_image). Avoids
        # the 2-3x bitshuffle decode that calling get_wavelength + get_spectrum +
        # get_spectrometer_image separately would incur on each frame.
        calib = calib or {}
        image_data = self._decode_image(message, image_shape, image_dtype)
        spectrometer = self._decode_spectrometer(message)
        spectrum = self._build_spectrum(
            spectrometer,
            calib.get("spectrum_eV_per_pixel"),
            calib.get("spectrum_eV_offset"),
        )
        # A 1D trace is fully preserved on NXbeam via the spectrum; only a 2D image
        # needs the separate second-stream archive at native dimensionality.
        spectrometer_image = (
            spectrometer
            if spectrometer is not None and spectrometer.ndim == 2
            else None
        )
        wavelength = self._resolve_wavelength(message, spectrum, calib)
        return image_data, wavelength, spectrum, spectrometer_image

    def _resolve_wavelength(self, message, spectrum, calib):
        # Wavelength priority (mirrors FormatXTC._beam), given an already-built
        # spectrum so the spectrometer is not decoded again:
        #   1. a calibrated per-shot spectrum (its intensity-weighted wavelength),
        #   2. the eBeam photon_energy (eV) carried on the image message,
        #   3. the photon_wavelength PV (often 0/unreliable, hence after energy),
        #   4. the configured wavelength_fallback.
        # wavelength_scale/offset are applied only when the value did NOT come from a
        # spectrum (a spectrum already carries the true per-shot value).
        if spectrum is not None:
            return spectrum.get_weighted_wavelength()

        photon_energy = message.get("photon_energy")
        if photon_energy:
            wavelength = factor_ev_angstrom / photon_energy
        else:
            photon_wavelength = message.get("photon_wavelength")
            wavelength = photon_wavelength if photon_wavelength else None
        if not wavelength or wavelength <= 0:
            wavelength = calib.get("wavelength_fallback")
        if wavelength and wavelength > 0:
            if calib.get("wavelength_scale") is not None:
                wavelength *= calib["wavelength_scale"]
            if calib.get("wavelength_offset") is not None:
                wavelength += calib["wavelength_offset"]
            return wavelength
        return None

    def get_wavelength_candidates(self, message, spectrum=None, **calib):
        # The raw per-source wavelength estimates (Angstrom) that _resolve_wavelength
        # picks between, computed every frame so the losers can be archived next to the
        # resolved incident_wavelength for diagnostics. These are the uncorrected
        # per-sensor readings: the global wavelength_scale/offset applied to the chosen
        # fallback is deliberately NOT applied here, so each candidate reflects what its
        # sensor actually reported. A source absent on this frame is omitted.
        del calib  # candidates are uncorrected; signature mirrors the base resolver
        candidates = {}
        if spectrum is not None:
            candidates["spectrum"] = spectrum.get_weighted_wavelength()
        photon_energy = message.get("photon_energy")
        if photon_energy:
            candidates["ebeam"] = factor_ev_angstrom / photon_energy
        photon_wavelength = message.get("photon_wavelength")
        if photon_wavelength and photon_wavelength > 0:
            candidates["pv"] = photon_wavelength
        return candidates

    def get_wavelength(self, message, **calib):
        spectrum = self.get_spectrum(
            message,
            spectrum_eV_per_pixel=calib.get("spectrum_eV_per_pixel"),
            spectrum_eV_offset=calib.get("spectrum_eV_offset"),
        )
        return self._resolve_wavelength(message, spectrum, calib)

    def _decode_spectrometer(self, message):
        # Decompress the optional per-frame spectrometer array (1D trace or 2D
        # image), compressed bitshuffle-lz4 exactly like the main detector frame.
        # Returns a numpy array, or None when no spectrometer rode this message.
        payload = message.get("spectrometer_compressed_data")
        shape = message.get("spectrometer_shape")
        dtype = message.get("spectrometer_dtype")
        if payload is None or shape is None or dtype is None:
            return None
        return bitshuffle.decompress_lz4(
            np.frombuffer(payload, dtype=np.uint8),
            shape=tuple(shape),
            dtype=np.dtype(dtype),
            block_size=2**12,
        )

    def _build_spectrum(self, array, spectrum_eV_per_pixel, spectrum_eV_offset):
        # Build a dxtbx Spectrum from an already-decoded spectrometer array. Needs the
        # eV calibration; without it the raw pixels cannot be placed on an energy
        # axis, so return None and let the wavelength resolver fall back.
        if array is None:
            return None
        if spectrum_eV_per_pixel is None or spectrum_eV_offset is None:
            return None
        # Collapse a 2D spectrometer image to a 1D trace by averaging the slow axis
        # (matching FormatXTC._spectrum_hutch); a 1D trace is used as-is. The full
        # native-dimensionality pattern is archived separately; here we only need
        # the 1D distribution to weight the wavelength.
        trace = array.mean(axis=0) if array.ndim == 2 else array
        trace = np.asarray(trace, dtype=float)
        energies_eV = spectrum_eV_offset + spectrum_eV_per_pixel * np.arange(len(trace))
        # Guard the dxtbx.model.Spectrum invariant (weighted_sum > 0 and summed
        # weights > 0): a non-finite trace or one whose weights sum to <= 0 (e.g. an
        # over-subtracted background, or a partially-dropped reading) would fire a
        # fatal CCTBX_ASSERT in get_weighted_wavelength and take down the worker. Treat
        # such a degenerate reading as absent so the resolver falls back instead.
        if not np.isfinite(trace).all() or float(trace.sum()) <= 0:
            return None
        if float((energies_eV * trace).sum()) <= 0:
            return None
        return Spectrum(flex.double(energies_eV), flex.double(trace))

    def get_spectrum(
        self, message, spectrum_eV_per_pixel=None, spectrum_eV_offset=None, **kwargs
    ):
        return self._build_spectrum(
            self._decode_spectrometer(message),
            spectrum_eV_per_pixel,
            spectrum_eV_offset,
        )

    def get_spectrometer_image(self, message):
        # Archive the full 2D spectrometer pattern as a second image stream. A 1D
        # trace returns None: its full distribution is already stored on NXbeam by
        # get_spectrum, so there is nothing extra to keep at native dimensionality.
        array = self._decode_spectrometer(message)
        if array is None or array.ndim != 2:
            return None
        return array

    def get_data_scale_factor(self, wavelength: float) -> float:
        # LCLStreamer Jungfrau pixels are deposited energy in keV. Photons are
        # energy_keV / photon_energy_keV, where photon_energy_keV =
        # factor_kev_angstrom / wavelength. So the NeXus data_scale_factor (the
        # multiplier applied on read) is its reciprocal. The caller passes the
        # resolved wavelength: the wire's per-frame photon_wavelength is unreliable
        # (often 0), so the wavelength comes from the PHIL override / reference beam.
        return wavelength / factor_kev_angstrom

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
            for col_idx in range(4):
                # Calculate source position (256x256 asic)
                src_row_start = row * 256
                src_row_end = (row + 1) * 256
                src_col_start = col_idx * 256
                src_col_end = (col_idx + 1) * 256

                # Calculate destination position (with 2 pixel gaps)
                dst_row_start = src_row_start + 2 * row
                dst_row_end = dst_row_start + 256
                dst_col_start = src_col_start + 2 * col_idx
                dst_col_end = dst_col_start + 256

                # Copy ASIC data
                output[:, dst_row_start:dst_row_end, dst_col_start:dst_col_end] = (
                    image_data[:, src_row_start:src_row_end, src_col_start:src_col_end]
                )

        for col_idx in range(1, 4):
            left_index = 256 * col_idx + 2 * (col_idx - 1) - 1
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
