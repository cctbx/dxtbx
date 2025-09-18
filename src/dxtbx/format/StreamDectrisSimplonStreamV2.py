import cbor2
import datetime
from dectris.compression import decompress
from dxtbx.format.Stream import StreamClass
import numpy as np


def decode_multi_dim_array(tag, column_major):
    dimensions, contents = tag.value
    if isinstance(contents, list):
        array = np.empty((len(contents),), dtype=object)
        array[:] = contents
    elif isinstance(contents, (np.ndarray, np.generic)):
        array = contents
    else:
        raise cbor2.CBORDecodeValueError("expected array or typed array")
    return array.reshape(dimensions, order="F" if column_major else "C")


def decode_typed_array(tag, dtype):
    if not isinstance(tag.value, bytes):
        raise cbor2.CBORDecodeValueError("expected byte string in typed array")
    return np.frombuffer(tag.value, dtype=dtype)


def decode_dectris_compression(tag):
    algorithm, elem_size, encoded = tag.value
    return decompress(encoded, algorithm, elem_size=elem_size)


tag_decoders = {
    40: lambda tag: decode_multi_dim_array(tag, column_major=False),
    64: lambda tag: decode_typed_array(tag, dtype="u1"),
    65: lambda tag: decode_typed_array(tag, dtype=">u2"),
    66: lambda tag: decode_typed_array(tag, dtype=">u4"),
    67: lambda tag: decode_typed_array(tag, dtype=">u8"),
    68: lambda tag: decode_typed_array(tag, dtype="u1"),
    69: lambda tag: decode_typed_array(tag, dtype="<u2"),
    70: lambda tag: decode_typed_array(tag, dtype="<u4"),
    71: lambda tag: decode_typed_array(tag, dtype="<u8"),
    72: lambda tag: decode_typed_array(tag, dtype="i1"),
    73: lambda tag: decode_typed_array(tag, dtype=">i2"),
    74: lambda tag: decode_typed_array(tag, dtype=">i4"),
    75: lambda tag: decode_typed_array(tag, dtype=">i8"),
    77: lambda tag: decode_typed_array(tag, dtype="<i2"),
    78: lambda tag: decode_typed_array(tag, dtype="<i4"),
    79: lambda tag: decode_typed_array(tag, dtype="<i8"),
    80: lambda tag: decode_typed_array(tag, dtype=">f2"),
    81: lambda tag: decode_typed_array(tag, dtype=">f4"),
    82: lambda tag: decode_typed_array(tag, dtype=">f8"),
    83: lambda tag: decode_typed_array(tag, dtype=">f16"),
    84: lambda tag: decode_typed_array(tag, dtype="<f2"),
    85: lambda tag: decode_typed_array(tag, dtype="<f4"),
    86: lambda tag: decode_typed_array(tag, dtype="<f8"),
    87: lambda tag: decode_typed_array(tag, dtype="<f16"),
    1040: lambda tag: decode_multi_dim_array(tag, column_major=True),
    56500: lambda tag: decode_dectris_compression(tag),
}


def tag_hook(decoder, tag):
    tag_decoder = tag_decoders.get(tag.tag)
    return tag_decoder(tag) if tag_decoder else tag


class StreamDectrisSimplonStreamV2(StreamClass):
    def __init__(
        self,
        port=None,
        ip_address=None,
        ipc_path=None,
        socket_library=None,
        socket_type=None,
        socket_mode=None,
        zmq_context=None,
        rcvhwm=None,
        rcvbuf=None,
    ):
        super().__init__(
            port,
            ip_address,
            ipc_path,
            socket_library,
            socket_type,
            socket_mode,
            zmq_context,
            rcvhwm,
            rcvbuf,
        )
        self.name = "DectrisSimplonStreamV2"

    def decode(self, encoded_message):
        return cbor2.loads(encoded_message, tag_hook=tag_hook)

    def _get_message_type(self, encoded_message=None, message=None):
        # len(encoded_message) gets the message size in bytes.
        # cbor2.loads takes ~0.5s with an Eiger 16M detector. len(encoded_message) is trivial.
        # The start and end messages are less than 1,000 bytes.
        # The image message is more than 1,000,000 bytes.
        # If the encoded message is more than 10,000 bytes, assume it is an image message and
        # do not decode it. Decoding the start and end messages is very fast.
        if message is None:
            type_section = encoded_message[:100]
            if b"image" in type_section:
                return "image"
            elif b"start" in type_section or b"end" in type_section:
                return "control"
            elif len(encoded_message) < 100000:
                return "control"
            else:
                return "image"
        else:
            return message["type"]

    def recv(self, copy=True, decode=True):
        encoded_message = self.socket.recv(copy=True)
        if decode:
            message = self.socket._decode(encoded_message)
            message_type = self._get_message_type(message)
            return message_type, message
        else:
            message_type = self._get_message_type(encoded_message)
            return message_type, encoded_message

    def handle_start_message(
        self, encoded_message=None, message=None, reference_experiment=None
    ):
        from dxtbx.format.nxmx_writer import phil_scope as nxmx_writer_phil_scope
        from dxtbx.format.nxmx_stream_writer import NXmxStreamWriter

        if message is None:
            message = self.decode(encoded_message)
        if isinstance(message["detector_description"], bytes):
            message["detector_description"] = message["detector_description"].decode()
        if isinstance(message["sensor_material"], bytes):
            message["sensor_material"] = message["sensor_material"].decode()
        if isinstance(message["image_dtype"], bytes):
            message["image_dtype"] = message["image_dtype"].decode()
        
        file_writer_params = nxmx_writer_phil_scope.extract()

        file_writer_params.dtype = message["image_dtype"]
        file_writer_params.nexus_details.instrument_name = message["detector_description"]
        file_writer_params.nexus_details.instrument_short_name = None
        file_writer_params.nexus_details.source_name = "ALS 2.0.1"
        file_writer_params.nexus_details.source_short_name = "ALS"

        if "start" in message.keys() and message["start"]:
            file_writer_params.nexus_details.start_time = message["start"]
        else:
            file_writer_params.nexus_details.start_time = datetime.datetime.now(
                datetime.UTC
            ).strftime("%Y-%m-%dT%H:%M:%SZ")
        collection_time_estimate = message["frame_time"] * message["number_of_images"]

        # This is wrong. The first term should be the params.nexus_details.start_time
        end_time_estimated = datetime.datetime.now(datetime.UTC) + datetime.timedelta(
            seconds=collection_time_estimate
        )
        file_writer_params.nexus_details.end_time_estimated = (
            end_time_estimated.strftime("%Y-%m-%dT%H:%M:%SZ")
        )
        file_writer_params.nexus_details.count_time = message["count_time"]
        file_writer_params.nexus_details.frame_time = message["frame_time"]
        file_writer_params.nexus_details.sample_name = "sample"
        file_writer_params.detector.sensor_material = message["sensor_material"]
        file_writer_params.detector.sensor_thickness = message["sensor_thickness"]

        if reference_experiment is None:
            from dxtbx.model.beam import beam_phil_scope
            from dxtbx.model.beam import BeamFactory
            from dxtbx.model.detector import detector_phil_scope
            from dxtbx.model.detector import DetectorFactory
            from dxtbx.model.experiment_list import Experiment

            # Construct beam
            beam_params = beam_phil_scope.extract()
            beam_params_ = beam_phil_scope.format(beam_params).as_str()
            beam_params.beam.direction = [0, 0, 1]
            beam_params.beam.divergence = None
            beam_params.beam.flux = None
            beam_params.beam.polarization_fraction = None
            beam_params.beam.polarization_normal = None
            beam_params.beam.probe = None
            beam_params.beam.sample_to_source_distance = 0
            beam_params.beam.sigma_divergence = None
            beam_params.beam.transmission = None
            beam_params.beam.type = None
            beam_params.beam.wavelength = message["incident_wavelength"]
            beam_params.beam.wavelength_range = None

            # Construct detector
            detector_params = detector_phil_scope.extract()
            # detector_params.detector.distance
            # detector_params.detector.fast_slow_beam_centre
            # detector_params.detector.mosflm_beam_centre
            # detector_params.detector.slow_fast_beam_centre

            detector_params.detector.hierarchy.fast_axis = (1, 0, 0)
            detector_params.detector.hierarchy.origin = (0, 0, 0)
            detector_params.detector.hierarchy.slow_axis = (0, -1, 0)

            # detector_params.detector.hierarchy.group.fast_axis
            # detector_params.detector.hierarchy.group.origin
            # detector_params.detector.hierarchy.group.id
            # detector_params.detector.hierarchy.group.panel
            # detector_params.detector.hierarchy.group.name
            # detector_params.detector.hierarchy.group.slow_axis

            detector_params.detector.panel[0].fast_axis = (1, 0, 0)
            detector_params.detector.panel[0].gain = 1.0
            detector_params.detector.panel[0].image_size = (
                message["image_size_x"],
                message["image_size_y"],
            )
            detector_params.detector.panel[0].material = message["sensor_material"]
            detector_params.detector.panel[0].origin = (
                -1000 * message["detector_translation"][0],
                1000 * message["detector_translation"][1],
                1000 * message["detector_translation"][2],
            )
            detector_params.detector.panel[0].pedestal = 0.0
            detector_params.detector.panel[0].pixel_size = (
                1000 * message["pixel_size_x"],
                1000 * message["pixel_size_y"],
            )
            detector_params.detector.panel[0].slow_axis = (0, -1, 0)
            detector_params.detector.panel[0].thickness = message["sensor_thickness"]
            detector_params.detector.panel[0].trusted_range = (0, 2147483647)

            reference_experiment = Experiment(
                beam=BeamFactory.from_phil(beam_params),
                detector=DetectorFactory.generate_from_phil(detector_params, ref_beam),
            )

        file_writer = NXmxStreamWriter(file_writer_params)
        file_writer(experiments=reference_experiment, in_memory=True)
        return file_writer

    def handle_end_message(self, encoded_message=None, message=None):
        if message is None:
            message = self.decode(encoded_message)
        return message

    def handle_image_message(self, encoded_message=None, message=None):
        if message is None:
            message = self.decode(encoded_message)
        return message

    def get_data(self, message, **kwargs):
        return message["data"]["threshold_1"]
