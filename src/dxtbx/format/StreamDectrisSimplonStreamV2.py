from __future__ import annotations

import datetime
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple, Union

import cbor2
import numpy as np
from dectris.compression import decompress

from dxtbx.format.Stream import StreamClass
from dxtbx.model.experiment_list import Experiment, ExperimentList

if TYPE_CHECKING:
    import zmq

    from libtbx.phil import scope_extract

    from dxtbx.imageset import StreamReader


def decode_multi_dim_array(tag: cbor2.CBORTag, column_major: bool) -> np.ndarray:
    dimensions, contents = tag.value
    if isinstance(contents, list):
        array = np.empty((len(contents),), dtype=object)
        array[:] = contents
    elif isinstance(contents, (np.ndarray, np.generic)):
        array = contents
    else:
        raise cbor2.CBORDecodeValueError("expected array or typed array")
    return array.reshape(dimensions, order="F" if column_major else "C")


def decode_typed_array(tag: cbor2.CBORTag, dtype: Any) -> np.ndarray:
    if not isinstance(tag.value, bytes):
        raise cbor2.CBORDecodeValueError("expected byte string in typed array")
    return np.frombuffer(tag.value, dtype=dtype)


def decode_dectris_compression(tag: cbor2.CBORTag) -> bytes:
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


def tag_hook(arg0: Any, arg1: Any = None) -> Any:
    # cbor2 >= 6 calls tag_hook(tag, immutable); older cbor2 calls
    # tag_hook(decoder, tag). Select whichever argument is the CBORTag.
    tag = arg0 if isinstance(arg0, cbor2.CBORTag) else arg1
    tag_decoder = tag_decoders.get(tag.tag)
    return tag_decoder(tag) if tag_decoder else tag


class StreamDectrisSimplonStreamV2(StreamClass):
    def __init__(
        self,
        port: Optional[int] = None,
        ip_address: Optional[str] = None,
        zmq_context: Optional[zmq.Context] = None,
        rcvhwm: Optional[int] = None,
        rcvbuf: Optional[int] = None,
        tcp_keepalive: bool = False,
    ) -> None:
        super().__init__(
            port=port,
            ip_address=ip_address,
            zmq_context=zmq_context,
            rcvhwm=rcvhwm,
            rcvbuf=rcvbuf,
            tcp_keepalive=tcp_keepalive,
        )
        self.name = "dectris"

    def decode(self, encoded_message: bytes) -> Dict[str, Any]:
        # cbor2 >= 6 returns an immutable frozendict when the message is wrapped
        # in the self-describe tag (55799); rebuild a mutable dict so the
        # in-place updates below (and downstream) work as they did before.
        message = dict(cbor2.loads(encoded_message, tag_hook=tag_hook))
        if "image_size_x" in message.keys():
            message["image_shape"] = (message["image_size_y"], message["image_size_x"])
        if "series_id" in message and "run_id" not in message:
            message["run_id"] = message.pop("series_id")
        # The dtype is supposedly a function of frame rate. This will be useful
        # once those boundaries are determined.
        if message["type"] == "start":
            if message["count_time"] >= 0.0145:
                message.setdefault("image_dtype", "uint32")
            elif message["count_time"] >= 0.00:
                message.setdefault("image_dtype", "uint16")
            else:
                message.setdefault("image_dtype", "uint8")

        return message

    def recv(self, copy: bool = True) -> Union[bytes, memoryview]:
        return self.socket.recv(copy=copy)

    def handle_start_message(
        self,
        message: Dict[str, Any],
        reference_experiment: Optional[ExperimentList] = None,
        sync_reference_geom: bool = True,
        wavelength: Optional[float] = None,
    ) -> Tuple[scope_extract, ExperimentList]:
        from dials.command_line.stills_process import sync_geometry

        from dxtbx.format.nxmx_writer import phil_scope as nxmx_writer_phil_scope
        from dxtbx.model.beam import BeamFactory, beam_phil_scope

        if isinstance(message["detector_description"], bytes):
            message["detector_description"] = message["detector_description"].decode()
        if isinstance(message["sensor_material"], bytes):
            message["sensor_material"] = message["sensor_material"].decode()
        if isinstance(message["image_dtype"], bytes):
            message["image_dtype"] = message["image_dtype"].decode()

        file_writer_params = nxmx_writer_phil_scope.extract()

        file_writer_params.dtype = message["image_dtype"]
        file_writer_params.nexus_details.instrument_name = message[
            "detector_description"
        ]
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
        file_writer_params.detector.sensor_thickness = (
            1000 * message["sensor_thickness"]
        )

        # Construct beam
        beam_params = beam_phil_scope.extract()
        beam_params.beam.direction = [0, 0, 1]
        beam_params.beam.divergence = None
        beam_params.beam.flux = None
        beam_params.beam.polarization_fraction = None
        beam_params.beam.polarization_normal = None
        beam_params.beam.probe = "x-ray"
        beam_params.beam.sample_to_source_distance = 0
        beam_params.beam.sigma_divergence = None
        beam_params.beam.transmission = None
        beam_params.beam.type = "monochromatic"
        beam_params.beam.wavelength = float(message["incident_wavelength"])
        beam_params.beam.wavelength_range = None
        beam = BeamFactory.from_phil(beam_params)

        if reference_experiment is None or sync_reference_geom:
            from cctbx.eltbx import attenuation_coefficient

            from dxtbx.model import ParallaxCorrectedPxMmStrategy
            from dxtbx.model.detector import Detector

            # Construct detector via Detector.from_dict rather than
            # DetectorFactory.generate_from_phil. With a hierarchy.group,
            # generate_from_phil builds a 3-level tree (root -> group -> panel),
            # but a refined reference .expt with one panel is a 2-level tree
            # (root -> panel). sync_geometry recursively zips src/dest children and
            # copies *local* frames assuming identical structure, so the depth
            # mismatch makes it pin the reference panel's frame onto our
            # intermediate group and skip the real panel; the untouched panel then
            # composes through the now-flipped group and its global slow axis comes
            # out (0, +1, 0) instead of (0, -1, 0). from_dict yields a 2-level
            # root -> panel tree matching the reference, so sync_geometry pairs
            # panel-with-panel and copies cleanly. (The no-reference path is
            # unaffected: the constructed global frame is correct either way.)
            detector = Detector.from_dict(
                {
                    "panels": [
                        {
                            "name": "panel0",
                            "type": "SENSOR_PAD",
                            "fast_axis": (1.0, 0.0, 0.0),
                            "slow_axis": (0.0, -1.0, 0.0),
                            "origin": (
                                -1000 * message["detector_translation"][0],
                                1000 * message["detector_translation"][1],
                                1000 * message["detector_translation"][2],
                            ),
                            "image_size": (
                                message["image_size_x"],
                                message["image_size_y"],
                            ),
                            "pixel_size": (
                                1000 * message["pixel_size_x"],
                                1000 * message["pixel_size_y"],
                            ),
                            "trusted_range": (0, 2147483647),
                            "thickness": 1000 * message["sensor_thickness"],
                            "material": message["sensor_material"],
                            "gain": 1.0,
                            "pedestal": 0.0,
                            "mask": [],
                        }
                    ],
                    "hierarchy": {
                        "fast_axis": (1.0, 0.0, 0.0),
                        "slow_axis": (0.0, 1.0, 0.0),
                        "origin": (0.0, 0.0, 0.0),
                        "children": [{"panel": 0}],
                    },
                }
            )
            # Apply parallax correction so the px<->mm mapping matches the
            # file-based path (FormatCBFMini / nexus.build_detector). Detector.from_dict
            # leaves the default SimplePxMmStrategy, which ignores the depth at which
            # X-rays are absorbed in the sensor and shifts spot positions.
            mu = (
                attenuation_coefficient.get_table(
                    detector[0].get_material()
                ).mu_at_angstrom(beam.get_wavelength())
                / 10.0
            )  # cm^-1 -> mm^-1
            detector[0].set_mu(mu)
            detector[0].set_px_mm_strategy(
                ParallaxCorrectedPxMmStrategy(mu, detector[0].get_thickness())
            )
            if reference_experiment and sync_reference_geom:
                sync_geometry(
                    reference_experiment[0].detector.hierarchy(),
                    detector.hierarchy(),
                )

            reference_experiment = ExperimentList(
                [Experiment(beam=beam, detector=detector)]
            )
        else:
            # If the reference_experiment has an imageset, it gets removed by
            # creating a new experiment without the imageset.
            reference_experiment = ExperimentList(
                [
                    Experiment(
                        beam=beam,
                        detector=reference_experiment[0].detector,
                    )
                ]
            )
        reference_experiment = ExperimentList(
            [
                Experiment(
                    beam=reference_experiment[0].beam,
                    detector=reference_experiment[0].detector,
                )
            ]
        )
        return file_writer_params, reference_experiment

    def get_data(
        self, message: Dict[str, Any], **kwargs: Any
    ) -> Tuple[np.ndarray, Optional[float]]:
        image_data = message["data"]["threshold_1"]
        return image_data, None

    def get_reader(self, image_data: np.ndarray, **kwargs: Any) -> StreamReader:
        from dials.array_family import flex

        from dxtbx.imageset import StreamReader

        # if 32 bit then it is a signed int, I think if 8, 16 then it is
        # unsigned with the highest two values assigned as masking values
        if image_data.dtype.name == "uint32":
            top = 2**32
        elif image_data.dtype.name == "uint16":
            top = 2**16
        elif image_data.dtype.name == "uint8":
            top = 2**8
        else:
            raise Exception(
                f"Unhandled data type {image_data.dtype.name} in StreamDectrisSimplonStreamV2.get_data"
            )
        image_data_int = image_data.astype(np.int32)
        image_data_int[image_data == (top - 1)] = -1
        image_data_int[image_data == (top - 2)] = -2

        image_data = flex.double(image_data_int)
        return StreamReader([image_data])
