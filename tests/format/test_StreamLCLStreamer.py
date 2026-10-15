"""Reader-boundary tests for the LCLStreamer stream class."""

from __future__ import annotations

import bitshuffle
import cbor
import numpy as np

from dxtbx.format.StreamLCLStreamer import LCLStreamer


def test_decode_normalizes_run_and_message_ids() -> None:
    """decode() casts the run_id string to an int and renames message_id to image_id."""
    reader = LCLStreamer()
    message = reader.decode(
        cbor.dumps({"type": "image", "run_id": "7", "message_id": 3})
    )
    assert message["type"] == "image"
    assert message["run_id"] == 7
    assert message["image_id"] == 3


def _frame() -> np.ndarray:
    # One Jungfrau module, 512 x 1024
    return np.arange(512 * 1024, dtype=np.float32).reshape(1, 512, 1024)


def test_uncompressed_frame_matches_compressed() -> None:
    """A frame sent uncompressed as "data" decodes to the same image as the same frame
    sent bitshuffle-lz4 compressed as "compressed_data"."""
    reader = LCLStreamer()
    frame = _frame()
    plain = reader.decode(cbor.dumps({"type": "image", "data": frame.tobytes()}))
    packed = reader.decode(
        cbor.dumps(
            {
                "type": "image",
                "compressed_data": bitshuffle.compress_lz4(
                    frame, block_size=2**12
                ).tobytes(),
            }
        )
    )
    image = reader._decode_image(plain, frame.shape, "float32")
    assert image.shape == (1, 514, 1030)
    assert np.array_equal(image, reader._decode_image(packed, frame.shape, "float32"))
    # Not the gap-marker frame
    assert image.max() > 0


def test_message_without_frame_is_gap_marker() -> None:
    """A message with neither "data" nor "compressed_data" is a fully-masked frame."""
    reader = LCLStreamer()
    image = reader._decode_image({"type": "image"}, (1, 512, 1024), "float32")
    assert image.max() < 0


def test_uncompressed_spectrometer() -> None:
    """The spectrum decodes from "spectrometer_data" (uncompressed) and from
    "spectrometer_compressed_data", and is absent when neither is sent."""
    reader = LCLStreamer()
    spectrum = np.arange(1, 9, dtype=np.float32)
    fields = {"spectrometer_dtype": "float32", "spectrometer_shape": "8"}
    plain = reader.decode(
        cbor.dumps({"type": "image", "spectrometer_data": spectrum.tobytes(), **fields})
    )
    packed = reader.decode(
        cbor.dumps(
            {
                "type": "image",
                "spectrometer_compressed_data": bitshuffle.compress_lz4(
                    spectrum, block_size=2**12
                ).tobytes(),
                **fields,
            }
        )
    )
    assert np.array_equal(reader._decode_spectrometer(plain), spectrum)
    assert np.array_equal(reader._decode_spectrometer(packed), spectrum)
    assert reader._decode_spectrometer(reader.decode(cbor.dumps(fields))) is None
