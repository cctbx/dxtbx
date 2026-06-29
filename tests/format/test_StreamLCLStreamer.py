"""Reader-boundary normalization tests for the LCLStreamer stream class."""

from __future__ import annotations

import cbor2

from dxtbx.format.StreamLCLStreamer import LCLStreamer


def test_decode_normalizes_stop_to_end() -> None:
    """LCLStreamer signals run-end with type "stop"; decode() must normalize it
    to the canonical "end" the rest of the system finalizes on."""
    reader = LCLStreamer()
    message = reader.decode(cbor2.dumps({"type": "stop", "run": 7}))
    assert message["type"] == "end"
    # Run-end messages carry "run", which decode() casts to an int run_id.
    assert message["run_id"] == 7


def test_decode_leaves_other_types_unchanged() -> None:
    """Only "stop" is rewritten; start/image types pass through verbatim."""
    reader = LCLStreamer()
    start = reader.decode(cbor2.dumps({"type": "start", "run_number": 7}))
    assert start["type"] == "start"
    assert start["run_id"] == 7

    image = reader.decode(cbor2.dumps({"type": "image", "run": 7, "message_id": 3}))
    assert image["type"] == "image"
    assert image["image_id"] == 3
