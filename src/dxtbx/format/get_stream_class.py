from __future__ import annotations


def get_stream_class(name):
    if name == "dectris":
        from dxtbx.format.StreamDectrisSimplonStreamV2 import (
            StreamDectrisSimplonStreamV2,
        )

        return StreamDectrisSimplonStreamV2
    elif name == "dectris_emulator":
        from dials_streaming.format_class.dectris_simplon_stream_v2 import (
            StreamDectrisSimplonStreamV2Emulator,
        )

        return StreamDectrisSimplonStreamV2Emulator
    elif name == "lcls":
        from dxtbx.format.StreamLCLStreamer import LCLStreamer

        return LCLStreamer
    elif name == "nxmx":
        from dials_streaming.format_class.nxmx_stream import NXmxEmulatorStreamer

        return NXmxEmulatorStreamer
    else:
        print(f"EXCEPTION: {name} IS NOT IMPLEMENTED")
        assert False
