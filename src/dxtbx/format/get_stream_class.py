def get_stream_class(name):
    if name == "DectrisSimplonStreamV2":
        from dxtbx.format.StreamDectrisSimplonStreamV2 import (
            StreamDectrisSimplonStreamV2,
        )

        return StreamDectrisSimplonStreamV2
    elif name == "DectrisSimplonStreamV2Emulator":
        from dials_streaming.format_class.dectris_simplon_stream_v2 import (
            StreamDectrisSimplonStreamV2Emulator,
        )

        return StreamDectrisSimplonStreamV2Emulator
    elif name == "DectrisSimplonStreamV2HitFiltered":
        from dials_streaming.format_class.dectris_simplon_stream_v2 import (
            StreamDectrisSimplonStreamV2HitFiltered,
        )

        return StreamDectrisSimplonStreamV2HitFiltered

    elif name == "LCLStreamer":
        from dxtbx.format.StreamLCLStreamer import LCLStreamer
        return LCLStreamer
    elif name == "NXmxEmulator":
        from dials_streaming.format_class.nxmx_stream import NXmxEmulatorStreamer

        return NXmxEmulatorStreamer
    else:
        print(f"EXCEPTION: {name} IS NOT IMPLEMENTED")
        assert False
