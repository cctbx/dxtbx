def get_stream_class(name):
    if name == "DectrisSimplonStreamV2":
        from dxtbx.format.StreamDectrisSimplonStreamV2 import (
            StreamDectrisSimplonStreamV2,
        )

        return StreamDectrisSimplonStreamV2
    elif name == "DectrisSimplonStreamV2Emulator":
        from format_class.dectris_simplon_stream_v2 import (
            StreamDectrisSimplonStreamV2Emulator,
        )

        return StreamDectrisSimplonStreamV2Emulator
    elif name == "DectrisSimplonStreamV2HitFiltered":
        from format_class.dectris_simplon_stream_v2 import (
            StreamDectrisSimplonStreamV2HitFiltered,
        )

        return StreamDectrisSimplonStreamV2HitFiltered()
    elif name == "DectrisSimplonStreamV2EmulatorHitFiltered":
        from format_class.dectris_simplon_stream_v2 import (
            StreamDectrisSimplonStreamV2EmulatorHitFiltered,
        )

        return StreamDectrisSimplonStreamV2EmulatorHitFiltered
    elif name == "DectrisSimplonStreamV1":
        print("EXCEPTION: DectrisSimplonStreamV1 IS NOT IMPLEMENTED")
        assert False
    elif name == "DectrisSimplonStreamV1Emulator":
        print("EXCEPTION: DectrisSimplonStreamV1Emulator IS NOT IMPLEMENTED")
        assert False
    elif name == "LCLStreamer":
        from dxtbx.format.StreamLCLStreamer import StreamLCLStreamer

    else:
        print(f"EXCEPTION: {name} IS NOT IMPLEMENTED")
        assert False
