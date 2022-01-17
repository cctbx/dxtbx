try:
    from importlib import metadata

    version = metadata.version("dxtbx")
except ImportError:
    try:
        # Running on pre-3.8 Python; use importlib-metadata backport package
        import importlib_metadata as metadata  # type: ignore

        version = metadata.version("dxtbx")
    except ImportError:
        # If we don't have importlib or the backport assume 3.8 series.
        # In this scenario we probably aren't going to be checking this anyway.
        version = "3.8.dev"
