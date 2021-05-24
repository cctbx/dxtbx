import h5py

from dxtbx.format.FormatNexus import FormatNexus


class FormatNexusEiger(FormatNexus):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            try:
                det = handle["/entry/instrument/detector/description"][()]
                if det.lower().startswith(b"eiger"):
                    return True
            except (KeyError, AttributeError):
                pass

        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super().__init__(image_file, **kwargs)

    def get_detector(self, index=None):
        # override base detector to assign correct value for minimum
        # trusted if it's -0x7ffffff

        detector = self._detector()

        for panel in detector:
            trusted = panel.get_trusted_range()
            if trusted[0] == -0x7FFFFFFF:
                panel.set_trusted_range((-1, trusted[1]))

        return detector
