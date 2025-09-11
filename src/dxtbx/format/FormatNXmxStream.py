from dials.array_family import flex
from dxtbx.format.FormatNXmx import FormatNXmx


class FormatNXmxStream(FormatNXmx):
    # FormatMultiImage, Format > FormatHDF5 > FormatNexus > FormatNXmx
    @staticmethod
    def understand(image_file):
        return True
        """
        with h5py.File(image_file) as handle:
            return bool(
                [
                    entry
                    for entry in nxmx.find_class(handle, "NXentry")
                    if "definition" in entry
                    and nxmx.h5str(entry["definition"][()]) == "NXmx"
                ]
            )
        """

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the image file."""
        super().__init__(image_file, **kwargs)

    def get_format_class(self):
        return FormatNXmx

    def set_image(self, image):
        self.image = image

    def get_raw_data(self, index=None):
        if self._bit_depth_readout:
            # if 32 bit then it is a signed int, I think if 8, 16 then it is
            # unsigned with the highest two values assigned as masking values
            if self._bit_depth_readout == 32:
                top = 2**31
            else:
                top = 2**self._bit_depth_readout

        if len(self.image.shape) == 2:
            raw_data = flex.double(self.image)
            if self._bit_depth_readout:
                d1d = raw_data.as_1d()
                d1d.set_selected(d1d == top - 1, -1)
                d1d.set_selected(d1d == top - 2, -2)
        else:
            raw_data = tuple([flex.double(panel) for panel in self.image])
            if self._bit_depth_readout:
                for data in raw_data:
                    d1d = data.as_1d()
                    d1d.set_selected(d1d == top - 1, -1)
                    d1d.set_selected(d1d == top - 2, -2)
        return raw_data
