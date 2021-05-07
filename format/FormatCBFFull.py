"""
Base implementation of fullCBF format - as used with Dectris detectors
amongst others - this will read the header and construct the full model,
but will allow for extension for specific implementations of CBF.
"""


import sys

import numpy as np
import pycbf

from iotbx.detectors.cbf import CBFImage
from scitbx.array_family import flex

from dxtbx.format.Format import bz2, gzip
from dxtbx.format.FormatCBF import FormatCBF
from dxtbx.format.FormatStill import FormatStill
from dxtbx.format.image import cbf_read_buffer


class FormatCBFFull(FormatCBF):
    """An image reading class for full CBF format images i.e. those from
    a variety of cameras which support this format."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        header = FormatCBF.get_cbf_header(image_file)

        if "_diffrn.id" not in header and "_diffrn_source" not in header:
            return False

        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        # It appears Pycbf can not handle unicode filenames (see dials/dials#256)
        self._raw_data = None
        super().__init__(str(image_file), **kwargs)

    def __del__(self):
        self._cbf_handle.__swig_destroy__(self._cbf_handle)

    def _get_cbf_handle(self):
        try:
            return self._cbf_handle
        except AttributeError:
            self._cbf_handle = pycbf.cbf_handle_struct()

            buffer = None
            # Reopen to tell if it's a gzip - this should be cached, so fast
            with FormatCBF.open_file(self._image_file, "rb") as fin:
                # If this was a gzip or bzip file, read as a buffer
                if (gzip and isinstance(fin._cache_object._file, gzip.GzipFile)) or (
                    bz2 and isinstance(fin._cache_object._file, bz2.BZ2File)
                ):
                    buffer = fin.read()

            if buffer:
                cbf_read_buffer(self._cbf_handle, buffer, pycbf.MSG_DIGEST)
            else:
                self._cbf_handle.read_widefile(
                    self._image_file.encode(), pycbf.MSG_DIGEST
                )

            return self._cbf_handle

    def _goniometer(self):
        """Return a working goniometer instance."""

        return self._goniometer_factory.imgCIF_H(self._get_cbf_handle())

    def _detector(self):
        """Return a working detector instance."""

        return self._detector_factory.imgCIF_H(self._get_cbf_handle(), "unknown")

    def _beam(self):
        """Return a working beam instance."""

        return self._beam_factory.imgCIF_H(self._get_cbf_handle())

    def _scan(self):
        """Return a working scan instance."""

        return self._scan_factory.imgCIF_H(self._image_file, self._get_cbf_handle())

    def detectorbase_start(self):
        self.detectorbase = CBFImage(self._image_file)
        self.detectorbase.readHeader()

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array."""
        if self._raw_data is not None:
            return self._raw_data

        self.detectorbase_start()
        try:
            image = self.detectorbase
            image.read()
            self._raw_data = image.get_raw_data()

            return self._raw_data
        except Exception:
            return None


class FormatCBFFullStill(FormatStill, FormatCBFFull):
    """An image reading class for full CBF format images i.e. those from
    a variety of cameras which support this format. Custom derived from
    the FormatStill to handle images without a gonimeter or scan"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        header = FormatCBF.get_cbf_header(image_file)

        if "_diffrn.id" not in header and "_diffrn_source" not in header:
            return False

        # According to ImageCIF, "Data items in the DIFFRN_MEASUREMENT_AXIS
        # category associate axes with goniometers."
        # http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Cdiffrn_measurement_axis.html
        if "diffrn_measurement_axis" in header:
            return False

        # This implementation only supports single panel.
        try:
            cbf_handle = pycbf.cbf_handle_struct()
            cbf_handle.read_widefile(image_file.encode(), pycbf.MSG_DIGEST)
        except Exception as e:
            if "CBFlib Error" in str(e):
                return False

        # check if multiple arrays
        try:
            return cbf_handle.count_elements() == 1
        except Exception as e:
            if "CBFlib Error" in str(e):
                return False

        return True

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array."""
        if self._raw_data is not None:
            return self._raw_data

        # Override parent's get_raw_data, which relies on iotbx, which in turn
        # relies on cbflib_adaptbx, which in turn expects a gonio
        cbf = self._get_cbf_handle()

        cbf.find_category(b"array_structure")
        cbf.find_column(b"encoding_type")
        cbf.select_row(0)
        types = []
        for i in range(cbf.count_rows()):
            types.append(cbf.get_value())
            cbf.next_row()
        assert (
            len(types) == cbf.count_rows() == 1
        )  # multi-tile data read by different class
        dtype = types[0]

        # find the data
        cbf.select_category(0)
        while cbf.category_name().lower() != b"array_data":
            try:
                cbf.next_category()
            except Exception:
                return None
        cbf.select_column(0)
        cbf.select_row(0)

        cbf.find_column(b"data")
        assert cbf.get_typeofvalue().find(b"bnry") > -1

        # handle floats vs ints
        if dtype == b"signed 32-bit integer":
            array_string = cbf.get_integerarray_as_string()
            self._raw_data = flex.int(np.fromstring(array_string, np.int32))
            parameters = cbf.get_integerarrayparameters_wdims_fs()
            slow, mid, fast = (parameters[11], parameters[10], parameters[9])
            assert slow == 1  # sections not supported
            array_size = mid, fast
        elif dtype == b"signed 64-bit real IEEE":
            array_string = cbf.get_realarray_as_string()
            self._raw_data = flex.double(np.fromstring(array_string, np.float))
            parameters = cbf.get_realarrayparameters_wdims_fs()
            slow, mid, fast = (parameters[7], parameters[6], parameters[5])
            assert slow == 1  # sections not supported
            array_size = mid, fast
        else:
            return None  # type not supported

        self._raw_data.reshape(flex.grid(*array_size))
        return self._raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFFull.understand(arg))
