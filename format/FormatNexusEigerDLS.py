from __future__ import absolute_import, division, print_function

import ast

import h5py

from dxtbx.format.FormatNexus import FormatNexus


def get_count_limit_from_meta(meta_file_name):
    with h5py.File(meta_file_name, "r") as f:

        config = f["/config"][()]
        config_data = ast.literal_eval(config.decode("utf-8"))

    return config_data["countrate_correction_count_cutoff"]


def get_bit_depth_from_meta(meta_file_name):
    with h5py.File(meta_file_name, "r") as f:

        config = f["/config"][()]
        config_data = ast.literal_eval(config.decode("utf-8"))

    return config_data["bit_depth_image"]


def find_meta_filename(master_like):
    meta_filename = None
    f = h5py.File(master_like, "r")

    def _local_visit(name):
        obj = f[name]
        if not hasattr(obj, "keys"):
            return None
        for k in obj.keys():
            kclass = obj.get(k, getlink=True, getclass=True)
            if kclass is h5py._hl.group.ExternalLink:
                kfile = obj.get(k, getlink=True).filename
                if kfile.split(".")[0].endswith("meta"):
                    return kfile

    meta_filename = f.visit(_local_visit)
    return meta_filename


class FormatNexusEigerDLS(FormatNexus):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            name = FormatNexusEigerDLS.get_instrument_name(handle)
            if name is None or name.lower() not in (b"i03", b"i04", b"vmxi"):
                return False

        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super(FormatNexusEigerDLS, self).__init__(image_file, **kwargs)
        self._meta = find_meta_filename(image_file)
        try:
            self._bit_depth_image = get_bit_depth_from_meta(self._meta)
        except Exception:
            self._bit_depth_image = 0

    def get_detector(self, index=None):
        # workaround for https://jira.diamond.ac.uk/browse/I03-365
        # read the count limit from the meta file - if anything goes
        # wrong, do nothing

        detector = self._detector()

        try:
            limit = get_count_limit_from_meta(self._meta)

            assert limit > 0

            for panel in detector:
                trusted = panel.get_trusted_range()
                panel.set_trusted_range((trusted[0], limit))

        except Exception:
            pass

        return detector

    def get_raw_data(self, index):
        data = self._raw_data[index]
        if self._bit_depth_image:
            # if 32 bit then it is a signed int, I think if 8, 16 then it is
            # unsigned with the highest two values assigned as masking values
            if self._bit_depth_image == 32:
                top = 2 ** 31
            else:
                top = 2 ** self._bit_depth_image
            d1d = data.as_1d()
            d1d.set_selected(d1d == top - 1, -1)
            d1d.set_selected(d1d == top - 2, -2)

        return data
