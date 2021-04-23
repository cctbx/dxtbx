import calendar

import h5py

from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
from dxtbx.format.FormatNexusEigerDLS16M import FormatNexusEigerDLS16M
from dxtbx.format.nexus import h5str


def get_start_time(filename):
    with h5py.File(filename, mode="r") as f:
        timestamp = h5str(f["/entry/start_time"][()]).rstrip("Z")
    return get_pilatus_timestamp(timestamp)


class FormatNexusEigerDLS16MI03(FormatNexusEigerDLS16M):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            name = h5str(FormatNexusEigerDLS16M.get_instrument_name(handle))
            return name and name.upper() in {"DIAMOND BEAMLINE I03", "DLS I03"}

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super().__init__(image_file, **kwargs)
        self._image_file = image_file

    def _start(self):
        super()._start()
        nx = 1028  # module pixels x
        ny = 512  # module pixels y
        dx = 12  # module gap size
        dy = 38  # module gap size
        timestamp = get_start_time(self._image_file)
        if timestamp > calendar.timegm(
            (2021, 4, 22, 0, 0, 0)
        ) and timestamp < calendar.timegm((2021, 5, 4, 0, 0, 0)):
            # 2021 run 2
            # Full module
            bad_modules = {(2, 5), (3, 5)}
            for i, j in bad_modules:
                self._detector_model[0].add_mask(
                    (nx + dx) * i, (ny + dy) * j, (nx + dx) * i + nx, (ny + dy) * j + ny
                )
