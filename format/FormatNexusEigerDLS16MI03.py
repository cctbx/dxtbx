import calendar

import h5py

from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
from dxtbx.format.FormatNexusEigerDLS16M import FormatNexusEigerDLS16M
from dxtbx.format.nexus import h5str


class FormatNexusEigerDLS16MI03(FormatNexusEigerDLS16M):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as fh:
            name = h5str(FormatNexusEigerDLS16M.get_instrument_name(fh))
            if name and name.upper() in {"DIAMOND BEAMLINE I03", "DLS I03"}:
                timestamp = get_pilatus_timestamp(
                    h5str(fh["/entry/start_time"][()]).rstrip("Z")
                )
                return timestamp > calendar.timegm(
                    (2021, 4, 22, 0, 0, 0)
                ) and timestamp < calendar.timegm((2021, 5, 4, 0, 0, 0))

    def _start(self):
        super()._start()
        nx = 1028  # module pixels x
        ny = 512  # module pixels y
        dx = 12  # module gap size
        dy = 38  # module gap size
        if self._detector_model[0].get_image_size() == (4148, 4362):
            # 2021 run 2
            # Full module
            bad_modules = {(2, 5), (3, 5)}
            for i_col, i_row in bad_modules:
                self._detector_model[0].add_mask(
                    (nx + dx) * i_col,
                    (ny + dy) * i_row,
                    (nx + dx) * i_col + nx,
                    (ny + dy) * i_row + ny,
                )
