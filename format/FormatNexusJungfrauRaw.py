from dxtbx.format.FormatNexus import FormatNexus
import h5py
import numpy as np
from dials.array_family import flex
from cctbx import factor_kev_angstrom


class FormatNexusJungfrauRaw(FormatNexus):
    @staticmethod
    def understand(image_file):
        with h5py.File(image_file, "r") as handle:
            try:
                data = handle["/entry/data"]
                if "pedestal" in data and "gains" in data:
                    return True
            except (KeyError, AttributeError):
                pass

        return False

    def _start(self):
        super(FormatNexusJungfrauRaw, self)._start()
        data_handle = self._reader.entries[0].data[0].handle

        slices = self._raw_data._datalists[0]._all_slices

        def slice_correction_array(array):
            return np.array(
                [array[:, slices[i][0], slices[i][1]] for i in range(len(slices))]
            )

        self._gain = slice_correction_array(data_handle["gains"][()])
        self._pedestal = slice_correction_array(data_handle["pedestal"][()])
        if "pedestalRMS" in data_handle:
            self._pedestalRMS = slice_correction_array(data_handle["pedestalRMS"][()])
        else:
            self._pedestalRMS = None

        self._gainRMS = [np.std(self._gain[:, i]) for i in (0, 1, 2)]

        self._current_index = None
        shape = self._gain[:, 0].shape
        self._current_raw = np.zeros(shape)
        self._current_pedestal = np.zeros(shape)
        self._current_pedestalRMS = np.zeros(shape)
        self._current_gain = np.zeros(shape)
        self._current_gainRMS = np.zeros(shape)

    def unpack_arrays(self, index):
        if index == self._current_index:
            return
        self._current_index = index
        raw = np.array([p.as_numpy_array() for p in self._raw_data[index]])
        gain_map = raw >> 14
        self._current_raw = raw & 0x3FFF  # 0011111111111111
        for gain_mode in (0, 1, 3):
            sel = gain_map == gain_mode
            gain_index = gain_mode if gain_mode != 3 else 2
            self._current_pedestal[sel] = self._pedestal[:, gain_index][sel]
            self._current_pedestalRMS[sel] = self._pedestalRMS[:, gain_index][sel]
            self._current_gain[sel] = self._gain[:, gain_index][sel]
            self._current_gainRMS[sel] = self._gainRMS[gain_index]

    def get_raw_data(self, index):
        self.unpack_arrays(index)
        beam = self.get_beam(index)
        keV = factor_kev_angstrom / beam.get_wavelength()
        corrected = (
            (self._current_raw - self._current_pedestal) / self._current_gain / keV
        )
        return tuple(flex.double(corrected[i]) for i in range(corrected.shape[0]))

    def get_pixel_variance(self, index):
        if self._pedestalRMS is None:
            return None
        self.unpack_arrays(index)
        beam = self.get_beam(index)
        keV = factor_kev_angstrom / beam.get_wavelength()
        spectrum_var = (
            self.get_spectrum(index).get_weighted_energy_variance() / 1000 / 1000
        )  # (sqrt(spectrum var) /1000)^2

        gsq = self._current_gain ** 2

        error = (
            (
                self._current_raw
                + self._current_pedestalRMS ** 2
                + (self._current_gainRMS ** 2 / gsq + spectrum_var / keV ** 2)
                * (self._current_raw - self._current_pedestal) ** 2
            )
            / gsq
            / keV ** 2
        )
        return tuple(flex.double(error[i]) for i in range(error.shape[0]))
