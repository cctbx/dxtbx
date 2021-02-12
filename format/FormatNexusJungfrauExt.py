import h5py
import numpy as np

from xfel.util.jungfrau import pad_stacked_format

from dials.array_family import flex

from dxtbx.format.FormatNexus import FormatNexus


class FormatNexusJungfrauResRaw(FormatNexus):
    @staticmethod
    def understand(image_file):
        with h5py.File(image_file, "r") as handle:
            try:
                data = handle["/entry/data"]
                if "pedestal" in data and "gains" in data and "raw" in data:
                    return True
            except (KeyError, AttributeError):
                pass
        return False

    @staticmethod
    def slice_correction_array(array, slices):
        return np.array(
            [array[:, slices[i][0], slices[i][1]] for i in range(len(slices))]
        )

    def _start(self):
        super(FormatNexusJungfrauResRaw, self)._start()
        data_handle = self._reader.entries[0].data[0].handle

        self.slices = self._raw_data._datalists[0]._all_slices

        self._gain = self.slice_correction_array(data_handle["gains"][()], self.slices)
        self._pedestal = self.slice_correction_array(
            data_handle["pedestal"][()], self.slices
        )
        if "pedestalRMS" in data_handle:
            self._pedestalRMS = self.slice_correction_array(
                data_handle["pedestalRMS"][()], self.slices
            )
        else:
            self._pedestalRMS = None

        self.dxtbx_detector_shape = self._gain[:, 0].shape

    def get_14bit_component(self, index, as_flex=False):
        data_handle = self._reader.entries[0].data[0].handle
        raw_16bit = data_handle["raw"][index]
        raw_14bit = pad_stacked_format(raw_16bit & 0x3FFF)
        raw_14bit = self.slice_correction_array(np.array([raw_14bit]), self.slices)[
            :, 0
        ]
        if as_flex:
            raw_14bit = tuple([flex.double(d.astype(np.float64)) for d in raw_14bit])
        return raw_14bit

    def get_pedestal_rms(
        self, index, as_flex=False, return_gain_modes=False, divide_by_gain=True
    ):
        if self._pedestalRMS is None:
            raise AttributeError("There is not pedestalRMS in the master file ... ")
        data_handle = self._reader.entries[0].data[0].handle
        raw = data_handle["raw"][index]
        shot_gain_modes = pad_stacked_format(raw >> 14)
        shot_gain_modes = self.slice_correction_array(
            np.array([shot_gain_modes]), self.slices
        )
        shot_gain_modes = shot_gain_modes[:, 0]
        shot_rms = np.zeros(self.dxtbx_detector_shape, np.float64)
        for mode in 0, 1, 3:
            gain_index = mode if mode != 3 else 2  # jungfrau quirk
            sel = shot_gain_modes == mode  # pixels in the gain mode
            pedestalRms = self._pedestalRMS[:, gain_index]
            if divide_by_gain:
                gain = self._gain[:, gain_index]
                shot_rms[sel] = pedestalRms[sel] / np.abs(gain[sel])
            else:
                shot_rms[sel] = pedestalRms[sel]
        if as_flex:
            shot_rms = tuple([flex.double(p) for p in shot_rms])
        return_packet = shot_rms
        if return_gain_modes:
            return_packet = shot_rms, shot_gain_modes
        return return_packet
