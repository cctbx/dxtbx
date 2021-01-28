import h5py
import numpy as np

from xfel.util.jungfrau import correct_panel

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

    def pad_raw_data(self, raw_img):
        # NOTE:  currently correct-panel is written for single panels.
        #  TODO speed this up,  do in one-go,
        padded = np.vstack(
            [
                correct_panel(raw_img[i * 512 : (i + 1) * 512], divide=False)
                for i in range(32)
            ]
        )
        return padded

    @staticmethod
    def slice_correction_array(array, slices):
        return np.array(
            [array[:, slices[i][0], slices[i][1]] for i in range(len(slices))]
        )

    def _pad_and_slice_raw_images(self, images):
        padded = np.array([self.pad_raw_data(img) for img in images])
        return self.slice_correction_array(padded, self.slices)

    def _start(self):
        super(FormatNexusJungfrauResRaw, self)._start()
        data_handle = self._reader.entries[0].data[0].handle

        self.slices = self._raw_data._datalists[0]._all_slices
        self._gain = self._pad_and_slice_raw_images(data_handle["gains"])
        self._pedestal = self._pad_and_slice_raw_images(data_handle["pedestal"])
        if "pedestalRMS" in data_handle:
            self._pedestalRMS = self._pad_and_slice_raw_images(
                data_handle["pedestalRMS"]
            )
        else:
            self._pedestalRMS = None

        self.dxtbx_detector_shape = self._gain[:, 0].shape

    def get_pedestal_rms(self, index, as_flex=True, return_gain_modes=False):
        if self._pedestalRMS is None:
            raise AttributeError("There is not pedestalRMS in the master file ... ")
        data_handle = self._reader.entries[0].data[0].handle
        raw = data_handle["raw"][index]
        shot_gain_modes = self.pad_raw_data(raw >> 14)
        shot_gain_modes = self.slice_correction_array(
            np.array([shot_gain_modes]), self.slices
        )
        shot_gain_modes = shot_gain_modes[:, 0]
        shot_rms = np.zeros(self.dxtbx_detector_shape, np.float64)
        for mode in 0, 1, 3:
            gain_index = mode if mode != 3 else 2  # jungfrau quirk
            sel = shot_gain_modes == mode  # pixels in the gain mode
            shot_rms[sel] = self._pedestalRMS[:, gain_index][sel]
        if as_flex:
            shot_rms = tuple([flex.double(p) for p in shot_rms])
        return_packet = shot_rms
        if return_gain_modes:
            return_packet = shot_rms, shot_gain_modes
        return return_packet
