import sys

import h5py
import numpy as np

try:
    from xfel.util.jungfrau import pad_stacked_format
except ImportError:
    pass

from dials.array_family import flex

from dxtbx.format.FormatNexus import FormatNexus

"""
This class can read a Jungfrau image file that includes, in addition to the corrected data,
the gain maps, and the dark (pedestal) averages and RMS, providing the user access
to that data.
This was written to process SWISSFEL image files created using the CCTBX script jf16m_cxigeom2nexus.py

  libtbx.python cctbx_project/xfel/swissfel/jf16m_cxigeom2nexus.py \
  unassembled_file=X \
  raw_file=Y \
  geom_file=Z \
  pedestal_file=P \
  gain_file=G \
  detector_distance=D \
  include_spectra=True \
  output_file=O \
  raw=False \
  beam_file=B

where X,Y are the corrected,raw hdf5 files obtained at SWISSFEL for a particular run, Z is the crystFEL geom file
provided by SWISSFEL for the jungfrau, P is the hdf5 file provided by swissfel containing the pedestal data,  G is
the gainmap hdf5 file provided by SWISSFEL, D is the detector distance in mm, O is an output filename (to be read by
this format class) and B is the hdf5 file provided by SWISSFEL containing the wavelength spectra measurements (optional)
"""


class FormatNexusJungfrauExt(FormatNexus):
    @staticmethod
    def understand(image_file):
        if "xfel" not in sys.modules:
            return False
        with h5py.File(image_file, "r") as handle:
            try:
                data = handle["/entry/data"]
                if "pedestal" in data and "gains" in data and "raw" in data:
                    return True
            except (KeyError, AttributeError):
                pass
        return False

    def _slice_correction_array(self, array):
        """
        :param array: numpy array where the slow axis is image number
        and e.g. array[i] corresponds to the ith nexus data entry
        :return: array sliced according to the nexus geometry specification
        """
        return np.array(
            [
                array[:, self._slices[i][0], self._slices[i][1]]
                for i in range(len(self._slices))
            ]
        )

    def _start(self):
        super()._start()
        data_handle = self._reader.entries[0].data[0].handle

        # nexus slices
        self._slices = self._raw_data._datalists[0]._all_slices

        self._gain = self._slice_correction_array(data_handle["gains"][()])
        self._pedestal = self._slice_correction_array(data_handle["pedestal"][()])
        self._pedestalRMS = self._slice_correction_array(data_handle["pedestalRMS"][()])

        self.dxtbx_detector_shape = self._gain[:, 0].shape

    def get_14bit_component(self, index, as_flex=False):
        """
        Jungfrau raw data are stored in 16-bit numbers, the 2 significant bits
        encode the gain mode (a number 0,1,2,3), and the other bits store the measurement
        (sum of a photon scattering term, and a pedestal term)
        This function returns the latter term (14-bit measurement)
        :param index: int, shot index
        :param as_flex: bool, return as tuple of flex arrays
        """
        data_handle = self._reader.entries[0].data[0].handle
        raw_16bit = data_handle["raw"][index]
        raw_14bit = pad_stacked_format(raw_16bit & 0x3FFF)
        raw_14bit = self._slice_correction_array(np.array([raw_14bit]))[:, 0]
        if as_flex:
            raw_14bit = tuple([flex.double(d.astype(np.float64)) for d in raw_14bit])
        return raw_14bit

    def get_pedestal_rms(
        self, index, as_flex=False, return_gain_modes=False, divide_by_gain=True
    ):
        """
        This method uses the dynamic (per-shot) pixel gain mode to lookup the per-shot pedestal RMS
        This is important for detailed error analysis of pixel measurements.
        :param index: image index (0-num_images)
        :param as_flex: if True, return tuple of flex arrays (like get_raw_data), else return as a numpy array
        :param return_gain_modes: whether to return the gain modes in addition to the pedestalRMS
        :param divide_by_gain: normalize the pedestalRMS such that the units are the same as get_raw_data (keV)
        """
        if self._pedestalRMS is None:
            raise AttributeError("There is not pedestalRMS in the master file ... ")
        data_handle = self._reader.entries[0].data[0].handle
        raw = data_handle["raw"][index]
        shot_gain_modes = pad_stacked_format(raw >> 14)
        shot_gain_modes = self._slice_correction_array(np.array([shot_gain_modes]))
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
