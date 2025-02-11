from __future__ import annotations

import weakref

import h5py
import nxmx

import scitbx.array_family.flex as flex

import dxtbx.nexus
from dxtbx.format.FormatNexus import FormatNexus


class _MaskCache:
    """A singleton to hold unique static_mask objects to avoid duplications"""

    def __init__(self):
        self.local_mask_cache = weakref.WeakValueDictionary()

    def _mask_hasher(self, mask: flex.bool) -> int:
        return hash(mask.as_numpy_array().tobytes())

    def store_unique_and_get(
        self, mask_tuple: tuple[flex.bool, ...] | None
    ) -> tuple[flex.bool, ...] | None:
        if mask_tuple is None:
            return None
        output = []
        for mask in mask_tuple:
            mask_hash = self._mask_hasher(mask)
            mask = self.local_mask_cache.setdefault(mask_hash, mask)
            output.append(mask)
        return tuple(output)


mask_cache = _MaskCache()


def detector_between_sample_and_source(detector, beam):
    """Check if the detector is perpendicular to beam and
    upstream of the sample."""

    if len(detector) != 1:
        return False
    p = detector[0]

    n = p.get_normal()
    o = p.get_origin()
    x = beam.get_sample_to_source_direction()

    dot = sum(_n * _x for _n, _x in zip(n, x))
    doto = sum(_o * _x for _o, _x in zip(o, x))

    if abs(dot) > 0.99 and doto > 0:
        return True

    return False


def inverted_distance_detector(detector):
    """Move the detector origin to x, y, -z"""

    origin = detector[0].get_origin()
    origin = origin[0], origin[1], -origin[2]
    fast = detector[0].get_fast_axis()
    slow = detector[0].get_slow_axis()
    detector[0].set_frame(fast, slow, origin)
    return detector


class FormatNXmx(FormatNexus):
    """Read NXmx-flavour NeXus-format HDF5 data."""

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            return (
                bool(
                    [
                        entry
                        for entry in nxmx.find_class(handle, "NXentry")
                        if "definition" in entry
                        and nxmx.h5str(entry["definition"][()]) == "NXmx"
                    ]
                )
                and "/entry/data/data" in handle
            )

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        super().__init__(image_file, **kwargs)

    def _start(self):
        self._static_mask = None

        self._cached_file_handle = h5py.File(self._image_file, swmr=True)
        nxmx_obj = self._get_nxmx(self._cached_file_handle)
        nxentry = nxmx_obj.entries[0]
        nxsample = nxentry.samples[0]
        nxinstrument = nxentry.instruments[0]
        nxdetector = nxinstrument.detectors[0]
        nxbeam = nxinstrument.beams[0]
        nxdata = nxmx_obj.entries[0].data[0]
        self._goniometer_model = dxtbx.nexus.get_dxtbx_goniometer(nxsample)
        self._beam_factory = dxtbx.nexus.CachedWavelengthBeamFactory(nxbeam)
        wavelength = self._beam_factory.make_beam(index=0).get_wavelength()
        self._detector_model = dxtbx.nexus.get_dxtbx_detector(
            nxdetector, wavelength, nxdata
        )

        # if the detector is between the sample and the source, and perpendicular
        # to the beam, then invert the distance vector, as this is probably wrong
        beam = self._beam()
        if detector_between_sample_and_source(self._detector_model, beam):
            self._detector_model = inverted_distance_detector(self._detector_model)

        self._scan_model = dxtbx.nexus.get_dxtbx_scan(nxsample, nxdetector)
        self._static_mask = mask_cache.store_unique_and_get(
            dxtbx.nexus.get_static_mask(nxdetector)
        )
        self._bit_depth_readout = nxdetector.bit_depth_readout

        if self._scan_model:
            self._num_images = len(self._scan_model)
        else:
            if nxdata.signal:
                data = nxdata[nxdata.signal]
            else:
                data = list(nxdata.values())[0]
            self._num_images, *_ = data.shape

        # handle reversed rotation by inverting both the rotation angles
        # and the rotation axis
        if self._scan_model:
            scan = self._scan_model
            oscillation = scan.get_oscillation()

            if oscillation[1] < 0:
                # invert the scan parameters so that the attitude of roation is positive
                self._scan_model.set_oscillation((-oscillation[0], -oscillation[1]))

                # invert the axis direction so that the universe is aligned
                axes = self._goniometer_model.get_axes()
                axis = self._goniometer_model.get_scan_axis()
                axes[axis] = -1 * axes[axis][0], -1 * axes[axis][1], -1 * axes[axis][2]
                self._goniometer_model.set_axes(axes)

    def _get_nxmx(self, fh: h5py.File):
        return nxmx.NXmx(fh)

    def _beam(self, index: int | None = None) -> dxtbx.model.Beam:
        return self._beam_factory.make_beam(index=index or 0)

    def get_spectrum(self, index=None):
        return self._beam_factory.make_spectrum(index=index or 0)

    def get_num_images(self) -> int:
        return self._num_images

    def get_static_mask(self, index=None, goniometer=None):
        return self._static_mask

    def get_raw_data(self, index):
        nxmx_obj = self._get_nxmx(self._cached_file_handle)
        nxdata = nxmx_obj.entries[0].data[0]
        nxdetector = nxmx_obj.entries[0].instruments[0].detectors[0]
        raw_data = dxtbx.nexus.get_raw_data(
            nxdata, nxdetector, index, bit_depth=self._bit_depth_readout
        )
        if self._bit_depth_readout:
            # if 32 bit then it is a signed int, I think if 8, 16 then it is
            # unsigned with the highest two values assigned as masking values
            if self._bit_depth_readout == 32:
                top = 2**31
            else:
                top = 2**self._bit_depth_readout
            for data in raw_data:
                d1d = data.as_1d()
                d1d.set_selected(d1d == top - 1, -1)
                d1d.set_selected(d1d == top - 2, -2)
        return raw_data
