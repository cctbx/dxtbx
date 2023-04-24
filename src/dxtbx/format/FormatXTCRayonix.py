from __future__ import annotations

import sys

try:
    import psana
except ImportError:
    psana = None
except TypeError:
    # Check if SIT_* environment variables are set
    import os

    if os.environ.get("SIT_ROOT"):
        # Variables are present, so must have been another error
        raise
    psana = None

from libtbx.phil import parse
from scitbx.array_family import flex

from dxtbx.format.FormatXTC import FormatXTC, locator_str

from serialtbx.detector import rayonix

rayonix_locator_str = """
  rayonix {
    bin_size = None
      .type = int
      .help = Detector binning mode
  }
"""

rayonix_locator_scope = parse(rayonix_locator_str + locator_str, process_includes=True)


class FormatXTCRayonix(FormatXTC):
    def __init__(self, image_file, **kwargs):
        super().__init__(image_file, locator_scope=rayonix_locator_scope, **kwargs)

        cfgs = self._ds.env().configStore()
        rayonix_cfg = cfgs.get(psana.Rayonix.ConfigV2, psana.Source("Rayonix"))
        assert rayonix_cfg.binning_f() == rayonix_cfg.binning_s()
        bin_size = rayonix_cfg.binning_f()
        if self.params.rayonix.bin_size is not None:
            assert bin_size == self.params.rayonix.bin_size
        self._pixel_size = rayonix.get_rayonix_pixel_size(bin_size)
        self._image_size = rayonix.get_rayonix_detector_dimensions(self._ds.env())

    @staticmethod
    def understand(image_file):
        try:
            params = FormatXTC.params_from_phil(rayonix_locator_scope, image_file)
        except Exception:
            return False
        return any(["rayonix" in src.lower() for src in params.detector_address])

    def get_raw_data(self, index=None):
        if index is None:
            index = 0
        assert len(self.params.detector_address) == 1
        data = rayonix.get_data_from_psana_event(
            self._get_event(index), self.params.detector_address[0]
        )
        return flex.double(data)

    def get_detector(self, index=None):
        return self._detector()

    def _detector(self):
        return self._detector_factory.simple(
            sensor="UNKNOWN",
            distance=100.0,
            beam_centre=(50.0, 50.0),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(self._pixel_size, self._pixel_size),
            image_size=self._image_size,
            trusted_range=(
                rayonix.rayonix_min_trusted_value,
                rayonix.rayonix_max_trusted_value,
            ),
            mask=[],
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatXTCRayonix.understand(arg))
