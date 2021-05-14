import sys

import numpy as np

from cctbx.eltbx import attenuation_coefficient
from libtbx.phil import parse
from scitbx.array_family import flex
from scitbx.matrix import col
from xfel.cftbx.detector.cspad_cbf_tbx import read_slac_metrology
from xfel.cxi.cspad_ana.cspad_tbx import env_distance

from dxtbx.format.FormatXTC import FormatXTC, locator_str
from dxtbx.model import Detector, ParallaxCorrectedPxMmStrategy

try:
    from xfel.cftbx.detector import cspad_cbf_tbx
    from xfel.cxi.cspad_ana import cspad_tbx
except ImportError:
    # xfel not configured
    pass

cspad_locator_str = """
  cspad {
    detz_offset = None
      .type = float
      .help = Distance from back of detector rail to sample interaction region (CXI) \
              or actual detector distance (XPP/MFX)
    apply_gain_mask = True
      .type = bool
      .help = flag to indicate if gain should be applied to cspad data
    dark_correction = True
      .type = bool
      .help = flag to decide if dark correction should be done
    use_psana_calib = False
      .type = bool
      .help = Use the psana calibration
    common_mode = default
      .type = str
      .help = Common mode correction "default", "cspad_default", or "unbonded"\
              see https://confluence.slac.stanford.edu/display/PSDM/Common+mode+correction+algorithms\
              default means no common mode corrections... the other two are psana corrections
    }
"""

cspad_locator_scope = parse(cspad_locator_str + locator_str, process_includes=True)


class FormatXTCCspad(FormatXTC):
    def __init__(self, image_file, locator_scope=cspad_locator_scope, **kwargs):
        super().__init__(image_file, locator_scope=locator_scope, **kwargs)
        assert (
            self.params.cspad.detz_offset is not None
        ), "Supply a detz_offset for the cspad"
        self._cache_psana_pedestals()  # NOTE: move to base FormatXTC class
        self._psana_gain_map_cache = {}

    @staticmethod
    def understand(image_file):
        try:
            params = FormatXTC.params_from_phil(cspad_locator_scope, image_file)
        except Exception:
            return False
        return any("cspad" in src.lower() for src in params.detector_address)

    def _get_psana_gain_map(self, run):
        """
        checks if user wants gain applied and caches a gain map per run
        """
        if run.run() not in self._psana_gain_map_cache:
            if self.params.cspad.apply_gain_mask:
                self._psana_gain_map_cache[run.run()] = (
                    self._get_psana_detector(run).gain_mask(run) > 0
                )
            else:
                self._psana_gain_map_cache[run.run()] = None

    def _cache_psana_pedestals(self):
        """Store a pedestal for each psana detector instance"""
        self._pedestals = {}
        for run_number, run in self._psana_runs.items():
            det = self._get_psana_detector(run)
            self._pedestals[run_number] = det.pedestals(run)

    def get_raw_data(self, index=None):
        if index is None:
            index = 0
        assert len(self.params.detector_address) == 1
        d = self.get_detector(index)
        event = self._get_event(index)
        run_number = event.run()
        run = self._psana_runs[run_number]
        det = self._get_psana_detector(run)
        data = cspad_cbf_tbx.get_psana_corrected_data(
            det,
            event,
            use_default=self.params.cspad.use_psana_calib,
            dark=self._pedestals[run_number],
            common_mode=self.params.cspad.common_mode,
            apply_gain_mask=self.params.cspad.apply_gain_mask,
            gain_mask_value=None,
            per_pixel_gain=False,
            gain_mask=self._get_psana_gain_map(run),
        )
        data = data.astype(np.float64)
        self._raw_data = []
        for quad_count, quad in enumerate(d.hierarchy()):
            for sensor_count, sensor in enumerate(quad):
                for asic_count, asic in enumerate(sensor):
                    fdim, sdim = asic.get_image_size()
                    asic_data = data[
                        sensor_count + quad_count * 8,
                        :,
                        asic_count * fdim : (asic_count + 1) * fdim,
                    ]  # 8 sensors per quad
                    self._raw_data.append(flex.double(np.array(asic_data)))
        assert len(d) == len(self._raw_data)
        return tuple(self._raw_data)

    def get_detector(self, index=None):
        return FormatXTCCspad._detector(self, index)

    # XXX Implement recursive version
    def _detector(self, index=None):
        if index is None:
            index = 0

        run = self.get_run_from_index(index)
        det = self._get_psana_detector(run)
        geom = det.pyda.geoaccess(run.run())
        cob = read_slac_metrology(geometry=geom, include_asic_offset=True)
        distance = env_distance(
            self.params.detector_address[0], run.env(), self.params.cspad.detz_offset
        )
        d = Detector()
        pg0 = d.hierarchy()
        # first deal with D0
        det_num = 0
        origin = col((cob[(0,)] * col((0, 0, 0, 1)))[0:3])
        fast = col((cob[(0,)] * col((1, 0, 0, 1)))[0:3]) - origin
        slow = col((cob[(0,)] * col((0, 1, 0, 1)))[0:3]) - origin
        origin += col((0.0, 0.0, -distance))
        pg0.set_local_frame(fast.elems, slow.elems, origin.elems)
        pg0.set_name("D%d" % (det_num))
        for quad_num in range(4):
            # Now deal with Qx
            pg1 = pg0.add_group()
            origin = col((cob[(0, quad_num)] * col((0, 0, 0, 1)))[0:3])
            fast = col((cob[(0, quad_num)] * col((1, 0, 0, 1)))[0:3]) - origin
            slow = col((cob[(0, quad_num)] * col((0, 1, 0, 1)))[0:3]) - origin
            pg1.set_local_frame(fast.elems, slow.elems, origin.elems)
            pg1.set_name("D%dQ%d" % (det_num, quad_num))
            for sensor_num in range(8):
                # Now deal with Sy
                pg2 = pg1.add_group()
                origin = col((cob[(0, quad_num, sensor_num)] * col((0, 0, 0, 1)))[0:3])
                fast = (
                    col((cob[(0, quad_num, sensor_num)] * col((1, 0, 0, 1)))[0:3])
                    - origin
                )
                slow = (
                    col((cob[(0, quad_num, sensor_num)] * col((0, 1, 0, 1)))[0:3])
                    - origin
                )
                pg2.set_local_frame(fast.elems, slow.elems, origin.elems)
                pg2.set_name("D%dQ%dS%d" % (det_num, quad_num, sensor_num))
                # Now deal with Az
                for asic_num in range(2):
                    val = "ARRAY_D0Q%dS%dA%d" % (quad_num, sensor_num, asic_num)
                    p = pg2.add_panel()
                    origin = col(
                        (cob[(0, quad_num, sensor_num, asic_num)] * col((0, 0, 0, 1)))[
                            0:3
                        ]
                    )
                    fast = (
                        col(
                            (
                                cob[(0, quad_num, sensor_num, asic_num)]
                                * col((1, 0, 0, 1))
                            )[0:3]
                        )
                        - origin
                    )
                    slow = (
                        col(
                            (
                                cob[(0, quad_num, sensor_num, asic_num)]
                                * col((0, 1, 0, 1))
                            )[0:3]
                        )
                        - origin
                    )
                    p.set_local_frame(fast.elems, slow.elems, origin.elems)
                    p.set_pixel_size(
                        (cspad_cbf_tbx.pixel_size, cspad_cbf_tbx.pixel_size)
                    )
                    p.set_image_size(cspad_cbf_tbx.asic_dimension)
                    p.set_trusted_range(
                        (
                            cspad_tbx.cspad_min_trusted_value,
                            cspad_tbx.cspad_saturated_value,
                        )
                    )
                    p.set_name(val)

        try:
            beam = self._beam(index)
        except Exception:
            print(
                "No beam object initialized. Returning CSPAD detector without parallax corrections"
            )
            return d

        # take into consideration here the thickness of the sensor also the
        # wavelength of the radiation (which we have in the same file...)
        wavelength = beam.get_wavelength()
        thickness = 0.5  # mm, see Hart et al. 2012

        table = attenuation_coefficient.get_table("Si")
        # mu_at_angstrom returns cm^-1
        mu = table.mu_at_angstrom(wavelength) / 10.0  # mu: mm^-1
        t0 = thickness
        for panel in d:
            panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
        return d


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatXTCCspad.understand(arg))
