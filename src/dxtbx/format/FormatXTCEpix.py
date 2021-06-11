import sys

import numpy as np

from cctbx import factor_kev_angstrom
from libtbx.phil import parse
from scitbx.array_family import flex
from scitbx.matrix import col

from dxtbx.format.FormatXTC import FormatXTC, locator_str
from dxtbx.model import Detector

epix_locator_str = """
"""

epix_locator_scope = parse(epix_locator_str + locator_str, process_includes=True)


class FormatXTCEpix(FormatXTC):
    def __init__(self, image_file, **kwargs):
        super().__init__(image_file, locator_scope=epix_locator_scope, **kwargs)
        self._cached_detector = {}

    @staticmethod
    def understand(image_file):
        try:
            params = FormatXTC.params_from_phil(epix_locator_scope, image_file)
        except Exception:
            return False
        return any(["epix" in src.lower() for src in params.detector_address])

    def get_raw_data(self, index=None):
        if index is None:
            index = 0
        d = FormatXTCEpix.get_detector(self, index)
        evt = self._get_event(index)
        run = self.get_run_from_index(index)
        det = self._get_psana_detector(run)
        data = det.calib(evt)
        data = data.astype(np.float64)
        # the shape of the epix 10k data is (16, 352, 384)

        self._raw_data = []
        for module_count, module in enumerate(d.hierarchy()):
            for asic_count, asic in enumerate(module):
                fdim, sdim = asic.get_image_size()
                sensor_id = asic_count // 2  # There are 2X2 asics per module
                asic_in_sensor_id = asic_count % 2  # this number will be 0 or 1
                asic_data = data[module_count][
                    sensor_id * sdim : (sensor_id + 1) * sdim,
                    asic_in_sensor_id * fdim : (asic_in_sensor_id + 1) * fdim,
                ]  # 8 sensors per module
                self._raw_data.append(flex.double(np.array(asic_data)))
        assert len(d) == len(self._raw_data), (len(d), len(self._raw_data))
        return tuple(self._raw_data)

    def get_detector(self, index=None):
        return FormatXTCEpix._detector(self, index)

    def _detector(self, index=None):
        from PSCalib.SegGeometryStore import sgs

        from xfel.cftbx.detector.cspad_cbf_tbx import basis_from_geo

        run = self.get_run_from_index(index)
        if run.run() in self._cached_detector:
            return self._cached_detector[run.run()]

        if index is None:
            index = 0
        assert len(self.params.detector_address) == 1

        wavelength = self.get_beam(index).get_wavelength()

        det = self._get_psana_detector(run)

        geom = det.pyda.geoaccess(self._get_event(index).run())
        pixel_size = det.pixel_size(self._get_event(index)) / 1000.0  # convert to mm
        d = Detector()
        pg0 = d.hierarchy()
        # first deal with D0
        det_num = 0
        root = geom.get_top_geo()
        root_basis = basis_from_geo(root)
        while len(root.get_list_of_children()) == 1:
            sub = root.get_list_of_children()[0]
            sub_basis = basis_from_geo(sub)
            root = sub
            root_basis = root_basis * sub_basis
        t = root_basis.translation
        root_basis.translation = col((t[0], t[1], -abs(t[2])))

        origin = col((root_basis * col((0, 0, 0, 1)))[0:3])
        fast = col((root_basis * col((1, 0, 0, 1)))[0:3]) - origin
        slow = col((root_basis * col((0, 1, 0, 1)))[0:3]) - origin
        pg0.set_local_frame(fast.elems, slow.elems, origin.elems)
        pg0.set_name("D%d" % (det_num))

        # Modules next
        for module_num in range(len(root.get_list_of_children())):
            pg1 = pg0.add_group()
            module = root.get_list_of_children()[module_num]
            module_basis = basis_from_geo(module)
            origin = col((module_basis * col((0, 0, 0, 1)))[0:3])
            fast = col((module_basis * col((1, 0, 0, 1)))[0:3]) - origin
            slow = col((module_basis * col((0, 1, 0, 1)))[0:3]) - origin
            pg1.set_local_frame(fast.elems, slow.elems, origin.elems)
            pg1.set_name("D%dM%d" % (det_num, module_num))

            # Read the known layout of the Epix 2x2 module
            sg = sgs.Create(segname=module.oname)
            xx, yy = sg.get_seg_xy_maps_um()
            xx = xx / 1000
            yy = yy / 1000

            # Now deal with ASICs
            for asic_num in range(4):
                val = "ARRAY_D0M%dA%d" % (module_num, asic_num)
                p = pg1.add_panel()
                dim_slow = xx.shape[0]
                dim_fast = xx.shape[1]
                sensor_id = asic_num // 2  # There are 2X2 asics per module
                asic_in_sensor_id = asic_num % 2  # this number will be 0 or 1
                id_slow = sensor_id * (dim_slow // 2)
                id_fast = asic_in_sensor_id * (dim_fast // 2)
                origin = col((xx[id_slow][id_fast], yy[id_slow][id_fast], 0))
                fp = col((xx[id_slow][id_fast + 1], yy[id_slow][id_fast + 1], 0))
                sp = col((xx[id_slow + 1][id_fast], yy[id_slow + 1][id_fast], 0))
                fast = (fp - origin).normalize()
                slow = (sp - origin).normalize()
                p.set_local_frame(fast.elems, slow.elems, origin.elems)
                p.set_pixel_size((pixel_size, pixel_size))
                p.set_image_size((dim_fast // 2, dim_slow // 2))
                p.set_trusted_range((-1, 2e6))
                p.set_gain(factor_kev_angstrom / wavelength)
                p.set_name(val)
        self._cached_detector[run.run()] = d
        return d


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        # Bug, should call this part differently for understand method to work
        print(FormatXTCEpix.understand(arg))
