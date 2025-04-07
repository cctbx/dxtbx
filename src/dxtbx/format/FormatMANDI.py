from __future__ import annotations

import sys
from multiprocessing import Pool, cpu_count
from os.path import join

import h5py
import numpy as np

from cctbx.array_family import flex

import dxtbx_flumpy as flumpy
from dxtbx import IncorrectFormatError
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.model import Detector
from dxtbx.model.beam import BeamFactory, PolychromaticBeam, Probe
from dxtbx.model.goniometer import Goniometer, GoniometerFactory
from dxtbx.model.scan import Scan, ScanFactory


class FormatMANDI(FormatHDF5):
    """
    Class to read NXTOFRAW from MaNDi
    (https://neutrons.ornl.gov/mandi)
    """

    def __init__(self, image_file: str, **kwargs) -> None:
        super().__init__(image_file)
        if not FormatMANDI.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self.image_file = image_file
        self.nxs_file = self.open_file(image_file)
        self._base_entry = self.get_base_entry_name(self.nxs_file)
        self.detector = None
        self.raw_data = None

    def open_file(self, image_file_path: str) -> h5py.File:
        return h5py.File(image_file_path, "r")

    def get_base_entry_name(self, nxs_file: h5py.File) -> str:
        return list(nxs_file.keys())[0]

    @staticmethod
    def understand(image_file: str) -> bool:
        try:
            return FormatMANDI.is_mandi_file(image_file)
        except OSError:
            return False

    @staticmethod
    def is_mandi_file(image_file: str) -> bool:
        def get_name(image_file: str) -> str:
            with h5py.File(image_file, "r") as handle:
                if len(handle) == 0:
                    return ""
                base_entry = list(handle.keys())[0]
                if f"{base_entry}/instrument/name" not in handle:
                    return ""
                try:
                    return handle[f"{base_entry}/instrument/name"][0].decode()
                except (ValueError, IndexError):
                    return ""

        return get_name(image_file) == "MANDI"

    def get_raw_data(self, index: int) -> tuple[flex.int]:
        raw_data = []
        panel_size = self._get_image_size()
        for panel_name in self._get_panel_names():
            spectra = np.reshape(
                self.nxs_file[self._base_entry][f"{panel_name}_events"]["spectra"][
                    :, index : index + 1
                ],
                panel_size,
            ).T
            # spectra = self.nxs_file[self._base_entry][f"{panel_name}_events"]["spectra"]
            raw_data.append(flumpy.from_numpy(np.ascontiguousarray(spectra)))

        return tuple(raw_data)

    def get_detector(self) -> Detector:
        num_panels = self._get_num_panels()
        panel_names = self._get_panel_names()
        panel_type = self._get_panel_type()
        image_size = self._get_image_size()
        trusted_range = self._get_panel_trusted_range()
        pixel_size = self._get_pixel_size()
        fast_axes = self._get_panel_fast_axes()
        slow_axes = self._get_panel_slow_axes()
        panel_origins = self._get_panel_origins()
        gain = self._get_panel_gain()
        panel_projections = self._get_panel_projections_2d()
        detector = Detector()
        root = detector.hierarchy()

        for i in range(num_panels):
            panel = root.add_panel()
            panel.set_type(panel_type)
            panel.set_name(panel_names[i])
            panel.set_image_size(image_size)
            panel.set_trusted_range(trusted_range)
            panel.set_pixel_size(pixel_size)
            panel.set_local_frame(fast_axes[i], slow_axes[i], panel_origins[i])
            panel.set_gain(gain)
            r, t = panel_projections[i + 1]
            r = tuple(map(int, r))
            t = tuple(map(int, t))
            panel.set_projection_2d(r, t)

        return detector

    def _get_time_channel_bins(self) -> list[float]:
        # (usec)
        return self.nxs_file[self._base_entry]["time_of_flight"][:]

    def _get_time_of_flight(self) -> list[float]:
        # (usec)
        bins = self._get_time_channel_bins()
        return tuple([float(bins[i] + bins[i + 1]) * 0.5 for i in range(len(bins) - 1)])

    def _get_sample_to_source_distance(self) -> float:
        # (mm)
        return 30000

    def _get_sample_to_source_direction(self) -> tuple[float, float, float]:
        return (0, 0, -1)

    def _get_num_panels(self) -> int:
        return 40

    def _get_panel_names(self) -> tuple[str]:
        return (
            "bank1",
            "bank2",
            "bank3",
            "bank5",
            "bank7",
            "bank8",
            "bank11",
            "bank12",
            "bank13",
            "bank17",
            "bank18",
            "bank19",
            "bank20",
            "bank21",
            "bank22",
            "bank23",
            "bank26",
            "bank27",
            "bank28",
            "bank29",
            "bank31",
            "bank32",
            "bank33",
            "bank37",
            "bank39",
            "bank40",
            "bank41",
            "bank42",
            "bank43",
            "bank46",
            "bank47",
            "bank48",
            "bank49",
            "bank50",
            "bank51",
            "bank52",
            "bank53",
            "bank57",
            "bank58",
            "bank59",
        )

    def _get_panel_gain(self) -> float:
        return 1.0

    def _get_panel_trusted_range(self) -> tuple[int, int]:
        return (-1, 100000)

    def _get_panel_origins(self) -> tuple[tuple[float, float, float], ...]:
        return (
            (71.59649688799595, -407.9875271017446, 91.63277006330976),
            (247.78726107849047, -370.86794223009696, 93.05022800601382),
            (87.97790138773688, -409.33375923673844, -82.03599884101118),
            (-85.64058504589039, -403.12747733399465, -97.65220987551567),
            (-102.35072309461592, -404.43783650393084, 75.318849784288),
            (-101.376307910567, -362.7034047105274, 249.52766339419702),
            (188.46366913374928, -340.8394179451527, 264.04352540644584),
            (304.0861943336393, -339.0817491399483, 97.69270882922524),
            (306.2755999937922, -340.057049383876, -104.79738499213931),
            (-317.25124354686346, -327.83699783860874, -104.13714348985233),
            (-317.64898298095443, -328.5375200345585, 98.10329861513794),
            (-199.4658588654215, -332.9638764406945, 263.6279591286127),
            (114.10208135428759, -227.42718157198428, 359.56893399219473),
            (303.4066657317839, -229.171259967358, 219.83460648731807),
            (373.7935811124756, -229.42481664792578, -3.4456330231885866),
            (307.8384913706013, -230.97676263539856, -228.302717512744),
            (-312.1677269773796, -218.49351412060085, -227.8285640453776),
            (-380.89278569326314, -215.05244003524544, -2.9178356187804786),
            (-309.1414237681691, -217.46267271820048, 219.27597207117535),
            (-121.94996995674012, -223.24489754625657, 361.1688851013957),
            (233.5580884765005, -115.43121703586225, 320.63922658665984),
            (378.45031177159586, -117.16549490064116, 120.41696423222004),
            (372.89800332945197, -115.93989744162306, -124.92982019508003),
            (-376.5169852579795, -102.17355196399106, -124.26626333114893),
            (-234.11393874355227, -106.04709661250322, 316.95386906093717),
            (137.13740529014822, 9.162292313720426, 417.1619791626718),
            (357.2696434729768, 4.7468832636133955, 256.1368430867005),
            (438.2478931969796, 3.3972037803596535, -3.291783715503191),
            (361.08469499737765, 6.145340038510798, -263.82143408424554),
            (-357.74327638098214, 18.735783076242843, -263.98028087869886),
            (-439.2712284275042, 20.01853918185089, -2.3600188451399426),
            (-355.6934447072983, 17.71984081232332, 255.6674822459931),
            (-136.15700784634026, 13.793781420194758, 419.92464356648065),
            (3.0211639276166506, 147.5028475462614, 442.29307467579065),
            (264.40926501753853, 142.959186239358, 357.0786864461248),
            (427.5813705436181, 141.24925820955698, 134.73307743251732),
            (428.2623879171177, 143.0955989020943, -141.14258845564459),
            (-418.3256035559564, 155.87101472044327, -139.2255203535932),
            (-421.1137777088633, 156.31099046243813, 134.80694824067064),
            (-259.44611581654755, 152.81413474372576, 357.96944889204525),
        )

    def _get_panel_fast_axes(self) -> tuple[tuple[float, float, float], ...]:
        return (
            (-0.9999583639962794, 0.004549995022478471, 0.007909982249012876),
            (-0.9270112280582247, -0.36646049406596276, 0.07973010311615258),
            (9.755968675898316e-17, 5.59387233726922e-17, 1.0),
            (0.9999739794108949, -0.005200025869334086, 0.0050000232096750075),
            (-0.0007900093059066182, -3.0005133535454756e-05, -0.9999996874924455),
            (-0.09582968262311213, -0.39675869523887103, -0.9129070104236807),
            (-0.3387401541674543, 0.6009102631780291, 0.7239903062623982),
            (0.163969987359823, 0.5861999502774039, 0.7933999379505827),
            (0.570700354137915, 0.6004903814689775, 0.5601003548919607),
            (-0.13169977290772172, 0.6018189486770399, -0.7876986243667796),
            (-0.5578201777633975, 0.6146602050286276, -0.5577001718075555),
            (-0.7780136583020894, 0.6167928979211943, -0.1194205533791195),
            (-0.6024392271286111, 0.6822191251868116, 0.41429946034966436),
            (-0.26297968340353617, 0.6894191631846113, 0.6749391850758157),
            (0.17966967519294516, 0.6836287651503686, 0.7073687293590546),
            (0.5733071659167379, 0.6797166377204127, 0.45749774416337297),
            (0.25794006072273956, 0.6797101771949186, -0.6866301770909884),
            (-0.18668894896436913, 0.6831461280707687, -0.7060160083429384),
            (-0.561269909625896, 0.6879698755373065, -0.4600799266450577),
            (-0.7216806782398583, 0.6916006515467016, -0.02942001759747731),
            (-0.5406021594612158, 0.7160628615705423, 0.44159176109306714),
            (-0.22815117249667427, 0.7086036629348159, 0.667703445673118),
            (0.24252895759587326, 0.7016769746800912, 0.669947108308748),
            (0.23261954038002947, 0.7028586076557112, -0.6722186601677009),
            (-0.5597197589436586, 0.717469677505232, -0.4146698123912455),
            (-0.7199483374110611, 0.6930083967891861, 0.03759991274897956),
            (-0.5756319823045508, 0.6837223456610776, 0.44852154350916884),
            (-0.18817009429166837, 0.6803903231242386, 0.7082803285515974),
            (0.25794006072273956, 0.6797101771949186, 0.6866301770909884),
            (0.5733071659167379, 0.6797166377204127, -0.45749774416337297),
            (0.20908971333938364, 0.6753790691234415, -0.7072090248049779),
            (-0.24702047603368216, 0.6870212856103326, -0.6833612789282247),
            (-0.5932592018008807, 0.6785590800586375, -0.43312941986034337),
            (-0.7000908028386412, 0.5988306911020218, -0.38894044682318896),
            (-0.7843286772076716, 0.6097189697574565, 0.11432980375899326),
            (-0.5632010367672651, 0.6076211250557385, 0.5600010362225023),
            (-0.14986972872339283, 0.5870889266408779, 0.7955285391662817),
            (0.5752912394628964, 0.5971712830921335, -0.5589512039948916),
            (0.1400298828594605, 0.605939509647955, -0.7830893579623988),
            (-0.3374699569875644, 0.606049917647629, -0.720289889870817),
        )

    def _get_panel_slow_axes(self) -> tuple[tuple[float, float, float], ...]:
        return (
            (0.009080058381285394, 0.4099593618136086, 0.9120585914299427),
            (-0.07844604621503816, 0.39736287118342395, 0.914302448010555),
            (0.9135486805116207, 0.4067294043162809, -1.1187744674538443e-16),
            (0.006681901394788637, 0.40642071283978043, -0.9136616202777456),
            (-0.9000874338744042, 0.4357087607239177, 0.0006980041674738889),
            (-0.9107931255862123, 0.40496046981334777, -0.08039216549795625),
            (0.7862610251238533, 0.6034193349109014, -0.13296129747739685),
            (0.5893302132278937, 0.5867709742995678, -0.555328482518501),
            (0.14371018354248974, 0.598521065045708, -0.7881116150918823),
            (-0.5542598612993068, 0.6140978131535265, 0.5618539686008542),
            (-0.14408786300718987, 0.590035000756182, 0.7944163805062643),
            (0.35367928103827984, 0.5871085430589494, 0.7281583102811113),
            (0.7350692860342056, 0.6764910100824224, -0.04508944451672663),
            (0.5937905665933713, 0.6670312021437351, -0.44998015333057373),
            (0.23556951108705004, 0.6682489972510377, -0.7056559240303951),
            (-0.2579398836750756, 0.6797104022878586, -0.6866300207763755),
            (-0.5733092269127125, 0.6797160990021477, 0.4574959618335998),
            (-0.18693934596038786, 0.680817634844769, 0.7081956149371988),
            (0.24370146751944194, 0.6686329796911721, 0.7025236887096229),
            (0.6312723758233735, 0.6749771526277607, 0.3819699346189881),
            (0.5598433004196994, 0.6980113089477805, -0.446492655657609),
            (0.2109292083765029, 0.7054769375484777, -0.6766174396517453),
            (-0.21138055008498535, 0.712193415198309, -0.669401824313273),
            (-0.20838062049588585, 0.7111470099447131, 0.671451746031272),
            (0.5778286402826893, 0.696587990671485, 0.4252989933227231),
            (0.6222313445360602, 0.6685195040380297, -0.4073202997612487),
            (0.2748273430210124, 0.678368281303838, -0.6813856517779691),
            (-0.15815963132467423, 0.6907582508433235, -0.705576764010207),
            (-0.5733092269127125, 0.6797160990021477, -0.4574959618335998),
            (-0.2579398836750756, 0.6797104022878586, 0.6866300207763755),
            (0.20974959544136526, 0.6753887094166133, 0.7070043835824111),
            (0.5889124037382343, 0.6664719408308167, 0.45716226091887696),
            (0.7351718509309583, 0.6758902572003087, 0.0519105944913895),
            (0.7138073326114943, 0.6012290545520337, -0.3591694807086639),
            (0.3550890340368667, 0.5923840954415561, -0.7231824537242688),
            (-0.13366751242833205, 0.6017868942416856, -0.7873915989138873),
            (-0.5813390730622492, 0.598508802800182, -0.5512096652833816),
            (0.14377848861475973, 0.5988879167456186, 0.7878204169655343),
            (0.5985395808907364, 0.5782215084773593, 0.5544459010952365),
            (0.8005222527712781, 0.587343547301513, 0.11912883891515713),
        )

    def get_num_images(self) -> int:
        return len(self._get_time_of_flight())

    def get_beam(self, index: int = None) -> PolychromaticBeam:
        direction = self._get_sample_to_source_direction()
        distance = self._get_sample_to_source_distance()
        return BeamFactory.make_polychromatic_beam(
            direction=direction, sample_to_source_distance=distance, probe=Probe.neutron
        )

    def get_scan(self, index=None) -> Scan:
        image_range = (1, self.get_num_images())
        properties = {"time_of_flight": self._get_time_of_flight()}
        return ScanFactory.make_scan_from_properties(
            image_range=image_range, properties=properties
        )

    def get_goniometer(self, idx: int = None) -> Goniometer:
        rotation_axis_phi = (0.0, 1.0, 0.0)
        rotation_axis_omega = (0.0, 1.0, 0.0)
        rotation_axis_chi = (0.0, 0.0, 1.0)
        fixed_rotation = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

        goniometer = GoniometerFactory.make_goniometer(
            rotation_axis_phi, fixed_rotation
        )
        phi = self.nxs_file["entry/DASlogs/phi/average_value"][0]
        omega = self.nxs_file["entry/DASlogs/omega/average_value"][0]
        chi = self.nxs_file["entry/DASlogs/chi/average_value"][0]

        goniometer.rotate_around_origin(rotation_axis_phi, phi)
        goniometer.rotate_around_origin(rotation_axis_chi, chi)
        goniometer.rotate_around_origin(rotation_axis_omega, omega)
        return goniometer

    def _get_image_size(self) -> tuple[int, int]:
        return (256, 256)

    def _get_pixel_size(self) -> tuple[float, float]:
        return (0.618, 0.618)

    def _get_panel_type(self) -> str:
        return "SENSOR_PAD"

    def _get_panel_projections_2d(self) -> dict:
        p_w, p_h = self._get_pixel_size()
        panel_pos = {}
        count = 1
        for i in range(8):
            for j in range(5):
                panel_pos[count] = ((1, 0, 0, 1), (p_h * i, p_w * j))
                count += 1
        return panel_pos

    @staticmethod
    def add_histogram_data_to_nxs_file(
        nxs_file_path: str,
        remove_event_data: bool,
        spectra_output_name: str = "spectra",
        write_tof_bins: bool = True,
        delta_tof: float = 5,  # (usec)
        tof_padding: float = 100,  # (usec)
        panel_size: tuple[int, int] = (256, 256),  # (px)
        nproc: int = 8,
    ) -> None:
        tof_bins = FormatMANDI.generate_tof_bins(
            nxs_file=nxs_file_path,
            panel_size=panel_size,
            delta_tof=delta_tof,
            padding=tof_padding,
        )
        FormatMANDI.write_histogram_data(
            nxs_file_path=nxs_file_path,
            tof_bins=tof_bins,
            panel_size=panel_size,
            remove_event_data=remove_event_data,
            spectra_output_name=spectra_output_name,
            write_tof_bins=write_tof_bins,
        )

    @staticmethod
    def write_histogram_data(
        nxs_file_path: str,
        tof_bins: np.array,
        panel_size: tuple[int, int],
        remove_event_data: bool,
        spectra_output_name: str = "spectra",
        write_tof_bins: bool = True,
        nproc: int = 8,
    ) -> None:
        """
        Generates histogram spectra for a given detector and writes it to nxs_file
        """

        def delete_event_data(nxs_file, base_dir, panel_name):
            del nxs_file[join(join(base_dir, panel_name), "event_index")]
            del nxs_file[join(join(base_dir, panel_name), "event_id")]
            del nxs_file[join(join(base_dir, panel_name), "event_time_zero")]
            del nxs_file[join(join(base_dir, panel_name), "event_time_offset")]

        print(f"Writing histrogram data in {nxs_file_path}")
        print(f"Remove event data: {remove_event_data}")
        nxs_file = h5py.File(nxs_file_path, "r+")
        base_dir = list(nxs_file.keys())[0]

        panel_names = FormatMANDI.get_panel_names(nxs_file)
        written_tof_bins = False
        for panel_name in panel_names:
            print(f"Processing panel {panel_name}")
            output_path = join(base_dir, panel_name)
            output_path = join(output_path, spectra_output_name)
            print(f"Writing spectra to {output_path}")
            if nxs_file.get(output_path):
                print(f"deleting {output_path}...")
                del nxs_file[output_path]
            panel_spectra = FormatMANDI.generate_histogram_data_for_panel(
                nxs_file, tof_bins, panel_size, panel_name, nproc
            )

            nxs_file.create_dataset(output_path, data=panel_spectra, compression="gzip")
            if remove_event_data:
                delete_event_data(nxs_file, base_dir, panel_name)
                print(f"Removed event data for {panel_name}")
            if write_tof_bins and not written_tof_bins:
                tof_path = join(base_dir, "time_of_flight")
                if nxs_file.get(tof_path):
                    del nxs_file[tof_path]
                print(f"Writing time of flight bins to {tof_path}")
                nxs_file.create_dataset(tof_path, data=tof_bins, compression="gzip")
                written_tof_bins = True

        nxs_file.close()

    @staticmethod
    def compute_event_histogram(
        args: tuple[int, np.array, np.array, np.array],
    ) -> np.array:
        pixel_idx, event_time_offset, corrected_event_id, tof_bins = args
        h, _ = np.histogram(
            event_time_offset[corrected_event_id == pixel_idx], tof_bins
        )
        return h

    @staticmethod
    def generate_histogram_data_for_panel(
        nxs_file: h5py.File,
        tof_bins: np.array,
        panel_size: tuple[int, int],
        panel_name: str,
        nproc=8,
    ) -> np.array:
        """
        Generates histogram data for a given panel
        """

        ## Get panel data
        panel_number = FormatMANDI.get_panel_number(panel_name)
        # Actual pixel ids, starting from bottom left and going up y axis first
        event_id = nxs_file[f"entry/{panel_name}/event_id"][:]
        # Time each event_id was triggered after event_time_zero (sec)
        event_time_offset = nxs_file[f"entry/{panel_name}/event_time_offset"][:]

        # event_ids are given with an offset
        event_id_offset = panel_number * panel_size[0] * panel_size[1]
        corrected_event_id = event_id - event_id_offset

        num_pixels = panel_size[0] * panel_size[1]
        spectra = np.zeros((num_pixels, len(tof_bins) - 1), dtype=np.int32)

        num_cpu = cpu_count()
        if nproc > num_cpu:
            nproc = num_cpu

        pool = Pool(processes=nproc)

        args_list = [
            (i, event_time_offset, corrected_event_id, tof_bins)
            for i in range(num_pixels)
        ]
        results = pool.map(FormatMANDI.compute_event_histogram, args_list)

        pool.close()
        pool.join()

        for i, h in enumerate(results):
            spectra[i] = h

        return spectra

    @staticmethod
    def get_time_range_for_panel(
        nxs_file: h5py.File, panel_size: tuple[float, float], panel_name: str
    ) -> tuple[float, float]:
        """
        Returns the range of event times for a given panel
        """

        def event_data_is_valid(event_id, event_time_offset):
            if len(event_id) == 0 or len(event_time_offset) == 0:
                return False
            return len(event_id) == len(event_time_offset)

        panel_number = FormatMANDI.get_panel_number(panel_name)
        event_index = nxs_file[f"entry/{panel_name}/event_index"]
        event_id = nxs_file[f"entry/{panel_name}/event_id"]
        event_time_zero = nxs_file[f"entry/{panel_name}/event_time_zero"]
        event_time_offset = nxs_file[f"entry/{panel_name}/event_time_offset"]

        if not event_data_is_valid(event_id, event_time_offset):
            return None, None

        num_pixels = panel_size[0] * panel_size[1]
        event_id_offset = panel_number * num_pixels - 1

        raw_event_id = event_id[event_index[0]]
        corrected_event_id = raw_event_id - event_id_offset
        min_event_time = event_time_zero[0] + event_time_offset[corrected_event_id]

        max_idx = int(event_index[-1] - 1)
        raw_event_id = event_id[max_idx]
        corrected_event_id = raw_event_id - event_id_offset
        max_event_time = (
            event_time_zero[max_idx] + event_time_offset[corrected_event_id]
        )

        return min_event_time, max_event_time

    @staticmethod
    def get_time_range_for_dataset(
        nxs_file_path: str, panel_size: tuple[int, int]
    ) -> tuple[float, float]:
        """
        Iterates over num_panels to find the overall min/max tof event recorded
        """

        # TODO get panel_size and panel_number from nxs_file xml

        nxs_file = h5py.File(nxs_file_path, "r")

        min_tof = -1
        max_tof = -1

        panel_names = FormatMANDI.get_panel_names(nxs_file)

        for panel_name in panel_names:
            try:
                min_event_time, max_event_time = FormatMANDI.get_time_range_for_panel(
                    nxs_file, panel_size, panel_name
                )
                if min_event_time is None or max_event_time is None:
                    # Some banks contain no data
                    continue
                if min_tof < 0 or min_event_time < min_tof:
                    min_tof = min_event_time
                if max_event_time > max_tof:
                    max_tof = max_event_time
            except KeyError:
                # Some banks not always present
                pass

        nxs_file.close()

        return min_tof, max_tof

    @staticmethod
    def generate_tof_bins(
        nxs_file: str,
        panel_size: tuple[float, float],
        delta_tof: float = 50,
        padding: float = 100,
    ) -> np.ndarray:
        """
        delta_tof: float (usec)
        padding: float (usec)
        """

        min_tof, max_tof = FormatMANDI.get_time_range_for_dataset(nxs_file, panel_size)
        min_tof = min_tof - padding
        max_tof = max_tof + padding
        print(
            f"Time of flight range for {nxs_file}: {round(min_tof, 3)} - {round(max_tof, 3)} (usec)"
        )
        num_bins = int((max_tof - min_tof) / delta_tof)
        return np.linspace(min_tof, max_tof, num_bins)

    @staticmethod
    def get_panel_names(nxs_file: h5py.File) -> list[str]:
        raw_names = [i for i in nxs_file[list(nxs_file.keys())[0]] if "bank" in i]
        names = []
        for name in raw_names:
            try:
                entry = (name, int("".join([j for j in name if j.isdigit()])))
                names.append(entry)
            except ValueError:  # Other fields with bank in them
                continue
        return [i[0] for i in sorted(names, key=lambda x: x[0])]

    @staticmethod
    def get_panel_number(panel_name: str) -> int:
        return int("".join([i for i in panel_name if i.isdigit()]))


if __name__ == "__main__":
    nxs_file = sys.argv[1]
    remove_event_data = False
    delta_tof = 50  # usec
    FormatMANDI.add_histogram_data_to_nxs_file(
        nxs_file_path=nxs_file, remove_event_data=remove_event_data, delta_tof=delta_tof
    )
