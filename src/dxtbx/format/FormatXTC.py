import functools
import sys
import time

import numpy as np

from cctbx import factor_ev_angstrom
from libtbx.phil import parse

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format, abstract
from dxtbx.format.FormatMultiImage import Reader
from dxtbx.format.FormatMultiImageLazy import FormatMultiImageLazy
from dxtbx.format.FormatStill import FormatStill

try:
    import psana

    from xfel.cxi.cspad_ana import cspad_tbx
except ImportError:
    psana = None
    cspad_tbx = None

locator_str = """
  experiment = None
    .type = str
    .help = Experiment identifier, e.g. mfxo1916
  run = None
    .type = ints
    .help = Run number or a list of runs to process
  mode = idx
    .type = str
    .help = Mode for reading the xtc data (see LCLS documentation)
  data_source = None
    .type = str
    .help = Complete LCLS data source.  Overrides experiment and run.  Example: \
            exp=mfxo1916:run=20:smd \
            More info at https://confluence.slac.stanford.edu/display/PSDM/Manual#Manual-Datasetspecification
  detector_address = None
    .type = str
    .multiple = True
    .help = detector used for collecting the data at LCLS
  calib_dir = None
    .type = str
    .help = Specify path to custom calib directory if needed
  use_ffb = False
    .type = bool
    .help = Run on the ffb if possible. Only for active users!
  wavelength_offset = None
    .type = float
    .help = Optional constant shift to apply to each wavelength
  spectrum_eV_per_pixel = None
    .type = float
    .help = If not None, use the FEE spectrometer to determine the wavelength. \
            spectrum_eV_offset should be also specified. A weighted average of \
            the horizontal projection of the per-shot FEE spectrometer is used.\
            The equation for each pixel eV is \
            eV = (spectrum_eV_per_pixel * pixel_number) + spectrum_eV_offset
  spectrum_eV_offset = None
    .type = float
    .help = See spectrum_eV_per_pixel
"""
locator_scope = parse(locator_str)


class XtcReader(Reader):
    def nullify_format_instance(self):
        """No-op for XTC streams. No issue with multiprocessing."""
        pass


@abstract
class FormatXTC(FormatMultiImageLazy, FormatStill, Format):
    def __init__(self, image_file, **kwargs):

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatMultiImageLazy.__init__(self, **kwargs)
        FormatStill.__init__(self, image_file, **kwargs)
        Format.__init__(self, image_file, **kwargs)
        self.current_index = None
        self.current_event = None
        self._psana_runs = {}  # empty container, to prevent breaking other formats
        if "locator_scope" in kwargs:
            self.params = FormatXTC.params_from_phil(
                master_phil=kwargs["locator_scope"], user_phil=image_file, strict=True
            )
        else:
            self.params = FormatXTC.params_from_phil(
                master_phil=locator_scope, user_phil=image_file, strict=True
            )
        assert self.params.mode == "idx", "idx mode should be used for analysis"

        self._ds = FormatXTC._get_datasource(image_file, self.params)
        self.populate_events()
        self.n_images = len(self.times)

        self._cached_psana_detectors = {}
        self._beam_index = None
        self._beam_cache = None
        self._initialized = True
        self._fee = None

    @staticmethod
    def understand(image_file):
        """Extracts the datasource and detector_address from the image_file and then feeds it to PSANA
        If PSANA fails to read it, then input may not be an xtc/smd file. If success, then OK.
        If detector_address is not provided, a command line promp will try to get the address
        from the user"""
        if not psana or not cspad_tbx:
            return False
        try:
            params = FormatXTC.params_from_phil(locator_scope, image_file)
        except Exception:
            return False
        if params is None:
            return False

        try:
            FormatXTC._get_datasource(image_file, params)
        except Exception:
            return False
        return True

    @staticmethod
    def params_from_phil(master_phil, user_phil, strict=False):
        """Read the locator file"""
        try:
            user_input = parse(file_name=user_phil)
            working_phil, unused = master_phil.fetch(
                sources=[user_input], track_unused_definitions=True
            )
            unused_args = ["%s=%s" % (u.path, u.object.words[0].value) for u in unused]
            if len(unused_args) > 0 and strict:
                for unused_arg in unused_args:
                    print(unused_arg)
                print(
                    "Incorrect of unused parameter in locator file. Please check and retry"
                )
                return None
            params = working_phil.extract()
            return params
        except Exception:
            return None

    @classmethod
    def get_reader(cls):
        """
        Return a reader class
        """
        return functools.partial(XtcReader, cls)

    def populate_events(self):
        """Read the timestamps from the XTC stream.  Assumes the psana idx mode of reading data.
        Handles multiple LCLS runs by concatenating the timestamps from multiple runs together
        in a single list and creating a mapping."""
        if hasattr(self, "times") and len(self.times) > 0:
            return

        if not self._psana_runs:
            self._psana_runs = self._get_psana_runs(self._ds)

        self.times = []
        self.run_mapping = {}
        for run in self._psana_runs.values():
            times = run.times()
            self.run_mapping[run.run()] = (
                len(self.times),
                len(self.times) + len(times),
                run,
            )
            self.times.extend(times)

    def get_run_from_index(self, index=None):
        """Look up the run number given an index"""
        if index is None:
            index = 0
        for run_number in self.run_mapping:
            start, stop, run = self.run_mapping[run_number]
            if index >= start and index < stop:
                return run
        raise IndexError("Index is not within bounds")

    def _get_event(self, index):
        """Retrieve a psana event given and index. This is the slow step for reading XTC streams,
        so implement a cache for the last read event."""
        if index == self.current_index:
            return self.current_event
        else:
            self.current_index = index
            self.current_event = self.get_run_from_index(index).event(self.times[index])
            return self.current_event

    @staticmethod
    def _get_datasource(image_file, params):
        """Construct a psana data source object given the locator parameters"""
        if params.calib_dir is not None:
            psana.setOption("psana.calib-dir", params.calib_dir)
        if params.data_source is None:
            if (
                params.experiment is None
                or params.run is None
                or params.mode is None
                or len(params.run) == 0
            ):
                return False
            img = "exp=%s:run=%s:%s" % (
                params.experiment,
                ",".join(["%d" % r for r in params.run]),
                params.mode,
            )

            if params.use_ffb:
                # as ffb is only at SLAC, ok to hardcode /reg/d here
                img += ":dir=/reg/d/ffb/%s/%s/xtc" % (
                    params.experiment[0:3],
                    params.experiment,
                )
        else:
            img = params.data_source
        return psana.DataSource(img)

    @staticmethod
    def _get_psana_runs(datasource):
        """
        Extracts the runs,
        These can only be extracted once,
        only call this method after datasource is set
        """
        # this is key,value = run_integer, psana.Run, e.g. {62: <psana.Run(@0x7fbd0e23c990)>}
        psana_runs = {r.run(): r for r in datasource.runs()}
        return psana_runs

    def _get_psana_detector(self, run):
        """Returns the psana detector for the given run"""
        if run.run() not in self._cached_psana_detectors:
            assert len(self.params.detector_address) == 1
            self._cached_psana_detectors[run.run()] = psana.Detector(
                self.params.detector_address[0], run.env()
            )
        return self._cached_psana_detectors[run.run()]

    def get_psana_timestamp(self, index):
        """Get the cctbx.xfel style event timestamp given an index"""
        evt = self._get_event(index)
        time = evt.get(psana.EventId).time()
        # fid = evt.get(psana.EventId).fiducials()

        sec = time[0]
        nsec = time[1]

        return cspad_tbx.evt_timestamp((sec, nsec / 1e6))

    def get_num_images(self):
        return self.n_images

    def get_beam(self, index=None):
        return self._beam(index)

    def _beam(self, index=None):
        """Returns a simple model for the beam"""
        if index is None:
            index = 0
        if self._beam_index != index:
            self._beam_index = index
            evt = self._get_event(index)
            if self.params.spectrum_eV_per_pixel is not None:
                if self._fee is None:
                    self._fee = psana.Detector("FEE-SPEC0")
                fee = self._fee.get(evt)
                if fee is None:
                    wavelength = cspad_tbx.evt_wavelength(evt)
                else:
                    x = (
                        self.params.spectrum_eV_per_pixel
                        * np.array(range(len(fee.hproj())))
                    ) + self.params.spectrum_eV_offset
                    wavelength = factor_ev_angstrom / np.average(x, weights=fee.hproj())
            else:
                wavelength = cspad_tbx.evt_wavelength(evt)
            if wavelength is None:
                self._beam_cache = None
            else:
                if self.params.wavelength_offset is not None:
                    wavelength += self.params.wavelength_offset
                self._beam_cache = self._beam_factory.simple(wavelength)
            s, nsec = evt.get(psana.EventId).time()
            evttime = time.gmtime(s)
            if (
                evttime.tm_year == 2020 and evttime.tm_mon >= 7
            ) or evttime.tm_year > 2020:
                self._beam_cache.set_polarization_normal((1, 0, 0))

        return self._beam_cache

    def get_goniometer(self, index=None):
        return None

    def get_scan(self, index=None):
        return None


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatXTC.understand(arg))
