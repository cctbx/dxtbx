from __future__ import annotations

import functools
import sys
import time

import numpy as np

from libtbx.phil import parse
from scitbx.array_family import flex

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format, abstract
from dxtbx.format.FormatMultiImage import FormatMultiImage, Reader
from dxtbx.format.FormatStill import FormatStill
from dxtbx.model import Spectrum
from dxtbx.util.rotate_and_average import rotate_and_average

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
  wavelength_delta_k = 0
    .type = float
    .help = Correction factor, needed during 2014
  wavelength_offset = None
    .type = float
    .help = Optional constant shift to apply to each wavelength
  spectrum_address = FEE-SPEC0
    .type = str
    .help = Address for incident beam spectrometer
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
  spectrum_rotation_angle = None
    .type = float
    .help = If present, first rotate the spectrum image a given amount in \
            degrees. Only applies to 2D spectrometers
  spectrum_pedestal = None
    .type = path
    .help = Path to pickled pedestal file to subtract from the pedestal
  filter {
    evr_address = evr1
      .type = str
      .help = Address for evr object which stores event codes. Should be evr0,\
              evr1, or evr2.
    required_present_codes = None
      .type = int
      .multiple = True
      .help = These codes must be present to keep the event
    required_absent_codes = None
      .type = int
      .multiple = True
      .help = These codes must be absent to keep the event
    pre_filter = False
      .type = bool
      .help = If True, read the event codes for all events up front, and \
              apply the filter then. Otherwise, apply when loading an \
              event.
  }
"""
locator_scope = parse(locator_str)


class XtcReader(Reader):
    def nullify_format_instance(self):
        """No-op for XTC streams. No issue with multiprocessing."""


@abstract
class FormatXTC(FormatMultiImage, FormatStill, Format):
    def __init__(self, image_file, **kwargs):

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self.lazy = kwargs.get("lazy", True)
        FormatMultiImage.__init__(self, **kwargs)
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
        assert self.params.mode in [
            "idx",
            "smd",
        ], "idx or smd mode should be used for analysis (idx is often faster)"

        self._ds = FormatXTC._get_datasource(image_file, self.params)
        self._evr = None
        self.populate_events()

        self._cached_psana_detectors = {}
        self._beam_index = None
        self._beam_cache = None
        self._initialized = True
        self._fee = None

        if self.params.spectrum_pedestal:
            from libtbx import easy_pickle

            self._spectrum_pedestal = easy_pickle.load(self.params.spectrum_pedestal)
        else:
            self._spectrum_pedestal = None

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
        if self.params.mode == "idx":
            if hasattr(self, "times") and len(self.times) > 0:
                return
        elif self.params.mode == "smd":
            if hasattr(self, "run_mapping") and self.run_mapping:
                return

        if not self._psana_runs:
            self._psana_runs = self._get_psana_runs(self._ds)

        self.times = []
        self.run_mapping = {}

        if self.params.mode == "idx":
            for run in self._psana_runs.values():
                times = run.times()
                if (
                    self.params.filter.required_present_codes
                    or self.params.filter.required_absent_codes
                ) and self.params.filter.pre_filter:
                    times = [t for t in times if self.filter_event(run.event(t))]
                self.run_mapping[run.run()] = (
                    len(self.times),
                    len(self.times) + len(times),
                    run,
                )
                self.times.extend(times)
            self.n_images = len(self.times)

        elif self.params.mode == "smd":
            self._ds = FormatXTC._get_datasource(self._image_file, self.params)
            for event in self._ds.events():
                run = event.run()
                if run not in self.run_mapping:
                    self.run_mapping[run] = []
                if (
                    self.params.filter.required_present_codes
                    or self.params.filter.required_absent_codes
                ) and self.params.filter.pre_filter:
                    if self.filter_event(event):
                        self.run_mapping[run].append(event)
                else:
                    self.run_mapping[run].append(event)
            total = 0
            remade_mapping = {}
            for run in sorted(self.run_mapping):
                start = total
                end = len(self.run_mapping[run]) + total
                total += len(self.run_mapping[run])
                events = self.run_mapping[run]
                remade_mapping[run] = start, end, run, events
            self.run_mapping = remade_mapping
            self.n_images = sum(
                [
                    self.run_mapping[r][1] - self.run_mapping[r][0]
                    for r in self.run_mapping
                ]
            )

    def filter_event(self, evt):
        """Return True to keep the event, False to reject it."""
        if not (
            self.params.filter.required_present_codes
            or self.params.filter.required_absent_codes
        ):
            return True
        if not self._evr:
            self._evr = psana.Detector(self.params.filter.evr_address)
        codes = self._evr.eventCodes(evt)

        if self.params.filter.required_present_codes and not all(
            [c in codes for c in self.params.filter.required_present_codes]
        ):
            return False
        if self.params.filter.required_absent_codes and any(
            [c in codes for c in self.params.filter.required_absent_codes]
        ):
            return False
        return True

    def get_run_from_index(self, index=None):
        """Look up the run number given an index"""
        if index is None:
            index = 0
        for run_number in self.run_mapping:
            start, stop, run = self.run_mapping[run_number][0:3]
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
            if self.params.mode == "idx":
                evt = self.get_run_from_index(index).event(self.times[index])
            elif self.params.mode == "smd":
                for run_number in self.run_mapping:
                    start, stop, run, events = self.run_mapping[run_number]
                    if index >= start and index < stop:
                        evt = events[index - start]
            if (
                (
                    self.params.filter.required_present_codes
                    or self.params.filter.required_absent_codes
                )
                and not self.params.filter.pre_filter
                and not self.filter_event(evt)
            ):
                evt = None
            self.current_event = evt
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
        if not evt:
            return None
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
            if not evt:
                self._beam_cache = None
                return None
            spectrum = self.get_spectrum(index)
            if spectrum:
                wavelength = spectrum.get_weighted_wavelength()
            else:
                wavelength = cspad_tbx.evt_wavelength(
                    evt, delta_k=self.params.wavelength_delta_k
                )
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
                if self._beam_cache is not None:
                    self._beam_cache.set_polarization_normal((1, 0, 0))

        return self._beam_cache

    def get_spectrum(self, index=None):
        if index is None:
            index = 0
        if self.params.spectrum_eV_per_pixel is None:
            return None

        evt = self._get_event(index)
        if not evt:
            return None
        if self._fee is None:
            self._fee = psana.Detector(self.params.spectrum_address)
        if self._fee is None:
            return None
        try:
            fee = self._fee.get(evt)
            y = fee.hproj()
            if self.params.spectrum_pedestal:
                y = y - self._spectrum_pedestal.as_numpy_array()

        except AttributeError:  # Handle older spectometers without the hproj method
            try:
                img = self._fee.image(evt)
            except AttributeError:
                return None
            if self.params.spectrum_rotation_angle is None:
                x = (
                    self.params.spectrum_eV_per_pixel * np.array(range(img.shape[1]))
                ) + self.params.spectrum_eV_offset
                y = img.mean(axis=0)  # Collapse 2D image to 1D trace
            else:
                mask = img == 2**16 - 1
                mask = np.invert(mask)

                x, y = rotate_and_average(
                    img, self.params.spectrum_rotation_angle, deg=True, mask=mask
                )
                x = (
                    self.params.spectrum_eV_per_pixel * x
                ) + self.params.spectrum_eV_offset

        else:
            x = (
                self.params.spectrum_eV_per_pixel * np.array(range(len(y)))
            ) + self.params.spectrum_eV_offset
        return Spectrum(flex.double(x), flex.double(y))

    def get_goniometer(self, index=None):
        return None

    def get_scan(self, index=None):
        return None


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatXTC.understand(arg))
