from __future__ import absolute_import, division, print_function
import numpy as np
import glob
import os
import re
from libtbx.phil import parse
from dxtbx.format.FormatTIFFRayonix import FormatTIFFRayonix
from dxtbx.format.FormatStill import FormatStill
from dxtbx.format.FormatMultiImage import Reader
from dxtbx.format.FormatMultiImageLazy import FormatMultiImageLazy
from dxtbx.format.FormatPhilLocator import FormatPhilLocator

"""
This format class only works on the queuing systems at SACLA
due to how the library dbpy is built.

The class uses a locator file with the below specification
to read TIFF images from the SACLA Rayonix detector and
combine them with per-shot wavelength data from the SACLA
dbpy stores.
"""

locator_str = """
  beamline = 2
    .type = int
    .help = Beamline number (1, 2 or 3)
  experiment = None
    .type = str
    .help = Experiment identifier, e.g. 2019A8088
  run = None
    .type = int
    .help = Run number
  rayonix_root = "/xustrg0"
    .type = path
    .help = Path to tiff files
"""
locator_scope = parse(locator_str)

class SaclaTiffReader(Reader):
    def nullify_format_instance(self):
        """ No-op. No issue with multiprocessing. """
        pass

def str2float(s):
    m = re.match("-?\d+(.\d+)?(e[+-]?\d+)?", s)
    if m is not None:
        return float(m.group(0))
    else:
        return float("nan")

class FormatTIFFRayonixSACLA(FormatPhilLocator, FormatMultiImageLazy, FormatStill, FormatTIFFRayonix):
    def __init__ (self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError
        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self.image_file=image_file
        if "locator_scope" in kwargs:
            self.params = FormatTIFFRayonixSACLA.params_from_phil(
                master_phil=kwargs["locator_scope"], user_phil=image_file, strict=True
            )
        else:
            self.params = FormatTIFFRayonixSACLA.params_from_phil(
                master_phil=locator_scope, user_phil=image_file, strict=True
            )
        self.img_files = sorted(glob.glob(os.path.join(self.params.rayonix_root, self.params.experiment, str(self.params.run), "data_*.img")))

        FormatMultiImageLazy.__init__(self, **kwargs)
        FormatStill.__init__(self, image_file, **kwargs)

        # Work arund to get FormatTIFFRayonix to understand this file
        # during its init method, since it cannot understand the locator
        tmp = self.understand
        self.understand = FormatTIFFRayonix.understand
        FormatTIFFRayonix.__init__(self, self.img_files[0], **kwargs) # init base class using first image in run
        self.understand = tmp

    @staticmethod
    def understand(image_file):
        try:
            # Following import only works on a compute node on sacla hpc
            import dbpy
        except Exception:
            return False
        try:
            params = FormatTIFFRayonixSACLA.params_from_phil(locator_scope, image_file, strict=True)
        except Exception:
            return False
        return params is not None

    @staticmethod
    def params_from_phil(master_phil, user_phil, strict=False):
        """ Read the locator file """
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
                    "Incorrect or unused parameter in locator file. Please check and retry"
                )
                return None
            params = working_phil.extract()
            return params
        except Exception:
            return None

    @classmethod
    def get_reader(Class):
        """
        Return a reader class

        """
        obj = SaclaTiffReader
        # Note, need to set this on the parent class since it's a scoped global variable
        Reader._format_class_ = Class
        return obj

    def _start(self):
        import dbpy

        sensor_spec = "xfel_bl_%d_tc_spec_1/energy"%self.params.beamline
        sensor_shutter = "xfel_bl_%d_shutter_1_open_valid/status"%self.params.beamline
        # Get run info
        try:
            run_info = dbpy.read_runinfo(self.params.beamline, self.params.run)
            high_tag = dbpy.read_hightagnumber(self.params.beamline, self.params.run)
        except Exception as e:
            return -1
        start_tag = run_info['start_tagnumber']
        end_tag = run_info['end_tagnumber']

        tag_list = np.array(dbpy.read_taglist_byrun(self.params.beamline, self.params.run))
        try:
            shutter = np.array(map(str2float, dbpy.read_syncdatalist(sensor_shutter, high_tag, tuple(tag_list))))
        except Exception as e:
            return -1

        high_tag=dbpy.read_hightagnumber(self.params.beamline, self.params.run)
        valid_tags=tag_list[shutter==1]
        if len(self.img_files)+1 != len(valid_tags): # last valid tag is not saved.
            print ("# WARNING!! img_files and valid_tag number mismatch")

            img_numbers = map(lambda x: int(x[x.rindex("_")+1:-4]), self.img_files)
            dropped_frames = sorted(set(range(1, len(valid_tags))).difference(img_numbers))
            print ("# Unsaved frame numbers =", tuple(dropped_frames))
            print ("# DEBUG::", len(self.img_files)-len(dropped_frames)+1, len(valid_tags))
            if len(self.img_files)+len(dropped_frames)+1 == len(valid_tags):
                print ("#  %d unsaved img files found, which explains number mismatch" % len(dropped_frames))
                valid_tags = np.delete(valid_tags, np.array(dropped_frames)-1)
                assert len(self.img_files)+1 == len(valid_tags)
            else:
                print ("# Assuming last %d img files are generated after stopping run.." % (len(self.img_files)-len(valid_tags)+1))
                self.img_files = self.img_files[:len(valid_tags)-1]
                assert len(self.img_files)+1 == len(valid_tags)

        self.photon_energies_in_keV=np.array([str2float(s) for s in dbpy.read_syncdatalist(sensor_spec, high_tag, tuple(valid_tags))])

        loc = self._image_file
        self._image_file = self.img_files[0]
        FormatTIFFRayonix._start(self)
        self._image_file = loc

    def get_num_images(self):
        return len(self.img_files)

    def get_detector(self, index=None):
        return self._detector(index)

    def get_beam(self, index=None):
        return self._beam(index)

    def _beam(self, index=None):
        if index is None: index = 0

        #wavelength=12.3984/self.photon_energies_in_keV[index] # cheetah uses this multiplier
        wavelength=12.3984187/self.photon_energies_in_keV[index] # this is more accurate multiplier
        return self._beam_factory.simple(wavelength)

    def _detector(self, index=None):
        if index is None: index = 0
        loc = self._image_file
        self._image_file = self.img_files[index]
        detector = FormatTIFFRayonix._detector(self)
        self._image_file = loc
        return detector

    def get_raw_data(self, index=None):
        if index is None: index = 0
        loc = self._image_file
        self._image_file = self.img_files[index]
        data = FormatTIFFRayonix.get_raw_data(self)
        self._image_file = loc
        return data

if __name__ == "__main__":
    import sys
    for arg in sys.argv[1:]:
        print(FormatTIFFRayonixSACLA.understand(arg))
