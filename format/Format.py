"""
A top-level class to represent image formats which does little else but
(i) establish an abstract class for what needs to be implemented and
(ii) include the format registration code for any image formats which
inherit from this. This will also contain links to the static methods
from the X(component)Factories which will allow construction of e.g.
goniometers etc. from the headers and hence a format specific factory.

isort:skip_file
"""

from __future__ import absolute_import, division, print_function

from future import standard_library

standard_library.install_aliases()

import functools
import sys
import urllib.parse
import urllib.request
from builtins import range
from os.path import abspath

import libtbx

import dxtbx.filecache_controller
from dxtbx.format.image import ImageBool
from dxtbx.model import MultiAxisGoniometer
from dxtbx.model.beam import BeamFactory
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.scan import ScanFactory
from dxtbx.sequence_filenames import template_regex


if sys.hexversion < 0x3040000:
    # try Python3.3 backport bz2 pypi module first.
    # this supports multiple compression streams.
    # to install run  libtbx.pip install bz2file
    try:
        import bz2file

        bz2 = bz2file
    except ImportError:
        bz2 = None
else:
    bz2 = None

if not bz2:
    try:
        import bz2
    except ImportError:
        bz2 = None

try:
    import gzip
except ImportError:
    gzip = None


# import access to all of the factories that we will be needing


_cache_controller = dxtbx.filecache_controller.simple_controller()


class Reader(object):
    def __init__(self, format_class, filenames, **kwargs):
        self._kwargs = kwargs
        self.format_class = format_class
        self._filenames = filenames

    def read(self, index):
        format_instance = self.format_class.get_instance(
            self._filenames[index], **self._kwargs
        )
        return format_instance.get_raw_data()

    def paths(self):
        return self._filenames

    def identifiers(self):
        return self._filenames

    def __len__(self):
        return len(self._filenames)

    def copy(self, filenames):
        return Reader(self.format_class, filenames)

    def is_single_file_reader(self):
        return False

    def master_path(self):
        return ""


class Format(object):
    """A base class for the representation and interrogation of diffraction
    image formats, from which all classes for reading the header should be
    inherited. This includes: autoregistration of implementation classes,
    stubs which need to be overridden and links to static factory methods
    which will prove to be useful in other implementations."""

    @staticmethod
    def understand(image_file):
        """Overload this to publish whether this class instance
        understands a given file.  N.B. to say that we understand it,
        return True.  If a subclass also understands the image
        (because, for example, its detector serial number takes a
        certain value), it will by definition understand it better
        than its superclass.  Thus, the preferred class will be the
        deepest subclass in the inheritance hierarchy.  Finally, if
        you are writing this subclass for a given instrument and you
        are given a different example return False.

        Implementing understand() in a subclass, one can safely assume
        that the superclasss understand() function returned True.
        The understand() function of two different classes directly
        derived from the same base should never both return True for
        the same input image."""

        return False

    @classmethod
    def ignore(cls):
        return False

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return False
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialize a class instance from an image file."""

        self._image_file = image_file

        self._goniometer_instance = None
        self._detector_instance = None
        self._beam_instance = None
        self._scan_instance = None

        self._goniometer_factory = GoniometerFactory
        self._detector_factory = DetectorFactory
        self._beam_factory = BeamFactory
        self._scan_factory = ScanFactory

        self.setup()

    def setup(self):
        """Read the image file, construct the information which we will be
        wanting about the experiment from this. N.B. in your implementation
        of this you will probably want to make use of the static methods
        below and probably add some format parsing code too. Please also keep
        in mind that your implementation may be further subclassed by
        someone else."""

        self._start()

        try:
            goniometer_instance = self._goniometer()
            self._goniometer_instance = goniometer_instance

            detector_instance = self._detector()
            self._detector_instance = detector_instance

            beam_instance = self._beam()
            self._beam_instance = beam_instance

            scan_instance = self._scan()
            self._scan_instance = scan_instance
        finally:
            self._end()

    @staticmethod
    def get_cache_controller():
        return _cache_controller

    def get_goniometer(self):
        """Get the standard goniometer instance which was derived from the
        image headers."""

        return self._goniometer_instance

    def get_detector(self):
        """Get the standard detector instance which was derived from the
        image headers."""

        return self._detector_instance

    def get_beam(self):
        """Get the standard beam instance which was derived from the image
        headers."""

        return self._beam_instance

    def get_scan(self):
        """Get the standard scan instance which was derived from the image
        headers."""

        return self._scan_instance

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array."""
        try:
            image = self.detectorbase
        except AttributeError:
            return None
        image.read()
        raw_data = image.get_raw_data()
        return raw_data

    def get_vendortype(self):
        return "no dxtbx Format vendortype"

    def detectorbase_start(self):
        raise NotImplementedError

    def get_detectorbase(self):
        """Return the instance of detector base."""
        self.detectorbase_start()

        # XXX Temporary proxy to aid in the transition to dxtbx.  Remove
        # once completed.
        class _detectorbase_proxy(object):
            def __init__(self, format_instance):
                self._fi = format_instance
                if not hasattr(self, "vendortype"):
                    self.vendortype = "generic"

            def __getattribute__(self, name):
                if name == "__class__":
                    return self._fi.__class__
                return object.__getattribute__(self, name)

            def __getattr__(self, name):
                try:
                    return self._fi.__getattribute__(name)
                except AttributeError:
                    # print >> sys.stderr, \
                    #    "requesting iotbx.detectors.detectorbase.%s" % name
                    try:
                        return self._fi.detectorbase.__getattribute__(name)
                    except AttributeError:
                        return self._fi.detectorbase.__getattr__(name)

        return _detectorbase_proxy(self)

    @classmethod
    def get_instance(Class, filename, **kwargs):
        if (
            not hasattr(Class, "_current_instance_")
            or Class._current_filename_ != filename
            or Class._current_kwargs_ != kwargs
        ):
            Class._current_instance_ = Class(filename, **kwargs)
            Class._current_filename_ = filename
            Class._current_kwargs_ = kwargs
        return Class._current_instance_

    @classmethod
    def get_reader(cls):
        """
        Return a reader class
        """
        return functools.partial(Reader, cls)

    def get_masker(self, goniometer=None):
        """
        Return a masker class
        """
        if (
            isinstance(goniometer, MultiAxisGoniometer)
            and hasattr(self, "_dynamic_shadowing")
            and self._dynamic_shadowing
        ):
            masker = self.get_goniometer_shadow_masker(goniometer=goniometer)
        else:
            masker = None
        return masker

    @classmethod
    def get_imageset(
        Class,
        filenames,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
        as_imageset=False,
        as_sequence=False,
        single_file_indices=None,
        format_kwargs=None,
        template=None,
        check_format=True,
    ):
        """
        Factory method to create an imageset

        """
        # Import here to avoid cyclic imports
        from dxtbx.imageset import ImageSet, ImageSetData, ImageSequence

        # Get filename absolute paths
        filenames = tuple(map(abspath, filenames))

        # Make it a dict
        if format_kwargs is None:
            format_kwargs = {}

        # Get some information from the format class
        reader = Class.get_reader()(filenames, **format_kwargs)

        # Get the format instance
        if check_format is True:
            format_instance = Class(filenames[0], **format_kwargs)
        else:
            format_instance = None

        # Read the vendor type
        if check_format is True:
            vendor = format_instance.get_vendortype()
        else:
            vendor = ""

        # Get the format kwargs
        params = format_kwargs

        # Make sure only 1 or none is set
        assert [as_imageset, as_sequence].count(True) < 2
        if as_imageset:
            is_sequence = False
        elif as_sequence:
            is_sequence = True
        else:
            if scan is None and format_instance is None:
                raise RuntimeError(
                    """
          One of the following needs to be set
            - as_imageset=True
            - as_sequence=True
            - scan
            - check_format=True
      """
                )
            if scan is None:
                test_scan = format_instance.get_scan()
            else:
                test_scan = scan
            if test_scan is not None and test_scan.get_oscillation()[1] != 0:
                is_sequence = True
            else:
                is_sequence = False

        # Create an imageset or sequence
        if not is_sequence:

            # Create the imageset
            iset = ImageSet(
                ImageSetData(
                    reader=reader,
                    masker=None,
                    vendor=vendor,
                    params=params,
                    format=Class,
                )
            )

            # If any are None then read from format
            if [beam, detector, goniometer, scan].count(None) != 0:

                # Get list of models
                beam = []
                detector = []
                goniometer = []
                scan = []
                for f in filenames:
                    format_instance = Class(f, **format_kwargs)
                    beam.append(format_instance.get_beam())
                    detector.append(format_instance.get_detector())
                    goniometer.append(format_instance.get_goniometer())
                    scan.append(format_instance.get_scan())

            # Set the list of models
            for i in range(len(filenames)):
                iset.set_beam(beam[i], i)
                iset.set_detector(detector[i], i)
                iset.set_goniometer(goniometer[i], i)
                iset.set_scan(scan[i], i)

        else:

            # Get the template
            if template is None:
                template = template_regex(filenames[0])[0]
            else:
                template = str(template)

            # Check scan makes sense
            if scan:
                if check_format is True:
                    assert scan.get_num_images() == len(filenames)

            # If any are None then read from format
            if beam is None and format_instance is not None:
                beam = format_instance.get_beam()
            if detector is None and format_instance is not None:
                detector = format_instance.get_detector()
            if goniometer is None and format_instance is not None:
                goniometer = format_instance.get_goniometer()
            if scan is None and format_instance is not None:
                scan = format_instance.get_scan()
                if scan is not None:
                    for f in filenames[1:]:
                        format_instance = Class(f, **format_kwargs)
                        scan += format_instance.get_scan()

            assert beam is not None, "Can't create Sequence without beam"
            assert detector is not None, "Can't create Sequence without detector"
            assert goniometer is not None, "Can't create Sequence without goniometer"
            assert scan is not None, "Can't create Sequence without scan"

            # Create the masker
            if format_instance is not None:
                masker = format_instance.get_masker(goniometer=goniometer)
            else:
                masker = None

            # Create the sequence
            iset = ImageSequence(
                ImageSetData(
                    reader=reader,
                    masker=masker,
                    vendor=vendor,
                    params=params,
                    format=Class,
                    template=template,
                ),
                beam=beam,
                detector=detector,
                goniometer=goniometer,
                scan=scan,
            )

        if format_instance is not None:
            static_mask = format_instance.get_static_mask()
            if static_mask is not None:
                if not iset.external_lookup.mask.data.empty():
                    for m1, m2 in zip(static_mask, iset.external_lookup.mask.data):
                        m1 &= m2.data()
                    iset.external_lookup.mask.data = ImageBool(static_mask)
                else:
                    iset.external_lookup.mask.data = ImageBool(static_mask)

        # Return the imageset
        return iset

    def get_image_file(self):
        """Get the image file provided to the constructor."""

        return self._image_file

    # methods which must be overloaded in order to produce a useful Format
    # class implementation

    def _start(self):
        """Start code for handling this image file, which may open a link
        to it once, say, and pass this around within the implementation.
        How you use this is up to you, though you probably want to overload
        it..."""

        return

    def _end(self):
        """Clean up things - keeping in mind that this should run even in the
        case of an exception being raised."""

        return

    def _goniometer(self):
        """Overload this method to read the image file however you like so
        long as the result is an goniometer."""
        return None

    def _detector(self):
        """Overload this method to read the image file however you like so
        long as the result is an detector."""
        return None

    def _beam(self):
        """Overload this method to read the image file however you like so
        long as the result is an beam."""
        return None

    def _scan(self):
        """Overload this method to read the image file however you like so
        long as the result is an scan."""
        return None

    def get_static_mask(self):
        """Overload this method to override the static mask."""
        return None

    def get_goniometer_shadow_masker(self, goniometer=None):
        """Overload this method to allow generation of dynamic goniometer shadow
        masks to be used during spotfinding or integration."""

        return None

    ####################################################################
    #                                                                  #
    # Helper functions for dealing with compressed images.             #
    #                                                                  #
    ####################################################################

    @staticmethod
    def is_url(path):
        """See if the file is a URL."""

        # Windows file paths can get caught up in this - check that the
        # first letter is one character (which I think should be safe: all
        # URL types are longer than this right?)
        scheme = urllib.parse.urlparse(path).scheme
        return scheme and len(scheme) != 1

    @classmethod
    def open_file(cls, filename, mode="rb", url=False):
        """Open file for reading, decompressing silently if necessary,
        caching transparently if possible."""

        if url and Format.is_url(filename):
            fh_func = functools.partial(urllib.request.urlopen, filename)

        elif filename.endswith(".bz2"):
            if not bz2:
                raise RuntimeError("bz2 file provided without bz2 module")
            fh_func = functools.partial(bz2.BZ2File, filename, mode=mode)

        elif filename.endswith(".gz"):
            if not gzip:
                raise RuntimeError("gz file provided without gzip module")
            fh_func = functools.partial(gzip.GzipFile, filename, mode=mode)

        else:
            fh_func = functools.partial(open, filename, mode=mode)

        ##  To disable caching logic:
        # return fh_func()
        return cls.get_cache_controller().check(filename, fh_func)
