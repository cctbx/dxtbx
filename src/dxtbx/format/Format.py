"""
A top-level class to represent image formats which does little else but
(i) establish an abstract class for what needs to be implemented and
(ii) include the format registration code for any image formats which
inherit from this. This will also contain links to the static methods
from the X(component)Factories which will allow construction of e.g.
goniometers etc. from the headers and hence a format specific factory.
"""

import bz2
import functools
import os
from pathlib import Path
from typing import ClassVar, List

import libtbx
from libtbx.utils import Sorry

import dxtbx.filecache_controller
from dxtbx.format.image import ImageBool
from dxtbx.model import MultiAxisGoniometer
from dxtbx.model.beam import MonochromaticBeamFactory
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.sequence import SequenceFactory
from dxtbx.sequence_filenames import template_regex
from dxtbx.util import get_url_scheme

try:
    import gzip
except ImportError:
    gzip = None


_cache_controller = dxtbx.filecache_controller.simple_controller()


def abstract(cls):
    """
    Mark a Format class as abstract, but not it's subclasses.

    This means that it should not be directly used to parse user images.
    """

    cls._abstract_base = cls
    return cls


class Reader:
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


@abstract
class Format:
    """A base class for the representation and interrogation of diffraction
    image formats, from which all classes for reading the header should be
    inherited. This includes: autoregistration of implementation classes,
    stubs which need to be overridden and links to static factory methods
    which will prove to be useful in other implementations.

    Attributes:
        schemes: List of schemes this Format class supports.
            An empty string ("") is treated as non-URL filenames, where
            "file" is treated as explicit file:// support. Format classes
            that do not support the listed schema will _not_ be called
            to understand() with the locator - note that this means if a
            class does not declare "file", it will not get called to
            understand "file://" URIs passed to it. An empty `schemes` will
            not ever get called in Format registry searches.
    """

    schemes = [""]  # type: List[str]

    # Which class is the abstract base. Assigned with the @abstract
    # decorator, and used to check initialization
    _abstract_base: "ClassVar[Format]" = None

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
    def is_abstract(cls) -> bool:
        """
        Is this Format class allowed to be instantiated?

        Format classes that describe an abstract class of formats but do
        not explicitly read a file-format should be marked as such with
        the dxtbx.format.abstract decorator, and they will not be
        instantiable or used as a registry match via load lookups.
        """
        return cls is cls._abstract_base

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return False
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialize a class instance from an image file."""

        # Don't allow abstract instantion
        # - except for Format, which is used as a placeholder in many
        # places (e.g. still, check_format=False) so needs to be allowed.
        if self.is_abstract() and not type(self) is Format:
            raise TypeError(
                f"Cannot instantiate: Format class '{type(self).__name__}' is marked abstract"
            )

        self._image_file = image_file

        self._goniometer_instance = None
        self._detector_instance = None
        self._beam_instance = None
        self._sequence_instance = None

        self._goniometer_factory = GoniometerFactory
        self._detector_factory = DetectorFactory
        self._beam_factory = MonochromaticBeamFactory
        self._sequence_factory = SequenceFactory

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

            sequence_instance = self._sequence()
            self._sequence_instance = sequence_instance
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

    def get_sequence(self):
        """Get the standard sequence instance which was derived from the image
        headers."""

        return self._sequence_instance

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

        # Proxy to link dxtbx to iotbx
        class _detectorbase_proxy:
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

    @staticmethod
    def identify_imageset_type(sequence, goniometer, beam, format_instance):

        from dxtbx.imageset import ImageSetType
        from dxtbx.model.beam import TOFBeam
        from dxtbx.model.sequence import Scan, TOFSequence

        def is_tof_imagesequence(sequence, beam, format_instance):
            if beam is None and format_instance is not None:
                beam = format_instance.get_beam()
            if sequence is None and format_instance is not None:
                sequence = format_instance.get_sequence()
            if isinstance(beam, TOFBeam):
                if not isinstance(sequence, TOFSequence):
                    raise NotImplementedError(
                        "ToFBeam only implemented to work with TOFSequence"
                    )
                else:
                    return True
            return False

        def is_rot_imagesequence(sequence, goniometer, format_instance):
            if sequence is not None and goniometer is not None:
                if isinstance(sequence, Scan):
                    return True
            if format_instance is not None:
                sequence = format_instance.get_sequence()
                goniometer = format_instance.get_goniometer()
                if sequence is not None and goniometer is not None:
                    if isinstance(sequence, Scan):
                        return True
            return False

        if is_tof_imagesequence(sequence, beam, format_instance):
            return ImageSetType.ImageSequence
        elif is_rot_imagesequence(sequence, goniometer, format_instance):
            return ImageSetType.ImageSequence
        return ImageSetType.ImageSet

    @staticmethod
    def format_instance_required(beam, detector, imageset_type):

        if beam is None or detector is None:
            return True
        return False

    @classmethod
    def get_imageset(
        cls,
        filenames,
        beam=None,
        detector=None,
        goniometer=None,
        sequence=None,
        imageset_type=None,
        format_kwargs=None,
        **kwargs,
    ):
        """
        Factory method to create an imageset

        """
        # Import here to avoid cyclic imports
        from dxtbx.imageset import ImageSequence, ImageSet, ImageSetData, ImageSetType

        def process_filenames(filenames):
            return [
                os.path.abspath(x) if not get_url_scheme(x) else x for x in filenames
            ]

        def process_format_kwargs(format_kwargs):
            if format_kwargs is None:
                format_kwargs = {}
            return format_kwargs

        def create_imagesequence(
            cls,
            filenames,
            beam,
            detector,
            goniometer,
            sequence,
            format_instance,
            format_kwargs,
            imageset_type,
            **kwargs,
        ):
            # Get the template
            if "template" not in kwargs or kwargs["template"] is None:
                template = template_regex(filenames[0])[0]
            else:
                template = str(kwargs["template"])

            # Check sequence makes sense
            if sequence:
                assert sequence.get_num_images() == len(filenames)

            # If any are None then read from format
            if beam is None and format_instance is not None:
                beam = format_instance.get_beam()
            if detector is None and format_instance is not None:
                detector = format_instance.get_detector()
            if goniometer is None and format_instance is not None:
                goniometer = format_instance.get_goniometer()
            if sequence is None and format_instance is not None:
                sequence = format_instance.get_sequence()
                if sequence is not None:
                    for f in filenames[1:]:
                        format_instance = cls(f, **format_kwargs)
                        sequence += format_instance.get_sequence()

            assert beam is not None, "Can't create ImageSequence without beam"
            assert detector is not None, "Can't create ImageSequence without detector"
            assert (
                goniometer is not None
            ), "Can't create ImageSequence without goniometer"
            assert sequence is not None, "Can't create ImageSequence without sequence"

            # Create the masker
            if format_instance is not None:
                masker = format_instance.get_masker(goniometer=goniometer)
            else:
                masker = None

            # Create the sequence
            imageset_type = Format.identify_imageset_type(
                sequence=sequence,
                goniometer=goniometer,
                beam=beam,
                format_instance=format_instance,
            )

            if imageset_type == ImageSetType.ImageSequence:
                iset_data = ImageSetData(
                    reader=reader,
                    masker=masker,
                    vendor=vendor,
                    params=format_kwargs,
                    format=cls,
                    template=template,
                )
                iset = ImageSequence(
                    iset_data,
                    beam=beam,
                    detector=detector,
                    goniometer=goniometer,
                    sequence=sequence,
                )
            else:
                raise (
                    NotImplementedError,
                    "Unknown ImageSequence type %s" % imageset_type,
                )

            return iset

        def create_imageset(
            cls,
            filenames,
            beam,
            detector,
            goniometer,
            sequence,
            format_instance,
            vendor,
            format_kwargs,
        ):
            # Create the imageset
            iset = ImageSet(
                ImageSetData(
                    reader=reader,
                    masker=None,
                    vendor=vendor,
                    params=format_kwargs,
                    format=cls,
                )
            )

            # If any are None then read from format
            if [beam, detector, goniometer, sequence].count(None) != 0:

                # Get list of models
                beam = []
                detector = []
                goniometer = []
                sequence = []
                for f in filenames:
                    format_instance = cls(f, **format_kwargs)
                    beam.append(format_instance.get_beam())
                    detector.append(format_instance.get_detector())
                    goniometer.append(format_instance.get_goniometer())
                    sequence.append(format_instance.get_sequence())

            # Set the list of models
            for i in range(len(filenames)):
                iset.set_beam(beam[i], i)
                iset.set_detector(detector[i], i)
                iset.set_goniometer(goniometer[i], i)
                iset.set_sequence(sequence[i], i)

            return iset

        # Process input
        filenames = process_filenames(filenames)
        format_kwargs = process_format_kwargs(format_kwargs)

        # Get some information from the format class
        reader = cls.get_reader()(filenames, **format_kwargs)
        if Format.format_instance_required(beam, detector, imageset_type):
            cls._current_filename_ = None
            cls._current_instance_ = None
            format_instance = cls.get_instance(filenames[0], **format_kwargs)
            vendor = format_instance.get_vendortype()
        else:
            format_instance = None
            vendor = ""

        # Attempt to identify imageset type from models if type is not given
        if imageset_type is None:
            imageset_type = Format.identify_imageset_type(
                sequence, goniometer, beam, format_instance
            )

        if imageset_type == ImageSetType.ImageSet:
            iset = create_imageset(
                cls,
                filenames,
                beam,
                detector,
                goniometer,
                sequence,
                format_instance,
                vendor,
                format_kwargs,
            )
        elif imageset_type == ImageSetType.ImageSequence:
            iset = create_imagesequence(
                cls,
                filenames,
                beam,
                detector,
                goniometer,
                sequence,
                format_instance,
                format_kwargs,
                imageset_type=imageset_type,
                **kwargs,
            )
        else:
            raise NotImplementedError(f"Imageset_type {imageset_type} not implemented")

        if format_instance is not None:
            static_mask = format_instance.get_static_mask()
            if static_mask is not None:
                if not iset.external_lookup.mask.data.empty():
                    for m1, m2 in zip(static_mask, iset.external_lookup.mask.data):
                        m1 &= m2.data()
                    iset.external_lookup.mask.data = ImageBool(static_mask)
                else:
                    iset.external_lookup.mask.data = ImageBool(static_mask)

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
        long as the result is a goniometer."""
        return None

    def _detector(self):
        """Overload this method to read the image file however you like so
        long as the result is a detector."""
        return None

    def _beam(self):
        """Overload this method to read the image file however you like so
        long as the result is a beam."""
        return None

    def _sequence(self):
        """Overload this method to read the image file however you like so
        long as the result is a sequence."""
        return None

    def get_static_mask(self):
        """Overload this method to override the static mask."""
        return None

    def get_spectrum(self):
        """Overload this method to read the image file however you like so
        long as the result is a spectrum
        """
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

    @classmethod
    def open_file(cls, filename, mode="rb"):
        """Open file for reading, decompressing silently if necessary,
        caching transparently if possible."""

        if not Path(filename).is_file():
            raise Sorry(f"File {filename} not found")

        if filename.endswith(".bz2"):
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
