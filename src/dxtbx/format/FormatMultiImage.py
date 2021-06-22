import functools
import os

from dxtbx.format.Format import Format, abstract
from dxtbx.format.image import ImageBool
from dxtbx.imageset import (
    ImageSequence,
    ImageSet,
    ImageSetData,
    ImageSetLazy,
    ImageSetType,
    TOFImageSet,
    TOFImageSetData,
)
from dxtbx.model import MultiAxisGoniometer


def _add_static_mask_to_iset(format_instance: Format, iset: ImageSet) -> None:
    """Combine any static mask from a Format with the ImageSet's mask"""
    if format_instance is not None:
        static_mask = format_instance.get_static_mask()
        if static_mask is not None:
            if not iset.external_lookup.mask.data.empty():
                for m1, m2 in zip(static_mask, iset.external_lookup.mask.data):
                    m1 &= m2.data()
                iset.external_lookup.mask.data = ImageBool(static_mask)
            else:
                iset.external_lookup.mask.data = ImageBool(static_mask)


class Reader:
    def __init__(self, format_class, filenames, num_images=None, **kwargs):
        self.kwargs = kwargs
        self.format_class = format_class
        assert len(filenames) == 1
        self._filename = filenames[0]
        if num_images is None:
            format_instance = self.format_class.get_instance(
                self._filename, **self.kwargs
            )
            self._num_images = format_instance.get_num_images()
        else:
            self._num_images = num_images

    def nullify_format_instance(self):
        self.format_class._current_instance_ = None
        self.format_class._current_filename_ = None

    def read(self, index):
        format_instance = self.format_class.get_instance(self._filename, **self.kwargs)
        return format_instance.get_raw_data(index)

    def paths(self):
        return [self._filename]

    def __len__(self):
        return self._num_images

    def copy(self, filenames, indices=None):
        return Reader(self.format_class, filenames, indices)

    def identifiers(self):
        return ["%s-%d" % (self._filename, index) for index in range(len(self))]

    def is_single_file_reader(self):
        return True

    def master_path(self):
        return self._filename


@abstract
class FormatMultiImage(Format):
    def __init__(self, **kwargs):
        pass

    def get_num_images(self):
        raise NotImplementedError

    def get_goniometer(self, index=None):
        return self._goniometer_instance

    def get_detector(self, index=None):
        return self._detector_instance

    def get_beam(self, index=None):
        return self._beam_instance

    def get_spectrum(self, index=None):
        raise NotImplementedError

    def get_scan(self, index=None):
        return self._scan_instance

    def get_raw_data(self, index=None):
        raise NotImplementedError

    def get_detectorbase(self, index=None):
        raise NotImplementedError

    @classmethod
    def get_reader(cls):
        """
        Return a reader class
        """
        return functools.partial(Reader, cls)

    def get_masker(self, goniometer=None):
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
        cls,
        filenames,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
        imageset_type=None,
        format_kwargs=None,
        **kwargs,
    ):
        """
        Factory method to create an imageset
        """

        def create_imagesequence(
            cls,
            filenames,
            beam,
            detector,
            goniometer,
            scan,
            single_file_indices,
            format_instance,
            format_kwargs,
        ):

            # Check indices are sequential
            if single_file_indices is not None:
                assert all(
                    i + 1 == j
                    for i, j in zip(single_file_indices[:-1], single_file_indices[1:])
                )
                num_images = len(single_file_indices)
            else:
                num_images = format_instance.get_num_images()

            reader = get_reader(cls, filenames, num_images, **format_kwargs)

            # Check the scan makes sense - we must want to use <= total images
            if scan is not None:
                assert scan.get_num_images() <= num_images

            if format_instance is None:
                vendor = ""
            else:
                vendor = format_instance.get_vendortype()

            # If any are None then read from format
            if beam is None:
                beam = format_instance.get_beam()
            if detector is None:
                detector = format_instance.get_detector()
            if goniometer is None:
                goniometer = format_instance.get_goniometer()
            if scan is None:
                scan = format_instance.get_scan()

            # Create the masker
            if format_instance is not None:
                masker = format_instance.get_masker(goniometer=goniometer)
            else:
                masker = None

            isetdata = ImageSetData(
                reader=reader,
                masker=masker,
                vendor=vendor,
                params=format_kwargs,
                format=cls,
                template=filenames[0],
            )

            # Create the sequence
            iset = ImageSequence(
                isetdata,
                beam=beam,
                detector=detector,
                goniometer=goniometer,
                scan=scan,
                indices=single_file_indices,
            )

            _add_static_mask_to_iset(format_instance, iset)

            return iset

        def create_imageset(
            cls,
            filenames,
            beam,
            detector,
            single_file_indices,
            format_instance,
            format_kwargs,
        ):

            num_images = get_num_images(single_file_indices, format_instance)
            reader = get_reader(cls, filenames, num_images, **format_kwargs)

            if format_instance is None:
                vendor = ""
            else:
                vendor = format_instance.get_vendortype()

            iset = ImageSet(
                ImageSetData(
                    reader=reader,
                    masker=None,
                    vendor=vendor,
                    params=format_kwargs,
                    format=cls,
                ),
                indices=single_file_indices,
            )

            if single_file_indices is None:
                single_file_indices = range(format_instance.get_num_images())

                # If any are None then read from format
                num_images = format_instance.get_num_images()
                beam = [None] * num_images
                detector = [None] * num_images
                goniometer = [None] * num_images
                scan = [None] * num_images
                for i in single_file_indices:
                    beam[i] = format_instance.get_beam(i)
                    detector[i] = format_instance.get_detector(i)
                    goniometer[i] = format_instance.get_goniometer(i)
                    scan[i] = format_instance.get_scan(i)

            # Set the list of models
            for i, index in enumerate(single_file_indices):
                iset.set_beam(beam[index], i)
                iset.set_detector(detector[index], i)
                iset.set_goniometer(goniometer[index], i)
                iset.set_scan(scan[index], i)

            _add_static_mask_to_iset(format_instance, iset)

            return iset

        def create_imageset_lazy(cls, filenames, single_file_indices, format_kwargs):

            format_instance = None
            num_images = get_num_images(single_file_indices, format_instance)
            reader = get_reader(cls, filenames, num_images, **format_kwargs)
            vendor = ""

            isd = ImageSetData(
                reader=reader,
                masker=None,
                vendor=vendor,
                params=format_kwargs,
                format=cls,
            )

            iset = ImageSetLazy(
                isd,
                indices=single_file_indices,
            )
            _add_static_mask_to_iset(format_instance, iset)
            return iset

        def create_tof_imageset(
            cls,
            filenames,
            beam,
            detector,
            single_file_indices,
            format_instance,
            format_kwargs,
        ):
            num_images = get_num_images(single_file_indices, format_instance)
            tof_in_seconds = format_instance.get_tof_in_seconds()
            reader = get_reader(cls, filenames, num_images, **format_kwargs)
            vendor = format_instance.get_vendortype()

            if "tof_indices" in format_kwargs:
                if single_file_indices is not None:
                    iset = TOFImageSet(
                        TOFImageSetData(
                            reader=reader,
                            masker=None,
                            vendor=vendor,
                            params=format_kwargs,
                            format=cls,
                        ),
                        tof_in_seconds=tof_in_seconds,
                        tof_indices=format_kwargs["tof_indices"],
                        indices=single_file_indices,
                    )
                else:
                    iset = TOFImageSet(
                        TOFImageSetData(
                            reader=reader,
                            masker=None,
                            vendor=vendor,
                            params=format_kwargs,
                            format=cls,
                        ),
                        tof_in_seconds=tof_in_seconds,
                        tof_indices=format_kwargs["tof_indices"],
                    )
            else:
                iset = TOFImageSet(
                    TOFImageSetData(
                        reader=reader,
                        masker=None,
                        vendor=vendor,
                        params=format_kwargs,
                        format=cls,
                    ),
                    tof_in_seconds=tof_in_seconds,
                )

            if single_file_indices is None:
                single_file_indices = range(format_instance.get_num_images())

                # If any are None then read from format
                num_images = format_instance.get_num_images()
                beam = [None] * num_images
                detector = [None] * num_images
                goniometer = [None] * num_images
                scan = [None] * num_images
                for i in single_file_indices:
                    beam[i] = format_instance.get_beam(i)
                    detector[i] = format_instance.get_detector(i)

            # Set the list of models
            for i, index in enumerate(single_file_indices):
                iset.set_beam(beam[index], i)
                iset.set_detector(detector[index], i)
                iset.set_goniometer(goniometer[index], i)
                iset.set_scan(scan[index], i)

            _add_static_mask_to_iset(format_instance, iset)
            return iset

        def process_filenames(filenames):

            # Ensure filenames is list of length 1
            if isinstance(filenames, str):
                filenames = [filenames]
            elif len(filenames) > 1:
                assert len(set(filenames)) == 1
                filenames = filenames[0:1]

            # Make filenames absolute
            filenames = [os.path.abspath(x) for x in filenames]

            return filenames

        def process_format_kwargs(format_kwargs):
            if format_kwargs is None:
                format_kwargs = {}
            return format_kwargs

        def get_single_file_indices(**kwargs):
            if "single_file_indices" in kwargs:
                return kwargs["single_file_indices"]
            else:
                return None

        def sanity_check_params(cls, single_file_indices):
            if cls.get_num_images == FormatMultiImage.get_num_images:
                assert single_file_indices is not None
                assert min(single_file_indices) >= 0

        def get_num_images(single_file_indices, format_instance):
            if single_file_indices is not None:
                return max(single_file_indices) + 1
            elif format_instance is not None:
                return format_instance.get_num_images()
            return None

        def get_reader(cls, filenames, num_images, **format_kwargs):
            return cls.get_reader()(filenames, num_images=num_images, **format_kwargs)

        def format_instance_required(
            beam, detector, single_file_indices, imageset_type
        ):
            if beam is None or detector is None:
                return True
            if single_file_indices is None:
                return True
            if imageset_type == ImageSetType.TOFImageSet:
                return True
            return False

        # Process input
        single_file_indices = get_single_file_indices(**kwargs)
        sanity_check_params(cls, single_file_indices)
        filenames = process_filenames(filenames)
        format_kwargs = process_format_kwargs(format_kwargs)

        # ImageSetLazy is unique in not needing model instances, and so is checked first
        if imageset_type == ImageSetType.ImageSetLazy:
            return create_imageset_lazy(
                cls, filenames, single_file_indices, format_kwargs
            )

        # A format instance is only created if there is missing information in params
        if format_instance_required(beam, detector, single_file_indices, imageset_type):
            cls._current_filename_ = None
            cls._current_instance_ = None
            format_instance = cls.get_instance(filenames[0], **format_kwargs)
        else:
            format_instance = None

        # Attempt to identify imageset type from models if type is not given
        if imageset_type is None:
            imageset_type = Format.identify_imageset_type(
                scan, goniometer, beam, format_instance
            )

        if imageset_type == ImageSetType.ImageSet:
            return create_imageset(
                cls,
                filenames,
                beam,
                detector,
                single_file_indices,
                format_instance,
                format_kwargs,
            )
        elif imageset_type == ImageSetType.ImageSequence:
            return create_imagesequence(
                cls,
                filenames,
                beam,
                detector,
                goniometer,
                scan,
                single_file_indices,
                format_instance,
                format_kwargs,
            )
        elif imageset_type == ImageSetType.TOFImageSet:
            return create_tof_imageset(
                cls,
                filenames,
                beam,
                detector,
                single_file_indices,
                format_instance,
                format_kwargs,
            )
        else:
            raise NotImplementedError(f"Imageset_type {imageset_type} not implemented")
