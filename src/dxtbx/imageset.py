from __future__ import annotations

import typing
import warnings
from collections.abc import Iterable
from typing import Sequence

import natsort

import boost_adaptbx.boost.python

import dxtbx.format.image  # noqa: F401, import dependency for unpickling
import dxtbx.format.Registry
from dxtbx.sequence_filenames import group_files_by_imageset, template_image_range

from . import model

try:
    from .dxtbx_imageset_ext import (
        ExternalLookup,
        ExternalLookupItemBool,
        ExternalLookupItemDouble,
        ImageGrid,
        ImageSequence,
        ImageSet,
        ImageSetData,
    )
except ModuleNotFoundError:
    from dxtbx_imageset_ext import (  # type: ignore
        ExternalLookup,
        ExternalLookupItemBool,
        ExternalLookupItemDouble,
        ImageGrid,
        ImageSequence,
        ImageSet,
        ImageSetData,
    )

if typing.TYPE_CHECKING:
    # Here because circular dependency
    from dxtbx.format.Format import Format

ext = boost_adaptbx.boost.python.import_ext("dxtbx_ext")

__all__ = (
    "ExternalLookup",
    "ExternalLookupItemBool",
    "ExternalLookupItemDouble",
    "ImageGrid",
    "ImageSet",
    "ImageSetData",
    "ImageSetFactory",
    "ImageSetLazy",
    "ImageSequence",
    "MemReader",
)


def _expand_template_to_sorted_filenames(
    template: str, indices: Iterable[int]
) -> list[str]:
    """Expand a template string to a list of filenames.

    Args:
        template: The template string, with a block of '#' to replace
        indices:  The numeric indices to insert
    """
    pfx = template.split("#")[0]
    sfx = template.split("#")[-1]
    count = template.count("#")
    if count == 1:
        # Special handling for a template with a single "#", which does not
        # assume a zero-padded index.
        filenames = [f"{pfx}{index}{sfx}" for index in indices]
    else:
        filenames = [f"{pfx}{index:0{count}}{sfx}" for index in indices]
    return natsort.natsorted(filenames)


class MemReader:
    """A reader for data already loaded in memory"""

    def __init__(self, images):
        self._images = images

    def copy(self, paths):
        """
        Experimental implementation where a copy of the reader also copies all
        the data
        """
        return MemReader(self._images)

    def paths(self):
        return ["" for im in self._images]

    def identifiers(self):
        return self.paths()

    def __len__(self):
        return len(self._images)

    def read(self, index):
        format_instance = self._images[index]
        return format_instance.get_raw_data()

    @staticmethod
    def is_single_file_reader():
        return False

    @staticmethod
    def master_path():
        return ""


@boost_adaptbx.boost.python.inject_into(ImageSet)
class _:
    """
    A class to inject additional methods into the imageset class
    """

    def __getitem__(self, item):
        """Get an item from the image set stream.

        If the item is an index, read and return the image at the given index.
        Otherwise, if the item is a slice, then create a new ImageSet object
        with the given number of array indices from the slice.

        Params:
            item The index or slice

        Returns:
            An image or new ImageSet object
        """
        if isinstance(item, slice):
            start = item.start or 0
            stop = item.stop or len(self)
            if item.step is not None and item.step != 1:
                raise IndexError("Step must be 1")
            if self.data().has_single_file_reader():
                reader = self.reader().copy(self.reader().paths(), stop - start)
            else:
                reader = self.reader().copy(self.reader().paths())
            return self.partial_set(reader, start, stop)
        else:
            return self.get_corrected_data(item)

    def __iter__(self):
        """Iterate over the array indices and read each image in turn."""
        for i in range(len(self)):
            yield self[i]

    def get_vendortype(self, index):
        """Get the vendor information."""
        return self.data().get_vendor()

    def get_format_class(self):
        """Get format class name"""
        return self.data().get_format_class()

    def get_spectrum(self, index):
        """Get the spectrum if available"""
        kwargs = self.params()
        if self.data().has_single_file_reader():
            format_instance = self.get_format_class().get_instance(
                self.data().get_master_path(), **kwargs
            )
        else:
            format_instance = self.get_format_class().get_instance(
                self.get_path(index), **kwargs
            )
        try:
            return format_instance.get_spectrum(self.indices()[index])
        except TypeError:
            return format_instance.get_spectrum()

    def params(self):
        """Get the parameters"""
        return self.data().get_params()

    def get_detectorbase(self, index):
        """
        A function to be injected into the imageset to get the detectorbase instance
        """
        kwargs = self.params()
        if self.data().has_single_file_reader():
            format_instance = self.get_format_class().get_instance(
                self.data().get_master_path(), **kwargs
            )
            return format_instance.get_detectorbase(self.indices()[index])
        else:
            format_instance = self.get_format_class().get_instance(
                self.get_path(index), **kwargs
            )
            return format_instance.get_detectorbase()

    def reader(self):
        """
        Return the reader
        """
        return self.data().reader()

    def masker(self):
        """
        Return the masker
        """
        return self.data().masker()

    def paths(self):
        """
        Return the list of paths
        """
        if self.data().has_single_file_reader():
            return [self.get_path(i) for i in range(len(self))]
        else:
            return [self.reader().paths()[i] for i in self.indices()]


class ImageSetLazy(ImageSet):
    """
    Lazy ImageSet class that doesn't necessitate setting the models ahead of time.
    Only when a particular model (like detector or beam) for an image is requested,
    it sets the model using the format class and then returns the model
    """

    def _get_item_from_parent_or_format(self, item_name, index):
        """
        Obtain an $item_name (eg. detector, beam, ...) of the given index from
        the parent class using get_detector, get_beam, ...
        If the parent class returns None then lookup the item using the format
        class (if defined) and store a local reference to the item using
        self.set_detector, set_beam, ..
        """
        if index is None:
            index = 0
        item = getattr(super(), "get_" + item_name)(index)
        if item is None:
            # If check_format=False was used, then _current_instance_ will not be set, so assume a None is correct
            format_class = self.get_format_class()
            if (
                hasattr(format_class, "_current_instance_")
                and format_class._current_instance_ is not None
            ):
                format_instance = format_class._current_instance_
                getter_function = getattr(format_instance, "get_" + item_name)
                item = getter_function(self.indices()[index])
                setter_function = getattr(self, "set_" + item_name)
                setter_function(item, index)
        return item

    def get_detector(self, index=None):
        return self._get_item_from_parent_or_format("detector", index)

    def get_beam(self, index=None):
        return self._get_item_from_parent_or_format("beam", index)

    def get_mask(self, index=None):
        """
        ImageSet::get_mask internally dereferences a pointer to the _detector
        member of ImageSetData, so we ensure the detector gets populated first.
        """
        if getattr(super(), "get_detector")(index) is None:
            self._load_models(index)
        return self._get_item_from_parent_or_format("mask", index)

    def get_goniometer(self, index=None):
        return self._get_item_from_parent_or_format("goniometer", index)

    def get_scan(self, index=None):
        return self._get_item_from_parent_or_format("scan", index)

    def _load_models(self, index):
        if index is None:
            index = 0
        # Sets the list for detector, beam etc before being accessed by functions in imageset.h
        self.get_detector(index)
        self.get_beam(index)
        self.get_goniometer(index)
        self.get_scan(index)

    def __getitem__(self, item):
        if isinstance(item, slice):
            start = item.start or 0
            stop = item.stop or len(self)
            if item.step is not None and item.step != 1:
                raise IndexError("Step must be 1")
            if self.data().has_single_file_reader():
                reader = self.reader().copy(self.reader().paths(), stop - start)
            else:
                reader = self.reader().copy(self.reader().paths())
            return ImageSetLazy(
                self.data().partial_data(reader, start, stop),
                indices=self.indices()[item],
            )
        self._load_models(item)
        return super().__getitem__(item)

    def get_corrected_data(self, index):
        self._load_models(index)
        return super().get_corrected_data(index)

    def get_gain(self, index):
        self._load_models(index)
        return super().get_gain(index)


@boost_adaptbx.boost.python.inject_into(ImageSequence)
class _imagesequence:
    def __getitem__(self, item):
        """Get an item from the sequence stream.

        If the item is an index, read and return the image at the given index.
        Otherwise, if the item is a slice, then create a new Sequence object
        with the given number of array indices from the slice.

        Params:
            item The index or slice

        Returns:
            An image or new Sequence object

        """
        if not isinstance(item, slice):
            return self.get_corrected_data(item)
        else:
            if item.step is not None:
                raise IndexError("Sequences must be sequential")

            start = item.start or 0
            stop = item.stop or len(self)

            if self.data().has_single_file_reader():
                reader = self.reader().copy(self.reader().paths(), stop - start)
            else:
                reader = self.reader().copy(self.reader().paths())
            return self.partial_set(reader, start, stop)

    def get_template(self):
        """Return the template"""
        return self.data().get_template()


def _analyse_files(filenames: list[str]) -> list[tuple[str, list[int | None], bool]]:
    """Group images by filename into image sets.

    Params:
        filenames The list of filenames

    Returns:
        A list of (template, [indices], is_sequence)

    """
    # Analyse filenames to figure out how many imagesets we have
    filelist_per_imageset = group_files_by_imageset(filenames)

    def _indices_sequential_ge_zero(indices):
        """Determine if indices are sequential."""
        prev = indices[0]
        if prev < 0:
            return False
        for curr in indices[1:]:
            if curr != prev + 1:
                return False
            prev = curr

        return True

    def _is_imageset_a_sequence(indices: Sequence[int | None]) -> bool:
        """
        Check if a set of indices is sequential or not.

        Where more than 1 image that follow sequential numbers are given,
        the images are categorised as belonging to a sequence; otherwise they
        belong to an image set.

        A single index on it's own is never a sequence.
        """
        if len(indices) <= 1:
            return False
        indices = sorted(indices)
        return _indices_sequential_ge_zero(indices)

    # Label each group as either an imageset or a sequence.
    file_groups = []
    for template, indices in filelist_per_imageset.items():
        # Check if this imageset is a sequence
        is_sequence = _is_imageset_a_sequence(indices)

        # Append the items to the group list
        file_groups.append((template, indices, is_sequence))

    return file_groups


# FIXME Lots of duplication in this class, need to tidy up
class ImageSetFactory:
    """Factory to create imagesets and sequences."""

    @staticmethod
    def new(
        filenames: list[str] | str,
        *args,
        **kwargs,
    ) -> list[ImageSet | ImageSequence]:
        """
        Create a list of imageset and/or sequences.

        Args:
            filenames: A list of filename templates. These will be expanded, and
                any that appear to be a sequence will be loaded into ImageSequences.
            check_headers: [OBSOLETE] Did nothing.
            ignore_unknown: [OBSOLETE] Blanket ignored all errors when processing.
        """
        if args or kwargs:
            # We used to have two parameters here - that did nothing, or hid errors.
            # These were always optional, so warn the user if they are setting it.
            warnings.warn(
                "check_headers and ignore_unknown arguments to ImageSetFactory::new are obsolete",
                DeprecationWarning,
                stacklevel=2,
            )
        # Ensure we have enough images
        if isinstance(filenames, list):
            assert filenames
        elif isinstance(filenames, str):
            filenames = [filenames]
        else:
            raise RuntimeError("unknown argument passed to ImageSetFactory")

        # Analyse the filenames and group the images into imagesets.
        #
        # For each file list denoting an image set, create the imageset
        # and return as a list of imagesets. N.B sequences and image sets are
        # returned in the same list.
        imagesetlist = []
        for template, indices, is_sequence in _analyse_files(filenames):
            if is_sequence:
                iset = ImageSetFactory._create_sequence(template, indices)
            else:
                iset = ImageSetFactory._create_imageset(template, indices)
            imagesetlist.append(iset)

        return imagesetlist

    @staticmethod
    def from_template(
        template,
        image_range=None,
        check_headers=False,
        check_format=True,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
    ):
        """Create a new sequence from a template.

        Params:
            template The template argument
            image_range The image range
            check_headers Check the headers to ensure all images are valid

        Returns:
            A list of sequences
        """
        if not check_format:
            assert not check_headers

        # Import here as Format and Imageset have cyclic dependencies
        from dxtbx.format.Format import Format
        from dxtbx.format.FormatMultiImage import FormatMultiImage

        # Check the template is valid
        if "#" in template:
            # Get the template image range
            if image_range is None:
                image_range = template_image_range(template)

            # Set the image range
            indices = range(image_range[0], image_range[1] + 1)
            filenames = _expand_template_to_sorted_filenames(template, indices)
            if check_format:
                format_class = dxtbx.format.Registry.get_format_class_for_file(
                    filenames[0]
                )
            else:
                format_class = Format
        else:
            # Note, this assumes image stacks can only be written by dectris detectors,
            # but doesn't account for other imagestack formats like MRC.
            if "master" not in template:
                raise ValueError("Invalid template")
            filenames = [template]
            indices = range(image_range[0], image_range[1] + 1)
            if check_format:
                format_class = dxtbx.format.Registry.get_format_class_for_file(
                    filenames[0]
                )
            else:
                format_class = FormatMultiImage

        # Create the sequence object
        sequence = format_class.get_imageset(
            filenames,
            single_file_indices=indices,
            template=template,
            as_sequence=True,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            check_format=check_format,
        )

        return [sequence]

    @staticmethod
    def _create_imageset(template: str, indices: list[int | None]):
        """Create an image set"""
        # Get the template format
        if "#" in template:
            filenames = _expand_template_to_sorted_filenames(template, indices)
        else:
            filenames = [template]

        # Get the format object
        format_class = dxtbx.format.Registry.get_format_class_for_file(filenames[0])

        # Create and return the imageset
        return format_class.get_imageset(filenames, as_imageset=True)

    @staticmethod
    def _create_sequence(template: str, indices: list[int | None]) -> ImageSequence:
        """Create a sequence"""
        # Expand the template if necessary
        if "#" in template:
            assert not any(x is None for x in indices), (
                "Only accept None-indices for non-template filenames"
            )
            filenames = _expand_template_to_sorted_filenames(template, indices)
        else:
            filenames = [template]

        # Get the format object
        format_class = dxtbx.format.Registry.get_format_class_for_file(filenames[0])

        return format_class.get_imageset(filenames, template=template, as_sequence=True)

    @staticmethod
    def make_imageset(
        filenames,
        format_class=None,
        check_format=True,
        single_file_indices=None,
        format_kwargs=None,
    ):
        """Create an image set"""
        # Import here as Format and Imageset have cyclic dependencies
        from dxtbx.format.Format import Format

        # So does FormatMultiImage
        from dxtbx.format.FormatMultiImage import FormatMultiImage

        # Get the format object
        if format_class is None:
            if check_format:
                format_class = dxtbx.format.Registry.get_format_class_for_file(
                    filenames[0]
                )
            else:
                if single_file_indices is None or len(single_file_indices) == 0:
                    format_class = Format
                else:
                    format_class = FormatMultiImage

        return format_class.get_imageset(
            filenames,
            single_file_indices=single_file_indices,
            as_imageset=True,
            format_kwargs=format_kwargs,
            check_format=check_format,
        )

    @staticmethod
    def make_sequence(
        template: str,
        indices: Sequence[int],
        format_class: Format | None = None,
        beam: model.Beam | None = None,
        detector: model.Detector | None = None,
        goniometer: model.Goniometer | None = None,
        scan: model.Scan | None = None,
        check_format: bool = True,
        format_kwargs: dict | None = None,
    ) -> ImageSequence:
        """
        Create a sequence

        Args:
            format_class: The Format class to use. If unspecified, then
                this will be either automatically detected via dxtbx (if
                `check_format` is False), or default to Format or
                FormatMultiImage depending on whether the template looks
                like a multiimage template or not.
        """
        indices = sorted(indices)

        # Import here as Format and Imageset have cyclic dependencies
        from dxtbx.format.Format import Format
        from dxtbx.format.FormatMultiImage import FormatMultiImage

        if "#" in template:
            filenames = _expand_template_to_sorted_filenames(template, indices)
            default_format_class = Format
        else:
            filenames = [template]
            default_format_class = FormatMultiImage

        if format_class is None:
            if check_format:
                format_class = dxtbx.format.Registry.get_format_class_for_file(
                    filenames[0]
                )
            else:
                format_class = default_format_class
        assert format_class is not None

        # Set the image range
        array_range = (min(indices) - 1, max(indices))
        if scan is not None:
            assert array_range == scan.get_array_range()
            scan.set_batch_offset(array_range[0])

        return format_class.get_imageset(
            filenames,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            format_kwargs=format_kwargs,
            template=template,
            as_sequence=True,
            check_format=check_format,
            single_file_indices=list(range(*array_range)),
        )

    @staticmethod
    def imageset_from_anyset(imageset):
        """Create a new ImageSet object from an imageset object. Converts ImageSequence to ImageSet."""
        if isinstance(imageset, ImageSetLazy):
            return ImageSetLazy(imageset.data(), imageset.indices())
        elif isinstance(imageset, ImageSequence) or isinstance(imageset, ImageSet):
            return ImageSet(imageset.data(), imageset.indices())
        else:
            raise ValueError("Unrecognized imageset type: %s" % str(type(imageset)))
