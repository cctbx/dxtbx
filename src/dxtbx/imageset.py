import boost_adaptbx.boost.python

import dxtbx.format.image  # noqa: F401, import dependency for unpickling
import dxtbx.format.Registry
from dxtbx.sequence_filenames import group_files_by_imageset, template_image_range
from dxtbx_imageset_ext import (
    ExternalLookup,
    ExternalLookupItemBool,
    ExternalLookupItemDouble,
    ImageGrid,
    ImageSequence,
    ImageSet,
    ImageSetData,
)

ext = boost_adaptbx.boost.python.import_ext("dxtbx_ext")

from typing import Iterable, List

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


def _expand_template(template: str, indices: Iterable[int]) -> List[str]:
    """Expand a template string to a list of filenames.

    Args:
        template: The template string, with a block of '#' to replace
        indices:  The numeric indices to insert
    """
    pfx = template.split("#")[0]
    sfx = template.split("#")[-1]
    count = template.count("#")
    return [f"{pfx}{index:0{count}}{sfx}" for index in indices]


class MemReader:
    """A reader for data already loaded in memory"""

    def __init__(self, images):
        self._images = images

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
            return self.partial_set(start, stop)
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
        return format_instance.get_spectrum(self.indices()[index])

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
        return [self.get_path(i) for i in range(len(self))]


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
            return ImageSetLazy(self.data(), indices=self.indices()[item])
        self._load_models(item)
        return super().__getitem__(item)

    def get_corrected_data(self, index):
        self._load_models(index)
        return super().get_corrected_data(index)

    def get_gain(self, index):
        self._load_models(index)
        return super().get_gain(index)


@boost_adaptbx.boost.python.inject_into(ImageSequence)
class _:
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
        if isinstance(item, slice):
            offset = self.get_scan().get_batch_offset()
            if item.step is not None:
                raise IndexError("Sequences must be sequential")

            # nasty workaround for https://github.com/dials/dials/issues/1153
            # slices with -1 in them are meaningful :-/ so grab the original
            # constructor arguments of the slice object.
            # item.start and item.stop may have been compromised at this point.
            if offset < 0:
                start, stop, step = item.__reduce__()[1]
                if start is None:
                    start = 0
                else:
                    start -= offset
                if stop is None:
                    stop = len(self)
                else:
                    stop -= offset
                return self.partial_set(start, stop)
            else:
                start = item.start or 0
                stop = item.stop or (len(self) + offset)
                return self.partial_set(start - offset, stop - offset)
        else:
            return self.get_corrected_data(item)

    def get_template(self):
        """Return the template"""
        return self.data().get_template()


def _analyse_files(filenames):
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

    def _is_imageset_a_sequence(template, indices):
        """Return True/False if the imageset is a sequence or not.

        Where more than 1 image that follow sequential numbers are given
        the images are catagorised as belonging to a sequence, otherwise they
        belong to an image set.

        """
        if len(indices) <= 1:
            return False
        indices = sorted(indices)
        return _indices_sequential_ge_zero(indices)

    # Label each group as either an imageset or a sequence.
    file_groups = []
    for template, indices in filelist_per_imageset.items():

        # Check if this imageset is a sequence
        is_sequence = _is_imageset_a_sequence(template, indices)

        # Append the items to the group list
        file_groups.append((template, indices, is_sequence))

    return file_groups


# FIXME Lots of duplication in this class, need to tidy up
class ImageSetFactory:
    """Factory to create imagesets and sequences."""

    @staticmethod
    def new(filenames, check_headers=False, ignore_unknown=False):
        """Create an imageset or sequence

        Params:
            filenames A list of filenames
            check_headers Check the headers to ensure all images are valid
            ignore_unknown Ignore unknown formats

        Returns:
            A list of imagesets

        """
        # Ensure we have enough images
        if isinstance(filenames, list):
            assert filenames
        elif isinstance(filenames, str):
            filenames = [filenames]
        else:
            raise RuntimeError("unknown argument passed to ImageSetFactory")

        # Analyse the filenames and group the images into imagesets.
        filelist_per_imageset = _analyse_files(filenames)

        # For each file list denoting an image set, create the imageset
        # and return as a list of imagesets. N.B sequences and image sets are
        # returned in the same list.
        imagesetlist = []
        for filelist in filelist_per_imageset:
            try:
                if filelist[2] is True:
                    iset = ImageSetFactory._create_sequence(filelist, check_headers)
                else:
                    iset = ImageSetFactory._create_imageset(filelist, check_headers)
                imagesetlist.append(iset)
            except Exception:
                if not ignore_unknown:
                    raise

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

        # Check the template is valid
        if "#" in template:
            # Get the template image range
            if image_range is None:
                image_range = template_image_range(template)

            # Set the image range
            indices = range(image_range[0], image_range[1] + 1)
            filenames = _expand_template(template, indices)
        else:
            if "master" not in template:
                raise ValueError("Invalid template")
            filenames = [template]

        # Import here as Format and Imageset have cyclic dependencies
        from dxtbx.format.Format import Format

        # Get the format class
        if check_format:
            format_class = dxtbx.format.Registry.get_format_class_for_file(filenames[0])
        else:
            format_class = Format

        # Create the sequence object
        sequence = format_class.get_imageset(
            filenames,
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
    def _create_imageset(filelist, check_headers):
        """Create an image set"""
        # Extract info from filelist
        template, indices, is_sequence = filelist

        # Get the template format
        if "#" in template:
            filenames = sorted(_expand_template(template, indices))
        else:
            filenames = [template]

        # Get the format object
        format_class = dxtbx.format.Registry.get_format_class_for_file(filenames[0])

        # Create and return the imageset
        return format_class.get_imageset(filenames, as_imageset=True)

    @staticmethod
    def _create_sequence(filelist, check_headers):
        """Create a sequence"""
        template, indices, is_sequence = filelist

        # Expand the template if necessary
        if "#" in template:
            filenames = sorted(_expand_template(template, indices))
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
        template,
        indices,
        format_class=None,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
        check_format=True,
        format_kwargs=None,
    ):
        """Create a sequence"""
        indices = sorted(indices)

        # Get the template format
        if "#" in template:
            filenames = sorted(_expand_template(template, indices))
        else:
            filenames = [template]

        # Set the image range
        array_range = (min(indices) - 1, max(indices))
        if scan is not None:
            assert array_range == scan.get_array_range()
            scan.set_batch_offset(array_range[0])

        # Get the format object and reader
        if format_class is None:
            # Import here as Format and Imageset have cyclic dependencies
            from dxtbx.format.Format import Format

            if check_format:
                format_class = dxtbx.format.Registry.get_format_class_for_file(
                    filenames[0]
                )
            else:
                format_class = Format

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
