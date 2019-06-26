from __future__ import absolute_import, division, print_function

from builtins import range
import boost.python
import dxtbx.format.image  # noqa: F401, import dependency for unpickling
import dxtbx.format.Registry

ext = boost.python.import_ext("dxtbx_ext")
from dxtbx_imageset_ext import *


class MemReader(object):
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


class MemMasker(object):
    def __init__(self, images):
        self._images = images

    def get(self, index, goniometer=None):
        format_instance = self._images[index]
        return format_instance.get_mask(goniometer=goniometer)

    def paths(self):
        return ["" for im in self._images]

    def identifiers(self):
        return self.paths()

    def __len__(self):
        return len(self._images)


class ImageSetAux(boost.python.injector, ImageSet):
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
        """ Iterate over the array indices and read each image in turn. """
        for i in range(len(self)):
            yield self[i]

    def get_vendortype(self, index):
        """ Get the vendor information. """
        return self.data().get_vendor()

    def get_format_class(self):
        """ Get format class name """
        return self.data().get_format_class()

    def params(self):
        """ Get the parameters """
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

    def get_detector(self, index=None):
        if index is None:
            index = 0
        detector = super(ImageSetLazy, self).get_detector(index)
        if detector is None:
            # If check_format=False was used, then _current_instance_ will not be set, so assume a None is correct
            format_class = self.get_format_class()
            if (
                hasattr(format_class, "_current_instance_")
                and format_class._current_instance_ is not None
            ):
                format_instance = format_class._current_instance_
                detector = format_instance.get_detector(self.indices()[index])
                self.set_detector(detector, index)
        return detector

    def get_beam(self, index=None):
        if index is None:
            index = 0
        beam = super(ImageSetLazy, self).get_beam(index)
        if beam is None:
            # If check_format=False was used, then _current_instance_ will not be set, so assume a None is correct
            format_class = self.get_format_class()
            if (
                hasattr(format_class, "_current_instance_")
                and format_class._current_instance_ is not None
            ):
                format_instance = format_class._current_instance_
                beam = format_instance.get_beam(self.indices()[index])
                self.set_beam(beam, index)
        return beam

    def get_goniometer(self, index=None):
        if index is None:
            index = 0
        goniometer = super(ImageSetLazy, self).get_goniometer(index)
        if goniometer is None:
            # If check_format=False was used, then _current_instance_ will not be set, so assume a None is correct
            format_class = self.get_format_class()
            if (
                hasattr(format_class, "_current_instance_")
                and format_class._current_instance_ is not None
            ):
                format_instance = format_class._current_instance_
                goniometer = format_instance.get_goniometer(self.indices()[index])
                self.set_goniometer(goniometer, index)
        return goniometer

    def get_scan(self, index=None):
        if index is None:
            index = 0
        scan = super(ImageSetLazy, self).get_scan(index)
        if scan is None:
            # If check_format=False was used, then _current_instance_ will not be set, so assume a None is correct
            format_class = self.get_format_class()
            if (
                hasattr(format_class, "_current_instance_")
                and format_class._current_instance_ is not None
            ):
                format_instance = format_class._current_instance_
                scan = format_instance.get_scan(self.indices()[index])
                self.set_scan(scan, index)
        return scan

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
        return super(ImageSetLazy, self).__getitem__(item)

    def get_corrected_data(self, index):
        self._load_models(index)
        return super(ImageSetLazy, self).get_corrected_data(index)

    def get_gain(self, index):
        self._load_models(index)
        return super(ImageSetLazy, self).get_gain(index)


class ImageSweepAux(boost.python.injector, ImageSweep):
    def __getitem__(self, item):
        """Get an item from the sweep stream.

        If the item is an index, read and return the image at the given index.
        Otherwise, if the item is a slice, then create a new Sweep object
        with the given number of array indices from the slice.

        Params:
            item The index or slice

        Returns:
            An image or new Sweep object

        """
        if isinstance(item, slice):
            start = item.start or 0
            stop = item.stop or len(self)
            if item.step is not None:
                raise IndexError("Sweeps must be sequential")
            return self.partial_set(start, stop)
        else:
            return self.get_corrected_data(item)

    def get_template(self):
        """ Return the template """
        return self.data().get_template()


def _analyse_files(filenames):
    """Group images by filename into image sets.

    Params:
        filenames The list of filenames

    Returns:
        A list of (template, [indices], is_sweep)

    """
    from dxtbx.sweep_filenames import group_files_by_imageset

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

    def _is_imageset_a_sweep(template, indices):
        """Return True/False if the imageset is a sweep or not.

        Where more than 1 image that follow sequential numbers are given
        the images are catagorised as belonging to a sweep, otherwise they
        belong to an image set.

        """
        if len(indices) <= 1:
            return False
        indices = sorted(indices)
        return _indices_sequential_ge_zero(indices)

    # Label each group as either an imageset or a sweep.
    file_groups = []
    for template, indices in filelist_per_imageset.items():

        # Check if this imageset is a sweep
        is_sweep = _is_imageset_a_sweep(template, indices)

        # Append the items to the group list
        file_groups.append((template, indices, is_sweep))

    # Return the groups of files
    return file_groups


# FIXME Lots of duplication in this class, need to tidy up
class ImageSetFactory(object):
    """ Factory to create imagesets and sweeps. """

    @staticmethod
    def new(filenames, check_headers=False, ignore_unknown=False):
        """Create an imageset or sweep

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
        # and return as a list of imagesets. N.B sweeps and image sets are
        # returned in the same list.
        imagesetlist = []
        for filelist in filelist_per_imageset:
            try:
                if filelist[2] is True:
                    iset = ImageSetFactory._create_sweep(filelist, check_headers)
                else:
                    iset = ImageSetFactory._create_imageset(filelist, check_headers)
                imagesetlist.append(iset)
            except Exception:
                if not ignore_unknown:
                    raise

        # Return the imageset list
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
        """Create a new sweep from a template.

        Params:
            template The template argument
            image_range The image range
            check_headers Check the headers to ensure all images are valid

        Returns:
            A list of sweeps

        """
        from dxtbx.sweep_filenames import template_image_range
        from dxtbx.format.Format import Format

        if not check_format:
            assert not check_headers

        # Check the template is valid
        if template.count("#") == 0:
            if "master" not in template:
                raise ValueError("Invalid template")
            filenames = [template]
        else:

            # Get the template format
            pfx = template.split("#")[0]
            sfx = template.split("#")[-1]
            template_format = "%s%%0%dd%s" % (pfx, template.count("#"), sfx)

            # Get the template image range
            if image_range is None:
                image_range = template_image_range(template)

            # Set the image range
            array_range = range(image_range[0] - 1, image_range[1])

            # Create the sweep file list
            filenames = [template_format % (i + 1) for i in array_range]

        # Get the format class
        if check_format:
            format_class = dxtbx.format.Registry.get_format_class_for_file(filenames[0])
        else:
            format_class = Format

        # Create the sweep object
        sweep = format_class.get_imageset(
            filenames,
            template=template,
            as_sweep=True,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            check_format=check_format,
        )

        # Return the sweep
        return [sweep]

    @staticmethod
    def _create_imageset(filelist, check_headers):
        """Create an image set"""
        # Extract info from filelist
        template, indices, is_sweep = filelist

        # Get the template format
        count = template.count("#")
        if count > 0:
            pfx = template.split("#")[0]
            sfx = template.split("#")[-1]
            template_format = "%s%%0%dd%s" % (pfx, template.count("#"), sfx)
            filenames = [template_format % index for index in indices]
        else:
            filenames = [template]

        # Sort the filenames
        filenames = sorted(filenames)

        # Get the format object
        format_class = dxtbx.format.Registry.get_format_class_for_file(filenames[0])

        # Create and return the imageset
        return format_class.get_imageset(filenames, as_imageset=True)

    @staticmethod
    def _create_sweep(filelist, check_headers):
        """Create a sweep"""
        # Extract info from filelist
        template, indices, is_sweep = filelist

        # Get the template format
        count = template.count("#")
        if count > 0:
            pfx = template.split("#")[0]
            sfx = template.split("#")[-1]
            template_format = "%s%%0%dd%s" % (pfx, template.count("#"), sfx)
            filenames = [template_format % index for index in indices]
        else:
            filenames = [template]

        # Sort the filenames
        filenames = sorted(filenames)

        # Get the format object
        format_class = dxtbx.format.Registry.get_format_class_for_file(filenames[0])

        # Get the template format
        pfx = template.split("#")[0]
        sfx = template.split("#")[-1]
        template_format = "%s%%0%dd%s" % (pfx, template.count("#"), sfx)

        # Set the image range
        array_range = range(min(indices) - 1, max(indices))

        # Create the sweep file list
        filenames = [template_format % (i + 1) for i in array_range]

        sweep = format_class.get_imageset(filenames, template=template, as_sweep=True)

        # Return the sweep
        return sweep

    @staticmethod
    def make_imageset(
        filenames,
        format_class=None,
        check_format=True,
        single_file_indices=None,
        format_kwargs=None,
    ):
        """Create an image set"""
        from dxtbx.format.Format import Format
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
    def make_sweep(
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
        """Create a sweep"""
        from dxtbx.format.Format import Format

        indices = sorted(indices)

        # Get the template format
        count = template.count("#")
        if count > 0:
            pfx = template.split("#")[0]
            sfx = template.split("#")[-1]
            template_format = "%s%%0%dd%s" % (pfx, template.count("#"), sfx)
            filenames = [template_format % index for index in indices]
        else:
            template_format = None
            filenames = [template]

        # Sort the filenames
        filenames = sorted(filenames)

        # Set the image range
        array_range = (min(indices) - 1, max(indices))
        if scan is not None:
            assert array_range == scan.get_array_range()

        # Get the format object and reader
        if format_class is None:
            if check_format:
                format_class = dxtbx.format.Registry.get_format_class_for_file(
                    filenames[0]
                )
            else:
                format_class = Format

        sweep = format_class.get_imageset(
            filenames,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            format_kwargs=format_kwargs,
            template=template,
            as_sweep=True,
            check_format=check_format,
            single_file_indices=list(range(*array_range)),
        )

        # Return the sweep
        return sweep

    @staticmethod
    def imageset_from_anyset(imageset):
        """ Create a new ImageSet object from an imageset object. Converts ImageSweep to ImageSet. """
        from dxtbx.imageset import ImageSet, ImageSweep, ImageSetLazy

        if isinstance(imageset, ImageSetLazy):
            return ImageSetLazy(imageset.data(), imageset.indices())
        elif isinstance(imageset, ImageSweep) or isinstance(imageset, ImageSet):
            return ImageSet(imageset.data(), imageset.indices())
        else:
            assert False, "Unrecognized imageset type: %s" % str(type(imageset))


try:
    boost.python.inject_into(ImageSet)(ImageSetAux)
    boost.python.inject_into(ImageSweep)(ImageSweepAux)
except AttributeError:
    pass
