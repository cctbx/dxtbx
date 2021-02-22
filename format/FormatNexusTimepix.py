"""Format class and assorted helper functions for the LATRD 'Tristan' detector."""

from typing import Sequence, Union

import h5py
import numpy as np

from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.model import Detector

Vector = Union[int, Sequence[int]]


def region_size(object_size: Vector, layout: Vector, stride: Vector) -> Vector:
    """Derive the total size of a grid of objects from their size, layout and stride."""
    layout = np.full((2,), layout)
    return object_size + (layout - 1) * stride


def object_layout(region_size: Vector, object_size: Vector, stride: Vector) -> Vector:
    """Derive the layout of a grid of objects from the size of the bounding region."""
    region_size = np.full((2,), region_size)
    return 1 + (region_size - object_size) / stride


# Each LATRD module consists of an arrangement of 8 × 2 Timepix3 chips.
chip_layout = (8, 2)  # Number of chips.
# Each Timepix3 chip is composed of 256 × 256 pixels, each 55µm × 55µm in size.
chip_size = 256  # Number of pixels.
# The edge pixels are actually slightly elongated, so there is effectively 3px-worth
# of space between adjacent chips in a LATRD module.
chip_stride = 259  # Number of pixels.

# From this, we can find the overall module size, including inter-chip spacing.
module_size = region_size(chip_size, chip_layout, chip_stride)  # Number of pixels.
# The LATRD comprises one or more such modules arranged in a grid,
# with some space between adjacent modules.
module_stride = (2114, 632)  # Number of pixels.


def gap_masks(object_size: Vector, layout: Vector, stride: Vector, inflate: Vector = 0):
    """
    Find the extents of the gaps between regularly spaced detector components.

    Area detectors often consist of a rectangular array of identical components
    (sensor chips, panels, etc.), with gaps in between that we wish to mask.  Given a
    region containing such an array of components, and knowing their size and
    spacing, the coordinates of the masked gap regions are computed.  Where the array
    consists of n_x × n_y components, there will be N = n_x + n_y - 2 such gap regions.

    Optionally, the masked region may be inflated by a specified number of pixels
    either side.  Sometimes this is necessary when the edge pixels behave in a
    different manner to the rest of the component, systematically over- or
    under-counting for example.  To illustrate, if adjacent components are separated
    by a gap three pixels wide and if a value inflate=1 has been specified,
    the masked region will be five pixels wide.

    Args:
        object_size:  The dimensions (size x, size y) in pixels of each component of
            the detector array.  If the components are square, the edge length can be
            passed as a scalar.
        layout:  The number of repeats (number x, number y) of detector components in
            the rectangular grid.  If there are the same number of components in each
            dimension, that number can be passed as a scalar.
        stride:  The spacing (stride x, stride y) of adjacent components in the grid.
            If the spacing is equal in each dimension, it can be passed as a scalar.
            The stride in each dimension will be the object size plus the width of the
            gap.
        inflate:  Width in pixels (number x, number y) to be added to each side of
            the masked gap region, in order to mask the edges of the components.  If
            the gaps in both the x and y directions are to be inflated by the same
            number of pixels, that number can be passed as a scalar.

    Returns:
        The parameters of the regions to be masked, a numpy.ndarray of shape (N, 2, 2),
        where N is the number of gaps between components.  Each of the N elements in
        the first dimension is a pair of two-vectors.  Each pair represents the
        bottom-left and top-right pixel, respectively, of a gap region to be masked.
        Each pixel is represented by a (x, y) vector of its position.
    """
    # Ensure all the arguments are cast to ndarrays of shape (2,).
    x_size, y_size = np.full((2,), object_size)
    x_layout, y_layout = np.full((2,), layout) - 1
    x_stride, y_stride = stride = np.full((2,), stride)
    x_inflate, y_inflate = np.full((2,), inflate)

    x_span, y_span = region_size(object_size, layout, stride)

    # The extents of the first gaps in x & y.
    x_gap_extent = np.array(
        (
            (x_size - x_inflate, 0),  # Bottom-left corner.
            (x_stride + x_inflate, y_span),  # Top-right corner.
        )
    )
    y_gap_extent = np.array(
        (
            (0, y_size - y_inflate),  # Bottom-left corner.
            (x_span, y_stride + y_inflate),  # Top-right corner.
        )
    )

    # Get the positions of the bottom-left corner of each gap region.
    x_gap_pos = x_stride * np.column_stack((np.arange(x_layout), np.zeros(x_layout)))
    y_gap_pos = y_stride * np.column_stack((np.zeros(y_layout), np.arange(y_layout)))

    # Compute the positions and extents of all the gaps in x & y.
    x_gap_masks = np.add(x_gap_extent, x_gap_pos[:, np.newaxis])
    y_gap_masks = np.add(y_gap_extent, y_gap_pos[:, np.newaxis])

    return np.concatenate((x_gap_masks, y_gap_masks))


def latrd_masks(module_layout: Vector, inflate: Vector = 0):
    """
    Find the regions to be masked on a Large Area Time Resolved Detector.

    A LATRD consists of a rectangular array of n_x × n_y LATRD modules, each in turn
    consisting of an array of 8 × 2 Timepix3 chips.  The gap regions between adjacent
    chips and between adjacent modules are calculated so that they can be masked.
    The masks for the gaps between adjacent chips can optionally be inflated by a
    specified number of pixels.

    There will be N = 8 × n_x × n_y + n_x + n_y - 2 rectangular regions to be masked.

    See Also:
        dxtbx.format.FormatNexusTimepix.gap_masks

    Args:
        module_layout:  The number of repeats (number x, number y) of LATRD modules
            in the detector.  If there are the same number of components in each
            dimension, that number can be passed as a scalar.
        inflate:  Width in pixels (number x, number y) to be added to each side of
            the masked gap region between Timepix3 chips within each LATRD module,
            in order to mask the edges of the chips.  If the gaps in both the x and y
            directions are to be inflated by the same number of pixels, that number
            can be passed as a scalar.

    Returns:
        The parameters of the regions to be masked, a numpy.ndarray of shape (N, 2, 2),
        where N is the number of gaps between components.  Each of the N elements in
        the first dimension is a pair of two-vectors.  Each pair represents the
        bottom-left and top-right pixel, respectively, of a gap region to be masked.
        Each pixel is represented by a (x, y) vector of its position.
    """
    # Generate the mask coordinates for the gaps between chips within the first
    # module.  Inflate the chip gap masks by a width of one pixel either side.
    chip_gap_mask_first_module = gap_masks(chip_size, chip_layout, chip_stride, inflate)

    # Find the positions of each of the modules.
    module_steps = np.concatenate(np.indices(module_layout).T)
    module_displacements = module_stride * module_steps

    # Place a copy of the chip gap masks over each module.
    chip_gap_mask_first_module = chip_gap_mask_first_module[np.newaxis, ...]
    module_displacements = module_displacements[:, np.newaxis, np.newaxis, :]
    chip_gap_masks = np.add(chip_gap_mask_first_module, module_displacements)
    chip_gap_masks = np.concatenate(chip_gap_masks)

    # Generate the mask coordinates for the gaps between modules.
    module_gap_mask = gap_masks(module_size, module_layout, module_stride)

    return np.concatenate((chip_gap_masks, module_gap_mask))


class FormatNexusTimepix(FormatNexus):
    """Format for binned images from a Timepix3 LATRD 'Tristan' detector."""

    @staticmethod
    def understand(image_file):
        """
        Check that the image format is recognised as binned images from a LATRD.

        The image format is understood if:
            • the NeXus tree contains an item
                "/entry/instrument/detector/module/data_size", and
            • the size recorded there corresponds to a known rectangular arrangement of
                LATRD modules.
        """
        known_module_layouts = (
            (1, 1),  # 1M detector
            (1, 2),  # 2M detector
            (2, 5),  # 10M detector
        )
        with h5py.File(image_file, "r") as handle:
            if "/entry/instrument/detector/module/data_size" in handle:
                size = handle["/entry/instrument/detector/module/data_size"]
                module_layout = object_layout(size[()], module_size, module_stride)
                if tuple(module_layout) in known_module_layouts:
                    return True

        return False

    def __init__(self, image_file, **kwargs):
        super().__init__(image_file, **kwargs)

        with h5py.File(image_file, "r") as handle:
            image_size = handle["/entry/instrument/detector/module/data_size"]
            self.module_layout = object_layout(image_size, module_size, module_stride)

    def _detector(self) -> Detector:
        """Ensure that the detector has appropriate masks for the non-sensitive area."""

        detector = super()._detector()

        masks = latrd_masks(self.module_layout, inflate=1)
        # Reshape each mask from a pair of vectors to the four-number list required
        # by dxtbx.model.Panel.add_mask()
        masks = masks.reshape(masks.shape[0], 4)
        # Add all the masks to the detector panel object.
        np.apply_along_axis(lambda mask: detector[0].add_mask(*mask), 1, masks)

        return detector
