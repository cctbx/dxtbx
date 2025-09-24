from dxtbx.format.nxmx_writer import NXmxWriter
from dxtbx.format.nxmx_writer import phil_scope as nxmx_writer_phil_scope
from dxtbx.format.nxmx_writer import get_compression
import h5py
import hdf5plugin
from libtbx.phil import parse


def compress_with_hdf5_filters(data, params):
    """Compress data using HDF5 filters without writing to disk"""
    # Create an in-memory HDF5 file using the core driver

    with h5py.File(name="in-memory", driver="core", backing_store=False, mode="w") as f:
        # Create dataset with the desired compression settings
        dset = f.create_dataset(
            "data",
            data=data,
            chunks=data.shape,
            compression=get_compression(params),
        )

        # Get the compressed chunk directly
        chunk_info = dset.id.get_chunk_info(0)
        address = chunk_info.byte_offset
        size = chunk_info.size

        # Read the raw compressed bytes from the in-memory file
        data_compressed_encoded = f.id.get_file_image()[address : address + size]
    return data_compressed_encoded


class NXmxStreamWriter(NXmxWriter):
    """
    https://www.dectris.com/en/support/downloads/header-docs/nexus/
    https://manual.nexusformat.org/index.html

    This object adds to NXmxWriter to be compatible with streaming.
    """

    def __init__(self, params, experiments=None, imageset=None):
        self.params = params
        self.detector = None
        self.handle = None

        self.image_count = 0
        self.data_group = None
        self.dset = None

    def __call__(self, experiments=None, imageset=None, in_memory=True):
        """
        This method is overwritten to remove the self.append_all_frames() call.
        Archivers use this method to setup the file writer before they have received
        data to save. They will then build up the h5 file incrementally.

        Because Archivers can archive compressed data, they need to know the image shape.
        It gets the image shape from the reference experiments
        """
        if experiments or imageset:
            self.setup(experiments, imageset, in_memory)
        self.validate()
        self.construct_detector()
        self.add_all_beams()
        self.add_scan_and_gonio()

        if experiments:
            n_panels = len(experiments[0].detector)
            panel_shape = experiments[0].detector[0].get_image_size()
            if n_panels == 1:
                self.image_shape = panel_shape
            else:
                self.image_shape = (n_panels, *panel_shape)

        if in_memory:
            # This ensures that the h5 file is written to the buffer.
            self.handle.flush()

    def write_master(self, data_file_names, sort_values=None):
        """
        This method is used by the control to link together the data files written
        by the archivers.
        Args:
            data_file_names: list of tuples. each individual tuple has the first
                element as the paths to the archived h5 files. The second element
                is the total number of images archived to that file.
            sort_values: list of lists. First level corresponds to data files,
                second level contains sort values for each image in that file.
        """
        total_images = sum(n_images for _, n_images in data_file_names)
        total_shape = (total_images, *self.image_shape)

        # Create a virtual layout for the combined dataset
        layout = h5py.VirtualLayout(shape=total_shape, dtype=self.params.dtype)
        if sort_values is None:
            start = 0
            for data_file_name, n_images in data_file_names:
                layout[start : start + n_images] = h5py.VirtualSource(
                    data_file_name,
                    "/entry/data/data_000001",
                    shape=(n_images, *self.image_shape),
                    dtype=self.params.dtype,
                )
                start += n_images
        else:
            # Create list of (file_path, local_index, sort_value) for all images
            all_images = []
            for file_index, ((data_file_name, n_images), file_sort_values) in enumerate(
                zip(data_file_names, sort_values)
            ):
                for local_index, sort_value in enumerate(file_sort_values):
                    all_images.append(
                        (data_file_name, local_index, sort_value, n_images)
                    )

            # Sort by the sort values
            all_images.sort(key=lambda x: x[2])

            # Map each image to its sorted position
            for virtual_index, (
                file_path,
                local_index,
                sort_value,
                n_images,
            ) in enumerate(all_images):
                layout[virtual_index] = h5py.VirtualSource(
                    file_path,
                    "/entry/data/data_000001",
                    shape=(n_images, *self.image_shape),
                    dtype=self.params.dtype,
                )[local_index]

        # Create the virtual dataset
        self.data_group = self.handle["entry"].create_group("data")
        self.data_group.attrs["NX_class"] = "NXdata"
        self.handle.create_virtual_dataset("/entry/data/data_000001", layout)
        self.handle.flush()

    def initialize_dataset(self):
        if self.data_group is None:
            self.data_group = self.handle["entry"].create_group("data")
            self.data_group.attrs["NX_class"] = "NXdata"

            self.dset = self.data_group.create_dataset(
                "data_000001",
                (0, *self.image_shape),
                maxshape=(None, *self.image_shape),
                chunks=(1, *self.image_shape),  # Each image is one chunk
                dtype=self.params.dtype,
                compression=get_compression(self.params.compression),
            )

    def append_compressed_image(self, compressed_data):
        """Append already-compressed image data directly"""
        # Make sure dataset is initialized
        if self.dset is None:
            self.initialize_dataset()

        # Resize dataset to accommodate new image
        current_size = self.dset.shape[0]
        self.dset.resize(current_size + 1, axis=0)

        # Calculate chunk index for this image
        # For a dataset with shape (N, height, width) chunked as (1, height, width)
        chunk_index = (current_size, 0, 0)

        # Write compressed data directly to chunk
        self.dset.id.write_direct_chunk(chunk_index, compressed_data, filter_mask=0)
        self.image_count += 1

    def append_image(self, image_data):
        """Method for uncompressed data"""
        # Make sure dataset is initialized
        if self.dset is None:
            self.initialize_dataset()
        # Resize and add image the traditional way
        current_size = self.dset.shape[0]
        self.dset.resize(current_size + 1, axis=0)
        self.dset[-1:] = image_data
        self.image_count += 1
