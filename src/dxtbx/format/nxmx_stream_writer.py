from __future__ import annotations

import os

import h5py
import numpy as np

from dxtbx.format.nxmx_writer import NXmxWriter, get_compression


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
        # Per-frame scaling metadata datasets, grown one element per appended image
        # (see append_frame_metadata) so the data file is self-consistent even if the
        # run is aborted before finalize.
        self.scale_factor_dset = None
        self.wavelength_dset = None

    def __call__(self, experiments=None, imageset=None):
        """
        This method is overwritten to remove the self.append_all_frames() call.
        Archivers use this method to setup the file writer before they have received
        data to save. They will then build up the h5 file incrementally.
        Because Archivers can archive compressed data, they need to know the image shape.
        It gets the image shape from the reference experiments
        """
        if experiments or imageset:
            self.setup(experiments, imageset)
        self.validate()
        self.construct_detector()
        self.add_all_beams()
        self.add_scan_and_gonio()

        value = np.dtype(self.params.dtype).itemsize * 8
        self.handle["/entry/instrument/detector"].create_dataset(
            "bit_depth_readout", (), "i"
        )
        self.handle["/entry/instrument/detector/bit_depth_readout"][()] = value

        if experiments:
            n_panels = len(experiments[0].detector)
            panel_shape = experiments[0].detector[0].get_image_size()
            # NOT SURE IF THE ORDER OF THE PANEL SHAPE IS CORRECT
            if n_panels == 1:
                self.image_shape = [panel_shape[1], panel_shape[0]]
            else:
                self.image_shape = (n_panels, panel_shape[1], panel_shape[0])

    def write_master(
        self,
        data_file_names,
        sort_values=None,
        data_scale_factors=None,
        incident_wavelengths=None,
    ):
        """
        This method is used by the control to link together the data files written
        by the archivers.
        Args:
            data_file_names: list of tuples. each individual tuple has the first
                element as the paths to the archived h5 files. The second element
                is the total number of images archived to that file.
            sort_values: list of lists. First level corresponds to data files,
                second level contains sort values for each image in that file.
            data_scale_factors: list of lists (same nesting as sort_values) of the
                per-frame NeXus data_scale_factor; written as a 1-D dataset under
                NXdata in the same order as the sorted VDS frames so the reader can
                recover photons from the archived (native-unit) pixels per frame.
            incident_wavelengths: list of lists (same nesting) of the per-frame
                incident wavelength; replaces the scalar incident_wavelength that
                add_all_beams wrote, so the re-read beam is per-frame correct.
        """
        total_images = sum(n_images for _, n_images in data_file_names)
        total_shape = (total_images, *self.image_shape)

        # Store each VDS source path relative to the master file's own directory (a
        # basename, since the master and its data files are co-located). HDF5 resolves a
        # relative VDS source against the master file's directory, so the master + data
        # files can be moved together to any directory (or machine) and still resolve,
        # regardless of the reader's working directory. An absolute path would break that
        # move, and an output-dir-relative path (e.g. rNNNN/rNNNN.h5) silently yields
        # fill-value (blank) frames when read from a different directory.
        master_dir = os.path.dirname(self.params.output_file)
        source_name = {
            data_file_name: os.path.relpath(data_file_name, master_dir)
            for data_file_name, _ in data_file_names
        }

        # Create a virtual layout for the combined dataset. The per-frame metadata
        # arrays are assembled in the SAME frame order as the VDS so element i of
        # each 1-D array describes virtual frame i.
        layout = h5py.VirtualLayout(shape=total_shape, dtype=self.params.dtype)
        ordered_scale_factors = []
        ordered_wavelengths = []
        if sort_values is None:
            start = 0
            for file_index, (data_file_name, n_images) in enumerate(data_file_names):
                layout[start : start + n_images] = h5py.VirtualSource(
                    source_name[data_file_name],
                    "/entry/data/data_000001",
                    shape=(n_images, *self.image_shape),
                    dtype=self.params.dtype,
                )
                start += n_images
                if data_scale_factors is not None:
                    ordered_scale_factors.extend(data_scale_factors[file_index])
                if incident_wavelengths is not None:
                    ordered_wavelengths.extend(incident_wavelengths[file_index])
        else:
            # Create list of (file_path, local_index, sort_value, n_images, scale,
            # wavelength) for all images so the per-frame metadata is carried through
            # the same sort that orders the pixel frames.
            all_images = []
            for file_index, ((data_file_name, n_images), file_sort_values) in enumerate(
                zip(data_file_names, sort_values)
            ):
                for local_index, sort_value in enumerate(file_sort_values):
                    scale = (
                        data_scale_factors[file_index][local_index]
                        if data_scale_factors is not None
                        else None
                    )
                    wavelength = (
                        incident_wavelengths[file_index][local_index]
                        if incident_wavelengths is not None
                        else None
                    )
                    all_images.append(
                        (
                            data_file_name,
                            local_index,
                            sort_value,
                            n_images,
                            scale,
                            wavelength,
                        )
                    )

            # Sort by the sort values
            all_images.sort(key=lambda x: x[2])

            # Map each image to its sorted position
            for virtual_index, (
                file_path,
                local_index,
                sort_value,
                n_images,
                scale,
                wavelength,
            ) in enumerate(all_images):
                layout[virtual_index] = h5py.VirtualSource(
                    source_name[file_path],
                    "/entry/data/data_000001",
                    shape=(n_images, *self.image_shape),
                    dtype=self.params.dtype,
                )[local_index]
                ordered_scale_factors.append(scale)
                ordered_wavelengths.append(wavelength)

        # Create the virtual dataset
        self.data_group = self.handle["entry"].create_group("data")
        self.data_group.attrs["NX_class"] = "NXdata"
        # Name the signal so the reader associates the group-level data_scale_factor
        # with the image VDS (not the "first member" fallback).
        self.data_group.attrs["signal"] = "data_000001"
        self.handle.create_virtual_dataset("/entry/data/data_000001", layout)

        self.write_per_frame_metadata(ordered_scale_factors, ordered_wavelengths)
        self.handle.flush()

    def write_per_frame_metadata(
        self, scale_factors: list[float], wavelengths: list[float]
    ) -> None:
        """Write the per-frame data_scale_factor and incident_wavelength datasets.

        The lists must be in the same frame order as this file's frames -- sorted
        VDS order for a master, arrival order for a per-archiver data file (the order
        append_image was called). Both the data files and the master carry these so a
        reader that opens either directly recovers photons per frame (the reader
        applies data_scale_factor in get_raw_data; the panel gain stays 1.0). The
        scale factor is omitted when no frame needs scaling (every factor == 1.0), so
        paths whose data is already in photons (e.g. Dectris) are read back
        byte-for-byte unchanged. The incident_wavelength scalar written by
        add_all_beams is replaced by the per-frame array.
        """
        if scale_factors and any(factor != 1.0 for factor in scale_factors):
            # NeXus: corrected = (data + offset) * data_scale_factor. Applied
            # per-frame by the reader to recover photons from native-unit pixels.
            self.data_group.create_dataset(
                "data_scale_factor", data=np.asarray(scale_factors, dtype="f8")
            )

        if wavelengths and "entry/instrument/beam" in self.handle:
            beam = self.handle["entry/instrument/beam"]
            # h5py will not overwrite a dataset in place, so drop the scalar that
            # add_all_beams wrote before recreating it as the per-frame array.
            if "incident_wavelength" in beam:
                del beam["incident_wavelength"]
            beam.create_dataset(
                "incident_wavelength", data=np.asarray(wavelengths, dtype="f8")
            )
            beam["incident_wavelength"].attrs["units"] = "angstrom"

    def initialize_dataset(self):
        if self.data_group is None:
            self.data_group = self.handle["entry"].create_group("data")
            self.data_group.attrs["NX_class"] = "NXdata"
            # Name the signal dataset explicitly so the reader associates the
            # group-level data_scale_factor / data_offset with the image data, rather
            # than relying on the "first member" fallback (which would otherwise pick
            # whichever NXdata child sorts first).
            self.data_group.attrs["signal"] = "data_000001"
            self.dset = self.data_group.create_dataset(
                "data_000001",
                (0, *self.image_shape),
                maxshape=(None, *self.image_shape),
                chunks=(1, *self.image_shape),  # Each image is one chunk
                dtype=self.params.dtype,
                compression=get_compression(self.params.compression),
            )
            # Resizable per-frame scaling metadata, grown one element per image by
            # append_frame_metadata. Written incrementally (not at finalize) so the
            # data file stays self-consistent if the run is aborted mid-stream.
            self.scale_factor_dset = self.data_group.create_dataset(
                "data_scale_factor", (0,), maxshape=(None,), dtype="f8"
            )
            beam = self.handle.get("entry/instrument/beam")
            if beam is not None:
                # Replace the scalar incident_wavelength add_all_beams wrote with a
                # resizable per-frame array (h5py cannot overwrite a dataset in place).
                if "incident_wavelength" in beam:
                    del beam["incident_wavelength"]
                self.wavelength_dset = beam.create_dataset(
                    "incident_wavelength", (0,), maxshape=(None,), dtype="f8"
                )
                self.wavelength_dset.attrs["units"] = "angstrom"

    def append_frame_metadata(
        self, data_scale_factor: float, incident_wavelength: float
    ) -> None:
        """Append one frame's data_scale_factor and incident_wavelength.

        Called once per image written to this data file, in lockstep with
        append_image, so element i of each 1-D dataset describes frame i. Writing it
        per frame (rather than batched at finalize) means a reader that opens the data
        file directly recovers photons per frame even if the run never finalized.
        """
        if self.scale_factor_dset is None:
            self.initialize_dataset()
        n = self.scale_factor_dset.shape[0]
        self.scale_factor_dset.resize(n + 1, axis=0)
        self.scale_factor_dset[n] = data_scale_factor
        if self.wavelength_dset is not None:
            m = self.wavelength_dset.shape[0]
            self.wavelength_dset.resize(m + 1, axis=0)
            self.wavelength_dset[m] = incident_wavelength

    def append_image(self, image_data, compressed=False):
        """Method for uncompressed data"""
        # Make sure dataset is initialized
        if self.dset is None:
            self.initialize_dataset()

        current_size = self.dset.shape[0]
        self.dset.resize(current_size + 1, axis=0)

        if compressed:
            # Calculate chunk index for this image
            # For a dataset with shape (N, height, width) chunked as (1, height, width)
            if len(self.image_shape) == 2:
                chunk_index = (current_size, 0, 0)
            # For a dataset with shape (N, n_panels, height, width) chunked as (1, n_panels, height, width)
            elif len(self.image_shape) == 3:
                chunk_index = (current_size, 0, 0, 0)
            # Write compressed data directly to chunk
            self.dset.id.write_direct_chunk(chunk_index, image_data, filter_mask=0)
        else:
            self.dset[-1:] = image_data

        self.image_count += 1
