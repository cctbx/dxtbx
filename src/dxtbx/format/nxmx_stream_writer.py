from __future__ import annotations

import logging
import os

import h5py
import numpy as np

from cctbx import factor_ev_angstrom

from dxtbx.format.nxmx_writer import NXmxWriter, get_compression

logger = logging.getLogger(__name__)


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
        # Optional SPARSE second image stream: the raw 2D spectrometer images, stored
        # under /entry/spectrum alongside a per-entry image_id index. Only frames that
        # carry a real reading append, so this stream is shorter than the detector
        # stream and is correlated back by image_id, not by position.
        self.spectrum_group = None
        self.spectrum_dset = None
        self.spectrum_index_dset = None
        self.spectrum_shape = None
        self.spectrum_dtype = None
        # 1D incident-beam spectrum on NXbeam: a shared (n_channels,) wavelength axis
        # plus resizable per-frame weights + present mask, grown by append_frame_metadata
        # in lockstep with the scalars so the data file is self-consistent even if the run
        # never finalizes. Distinct from the 2D /entry/spectrum image stream above; the
        # master VDSes these per-frame weights/present from the data files.
        self.incident_spectrum_axis_dset = None
        self.incident_spectrum_weights_dset = None
        self.incident_spectrum_present_dset = None
        self.incident_spectrum_channels = None
        # Per-source candidate wavelengths on NXbeam (incident_wavelength_<source>):
        # the raw per-sensor readings the resolver chose between, grown per frame
        # alongside the scalar incident_wavelength. One resizable dataset per source,
        # created lazily the first time that source appears; absent readings are NaN.
        self.candidate_wl_dsets: dict = {}

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
        spectrum_has_weights=None,
        spectrum_n_channels=None,
        wavelength_candidate_names=None,
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
            spectrum_has_weights: list (aligned with data_file_names) of bool -- did
                that data file archive any 1D incident spectrum (so it carries the
                per-frame weights/present datasets the master VDSes from).
            spectrum_n_channels: shared spectrometer channel count, or None.
            wavelength_candidate_names: list (aligned with data_file_names) of the
                candidate source names each data file archived (incident_wavelength_<name>
                datasets), used to VDS the per-source candidate wavelengths from the data
                files in the same sorted frame order. None when no API reported them.
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
        # (file_path, local_index, n_images) per virtual frame, in master order -- used
        # to build the NXbeam 1D-spectrum weights/present VDS over the data files.
        ordered_frames = []
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
                ordered_frames.extend(
                    (data_file_name, local_index, n_images)
                    for local_index in range(n_images)
                )
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
                ordered_frames.append((file_path, local_index, n_images))

        # Create the virtual dataset
        self.data_group = self.handle["entry"].create_group("data")
        self.data_group.attrs["NX_class"] = "NXdata"
        # Name the signal so the reader associates the group-level data_scale_factor
        # with the image VDS (not the "first member" fallback).
        self.data_group.attrs["signal"] = "data_000001"
        self.handle.create_virtual_dataset("/entry/data/data_000001", layout)

        self.write_per_frame_metadata(ordered_scale_factors, ordered_wavelengths)
        # The 1D incident spectrum is VDSed from the data files' per-frame weights/present
        # datasets (not copied), in the same sorted order, after the scalar
        # incident_wavelength exists (the variant attr points back to it).
        if spectrum_has_weights is not None:
            has_weights = {
                data_file_name: bool(flag)
                for (data_file_name, _), flag in zip(
                    data_file_names, spectrum_has_weights
                )
            }
            self._write_incident_spectrum_vds(
                ordered_frames, source_name, has_weights, spectrum_n_channels
            )
        if wavelength_candidate_names is not None:
            candidate_names_by_file = {
                data_file_name: set(names)
                for (data_file_name, _), names in zip(
                    data_file_names, wavelength_candidate_names
                )
            }
            self._write_candidate_wavelengths_vds(
                ordered_frames, source_name, candidate_names_by_file
            )
        self.handle.flush()

    def write_spectrum_master(
        self, spectrum_file_names, spectrum_image_ids, spectrum_shape, spectrum_dtype
    ):
        """Build the master's sparse /entry/spectrum VDS over the per-archiver files.

        spectrum_file_names: list of (data_file_path, n_spectrum_images) per archiver.
        spectrum_image_ids: list-of-lists of source image_ids, aligned with
            spectrum_file_names and concatenated in the same order as the VDS so the
            image_id index correlates each virtual spectrum frame to its detector frame.

        Sparse and unsorted by design: nothing reads this stream positionally against
        the detector frames, so it is concatenated in file order and correlated by
        image_id. Skipped entirely when no spectrometer images were archived.
        """
        total = sum(n for _, n in spectrum_file_names)
        if total == 0:
            return
        spectrum_shape = tuple(int(n) for n in spectrum_shape)
        dtype = np.dtype(spectrum_dtype)
        master_dir = os.path.dirname(self.params.output_file)
        layout = h5py.VirtualLayout(shape=(total, *spectrum_shape), dtype=dtype)
        start = 0
        for data_file_name, n_images in spectrum_file_names:
            if n_images == 0:
                continue
            # Basename relative to the master dir, same rule as the primary VDS, so the
            # master + data files relocate together.
            layout[start : start + n_images] = h5py.VirtualSource(
                os.path.relpath(data_file_name, master_dir),
                "/entry/spectrum/data_000001",
                shape=(n_images, *spectrum_shape),
                dtype=dtype,
            )
            start += n_images
        group = self.handle["entry"].create_group("spectrum")
        group.attrs["NX_class"] = "NXdata"
        group.attrs["signal"] = "data_000001"
        self.handle.create_virtual_dataset("/entry/spectrum/data_000001", layout)
        all_ids = [image_id for ids in spectrum_image_ids for image_id in ids]
        group.create_dataset("image_id", data=np.asarray(all_ids, dtype="i8"))
        self.handle.flush()

    def write_per_frame_metadata(
        self,
        scale_factors: list[float],
        wavelengths: list[float],
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

    def _append_incident_spectrum(self, frame_index, spectrum) -> None:
        """Append one frame's 1D incident spectrum to this data file's NXbeam, in
        lockstep with the per-frame scalars.

        On the first reading, create the shared 1-D wavelength axis
        (incident_wavelength_1Dspectrum, stored once because a fixed eV calibration
        gives every frame the same axis), the resizable per-frame weights
        (incident_wavelength_1Dspectrum_weights) and the present mask
        (incident_wavelength_1Dspectrum_present), back-filling flat placeholder rows for
        any frames already written before the first reading arrived. A missing or
        degenerate reading gets a flat placeholder row flagged present=0; the per-shot
        beam wavelength still comes from the scalar incident_wavelength, so the
        placeholder never perturbs it. Written per frame (not at finalize) so a reader
        that opens the data file directly sees the spectrum even if the run is aborted.
        """
        beam = self.handle.get("entry/instrument/beam")
        if beam is None:
            return
        if self.incident_spectrum_weights_dset is None:
            if spectrum is None:
                return  # no spectrum stream yet; create it on the first real reading
            n_channels = len(spectrum["weights"])
            self.incident_spectrum_channels = n_channels
            # The reader stores wavelengths and converts back to energies, so write the
            # axis as wavelengths (Angstrom), as nxmx_writer.add_beams does.
            axis = factor_ev_angstrom / np.asarray(spectrum["energies_eV"], dtype="f8")
            self.incident_spectrum_axis_dset = beam.create_dataset(
                "incident_wavelength_1Dspectrum", data=axis, dtype="f8"
            )
            self.incident_spectrum_axis_dset.attrs["units"] = "angstrom"
            if "incident_wavelength" in beam:
                beam["incident_wavelength"].attrs["variant"] = (
                    "incident_wavelength_1Dspectrum"
                )
            self.incident_spectrum_weights_dset = beam.create_dataset(
                "incident_wavelength_1Dspectrum_weights",
                (frame_index, n_channels),
                maxshape=(None, n_channels),
                chunks=(1, n_channels),
                dtype="f8",
            )
            self.incident_spectrum_present_dset = beam.create_dataset(
                "incident_wavelength_1Dspectrum_present",
                (frame_index,),
                maxshape=(None,),
                dtype="u1",
            )
            if frame_index > 0:
                self.incident_spectrum_weights_dset[...] = 1.0  # flat placeholder rows
                self.incident_spectrum_present_dset[...] = 0
        i = self.incident_spectrum_weights_dset.shape[0]
        self.incident_spectrum_weights_dset.resize(i + 1, axis=0)
        self.incident_spectrum_present_dset.resize(i + 1, axis=0)
        if (
            spectrum is not None
            and len(spectrum["weights"]) == self.incident_spectrum_channels
        ):
            self.incident_spectrum_weights_dset[i] = np.asarray(
                spectrum["weights"], dtype="f8"
            )
            self.incident_spectrum_present_dset[i] = 1
        else:
            self.incident_spectrum_weights_dset[i] = 1.0  # flat placeholder
            self.incident_spectrum_present_dset[i] = 0

    def _write_incident_spectrum_vds(
        self, ordered_frames, source_name, has_weights, n_channels
    ) -> None:
        """Build the master's NXbeam 1D-spectrum VDS over the data files' per-frame
        weights/present datasets, in the same sorted order as the image VDS.

        ordered_frames: list of (data_file_path, local_index, n_images) in master order.
        has_weights: dict data_file_path -> bool (did that file archive any spectrum, so
            it carries the per-frame weights/present datasets to source from).
        n_channels: shared channel count.

        The shared 1-D wavelength axis is copied once from the first contributing data
        file. Frames from a file that archived no spectrum at all (no weights dataset)
        fall to the layout fill value (flat weights, present=0) -- never expected at a
        non-trivial FEE rate, but kept robust.
        """
        if not n_channels or not any(has_weights.values()):
            return
        weights_path = "/entry/instrument/beam/incident_wavelength_1Dspectrum_weights"
        present_path = "/entry/instrument/beam/incident_wavelength_1Dspectrum_present"
        total = len(ordered_frames)
        weights_layout = h5py.VirtualLayout(shape=(total, n_channels), dtype="f8")
        present_layout = h5py.VirtualLayout(shape=(total,), dtype="u1")
        for virtual_index, (file_path, local_index, n_images) in enumerate(
            ordered_frames
        ):
            if not has_weights.get(file_path):
                continue  # no weights dataset in this file -> layout fill value
            rel = source_name[file_path]
            weights_layout[virtual_index] = h5py.VirtualSource(
                rel, weights_path, shape=(n_images, n_channels), dtype="f8"
            )[local_index]
            present_layout[virtual_index] = h5py.VirtualSource(
                rel, present_path, shape=(n_images,), dtype="u1"
            )[local_index]
        beam = self.handle["entry/instrument/beam"]
        # Copy the shared 1-D wavelength axis once from a contributing data file.
        axis = None
        for file_path, _, _ in ordered_frames:
            if has_weights.get(file_path):
                with h5py.File(file_path, "r") as data_file:
                    axis = data_file[
                        "entry/instrument/beam/incident_wavelength_1Dspectrum"
                    ][()]
                break
        if axis is None:
            return
        beam.create_dataset("incident_wavelength_1Dspectrum", data=axis, dtype="f8")
        beam["incident_wavelength_1Dspectrum"].attrs["units"] = "angstrom"
        # fillvalue covers frames from a data file that archived no spectrum (no weights
        # dataset to source): flat weights (1.0) marked absent (present=0).
        self.handle.create_virtual_dataset(weights_path, weights_layout, fillvalue=1.0)
        self.handle.create_virtual_dataset(present_path, present_layout, fillvalue=0)
        if "incident_wavelength" in beam:
            beam["incident_wavelength"].attrs["variant"] = (
                "incident_wavelength_1Dspectrum"
            )

    def _write_candidate_wavelengths_vds(
        self, ordered_frames, source_name, candidate_names_by_file
    ) -> None:
        """Build the master's NXbeam per-source candidate-wavelength VDS over the data
        files' per-frame incident_wavelength_<source> datasets, in the same sorted order
        as the image VDS.

        ordered_frames: list of (data_file_path, local_index, n_images) in master order.
        candidate_names_by_file: dict data_file_path -> set of candidate source names
            that file archived. Frames from a file that did not archive a given source
            fall to the layout NaN fill value.
        """
        all_names = sorted(
            {name for names in candidate_names_by_file.values() for name in names}
        )
        if not all_names:
            return
        beam = self.handle["entry/instrument/beam"]
        total = len(ordered_frames)
        for name in all_names:
            path = f"/entry/instrument/beam/incident_wavelength_{name}"
            layout = h5py.VirtualLayout(shape=(total,), dtype="f8")
            for virtual_index, (file_path, local_index, n_images) in enumerate(
                ordered_frames
            ):
                if name not in candidate_names_by_file.get(file_path, ()):
                    continue  # source absent in this file -> layout fill value (NaN)
                layout[virtual_index] = h5py.VirtualSource(
                    source_name[file_path], path, shape=(n_images,), dtype="f8"
                )[local_index]
            self.handle.create_virtual_dataset(path, layout, fillvalue=np.nan)
            beam[f"incident_wavelength_{name}"].attrs["units"] = "angstrom"

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
        self,
        data_scale_factor: float,
        incident_wavelength: float,
        incident_spectrum=None,
        wavelength_candidates=None,
    ) -> None:
        """Append one frame's data_scale_factor, incident_wavelength, the (optional) 1D
        incident spectrum, and the (optional) per-source candidate wavelengths.

        Called once per image written to this data file, in lockstep with
        append_image, so element i of each 1-D dataset describes frame i. Writing it
        per frame (rather than batched at finalize) means a reader that opens the data
        file directly recovers photons per frame even if the run never finalized.
        incident_spectrum is a dict with energies_eV + weights, or None when this frame
        carried no calibrated spectrum. wavelength_candidates is a dict source-name ->
        wavelength (Angstrom) of the raw per-sensor readings the resolver chose between,
        or None/empty when the API carries a single wavelength source.
        """
        if self.scale_factor_dset is None:
            self.initialize_dataset()
        n = self.scale_factor_dset.shape[0]
        self.scale_factor_dset.resize(n + 1, axis=0)
        self.scale_factor_dset[n] = data_scale_factor
        if self.wavelength_dset is not None:
            self.wavelength_dset.resize(n + 1, axis=0)
            self.wavelength_dset[n] = incident_wavelength
        self._append_incident_spectrum(n, incident_spectrum)
        self._append_wavelength_candidates(n, wavelength_candidates)

    def _append_wavelength_candidates(self, frame_index, candidates) -> None:
        """Append one frame's per-source candidate wavelengths to NXbeam, in lockstep
        with the scalar incident_wavelength.

        Each source gets its own resizable dataset incident_wavelength_<source>,
        created the first time that source appears and back-filled with NaN for the
        frames already written (a fixed fill value, since a missing reading carries no
        wavelength). Every known dataset is advanced one row per frame -- with NaN when
        that source had no reading this frame -- so all stay aligned with
        incident_wavelength at index frame_index. A no-op for APIs that report a single
        wavelength source (candidates always empty).
        """
        beam = self.handle.get("entry/instrument/beam")
        if beam is None:
            return
        candidates = candidates or {}
        for name in set(self.candidate_wl_dsets) | set(candidates):
            dset = self.candidate_wl_dsets.get(name)
            if dset is None:
                dset = beam.create_dataset(
                    f"incident_wavelength_{name}",
                    (frame_index,),
                    maxshape=(None,),
                    dtype="f8",
                    fillvalue=np.nan,
                )
                dset.attrs["units"] = "angstrom"
                self.candidate_wl_dsets[name] = dset
            dset.resize(frame_index + 1, axis=0)
            value = candidates.get(name)
            dset[frame_index] = np.nan if value is None else value

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

    def initialize_spectrum_dataset(self, spectrum_shape, spectrum_dtype):
        """Lazily create the sparse /entry/spectrum NXdata group on first append.

        Created only when a run actually carries 2D spectrometer images, so a
        spectrum-free run (e.g. Dectris) never grows an /entry/spectrum group and its
        file stays byte-for-byte identical to the pre-feature output. Sibling image_id
        dataset records which detector frame each stored spectrum belongs to.
        """
        self.spectrum_shape = tuple(int(n) for n in spectrum_shape)
        self.spectrum_dtype = np.dtype(spectrum_dtype)
        group = self.handle["entry"].create_group("spectrum")
        group.attrs["NX_class"] = "NXdata"
        group.attrs["signal"] = "data_000001"
        self.spectrum_group = group
        self.spectrum_dset = group.create_dataset(
            "data_000001",
            (0, *self.spectrum_shape),
            maxshape=(None, *self.spectrum_shape),
            chunks=(1, *self.spectrum_shape),
            dtype=self.spectrum_dtype,
            compression=get_compression(self.params.compression),
        )
        self.spectrum_index_dset = group.create_dataset(
            "image_id", (0,), maxshape=(None,), dtype="i8"
        )

    def append_spectrum_image(
        self, image_data, image_id, spectrum_shape, spectrum_dtype, compressed=True
    ):
        """Append one sparse 2D spectrometer image (and its source image_id).

        Called only for frames that carry a real reading, so this stream is shorter
        than the detector stream; the image_id index lets a reader correlate each
        stored spectrum back to its detector frame without positional alignment.
        """
        if self.spectrum_dset is None:
            self.initialize_spectrum_dataset(spectrum_shape, spectrum_dtype)
        current_size = self.spectrum_dset.shape[0]
        self.spectrum_dset.resize(current_size + 1, axis=0)
        if compressed:
            chunk_index = (current_size,) + (0,) * len(self.spectrum_shape)
            self.spectrum_dset.id.write_direct_chunk(
                chunk_index, image_data, filter_mask=0
            )
        else:
            self.spectrum_dset[-1:] = image_data
        self.spectrum_index_dset.resize(current_size + 1, axis=0)
        self.spectrum_index_dset[current_size] = image_id
