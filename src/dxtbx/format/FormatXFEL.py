"""Mixin for XFEL format classes that produce a per-wavelength ImageSequence."""

from __future__ import annotations

import os


class FormatXFEL:
    """Mixin for XFEL format classes.

    Subclasses must implement get_wavelengths() returning a Python list of
    per-frame wavelengths in Angstrom.

    get_imageset() builds a regular ImageSequence with:
      - XFELBeam (direction + divergence, no wavelength or wavelength_range)
      - zero-oscillation Scan with "wavelength" property (one entry per frame)
    """

    def get_wavelengths(self):
        """Return list of per-frame wavelengths (Angstrom). Override in subclass."""
        raise NotImplementedError

    @classmethod
    def get_imageset(
        cls,
        filenames,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
        format_kwargs=None,
        single_file_indices=None,
        **kwargs,
    ):
        from scitbx.array_family import flex

        from dxtbx.format.FormatMultiImage import _add_static_mask_to_iset
        from dxtbx.imageset import ImageSequence, ImageSetData
        from dxtbx.model import Scan
        from dxtbx.model.beam import BeamFactory

        # The parent's get_imageset switches on these flags; for XFEL we always
        # build an ImageSequence directly so they're irrelevant here.
        kwargs.pop("as_imageset", None)
        kwargs.pop("as_sequence", None)
        kwargs.pop("lazy", None)
        kwargs.pop("check_format", None)
        kwargs.pop("template", None)

        if isinstance(filenames, str):
            filenames = [filenames]
        elif len(filenames) > 1:
            assert len(set(filenames)) == 1
            filenames = filenames[0:1]
        filenames = [os.path.abspath(x) for x in filenames]

        if format_kwargs is None:
            format_kwargs = {}

        # Single shared format instance — its _start cache holds the detector
        # and beam, so subsequent get_beam(0)/get_detector(0) calls are free.
        # This replaces FormatMultiImage.get_imageset's per-frame loop, which
        # called format_instance.get_beam(i) for every frame in
        # single_file_indices only to discard the results once the XFEL mixin
        # overwrote them with a single shared XFELBeam.
        fmt = cls.get_instance(filenames[0], **format_kwargs)
        num_images = fmt.get_num_images()

        if single_file_indices is not None:
            assert len(single_file_indices)
            single_file_indices = flex.size_t(single_file_indices)
            indices = list(single_file_indices)
        else:
            indices = list(range(num_images))
            single_file_indices = flex.size_t(indices)

        n = len(indices)

        reader = cls.get_reader()(filenames, num_images=num_images, **format_kwargs)
        # XFEL stills have goniometer=None, so the multi-axis shadow path
        # in Format.get_masker() is never taken; pass None through.
        masker = fmt.get_masker(goniometer=None)
        isetdata = ImageSetData(
            reader=reader,
            masker=masker,
            vendor=fmt.get_vendortype(),
            params=format_kwargs,
            format=cls,
            template=filenames[0],
        )

        ref_beam = fmt.get_beam(0)
        if detector is None:
            detector = fmt.get_detector(0)

        wavelengths_all = fmt.get_wavelengths()  # list[float], Angstrom, full file
        if len(indices) != len(wavelengths_all):
            # Subset of the file (e.g. loading a composite stills output).
            # single_file_indices holds the absolute 0-based frame positions;
            # select the matching wavelengths so the property length equals n.
            wavelengths = [wavelengths_all[i] for i in indices]
        else:
            wavelengths = wavelengths_all

        # Carry the instrument beam's polarization (and flux/transmission/probe
        # /etc.) onto the XFELBeam so that monochromatic beams derived from it
        # via get_monochromatic_beam() retain the correct values for the
        # Lorentz-polarization correction downstream.
        xfel_beam = BeamFactory.make_xfel_beam(
            direction=ref_beam.get_sample_to_source_direction(),
            divergence=ref_beam.get_divergence(),  # degrees (Python default)
            sigma_divergence=ref_beam.get_sigma_divergence(),
            polarization_normal=ref_beam.get_polarization_normal(),
            polarization_fraction=ref_beam.get_polarization_fraction(),
            flux=ref_beam.get_flux(),
            transmission=ref_beam.get_transmission(),
            probe=ref_beam.get_probe(),
            sample_to_source_distance=ref_beam.get_sample_to_source_distance(),
        )
        seq_scan = Scan((1, n), (0.0, 0.0))
        seq_scan.set_property("wavelength", flex.double(wavelengths))

        iset = ImageSequence(
            isetdata,
            beam=xfel_beam,
            detector=detector,
            goniometer=None,  # no goniometer for XFEL stills
            scan=seq_scan,
            indices=single_file_indices,
        )

        _add_static_mask_to_iset(fmt, iset)

        return iset
