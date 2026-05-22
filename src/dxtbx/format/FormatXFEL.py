"""Mixin for XFEL format classes that produce a per-wavelength ImageSequence."""


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
        **kwargs,
    ):
        from dxtbx.imageset import ImageSequence
        from dxtbx.model import Scan
        from dxtbx.model.beam import BeamFactory
        from scitbx.array_family import flex

        kwargs.pop("as_imageset", None)
        kwargs.pop("as_sequence", None)

        raw_iset = super().get_imageset(
            filenames,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            as_imageset=True,
            format_kwargs=format_kwargs,
            **kwargs,
        )

        fmt = cls.get_instance(filenames[0], **(format_kwargs or {}))
        wavelengths_all = fmt.get_wavelengths()  # list[float], Angstrom, full file

        ref_beam = raw_iset.get_beam(0)
        xfel_beam = BeamFactory.make_xfel_beam(
            direction=ref_beam.get_sample_to_source_direction(),
            divergence=ref_beam.get_divergence(),  # degrees (Python default)
            sigma_divergence=ref_beam.get_sigma_divergence(),
        )

        n = len(raw_iset)
        indices = list(raw_iset.indices())
        if len(indices) != len(wavelengths_all):
            # Subset of the file (e.g. loading a composite stills output).
            # raw_iset.indices() holds the absolute 0-based frame positions;
            # select the matching wavelengths so the property length equals n.
            wavelengths = [wavelengths_all[i] for i in indices]
        else:
            wavelengths = wavelengths_all
        seq_scan = Scan((1, n), (0.0, 0.0))
        seq_scan.set_property("wavelength", flex.double(wavelengths))

        return ImageSequence(
            raw_iset.data(),
            raw_iset.indices(),
            xfel_beam,
            raw_iset.get_detector(0),
            None,       # no goniometer for XFEL stills
            seq_scan,
        )
