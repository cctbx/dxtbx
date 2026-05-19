"""Mixin for XFEL format classes that produce XFELImageSequence."""


class FormatXFEL:
    """Mixin for XFEL format classes that produce XFELImageSequence.

    The parent format class (e.g. FormatNXmx) handles detector/beam/geometry.
    This mixin wraps that imageset in XFELImageSequence with per-frame wavelengths.
    Subclasses must implement get_wavelengths().
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
        from dxtbx.imageset import XFELImageSequence
        from dxtbx.model import Scan

        # Let the parent class build the raw ImageSet (populates beam/detector models)
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

        # Read per-frame wavelengths via the subclass implementation
        fmt = cls.get_instance(filenames[0], **(format_kwargs or {}))
        wavelengths = fmt.get_wavelengths()

        n = len(raw_iset)
        return XFELImageSequence(
            raw_iset.data(),
            raw_iset.indices(),
            beam=raw_iset.get_beam(0),
            detector=raw_iset.get_detector(0),
            goniometer=None,  # XFEL stills have no goniometer
            scan=Scan((1, n), (0.0, 0.0)),  # zero-osc: is_still()=True
            wavelengths=wavelengths,
        )
