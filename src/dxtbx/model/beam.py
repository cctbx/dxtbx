import math

import pycbf

import libtbx.phil

from dxtbx_model_ext import MonochromaticBeam, TOFBeam

beam_phil_scope = libtbx.phil.parse(
    """
  beam
    .expert_level = 1
    .short_caption = "Beam overrides"
  {

    wavelength = None
      .type = float
      .help = "Override the beam wavelength"

    direction = None
      .type = floats(size=3)
      .help = "Override the sample to source direction"
      .short_caption = "Sample to source direction"

    polarization_normal = None
      .type = floats(size=3)
      .help = "Override the polarization normal"
      .short_caption = "Polarization normal"

    polarization_fraction = None
      .type = float(value_min=0.0, value_max=1.0)
      .help = "Override the polarization fraction"
      .short_caption = "Polarization fraction"

    sample_to_moderator_distance = None
        .type = float
        .help = "Override sample to moderator distance"

  }
"""
)


class BeamFactory:
    """A factory class for beam objects, which encapsulate standard beam
    models. In cases where a full cbf description is available this
    will be used, otherwise simplified descriptions can be applied."""

    @staticmethod
    def from_phil(params, reference=None):
        """
        Convert the phil parameters into a beam model

        """

        def check_for_required_basic_params(params):
            if params.beam.direction is None:
                raise RuntimeError("Cannot create beam: direction not set")

        def check_for_required_tof_params(params):
            if params.beam.sample_to_moderator_distance is None:
                raise RuntimeError(
                    "Cannot create ToF beam: sample_to_moderator_distance not set"
                )

        def check_for_required_monochromatic_params(params):
            if params.beam.wavelength is None:
                raise RuntimeError(
                    "Cannot create monochromatic beam: wavelength not set"
                )

        def tof_from_phil(params, beam):

            check_for_required_basic_params(params)
            check_for_required_tof_params(params)

            if beam is None:
                beam = TOFBeam()

            beam.set_sample_to_source_direction(params.beam.direction)
            beam.set_sample_to_moderator_distance(params.sample_to_moderator_distance)

            return beam

        def monochromatic_from_phil(params, beam):

            check_for_required_basic_params(params)
            check_for_required_monochromatic_params(params)

            if beam is None:
                beam = MonochromaticBeam()

            beam.set_wavelength(params.beam.wavelength)
            if params.beam.polarization_normal is not None:
                beam.set_polarization_normal(params.beam.polarization_normal)
            if params.beam.polarization_fraction is not None:
                beam.set_polarization_fraction(params.beam.polarization_fraction)
            return beam

        try:
            return monochromatic_from_phil(reference, params)
        except RuntimeError:
            return tof_from_phil(reference, params)

    @staticmethod
    def from_dict(d, t=None):
        """Convert the dictionary to a beam model

        Params:
            d The dictionary of parameters
            t The template dictionary to use

        Returns:
            The beam model
        """

        def check_for_required_keys(d, required_keys, beam_type):
            for i in required_keys:
                if i not in d:
                    raise RuntimeError(
                        f"Cannot create {beam_type}: {i} not in dictionary"
                    )

        basic_required_keys = ["direction"]
        monochromatic_required_keys = ["wavelength"]
        tof_required_keys = ["sample_to_moderator_distance"]

        if d is None and t is None:
            return None

        joint = t.copy() if t else {}
        joint.update(d)

        check_for_required_keys(joint, basic_required_keys, "beam")
        try:
            check_for_required_keys(
                joint, monochromatic_required_keys, "monochromatic beam"
            )
            return MonochromaticBeam.from_dict(joint)
        except RuntimeError:
            check_for_required_keys(joint, tof_required_keys, "ToF beam")
            return TOFBeam.from_dict(joint)

    @staticmethod
    def make_tof_beam(sample_to_source_direction, sample_to_moderator_distance):
        return TOFBeam(
            tuple(map(float, sample_to_source_direction)),
            float(sample_to_moderator_distance),
        )

    @staticmethod
    def make_monochromatic_beam(
        sample_to_source=None,
        wavelength=None,
        s0=None,
        unit_s0=None,
        divergence=None,
        sigma_divergence=None,
    ):

        if divergence is None or sigma_divergence is None:
            divergence = 0.0
            sigma_divergence = 0.0

        if sample_to_source:
            assert wavelength
            return MonochromaticBeam(
                tuple(map(float, sample_to_source)),
                float(wavelength),
                float(divergence),
                float(sigma_divergence),
            )
        elif unit_s0:
            assert wavelength
            return MonochromaticBeam(
                tuple(-float(x) for x in unit_s0),
                float(wavelength),
                float(divergence),
                float(sigma_divergence),
            )
        else:
            assert s0
            return MonochromaticBeam(tuple(map(float, s0)))

    @staticmethod
    def make_polarized_beam(
        sample_to_source=None,
        wavelength=None,
        s0=None,
        unit_s0=None,
        polarization=None,
        polarization_fraction=None,
        divergence=None,
        sigma_divergence=None,
        flux=None,
        transmission=None,
    ):
        assert polarization
        assert 0.0 <= polarization_fraction <= 1.0

        if divergence is None or sigma_divergence is None:
            divergence = 0.0
            sigma_divergence = 0.0

        if flux is None:
            flux = 0
        if transmission is None:
            transmission = 1.0

        if sample_to_source:
            assert wavelength
            return MonochromaticBeam(
                tuple(map(float, sample_to_source)),
                float(wavelength),
                float(divergence),
                float(sigma_divergence),
                tuple(map(float, polarization)),
                float(polarization_fraction),
                float(flux),
                float(transmission),
            )
        elif unit_s0:
            assert wavelength
            return MonochromaticBeam(
                tuple(-float(x) for x in unit_s0),
                float(wavelength),
                float(divergence),
                float(sigma_divergence),
                tuple(map(float, polarization)),
                float(polarization_fraction),
                float(flux),
                float(transmission),
            )
        else:
            assert s0
            return MonochromaticBeam(
                tuple(map(float, s0)),
                float(divergence),
                float(sigma_divergence),
                tuple(map(float, polarization)),
                float(polarization_fraction),
                float(flux),
                float(transmission),
            )

    @staticmethod
    def simple_monochromatic(wavelength):
        """Construct a beam object on the principle that the beam is aligned
        with the +z axis, as is quite normal. Also assume the beam has
        polarization fraction 0.999 and is polarized in the x-z plane, unless
        it has a wavelength shorter than 0.05 Å in which case we assume
        electron diffraction and return an unpolarized beam model."""

        if wavelength > 0.05:
            return BeamFactory.make_monochromatic_beam(
                sample_to_source=(0.0, 0.0, 1.0), wavelength=wavelength
            )
        else:
            return BeamFactory.make_polarized_beam(
                sample_to_source=(0.0, 0.0, 1.0),
                wavelength=wavelength,
                polarization=(0, 1, 0),
                polarization_fraction=0.5,
            )

    @staticmethod
    def simple_directional_monochromatic(sample_to_source, wavelength):
        """Construct a beam with direction and wavelength."""

        if wavelength > 0.05:
            return BeamFactory.make_monochromatic_beam(
                sample_to_source=sample_to_source, wavelength=wavelength
            )
        else:
            return BeamFactory.make_polarized_beam(
                sample_to_source=sample_to_source,
                wavelength=wavelength,
                polarization=(0, 1, 0),
                polarization_fraction=0.5,
            )

    @staticmethod
    def complex_monochromatic(
        sample_to_source, polarization_fraction, polarization_plane_normal, wavelength
    ):
        """Full access to the constructor for cases where we do know everything
        that we need..."""

        return BeamFactory.make_polarized_beam(
            sample_to_source=sample_to_source,
            wavelength=wavelength,
            polarization=polarization_plane_normal,
            polarization_fraction=polarization_fraction,
        )

    @staticmethod
    def imgCIF(cif_file):
        """Initialize a detector model from an imgCIF file. N.B. the
        definition of the polarization plane is not completely helpful
        in this - it is the angle between the polarization plane and the
        +Y laboratory frame vector."""

        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_widefile(cif_file.encode(), pycbf.MSG_DIGEST)

        result = BeamFactory.imgCIF_H(cbf_handle)

        return result

    @staticmethod
    def imgCIF_H(cbf_handle):
        """Initialize a detector model from an imgCIF file. N.B. the
        definition of the polarization plane is not completely helpful
        in this - it is the angle between the polarization plane and the
        +Y laboratory frame vector. This example works from a cbf_handle,
        which is already configured."""

        d2r = math.pi / 180.0

        cbf_handle.find_category(b"axis")

        # find record with equipment = source

        try:
            cbf_handle.find_column(b"equipment")
            cbf_handle.find_row(b"source")

            # then get the vector and offset from this
            direction = []

            for j in range(3):
                cbf_handle.find_column(b"vector[%d]" % (j + 1))
                direction.append(cbf_handle.get_doublevalue())
        except Exception as e:
            if str(e).split()[-1] != "CBF_NOTFOUND":
                raise
            direction = [0, 0, 1]

        # and the wavelength
        wavelength = cbf_handle.get_wavelength()

        # and information about the polarization - FIXME this should probably
        # be a rotation about the beam not about the Z axis. Should also check
        # to see if this is Cu K-alpha wavelength (i.e. lab source...)

        try:
            polar_fraction, polar_angle = cbf_handle.get_polarization()
        except Exception:
            polar_fraction = 0.999
            polar_angle = 0.0

        polar_plane_normal = (
            math.sin(polar_angle * d2r),
            math.cos(polar_angle * d2r),
            0.0,
        )

        return BeamFactory.make_polarized_beam(
            sample_to_source=direction,
            wavelength=wavelength,
            polarization=polar_plane_normal,
            polarization_fraction=polar_fraction,
        )
