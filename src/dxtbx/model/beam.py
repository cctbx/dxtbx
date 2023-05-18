from __future__ import annotations

import math
from typing import Dict, Tuple, Union

import pycbf

import libtbx.phil

try:
    from ..dxtbx_model_ext import Beam, PolychromaticBeam
except ModuleNotFoundError:
    from dxtbx_model_ext import Beam, PolychromaticBeam  # type: ignore

Vec3Float = Tuple[float, float, float]

beam_phil_scope = libtbx.phil.parse(
    """
  beam
    .expert_level = 1
    .short_caption = "Beam overrides"
  {
    type = *monochromatic polychromatic
      .type = choice
      .help = "Override the beam type"
      .short_caption = "beam_type"

    wavelength = None
      .type = float
      .help = "Override the beam wavelength"

    direction = None
      .type = floats(size=3)
      .help = "Override the sample to source direction"
      .short_caption = "Sample to source direction"

    divergence = None
        .type = float
        .help = "Override the beam divergence"

    sigma_divergence = None
        .type = float
        .help = "Override the beam sigma divergence"

    polarization_normal = None
      .type = floats(size=3)
      .help = "Override the polarization normal"
      .short_caption = "Polarization normal"

    polarization_fraction = None
      .type = float(value_min=0.0, value_max=1.0)
      .help = "Override the polarization fraction"
      .short_caption = "Polarization fraction"

    transmission = None
        .type = float
        .help = "Override the transmission"
        .short_caption = "transmission"

    flux = None
        .type = float
        .help = "Override the flux"
        .short_caption = "flux"
  }
"""
)


class BeamFactory:
    """A factory class for beam objects, which encapsulate standard beam
    models. In cases where a full cbf description is available this
    will be used, otherwise simplified descriptions can be applied."""

    @staticmethod
    def from_phil(
        params: libtbx.phil.scope_extract,
        reference: Union[Beam, PolychromaticBeam] = None,
    ) -> Union[Beam, PolychromaticBeam]:
        """
        Convert the phil parameters into a beam model
        """

        # Check the input
        if reference is None:
            beam = (
                PolychromaticBeam() if params.beam.type == "polychromatic" else Beam()
            )
        else:
            beam = reference

        # Set the parameters
        if params.beam.type == "monochromatic":
            if params.beam.wavelength is not None:
                beam.set_wavelength(params.beam.wavelength)
            elif reference is None:
                raise RuntimeError("No wavelength set")
        if params.beam.direction is not None:
            beam.set_direction(params.beam.direction)
        elif reference is None:
            raise RuntimeError("No beam direction set")

        if params.beam.divergence is not None:
            beam.set_divergence(params.beam.divergence)
        if params.beam.sigma_divergence is not None:
            beam.set_sigma_divergence(params.beam.sigma_divergence)
        if params.beam.polarization_normal is not None:
            beam.set_polarization_normal(params.beam.polarization_normal)
        if params.beam.polarization_fraction is not None:
            beam.set_polarization_fraction(params.beam.polarization_fraction)
        if params.beam.transmission is not None:
            beam.set_transmission(params.beam.transmission)
        if params.beam.flux is not None:
            beam.set_flux(params.beam.flux)

        return beam

    @staticmethod
    def from_dict(dict: Dict, template: Dict = None) -> Union[Beam, PolychromaticBeam]:
        """Convert the dictionary to a beam model"""

        if template is not None:
            if "__id__" in dict and "__id__" in template:
                assert (
                    dict["__id__"] == template["__id__"]
                ), "Beam and template dictionaries are not the same type."

        if dict is None and template is None:
            return None
        joint = template.copy() if template else {}
        joint.update(dict)

        # Create the model from the joint dictionary
        if joint.get("__id__", "") == "polychromatic":
            return PolychromaticBeam.from_dict(joint)
        return Beam.from_dict(joint)

    @staticmethod
    def make_beam(
        sample_to_source: Vec3Float = None,
        wavelength: float = None,
        s0: Vec3Float = None,
        unit_s0: Vec3Float = None,
        divergence: float = None,
        sigma_divergence: float = None,
    ) -> Beam:

        if divergence is None or sigma_divergence is None:
            divergence = 0.0
            sigma_divergence = 0.0

        if sample_to_source:
            assert wavelength
            return Beam(
                tuple(map(float, sample_to_source)),
                float(wavelength),
                float(divergence),
                float(sigma_divergence),
            )
        elif unit_s0:
            assert wavelength
            return Beam(
                tuple(-float(x) for x in unit_s0),
                float(wavelength),
                float(divergence),
                float(sigma_divergence),
            )
        else:
            assert s0
            return Beam(tuple(map(float, s0)))

    @staticmethod
    def make_polychromatic_beam(
        direction: Vec3Float,
        divergence: float = 0.0,
        sigma_divergence: float = 0.0,
        polarization_normal: Vec3Float = (0.0, 1.0, 0.0),
        polarization_fraction: float = 0.5,
        flux: float = 0.0,
        transmission: float = 1.0,
        deg: bool = True,
    ) -> PolychromaticBeam:

        return PolychromaticBeam(
            tuple(map(float, direction)),
            float(divergence),
            float(sigma_divergence),
            tuple(map(float, polarization_normal)),
            float(polarization_fraction),
            float(flux),
            float(transmission),
            bool(deg),
        )

    @staticmethod
    def make_polarized_beam(
        sample_to_source: Vec3Float = None,
        wavelength: float = None,
        s0: Vec3Float = None,
        unit_s0: Vec3Float = None,
        polarization: Vec3Float = None,
        polarization_fraction: float = None,
        divergence: float = None,
        sigma_divergence: float = None,
        flux: float = None,
        transmission: float = None,
    ) -> Beam:
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
            return Beam(
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
            return Beam(
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
            return Beam(
                tuple(map(float, s0)),
                float(divergence),
                float(sigma_divergence),
                tuple(map(float, polarization)),
                float(polarization_fraction),
                float(flux),
                float(transmission),
            )

    @staticmethod
    def simple(wavelength: float) -> Beam:
        """Construct a beam object on the principle that the beam is aligned
        with the +z axis, as is quite normal. Also assume the beam has
        polarization fraction 0.999 and is polarized in the x-z plane, unless
        it has a wavelength shorter than 0.05 Å in which case we assume
        electron diffraction and return an unpolarized beam model."""

        if wavelength > 0.05:
            return BeamFactory.make_beam(
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
    def simple_directional(sample_to_source: Vec3Float, wavelength: float) -> Beam:
        """Construct a beam with direction and wavelength."""

        if wavelength > 0.05:
            return BeamFactory.make_beam(
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
    def complex(
        sample_to_source: Vec3Float,
        polarization_fraction: float,
        polarization_plane_normal: Vec3Float,
        wavelength: float,
    ) -> Beam:
        """Full access to the constructor for cases where we do know everything
        that we need..."""

        return BeamFactory.make_polarized_beam(
            sample_to_source=sample_to_source,
            wavelength=wavelength,
            polarization=polarization_plane_normal,
            polarization_fraction=polarization_fraction,
        )

    @staticmethod
    def imgCIF(cif_file: str) -> Beam:
        """Initialize a detector model from an imgCIF file. N.B. the
        definition of the polarization plane is not completely helpful
        in this - it is the angle between the polarization plane and the
        +Y laboratory frame vector."""

        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_widefile(cif_file.encode(), pycbf.MSG_DIGEST)

        result = BeamFactory.imgCIF_H(cbf_handle)

        return result

    @staticmethod
    def imgCIF_H(cbf_handle: pycbf.cbf_handle_struct) -> Beam:
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

        # and the flux if available
        try:
            cbf_handle.find_category(b"diffrn_radiation")
            cbf_handle.find_column(b"beam_flux")
            flux = cbf_handle.get_value()
        except Exception as e:
            if str(e).split()[-1] != "CBF_NOTFOUND":
                raise
            flux = None

        return BeamFactory.make_polarized_beam(
            sample_to_source=direction,
            wavelength=wavelength,
            polarization=polar_plane_normal,
            polarization_fraction=polar_fraction,
            flux=flux,
        )
