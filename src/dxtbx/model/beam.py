from __future__ import annotations

import math
from abc import ABC, abstractmethod
from typing import Dict, Tuple

import pycbf

import libtbx.phil

try:
    from ..dxtbx_model_ext import Beam, BeamBase, TOFBeam
except ModuleNotFoundError:
    from dxtbx_model_ext import BeamBase, Beam, TOFBeam  # type: ignore

from enum import Enum


class BeamType(Enum):
    Monochromatic = 1
    TOF = 2


Vec3Float = Tuple[float, float, float]

beam_phil_scope = libtbx.phil.parse(
    """
  beam
    .expert_level = 1
    .short_caption = "Beam overrides"
  {
    type = *monochromatic tof
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

    sample_to_moderator_distance = None
        .type = float
        .help = "Override sample to moderator distance"
        .short_caption = "Sample to moderator distance"
  }
"""
)


class AbstractBeamFactory(ABC):

    """
    Factory interface for all beam factories.
    """

    @staticmethod
    @abstractmethod
    def from_phil(params, reference: BeamBase = None) -> BeamBase:
        """Convert the phil parameters into a beam model"""
        ...

    @staticmethod
    @abstractmethod
    def from_dict(dict: Dict, template: Dict) -> BeamBase:
        """Convert the dictionary to a beam model"""
        ...

    @staticmethod
    @abstractmethod
    def make_beam(**kwargs) -> BeamBase:
        """Convert params into a beam model"""
        ...


class BeamBaseFactory(AbstractBeamFactory):

    """
    Factory for selecting which beam factory to use.
    """

    @staticmethod
    def from_phil(params, reference: BeamBase = None) -> BeamBase:
        """Convert the phil parameters into a beam model"""

        if params.beam.type == "tof":
            return TOFBeamFactory.from_phil(params=params, reference=reference)
        else:  # Default to monochromatic for back compatibility
            return BeamFactory.from_phil(params=params, reference=reference)

    @staticmethod
    def from_dict(dict: Dict, template: Dict = None) -> BeamBase:
        """Convert the dictionary to a beam model"""

        if template is not None:
            if "__id__" in dict and "__id__" in template:
                assert (
                    dict["__id__"] == template["__id__"]
                ), "Beam and template dictionaries are not the same type."

        # Assume dictionaries without "__id__" are for Beam objects
        if "__id__" not in dict or dict["__id__"] == "Monochromatic":
            return BeamFactory.from_dict(dict=dict, template=template)
        elif dict["__id__"] == "TOF":
            return TOFBeamFactory.from_dict(dict=dict, template=template)
        else:
            raise NotImplementedError(f"Unknown beam type {dict['__id__']}")

    @staticmethod
    def make_beam(**kwargs) -> BeamBase:

        """
        Convert params into a beam model. Any missing params default to None.
        """

        beam_type = kwargs.get("beam_type")

        # Default to Beam for back compatability
        if beam_type is None or beam_type == BeamType.Monochromatic:
            return BeamFactory.make_beam(**kwargs)
        elif beam_type == BeamType.TOF:
            return TOFBeamFactory.make_beam(**kwargs)
        else:
            raise NotImplementedError(f"Unknown beam type {beam_type}")


class BeamFactory(AbstractBeamFactory):

    """
    A factory class for Beam objects, which encapsulate standard
    monochromatic beam models. In cases where a full cbf description
    is available this will be used, otherwise simplified descriptions
    can be applied.
    """

    @staticmethod
    def from_phil(params, reference: Beam = None) -> Beam:

        """
        Convert the phil parameters into a beam model
        """

        # Check the input
        if reference is None:
            beam = Beam()
        else:
            beam = reference

        # Set the parameters
        if params.beam.wavelength is not None:
            beam.set_wavelength(params.beam.wavelength)
        elif reference is None:
            raise RuntimeError("No wavelength set")
        if params.beam.direction is not None:
            beam.set_direction(params.beam.direction)
        elif reference is None:
            raise RuntimeError("No beam direction set")
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
    def from_dict(dict: Dict, template: Dict = None) -> Beam:

        """Convert the dictionary to a beam model"""

        if dict is None and template is None:
            return None
        joint = template.copy() if template else {}
        joint.update(dict)

        # Create the model from the joint dictionary
        return Beam.from_dict(joint)

    @staticmethod
    def make_beam(**kwargs) -> Beam:

        """Convert params into a beam model"""

        sample_to_source = kwargs.get("sample_to_source")
        wavelength = kwargs.get("wavelength")
        s0 = kwargs.get("s0")
        unit_s0 = kwargs.get("unit_s0")
        divergence = kwargs.get("divergence")
        sigma_divergence = kwargs.get("sigma_divergence")

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

        """
        Construct a beam object on the principle that the beam is aligned
        with the +z axis, as is quite normal. Also assume the beam has
        polarization fraction 0.999 and is polarized in the x-z plane, unless
        it has a wavelength shorter than 0.05 Å in which case we assume
        electron diffraction and return an unpolarized beam model.
        """

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

        """
        Full access to the constructor for cases where we do know everything
        that we need.
        """

        return BeamFactory.make_polarized_beam(
            sample_to_source=sample_to_source,
            wavelength=wavelength,
            polarization=polarization_plane_normal,
            polarization_fraction=polarization_fraction,
        )

    @staticmethod
    def imgCIF(cif_file: str) -> Beam:

        """
        Initialize a Beam model from an imgCIF file. N.B. the
        definition of the polarization plane is not completely helpful
        in this - it is the angle between the polarization plane and the
        +Y laboratory frame vector.
        """

        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_widefile(cif_file.encode(), pycbf.MSG_DIGEST)

        result = BeamFactory.imgCIF_H(cbf_handle)

        return result

    @staticmethod
    def imgCIF_H(cbf_handle: pycbf.cbf_handle_struct) -> Beam:

        """
        Initialize a Beam model from an imgCIF file. N.B. the
        definition of the polarization plane is not completely helpful
        in this - it is the angle between the polarization plane and the
        +Y laboratory frame vector. This example works from a cbf_handle,
        which is already configured.
        """

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


class TOFBeamFactory(AbstractBeamFactory):

    """
    A factory for creating TOFBeam objects.
    """

    @staticmethod
    def from_phil(params, reference: TOFBeam = None) -> TOFBeam:

        """
        Convert the phil parameters into a TOFBeam model
        """

        def check_for_required_params(params, reference):
            if params.beam.direction is None and reference is None:
                raise RuntimeError("Cannot create TOFBeam: direction not set")
            if params.beam.sample_to_moderator_distance is None and reference is None:
                raise RuntimeError(
                    "Cannot create ToF beam: sample_to_moderator_distance not set"
                )

        check_for_required_params(params=params, reference=reference)

        if reference is None:
            beam = TOFBeam()
        else:
            beam = reference

        if params.beam.direction is not None:
            beam.set_direction(params.beam.direction)
        if params.beam.sample_to_moderator_distance is not None:
            beam.set_sample_to_moderator_distance(
                params.beam.sample_to_moderator_distance
            )
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
    def from_dict(dict: Dict, template: Dict = None) -> TOFBeam:

        """Convert the dictionary to a TOFBeam model"""

        def check_for_required_keys(dict, required_keys):
            for i in required_keys:
                if i not in dict:
                    raise RuntimeError(f"Cannot create TOFBeam: {i} not in dictionary")

        required_keys = [
            "direction",
            "sample_to_moderator_distance",
        ]
        if dict is None and template is None:
            return None
        # Use the template as the initial dictionary,
        # and update/replace fields with dict
        beam_dict = template.copy() if template else {}
        beam_dict.update(dict)
        check_for_required_keys(dict=beam_dict, required_keys=required_keys)

        return TOFBeam.from_dict(beam_dict)

    @staticmethod
    def make_beam(**kwargs) -> TOFBeam:

        """Convert params into a TOFBeam model"""

        sample_to_source_direction = kwargs.get("sample_to_source_direction")
        if not sample_to_source_direction:
            raise RuntimeError(
                "Cannot create TOFBeam: sample_to_source_direction not set"
            )

        sample_to_moderator_distance = kwargs.get("sample_to_moderator_distance")
        if not sample_to_moderator_distance:
            raise RuntimeError(
                "Cannot create TOFBeam: sample_to_moderator_distance not set"
            )

        return TOFBeam(
            tuple(map(float, sample_to_source_direction)),
            float(sample_to_moderator_distance),
        )
