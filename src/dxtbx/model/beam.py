import math
from abc import ABC
from enum import Enum
from typing import Dict

import pycbf

import libtbx.phil

from dxtbx_model_ext import Beam, MonochromaticBeam, TOFBeam


class BeamType(Enum):
    MonochromaticBeam = 1
    TOFBeam = 2


beam_phil_scope = libtbx.phil.parse(
    """
  beam
    .expert_level = 1
    .short_caption = "Beam overrides"
  {
    type = "Monochromatic"
      .type = str
      .help = "Override the beam type"
      .short_caption = "beam type"

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


class BeamFactoryBase(ABC):
    """A factory class for beam objects, which encapsulate standard beam
    models. In cases where a full cbf description is available this
    will be used, otherwise simplified descriptions can be applied."""

    @staticmethod
    def from_phil(params, reference: Beam = None) -> Beam:
        """Convert the phil parameters into a beam model

        :param params: phil params
        :type params: [type]
        :param reference: beam model, defaults to None
        :type reference: Beam, optional
        :return: beam model
        :rtype: Beam
        """
        raise NotImplementedError

    @staticmethod
    def from_dict(dict: Dict, template: Dict) -> Beam:
        """Convert the dictionary to a beam model

        :param dict: beam parameters
        :type dict: Dictionary
        :param template: Starting parameters that can be overriden by dict
        :type template: Dictionary
        :return: beam model
        :rtype: Beam
        """
        raise NotImplementedError

    @staticmethod
    def make_beam(**kwargs) -> Beam:
        """Convert params into a beam model. Any missing params default to None.

        :return: beam model
        :rtype: Beam
        """
        raise NotImplementedError


class BeamFactory(BeamFactoryBase):

    """Class to select the correct factory to use based on the input.
    (i.e the class actually used in the code, unless a method unique to a
    derived class is required.)
    """

    @staticmethod
    def from_phil(params, reference: Beam = None) -> Beam:
        """Convert the phil parameters into a beam model

        :param params: phil params
        :type params: [type]
        :param reference: beam model, defaults to None
        :type reference: Beam, optional
        :return: beam model
        :rtype: Beam
        """

        if params.beam.type == "Monochromatic":
            return MonochromaticBeamFactory.from_phil(
                params=params, reference=reference
            )
        elif params.beam.type == "TOFBeam":
            return TOFBeamFactory.from_phil(params=params, reference=reference)
        else:
            raise NotImplementedError("Unknown beam type {params.beam.type}")

    @staticmethod
    def from_dict(dict: Dict, template: Dict = None) -> Beam:
        """Convert the dictionary to a beam model

        :param dict: Beam parameters
        :type dict: Dictionary
        :param template: Starting parameters that can be overriden by dict
        :type template: Dictionary
        :return: beam model
        :rtype: Beam
        """

        if "__id__" in dict:
            if template:
                assert (
                    dict["__id__"] == template["__id__"]
                ), "Beam and template dictionaries are not the same type."

            if "__id__" not in dict or dict["__id__"] == "MonochromaticBeam":
                return MonochromaticBeamFactory.from_dict(dict=dict, template=template)
            elif dict["__id__"] == "TOFBeam":
                return TOFBeamFactory.from_dict(dict=dict, template=template)
            else:
                raise NotImplementedError(f"Unknown beam type {dict['__id__']}")
        else:
            # Legacy beams without __id__ will all be monochromatic
            return MonochromaticBeamFactory.from_dict(dict=dict, template=template)

    @staticmethod
    def make_beam(beam_type: BeamType, **kwargs) -> Beam:
        """Convert params into a beam model. Any missing params default to None.

        :param beam_type: Which beam type to make
        :type beam_type: BeamType
        :return: beam model
        :rtype: Beam
        """

        beam_type = kwargs.get("beam_type")

        if beam_type == BeamType.MonochromaticBeam:
            return MonochromaticBeamFactory.make_beam(kwargs=kwargs)
        elif beam_type == BeamType.TOFBeam:
            return TOFBeamFactory.make_beam(kwargs=kwargs)
        else:
            raise NotImplementedError(f"Unknown beam type {beam_type}")


class MonochromaticBeamFactory(BeamFactoryBase):
    @staticmethod
    def from_phil(params, reference: Beam = None) -> Beam:
        def check_for_required_params(params, reference):
            if params.beam.direction is None and reference is None:
                raise RuntimeError("Cannot create MonochromaticBeam: direction not set")
            if params.beam.wavelength is None and reference is None:
                raise RuntimeError(
                    "Cannot create MonochromaticBeam: wavelength not set"
                )

        check_for_required_params(params=params, reference=reference)

        if reference is None:
            beam = MonochromaticBeam()
        else:
            beam = reference

        if params.beam.wavelength is not None:
            beam.set_wavelength(params.beam.wavelength)
        if params.beam.polarization_normal is not None:
            beam.set_polarization_normal(params.beam.polarization_normal)
        if params.beam.polarization_fraction is not None:
            beam.set_polarization_fraction(params.beam.polarization_fraction)
        return beam

    @staticmethod
    def from_dict(dict: Dict, template: Dict = None) -> Beam:
        def check_for_required_keys(dict, required_keys):
            for i in required_keys:
                if i not in dict:
                    raise RuntimeError(
                        f"Cannot create MonochromaticBeam: {i} not in dictionary"
                    )

        required_keys = ["direction", "wavelength"]
        check_for_required_keys(dict=dict, required_keys=required_keys)

        if dict is None and template is None:
            return None

        # Use the template as the initial dictionary,
        # and update/replace fields with dict
        beam_dict = template.copy() if template else {}
        beam_dict.update(dict)

        return MonochromaticBeam.from_dict(beam_dict)

    @staticmethod
    def make_beam(**kwargs) -> Beam:

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
    def make_simple_beam(wavelength):
        """Construct a beam object on the principle that the beam is aligned
        with the +z axis, as is quite normal. Also assume the beam has
        polarization fraction 0.999 and is polarized in the x-z plane, unless
        it has a wavelength shorter than 0.05 Å in which case we assume
        electron diffraction and return an unpolarized beam model."""

        if wavelength > 0.05:
            return MonochromaticBeamFactory.make_beam(
                sample_to_source=(0.0, 0.0, 1.0), wavelength=wavelength
            )
        else:
            return MonochromaticBeamFactory.make_polarized_beam(
                sample_to_source=(0.0, 0.0, 1.0),
                wavelength=wavelength,
                polarization=(0, 1, 0),
                polarization_fraction=0.5,
            )

    @staticmethod
    def make_simple_directional_beam(sample_to_source, wavelength):
        """Construct a beam with direction and wavelength."""

        if wavelength > 0.05:
            return MonochromaticBeamFactory.make_beam(
                sample_to_source=sample_to_source, wavelength=wavelength
            )
        else:
            return MonochromaticBeamFactory.make_polarized_beam(
                sample_to_source=sample_to_source,
                wavelength=wavelength,
                polarization=(0, 1, 0),
                polarization_fraction=0.5,
            )

    @staticmethod
    def make_complex_beam(
        sample_to_source, polarization_fraction, polarization_plane_normal, wavelength
    ):
        """Full access to the constructor for cases where we do know everything
        that we need..."""

        return MonochromaticBeamFactory.make_polarized_beam(
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

        result = MonochromaticBeamFactory.imgCIF_H(cbf_handle)

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

        return MonochromaticBeamFactory.make_polarized_beam(
            sample_to_source=direction,
            wavelength=wavelength,
            polarization=polar_plane_normal,
            polarization_fraction=polar_fraction,
        )


class TOFBeamFactory(BeamFactoryBase):
    @staticmethod
    def from_phil(params, reference: Beam = None) -> Beam:
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

        beam.set_sample_to_source_direction(params.beam.direction)
        beam.set_sample_to_moderator_distance(params.sample_to_moderator_distance)

        return beam

    @staticmethod
    def from_dict(dict: Dict, template: Dict = None) -> Beam:
        def check_for_required_keys(dict, required_keys):
            for i in required_keys:
                if i not in dict:
                    raise RuntimeError(f"Cannot create TOFBeam: {i} not in dictionary")

        required_keys = ["direction, sample_to_moderator_distance"]
        check_for_required_keys(dict=dict, required_keys=required_keys)
        if dict is None and template is None:
            return None

        # Use the template as the initial dictionary,
        # and update/replace fields with dict
        beam_dict = template.copy() if template else {}
        beam_dict.update(dict)

        return TOFBeam.from_dict(beam_dict)

    @staticmethod
    def make_beam(**kwargs) -> Beam:

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
