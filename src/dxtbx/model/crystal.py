from __future__ import annotations

import iotbx.phil
from cctbx.sgtbx import space_group as SG
from scitbx import matrix

try:
    from ..dxtbx_model_ext import (
        Crystal,
        MosaicCrystalKabsch2010,
        MosaicCrystalSauter2014,
    )
except ModuleNotFoundError:
    from dxtbx_model_ext import (  # type: ignore
        Crystal,
        MosaicCrystalKabsch2010,
        MosaicCrystalSauter2014,
    )

crystal_phil_scope = iotbx.phil.parse(
    """
  crystal
    .expert_level = 2
    .short_caption = "Crystal overrides"
  {
    unit_cell = None
      .type = unit_cell

    A_matrix = None
      .type = floats(size=9)
      .help = "The crystal setting A=UB matrix. If set, this will override the unit cell."

    space_group = None
      .type = space_group
  }
"""
)


class CrystalFactory:
    @staticmethod
    def from_dict(d, t=None):
        """Convert the dictionary to a crystal model

        Params:
            d The dictionary of parameters
            t The template dictionary to use

        Returns:
            The crystal model

        """
        if d is None and t is None:
            return None
        joint = t.copy() if t else {}
        joint.update(d)

        # Create the model from the dictionary
        if "ML_half_mosaicity_deg" in joint:
            assert "ML_domain_size_ang" in joint
            if (
                joint["ML_half_mosaicity_deg"] is None
                or joint["ML_domain_size_ang"] is None
            ):
                assert (
                    joint["ML_half_mosaicity_deg"] is None
                    and joint["ML_domain_size_ang"] is None
                )
            else:
                if joint.get("mosaicity", 0) > 0:
                    print(
                        "Warning, two kinds of mosaicity found. Using Sauter2014 model"
                    )
                return MosaicCrystalSauter2014.from_dict(joint)
        if "mosaicity" in joint:
            return MosaicCrystalKabsch2010.from_dict(joint)
        else:
            return Crystal.from_dict(joint)

    @staticmethod
    def from_mosflm_matrix(
        mosflm_A_matrix, unit_cell=None, wavelength=None, space_group=None
    ):
        """
        Create a crystal_model from a Mosflm A matrix (a*, b*, c*); N.B. assumes
        the mosflm coordinate frame::

                                                         /!
                            Y-axis                      / !
                              ^                        /  !
                              !                       /   !
                              !                      /    !
                              !   /                 /  Xd !
                              !  /                 / * ^  !
                              ! /                  ! 3 !  !
                              !/      X-ray beam   !   !  !
                              /------------------------/--!---->X-axis
                             /                     !  / *1!
                          <-/-                     ! /    !
                           / \\+ve phi             ! Yd  /
                          /   /                    ! 2  /
                         /                         ! * /
                        Z-axis                  Ys ^ _/
                      Rotation                     ! /| Xs
                       axis                        !/
                                                   O

        Also assume that the mosaic spread is 0.

        :param mosflm_A_matrix: The A matrix in Mosflm convention.
        :type mosflm_A_matrix: tuple of floats
        :param unit_cell: The unit cell parameters which are used to determine the
                          wavelength from the Mosflm A matrix.
        :type unit_cell: cctbx.uctbx.unit_cell
        :param wavelength: The wavelength to scale the A matrix
        :type wavelength: float
        :param space_group: If the space group is None then the space_group will
                            be assigned as P1
        :type space_group: cctbx.sgtbx.space_group
        :returns: A crystal model derived from the given Mosflm A matrix
        :rtype: :py:class:`crystal_model`
        """
        if not space_group:
            space_group = SG("P1")

        A_star = matrix.sqr(mosflm_A_matrix)
        A = A_star.inverse()

        if unit_cell:
            unit_cell_constants = unit_cell.parameters()
            a = matrix.col(A.elems[0:3])
            wavelength = unit_cell_constants[0] / a.length()
            A *= wavelength
        elif wavelength:
            A *= wavelength
        else:
            # assume user has pre-scaled
            pass

        a = A.elems[0:3]
        b = A.elems[3:6]
        c = A.elems[6:9]
        rotate_mosflm_to_imgCIF = matrix.sqr((0, 0, 1, 0, 1, 0, -1, 0, 0))
        _a = rotate_mosflm_to_imgCIF * a
        _b = rotate_mosflm_to_imgCIF * b
        _c = rotate_mosflm_to_imgCIF * c

        return Crystal(_a, _b, _c, space_group=space_group)

    @staticmethod
    def from_phil(
        params: iotbx.phil.scope_extract,
        reference: Crystal | None = None,
    ) -> Crystal:
        """
        Convert the phil parameters into a crystal model
        """

        all_params = [
            params.crystal.unit_cell,
            params.crystal.A_matrix,
            params.crystal.space_group,
        ]
        if all_params.count(None) == len(all_params):
            return reference

        if reference is None:
            crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), "P1")
        else:
            crystal = reference
            crystal.reset_scan_points()

        if params.crystal.unit_cell is not None:
            crystal.set_unit_cell(params.crystal.unit_cell)
        if params.crystal.A_matrix is not None:
            crystal.set_A(params.crystal.A_matrix)
        if params.crystal.space_group is not None:
            crystal.set_space_group(params.crystal.space_group.group())

        return crystal
