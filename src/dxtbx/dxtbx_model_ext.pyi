from __future__ import annotations

from enum import Enum
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
    overload,
)

from cctbx.sgtbx import space_group
from cctbx.uctbx import unit_cell
from scitbx.array_family import flex
from scitbx.array_family import shared as flex_shared

# Attempt to use the stub typing for flex-inheritance
from scitbx.array_family.flex import FlexPlain

# TypeVar for the set of Experiment models that can be joint-accepted
# - profile, imageset and scalingmodel are handled as 'object'
TExperimentModel = TypeVar(
    "TExperimentModel", BeamBase, Detector, Goniometer, Scan, CrystalBase, object
)

Vec2Float = Tuple[float, float]
Vec3Float = Tuple[float, float, float]
Vec6Float = Tuple[float, float, float, float, float, float]
Vec9Float = Tuple[float, float, float, float, float, float, float, float, float]
Vec2Int = Tuple[int, int]
Vec4Int = Tuple[int, int, int, int]

class BeamBase:
    @property
    def num_scan_points(self) -> int: ...
    def get_divergence(self) -> float: ...
    def set_divergence(self, divergence: float) -> None: ...
    def get_flux(self) -> float: ...
    def set_flux(self, flux: float) -> None: ...
    def get_num_scan_points(self) -> int: ...
    def get_polarization_fraction(self) -> float: ...
    def set_polarization_fraction(self, polarization_fraction: float) -> None: ...
    def get_polarization_normal(self) -> Vec3Float: ...
    def set_polarization_normal(self, polarization_normal: Vec3Float) -> None: ...
    def get_s0(self) -> Vec3Float: ...
    def set_s0(self, s0: Vec3Float) -> None: ...
    def get_s0_at_scan_point(self, index: int) -> Vec3Float: ...
    def get_s0_at_scan_points(self) -> flex.vec3_double: ...
    def set_s0_at_scan_points(
        self, points: Union[Tuple[Vec3Float], List[Vec3Float]]
    ) -> None: ...
    def get_sample_to_source_direction(self) -> Vec3Float: ...
    def get_sigma_divergence(self) -> float: ...
    def set_sigma_divergence(self, sigma_divergence: float) -> None: ...
    def get_transmission(self) -> float: ...
    def set_transmission(self, transmission: float) -> None: ...
    def get_unit_s0(self) -> Vec3Float: ...
    def set_unit_s0(self, unit_s0: Vec3Float) -> None: ...
    def get_wavelength(self) -> float: ...
    def set_wavelength(self, wavelength: float) -> None: ...
    def is_similar_to(
        self,
        other: BeamBase,
        wavelength_tolerance: float = ...,
        direction_tolerance: float = ...,
        polarization_normal_tolerance: float = ...,
        polarization_fraction_tolerance: float = ...,
    ) -> bool: ...
    def reset_scan_points(self) -> None: ...
    def rotate_around_origin(
        self, axis: Vec3Float, angle: float, deg: bool = ...
    ) -> Any: ...
    def set_direction(self, direction: Vec3Float) -> None: ...

class Beam(BeamBase):
    @overload
    def __init__(self, beam: Beam) -> None: ...
    @overload
    def __init__(self, direction: Vec3Float, wavelength: float) -> None: ...
    @overload
    def __init__(self, s0: Vec3Float) -> None: ...
    @overload
    def __init__(
        self,
        direction: Vec3Float,
        wavelength: float,
        divergence: float,
        sigma_divergence: float,
        deg: bool = ...,
    ) -> None: ...
    @overload
    def __init__(
        self, s0: Vec3Float, divergence: float, sigma_divergence: float, deg: bool = ...
    ) -> None: ...
    @overload
    def __init__(
        self,
        direction: Vec3Float,
        wavelength: float,
        divergence: float,
        sigma_divergence: float,
        polarization_normal: Vec3Float,
        polarization_fraction: float,
        flux: float,
        transmission: float,
        deg: bool = ...,
    ) -> None: ...
    @staticmethod
    def from_dict(data: Dict) -> Beam: ...
    def to_dict(self) -> Dict: ...

class CrystalBase:
    @property
    def num_scan_points(self) -> int: ...
    def change_basis(self, cctbx) -> Any: ...
    def get_A(self) -> Vec9Float: ...
    def get_B(self) -> Vec9Float: ...
    def get_U(self) -> Vec9Float: ...
    def set_A(self, A: Vec9Float) -> None: ...
    def set_B(self, B: Vec9Float) -> None: ...
    def set_U(self, U: Vec9Float) -> None: ...
    def get_A_at_scan_point(self, index: int) -> Vec9Float: ...
    def get_B_at_scan_point(self, index: int) -> Vec9Float: ...
    def get_B_covariance(self) -> flex.double: ...
    def get_B_covariance_at_scan_point(self, index: int) -> flex.double: ...
    def get_B_covariance_at_scan_points(self) -> flex.double: ...
    def get_U_at_scan_point(self, index: int) -> Vec9Float: ...
    def get_cell_parameter_sd(self) -> Vec6Float: ...
    def get_cell_parameter_sd_at_scan_point(self, index: int) -> Vec6Float: ...
    def get_cell_volume_sd(self) -> float: ...
    def get_num_scan_points(self) -> int: ...
    def get_real_space_vectors(self) -> flex.vec3_double: ...
    def get_recalculated_cell_parameter_sd(self) -> Vec6Float: ...
    def get_recalculated_cell_volume_sd(self) -> float: ...
    def get_recalculated_unit_cell(self) -> Optional[unit_cell]: ...
    def get_unit_cell(self) -> unit_cell: ...
    def get_space_group(self) -> space_group: ...
    def get_unit_cell_at_scan_point(self, index: int) -> unit_cell: ...
    def is_similar_to(
        self,
        other: CrystalBase,
        angle_tolerance: float = ...,
        uc_rel_length_tolerance: float = ...,
        uc_abs_angle_tolerance: float = ...,
    ) -> bool: ...
    def reset_scan_points(self) -> None: ...
    def reset_unit_cell_errors(self) -> None: ...
    def rotate_around_origin(
        self, axis: Vec3Float, angle: float, deg: bool = ...
    ) -> Any: ...
    def set_A_at_scan_points(
        self, value: Union[List[Vec9Float], Tuple[Vec9Float], flex.mat3_double]
    ) -> None: ...
    def set_B_covariance(
        self, covariance: Union[Tuple[float, ...], flex.double]
    ) -> None: ...
    def set_B_covariance_at_scan_points(self, data: flex.double) -> None: ...
    def set_recalculated_cell_parameter_sd(self, value: Vec6Float) -> None: ...
    def set_recalculated_cell_volume_sd(self, volume: float) -> None: ...
    def set_recalculated_unit_cell(self, unit_cell: unit_cell) -> None: ...
    def set_space_group(self, space_group: space_group) -> None: ...
    def set_unit_cell(self, unit_cell: unit_cell) -> None: ...
    def update(self, other: CrystalBase) -> None: ...
    def update_B(self) -> None: ...

class Crystal(CrystalBase):
    @overload
    def __init__(self, other: Crystal) -> None: ...
    @overload
    def __init__(
        self,
        real_space_a: float,
        real_space_b: float,
        real_space_c: float,
        space_group: space_group,
    ) -> None: ...
    @overload
    def __init__(
        self,
        real_space_a: float,
        real_space_b: float,
        real_space_c: float,
        space_group_symbol: str,
    ) -> None: ...
    @overload
    def __init__(
        self, A: Vec9Float, space_group: space_group, reciprocal: bool = ...
    ) -> None: ...
    @overload
    def __init__(
        self, A: Vec9Float, space_group_symbol: str, reciprocal: bool = ...
    ) -> None: ...

class MosaicCrystalKabsch2010(Crystal):
    def is_similar_to(
        self,
        other: CrystalBase,
        angle_tolerance: float = ...,
        uc_rel_length_tolerance: float = ...,
        uc_abs_angle_tolerance: float = ...,
        mosaicity_tolerance: float = ...,
    ) -> bool: ...
    def get_mosaicity(self, deg: bool = True) -> float: ...
    def set_mosaicity(self, mosaicity: float, deg: bool = True) -> Any: ...

class MosaicCrystalSauter2014(Crystal):
    def is_similar_to(
        self,
        other: CrystalBase,
        angle_tolerance: float = ...,
        uc_rel_length_tolerance: float = ...,
        uc_abs_angle_tolerance: float = ...,
        half_mosaicity_tolerance: float = ...,
        domain_size_tolerance: float = ...,
    ) -> bool: ...
    def get_half_mosaicity_deg(self) -> float: ...
    def set_half_mosaicity_deg(self, half_mosaicity_deg: float) -> None: ...
    def get_domain_size_ang(self) -> float: ...
    def set_domain_size_ang(self, domain_size_ang: float) -> None: ...

class Detector:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, panel: Panel) -> None: ...
    def hierarchy(self) -> DetectorNode: ...
    @overload
    def add_group(self) -> DetectorNode: ...
    @overload
    def add_group(self, group: Panel) -> DetectorNode: ...
    @overload
    def add_panel(self) -> Panel: ...
    @overload
    def add_panel(self, panel: Panel) -> Panel: ...
    @staticmethod
    def from_dict(data: Dict) -> Detector: ...
    def to_dict(self) -> Dict: ...
    def get_max_inscribed_resolution(self, s0: Vec3Float) -> float: ...
    def get_max_resolution(self, s0: Vec3Float) -> float: ...
    def get_names(self) -> flex.std_string: ...
    def get_panel_intersection(self, s1: Vec3Float) -> int: ...
    def get_ray_intersection(self, s1: Vec3Float) -> Tuple[int, Vec2Float]: ...
    def has_projection_2d(self) -> bool: ...
    def is_similar_to(
        self,
        other: Detector,
        fast_axis_tolerance: float = ...,
        slow_axis_tolerance: float = ...,
        origin_tolerance: float = ...,
        static_only: bool = ...,
        ignore_trusted_range: bool = ...,
    ) -> bool: ...
    def rotate_around_origin(
        self, axis: Vec3Float, angle: float, deg: bool = ...
    ) -> Any: ...
    def __getitem__(self, index: int) -> Panel: ...
    def __setitem__(self, index: int, panel: Panel) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(boost) -> Iterator[DetectorNode]: ...

class Experiment:
    beam: BeamBase
    crystal: CrystalBase
    detector: Detector
    goniometer: Goniometer
    identifier: str
    imageset: object
    profile: object
    scaling_model: object
    scan: Scan
    def __init__(
        self,
        beam: BeamBase,
        detector: Detector,
        goniometer: Goniometer,
        scan: Scan,
        crystal: CrystalBase,
        profile: object,
        imageset: object,
        scaling_model: object,
        identifier: str = ...,
    ) -> None: ...
    def is_consistent(self) -> bool: ...
    def is_sequence(self) -> bool: ...
    def is_still(self) -> bool: ...
    def __contains__(self, obj: TExperimentModel) -> bool: ...

class ExperimentList:
    @overload
    def __init__(self, experiments: Sequence[Experiment]) -> None: ...
    @overload
    def __init__(self) -> None: ...
    def identifiers(self) -> flex.std_string: ...
    def find(self, identifier: str) -> int: ...
    def append(self, experiment: Experiment) -> None: ...
    def extend(self, other: ExperimentList) -> None: ...
    def clear(self) -> None: ...
    def empty(self) -> bool: ...
    @overload
    def __getitem__(self, index: int) -> Experiment: ...
    @overload
    def __getitem__(self, indices: slice) -> ExperimentList: ...
    def __setitem__(self, index: int, experiment: Experiment) -> None: ...
    def __delitem__(self, index: int) -> None: ...
    def __iter__(self) -> Iterator[Experiment]: ...
    def __contains__(self, obj: TExperimentModel) -> bool: ...
    def replace(self, obj: TExperimentModel, withobj: TExperimentModel) -> None: ...
    def indices(
        self, obj: Union[TExperimentModel, FlexPlain[TExperimentModel]]
    ) -> flex.size_t: ...
    def remove_on_experiment_identifiers(self, identifiers: List[str]) -> None: ...
    def select_on_experiment_identifiers(self, identifiers: List[str]) -> None: ...
    def where(
        self,
        beam: BeamBase,
        detector: Detector,
        goniometer: Goniometer,
        scan: Scan,
        profile: object,
        imageset: object,
        scaling_model: object,
    ) -> flex.size_t: ...
    def is_consistent(self) -> bool: ...
    def __len__(self) -> int: ...

class GoniometerBase:
    pass

class Goniometer(GoniometerBase):
    @overload
    def __init__(self, other: Goniometer) -> None: ...
    @overload
    def __init__(self, rotation_axis: Vec3Float) -> None: ...
    @overload
    def __init__(
        self, rotation_axis: Vec3Float, fixed_rotation_matrix: Vec9Float
    ) -> None: ...
    @overload
    def __init__(
        self,
        rotation_axis: Vec3Float,
        fixed_rotation_matrix: Vec9Float,
        setting_rotation_matrix: Vec9Float,
    ) -> None: ...
    @property
    def num_scan_points(self) -> int: ...
    @staticmethod
    def from_dict(data: Dict) -> Goniometer: ...
    def get_fixed_rotation(self) -> Vec9Float: ...
    def get_num_scan_points(self) -> int: ...
    def get_rotation_axis(self) -> Vec3Float: ...
    def get_rotation_axis_datum(self) -> Vec3Float: ...
    def get_setting_rotation(self) -> Vec9Float: ...
    def get_setting_rotation_at_scan_point(self, index: int) -> Vec9Float: ...
    def get_setting_rotation_at_scan_points(self) -> flex.mat3_double: ...
    def is_similar_to(
        self,
        other: Goniometer,
        rotation_axis_tolerance: float = ...,
        fixed_rotation_tolerance: float = ...,
        setting_rotation_tolerance: float = ...,
    ) -> bool: ...
    def reset_scan_points(self) -> None: ...
    def rotate_around_origin(
        self, axis: Vec3Float, angle: float, deg: bool = ...
    ) -> None: ...
    def set_fixed_rotation(self, rotation: Vec9Float) -> None: ...
    def set_rotation_axis(self, axis: Vec3Float) -> None: ...
    def set_rotation_axis_datum(self, axis_datum: Vec3Float) -> None: ...
    def set_setting_rotation(self, rotation: Vec9Float) -> None: ...
    def set_setting_rotation_at_scan_points(
        self,
        setting_rotations: Union[List[Vec9Float], Tuple[Vec9Float], flex.mat3_double],
    ) -> Any: ...
    def to_dict(self) -> Dict: ...

class KappaDirection(Enum):
    PlusY = ...
    PlusZ = ...
    MinusY = ...
    MinusZ = ...

class KappaScanAxis(Enum):
    Omega = ...
    Phi = ...

class KappaGoniometer(Goniometer):
    @overload
    def __init__(
        self,
        alpha: float,
        omega: float,
        kappa: float,
        phi: float,
        direction: KappaDirection,
        scan_axis: KappaScanAxis,
    ) -> None: ...
    @overload
    def __init__(
        self,
        alpha: float,
        omega: float,
        kappa: float,
        phi: float,
        direction: str,
        scan_axis: str,
    ) -> None: ...
    def get_alpha_angle(self) -> float: ...
    def get_direction(self) -> KappaDirection: ...
    def get_kappa_angle(self) -> float: ...
    def get_kappa_axis(self) -> Vec3Float: ...
    def get_omega_angle(self) -> float: ...
    def get_omega_axis(self) -> Vec3Float: ...
    def get_phi_angle(self) -> float: ...
    def get_phi_axis(self) -> Vec3Float: ...
    def get_scan_axis(self) -> Vec3Float: ...

class MultiAxisGoniometer(Goniometer):
    def __init__(
        self,
        axes: flex.vec3_double,
        angles: flex.double,
        names: flex.std_string,
        scan_axis: int,
    ) -> None: ...
    @staticmethod
    def from_dict(data: Dict) -> MultiAxisGoniometer: ...
    def to_dict(self) -> Dict: ...
    def get_angles(self) -> flex.double: ...
    def set_angles(self, angles: flex.double) -> None: ...
    def get_axes(self) -> flex.vec3_double: ...
    def set_axes(self, axes: flex.vec3_double) -> None: ...
    def get_names(self) -> flex.std_string: ...
    def get_scan_axis(self) -> int: ...

class VirtualPanelFrame:
    def get_beam_centre_lab(self, s0: Vec3Float) -> Vec3Float: ...
    def get_beam_centre(self, s0: Vec3Float) -> Vec2Float: ...
    def get_bidirectional_ray_intersection(self, s1: Vec3Float) -> Vec2Float: ...
    def get_d_matrix(self) -> Vec3Float: ...
    def get_D_matrix(self) -> Vec3Float: ...
    def get_directed_distance(self) -> float: ...
    def get_distance(self) -> float: ...
    def get_fast_axis(self) -> Vec3Float: ...
    @overload
    def get_lab_coord(self, xy: Vec2Float) -> Vec3Float: ...
    @overload
    def get_lab_coord(self, xy: flex.vec2_double) -> flex.vec3_double: ...
    def get_local_d_matrix(self) -> Vec9Float: ...
    def get_local_fast_axis(self) -> Vec3Float: ...
    def get_local_origin(self) -> Vec3Float: ...
    def get_local_slow_axis(self) -> Vec3Float: ...
    def get_normal_origin(self) -> Vec2Float: ...
    def get_normal(self) -> Vec3Float: ...
    def get_origin(self) -> Vec3Float: ...
    def get_parent_d_matrix(self) -> Vec9Float: ...
    def get_parent_fast_axis(self) -> Vec3Float: ...
    def get_parent_origin(self) -> Vec3Float: ...
    def get_parent_slow_axis(self) -> Vec3Float: ...
    def get_ray_intersection(self, s1: Vec3Float) -> Vec2Float: ...
    def get_slow_axis(self) -> Vec3Float: ...
    def set_frame(
        self, fast_axis: Vec3Float, slow_axis: Vec3Float, origin: Vec3Float
    ) -> None: ...
    def set_local_frame(
        self, fast_axis: Vec3Float, slow_axis: Vec3Float, origin: Vec3Float
    ) -> None: ...
    def set_parent_frame(
        self, fast_axis: Vec3Float, slow_axis: Vec3Float, origin: Vec3Float
    ) -> None: ...

class VirtualPanel(VirtualPanelFrame):
    def get_name(self) -> str: ...
    def get_type(self) -> str: ...
    def set_name(self, name: str) -> None: ...
    def set_type(self, typename: str) -> None: ...

class PanelData(VirtualPanel):
    def __init__(
        self,
        type: str,
        name: str,
        fast_axis: Vec3Float,
        slow_axis: Vec3Float,
        origin: Vec3Float,
        pixel_size: Vec2Float,
        image_size: Vec2Int,
        trusted_range: Vec2Float,
        thickness: float,
        material: str,
        mu: float = ...,
    ) -> None: ...
    def get_pixel_size(self) -> Vec2Float: ...
    def set_pixel_size(self, size: Vec2Float) -> None: ...
    def get_image_size(self) -> Vec2Int: ...
    def set_image_size(self, size: Vec2Int) -> None: ...
    def get_trusted_range(self) -> Vec2Int: ...
    def set_trusted_range(self, range: Vec2Float) -> None: ...
    def get_thickness(self) -> float: ...
    def set_thickness(self, thickness: float) -> None: ...
    def get_material(self) -> str: ...
    def set_material(self, material: str) -> None: ...
    def get_mu(self) -> float: ...
    def set_mu(self, mu: float) -> None: ...
    def get_raw_image_offset(self) -> Vec2Int: ...
    def set_raw_image_offset(self, offset: Vec2Int) -> None: ...
    def get_mask(self) -> flex_shared.tiny_int_4: ...
    def set_mask(self, mask: flex_shared.tiny_int_4): ...
    def add_mask(self, f0: int, s0: int, f1: int, s1: int) -> None: ...
    def is_similar_to(
        self,
        other: PanelData,
        fast_axis_tolerance: float = ...,
        slow_axis_tolerance: float = ...,
        origin_tolerance: float = ...,
        static_only: bool = ...,
        ignore_trusted_range: bool = ...,
    ) -> bool: ...

class Panel(PanelData):
    @overload
    def __init__(
        self,
        type: str,
        name: str,
        fast_axis: Vec3Float,
        slow_axis: Vec3Float,
        origin: Vec3Float,
        pixel_size: Vec2Float,
        image_size: Tuple[int, int],
        trusted_range: Vec2Float,
        thickness: float,
        material: str,
        mu: float = ...,
        identifier: str = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        type: str,
        name: str,
        fast_axis: Vec3Float,
        slow_axis: Vec3Float,
        origin: Vec3Float,
        pixel_size: Vec2Float,
        image_size: Tuple[int, int],
        trusted_range: Vec2Float,
        thickness: float,
        material: str,
        px_mm: PxMmStrategy,
        mu: float = ...,
        identifier: str = ...,
    ) -> None: ...
    @staticmethod
    @overload
    def from_dict(data: Dict) -> Panel: ...
    @staticmethod
    @overload
    def from_dict(data: Dict, dx: flex.double, dy: flex.double) -> Panel: ...
    def get_beam_centre_px(self, s0: Vec3Float) -> Vec2Float: ...
    def get_bidirectional_ray_intersection_px(self, s1: Vec3Float) -> Vec2Float: ...
    def get_cos2_two_theta_array(self, s0: Vec3Float) -> flex.double: ...
    def get_gain(self) -> float: ...
    def get_identifier(self) -> str: ...
    def get_image_size_mm(self) -> Vec2Float: ...
    def get_max_resolution_at_corners(self, s0: Vec3Float) -> float: ...
    def get_max_resolution_ellipse(self, s0: Vec3Float) -> float: ...
    def get_normal_origin_px(self) -> Vec2Float: ...
    def get_pedestal(self) -> float: ...
    def get_pixel_lab_coord(self, px: Vec2Float) -> Vec3Float: ...
    def get_projection_2d(self) -> Union[Tuple[Vec4Int, Vec2Int], Tuple[()]]: ...
    def get_px_mm_strategy(self) -> PxMmStrategy: ...
    def get_ray_intersection_px(self, s1: Vec3Float) -> Vec2Float: ...
    def get_resolution_at_pixel(self, s0: Vec3Float, xy: Vec2Float) -> float: ...
    def get_trusted_range_mask(
        self, image: Union[flex.int, flex.double]
    ) -> flex.bool: ...
    def get_two_theta_array(self, s0: Vec3Float) -> flex.double: ...
    def get_two_theta_at_pixel(self, s0: Vec3Float, xy: Vec2Float) -> float: ...
    def get_untrusted_rectangle_mask(self) -> flex.bool: ...
    def is_(self, other: Panel) -> bool: ...
    def is_coord_valid(self, xy: Vec2Float) -> bool: ...
    def is_coord_valid_mm(self, xy: Vec2Float) -> bool: ...
    def is_value_in_trusted_range(self, value: float) -> bool: ...
    @overload
    def millimeter_to_pixel(self, xy: Vec2Float) -> Vec2Float: ...
    @overload
    def millimeter_to_pixel(self, xy: flex.vec2_double) -> flex.vec2_double: ...
    @overload
    def pixel_to_millimeter(self, xy: Vec2Float) -> Vec2Float: ...
    @overload
    def pixel_to_millimeter(self, xy: flex.vec2_double) -> flex.vec2_double: ...
    def rotate_around_origin(
        self, axis: Vec3Float, angle: float, deg: bool = ...
    ) -> None: ...
    def set_gain(self, gain: float) -> None: ...
    def set_identifier(self, identifier: str) -> None: ...
    def set_pedestal(self, pedestal: float) -> None: ...
    def set_projection_2d(self, rotation: Vec4Int, translation: Vec2Int) -> None: ...
    def set_px_mm_strategy(self, strategy: PxMmStrategy) -> None: ...
    def to_dict(self) -> Dict: ...

class DetectorNode(Panel):
    @overload
    def add_group(self) -> DetectorNode: ...
    @overload
    def add_group(self, group: Panel) -> DetectorNode: ...
    @overload
    def add_panel(self) -> Panel: ...
    @overload
    def add_panel(self, panel: Panel) -> Panel: ...
    def children(self) -> Iterator[DetectorNode]: ...
    def empty(self) -> bool: ...
    def index(self) -> int: ...
    def is_group(self) -> bool: ...
    def is_panel(self) -> bool: ...
    def is_similar_to(
        self,
        other: PanelData,
        fast_axis_tolerance: float = ...,
        slow_axis_tolerance: float = ...,
        origin_tolerance: float = ...,
        static_only: bool = ...,
        ignore_trusted_range: bool = ...,
    ) -> bool: ...
    def parent(self) -> DetectorNode: ...
    def root(self) -> DetectorNode: ...
    def set_frame(
        self, fast_axis: Vec3Float, slow_axis: Vec3Float, origin: Vec3Float
    ) -> None: ...
    def set_local_frame(
        self, fast_axis: Vec3Float, slow_axis: Vec3Float, origin: Vec3Float
    ) -> None: ...
    def set_parent_frame(
        self, fast_axis: Vec3Float, slow_axis: Vec3Float, origin: Vec3Float
    ) -> None: ...
    def __getitem__(self, index: int) -> DetectorNode: ...
    def __iter__(boost) -> Iterator[DetectorNode]: ...
    def __len__(self) -> int: ...

class PxMmStrategy:
    def name(self) -> str: ...
    def to_millimeter(self, panel: Panel, xy: Vec2Float) -> Vec2Float: ...
    def to_pixel(self, panel: Panel, xy: Vec2Float) -> Vec2Float: ...

class SimplePxMmStrategy(PxMmStrategy):
    pass

class ParallaxCorrectedPxMmStrategy(PxMmStrategy):
    def __init__(self, mu: float, t0: float) -> None: ...
    def mu(self) -> float: ...
    def t0(self) -> float: ...

class OffsetPxMmStrategy(PxMmStrategy):
    def __init__(self, dx: flex.double, dy: flex.double) -> None: ...
    def dx(self) -> flex.double: ...
    def dy(self) -> flex.double: ...

class OffsetParallaxCorrectedPxMmStrategy(ParallaxCorrectedPxMmStrategy):
    def __init__(
        self, mu: float, t0: float, dx: flex.double, dy: flex.double
    ) -> None: ...
    def dx(self) -> flex.double: ...
    def dy(self) -> flex.double: ...

class ScanBase:
    def get_array_range(self) -> Vec2Int: ...
    def get_batch_for_array_index(self, index: int) -> int: ...
    def get_batch_for_image_index(self, index: int) -> int: ...
    def get_batch_offset(self) -> int: ...
    def set_batch_offset(self, batch_offset: int) -> int: ...
    def get_batch_range(self) -> Vec2Int: ...
    def get_image_range(self) -> Vec2Int: ...
    def set_image_range(self, image_range: Vec2Int) -> None: ...
    def get_num_images(self) -> int: ...
    def get_valid_image_ranges(self, i: str) -> List[Vec2Int]: ...
    def set_valid_image_ranges(self, i: str, ranges: List[Vec2Int]) -> None: ...
    def is_array_index_valid(self, index: int) -> bool: ...
    def is_batch_valid(self, batch: int) -> bool: ...
    def is_image_index_valid(self, index: int) -> bool: ...
    def is_still(self) -> bool: ...

class Scan(ScanBase):
    @overload
    def __init__(
        self,
        image_range: Vec2Int,
        oscillation: Vec2Float,
        batch_offset: int,
        deg: bool = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        image_range: Vec2Int,
        oscillation: Vec2Float,
        exposure_times: flex.double,
        epochs: flex.double,
        batch_offset: int,
        deg: bool = ...,
    ) -> None: ...
    def append(self, other: Scan, scan_tolerance: float) -> None: ...
    @staticmethod
    def from_dict(data: Dict) -> Scan: ...
    @overload
    def get_angle_from_array_index(self, index: float, deg: bool = ...) -> float: ...
    @overload
    def get_angle_from_array_index(
        self, index: flex.double, deg: bool = ...
    ) -> flex.double: ...
    def get_angle_from_image_index(self, index: float, deg: bool = ...) -> float: ...
    @overload
    def get_array_index_from_angle(self, angle: float, deg: bool = ...) -> float: ...
    @overload
    def get_array_index_from_angle(
        self, angle: flex.double, deg: bool = ...
    ) -> flex.double: ...
    def get_array_indices_with_angle(
        self, angle: float, deg: bool = ...
    ) -> flex.vec2_double: ...
    def get_epochs(self) -> flex.double: ...
    def set_epochs(self, epochs: flex.double) -> None: ...
    def get_exposure_times(self) -> flex.double: ...
    def set_exposure_times(self, times: flex.double) -> None: ...
    def get_image_epoch(self, index: int) -> float: ...
    def get_image_index_from_angle(self, angle: float, deg: bool = ...) -> float: ...
    def get_image_indices_with_angle(
        self, angle: float, deg: bool = ...
    ) -> flex.vec2_double: ...
    def get_image_oscillation(self, index: int, deg: bool = ...) -> Vec2Float: ...
    def set_image_range(self, image_range: Vec2Int) -> None: ...
    def get_oscillation(self, deg: bool = ...) -> Vec2Float: ...
    def set_oscillation(self, oscillation: Vec2Float, deg: bool = ...) -> None: ...
    def get_oscillation_range(self, deg: bool = ...) -> Vec2Float: ...
    @overload
    def is_angle_valid(self, angle: float, deg: bool = ...) -> bool: ...
    @overload
    def is_angle_valid(self, angle: flex.double, deg: bool = ...) -> flex.bool: ...
    def swap(self, other: Scan) -> None: ...
    def to_dict(self) -> Dict: ...
    def __getitem__(self, index: Union[int, slice]) -> Scan: ...
    def __len__(self) -> int: ...
    def __add__(self, other: Scan) -> Scan: ...
    def __iadd__(self, other: Scan) -> Scan: ...
    def __ge__(self, other: Scan) -> bool: ...
    def __gt__(self, other: Scan) -> bool: ...
    def __le__(self, other: Scan) -> bool: ...
    def __lt__(self, other: Scan) -> bool: ...

class Spectrum:
    @overload
    def __init__(self, other: Spectrum) -> None: ...
    @overload
    def __init__(self, energies: flex.double, weights: flex.double) -> None: ...
    @staticmethod
    def from_dict(data: Dict) -> Spectrum: ...
    def get_emax_eV(self) -> float: ...
    def get_emin_eV(self) -> float: ...
    def get_energies_eV(self) -> flex.double: ...
    def get_weighted_energy_eV(self) -> float: ...
    def get_weighted_energy_variance(self) -> float: ...
    def get_weighted_wavelength(self) -> float: ...
    def get_weights(self) -> flex.double: ...
    def to_dict(self) -> Dict: ...

class flex_Beam(FlexPlain[Beam]):
    pass

class flex_Spectrum(FlexPlain[Spectrum]):
    pass

def get_mod2pi_angles_in_range(
    range: Vec2Float, angle: float, deg: bool = False
) -> flex.double: ...
def is_angle_in_range(range: Vec2Float, angle: float, deg: bool = False) -> bool: ...
def get_range_of_mod2pi_angles(
    range: Vec2Float, angle: float, deg: bool = False
) -> Vec2Float: ...
@overload
def parallax_correction(
    d: float, la: float, xy0: Vec2Float, xy: Vec2Float
) -> Vec2Float: ...
@overload
def parallax_correction(
    mu: float,
    t0: float,
    xy: Vec2Float,
    fast_axis: Vec3Float,
    slow_axis: Vec3Float,
    origin: Vec3Float,
) -> Vec2Float: ...
@overload
def parallax_correction_inv(
    d: float, la: float, xy0: Vec2Float, xy: Vec2Float
) -> Vec2Float: ...
@overload
def parallax_correction_inv(
    mu: float,
    t0: float,
    xy: Vec2Float,
    fast_axis: Vec3Float,
    slow_axis: Vec3Float,
    origin: Vec3Float,
) -> Vec2Float: ...
