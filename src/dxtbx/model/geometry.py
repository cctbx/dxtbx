"""
Support for experimental metadata that a Format class cannot read from the
image files, and which must therefore be supplied by the user.

A Format class declares the geometry it cannot determine by listing PHIL paths
(relative to the ``geometry`` scope) in :attr:`dxtbx.format.Format.missing_metadata`,
or by calling :meth:`dxtbx.format.Format.declare_missing_metadata` from
``_start()``. The user's overrides are then applied to those items alone, as the
Format instance is constructed, and a :exc:`MissingMetadataError` is raised for
anything left unsatisfied. This replaces the previous situation in which a
Format had to invent dummy values that the user might silently fail to override.
"""

from __future__ import annotations

import copy

import libtbx.phil

from dxtbx.model.beam import BeamFactory, beam_phil_scope
from dxtbx.model.detector import DetectorFactory, detector_phil_scope
from dxtbx.model.goniometer import GoniometerFactory, goniometer_phil_scope
from dxtbx.model.scan import ScanFactory, scan_phil_scope


def _build_geometry_phil_scope():
    # Composed by adopting the four model scopes rather than with
    # "include scope dxtbx.model.beam.beam_phil_scope", because this module is
    # imported from dxtbx/model/__init__.py, and phil resolves an include by
    # attribute lookup on the dxtbx package - which cannot yet see dxtbx.model
    # while its __init__ is still executing.
    scope = libtbx.phil.parse(
        """
geometry
  .help = "Allow overrides of experimental geometry"
  .expert_level = 2
{
}
"""
    )
    for model_scope in (
        beam_phil_scope,
        detector_phil_scope,
        goniometer_phil_scope,
        scan_phil_scope,
    ):
        scope.objects[0].adopt_scope(model_scope)
    return scope


geometry_phil_scope = _build_geometry_phil_scope()

_PREFIX = "geometry."

# Parameters that modify an existing model rather than describe metadata. These
# must not be declarable: they are the only overrides that are not idempotent,
# and dials.import's ManualGeometryUpdater applies the user's geometry scope a
# second time after the Format has been constructed. scan.image_range is a
# property of a whole sequence, which a per-image Format instance cannot express.
_NON_METADATA_PATHS = frozenset(
    {
        "goniometer.invert_rotation_axis",
        "scan.extrapolate_scan",
        "scan.image_range",
    }
)

# Declaring either of these is satisfied by the user supplying either of them.
# GoniometerFactory.from_phil already rejects both being set at once.
_ALIASES = {
    "goniometer.axis": ("goniometer.axis", "goniometer.axes"),
    "goniometer.axes": ("goniometer.axis", "goniometer.axes"),
}


class MissingMetadataError(RuntimeError):
    """
    Raised when a Format class declares metadata it cannot read from the image
    file and the user has not supplied an override for it.
    """

    def __init__(self, format_name, unsatisfied, image_file="", hint=""):
        self.format_name = format_name
        self.unsatisfied = tuple(sorted(unsatisfied))
        self.image_file = image_file
        self.hint = hint
        super().__init__(self._message())

    def _message(self):
        where = f" from {self.image_file}" if self.image_file else ""
        listing = "\n".join(f"    geometry.{path}" for path in self.unsatisfied)
        example = " ".join(f"geometry.{path}=<value>" for path in self.unsatisfied)
        message = (
            f"The format class {self.format_name} cannot determine the following "
            f"experimental metadata{where}:\n\n"
            f"{listing}\n\n"
            "These values are not in the image headers, so they must be supplied "
            "explicitly, either on the command line\n\n"
            f"    dials.import <images> {example}\n\n"
            "or in a site.phil file."
        )
        if self.hint:
            message += f"\n\n{self.hint}"
        return message


def declarable_metadata_paths():
    """
    The PHIL paths, relative to the ``geometry`` scope, that a Format class may
    declare as missing.

    Excludes the modifiers in :data:`_NON_METADATA_PATHS`, and the repeatable
    ``detector.panel`` and ``detector.hierarchy`` sub-scopes, which have no
    unambiguous dotted path.
    """
    try:
        return declarable_metadata_paths._cache
    except AttributeError:
        pass
    paths = set()
    for definition in geometry_phil_scope.all_definitions():
        path = definition.path[len(_PREFIX) :]
        if ".panel" in path or ".hierarchy" in path:
            continue
        if path in _NON_METADATA_PATHS:
            continue
        paths.add(path)
    declarable_metadata_paths._cache = frozenset(paths)
    return declarable_metadata_paths._cache


def validate_metadata_path(path):
    """Raise ValueError unless ``path`` may be declared as missing metadata."""
    if path in declarable_metadata_paths():
        return
    if path in _NON_METADATA_PATHS:
        raise ValueError(
            f"{path!r} modifies an existing model rather than describing "
            "metadata, so it cannot be declared as missing."
        )
    if ".panel" in path or ".hierarchy" in path:
        raise ValueError(
            f"{path!r} is part of a repeatable sub-scope, which cannot yet be "
            "declared as missing metadata. Build a placeholder detector and "
            "declare detector.distance or a beam centre instead."
        )
    raise ValueError(
        f"{path!r} is not a geometry metadata path. Valid paths are: "
        + ", ".join(sorted(declarable_metadata_paths()))
    )


def _graft(root, path, definition):
    """Add a PHIL definition to ``root`` under scopes named by its dotted path,
    reusing any scope already there so the result is a single merged tree."""
    parts = path.split(".")
    parent = root
    for name in parts[:-1]:
        for obj in parent.objects:
            if obj.is_scope and obj.name == name:
                parent = obj
                break
        else:
            child = libtbx.phil.scope(name=name)
            parent.adopt(child)
            parent = child
    parent.adopt(definition.customized_copy(name=parts[-1]))


def select_geometry_phil(source, paths):
    """
    Mask ``source`` down to the requested paths.

    Args:
        source: a PHIL scope whose definitions live under ``geometry.``.
        paths: an iterable of paths relative to the ``geometry`` scope.

    Returns:
        (masked_source, supplied) where masked_source is a nameless scope
        suitable as a ``fetch()`` source for :data:`geometry_phil_scope`, and
        supplied is the set of requested paths actually present in source.
    """
    wanted = set()
    for path in paths:
        wanted.update(_ALIASES.get(path, (path,)))
    root = libtbx.phil.scope(name="")
    supplied = set()
    for definition in source.all_definitions():
        if not definition.path.startswith(_PREFIX):
            continue
        path = definition.path[len(_PREFIX) :]
        if path not in wanted:
            continue
        # An explicit None is how phil spells "unset", so it is not an override.
        # Skipping it also keeps the alias expansion above from reporting
        # goniometer.axes as supplied when only goniometer.axis was given.
        if definition.object.extract() is None:
            continue
        _graft(root, definition.path, definition.object)
        supplied.add(path)
    return root, supplied


def geometry_params_for_paths(overrides, paths):
    """
    Parse ``overrides`` and mask it down to ``paths``.

    Args:
        overrides: a PHIL string (or scope) rooted at ``geometry``, or None.
        paths: an iterable of paths relative to the ``geometry`` scope.

    Returns:
        (params, supplied) where params is a fresh extract of the ``geometry``
        scope carrying only the requested overrides.
    """
    if overrides is None or (isinstance(overrides, str) and not overrides.strip()):
        source = libtbx.phil.scope(name="", objects=[])
        supplied = set()
    else:
        if isinstance(overrides, str):
            overrides = libtbx.phil.parse(overrides)
        source, supplied = select_geometry_phil(overrides, paths)
    # Note fetch(), never geometry_phil_scope.extract(): the latter materialises
    # a bogus detector.panel and hierarchy.group entry with id=None, which makes
    # DetectorFactory.from_phil fail. A fetch of an empty source leaves both
    # repeatable scopes as empty lists, which is what the factories expect.
    return geometry_phil_scope.fetch(source=source).extract().geometry, supplied


def _satisfied(path, supplied):
    return bool(supplied.intersection(_ALIASES.get(path, (path,))))


def _apply_scan_overrides(params, scan):
    """
    Apply scan overrides, accounting for one-image scans.

    Format.get_imageset builds the scan of a multi-file sequence by summing the
    single-image scans of each Format instance. Giving every one of them the
    same oscillation start would make that sum fail, so offset the start by the
    image index. Multi-image formats (a single stack file) are unaffected and
    behave exactly as ManualGeometryUpdater does.
    """
    oscillation = params.scan.oscillation
    if oscillation is not None and scan is not None and scan.get_num_images() == 1:
        first = scan.get_image_range()[0]
        params = copy.deepcopy(params)
        params.scan.oscillation = (
            oscillation[0] + (first - 1) * oscillation[1],
            oscillation[1],
        )
    return ScanFactory.from_phil(params, scan)


def apply_geometry_overrides(
    declared,
    overrides,
    beam=None,
    detector=None,
    goniometer=None,
    scan=None,
    format_name="",
    image_file="",
):
    """
    Fill in declared-missing metadata from the user's overrides.

    Only models for which a path was both declared and supplied are touched, so
    a Format declaring nothing gets back exactly the models it passed in.

    Returns:
        (beam, detector, goniometer, scan)

    Raises:
        MissingMetadataError: if any declared path has no user-supplied value.
    """
    declared = set(declared)
    if not declared:
        return beam, detector, goniometer, scan

    params, supplied = geometry_params_for_paths(overrides, declared)

    unsatisfied = {path for path in declared if not _satisfied(path, supplied)}
    if unsatisfied:
        raise MissingMetadataError(format_name, unsatisfied, image_file)

    by_model = {"beam": set(), "detector": set(), "goniometer": set(), "scan": set()}
    for path in supplied:
        by_model[path.split(".", 1)[0]].add(path)

    # A model that the Format left as None can only be built from scratch if the
    # overrides carry everything the factory needs. Check up front so the user
    # gets an actionable message rather than a bare RuntimeError from the factory.
    if beam is None and by_model["beam"]:
        needed = {"beam.wavelength", "beam.direction"}
        if not needed <= supplied:
            raise MissingMetadataError(
                format_name,
                needed - supplied,
                image_file,
                hint=(
                    f"{format_name} provides no beam model at all, so building "
                    "one requires both a wavelength and a direction."
                ),
            )
    if detector is None and by_model["detector"]:
        raise MissingMetadataError(
            format_name,
            by_model["detector"],
            image_file,
            hint=(
                f"{format_name} provides no detector model at all. A detector "
                "cannot be built from these overrides alone; the format class "
                "should construct a placeholder detector and declare only the "
                "items it cannot determine."
            ),
        )

    # Beam first: the detector beam-centre and parallax branches need it. This
    # is the same order that dials.import's ManualGeometryUpdater uses.
    if by_model["beam"]:
        beam = BeamFactory.from_phil(params, beam)
    if by_model["detector"]:
        detector = DetectorFactory.from_phil(params, detector, beam)
    if by_model["goniometer"]:
        goniometer = GoniometerFactory.from_phil(params, goniometer)
    if by_model["scan"]:
        scan = _apply_scan_overrides(params, scan)

    return beam, detector, goniometer, scan


def check_getters_not_bypassed(format_class, missing):
    """
    Raise TypeError if ``format_class`` declares missing metadata for a model
    whose getter bypasses the corresponding ``_<model>_instance`` attribute.

    Overrides are applied by writing those attributes, so a format that
    overrides e.g. ``get_beam`` to rebuild the model on every call would discard
    them silently. Fail loudly at the point the author opts in instead.
    """
    from dxtbx.format.Format import Format
    from dxtbx.format.FormatMultiImage import FormatMultiImage

    bypassed = []
    for model in sorted({path.split(".", 1)[0] for path in missing}):
        name = f"get_{model}"
        safe = (getattr(Format, name), getattr(FormatMultiImage, name))
        if getattr(format_class, name) not in safe:
            bypassed.append(name)
    if bypassed:
        raise TypeError(
            f"{format_class.__name__} declares missing metadata but overrides "
            f"{', '.join(bypassed)} so that the _<model>_instance attributes are "
            "bypassed, which would silently discard the user's overrides. Route "
            "these getters through self._<model>_instance, as Format and "
            "FormatMultiImage do, before declaring missing metadata."
        )
