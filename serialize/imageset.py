import collections
import os

import six.moves.cPickle as pickle

from dxtbx.format.image import ImageBool, ImageDouble
from dxtbx.imageset import ImageSequence, ImageSet, ImageSetFactory
from dxtbx.model import BeamFactory, DetectorFactory, GoniometerFactory, ScanFactory
from dxtbx.serialize.filename import resolve_path


def filename_to_absolute(filename):
    """Convert filenames to absolute form."""

    if isinstance(filename, list):
        return [os.path.abspath(f) for f in filename]

    return os.path.abspath(filename)


def filename_or_none(filename):
    if filename is None or filename == "":
        return None
    return filename_to_absolute(filename)


def basic_imageset_to_dict(imageset):
    """Convert an imageset to a dictionary

    Params:
        imageset The imageset

    Returns:
        A dictionary of the parameters

    """

    return collections.OrderedDict(
        [
            ("__id__", "imageset"),
            ("filenames", filename_to_absolute(imageset.paths())),
            ("mask", filename_or_none(imageset.external_lookup.mask.filename)),
            ("gain", filename_or_none(imageset.external_lookup.gain.filename)),
            ("pedestal", filename_or_none(imageset.external_lookup.pedestal.filename)),
            ("beam", imageset.get_beam(0).to_dict()),
            ("detector", imageset.get_detector(0).to_dict()),
        ]
    )


def imagesequence_to_dict(sequence):
    """Convert a sequence to a dictionary

    Params:
        sequence The sequence

    Returns:
        A dictionary of the parameters

    """

    return collections.OrderedDict(
        [
            ("__id__", "imageset"),
            ("template", filename_to_absolute(sequence.get_template())),
            ("mask", filename_or_none(sequence.external_lookup.mask.filename)),
            ("gain", filename_or_none(sequence.external_lookup.gain.filename)),
            ("pedestal", filename_or_none(sequence.external_lookup.pedestal.filename)),
            ("beam", sequence.get_beam().to_dict()),
            ("detector", sequence.get_detector().to_dict()),
            ("goniometer", sequence.get_goniometer().to_dict()),
            ("scan", sequence.get_scan().to_dict()),
        ]
    )


def imageset_to_dict(imageset):
    """Convert the imageset to a dictionary

    Params:
        imageset The imageset

    Returns:
        A dictionary of the parameters

    """
    # If this is an imageset then return a list of filenames
    if isinstance(imageset, ImageSequence):
        return imagesequence_to_dict(imageset)
    elif isinstance(imageset, ImageSet):
        return basic_imageset_to_dict(imageset)
    else:
        raise TypeError("Unknown ImageSet Type")


def basic_imageset_from_dict(d, directory=None):
    """Construct an ImageSet class from the dictionary."""
    # Get the filename list and create the imageset
    filenames = [resolve_path(str(p), directory=directory) for p in d["filenames"]]
    imageset = ImageSetFactory.new(filenames)[0]

    # Set some external lookups
    if "mask" in d and d["mask"] is not None and d["mask"] != "":
        path = resolve_path(d["mask"], directory=directory)
        with open(path) as infile:
            imageset.external_lookup.mask.filename = path
            imageset.external_lookup.mask.data = ImageBool(pickle.load(infile))
    if "gain" in d and d["gain"] is not None and d["gain"] != "":
        path = resolve_path(d["gain"], directory=directory)
        with open(path) as infile:
            imageset.external_lookup.gain.filename = path
            imageset.external_lookup.gain.data = ImageDouble(pickle.load(infile))
    if "pedestal" in d and d["pedestal"] is not None and d["pedestal"] != "":
        path = resolve_path(d["pedestal"], directory=directory)
        with open(path) as infile:
            imageset.external_lookup.pedestal.filename = path
            imageset.external_lookup.pedestal.data = ImageDouble(pickle.load(infile))

    # Get the existing models as dictionaries
    beam_dict = imageset.get_beam(0).to_dict()
    detector_dict = imageset.get_detector(0).to_dict()

    # Set models
    imageset.set_beam(BeamFactory.from_dict(d.get("beam"), beam_dict))
    imageset.set_detector(DetectorFactory.from_dict(d.get("detector"), detector_dict))

    return imageset


def imagesequence_from_dict(d, check_format=True, directory=None):
    """Construct and image sequence from the dictionary."""
    # Get the template (required)
    template = resolve_path(str(d["template"]), directory=directory)

    # If the scan isn't set, find all available files
    scan_dict = d.get("scan")
    if scan_dict is None:
        image_range = None
    else:
        image_range = scan_dict.get("image_range")

    # Set the models with the exisiting models as templates
    beam = BeamFactory.from_dict(d.get("beam"))
    goniometer = GoniometerFactory.from_dict(d.get("goniometer"))
    detector = DetectorFactory.from_dict(d.get("detector"))
    scan = ScanFactory.from_dict(d.get("scan"))

    # Construct the sequence
    try:
        sequence = ImageSetFactory.from_template(
            template,
            image_range,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            check_format=check_format,
        )[0]
    except Exception:
        indices = list(range(image_range[0], image_range[1] + 1))
        sequence = ImageSetFactory.make_sequence(
            template,
            indices,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            check_format=check_format,
        )

    # Set some external lookups
    if "mask" in d and d["mask"] is not None and d["mask"] != "":
        path = resolve_path(d["mask"], directory=directory)
        with open(path) as infile:
            sequence.external_lookup.mask.filename = path
            sequence.external_lookup.mask.data = ImageBool(pickle.load(infile))
    if "gain" in d and d["gain"] is not None and d["gain"] != "":
        path = resolve_path(d["gain"], directory=directory)
        with open(path) as infile:
            sequence.external_lookup.gain.filename = path
            sequence.external_lookup.gain.data = ImageDouble(pickle.load(infile))
    if "pedestal" in d and d["pedestal"] is not None and d["pedestal"] != "":
        path = resolve_path(d["pedestal"], directory=directory)
        with open(path) as infile:
            sequence.external_lookup.pedestal.filename = path
            sequence.external_lookup.pedestal.data = ImageDouble(pickle.load(infile))

    return sequence


def imageset_from_dict(d, check_format=True, directory=None):
    """Convert the dictionary to a sequence

    Params:
        d The dictionary of parameters

    Returns:
        The sequence

    """
    # Check the input
    if d is None:
        return None

    # Check the version and id
    if str(d["__id__"]) != "imageset":
        raise ValueError('"__id__" does not equal "imageset"')

    if "filenames" in d:
        return basic_imageset_from_dict(d, directory=directory)
    elif "template" in d:
        return imagesequence_from_dict(
            d, check_format=check_format, directory=directory
        )
    else:
        raise TypeError("Unable to deserialize given imageset")
