from __future__ import annotations

from pathlib import Path

import setuptools

from build import build

# Static version number which is updated by bump2version
# Do not change this manually - use 'bump2version <major/minor/patch/release>'
__version_tag__ = "3.19.dev"

setup_kwargs = {
    "name": "dxtbx",
    "version": __version_tag__,
    "long_description": Path(__file__).parent.joinpath("README.md").read_text(),
    "description": "Diffraction Experiment Toolbox",
    "author": "Diamond Light Source",
    "license": "BSD-3-Clause",
    "author_email": "dials-support@lists.sourceforge.net",
    "project_urls": {
        "homepage": "https://dials.github.io",
        "repository": "https://github.com/cctbx/dxtbx",
    },
    "packages": setuptools.find_packages(where="src"),
    "package_dir": {"": "src"},
    "package_data": {
        "": ["*"],
        "dxtbx": ["array_family/*", "boost_python/*", "example/*", "py.typed"],
        "dxtbx.format": ["boost_python/*"],
        "dxtbx.masking": ["boost_python/*"],
        "dxtbx.model": ["boost_python/*"],
    },
    "entry_points": {
        "libtbx.precommit": ["dxtbx=dxtbx"],
        "libtbx.dispatcher.script": ["pytest=pytest"],
    },
}

build(setup_kwargs)
setuptools.setup(**setup_kwargs)
