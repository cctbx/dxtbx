# coding: utf-8
from __future__ import absolute_import, division, print_function

import ast
import imp
import pkgutil
import sys

import past.builtins  # noqa; lgtm; Ensure package grammar files are generated

import libtbx.pkg_utils

import dxtbx.format

try:
    import dials.precommitbx.nagger

    dials.precommitbx.nagger.nag()
except ImportError:
    pass

if sys.version_info.major == 2:
    sys.exit("Python 2 is no longer supported")


# --- format class registration ---

print("Enumerating format classes:")
format_classes = []
for _, name, _ in pkgutil.iter_modules(dxtbx.format.__path__):
    if not name.startswith("Format"):
        continue
    try:
        fid, pathname, desc = imp.find_module(name, dxtbx.format.__path__)
    except Exception:
        fid = None
    if not fid:
        print("  *** Could not read %s" % name)
        continue
    if desc[0] == ".pyc":
        print("  *** %s only present in compiled form, ignoring" % name)
        continue
    content = fid.read()
    fid.close()
    try:
        parsetree = ast.parse(content)
    except SyntaxError:
        print("  *** Could not parse %s" % name)
        continue
    for top_level_def in parsetree.body:
        if not isinstance(top_level_def, ast.ClassDef):
            continue
        base_names = [
            baseclass.id
            for baseclass in top_level_def.bases
            if isinstance(baseclass, ast.Name) and baseclass.id.startswith("Format")
        ]
        if base_names:
            classname = top_level_def.name
            format_classes.append(
                "{classname}:{baseclasses} = dxtbx.format.{modulename}:{classname}".format(
                    classname=classname,
                    modulename=name,
                    baseclasses=",".join(base_names),
                )
            )
            print("  found", classname, " based on ", str(base_names))

# Sanity check to catch old configurations
assert not any(
    x.endswith("cctbx_project/dxtbx") for x in dxtbx.__path__
), "dxtbx found inside cctbx_project. Please remove modules/cctbx_project/dxtbx contents"
assert format_classes, "No format classes found; something went wrong ¯\\_(ツ)_/¯"

libtbx.pkg_utils.define_entry_points({"dxtbx.format": sorted(format_classes)})
