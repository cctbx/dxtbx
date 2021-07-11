"""
Handle dynamic aspects of setup.py and building.

This is separate because the non-dynamic stuff can generally be moved
out of a setup.py, but mainly because at the moment it's how poetry
offloads the unresolved build phases.
"""

import ast
import itertools
import re
import sys
from pathlib import Path
from typing import Any, Dict, List


def get_entry_point(filename: Path, prefix: str, import_path: str) -> List[str]:
    """Returns the entry point string for a given path.

    This looks for LIBTBX_SET_DISPATCHER_NAME, and a root function
    named 'run'. It can return multiple results for each file, if more
    than one dispatcher name is bound.

    Args:
        filename:
            The python file to parse. Will look for a run() function
            and any number of LIBTBX_SET_DISPATCHER_NAME.
        prefix: The prefix to output the entry point console script with
        import_path: The import path to get to the package the file is in

    Returns:
        A list of entry_point specifications
    """
    contents = filename.read_text()
    tree = ast.parse(contents)
    # Find root functions named "run"
    has_run = any(
        x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "run"
    )
    if not has_run:
        return []
    # Find if we need an alternate name via LIBTBX_SET_DISPATCHER_NAME
    alternate_names = re.findall(
        r"^#\s*LIBTBX_SET_DISPATCHER_NAME\s+(.*)$", contents, re.M
    )
    if alternate_names:
        return [f"{name}={import_path}.{filename.stem}:run" for name in alternate_names]

    return [f"{prefix}.{filename.stem}={import_path}.{filename.stem}:run"]


def enumerate_format_classes(path: Path) -> List[str]:
    """Find all Format*.py files and contained Format classes in a path"""
    format_classes = []
    for filename in path.glob("Format*.py"):
        content = filename.read_bytes()
        try:
            parsetree = ast.parse(content)
        except SyntaxError:
            print(f"  *** Could not parse {filename.name}")
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
                    f"{classname}:{','.join(base_names)} = dxtbx.format.{filename.stem}:{classname}"
                )
                # print("  found", classname, " based on ", str(base_names))
    return format_classes


def build(setup_kwargs: Dict[str, Any]) -> None:
    """Called by setup.py to inject any dynamic configuration"""
    package_path = Path(__file__).parent / "src" / "dxtbx"
    entry_points = setup_kwargs.setdefault("entry_points", {})
    console_scripts = entry_points.setdefault("console_scripts", [])
    # Work out what dispatchers to add
    all_dispatchers = sorted(
        itertools.chain(
            *[
                get_entry_point(f, "dxtbx", "dxtbx.command_line")
                for f in (package_path / "command_line").glob("*.py")
            ]
        )
    )
    console_scripts.extend(x for x in all_dispatchers if x not in console_scripts)
    libtbx_dispatchers = entry_points.setdefault("libtbx.dispatcher.script", [])
    libtbx_dispatchers.extend(
        "{name}={name}".format(name=x.split("=")[0]) for x in console_scripts
    )

    dxtbx_format = entry_points.setdefault("dxtbx.format", [])
    format_classes = sorted(enumerate_format_classes(package_path / "format"))
    dxtbx_format.extend([x for x in format_classes if x not in dxtbx_format])

    print(f"Found {len(entry_points['console_scripts'])} dxtbx dispatchers")
    print(f"Found {len(entry_points['dxtbx.format'])} Format classes")


if __name__ == "__main__":
    sys.exit("Cannot call build.py directly, please use setup.py instead")
