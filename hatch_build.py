"""
Dynamically generate the list of console_scripts dxtbx.format entry-points.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

from hatchling.metadata.plugin.interface import MetadataHookInterface


def get_entry_point(
    filename: Path, prefix: str, import_path: str
) -> list[tuple[str, str]]:
    """Returns any entry point strings for a given path.

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
        x
        for x in tree.body
        if (isinstance(x, ast.FunctionDef) and x.name == "run")
        or (isinstance(x, ast.ImportFrom) and "run" in [a.name for a in x.names])
    )
    if not has_run:
        return []
    # Find if we need an alternate name via LIBTBX_SET_DISPATCHER_NAME
    alternate_names = re.findall(
        r"^#\s*LIBTBX_SET_DISPATCHER_NAME\s+(.*)$", contents, re.M
    )
    if alternate_names:
        return [
            (name, f"{import_path}.{filename.stem}:run") for name in alternate_names
        ]

    return [(f"{prefix}.{filename.stem}", f"{import_path}.{filename.stem}:run")]


def enumerate_format_classes(path: Path) -> list[tuple(str, str)]:
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
                    (
                        f"{classname}:{','.join(base_names)}",
                        f"dxtbx.format.{filename.stem}:{classname}",
                    )
                )
    return format_classes


class CustomMetadataHook(MetadataHookInterface):
    def update(self, metadata):
        scripts = metadata.setdefault("scripts", {})
        package_path = Path(self.root) / "src" / "dxtbx"
        for file in package_path.joinpath("command_line").glob("*.py"):
            for name, symbol in get_entry_point(file, "dxtbx", "dxtbx.command_line"):
                if name not in scripts:
                    scripts[name] = symbol

        plugins = metadata.setdefault("entry-points", {})
        formats = plugins.setdefault("dxtbx.format", {})
        for name, symbol in sorted(enumerate_format_classes(package_path / "format")):
            if name not in formats:
                formats[name] = symbol
