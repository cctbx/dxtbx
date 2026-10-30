#!/usr/bin/env python
from __future__ import annotations

import shlex
import subprocess
import sys
import textwrap
from pathlib import Path

template = Path(__file__).parent.joinpath("meta.yaml.template").read_text()


# Get the indents
def _get_indent_on_line_containing(text: str):
    lines = template.splitlines()
    for num, line in enumerate(lines):
        if text in lines[num]:
            break
    else:
        raise ValueError(f"Could not find {text} in template")

    return len(line) - len(line.lstrip())


# Work out how much we need to indent the requirement sections
DEPS_INDENT = _get_indent_on_line_containing("$REQUIRES")
TEST_INDENT = _get_indent_on_line_containing("$TEST_REQUIRES")


def _run_parser(*args):
    # Generate the dependency sections
    cmd = [
        sys.executable,
        Path(__file__).parent / "parse_dependencies.py",
        Path(__file__).parent.parent / "dependencies.yaml",
        "--conda-build",
        *args,
    ]
    cmd = [str(x) for x in cmd]
    try:
        proc = subprocess.run(
            cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
    except subprocess.CalledProcessError as p:
        print("Error getting dependencies " + " ".join(args) + ":")
        print(" + " + shlex.join(cmd))
        print(p.stdout)
        sys.exit(1)

    return proc.stdout


deps = _run_parser("--build", "--host", "--run")
deps_test = textwrap.dedent("\n".join(_run_parser("--test").splitlines()[1:]))

template = template.replace(
    "$REQUIRES", textwrap.indent(deps, " " * DEPS_INDENT).lstrip()
)
template = template.replace(
    "$TEST_REQUIRES", textwrap.indent(deps_test, " " * TEST_INDENT).lstrip()
)
out_file = Path(__file__).parent.parent / "meta.yaml"
out_file.write_text(template)
print(f"Written to {out_file}")
