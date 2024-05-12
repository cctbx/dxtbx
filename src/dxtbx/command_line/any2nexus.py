from __future__ import annotations

import sys

from dxtbx.format.nxmx_writer import run

if __name__ == "__main__":
    run(sys.argv[1:])
