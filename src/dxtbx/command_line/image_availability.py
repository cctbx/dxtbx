# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.image_availability
"""
Live-monitor how much of an NXmx master file has actually been written to disk.

During serial (SSX) data collection at DLS the master ``.h5``/``.nxs`` file and its
underlying ``data_*.h5`` sources are created up front (the virtual dataset declares the
full planned number of frames) and then filled incrementally, typically via SWMR. File
presence therefore does not mean the image data is present.

This utility does *no* processing: it simply polls the availability tools in
``dxtbx.nexus`` and shows, live, how many frames of each given master are genuinely
written so far. Handy for watching a collection fill up, or for debugging the
``xia2.ssx wait_for_images`` live-processing mode.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

from tqdm import tqdm

import dxtbx.util
from dxtbx.nexus import get_frame_counts

UNIT = " images"


def _snapshot(label: str, avail: int, total: int) -> str:
    """Render a single, static tqdm-styled line (no rate or ETA available)."""
    if total <= 0:
        return f"{label}  unreadable / no frames declared"
    return tqdm.format_meter(n=avail, total=total, elapsed=0, prefix=label, unit=UNIT)


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser(
        description=(
            "Live-monitor how many image frames of NXmx master file(s) have been "
            "written to disk. Does no processing — only inspects data availability."
        )
    )
    parser.add_argument(
        "master", nargs="+", metavar="MASTER.h5", help="NXmx master file(s) to watch"
    )
    parser.add_argument(
        "--method",
        choices=["swmr", "file_close"],
        default="swmr",
        help=(
            "How to test source-file readiness: 'swmr' reads the current written extent "
            "(per-image); 'file_close' waits until the writer closes each source file "
            "(per-block). Default: swmr."
        ),
    )
    parser.add_argument(
        "-n",
        "--interval",
        type=float,
        default=2.0,
        metavar="SECONDS",
        help="Polling interval in seconds (default: 2).",
    )
    parser.add_argument(
        "-1",
        "--once",
        action="store_true",
        help="Print a single snapshot and exit (no live updating).",
    )
    parser.add_argument(
        "--watch-complete",
        action="store_true",
        help=(
            "Keep polling even after every file is complete (default: exit once all "
            "files reach their full frame count)."
        ),
    )
    options = parser.parse_args(args)

    masters = options.master
    labels = [os.path.basename(m) for m in masters]
    width = max(len(label) for label in labels)
    labels = [label.ljust(width) for label in labels]

    # tqdm can only redraw in place on an interactive terminal; anywhere else fall
    # back to plain snapshot lines, one block per poll.
    live = not options.once and sys.stdout.isatty()

    def snapshot() -> bool:
        """Print one block of static lines; return True when every file is complete."""
        lines = []
        all_complete = True
        for master, label in zip(masters, labels):
            avail, total = get_frame_counts(master, options.method)
            lines.append(_snapshot(label, avail, total))
            if total <= 0 or avail < total:
                all_complete = False
        print("\n".join(lines))
        return all_complete

    if not live:
        try:
            while True:
                if snapshot() and not options.watch_complete:
                    break
                if options.once:
                    break
                time.sleep(options.interval)
        except KeyboardInterrupt:
            pass
        return

    # One bar per master. The declared total is unknown until the first successful
    # read, so start unbounded and fill it in as soon as we know it.
    bars = [
        tqdm(
            total=None,
            desc=label,
            unit=UNIT,
            position=index,
            leave=True,
            dynamic_ncols=True,
            file=sys.stdout,
        )
        for index, label in enumerate(labels)
    ]
    try:
        while True:
            all_complete = True
            for master, bar in zip(masters, bars):
                avail, total = get_frame_counts(master, options.method)
                if total > 0 and bar.total != total:
                    bar.total = total
                if avail >= bar.n:
                    bar.update(avail - bar.n)
                else:
                    # Shouldn't happen, but never let a bar run backwards silently.
                    bar.reset(total=bar.total)
                    bar.update(avail)
                bar.refresh()
                if total <= 0 or avail < total:
                    all_complete = False
            if all_complete and not options.watch_complete:
                break
            time.sleep(options.interval)
    except KeyboardInterrupt:
        # Leave the last rendered frame visible and exit quietly.
        pass
    finally:
        for bar in bars:
            bar.close()


if __name__ == "__main__":
    run()
