from __future__ import absolute_import, division, print_function

import sys

from scitbx import matrix

from dxtbx.imageset import ImageSetFactory


def print_sequence(list_of_images):

    sequences = ImageSetFactory.new(list_of_images)

    for sequence in sequences:
        print(sequence.get_detector())
        print(sequence.get_beam())
        print(sequence.get_goniometer())
        print(sequence.get_scan())

        # compute the beam centre... in mm... w.r.t. fast, slow axis

        print("Derived quantities:")

        d = sequence.get_detector()[0]
        b = sequence.get_beam()

        o = matrix.col(d.get_origin())
        f = matrix.col(d.get_fast_axis())
        s = matrix.col(d.get_slow_axis())
        s0 = matrix.col(b.get_sample_to_source_direction())

        beam_offset = o - o.dot(s0) * s0
        print(
            "    beam centre (mm, fast, slow): %.2f %.2f"
            % (-beam_offset.dot(f), -beam_offset.dot(s))
        )


if __name__ == "__main__":
    if len(sys.argv) == 2:
        print_sequence(sys.argv[1])
    else:
        print_sequence(sys.argv[1:])
