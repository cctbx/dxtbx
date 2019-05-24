from __future__ import absolute_import, division, print_function

try:
    import dials.precommitbx.nagger

    dials.precommitbx.nagger.nag()
except ImportError:
    pass
