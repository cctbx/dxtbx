from __future__ import annotations

import importlib.metadata
import logging

try:
    profile_entry_points = importlib.metadata.entry_points()["dxtbx.profile_model"]
except KeyError:
    profile_entry_points = []


class ProfileModelFactory:
    """
    A factory to create a profile model
    """

    @staticmethod
    def from_dict(obj):
        """
        Given a dictionary, convert to a profile model
        """
        if obj is None:
            return None
        for entry_point in profile_entry_points:
            if entry_point.name == obj["__id__"]:
                return entry_point.load().from_dict(obj)
        logging.getLogger("dxtbx.model.profile").warn(
            "No profile class %s registered" % obj["__id__"]
        )
        print(
            "dxtbx.model.profile: WARNING: No profile class %s registered"
            % obj["__id__"]
        )
        return None
