from __future__ import annotations

from typing import Any, Union

from scitbx.array_family import flex

ImageData = Union[flex.int, flex.double, flex.float]

class ExternalLookup:
    def __init__(self, *args, **kwargs) -> None: ...
    def __reduce__(self) -> Any: ...
    @property
    def dx(self) -> Any: ...
    @property
    def dy(self) -> Any: ...
    @property
    def gain(self) -> Any: ...
    @property
    def mask(self) -> Any: ...
    @property
    def pedestal(self) -> Any: ...

class ExternalLookupItemBool:
    data: Any
    filename: Any
    def __init__(self, *args, **kwargs) -> None: ...
    def __reduce__(self) -> Any: ...

class ExternalLookupItemDouble:
    data: Any
    filename: Any
    def __init__(self, *args, **kwargs) -> None: ...
    def __reduce__(self) -> Any: ...

class ImageGrid(ImageSet):
    def __init__(self, *args, **kwargs) -> None: ...
    def from_imageset(self, *args, **kwargs) -> Any: ...
    def get_grid_size(self) -> Any: ...
    def __getinitargs__(self) -> Any: ...
    def __reduce__(self) -> Any: ...

class ImageSequence(ImageSet):
    def __init__(self, *args, **kwargs) -> None: ...
    def complete_set(self) -> Any: ...
    def get_array_range(self) -> Any: ...
    def get_beam(self) -> Any: ...
    def get_detector(self) -> Any: ...
    def get_goniometer(self) -> Any: ...
    def get_scan(self) -> Any: ...
    def partial_set(self, *args, **kwargs) -> Any: ...
    def set_beam(self, boost) -> Any: ...
    def set_detector(self, boost) -> Any: ...
    def set_goniometer(self, boost) -> Any: ...
    def set_scan(self, boost) -> Any: ...
    def update_detector_px_mm_data(self) -> Any: ...
    def __getinitargs__(self) -> Any: ...
    def __reduce__(self) -> Any: ...

class ImageSet:
    def __init__(boost, dxtbx) -> None: ...
    def as_imageset(self) -> Any: ...
    def clear_cache(self) -> Any: ...
    def complete_set(self) -> Any: ...
    def data(self) -> Any: ...
    def get_beam(self) -> Any: ...
    def get_corrected_data(self, int) -> Any: ...
    def get_detector(self) -> Any: ...
    def get_gain(self, int) -> Any: ...
    def get_goniometer(self) -> Any: ...
    def get_image_identifier(self, int) -> Any: ...
    def get_mask(self, int) -> Any: ...
    def get_path(self, int) -> Any: ...
    def get_pedestal(self, int) -> Any: ...
    def get_raw_data(self, index: int) -> tuple[ImageData]: ...
    def get_scan(self) -> Any: ...
    def has_dynamic_mask(self) -> Any: ...
    def indices(self) -> Any: ...
    def is_marked_for_rejection(self, int) -> Any: ...
    def mark_for_rejection(self, int, bool) -> Any: ...
    def partial_set(self, *args, **kwargs) -> Any: ...
    def set_beam(self, boost) -> Any: ...
    def set_detector(self, boost) -> Any: ...
    def set_goniometer(self, boost) -> Any: ...
    def set_scan(self, boost) -> Any: ...
    def size(self) -> Any: ...
    def update_detector_px_mm_data(self) -> Any: ...
    def __eq__(self, other) -> Any: ...
    def __getinitargs__(self) -> Any: ...
    def __len__(self) -> Any: ...
    def __ne__(self, other) -> Any: ...
    def __reduce__(self) -> Any: ...
    @property
    def external_lookup(self) -> Any: ...

class ImageSetData:
    def __init__(self, *args, **kwargs) -> None: ...
    def get_beam(self, int) -> Any: ...
    def get_data(self, int) -> Any: ...
    def get_detector(self, int) -> Any: ...
    def get_format_class(self) -> Any: ...
    def get_goniometer(self, int) -> Any: ...
    def get_image_identifier(self, int) -> Any: ...
    def get_master_path(self) -> Any: ...
    def get_params(self) -> Any: ...
    def get_path(self, int) -> Any: ...
    def get_scan(self, int) -> Any: ...
    def get_template(self) -> Any: ...
    def get_vendor(self) -> Any: ...
    def has_single_file_reader(self) -> Any: ...
    def is_marked_for_rejection(self, int) -> Any: ...
    def mark_for_rejection(self, int, bool) -> Any: ...
    def masker(self) -> Any: ...
    def reader(self) -> Any: ...
    def set_beam(self, boost, int) -> Any: ...
    def set_detector(self, boost, int) -> Any: ...
    def set_format_class(self, boost) -> Any: ...
    def set_goniometer(self, boost, int) -> Any: ...
    def set_params(self, boost) -> Any: ...
    def set_scan(self, boost, int) -> Any: ...
    def set_template(self, *args, **kwargs) -> Any: ...
    def set_vendor(self, *args, **kwargs) -> Any: ...
    def __getinitargs__(self) -> Any: ...
    def __getstate__(self) -> Any: ...
    def __reduce__(self) -> Any: ...
    def __setstate__(self, boost) -> Any: ...
    @property
    def external_lookup(self) -> Any: ...
