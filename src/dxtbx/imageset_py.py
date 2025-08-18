from __future__ import annotations

from collections.abc import Callable
from typing import Any, Iterable, Protocol, Type

from dxtbx.format.Format import Format
from dxtbx.imageset import ExternalLookup
from dxtbx.masking import GoniometerShadowMasker


class Reader(Protocol):
    def read(self, index: int | None) -> Any: ...
    def is_single_file_reader(self) -> bool: ...
    def paths(self) -> list[str]: ...
    def master_path(self) -> str: ...
    def identifiers(self) -> list[str]: ...
    def __len__(self) -> int: ...


class ImageSetData:
    def __init__(
        self,
        reader: Reader,
        masker: GoniometerShadowMasker | Callable[[], GoniometerShadowMasker | None],
        template: str | None = None,
        vendor: str | None = None,
        params: dict | None = None,
        format: Type[Format] | None = None,
    ):
        self.reader = reader
        self._masker = masker
        self.template = template or None
        self.vendor = vendor or None
        self.format = format
        self.params = params or {}

        self._reject = [False for _ in range(len(self.reader))]
        self._beams: list[int | None] = [None for _ in range(len(self.reader))]
        self._detectors: list[int | None] = [None for _ in range(len(self.reader))]
        self._goniometers: list[int | None] = [None for _ in range(len(self.reader))]
        self._scans: list[int | None] = [None for _ in range(len(self.reader))]

        self.external_lookup = ExternalLookup()

    def __len__(self) -> int:
        return len(self.reader)

    @property
    def masker(self) -> GoniometerShadowMasker | None:
        if callable(self._masker):
            self._masker = self._masker()
        return self._masker

    def has_dynamic_mask(self) -> bool:
        return self.masker is not None

    def get_data(self, index: int | None):
        # Unsure if | None is needed?
        raise NotImplementedError

    def has_single_file_reader(self) -> bool:
        return self.reader.is_single_file_reader()

    def get_path(self):
        return self.reader.paths()[0]

    def get_master_path(self):
        return self.reader.master_path()

    def get_image_identifier(self, index: int) -> str:
        return self.reader.identifiers()[index]

    def mark_for_rejection(self, index: int, reject: bool) -> None:
        self._reject[index] = reject

    def is_marked_for_rejection(self, index: int) -> bool:
        return self._reject[index]

    def get_reject_list(self) -> list[bool]:
        return self._reject

    def set_reject_list(self, rejects: Iterable[bool]):
        if len(rejects) != len(self._reject):
            raise ValueError(
                f"Got reject list of {len(rejects)} items instead of {len(self._reject)}"
            )
        self._reject = list(rejects)

    def get_beam(self, index: int) -> int:
        item = self._beams[index]
        assert item is not None
        return item

    def set_beam(self, value: int, index: int):
        self._beams[index] = value

    def get_detector(self, index: int) -> int:
        item = self._detectors[index]
        assert item is not None
        return item

    def set_detector(self, value, index):
        self._detectors[index] = value

    def get_scan(self, index: int) -> int:
        item = self._scans[index]
        assert item is not None
        return item

    def set_scan(self, value, index):
        self._scans[index] = value

    def get_goniometer(self, index: int) -> int:
        item = self._goniometers[index]
        assert item is not None
        return item

    def set_goniometer(self, value, index):
        self._goniometers[index] = value

    def get_template(self) -> str | None:
        return self.template

    def set_template(self, value: str | None) -> None:
        self.template = value

    def get_vendor(self) -> str | None:
        return self.vendor

    def set_vendor(self, value: str | None) -> None:
        self.vendor = value

    def get_params(self) -> dict:
        return self.params

    def set_params(self, value: dict) -> None:
        self.params = value

    def get_format_class(self) -> Type[Format] | None:
        return self.format

    def set_format_class(self, value: Type[Format] | None) -> None:
        self.format = value

    def partial_data(self, reader: Reader, first: int, last: int) -> ImageSetData:
        raise NotImplementedError
