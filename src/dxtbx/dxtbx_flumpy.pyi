from typing import Union

import numpy as np

from scitbx.array_family import flex

# Simple, flat array types
FlexSimpleArray = Union[
    flex.bool,
    flex.uint8,
    flex.uint16,
    flex.uint32,
    flex.uint64,
    flex.size_t,
    flex.int,
    flex.int8,
    flex.int16,
    flex.int32,
    flex.int64,
    flex.double,
    flex.float,
    flex.complex_double,
]
FlexVecArray = Union[flex.vec2_double, flex.vec3_double, flex.vec3_int]
FlexArray = Union[FlexSimpleArray, FlexVecArray]

class Scuffer:
    def __init__(self, arg0: FlexArray) -> None: ...
    @property
    def base(self) -> FlexArray: ...

def from_numpy(arg0: np.ndarray) -> FlexArray: ...
def mat3_from_numpy(arg0: np.ndarray) -> flex.mat3_double: ...
def to_numpy(arg0: FlexArray) -> np.ndarray: ...
def vec_from_numpy(arg0: np.ndarray) -> FlexVecArray: ...
