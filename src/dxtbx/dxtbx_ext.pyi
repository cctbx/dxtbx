from boost_adaptbx.boost.python import streambuf
from scitbx.array_family import flex

def compress(array: flex.int) -> bytes: ...
def is_big_endian() -> bool: ...
def read_float32(file: streambuf, count: int) -> float: ...
def read_int16(file: streambuf, count: int) -> int: ...
def read_int32(file: streambuf, count: int) -> int: ...
def read_uint16(file: streambuf, count: int) -> int: ...
def read_uint16_bs(file: streambuf, count: int) -> int: ...
def read_uint32(file: streambuf, count: int) -> int: ...
def read_uint32_bs(file: streambuf, count: int) -> int: ...
def read_uint8(file: streambuf, count: int) -> int: ...
def uncompress(packed: bytes, slow: int, fast: int) -> flex.int: ...
