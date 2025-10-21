from abc import ABC, abstractmethod
import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Union
import zmq


def format_address(
    ip_address: Optional[str] = None,
    port: Optional[int] = None,
    ipc_path: Optional[Union[str, Path]] = None,
) -> str:
    """
    Format a ZeroMQ address string, automatically detecting TCP or IPC protocol.

    Args:
        ip_address: IP address for TCP connections. If None, uses "*" for binding
        port: Port number for TCP connections
        ipc_path: File path for IPC connections

    Returns:
        Formatted ZeroMQ address string

    Raises:
        ValueError: If invalid parameters are provided

    Examples:
        >>> format_address("192.168.1.1", 5555)
        "tcp://192.168.1.1:5555"

        >>> format_address(port=5555)  # Binding
        "tcp://*:5555"

        >>> format_address(ipc_path="/tmp/socket")
        "ipc:///tmp/socket"

        >>> format_address(ipc_path="socket")  # Relative path
        "ipc://socket"
    """
    # Check for conflicting parameters
    if ipc_path is not None and port is not None:
        raise ValueError("Cannot specify both ipc_path and port - choose TCP or IPC")

    if ipc_path is not None:
        # IPC protocol
        if not ipc_path:  # Empty string check
            raise ValueError("ipc_path cannot be empty")

        path_str = str(ipc_path).strip()
        return f"ipc://{path_str}"

    elif port is not None:
        # TCP protocol
        if not isinstance(port, int):
            raise ValueError(f"port must be an integer, got: {type(port)}")
            if port <= 0 or port > 65535:
                raise ValueError(
                    f"Port must be a positive integer between 1-65535, got: {port}"
                )

        host = (ip_address or "*").strip() if ip_address else "*"
        return f"tcp://{host}:{port}"

    else:
        raise ValueError("Must provide either port (for TCP) or ipc_path (for IPC)")


class StreamClass(ABC):
    def __init__(
        self,
        port=None,
        ports=None,
        ip_address=None,
        socket_library=None,
        socket_type=None,
        socket_mode=None,
        zmq_context=None,
        rcvhwm=None,
        rcvbuf=None,
    ):
        if socket_mode == "connect":
            if ports is None:
                assert port and ip_address
                self._address = format_address(ip_address, port)
            elif ports:
                assert ip_address
                self._addresses = [
                    format_address(ip_address, port_i) for port_i in ports
                ]
        elif socket_mode == "bind":
            if ports is None:
                self._address = format_address(port=port)
            else:
                self._addresses = [format_address(port=port_i) for port_i in ports]
        else:
            assert socket_mode is None
            self._address = None

        if socket_library == "zeromq":
            self.socket = zmq_context.socket(socket_type)
            if rcvhwm:
                self.socket.setsockopt(zmq.RCVHWM, rcvhwm)
            if rcvbuf:
                self.socket.setsockopt(zmq.RCVBUF, rcvbuf)
            if socket_mode == "connect":
                if ports is None:
                    self.socket.connect(self._address)
                elif ports:
                    for address in self._addresses:
                        self.socket.connect(address)
            elif socket_mode == "bind":
                if ports is None:
                    self.socket.bind(self._address)
                else:
                    for address in self._addresses:
                        self.socket.bind(address)
        elif socket_library == "nng":
            assert False
        else:
            assert socket_library is None
            self.socket = None

    @abstractmethod
    def recv(self, copy: bool = True) -> bytes:
        """Receive a message from the stream"""
        pass

    @abstractmethod
    def decode(self, encoded_message) -> dict:
        """Decode the message"""
        pass

    @abstractmethod
    def handle_start_message(
        self, encoded_message=None, message=None, reference_experiment=None
    ):
        """
        Convert a start message into
            1: file writer phil parameters
                dxtbx.format.nxmx_writer.phil_scope
            2: An ExperimentList
        """
        pass

    @abstractmethod
    def get_data(self, message, **kwargs):
        """Convert an image message to a numpy array"""
        pass

    @abstractmethod
    def get_reader(self, image_data, **kwargs):
        from dials.array_family import flex
        from dxtbx.imageset import StreamReader

        image_data = flex.double(image_data)
        return StreamReader([image_data])
