from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import zmq

if TYPE_CHECKING:
    from libtbx.phil import scope_extract

    from dxtbx.imageset import StreamReader
    from dxtbx.model import ExperimentList


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
        port: Optional[int] = None,
        ports: Optional[List[int]] = None,
        ip_address: Optional[str] = None,
        zmq_context: Optional[zmq.Context] = None,
        rcvhwm: Optional[int] = None,
        rcvbuf: Optional[int] = None,
        tcp_keepalive: bool = False,
    ) -> None:
        # A reader constructed without a ZeroMQ context is decoder-only: it
        # normalizes messages and holds no socket (port/ip are ignored). Live
        # readers always PULL-connect to the source; ``ports`` connects to
        # several source endpoints at once.
        if zmq_context is None:
            self._address = None
            self.socket = None
            return

        if ports is None:
            assert port and ip_address
            self._address = format_address(ip_address, port)
            addresses = [self._address]
        else:
            assert ip_address
            self._addresses = [format_address(ip_address, port_i) for port_i in ports]
            addresses = self._addresses

        self.socket = zmq_context.socket(zmq.PULL)
        self.socket.setsockopt(zmq.LINGER, 0)
        if rcvhwm:
            self.socket.setsockopt(zmq.RCVHWM, rcvhwm)
        if rcvbuf:
            self.socket.setsockopt(zmq.RCVBUF, rcvbuf)
        if tcp_keepalive:
            # Set before connect so the options apply to the connection: the kernel
            # then reaps a dead/orphaned peer (e.g. a consumer killed after an SSH
            # drop) instead of leaving it ESTABLISHED and silently stealing the
            # detector's round-robined stream.
            self.socket.setsockopt(zmq.TCP_KEEPALIVE, 1)
            self.socket.setsockopt(zmq.TCP_KEEPALIVE_IDLE, 30)
            self.socket.setsockopt(zmq.TCP_KEEPALIVE_INTVL, 10)
            self.socket.setsockopt(zmq.TCP_KEEPALIVE_CNT, 3)
        for address in addresses:
            self.socket.connect(address)

    def close_socket(self) -> None:
        self.socket.close()

    @abstractmethod
    def recv(self, copy: bool = True) -> bytes:
        """Receive a message from the stream"""
        pass

    @abstractmethod
    def decode(self, encoded_message: bytes) -> Dict[str, Any]:
        """Decode the message"""
        pass

    @abstractmethod
    def handle_start_message(
        self,
        message: Dict[str, Any],
        reference_experiment: Optional[ExperimentList] = None,
        sync_reference_geom: bool = True,
        wavelength: Optional[float] = None,
    ) -> Tuple[scope_extract, ExperimentList]:
        """
        Convert a start message into
            1: file writer phil parameters
                dxtbx.format.nxmx_writer.phil_scope
            2: An ExperimentList
        """
        pass

    @abstractmethod
    def get_data(
        self, message: Dict[str, Any], **kwargs: Any
    ) -> Tuple[Any, Optional[float]]:
        """Convert an image message to a numpy array"""
        pass

    def get_reader(self, image_data: Any, **kwargs: Any) -> StreamReader:
        from dials.array_family import flex

        from dxtbx.imageset import StreamReader

        image_data = flex.double(image_data)
        return StreamReader([image_data])

    def get_data_scale_factor(self, wavelength: float) -> float:
        """NeXus data_scale_factor for a frame at the given wavelength (Angstrom).

        The archived pixel data is stored in the detector's native units; the
        scale factor records what to multiply by on read to obtain photons
        (NeXus: corrected = (data + offset) * scaling_factor). The base class
        returns 1.0 (data is already in photons); detector APIs whose data is
        not in photons (e.g. LCLStreamer, keV) override this. The caller supplies
        the resolved wavelength so components never branch on API.
        """
        return 1.0

    def get_wavelength(self, message: dict, **calib: Any) -> Optional[float]:
        """Resolve the per-frame wavelength (Angstrom) from a decoded image message.

        This is the single authoritative wavelength resolver: every API-specific
        case (a per-shot spectrum, an eBeam photon energy, a wavelength PV, a
        configured fallback) is handled here so components make one call and never
        branch on which detector API is in use. ``calib`` carries the optional PHIL
        knobs (``spectrum_eV_per_pixel``, ``spectrum_eV_offset``,
        ``wavelength_offset``, ``wavelength_scale``, ``wavelength_fallback``).

        Base behavior: no per-frame wavelength (APIs that carry wavelength on the
        start message return None here, leaving the reference-beam wavelength in
        place). Readers whose image messages carry per-frame energy/spectrum
        override this with their fallback chain.
        """
        return None

    def get_wavelength_candidates(
        self, message: dict, spectrum: Optional[Any] = None, **calib: Any
    ) -> dict:
        """The per-source candidate wavelengths (Angstrom) that ``get_wavelength``
        chooses between, computed for every frame regardless of which one wins so the
        alternatives can be archived alongside the resolved ``incident_wavelength`` for
        diagnostics.

        Returns a dict mapping a short source name (e.g. ``"spectrum"``, ``"ebeam"``,
        ``"pv"``) to that source's raw wavelength in Angstrom; a source absent for this
        frame is omitted. ``spectrum`` is the already-built per-shot Spectrum (passed in
        so it is not decoded again). Base behavior: no candidates (the resolver returns
        the start-message reference). Readers whose image messages carry multiple
        wavelength sources override this.
        """
        return {}

    def get_spectrum(self, message: dict, **calib: Any) -> Optional[Any]:
        """Build a ``dxtbx.model.Spectrum`` for this frame, or None if unavailable.

        Base behavior: no per-frame spectrum. Readers whose image messages carry a
        spectrometer trace (or a 2D spectrometer image to collapse) override this.
        ``calib`` carries the spectrometer calibration PHIL (``spectrum_eV_per_pixel``,
        ``spectrum_eV_offset``).
        """
        return None

    def get_spectrometer_image(self, message: dict) -> Optional[Any]:
        """Return the raw per-frame spectrometer array (a 2D detector image) to archive
        as a second image stream, or None.

        Base behavior: none. Readers override this when the API carries a 2D
        spectrometer image whose full native-dimensionality pattern must be stored
        alongside the main detector frame (the 1D collapse used for the wavelength is
        a separate concern handled by get_spectrum). A 1D spectrometer trace returns
        None here because it is already preserved on NXbeam by the get_spectrum path.
        """
        return None

    def get_frame_data(
        self,
        message: dict,
        image_shape: Any = None,
        image_dtype: Any = None,
        calib: Optional[dict] = None,
    ) -> Tuple[Any, Optional[float], Optional[Any], Optional[Any]]:
        """Single per-frame entry point returning
        ``(image_data, wavelength, spectrum, spectrometer_image)``.

        Components call this once per frame so a reader can decode its payloads a
        single time and return mutually-consistent models. ``calib`` is the optional
        wavelength/spectrum PHIL dict. Base composition (for readers with no
        spectrometer) simply composes the individual accessors; readers that carry a
        spectrometer override this to decode it once.
        """
        calib = calib or {}
        image_data, _ = self.get_data(
            message, image_shape=image_shape, image_dtype=image_dtype
        )
        wavelength = self.get_wavelength(message, **calib)
        spectrum = self.get_spectrum(message, **calib)
        spectrometer_image = self.get_spectrometer_image(message)
        return image_data, wavelength, spectrum, spectrometer_image
