from __future__ import annotations

import binascii

import h5py
import hdf5plugin  # noqa; F401
import numpy as np

import dxtbx.model
import dxtbx.nexus.nxmx
from dxtbx.ext import compress


def compute_cbf_header(nxmx: dxtbx.nexus.nxmx.NXmx, nn: int):
    nxentry = nxmx.entries[0]
    nxsample = nxentry.samples[0]
    nxinstrument = nxentry.instruments[0]
    nxdetector = nxinstrument.detectors[0]
    nxbeam = nxinstrument.beams[0]

    distance = nxdetector.distance
    if distance is None:
        distance = dxtbx.nexus.get_dxtbx_detector(nxdetector, nxbeam)[0].get_distance()

    result = []

    A = nxinstrument.attenuators[0]["attenuator_transmission"][()]

    # this timestamp _should_ be in UTC - so can ignore timezone info - at
    # least for purposes here (causes errors for old versions of dxtbx reading
    # the data) - 19 chars needed
    timestamp = nxentry.start_time

    dependency_chain = dxtbx.nexus.nxmx.get_dependency_chain(nxsample.depends_on)

    result.append("###CBF: VERSION 1.5, CBFlib v0.7.8 - Eiger detectors")
    result.append("")
    result.append(f"data_{nn + 1:06d}")
    result.append("")
    result.append(
        """_array_data.header_convention "PILATUS_1.2"
_array_data.header_contents
;"""
    )
    result.append(
        f"# Detector: {nxdetector.description} S/N {nxdetector.serial_number} "
        f"{nxinstrument.short_name}"
    )
    result.append(f"# {timestamp}")
    result.append(
        f"# Pixel_size {nxdetector.modules[0].fast_pixel_direction[()]:~} "
        f"x {nxdetector.modules[0].slow_pixel_direction[()]:~}"
    )
    result.append(
        f"# {nxdetector.sensor_material} sensor, "
        f"thickness {nxdetector.sensor_thickness:~}"
    )
    result.append(f"# Exposure_time {nxdetector.count_time:.5f} s")
    result.append(f"# Exposure_period {nxdetector.count_time:.5f} s")
    result.append("# Tau = 1e-9 s")
    result.append(f"# Count_cutoff {nxdetector.saturation_value} counts")
    result.append("# Threshold_setting: 0 eV")
    result.append("# Gain_setting: mid gain (vrf = -0.200)")
    result.append("# N_excluded_pixels = 0")
    result.append(f"# Wavelength {nxbeam.incident_wavelength:.5f} A")
    result.append(f"# Detector_distance {distance / 1000.0:.5f} m")
    result.append(
        f"# Beam_xy ({nxdetector.beam_center_x:.2f}, {nxdetector.beam_center_y:.2f}) "
        "pixels"
    )
    result.append("# Flux 0.000000")
    result.append(f"# Filter_transmission {A:.3f}")

    transformations = {
        t.path.split("/")[-1].capitalize(): t for t in nxinstrument.transformations
    }
    two_theta = transformations.get("TWO_THETA", np.array(0))[()]
    result.append(f"# Detector_2theta {two_theta:.3f} deg.")

    result.append("# Polarization 0.990")

    rotations = {
        t.path.split("/")[-1].capitalize(): t
        for t in dependency_chain
        if t.transformation_type == "rotation"
    }

    if "PHI" in rotations and "KAPPA" in rotations:
        kappa_axis = rotations["KAPPA"].attrs["vector"]
        phi_axis = rotations["PHI"].attrs["vector"]
        alpha = np.degrees(np.arccos(np.dot(phi_axis, kappa_axis)))
    else:
        alpha = 0
    result.append(f"# Alpha {alpha:.3f} deg.")

    for name, r in rotations.items():
        # Find the first varying rotation axis
        if len(r) > 1 and not np.all(r[()] == r[0]):
            result.append(f"# {name} {r[nn].to('degree'):.4f} deg.")
            result.append(f"# {name}_increment {(r[1] - r[0]).to('degree'):.4f} deg")
            result.append(f"# Start_angle {r[nn].to('degree'):.4f} deg.")
            result.append(f"# Angle_increment {(r[1] - r[0]).to('degree'):.4f} deg")
        else:
            result.append(f"# {name} {r[0].to('degree'):.4f} deg.")
            result.append(f"# {name}_increment 0.0000 deg")

    result.append("# Oscillation_axis X.CW")
    result.append("# N_oscillations 1")
    result.append(";")

    return "\n".join(result)


def make_cbf(in_name, template):
    with h5py.File(in_name, "r") as f:
        start_tag = binascii.unhexlify("0c1a04d5")

        nxmx = dxtbx.nexus.nxmx.NXmx(f)
        nxsample = nxmx.entries[0].samples[0]
        nxinstrument = nxmx.entries[0].instruments[0]
        nxdetector = nxinstrument.detectors[0]
        nxdata = nxmx.entries[0].data[0]

        dependency_chain = dxtbx.nexus.nxmx.get_dependency_chain(nxsample.depends_on)
        scan_axis = None
        for t in dependency_chain:
            # Find the first varying rotation axis
            if (
                t.transformation_type == "rotation"
                and len(t) > 1
                and not np.all(t[()] == t[0])
            ):
                scan_axis = t
                break

        if scan_axis is None:
            # Fall back on the first varying axis of any type
            for t in dependency_chain:
                if len(t) > 1 and not np.all(t[()] == t[0]):
                    scan_axis = t
                    break

        if scan_axis is None:
            scan_axis = nxsample.depends_on

        num_images = len(scan_axis)
        static_mask = dxtbx.nexus.get_static_mask(nxdetector)
        bit_depth_readout = nxdetector.bit_depth_readout

        for j in range(num_images):
            header = compute_cbf_header(nxmx, j)
            data = dxtbx.nexus.get_raw_data(nxdata, nxdetector, j)[0]
            if bit_depth_readout:
                # if 32 bit then it is a signed int, I think if 8, 16 then it is
                # unsigned with the highest two values assigned as masking values
                if bit_depth_readout == 32:
                    top = 2 ** 31
                else:
                    top = 2 ** bit_depth_readout
                d1d = data.as_1d()
                d1d.set_selected(d1d == top - 1, -1)
                d1d.set_selected(d1d == top - 2, -2)

            # set the tile join regions to -1 - MOSFLM cares about this apparently
            if static_mask:
                data.as_1d().set_selected(static_mask[0].as_1d(), -1)
            compressed = compress(data)

            mime = f"""

_array_data.data
;
--CIF-BINARY-FORMAT-SECTION--
Content-Type: application/octet-stream;
     conversions="x-CBF_BYTE_OFFSET"
Content-Transfer-Encoding: BINARY
X-Binary-Size: {len(compressed)}
X-Binary-ID: 1
X-Binary-Element-Type: "signed 32-bit integer"
X-Binary-Element-Byte-Order: LITTLE_ENDIAN
X-Binary-Number-of-Elements: {data.size()}
X-Binary-Size-Fastest-Dimension: {data.focus()[1]}
X-Binary-Size-Second-Dimension: {data.focus()[0]}
X-Binary-Size-Padding: 4095

"""

            padding = (
                bytearray(4095)
                + b"""--CIF-BINARY-FORMAT-SECTION----
;"""
            )

            with open(template % (j + 1), "wb") as fout:
                print(template % (j + 1))
                fout.write(
                    ("".join(header) + mime).replace("\n", "\r\n").encode("latin-1")
                )
                fout.write(start_tag)
                fout.write(compressed)
                fout.write(padding)
