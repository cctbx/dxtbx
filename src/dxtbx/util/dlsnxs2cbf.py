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
    result.append("data_%06d" % (nn + 1))
    result.append("")
    result.append(
        """_array_data.header_convention "PILATUS_1.2"
_array_data.header_contents
;"""
    )
    result.append(
        f"# Detector: EIGER 2XE 16M S/N {nxdetector.serial_number} {nxinstrument.name}"
    )
    result.append(f"# {timestamp}")
    result.append("# Pixel_size 75e-6 m x 75e-6 m")
    result.append("# Silicon sensor, thickness 0.000450 m")
    result.append(f"# Exposure_time {nxdetector.count_time:.5f} s")
    result.append(f"# Exposure_period {nxdetector.count_time:.5f} s")
    result.append("# Tau = 1e-9 s")
    result.append("# Count_cutoff 65535 counts")
    result.append("# Threshold_setting: 0 eV")
    result.append("# Gain_setting: mid gain (vrf = -0.200)")
    result.append("# N_excluded_pixels = 0")
    result.append("# Excluded_pixels: badpix_mask.tif")
    result.append("# Flat_field: (nil)")
    result.append(f"# Wavelength {nxbeam.incident_wavelength:.5f} A")
    result.append(f"# Detector_distance {distance / 1000.0:.5f} m")
    result.append(
        f"# Beam_xy ({nxdetector.beam_center_x:.2f}, {nxdetector.beam_center_y:.2f}) pixels"
    )
    result.append("# Flux 0.000000")
    result.append("# Filter_transmission %.3f" % A)
    result.append("# Detector_2theta 0.0000 deg.")
    result.append("# Polarization 0.990")
    result.append("# Alpha 0.0000 deg.")
    for t in dependency_chain:
        name = t.path.split("/")[-1].capitalize()
        # Find the first varying rotation axis
        if (
            t.transformation_type == "rotation"
            and len(t) > 1
            and not np.all(t[()] == t[0])
        ):
            result.append(f"# {name} {t[nn].to('degree'):.4f} deg.")
            result.append(f"# {name}_increment {(t[1] - t[0]).to('degree'):.4f} deg")
            result.append(f"# Start_angle {t[nn].to('degree'):.4f} deg.")
            result.append(f"# Angle_increment {(t[1] - t[0]).to('degree'):.4f} deg")
        else:
            result.append(f"# {name} {t[0].to('degree'):.4f} deg.")
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
        static_mask = dxtbx.nexus.get_static_mask(nxdetector)[0]
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
            data.as_1d().set_selected(static_mask.as_1d(), -1)
            compressed = compress(data)

            mime = """

_array_data.data
;
--CIF-BINARY-FORMAT-SECTION--
Content-Type: application/octet-stream;
     conversions="x-CBF_BYTE_OFFSET"
Content-Transfer-Encoding: BINARY
X-Binary-Size: %d
X-Binary-ID: 1
X-Binary-Element-Type: "signed 32-bit integer"
X-Binary-Element-Byte-Order: LITTLE_ENDIAN
X-Binary-Number-of-Elements: %d
X-Binary-Size-Fastest-Dimension: %d
X-Binary-Size-Second-Dimension: %d
X-Binary-Size-Padding: 4095

""" % (
                len(compressed),
                data.size(),
                data.focus()[1],
                data.focus()[0],
            )

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
