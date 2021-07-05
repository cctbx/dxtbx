import sys

from iotbx.detectors.bruker import BrukerImage

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format


class FormatBruker(Format):
    """cctbx has no authoritative sources describing the Bruker format.
    Positive identification:  when listing the header out in 16-byte chunks,
    there are numerous chunks of the form
    KEYNAME:xxxxxxxx
    KEYNAME:xxxxxxxx
    ...in other words a 7-upper case character (or space) keyname followed by
    a colon, starting at the beginning of every fifth chunk.
    Here we will take a series of 12 of these in the first 1024 characters
    as a positive fingerprint of Bruker format.
    """

    @staticmethod
    def read_header_lines(image_path, max_bytes=512 * 50):
        """Attempt to read the whole Bruker header into a list of lines of 80 char
        length. Finish reading lines when either the amount of data determined by
        HDRBLKS is consumed, or if max_bytes is reached, or a line contains a
        non-ASCII character. HDRBLKS is interpreted as being an integer number
        of 512 byte blocks. HDRBLKS must be divisible by 5, so that the total
        number of bytes in the header is divisible by 80."""

        hdr_lines = []
        max_bytes = (max_bytes // 80) * 80
        with FormatBruker.open_file(image_path, "rb") as f:
            while True:

                # read a line and look for HDRBLKS
                line = f.read(80)
                if not line:
                    break
                try:
                    line = line.decode("ascii")
                except UnicodeDecodeError:
                    break
                if line.startswith("HDRBLKS"):
                    try:
                        val = int(line.split(":", 1)[1])
                        if val % 5 == 0:
                            max_bytes = min(val * 512, max_bytes)
                    except ValueError:
                        pass

                # looks like a valid line
                hdr_lines.append(line)

                if max_bytes is not None and f.tell() >= max_bytes:
                    break

        return hdr_lines

    @staticmethod
    def parse_header(header_lines):
        header_dic = {}

        for l in header_lines:
            separator = l.find(":")
            if separator == -1:
                continue
            k = l[:separator].strip()
            v = l[separator + 1 :].strip()
            if k in header_dic:
                header_dic[k] = header_dic[k] + "\n" + v
            else:
                header_dic[k] = v

        return header_dic

    @staticmethod
    def understand(image_file):
        try:
            with FormatBruker.open_file(image_file, "rb") as fh:
                tag = fh.read(1024).decode("latin-1", "replace")
        except OSError:
            return False
        matches = 0
        for x in range(0, 1024, 80):
            word = tag[x : x + 16]
            if len(word) != 16:
                return False
            if word[7] == ":" and word[0:7].isupper():
                matches += 1
                if matches >= 12:
                    return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        super().__init__(str(image_file), **kwargs)

    def detectorbase_start(self):
        pass

    def _start(self):
        """Open the image file, read the image header, copy the key / value
        pairs into an internal dictionary self._header_dictionary along with
        the length of the header in bytes self._header_size."""
        self.detectorbase = BrukerImage(self._image_file)
        # self.detectorbase.readHeader() #unnecessary for the Bruker specialization

    def _goniometer(self):
        if self.detectorbase.parameters["OSC_RANGE"] > 0:
            return self._goniometer_factory.single_axis()
        else:
            return self._goniometer_factory.single_axis_reverse()

    def _detector(self):
        """Return a model for a simple detector"""

        # At present, ignore non-zero two theta for the dxtbx model
        # XXX Return to this issue later.
        # twotheta = self.detectorbase.parameters["TWOTHETA"]

        return self._detector_factory.simple(
            sensor="CCD",
            distance=self.detectorbase.distance,
            beam_centre=(self.detectorbase.beamx, self.detectorbase.beamy),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(self.detectorbase.pixel_size, self.detectorbase.pixel_size),
            image_size=(self.detectorbase.size1, self.detectorbase.size2),
            trusted_range=(0, self.detectorbase.saturation),
            mask=[],
        )  # a list of dead rectangles

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(self.detectorbase.wavelength)

    def _scan(self):
        """Return the scan information for this image."""

        return self._scan_factory.single_file(
            filename=self._image_file,
            # It's not at all clear how to recover the exposure time from the header
            # or even whether it is recorded.
            # XXX Here it will simply be set to a default number.
            exposure_times=1,
            osc_start=self.detectorbase.osc_start,
            osc_width=abs(self.detectorbase.deltaphi),
            epoch=None,
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatBruker.understand(arg))
