#include "compression.h"
#include "dxtbx/error.h"
#include <assert.h>

typedef union {
  char b[2];
  short s;
} union_short;

typedef union {
  char b[4];
  int i;
} union_int;

void byte_swap_short(char *b) {
  char c;
  c = b[0];
  b[0] = b[1];
  b[1] = c;
  return;
}

void byte_swap_int(char *b) {
  char c;
  c = b[0];
  b[0] = b[3];
  b[3] = c;
  c = b[1];
  b[1] = b[2];
  b[2] = c;
  return;
}

bool little_endian() {
  int i = 0x1;
  char b = ((union_int *)&i)[0].b[0];
  if (b == 0) {
    return false;
  } else {
    return true;
  }
}

std::vector<char> dxtbx::boost_python::cbf_compress(const int *values,
                                                    const std::size_t &sz) {
  std::vector<char> packed(0);
  int current = 0;
  int delta, i;
  unsigned int j;
  bool le = little_endian();
  short s;
  char c;
  char *b;

  for (j = 0; j < sz; j++) {
    delta = values[j] - current;

    if ((-0x7f <= delta) && (delta < 0x80)) {
      c = (char)delta;
      packed.push_back(c);
      current += delta;
      continue;
    }

    packed.push_back(-0x80);

    if ((-0x7fff <= delta) && (delta < 0x8000)) {
      s = (short)delta;
      b = ((union_short *)&s)[0].b;

      if (!le) {
        byte_swap_short(b);
      }

      packed.push_back(b[0]);
      packed.push_back(b[1]);
      current += delta;
      continue;
    }

    s = -0x8000;
    b = ((union_short *)&s)[0].b;

    if (!le) {
      byte_swap_short(b);
    }

    packed.push_back(b[0]);
    packed.push_back(b[1]);

    assert(delta != -0x8000000);

    i = delta;
    b = ((union_int *)&i)[0].b;

    if (!le) {
      byte_swap_int(b);
    }

    packed.push_back(b[0]);
    packed.push_back(b[1]);
    packed.push_back(b[2]);
    packed.push_back(b[3]);
    current += delta;
  }

  return packed;
}

unsigned int dxtbx::boost_python::cbf_decompress(const char *packed,
                                                 std::size_t packed_sz,
                                                 int *values,
                                                 std::size_t values_sz) {
  int current = 0;
  int *original = values;
  unsigned int j = 0;
  short s;
  char c;
  int i;
  bool le = little_endian();

  while ((j < packed_sz) && ((values - original) < values_sz)) {
    c = packed[j];
    j += 1;

    if (c != -0x80) {
      current += c;
      *values = current;
      values++;
      continue;
    }

    DXTBX_ASSERT(j + 1 < packed_sz);
    ((union_short *)&s)[0].b[0] = packed[j];
    ((union_short *)&s)[0].b[1] = packed[j + 1];
    j += 2;

    if (!le) {
      byte_swap_short((char *)&s);
    }

    if (s != -0x8000) {
      current += s;
      *values = current;
      values++;
      continue;
    }

    DXTBX_ASSERT(j + 3 < packed_sz);
    ((union_int *)&i)[0].b[0] = packed[j];
    ((union_int *)&i)[0].b[1] = packed[j + 1];
    ((union_int *)&i)[0].b[2] = packed[j + 2];
    ((union_int *)&i)[0].b[3] = packed[j + 3];
    j += 4;

    if (!le) {
      byte_swap_int((char *)&i);
    }

    current += i;
    *values = current;
    values++;
  }

  return values - original;
}

inline uint32_t read_uint32_from_bytearray(const char *buf) {
  // `char` can be signed or unsigned depending on the platform.
  // For bit shift operations, we need unsigned values.
  // If `char` on the platform is signed, converting directly to "unsigned int" can
  // produce huge numbers because modulo 2^n is taken by the integral conversion
  // rules. Thus, we have to explicitly cast to `unsigned char` first (but is the
  // result valid in platforms that don't use two's complement?). Then the automatic
  // integral promotion converts them to `int`.

  return ((unsigned char)buf[0]) | (((unsigned char)buf[1]) << 8)
         | (((unsigned char)buf[2]) << 16) | (((unsigned char)buf[3]) << 24);
}

inline uint16_t read_uint16_from_bytearray(const char *buf) {
  return ((unsigned char)buf[0]) | ((unsigned char)buf[1] << 8);
}

void dxtbx::boost_python::TY6_decompress(int *const ret,
                                         const char *const buf_data,
                                         const char *const buf_offsets,
                                         const int slow,
                                         const int fast) {
  const size_t BLOCKSIZE = 8;             // Codes below assume this is at most 8
  const signed int SHORT_OVERFLOW = 127;  // after 127 is subtracted
  const signed int LONG_OVERFLOW = 128;

  const size_t nblock = (fast - 1) / (BLOCKSIZE * 2);
  const size_t nrest = (fast - 1) % (BLOCKSIZE * 2);

  for (size_t iy = 0; iy < slow; iy++) {
    size_t ipos = read_uint32_from_bytearray(buf_offsets + iy * sizeof(uint32_t));
    size_t opos = fast * iy;

    // Values from -127 to +126 (inclusive) are stored with an offset of 127
    // as 0 to 253. 254 and 255 mark short and long overflows.
    // Other values ("overflows") are represented in two's complement.

    int firstpx = (unsigned char)buf_data[ipos++] - 127;
    if (firstpx == LONG_OVERFLOW) {
      // FIXME: this is an unsigned to signed conversion and implementation dependent.
      firstpx = (signed int)read_uint32_from_bytearray(buf_data + ipos);
      ipos += 4;
    } else if (firstpx == SHORT_OVERFLOW) {
      firstpx = (signed short)read_uint16_from_bytearray(buf_data + ipos);
      ipos += 2;
    }
    ret[opos++] = firstpx;

    // For every two blocks
    for (int k = 0; k < nblock; k++) {
      const size_t bittypes = buf_data[ipos++];
      const size_t nbits[2] = {bittypes & 15, (bittypes >> 4) & 15};

      // One pixel is stored using `nbit` bits.
      // Although `nbit` itself is stored using 4 bits,
      // only values 1 (0001b) to 8 (1000b) are allowed.
      // Negative values are encoded as follows. (Not 2's complement!)
      // - When nbit = 1, the pixel value is 0 or 1
      // - When nbit = 2, the pixel value is -1, 0, 1, 2
      // - When nbit = 3, the pixel value is -3, -2, 1, 0, 1, 2, 3, 4
      // - When nbit - 8, the pixel value is -127, -126, ...,
      //   127 (== // SHORT_OVERFLOW_SIGNED), 128 (== LONG_OVERFLOW_SIGNED)

      // Load values
      for (int i = 0; i < 2; i++) {
        const size_t nbit = nbits[i];
        assert(nbit >= 0 && nbit <= 8);

        int zero_at = 0;
        if (nbit > 1) {
          zero_at = (1 << (nbit - 1)) - 1;
        }

        // Since nbit is at most 8, 8 * 8 (= BLOCKSIZE) = 64 bits are sufficient.
        unsigned long long v = 0;
        for (int j = 0; j < nbit; j++) {
          // Implicit promotion is only up to 32 bits, not 64 bits so we have to be
          // explicit.
          v |= (long long)((unsigned char)buf_data[ipos++]) << (BLOCKSIZE * j);
        }

        const unsigned long long mask = (1 << nbit) - 1;
        for (int j = 0; j < BLOCKSIZE; j++) {
          ret[opos++] = ((v >> (nbit * j)) & mask) - zero_at;
        }
      }

      // Apply differences. Load more values when overflown.
      for (size_t i = opos - 2 * BLOCKSIZE; i < opos; i++) {
        int offset = ret[i];

        if (offset == LONG_OVERFLOW) {
          offset = (signed int)read_uint32_from_bytearray(buf_data + ipos);
          ipos += 4;
        } else if (offset == SHORT_OVERFLOW) {
          offset = (signed short)read_uint16_from_bytearray(buf_data + ipos);
          ipos += 2;
        }

        ret[i] = offset + ret[i - 1];
      }
    }

    for (int i = 0; i < nrest; i++) {
      int offset = (unsigned char)buf_data[ipos++] - 127;

      if (offset == LONG_OVERFLOW) {
        offset = (signed int)read_uint32_from_bytearray(buf_data + ipos);
        ipos += 4;
      } else if (offset == SHORT_OVERFLOW) {
        offset = (signed short)read_uint16_from_bytearray(buf_data + ipos);
        ipos += 2;
      }

      ret[opos] = ret[opos - 1] + offset;
      opos++;
    }
  }
}
