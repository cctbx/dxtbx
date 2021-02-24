#include "compression.h"
#include "dxtbx/error.h"

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

    if ((-0xf777 <= delta) && (delta < 0x8000)) {
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

    if ((-0x7fffffff <= delta) && (delta < 0x80000000)) {
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
      continue;
    }
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

    if (c != -128) {
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

    if (s != -32768) {
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
