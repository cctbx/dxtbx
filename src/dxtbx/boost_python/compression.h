#ifndef DXTBX_COMPRESSION
#define DXTBX_COMPRESSION
#include <cstddef>
#include <vector>

namespace dxtbx { namespace boost_python {
  unsigned int cbf_decompress(const char*, std::size_t, int*, const std::size_t);
  std::vector<char> cbf_compress(const int*, const std::size_t&);
  void TY6_decompress(int* const,
                      const char* const,
                      const char* const,
                      const int,
                      const int);
}}  // namespace dxtbx::boost_python

#endif
