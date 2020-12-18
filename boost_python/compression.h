#ifndef DXTBX_COMPRESSION
#define DXTBX_COMPRESSION
#include <cstddef>
#include <vector>

namespace dxtbx { namespace boost_python {
  void cbf_decompress(const char*, std::size_t, int*);
  std::vector<char> cbf_compress(const int*, const std::size_t&);
}}  // namespace dxtbx::boost_python

#endif
