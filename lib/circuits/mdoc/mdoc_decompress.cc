// Copyright 2026 Google LLC.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "circuits/mdoc/mdoc_decompress.h"

#include <cstddef>
#include <cstdint>
#include <vector>

#include "util/log.h"
#include "zstd.h"

#if defined(__APPLE__)
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace proofs {

#if defined(__APPLE__)
// Decompress into an mmap-backed temporary file so the 150 MB output buffer
// is backed by the filesystem instead of anonymous heap memory.  iOS will
// page-out the mapped region under memory pressure rather than killing the
// process.
size_t decompress(std::vector<uint8_t>& /*unused*/, const uint8_t* compressed,
                  size_t compressed_len, int& fd_out, uint8_t*& map_out,
                  size_t capacity) {
  FILE* f = tmpfile();
  if (!f) {
    log(ERROR, "tmpfile() failed");
    return 0;
  }
  int fd = dup(fileno(f));
  fclose(f);
  if (fd < 0) {
    log(ERROR, "dup(fileno) failed");
    return 0;
  }
  if (ftruncate(fd, static_cast<off_t>(capacity)) != 0) {
    close(fd);
    log(ERROR, "ftruncate failed");
    return 0;
  }
  void* map =
      mmap(nullptr, capacity, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (map == MAP_FAILED) {
    close(fd);
    log(ERROR, "mmap failed");
    return 0;
  }
  size_t res = ZSTD_decompress(map, capacity, compressed, compressed_len);
  if (ZSTD_isError(res)) {
    munmap(map, capacity);
    close(fd);
    log(ERROR, "ZSTD_decompress failed: %s", ZSTD_getErrorName(res));
    return 0;
  }
  fd_out = fd;
  map_out = static_cast<uint8_t*>(map);
  return res;
}
#endif  // __APPLE__

// Decompress a circuit representation into a vector that has been reserved
// with size len.  The value len needs to be a good upper-bound estimate on
// the size of the uncompressed string.
size_t decompress(std::vector<uint8_t>& bytes, const uint8_t* compressed,
                  size_t compressed_len) {
  size_t res =
      ZSTD_decompress(bytes.data(), bytes.size(), compressed, compressed_len);

  if (ZSTD_isError(res)) {
    log(ERROR, "zlib.UncompressAtMost failed: %s", ZSTD_getErrorName(res));
    return 0;
  }
  return res;
}

}  // namespace proofs
