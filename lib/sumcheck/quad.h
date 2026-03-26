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

#ifndef PRIVACY_PROOFS_ZK_LIB_SUMCHECK_QUAD_H_
#define PRIVACY_PROOFS_ZK_LIB_SUMCHECK_QUAD_H_

#include <stddef.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

#include "algebra/compare.h"
#include "algebra/poly.h"
#include "arrays/affine.h"
#include "arrays/eqs.h"
#include "util/ceildiv.h"
#include "util/panic.h"
#define DEFINE_STRONG_INT_TYPE(a, b) using a = b

#if defined(__APPLE__)
#include <sys/mman.h>
#include <unistd.h>

// File-backed mmap arena.  All allocations come from a single temp file so iOS
// can page-out the data under memory pressure (clean file-backed pages).
// The arena hands out slices from one contiguous mapping; slices outlive the
// arena because the mmap stays valid until the destructor.
class MmapArena {
  int fd_ = -1;
  uint8_t *base_ = nullptr;
  size_t capacity_ = 0;
  size_t used_ = 0;

 public:
  MmapArena() = default;
  ~MmapArena() {
    if (base_) munmap(base_, capacity_);
    if (fd_ >= 0) ::close(fd_);
  }
  MmapArena(const MmapArena &) = delete;
  MmapArena &operator=(const MmapArena &) = delete;

  bool init(size_t cap) {
    FILE *f = tmpfile();
    if (!f) return false;
    fd_ = dup(fileno(f));
    fclose(f);
    if (fd_ < 0) return false;
    capacity_ = cap;
    if (ftruncate(fd_, static_cast<off_t>(cap)) != 0) return false;
    base_ = static_cast<uint8_t *>(
        mmap(nullptr, cap, PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0));
    if (base_ == MAP_FAILED) {
      base_ = nullptr;
      return false;
    }
    return true;
  }

  // Return aligned pointer into the arena, or nullptr on overflow.
  void *alloc(size_t bytes, size_t align = 16) {
    size_t off = (used_ + align - 1) & ~(align - 1);
    if (off + bytes > capacity_) return nullptr;
    used_ = off + bytes;
    return base_ + off;
  }

  bool valid() const { return base_ != nullptr; }
};
#endif  // __APPLE__

// ------------------------------------------------------------
// Special-purpose sparse array for use with sumcheck
namespace proofs {
template <class Field>
class Quad {
  using Elt = typename Field::Elt;
  using T2 = Poly<2, Field>;

 public:
  // To save space when representing large circuits, quad_corner_t
  // is defined as uint32_t.  (Note that Elt probably imposes uint64_t
  // alignment, so struct corner has holes.)
  //
  // To make the narrowing explicit, define corner_t as a
  // Google-specific strong int.  Outside of Google, replace
  // this definition with a typedef.
  DEFINE_STRONG_INT_TYPE(quad_corner_t, uint32_t);

  struct corner {
    quad_corner_t g;     // "gate" variable
    quad_corner_t h[2];  // two "hand" variables
    Elt v;

    bool operator==(const corner& y) const {
      return g == y.g &&
             morton::eq(size_t(h[0]), size_t(h[1]), size_t(y.h[0]),
                        size_t(y.h[1])) &&
             v == y.v;
    }

    bool eqndx(const corner& y) const {
      return (g == y.g && h[0] == y.h[0] && h[1] == y.h[1]);
    }

    void canonicalize() {
      quad_corner_t h0 = h[0], h1 = h[1];
      h[0] = std::min<quad_corner_t>(h0, h1);
      h[1] = std::max<quad_corner_t>(h0, h1);
    }

    static bool compare(const corner& x, const corner& y, const Field& F) {
      if (morton::lt(size_t(x.h[0]), size_t(x.h[1]), size_t(y.h[0]),
                     size_t(y.h[1]))) {
        return true;
      } else if (morton::eq(size_t(x.h[0]), size_t(x.h[1]), size_t(y.h[0]),
                            size_t(y.h[1]))) {
        if (x.g < y.g) return true;
        if (x.g > y.g) return false;
        return elt_less_than(x.v, y.v, F);
      } else {
        return false;
      }
    }
  };

  using index_t = size_t;
  index_t n_;

  // On Apple/iOS, corner storage is a raw pointer that can be backed by either
  // heap memory (c_heap_) or a file-backed MmapArena.  The arena path lets iOS
  // page out circuit data under memory pressure instead of jetsaming the app.
  // On other platforms, c_ is a plain std::vector.
#if defined(__APPLE__)
  corner *c_;
  std::vector<corner> c_heap_;  // owns memory when arena is not used

  void move_to_arena(MmapArena &arena) {
    if (c_heap_.empty()) return;
    void *p = arena.alloc(c_heap_.size() * sizeof(corner), alignof(corner));
    if (!p) return;
    memcpy(p, c_heap_.data(), c_heap_.size() * sizeof(corner));
    c_ = static_cast<corner *>(p);
    c_heap_.clear();
    c_heap_.shrink_to_fit();
  }
#else
  std::vector<corner> c_;
#endif

  bool operator==(const Quad& y) const {
    return n_ == y.n_ && std::equal(c_, c_ + n_, y.c_);
  }

  explicit Quad(index_t n) : n_(n) {
#if defined(__APPLE__)
    c_heap_.resize(n);
    c_ = c_heap_.data();
#else
    c_.resize(n);
#endif
  }

#if defined(__APPLE__)
  // Arena-backed constructor: corners are allocated directly in the mmap arena,
  // bypassing the heap entirely.  Falls back to heap if arena is exhausted.
  Quad(index_t n, MmapArena &arena) : n_(n) {
    void *p = arena.alloc(n * sizeof(corner), alignof(corner));
    if (p) {
      c_ = static_cast<corner *>(p);
      memset(c_, 0, n * sizeof(corner));
    } else {
      c_heap_.resize(n);
      c_ = c_heap_.data();
    }
  }
#endif

  // no copies, but see clone() below
  Quad(const Quad& y) = delete;
  Quad(const Quad&& y) = delete;
  Quad operator=(const Quad& y) = delete;

  std::unique_ptr<Quad> clone() const {
    auto s = std::make_unique<Quad>(n_);
    for (index_t i = 0; i < n_; ++i) {
      s->c_[i] = c_[i];
    }
    return s;
  }

  void bind_h(const Elt& r, size_t hand, const Field& F) {
    index_t rd = 0, wr = 0;
    while (rd < n_) {
      corner cc;
      cc.g = quad_corner_t(0);
      cc.h[hand] = c_[rd].h[hand] >> 1;
      cc.h[1 - hand] = c_[rd].h[1 - hand];

      size_t rd1 = rd + 1;
      if (rd1 < n_ &&                                         //
          c_[rd].h[1 - hand] == c_[rd1].h[1 - hand] &&        //
          (c_[rd].h[hand] >> 1) == (c_[rd1].h[hand] >> 1) &&  //
          c_[rd1].h[hand] == c_[rd].h[hand] + quad_corner_t(1)) {
        // we have two corners.
        cc.v = affine_interpolation(r, c_[rd].v, c_[rd1].v, F);
        rd += 2;
      } else {
        // we have one corner and the other one is zero.
        if ((c_[rd].h[hand] & quad_corner_t(1)) == quad_corner_t(0)) {
          cc.v = affine_interpolation_nz_z(r, c_[rd].v, F);
        } else {
          cc.v = affine_interpolation_z_nz(r, c_[rd].v, F);
        }
        rd = rd1;
      }

      c_[wr++] = cc;
    }

    // shrink the array
    n_ = wr;
  }

  // Set zero coefficients to BETA, then bind to both
  // G0 and G1 and take the linear combination bind(G0) + alpha*bind(G1)
  void bind_g(size_t logv, const Elt* G0, const Elt* G1, const Elt& alpha,
              const Elt& beta, const Field& F) {
    size_t nv = size_t(1) << logv;
    auto dot = Eqs<Field>::raw_eq2(logv, nv, G0, G1, alpha, F);
    for (index_t i = 0; i < n_; ++i) {
      if (c_[i].v == F.zero()) {
        c_[i].v = beta;
      }
      F.mul(c_[i].v, dot[corner_t(c_[i].g)]);
      c_[i].g = quad_corner_t(0);
    }

    // coalesce any duplicates that we may have created
    coalesce(F);
  }

  // Optimized combined bind_g + bind_h, nondestructive
  Elt bind_gh_all(
      // G bindings
      size_t logv, const Elt G0[/*logv*/], const Elt G1[/*logv*/],
      const Elt& alpha, const Elt& beta,
      // H bindings
      size_t logw, const Elt H0[/*logw*/], const Elt H1[/*logw*/],
      // field
      const Field& F) const {
    size_t nv = size_t(1) << logv;
    auto eqg = Eqs<Field>::raw_eq2(logv, nv, G0, G1, alpha, F);

    size_t nw = size_t(1) << logw;
    Eqs<Field> eqh0(logw, nw, H0, F);
    Eqs<Field> eqh1(logw, nw, H1, F);

    Elt s{};

    for (index_t i = 0; i < n_; ++i) {
      Elt q(c_[i].v);
      if (q == F.zero()) {
        q = beta;
      }
      F.mul(q, eqg[corner_t(c_[i].g)]);
      F.mul(q, eqh0.at(corner_t(c_[i].h[0])));
      F.mul(q, eqh1.at(corner_t(c_[i].h[1])));
      F.add(s, q);
    }
    return s;
  }

  Elt scalar() {
    check(n_ == 1, "n_ == 1");
    check(c_[0].g == quad_corner_t(0), "c_[0].g == 0");
    check(c_[0].h[0] == quad_corner_t(0), "c_[0].h[0] == 0");
    check(c_[0].h[1] == quad_corner_t(0), "c_[0].h[1] == 0");
    return c_[0].v;
  }

  void canonicalize(const Field& F) {
    for (index_t i = 0; i < n_; ++i) {
      c_[i].canonicalize();
    }
    std::sort(c_, c_ + n_, [&F](const corner& x, const corner& y) {
      return corner::compare(x, y, F);
    });
    coalesce(F);
  }

 private:
  void coalesce(const Field& F) {
    // Coalesce duplicates.
    // The (rd,wr)=(0,0) iteration executes the else{} branch and
    // continues with (1,1), so we start at (1,1) and avoid the
    // special case for wr-1 at wr=0.
    index_t wr = 1;
    for (index_t rd = 1; rd < n_; ++rd) {
      if (c_[rd].eqndx(c_[wr - 1])) {
        F.add(c_[wr - 1].v, c_[rd].v);
      } else {
        c_[wr] = c_[rd];
        wr++;
      }
    }
    n_ = wr;
  }
};
}  // namespace proofs

#endif  // PRIVACY_PROOFS_ZK_LIB_SUMCHECK_QUAD_H_
