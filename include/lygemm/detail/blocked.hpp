// lygemm -- blocked driver for runtime shapes, plus dispatch and BLAS hook.
//
// Same tile machinery as the compile time grid, with a runtime grid: full
// panels are consumed straight from memory, partial panels are staged through
// zero padded scratch buffers so that no access has to be masked or guarded.
// The fragment is written back through a scratch tile, which keeps the column
// tail inside the tile instead of in the caller's matrix.
//
// SPDX-License-Identifier: MIT
#ifndef LYGEMM_BLOCKED_HPP
#define LYGEMM_BLOCKED_HPP

#include <algorithm>
#include <cstddef>

#include "lygemm/config.hpp"
#include "lygemm/detail/arch.hpp"
#include "lygemm/detail/meta.hpp"
#include "lygemm/detail/tile.hpp"

/* Shapes answered by the compile time tile grid through the runtime
 * interface.  Cost: one instantiation per entry and scalar type. */
#ifndef LYGEMM_DISPATCH_LIST
#define LYGEMM_DISPATCH_LIST(X) \
  X(1, 1, 1)                    \
  X(2, 2, 2)                    \
  X(3, 3, 3)                    \
  X(4, 4, 4)                    \
  X(5, 5, 5)                    \
  X(6, 6, 6)                    \
  X(7, 7, 7)                    \
  X(8, 8, 8)                    \
  X(12, 12, 12)                 \
  X(16, 16, 16)
#endif

namespace lygemm {
namespace detail {

/* Zero padded copy of a runtime sized block into a R x C scratch panel. */
template <class T, std::size_t R, std::size_t C>
inline void stage_panel(const T* LYGEMM_RESTRICT src, std::ptrdiff_t rs,
                        T* LYGEMM_RESTRICT dst, std::size_t rows,
                        std::size_t cols) noexcept {
  for (std::size_t i = 0; i < rows; ++i) {
    const T* s = src + i * rs;
    T* d = dst + i * C;
    std::size_t j = 0;
    for (; j < cols; ++j) d[j] = s[j];
    for (; j < C; ++j) d[j] = T(0);
  }
  for (std::size_t i = rows; i < R; ++i) {
    T* d = dst + i * C;
    for (std::size_t j = 0; j < C; ++j) d[j] = T(0);
  }
}

template <class Ops, class T>
inline void gemm_blocked(std::size_t m, std::size_t n, std::size_t k,
                         const T* LYGEMM_RESTRICT A, std::ptrdiff_t ra,
                         const T* LYGEMM_RESTRICT B, std::ptrdiff_t rb,
                         T* LYGEMM_RESTRICT C, std::ptrdiff_t rc, T alpha,
                         T beta) noexcept {
  constexpr std::size_t BM = default_tile<Ops>::bm;
  constexpr std::size_t BN = default_tile<Ops>::bn;
  constexpr std::size_t BK = default_tile<Ops>::bk;

  T a_panel[BM * BK];
  T b_panel[BK * BN];
  T c_tile[BM * BN];

  for (std::size_t i0 = 0; i0 < m; i0 += BM) {
    const std::size_t bm = std::min(BM, m - i0);
    for (std::size_t j0 = 0; j0 < n; j0 += BN) {
      const std::size_t bn = std::min(BN, n - j0);
      fragment<Ops, BM, BN> acc;
      acc.clear();
      for (std::size_t k0 = 0; k0 < k; k0 += BK) {
        const std::size_t bk = std::min(BK, k - k0);
        const T* ap = A + i0 * ra + k0;
        std::ptrdiff_t sa = ra;
        if (bm != BM || bk != BK) {
          stage_panel<T, BM, BK>(ap, ra, a_panel, bm, bk);
          ap = a_panel;
          sa = BK;
        }
        const T* bp = B + k0 * rb + j0;
        std::ptrdiff_t sb = rb;
        if (bn != BN || bk != BK) {
          stage_panel<T, BK, BN>(bp, rb, b_panel, bk, bn);
          bp = b_panel;
          sb = BN;
        }
        mma<BK>(ap, sa, bp, sb, acc);
      }
      acc.store(c_tile, BN, alpha, T(0));
      for (std::size_t r = 0; r < bm; ++r) {
        const T* s = c_tile + r * BN;
        T* p = C + (i0 + r) * rc + j0;
        if (beta != T(0))
          for (std::size_t c = 0; c < bn; ++c) p[c] = s[c] + beta * p[c];
        else
          for (std::size_t c = 0; c < bn; ++c) p[c] = s[c];
      }
    }
  }
}

/* Runtime -> compile time dispatch.  True when the tile grid ran. */
template <class T>
LYGEMM_ALWAYS_INLINE bool gemm_static_if(std::size_t m, std::size_t n,
                                         std::size_t k,
                                         const T* LYGEMM_RESTRICT A,
                                         const T* LYGEMM_RESTRICT B,
                                         T* LYGEMM_RESTRICT C, T alpha,
                                         T beta) noexcept {
#define LYGEMM_DISPATCH_CASE(M, N, K)                                    \
  if (m == (M) && n == (N) && k == (K)) {                                \
    tile_gemm<typename shape_ops<T, N>::type, M, N, K>::run(A, K, B, N,  \
                                                             C, N, alpha, \
                                                             beta);       \
    return true;                                                         \
  }
  LYGEMM_DISPATCH_LIST(LYGEMM_DISPATCH_CASE)
#undef LYGEMM_DISPATCH_CASE
  return false;
}

#if LYGEMM_USE_BLAS
template <class T>
struct blas_gemm;

template <>
struct blas_gemm<double> {
  static inline void run(std::size_t m, std::size_t n, std::size_t k,
                         double alpha, const double* A, const double* B,
                         double beta, double* C) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, static_cast<int>(m),
                static_cast<int>(n), static_cast<int>(k), alpha, A,
                static_cast<int>(k), B, static_cast<int>(n), beta, C,
                static_cast<int>(n));
  }
};

template <>
struct blas_gemm<float> {
  static inline void run(std::size_t m, std::size_t n, std::size_t k,
                         float alpha, const float* A, const float* B,
                         float beta, float* C) {
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, static_cast<int>(m),
                static_cast<int>(n), static_cast<int>(k), alpha, A,
                static_cast<int>(k), B, static_cast<int>(n), beta, C,
                static_cast<int>(n));
  }
};
#endif /* LYGEMM_USE_BLAS */

} // namespace detail
} // namespace lygemm

#endif /* LYGEMM_BLOCKED_HPP */
