// lygemm -- header-only small matrix multiplication.
//
//   #include <lygemm.hpp>
//   double A[16] = {...}, B[16] = {...}, C[16];
//   lygemm::gemm<4, 4, 4>(A, B, C);          // compile time shape, tiled + unrolled
//   lygemm::gemm(m, n, k, A, B, C);          // runtime shape
//   lygemm::gemm_batch<8, 8, 8>(n, A, B, C); // n independent products
//
// Row major everywhere.  A tile is (rows x columns) of C; the K loop feeds
// panels of A and B into the tile fragment.
//
//   lygemm::tile::fragment<T, BM, BN>   register tile
//   lygemm::tile::mma<BK>(A, rsa, B, rsb, frag)   frag += A_panel * B_panel
//   lygemm::tile::gemm<M, N, K, T>(...)  run the tile grid for a shape
//
// C must not overlap A or B.  Shapes must be non zero.
//
// SPDX-License-Identifier: MIT
#ifndef LYGEMM_HPP
#define LYGEMM_HPP

#include <cstddef>

#include "lygemm/config.hpp"
#include "lygemm/detail/arch.hpp"
#include "lygemm/detail/blocked.hpp"
#include "lygemm/detail/jit.hpp"
#include "lygemm/detail/meta.hpp"
#include "lygemm/detail/tile.hpp"

namespace lygemm {

/* =========================================================================
 * compile time shapes
 * ====================================================================== */

/* C = alpha * A * B + beta * C.  beta == 0 leaves C untouched before the
 * store, so C does not have to be initialised. */
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
LYGEMM_ALWAYS_INLINE void gemm(T alpha, const T* LYGEMM_RESTRICT A,
                               const T* LYGEMM_RESTRICT B, T beta,
                               T* LYGEMM_RESTRICT C) noexcept {
  detail::tile_gemm<typename detail::shape_ops<T, N>::type, M, N, K>::run(
      A, K, B, N, C, N, alpha, beta);
}

/* C = A * B. */
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
LYGEMM_ALWAYS_INLINE void gemm(const T* LYGEMM_RESTRICT A,
                               const T* LYGEMM_RESTRICT B,
                               T* LYGEMM_RESTRICT C) noexcept {
  gemm<M, N, K, T>(T(1), A, B, T(0), C);
}

/* =========================================================================
 * runtime shapes
 * ====================================================================== */

/* Row major C[m][n] = alpha * A[m][k] * B[k][n] + beta * C[m][n].
 * Shapes in LYGEMM_DISPATCH_LIST run the compile time tile grid, everything
 * else the blocked driver (or BLAS, when built with LYGEMM_USE_BLAS and
 * m*n*k >= LYGEMM_BLAS_THRESHOLD). */
/* Everything the generated kernel cannot handle, kept out of the hot path so
 * that a call through it stays a handful of instructions. */
template <class T>
LYGEMM_NOINLINE void gemm_fallback(std::size_t m, std::size_t n, std::size_t k,
                                   const T* LYGEMM_RESTRICT A,
                                   const T* LYGEMM_RESTRICT B,
                                   T* LYGEMM_RESTRICT C, T alpha,
                                   T beta) noexcept {
  if (detail::gemm_static_if<T>(m, n, k, A, B, C, alpha, beta)) return;
#if LYGEMM_USE_BLAS
  if (m * n * k >= static_cast<std::size_t>(LYGEMM_BLAS_THRESHOLD)) {
    detail::blas_gemm<T>::run(m, n, k, alpha, A, B, beta, C);
    return;
  }
#endif
  using narrow = typename detail::narrow_ops<T>::type;
  if (n <= narrow::lanes)
    detail::gemm_blocked<narrow, T>(m, n, k, A, k, B, n, C, n, alpha, beta);
  else
    detail::gemm_blocked<typename detail::default_ops<T>::type, T>(
        m, n, k, A, k, B, n, C, n, alpha, beta);
}

/* The generated kernel comes first: it skips the dispatch chain and the blocked
 * driver, and it is emitted at the build's own vector width (512 bit when
 * compiled for AVX-512, with the A operand broadcast straight out of memory),
 * so it beats the compile time tile grid on the shapes it covers. */
template <class T>
inline void gemm(std::size_t m, std::size_t n, std::size_t k,
                 const T* LYGEMM_RESTRICT A, const T* LYGEMM_RESTRICT B,
                 T* LYGEMM_RESTRICT C, T alpha = T(1),
                 T beta = T(0)) noexcept {
  if (alpha == T(1) && beta == T(0)) {
    const auto fn = detail::jit::entry<T>::get(m, n, k);
    if (fn != nullptr) {
      fn(A, B, C);
      return;
    }
  }
  gemm_fallback<T>(m, n, k, A, B, C, alpha, beta);
}

/* =========================================================================
 * batches of compile time shapes
 *
 * A batch element is M*K / K*N / M*N elements apart unless a stride is given.
 * A stride of 0 broadcasts one operand across the batch.
 * ====================================================================== */

template <std::size_t M, std::size_t N, std::size_t K, class T = double>
inline void gemm_batch(std::size_t batch, const T* LYGEMM_RESTRICT A,
                       std::ptrdiff_t sa, const T* LYGEMM_RESTRICT B,
                       std::ptrdiff_t sb, T* LYGEMM_RESTRICT C,
                       std::ptrdiff_t sc, T alpha = T(1),
                       T beta = T(0)) noexcept {
  using Ops = typename detail::shape_ops<T, N>::type;
  for (std::size_t i = 0; i < batch; ++i)
    detail::tile_gemm<Ops, M, N, K>::run(A + i * sa, K, B + i * sb, N,
                                         C + i * sc, N, alpha, beta);
}

template <std::size_t M, std::size_t N, std::size_t K, class T = double>
inline void gemm_batch(std::size_t batch, const T* LYGEMM_RESTRICT A,
                       const T* LYGEMM_RESTRICT B, T* LYGEMM_RESTRICT C,
                       T alpha = T(1), T beta = T(0)) noexcept {
  gemm_batch<M, N, K, T>(batch, A, static_cast<std::ptrdiff_t>(M * K), B,
                         static_cast<std::ptrdiff_t>(K * N), C,
                         static_cast<std::ptrdiff_t>(M * N), alpha, beta);
}

/* =========================================================================
 * batches of runtime shapes
 * ====================================================================== */

template <class T>
inline void gemm_batch(std::size_t batch, std::size_t m, std::size_t n,
                       std::size_t k, const T* LYGEMM_RESTRICT A,
                       std::ptrdiff_t sa, const T* LYGEMM_RESTRICT B,
                       std::ptrdiff_t sb, T* LYGEMM_RESTRICT C,
                       std::ptrdiff_t sc, T alpha = T(1),
                       T beta = T(0)) noexcept {
  for (std::size_t i = 0; i < batch; ++i)
    gemm<T>(m, n, k, A + i * sa, B + i * sb, C + i * sc, alpha, beta);
}

template <class T>
inline void gemm_batch(std::size_t batch, std::size_t m, std::size_t n,
                       std::size_t k, const T* LYGEMM_RESTRICT A,
                       const T* LYGEMM_RESTRICT B, T* LYGEMM_RESTRICT C,
                       T alpha = T(1), T beta = T(0)) noexcept {
  gemm_batch<T>(batch, m, n, k, A, static_cast<std::ptrdiff_t>(m * k), B,
                static_cast<std::ptrdiff_t>(k * n), C,
                static_cast<std::ptrdiff_t>(m * n), alpha, beta);
}

/* =========================================================================
 * the generated kernel for a shape
 *
 * Fetch it once and call it in a loop: a fixed shape used repeatedly is then
 * as fast as a hand written kernel, with none of the per call dispatch and
 * lookup that lygemm::gemm(m, n, k, ...) pays.  Returns null when the shapes
 * are outside the generator's range (see the README) or the JIT is off.
 * ====================================================================== */

template <class T>
using kernel_fn = void (*)(const T*, const T*, T*);

template <class T = double>
inline kernel_fn<T> kernel(std::size_t m, std::size_t n,
                           std::size_t k) noexcept {
  return detail::jit::entry<T>::get(m, n, k);
}

/* =========================================================================
 * batches with one shared operand
 *
 *   C[b] = alpha * A * B[b] + beta * C[b]      (A shared)
 *   C[b] = alpha * A[b] * B + beta * C[b]      (B shared)
 *
 * The tile grid stays outside the batch loop here, so the tile offsets of the
 * shared operand are computed once for the whole batch instead of per element.
 * ====================================================================== */

template <std::size_t M, std::size_t N, std::size_t K, class T = double>
inline void gemm_batch_shared_a(std::size_t batch, const T* LYGEMM_RESTRICT A,
                                const T* LYGEMM_RESTRICT B, std::ptrdiff_t sb,
                                T* LYGEMM_RESTRICT C, std::ptrdiff_t sc,
                                T alpha = T(1), T beta = T(0)) noexcept {
  detail::shared_batch_gemm<typename detail::shape_ops<T, N>::type, M, N, K,
                            true>::run(batch, A, 0, B, sb, C, sc, alpha, beta);
}

template <std::size_t M, std::size_t N, std::size_t K, class T = double>
inline void gemm_batch_shared_b(std::size_t batch, const T* LYGEMM_RESTRICT A,
                                std::ptrdiff_t sa, const T* LYGEMM_RESTRICT B,
                                T* LYGEMM_RESTRICT C, std::ptrdiff_t sc,
                                T alpha = T(1), T beta = T(0)) noexcept {
  detail::shared_batch_gemm<typename detail::shape_ops<T, N>::type, M, N, K,
                            false>::run(batch, A, sa, B, 0, C, sc, alpha, beta);
}

/* =========================================================================
 * tile level API: compose your own kernel
 * ====================================================================== */

namespace tile {

/* Register tile of the output, BM rows x BN columns. */
template <class T, std::size_t BM, std::size_t BN>
using fragment = detail::fragment<typename detail::shape_ops<T, BN>::type, BM, BN>;

/* frag += A[BM][BK] * B[BK][BN]. */
template <std::size_t BK, class T, std::size_t BM, std::size_t BN>
LYGEMM_ALWAYS_INLINE void mma(const T* LYGEMM_RESTRICT A, std::ptrdiff_t rsa,
                              const T* LYGEMM_RESTRICT B, std::ptrdiff_t rsb,
                              fragment<T, BM, BN>& frag) noexcept {
  detail::mma<BK>(A, rsa, B, rsb, frag);
}

/* The tile grid of a compile time shape, with default tile sizes. */
template <std::size_t M, std::size_t N, std::size_t K, class T = double>
LYGEMM_ALWAYS_INLINE void gemm(const T* LYGEMM_RESTRICT A,
                               const T* LYGEMM_RESTRICT B,
                               T* LYGEMM_RESTRICT C) noexcept {
  detail::tile_gemm<typename detail::shape_ops<T, N>::type, M, N, K>::run(
      A, K, B, N, C, N, T(1), T(0));
}

} // namespace tile

} // namespace lygemm

#endif /* LYGEMM_HPP */
