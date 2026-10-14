// lygemm -- tile level primitives.
//
// The library is organised the way tile DSLs are: a GEMM is a grid of output
// tiles, every tile owns a register fragment, and the K loop feeds panels of A
// and B into that fragment.
//
//   fragment<Ops, BM, BN>   the BM x BN output tile in registers
//   mma<BK>(A, rsa, B, rsb, c)   c += A_panel * B_panel, K panel of BK
//   tile_gemm<...>::run(...)     the tile grid over a compile time shape
//
// Shapes are compile time constants, so the grid, the panel offsets and the
// register blocking are all resolved by the compiler; nothing is staged or
// predicated, edges included.
//
// SPDX-License-Identifier: MIT
#ifndef LYGEMM_TILE_HPP
#define LYGEMM_TILE_HPP

#include <cstddef>

#include "lygemm/config.hpp"
#include "lygemm/detail/arch.hpp"
#include "lygemm/detail/meta.hpp"

namespace lygemm {
namespace detail {

/* K is unrolled while it is short enough, looped afterwards. */
template <std::size_t K, class Body>
LYGEMM_ALWAYS_INLINE void k_loop(Body&& body) noexcept {
  if constexpr (K <= LYGEMM_MAX_UNROLL_K) {
    meta::static_for<K>([&](auto k) { body(k.value); });
  } else {
    for (std::size_t k = 0; k < K; ++k) body(k);
  }
}

/* Default tile: 4 rows against 3 column registers, the classic x86 register
 * blocking (12 accumulators).  bn is a whole number of vectors so the blocked
 * driver below never needs a masked store. */
template <class Ops>
struct default_tile {
  static constexpr std::size_t bm = 4;
  static constexpr std::size_t bn = 3 * Ops::lanes;
  static constexpr std::size_t bk = 2 * Ops::lanes;
};

/* Tile for a compile time shape: three column registers wide, and as many rows
 * as the register file allows.
 *
 * Every extra row tile re-reads the B panel and repeats the fragment epilogue,
 * so a shape that does not divide the tile pays for it twice.  Covering M in
 * one tile where the registers allow is worth 33% on 6x6x6, 36% on 5x5x5 and
 * 32% on 8x1x8 against stopping at four rows. */
template <class Ops, std::size_t M, std::size_t N>
struct shape_tile {
  static constexpr std::size_t lanes = Ops::lanes;
  static constexpr std::size_t nv =
      meta::clamp<(N + lanes - 1) / lanes, 1, 3>::value;
  static constexpr std::size_t bm =
      meta::clamp<M, 1, meta::max<1, LYGEMM_MAX_ACC / nv>::value>::value;
  static constexpr std::size_t bn = 3 * lanes;
};

/* The output tile: BM x BN accumulators held in registers.
 *
 * Columns are covered by whole column registers; when BN is not a multiple of
 * the vector width the last register of a row is accessed with a mask, so the
 * tail costs one masked load and one masked store instead of a scalar loop.
 * A scalar tail looks cheaper but the vectoriser rewrites it into a chain of
 * shuffles that costs several times what it saves. */
template <class Ops, std::size_t BM, std::size_t BN>
struct fragment {
  static_assert(BM > 0 && BN > 0, "lygemm: empty tile");
  using T = typename Ops::element;
  using V = typename Ops::vec;
  static constexpr std::size_t lanes = Ops::lanes;
  static constexpr std::size_t NV = (BN + lanes - 1) / lanes;
  static constexpr std::size_t TAIL = BN - (NV - 1) * lanes;

  V reg[BM][NV];

  LYGEMM_ALWAYS_INLINE void clear() noexcept {
    meta::static_for<BM>([&](auto i) {
      meta::static_for<NV>([&](auto v) { reg[i.value][v.value] = Ops::zero(); });
    });
  }

  /* Column register VI: whole vectors load and store plainly, the last one is
   * masked when the tile width is not a multiple of the vector width. */
  /* The mask is built once by the caller, not per column and not per k. */
  static LYGEMM_ALWAYS_INLINE typename Ops::mask_type tail_mask() noexcept {
    return Ops::mask_for(TAIL);
  }

  template <std::size_t VI>
  static LYGEMM_ALWAYS_INLINE V load_col(const T* p,
                                         typename Ops::mask_type m) noexcept {
    if constexpr (VI + 1 == NV && TAIL != lanes)
      return Ops::load_masked(p, m);
    else
      return Ops::load(p);
  }

  template <std::size_t VI>
  static LYGEMM_ALWAYS_INLINE void store_col(T* p, V v,
                                             typename Ops::mask_type m) noexcept {
    if constexpr (VI + 1 == NV && TAIL != lanes)
      Ops::store_masked(p, v, m);
    else
      Ops::store(p, v);
  }

  /* C = alpha * fragment + beta * C */
  LYGEMM_ALWAYS_INLINE void store(T* LYGEMM_RESTRICT C, std::ptrdiff_t rc,
                                  T alpha, T beta) const noexcept {
    const V va = Ops::set1(alpha);
    const V vb = Ops::set1(beta);
    const typename Ops::mask_type m = tail_mask();
    meta::static_for<BM>([&](auto i) {
      T* row = C + i.value * rc;
      meta::static_for<NV>([&](auto v) {
        constexpr std::size_t vi = decltype(v)::value;
        T* p = row + vi * lanes;
        V r = Ops::mul(va, reg[i.value][vi]);
        if (beta != T(0)) r = Ops::fma(vb, load_col<vi>(p, m), r);
        store_col<vi>(p, r, m);
      });
    });
  }
};

/* c += A[BM][BK] * B[BK][BN]. */
template <std::size_t BK, class Ops, std::size_t BM, std::size_t BN>
LYGEMM_ALWAYS_INLINE void mma(const typename Ops::element* LYGEMM_RESTRICT A,
                              std::ptrdiff_t rsa,
                              const typename Ops::element* LYGEMM_RESTRICT B,
                              std::ptrdiff_t rsb,
                              fragment<Ops, BM, BN>& c) noexcept {
  using F = fragment<Ops, BM, BN>;
  using V = typename Ops::vec;
  constexpr std::size_t NV = F::NV;

  const typename Ops::mask_type m = F::tail_mask();
  k_loop<BK>([&](std::size_t k) {
    const typename Ops::element* brow = B + k * rsb;
    V bv[NV];
    meta::static_for<NV>([&](auto v) {
      constexpr std::size_t vi = decltype(v)::value;
      bv[vi] = F::template load_col<vi>(brow + vi * Ops::lanes, m);
    });
    meta::static_for<BM>([&](auto i) {
      const V a = Ops::set1(A[i.value * rsa + k]);
      meta::static_for<NV>([&](auto v) {
        c.reg[i.value][v.value] =
            Ops::fma(a, bv[v.value], c.reg[i.value][v.value]);
      });
    });
  });
}

/* The tile grid.
 *
 * The grid itself is a runtime loop and only the tile shape is a template
 * argument.  Unrolling a grid of tiles instead multiplies the code size by the
 * number of tiles, and past a few kilobytes of straight line code the front
 * end, not the FMA units, sets the pace: 16x16x16 measures 4x faster with a
 * runtime grid than with the fully unrolled one.
 *
 * K is folded in panels of BK; the last panel is unrolled to the exact
 * remainder, so no access ever runs past the end of a panel. */
template <class Ops, std::size_t M, std::size_t N, std::size_t K,
          std::size_t BM = shape_tile<Ops, M, N>::bm,
          std::size_t BN = shape_tile<Ops, M, N>::bn,
          std::size_t BK = default_tile<Ops>::bk>
struct tile_gemm {
  using T = typename Ops::element;

  static_assert(M > 0 && N > 0 && K > 0,
                "lygemm: zero sized shapes are not supported");
  static_assert(BK <= LYGEMM_MAX_UNROLL_K,
                "lygemm: LYGEMM_MAX_UNROLL_K must cover a K panel");

  static constexpr std::size_t m_full = M / BM;
  static constexpr std::size_t n_full = N / BN;
  static constexpr std::size_t m_tail = M % BM;
  static constexpr std::size_t n_tail = N % BN;

  /* One BM x BN tile: clear, fold every K panel in, write back. */
  template <std::size_t TM, std::size_t TN>
  static LYGEMM_ALWAYS_INLINE void body(const T* LYGEMM_RESTRICT A,
                                        std::ptrdiff_t ra,
                                        const T* LYGEMM_RESTRICT B,
                                        std::ptrdiff_t rb,
                                        T* LYGEMM_RESTRICT C, std::ptrdiff_t rc,
                                        T alpha, T beta) noexcept {
    fragment<Ops, TM, TN> acc;
    acc.clear();
    std::size_t k = 0;
    for (; k + BK <= K; k += BK) mma<BK>(A + k, ra, B + k * rb, rb, acc);
    if constexpr (K % BK != 0)
      mma<K % BK>(A + k, ra, B + k * rb, rb, acc);
    acc.store(C, rc, alpha, beta);
  }

  static LYGEMM_ALWAYS_INLINE void run(const T* LYGEMM_RESTRICT A,
                                       std::ptrdiff_t ra,
                                       const T* LYGEMM_RESTRICT B,
                                       std::ptrdiff_t rb,
                                       T* LYGEMM_RESTRICT C, std::ptrdiff_t rc,
                                       T alpha, T beta) noexcept {
    constexpr std::size_t rows = m_full * BM;
    constexpr std::size_t cols = n_full * BN;
    for (std::size_t i = 0; i < rows; i += BM) {
      for (std::size_t j = 0; j < cols; j += BN)
        body<BM, BN>(A + i * ra, ra, B + j, rb, C + i * rc + j, rc, alpha,
                     beta);
      if constexpr (n_tail)
        body<BM, n_tail>(A + i * ra, ra, B + cols, rb, C + i * rc + cols, rc,
                         alpha, beta);
    }
    if constexpr (m_tail) {
      for (std::size_t j = 0; j < cols; j += BN)
        body<m_tail, BN>(A + rows * ra, ra, B + j, rb, C + rows * rc + j, rc,
                         alpha, beta);
      if constexpr (n_tail)
        body<m_tail, n_tail>(A + rows * ra, ra, B + cols, rb,
                             C + rows * rc + cols, rc, alpha, beta);
    }
  }
};

/* Grid of tiles for a batch that shares one operand: the shared operand keeps
 * tile invariant offsets, the batch index is the innermost loop.  ShareA picks
 * which of the two operands is the reused one; the other carries the batch
 * stride. */
template <class Ops, std::size_t M, std::size_t N, std::size_t K, bool ShareA,
          std::size_t BM = shape_tile<Ops, M, N>::bm,
          std::size_t BN = shape_tile<Ops, M, N>::bn,
          std::size_t BK = default_tile<Ops>::bk>
struct shared_batch_gemm {
  using T = typename Ops::element;
  using G = tile_gemm<Ops, M, N, K, BM, BN, BK>;

  static LYGEMM_ALWAYS_INLINE void run(std::size_t batch,
                                       const T* LYGEMM_RESTRICT A,
                                       std::ptrdiff_t sa,
                                       const T* LYGEMM_RESTRICT B,
                                       std::ptrdiff_t sb, T* LYGEMM_RESTRICT C,
                                       std::ptrdiff_t sc, T alpha,
                                       T beta) noexcept {
    constexpr std::size_t rows = G::m_full * BM;
    constexpr std::size_t cols = G::n_full * BN;
    const std::ptrdiff_t oa = ShareA ? 0 : sa; /* batch offset of A */
    const std::ptrdiff_t ob = ShareA ? sb : 0; /* batch offset of B */
    for (std::size_t i = 0; i < rows; i += BM) {
      for (std::size_t j = 0; j < cols; j += BN)
        for (std::size_t b = 0; b < batch; ++b)
          G::template body<BM, BN>(A + b * oa + i * K, K, B + b * ob + j, N,
                                   C + b * sc + i * N + j, N, alpha, beta);
      if constexpr (G::n_tail)
        for (std::size_t b = 0; b < batch; ++b)
          G::template body<BM, G::n_tail>(A + b * oa + i * K, K,
                                          B + b * ob + cols, N,
                                          C + b * sc + i * N + cols, N, alpha,
                                          beta);
    }
    if constexpr (G::m_tail) {
      for (std::size_t j = 0; j < cols; j += BN)
        for (std::size_t b = 0; b < batch; ++b)
          G::template body<G::m_tail, BN>(A + b * oa + rows * K, K,
                                          B + b * ob + j, N,
                                          C + b * sc + rows * N + j, N, alpha,
                                          beta);
      if constexpr (G::n_tail)
        for (std::size_t b = 0; b < batch; ++b)
          G::template body<G::m_tail, G::n_tail>(
              A + b * oa + rows * K, K, B + b * ob + cols, N,
              C + b * sc + rows * N + cols, N, alpha, beta);
    }
  }
};

} // namespace detail
} // namespace lygemm

#endif /* LYGEMM_TILE_HPP */
