// lygemm -- vector backend abstraction.
//
// One micro kernel is written against these members, so a backend is a thin
// wrapper over one instruction set:
//
//   element, vec, lanes
//   zero / load / store / set1 / mul / fma
//   fma(a, b, c) == a * b + c
//
// SPDX-License-Identifier: MIT
#ifndef LYGEMM_ARCH_HPP
#define LYGEMM_ARCH_HPP

#include <cstddef>
#include <type_traits>

#include "lygemm/config.hpp"

namespace lygemm {
namespace detail {

/* ------------------------------------------------------- portable scalar */
template <class T>
struct scalar_ops {
  using element = T;
  using vec = T;
  using mask_type = std::size_t; /* number of live lanes */
  static constexpr std::size_t lanes = 1;

  static LYGEMM_ALWAYS_INLINE vec zero() noexcept { return T(0); }
  static LYGEMM_ALWAYS_INLINE vec load(const T* p) noexcept { return *p; }
  static LYGEMM_ALWAYS_INLINE void store(T* p, vec v) noexcept { *p = v; }
  static LYGEMM_ALWAYS_INLINE mask_type mask_for(std::size_t c) noexcept {
    return c;
  }
  static LYGEMM_ALWAYS_INLINE vec load_masked(const T* p,
                                              mask_type m) noexcept {
    return m ? *p : T(0);
  }
  static LYGEMM_ALWAYS_INLINE void store_masked(T* p, vec v,
                                                mask_type m) noexcept {
    if (m) *p = v;
  }
  static LYGEMM_ALWAYS_INLINE vec set1(T a) noexcept { return a; }
  static LYGEMM_ALWAYS_INLINE vec mul(vec a, vec b) noexcept { return a * b; }
  static LYGEMM_ALWAYS_INLINE vec fma(vec a, vec b, vec c) noexcept {
    return a * b + c;
  }
};

#if LYGEMM_ARCH >= 2
/* ------------------------------------------------------------------ AVX2 */
struct avx2_ops_f64 {
  using element = double;
  using vec = __m256d;
  using mask_type = __m256i;
  static constexpr std::size_t lanes = 4;

  static LYGEMM_ALWAYS_INLINE vec zero() noexcept { return _mm256_setzero_pd(); }
  static LYGEMM_ALWAYS_INLINE vec load(const double* p) noexcept {
    return _mm256_loadu_pd(p);
  }
  static LYGEMM_ALWAYS_INLINE void store(double* p, vec v) noexcept {
    _mm256_storeu_pd(p, v);
  }
  static LYGEMM_ALWAYS_INLINE mask_type mask_for(std::size_t c) noexcept {
    const __m256i idx = _mm256_set_epi64x(3, 2, 1, 0);
    const __m256i cnt = _mm256_set1_epi64x(static_cast<long long>(c));
    return _mm256_cmpgt_epi64(cnt, idx);
  }
  static LYGEMM_ALWAYS_INLINE vec load_masked(const double* p,
                                              mask_type m) noexcept {
    return _mm256_maskload_pd(p, m);
  }
  static LYGEMM_ALWAYS_INLINE void store_masked(double* p, vec v,
                                                mask_type m) noexcept {
    _mm256_maskstore_pd(p, m, v);
  }
  static LYGEMM_ALWAYS_INLINE vec set1(double a) noexcept {
    return _mm256_set1_pd(a);
  }
  static LYGEMM_ALWAYS_INLINE vec mul(vec a, vec b) noexcept {
    return _mm256_mul_pd(a, b);
  }
  static LYGEMM_ALWAYS_INLINE vec fma(vec a, vec b, vec c) noexcept {
#if LYGEMM_HAS_FMA
    return _mm256_fmadd_pd(a, b, c);
#else
    return _mm256_add_pd(_mm256_mul_pd(a, b), c);
#endif
  }
};

struct avx2_ops_f32 {
  using element = float;
  using vec = __m256;
  using mask_type = __m256i;
  static constexpr std::size_t lanes = 8;

  static LYGEMM_ALWAYS_INLINE vec zero() noexcept { return _mm256_setzero_ps(); }
  static LYGEMM_ALWAYS_INLINE vec load(const float* p) noexcept {
    return _mm256_loadu_ps(p);
  }
  static LYGEMM_ALWAYS_INLINE void store(float* p, vec v) noexcept {
    _mm256_storeu_ps(p, v);
  }
  static LYGEMM_ALWAYS_INLINE mask_type mask_for(std::size_t c) noexcept {
    const __m256i idx = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    const __m256i cnt = _mm256_set1_epi32(static_cast<int>(c));
    return _mm256_cmpgt_epi32(cnt, idx);
  }
  static LYGEMM_ALWAYS_INLINE vec load_masked(const float* p,
                                              mask_type m) noexcept {
    return _mm256_maskload_ps(p, m);
  }
  static LYGEMM_ALWAYS_INLINE void store_masked(float* p, vec v,
                                                mask_type m) noexcept {
    _mm256_maskstore_ps(p, m, v);
  }
  static LYGEMM_ALWAYS_INLINE vec set1(float a) noexcept {
    return _mm256_set1_ps(a);
  }
  static LYGEMM_ALWAYS_INLINE vec mul(vec a, vec b) noexcept {
    return _mm256_mul_ps(a, b);
  }
  static LYGEMM_ALWAYS_INLINE vec fma(vec a, vec b, vec c) noexcept {
#if LYGEMM_HAS_FMA
    return _mm256_fmadd_ps(a, b, c);
#else
    return _mm256_add_ps(_mm256_mul_ps(a, b), c);
#endif
  }
};
#endif /* LYGEMM_ARCH >= 2 */

#if LYGEMM_ARCH >= 3
/* --------------------------------------------------------------- AVX-512 */
struct avx512_ops_f64 {
  using element = double;
  using vec = __m512d;
  using mask_type = __mmask8;
  static constexpr std::size_t lanes = 8;

  static LYGEMM_ALWAYS_INLINE vec zero() noexcept { return _mm512_setzero_pd(); }
  static LYGEMM_ALWAYS_INLINE vec load(const double* p) noexcept {
    return _mm512_loadu_pd(p);
  }
  static LYGEMM_ALWAYS_INLINE void store(double* p, vec v) noexcept {
    _mm512_storeu_pd(p, v);
  }
  static LYGEMM_ALWAYS_INLINE mask_type mask_for(std::size_t c) noexcept {
    return static_cast<__mmask8>((1u << c) - 1u);
  }
  static LYGEMM_ALWAYS_INLINE vec load_masked(const double* p,
                                              mask_type m) noexcept {
    return _mm512_maskz_loadu_pd(m, p);
  }
  static LYGEMM_ALWAYS_INLINE void store_masked(double* p, vec v,
                                                mask_type m) noexcept {
    _mm512_mask_storeu_pd(p, m, v);
  }
  static LYGEMM_ALWAYS_INLINE vec set1(double a) noexcept {
    return _mm512_set1_pd(a);
  }
  static LYGEMM_ALWAYS_INLINE vec mul(vec a, vec b) noexcept {
    return _mm512_mul_pd(a, b);
  }
  static LYGEMM_ALWAYS_INLINE vec fma(vec a, vec b, vec c) noexcept {
    return _mm512_fmadd_pd(a, b, c);
  }
};

struct avx512_ops_f32 {
  using element = float;
  using vec = __m512;
  using mask_type = __mmask16;
  static constexpr std::size_t lanes = 16;

  static LYGEMM_ALWAYS_INLINE vec zero() noexcept { return _mm512_setzero_ps(); }
  static LYGEMM_ALWAYS_INLINE vec load(const float* p) noexcept {
    return _mm512_loadu_ps(p);
  }
  static LYGEMM_ALWAYS_INLINE void store(float* p, vec v) noexcept {
    _mm512_storeu_ps(p, v);
  }
  static LYGEMM_ALWAYS_INLINE mask_type mask_for(std::size_t c) noexcept {
    return static_cast<__mmask16>((1u << c) - 1u);
  }
  static LYGEMM_ALWAYS_INLINE vec load_masked(const float* p,
                                              mask_type m) noexcept {
    return _mm512_maskz_loadu_ps(m, p);
  }
  static LYGEMM_ALWAYS_INLINE void store_masked(float* p, vec v,
                                                mask_type m) noexcept {
    _mm512_mask_storeu_ps(p, m, v);
  }
  static LYGEMM_ALWAYS_INLINE vec set1(float a) noexcept {
    return _mm512_set1_ps(a);
  }
  static LYGEMM_ALWAYS_INLINE vec mul(vec a, vec b) noexcept {
    return _mm512_mul_ps(a, b);
  }
  static LYGEMM_ALWAYS_INLINE vec fma(vec a, vec b, vec c) noexcept {
    return _mm512_fmadd_ps(a, b, c);
  }
};
#endif /* LYGEMM_ARCH >= 3 */

template <class T>
struct default_ops {
  using type = scalar_ops<T>;
};
#if LYGEMM_ARCH >= 3
template <>
struct default_ops<double> {
  using type = avx512_ops_f64;
};
template <>
struct default_ops<float> {
  using type = avx512_ops_f32;
};
#elif LYGEMM_ARCH >= 2
template <>
struct default_ops<double> {
  using type = avx2_ops_f64;
};
template <>
struct default_ops<float> {
  using type = avx2_ops_f32;
};
#endif

/* A wide backend is only worth it when the columns fill it.  With eight lanes
 * a four column problem leaves half of every FMA register idle, and the parts
 * that have AVX-512 also drop their clock for 512-bit work: on a Skylake-SP a
 * 4x4x4 runs in 20.5 ns with the 256-bit backend and 43.9 ns with the 512-bit
 * one.  shape_ops picks the widest backend whose vector width the column count
 * actually fills, the way a per-shape code generator would. */
template <class T, std::size_t N>
struct shape_ops {
  using type = typename default_ops<T>::type;
};
#if LYGEMM_ARCH >= 3
template <std::size_t N>
struct shape_ops<double, N> {
  using type = typename std::conditional<(N <= avx2_ops_f64::lanes),
                                         avx2_ops_f64, avx512_ops_f64>::type;
};
template <std::size_t N>
struct shape_ops<float, N> {
  using type = typename std::conditional<(N <= avx2_ops_f32::lanes),
                                         avx2_ops_f32, avx512_ops_f32>::type;
};
#endif

/* The narrow backend, for the cases where the column count is only known at
 * run time. */
template <class T>
struct narrow_ops {
  using type = typename default_ops<T>::type;
};
#if LYGEMM_ARCH >= 3
template <>
struct narrow_ops<double> {
  using type = avx2_ops_f64;
};
template <>
struct narrow_ops<float> {
  using type = avx2_ops_f32;
};
#endif

} // namespace detail
} // namespace lygemm

#endif /* LYGEMM_ARCH_HPP */
