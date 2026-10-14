// lygemm -- header-only small matrix multiplication kernels.
// SPDX-License-Identifier: MIT
#ifndef LYGEMM_CONFIG_HPP
#define LYGEMM_CONFIG_HPP

#include <cstddef>

#if defined(_MSVC_LANG)
#  define LYGEMM_CPLUSPLUS _MSVC_LANG
#else
#  define LYGEMM_CPLUSPLUS __cplusplus
#endif

#if LYGEMM_CPLUSPLUS < 201703L
#  error "lygemm requires C++17 or later (-std=c++17)"
#endif

/* -------------------------------------------------------------------------
 * Instruction set selection.
 *
 *   LYGEMM_ARCH == 0   portable scalar (compiler auto-vectorisation only)
 *   LYGEMM_ARCH == 2   AVX2, 4 doubles / 8 floats, FMA when available
 *   LYGEMM_ARCH == 3   AVX-512F + AVX-512DQ, 8 doubles / 16 floats
 *
 * The default follows the flags the translation unit is compiled with, so
 * `-mavx2 -mfma` (or `-march=native`) is enough to get the wide kernels.
 * ------------------------------------------------------------------------- */
#ifndef LYGEMM_ARCH
#  if defined(__AVX512F__) && defined(__AVX512DQ__)
#    define LYGEMM_ARCH 3
#  elif defined(__AVX2__)
#    define LYGEMM_ARCH 2
#  else
#    define LYGEMM_ARCH 0
#  endif
#endif

#if LYGEMM_ARCH == 2 && !defined(__AVX2__)
#  error "LYGEMM_ARCH=2 requires AVX2: compile with -mavx2"
#elif LYGEMM_ARCH == 3 && !(defined(__AVX512F__) && defined(__AVX512DQ__))
#  error "LYGEMM_ARCH=3 requires AVX-512F/DQ: compile with -mavx512f -mavx512dq"
#elif LYGEMM_ARCH != 0 && LYGEMM_ARCH != 2 && LYGEMM_ARCH != 3
#  error "LYGEMM_ARCH must be 0, 2 or 3"
#endif

#if LYGEMM_ARCH > 0
#  include <immintrin.h>
#endif

#if defined(__FMA__)
#  define LYGEMM_HAS_FMA 1
#else
#  define LYGEMM_HAS_FMA 0
#endif

/* -------------------------------------------------------------------------
 * Tuning knobs.  All of them can be overridden from the command line.
 * ------------------------------------------------------------------------- */
/* Accumulator registers a tile may occupy: on x86 with 16 vector registers a
 * tile of 12 leaves room for the operands and the loop state. */
#ifndef LYGEMM_MAX_ACC
#  define LYGEMM_MAX_ACC 12
#endif
/* K panels are unrolled completely up to this value, then looped. */
#ifndef LYGEMM_MAX_UNROLL_K
#  define LYGEMM_MAX_UNROLL_K 32
#endif

/* -------------------------------------------------------------------------
 * Optional BLAS fallback for the runtime interface.  Header-only by default:
 * define LYGEMM_USE_BLAS=1 and provide <cblas.h> plus a CBLAS library to make
 * the runtime entry point hand large problems over to BLAS.
 * ------------------------------------------------------------------------- */
#ifndef LYGEMM_USE_BLAS
#  define LYGEMM_USE_BLAS 0
#endif
#if LYGEMM_USE_BLAS
#  include <cblas.h>
#  ifndef LYGEMM_BLAS_THRESHOLD
#    define LYGEMM_BLAS_THRESHOLD 32768 /* m*n*k above which BLAS is used */
#  endif
#endif

/* ------------------------------------------------------------------------- */
#if defined(__GNUC__) || defined(__clang__)
#  define LYGEMM_ALWAYS_INLINE inline __attribute__((always_inline))
#  define LYGEMM_RESTRICT __restrict__
#  define LYGEMM_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#  define LYGEMM_ALWAYS_INLINE __forceinline
#  define LYGEMM_RESTRICT __restrict
#  define LYGEMM_NOINLINE __declspec(noinline)
#else
#  define LYGEMM_ALWAYS_INLINE inline
#  define LYGEMM_RESTRICT
#  define LYGEMM_NOINLINE
#endif

/* Compilers stop inlining once a translation unit grows too much, which would
 * leave an unrolled tile kernel outlined and its register tile on the stack. */
#if defined(__GNUC__) && !defined(__clang__)
#  define LYGEMM_FLATTEN __attribute__((flatten))
#else
#  define LYGEMM_FLATTEN
#endif

#define LYGEMM_VERSION_MAJOR 2
#define LYGEMM_VERSION_MINOR 0
#define LYGEMM_VERSION_PATCH 0

#endif /* LYGEMM_CONFIG_HPP */
