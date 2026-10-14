// lygemm benchmark: small GEMM against the libraries people reach for first.
//
//   make bench
//
// Eigen, OpenBLAS (CBLAS) and libxsmm are used when their headers and
// libraries are on the include / link path; a contender is compiled in only
// when its marker macro is defined, so the benchmark also builds standalone.
//
// Two things make the numbers trustworthy on a machine that is not quiet:
//
//   * every measured call works on a different one of eight buffer sets, so
//     the compiler cannot hoist an invariant product out of the loop, and each
//     result is read once, so it cannot delete the work either;
//   * each table re-measures the FMA ceiling next to the contenders and reports
//     the fraction of it.  This laptop drops its clock by ~3x under sustained
//     AVX2 load, so an absolute GFLOP/s figure from one moment says little
//     about the next, while the fraction of a simultaneous ceiling is stable.
//
// SPDX-License-Identifier: MIT
#include <immintrin.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "lygemm.hpp"

#ifndef LYGEMM_BENCH_BLAS
#define LYGEMM_BENCH_BLAS 0
#endif
#ifndef LYGEMM_BENCH_XSMM
#define LYGEMM_BENCH_XSMM 0
#endif

#if defined(__has_include)
#if __has_include(<Eigen/Dense>)
#define LYGEMM_BENCH_EIGEN 1
#include <Eigen/Dense>
#endif
#endif
#ifndef LYGEMM_BENCH_EIGEN
#define LYGEMM_BENCH_EIGEN 0
#endif

#if LYGEMM_BENCH_BLAS
#include <cblas.h>
#endif
#if LYGEMM_BENCH_XSMM
#include <libxsmm.h>
#endif

namespace {

using clock_type = std::chrono::steady_clock;

/* Slots per measurement: enough that the inputs differ on every call. */
constexpr std::size_t slots = 8;

/* Stops the compiler from caching values across iterations without costing a
 * memory access itself. */
inline void barrier() { __asm__ __volatile__("" ::: "memory"); }

inline double ns_since(clock_type::time_point t0) {
  return std::chrono::duration<double, std::nano>(clock_type::now() - t0).count();
}

double volatile g_sink = 0.0;

/* Accumulates one element of the output so the compiler cannot delete the work
 * it was just asked to do.  A plain per-slot array rather than a volatile
 * store per call: a volatile store adds a store-to-load chain that dominates
 * the small shapes.  The sum reaches the sink once per window. */
double g_acc_d[slots] = {};
float g_acc_f[slots] = {};
inline void consume(const double* c, std::size_t s) { g_acc_d[s] += c[0]; }
inline void consume(const float* c, std::size_t s) { g_acc_f[s] += c[0]; }
inline void flush_consumed() {
  for (std::size_t s = 0; s < slots; ++s) {
    g_sink = g_sink + g_acc_d[s] + double(g_acc_f[s]);
    g_acc_d[s] = 0.0;
    g_acc_f[s] = 0.0;
  }
}

/* Nanoseconds per call of f, where f takes the slot index of this call.  The
 * batch length is calibrated to ~50 us and replayed in many windows; the
 * fastest window is reported. */
template <class F>
__attribute__((always_inline)) inline double measure(F&& f) {
  constexpr double window_ns = 5.0e4;
  constexpr int windows = 200;
  std::size_t reps = 1;
  double elapsed = 0.0;
  for (int i = 0; i < 40; ++i) {
    const auto t0 = clock_type::now();
    for (std::size_t j = 0; j < reps; ++j) {
      std::size_t idx = j;
      __asm__ __volatile__("" : "+r"(idx));
      f(idx);
      barrier();
    }
    elapsed = ns_since(t0);
    flush_consumed();
    if (elapsed >= window_ns) break;
    const double grow = elapsed > 0.0 ? window_ns / elapsed : 64.0;
    reps = std::size_t(std::min(double(reps) * std::max(2.0, grow), 1e7));
  }
  double best = elapsed / double(reps);
  for (int w = 0; w < windows; ++w) {
    const auto t0 = clock_type::now();
    for (std::size_t j = 0; j < reps; ++j) {
      /* The index is opaque to the optimiser: without this it unrolls the loop
       * by the number of slots and precomputes every result once. */
      std::size_t idx = j;
      __asm__ __volatile__("" : "+r"(idx));
      f(idx);
      barrier();
    }
    best = std::min(best, ns_since(t0) / double(reps));
    flush_consumed();
  }
  return best;
}

/* Eight independent FMA chains: what the core achieves right now. */
__attribute__((noinline)) double fma_chain(std::size_t iters) {
  __m256d a0 = _mm256_set1_pd(1.0), a1 = _mm256_set1_pd(1.1);
  __m256d a2 = _mm256_set1_pd(1.2), a3 = _mm256_set1_pd(1.3);
  __m256d a4 = _mm256_set1_pd(1.4), a5 = _mm256_set1_pd(1.5);
  __m256d a6 = _mm256_set1_pd(1.6), a7 = _mm256_set1_pd(1.7);
  const __m256d b = _mm256_set1_pd(0.9999999);
  const __m256d c = _mm256_set1_pd(1e-12);
  for (std::size_t i = 0; i < iters; ++i) {
    a0 = _mm256_fmadd_pd(b, a0, c);
    a1 = _mm256_fmadd_pd(b, a1, c);
    a2 = _mm256_fmadd_pd(b, a2, c);
    a3 = _mm256_fmadd_pd(b, a3, c);
    a4 = _mm256_fmadd_pd(b, a4, c);
    a5 = _mm256_fmadd_pd(b, a5, c);
    a6 = _mm256_fmadd_pd(b, a6, c);
    a7 = _mm256_fmadd_pd(b, a7, c);
  }
  const __m256d s = _mm256_add_pd(_mm256_add_pd(a0, a1), _mm256_add_pd(a2, a3));
  g_sink = g_sink + _mm256_cvtsd_f64(s);
  (void)a4; (void)a5; (void)a6; (void)a7;
  return double(iters) * 8.0 * 4.0 * 2.0; /* chains, lanes, flops per FMA */
}

double peak_gflops() {
  const double ns = measure([](std::size_t) { fma_chain(2000); });
  return (2000.0 * 8.0 * 4.0 * 2.0) / ns;
}

/* Effective core clock, from a chain of dependent adds (one per cycle). */
double core_clock_ghz() {
  constexpr std::uint64_t n = 5e7;
  std::uint64_t acc = 0;
  const auto t0 = clock_type::now();
  for (std::uint64_t i = 0; i < n; ++i)
    __asm__ __volatile__("addq $1, %0" : "+r"(acc));
  return double(n) / ns_since(t0);
}

void print_row(const char* shape, const char* operand, const char* method,
               double ns, double flops, double peak) {
  std::printf("%-12s %-8s %-26s %10.2f %12.2f %7.1f%%\n", shape, operand,
              method, ns, flops / ns, 100.0 * (flops / ns) / peak);
}

template <class T>
void fill_random(std::vector<T>& v, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (auto& x : v) x = static_cast<T>(dist(gen));
}

template <class T>
void naive_single(std::size_t m, std::size_t n, std::size_t k, const T* A,
                  const T* B, T* C) {
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      T acc = T(0);
      for (std::size_t l = 0; l < k; ++l) acc += A[i * k + l] * B[l * n + j];
      C[i * n + j] = acc;
    }
}

/* ------------------------------------------------------------ single GEMM */
template <std::size_t M, std::size_t N, std::size_t K, class T>
__attribute__((noinline)) void bench_single(const char* shape) {
  const double flops = 2.0 * double(M) * double(N) * double(K);
  std::vector<T> A(slots * M * K), B(slots * K * N), C(slots * M * N, T(0));
  fill_random(A, 1u + M);
  fill_random(B, 2u + N * 3 + K);
  const T* a = A.data();
  const T* b = B.data();
  T* c = C.data();

  print_row(shape, "-", "naive",
            measure([&](std::size_t i) {
              const std::size_t s = i & (slots - 1);
              naive_single(M, N, K, a + s * M * K, b + s * K * N, c + s * M * N);
            }),
            flops, peak_gflops());
  print_row(shape, "-", "lygemm gemm<M,N,K>",
            measure([&](std::size_t i) {
              const std::size_t s = i & (slots - 1);
              lygemm::gemm<M, N, K, T>(a + s * M * K, b + s * K * N,
                                       c + s * M * N);
              consume(c + s * M * N, s);
            }),
            flops, peak_gflops());
  print_row(shape, "-", "lygemm gemm(m,n,k)",
            measure([&](std::size_t i) {
              const std::size_t s = i & (slots - 1);
              lygemm::gemm(std::size_t(M), std::size_t(N), std::size_t(K),
                           a + s * M * K, b + s * K * N, c + s * M * N);
              consume(c + s * M * N, s);
            }),
            flops, peak_gflops());

#if LYGEMM_BENCH_EIGEN
  {
    Eigen::Matrix<T, M, K> ea[slots];
    Eigen::Matrix<T, K, N> eb[slots];
    Eigen::Matrix<T, M, N> ec[slots];
    for (std::size_t s = 0; s < slots; ++s) {
      for (std::size_t e = 0; e < M * K; ++e) ea[s].data()[e] = a[s * M * K + e];
      for (std::size_t e = 0; e < K * N; ++e) eb[s].data()[e] = b[s * K * N + e];
    }
    print_row(shape, "-", "Eigen fixed",
              measure([&](std::size_t i) {
                const std::size_t s = i & (slots - 1);
                ec[s].noalias() = ea[s] * eb[s];
                consume(ec[s].data(), s);
              }),
              flops, peak_gflops());
    std::vector<Eigen::MatrixXd> da(slots, Eigen::MatrixXd(M, K));
    std::vector<Eigen::MatrixXd> db(slots, Eigen::MatrixXd(K, N));
    std::vector<Eigen::MatrixXd> dc(slots, Eigen::MatrixXd(M, N));
    for (std::size_t s = 0; s < slots; ++s) {
      for (std::size_t e = 0; e < M * K; ++e) da[s].data()[e] = a[s * M * K + e];
      for (std::size_t e = 0; e < K * N; ++e) db[s].data()[e] = b[s * K * N + e];
    }
    print_row(shape, "-", "Eigen dynamic",
              measure([&](std::size_t i) {
                const std::size_t s = i & (slots - 1);
                dc[s].noalias() = da[s] * db[s];
                consume(dc[s].data(), s);
              }),
              flops, peak_gflops());
  }
#endif

#if LYGEMM_BENCH_BLAS
  {
    const T one = T(1), zero = T(0);
    print_row(shape, "-", "OpenBLAS cblas",
              measure([&](std::size_t i) {
                const std::size_t s = i & (slots - 1);
                if constexpr (std::is_same<T, double>::value)
                  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, int(M),
                              int(N), int(K), one, a + s * M * K, int(K),
                              b + s * K * N, int(N), zero, c + s * M * N,
                              int(N));
                else
                  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, int(M),
                              int(N), int(K), one, a + s * M * K, int(K),
                              b + s * K * N, int(N), zero, c + s * M * N,
                              int(N));
                consume(c + s * M * N, s);
              }),
              flops, peak_gflops());
  }
#endif

#if LYGEMM_BENCH_XSMM
  {
    libxsmm_mmfunction<T> kernel(LIBXSMM_GEMM_FLAG_NONE, int(M), int(N), int(K),
                                 1.0, 0.0);
    if (kernel)
      print_row(shape, "-", "libxsmm",
                measure([&](std::size_t i) {
                  const std::size_t s = i & (slots - 1);
                  kernel(a + s * M * K, b + s * K * N, c + s * M * N);
                  consume(c + s * M * N, s);
                }),
                flops, peak_gflops());
    else
      std::printf("%-12s %-8s %-26s %10s %12s %8s\n", shape, "-", "libxsmm",
                  "n/a", "n/a", "n/a");
  }
#endif
  std::printf("\n");
}

/* ------------------------------------------------------------- batched */
template <std::size_t M, std::size_t N, std::size_t K, class T>
__attribute__((noinline)) void bench_batch(const char* shape,
                                           std::size_t batch) {
  const double flops = 2.0 * double(M) * double(N) * double(K);
  std::vector<T> A(batch * M * K), B(batch * K * N), C(batch * M * N, T(0));
  fill_random(A, 7u + M);
  fill_random(B, 9u + N * 5 + K);
  const T* a = A.data();
  const T* b = B.data();
  T* c = C.data();
  const auto per_item = [&](double ns) { return ns / double(batch); };

  /* the batch itself is the rotation: every item touches different matrices */
  print_row(shape, "indep", "naive loop",
            per_item(measure([&](std::size_t) {
              for (std::size_t i = 0; i < batch; ++i)
                naive_single(M, N, K, a + i * M * K, b + i * K * N,
                             c + i * M * N);
              consume(c, 0);
            })),
            flops, peak_gflops());
  print_row(shape, "indep", "lygemm gemm_batch<M,N,K>",
            per_item(measure([&](std::size_t) {
              lygemm::gemm_batch<M, N, K, T>(batch, a, b, c);
              consume(c, 0);
            })),
            flops, peak_gflops());
  print_row(shape, "indep", "lygemm gemm_batch(m,n,k)",
            per_item(measure([&](std::size_t) {
              lygemm::gemm_batch(batch, std::size_t(M), std::size_t(N),
                                 std::size_t(K), a, b, c);
              consume(c, 0);
            })),
            flops, peak_gflops());

  /* A shared by the whole batch: one M x K operand, `batch` different B */
  print_row(shape, "A", "gemm_batch(sa=0)",
            per_item(measure([&](std::size_t) {
              lygemm::gemm_batch<M, N, K, T>(batch, a, std::ptrdiff_t(0), b,
                                             std::ptrdiff_t(K * N), c,
                                             std::ptrdiff_t(M * N));
              consume(c, 0);
            })),
            flops, peak_gflops());
  print_row(shape, "A", "gemm_batch_shared_a",
            per_item(measure([&](std::size_t) {
              lygemm::gemm_batch_shared_a<M, N, K, T>(
                  batch, a, b, std::ptrdiff_t(K * N), c,
                  std::ptrdiff_t(M * N));
              consume(c, 0);
            })),
            flops, peak_gflops());

  /* B shared by the whole batch: one K x N operand, `batch` different A */
  print_row(shape, "B", "gemm_batch(sb=0)",
            per_item(measure([&](std::size_t) {
              lygemm::gemm_batch<M, N, K, T>(batch, a, std::ptrdiff_t(M * K), b,
                                             std::ptrdiff_t(0), c,
                                             std::ptrdiff_t(M * N));
              consume(c, 0);
            })),
            flops, peak_gflops());
  print_row(shape, "B", "gemm_batch_shared_b",
            per_item(measure([&](std::size_t) {
              lygemm::gemm_batch_shared_b<M, N, K, T>(
                  batch, a, std::ptrdiff_t(M * K), b, c,
                  std::ptrdiff_t(M * N));
              consume(c, 0);
            })),
            flops, peak_gflops());

#if LYGEMM_BENCH_BLAS
  {
    const T one = T(1), zero = T(0);
    print_row(shape, "indep", "OpenBLAS loop",
              per_item(measure([&](std::size_t) {
                for (std::size_t i = 0; i < batch; ++i) {
                  if constexpr (std::is_same<T, double>::value)
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                int(M), int(N), int(K), one, a + i * M * K,
                                int(K), b + i * K * N, int(N), zero,
                                c + i * M * N, int(N));
                  else
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                int(M), int(N), int(K), one, a + i * M * K,
                                int(K), b + i * K * N, int(N), zero,
                                c + i * M * N, int(N));
                }
                consume(c, 0);
              })),
              flops, peak_gflops());
  }
#endif

#if LYGEMM_BENCH_XSMM
  {
    libxsmm_mmfunction<T> kernel(LIBXSMM_GEMM_FLAG_NONE, int(M), int(N), int(K),
                                 1.0, 0.0);
    if (kernel)
      print_row(shape, "indep", "libxsmm loop",
                per_item(measure([&](std::size_t) {
                  for (std::size_t i = 0; i < batch; ++i)
                    kernel(a + i * M * K, b + i * K * N, c + i * M * N);
                  consume(c, 0);
                })),
                flops, peak_gflops());
  }
#endif
  std::printf("\n");
}

template <class T>
void run_single(const char* tag) {
  std::printf("== %s, single GEMM, ns/call  (core %.2f GHz) ==\n", tag,
              core_clock_ghz());
  std::printf("%-12s %-8s %-26s %10s %12s %8s\n", "shape", "operand", "method",
              "ns/call", "GFLOP/s", "of peak");
#define LYGEMM_BENCH_SHAPE(M, N, K) bench_single<M, N, K, T>(#M "x" #N "x" #K);
  LYGEMM_BENCH_SHAPE(2, 2, 2)
  LYGEMM_BENCH_SHAPE(3, 3, 3)
  LYGEMM_BENCH_SHAPE(4, 4, 4)
  LYGEMM_BENCH_SHAPE(4, 4, 16)
  LYGEMM_BENCH_SHAPE(6, 6, 6)
  LYGEMM_BENCH_SHAPE(8, 8, 8)
  LYGEMM_BENCH_SHAPE(8, 8, 64)
  LYGEMM_BENCH_SHAPE(12, 12, 12)
  LYGEMM_BENCH_SHAPE(16, 16, 16)
  LYGEMM_BENCH_SHAPE(16, 16, 64)
  LYGEMM_BENCH_SHAPE(24, 24, 24)
  LYGEMM_BENCH_SHAPE(32, 32, 32)
  LYGEMM_BENCH_SHAPE(64, 64, 64)
#undef LYGEMM_BENCH_SHAPE
}

template <class T>
void run_batch(const char* tag, std::size_t batch) {
  std::printf("== %s, batch of %zu, ns per item  (core %.2f GHz) ==\n", tag,
              batch, core_clock_ghz());
  std::printf("%-12s %-8s %-26s %10s %12s %8s\n", "shape", "operand", "method",
              "ns/item", "GFLOP/s", "of peak");
#define LYGEMM_BENCH_SHAPE(M, N, K) bench_batch<M, N, K, T>(#M "x" #N "x" #K, batch);
  LYGEMM_BENCH_SHAPE(2, 2, 2)
  LYGEMM_BENCH_SHAPE(4, 4, 4)
  LYGEMM_BENCH_SHAPE(8, 8, 8)
  LYGEMM_BENCH_SHAPE(12, 12, 12)
  LYGEMM_BENCH_SHAPE(16, 16, 16)
#undef LYGEMM_BENCH_SHAPE
}

} // namespace

int main(int argc, char** argv) {
  std::size_t batch = 4000;
  if (argc > 1) batch = std::size_t(std::atol(argv[1]));
  std::printf("lygemm %d.%d.%d, LYGEMM_ARCH=%d  (Eigen=%d BLAS=%d libxsmm=%d)\n",
              LYGEMM_VERSION_MAJOR, LYGEMM_VERSION_MINOR,
              LYGEMM_VERSION_PATCH, LYGEMM_ARCH, LYGEMM_BENCH_EIGEN,
              LYGEMM_BENCH_BLAS, LYGEMM_BENCH_XSMM);
  std::printf("ceiling now: %.1f GFLOP/s\n\n", peak_gflops());
  run_single<double>("double");
  run_batch<double>("double", batch);
  run_single<float>("float");
  return 0;
}
