// Correctness tests for the lygemm kernels.
//
// Every buffer is mmap'ed so that its last element ends exactly at a PROT_NONE
// page.  Any out of bounds access - including one performed by a SIMD load or
// store - faults immediately instead of silently corrupting memory, which is
// what the previous hand written kernels did for shapes they did not cover.
//
// SPDX-License-Identifier: MIT
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include "lygemm.hpp"

namespace {

int g_checks = 0;
int g_failures = 0;

void report(bool ok, const std::string& what) {
  ++g_checks;
  if (!ok) {
    ++g_failures;
    std::printf("FAIL %s\n", what.c_str());
  }
}

template <class T>
double tolerance() {
  return std::is_same<T, float>::value ? 1e-5 : 1e-12;
}

template <class T>
bool close(T a, T b) {
  const double d = std::fabs(static_cast<double>(a) - static_cast<double>(b));
  return d <= tolerance<T>() * (1.0 + std::fabs(static_cast<double>(b)));
}

/* Buffer whose payload ends exactly at a guard page, so over-reads and
 * over-writes fault.  When the payload spans several pages the payload also
 * starts right after a guard page, which catches under-reads. */
template <class T>
class guarded {
 public:
  explicit guarded(std::size_t n) : n_(n) {
    const std::size_t ps = page_size();
    const std::size_t need = n * sizeof(T);
    pages_ = (need + ps - 1) / ps + 2;
    byte_len_ = pages_ * ps;
    base_ = mmap(nullptr, byte_len_, PROT_READ | PROT_WRITE,
                 MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (base_ == MAP_FAILED) {
      std::perror("mmap");
      std::abort();
    }
    char* p = static_cast<char*>(base_);
    mprotect(p, ps, PROT_NONE);
    mprotect(p + (pages_ - 1) * ps, ps, PROT_NONE);
    data_ = reinterpret_cast<T*>(p + (pages_ - 1) * ps - need);
  }
  ~guarded() { munmap(base_, byte_len_); }
  guarded(const guarded&) = delete;
  guarded& operator=(const guarded&) = delete;

  T* data() { return data_; }
  const T* data() const { return data_; }
  std::size_t size() const { return n_; }
  void fill(T v) {
    for (std::size_t i = 0; i < n_; ++i) data_[i] = v;
  }
  void fill_rand(std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t i = 0; i < n_; ++i) data_[i] = static_cast<T>(dist(gen));
  }

 private:
  static std::size_t page_size() {
    return static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
  }
  void* base_ = nullptr;
  T* data_ = nullptr;
  std::size_t n_ = 0;
  std::size_t pages_ = 0;
  std::size_t byte_len_ = 0;
};

template <class T>
void reference(std::size_t m, std::size_t n, std::size_t k, T alpha, const T* A,
               const T* B, T beta, const T* C0, T* C) {
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      T acc = T(0);
      for (std::size_t l = 0; l < k; ++l) acc += A[i * k + l] * B[l * n + j];
      C[i * n + j] = alpha * acc + (beta != T(0) ? beta * C0[i * n + j] : T(0));
    }
}

template <class T>
void check_equal(const guarded<T>& got, const guarded<T>& want,
                 const std::string& what) {
  bool ok = true;
  for (std::size_t i = 0; i < got.size(); ++i)
    if (!close<T>(got.data()[i], want.data()[i])) {
      ok = false;
      break;
    }
  report(ok, what);
}

/* ---------------------------------------------------------------- static */
template <std::size_t M, std::size_t N, std::size_t K, class T>
void check_static(T alpha, T beta, const std::string& tag) {
  std::mt19937 gen(1234u + M * 7 + N * 13 + K * 17);
  guarded<T> A(M * K), B(K * N), C(M * N), C0(M * N), R(M * N);
  A.fill_rand(gen);
  B.fill_rand(gen);
  C0.fill_rand(gen);
  C.fill(T(-7));
  for (std::size_t i = 0; i < M * N; ++i) C.data()[i] = C0.data()[i];
  lygemm::gemm<M, N, K, T>(alpha, A.data(), B.data(), beta, C.data());
  reference<T>(M, N, K, alpha, A.data(), B.data(), beta, C0.data(), R.data());
  char buf[160];
  std::snprintf(buf, sizeof buf, "%s static %zux%zux%zu alpha=%g beta=%g",
                tag.c_str(), M, N, K, double(alpha), double(beta));
  check_equal<T>(C, R, buf);
}

template <class T>
void static_suite(const std::string& tag) {
#define LYGEMM_TEST_SHAPE(M, N, K)              \
  check_static<M, N, K, T>(T(1), T(0), tag);    \
  check_static<M, N, K, T>(T(2.5), T(0.5), tag);
  LYGEMM_TEST_SHAPE(1, 1, 1)
  LYGEMM_TEST_SHAPE(2, 2, 2)
  LYGEMM_TEST_SHAPE(3, 3, 3)
  LYGEMM_TEST_SHAPE(4, 4, 4)
  LYGEMM_TEST_SHAPE(5, 5, 5)
  LYGEMM_TEST_SHAPE(6, 6, 6)
  LYGEMM_TEST_SHAPE(7, 7, 7)
  LYGEMM_TEST_SHAPE(8, 8, 8)
  LYGEMM_TEST_SHAPE(9, 9, 9)
  LYGEMM_TEST_SHAPE(12, 12, 12)
  LYGEMM_TEST_SHAPE(16, 16, 16)
  LYGEMM_TEST_SHAPE(1, 4, 4)
  LYGEMM_TEST_SHAPE(4, 1, 4)
  LYGEMM_TEST_SHAPE(4, 4, 1)
  LYGEMM_TEST_SHAPE(3, 5, 7)
  LYGEMM_TEST_SHAPE(8, 3, 5)
  LYGEMM_TEST_SHAPE(2, 8, 16)
  LYGEMM_TEST_SHAPE(16, 2, 3)
  LYGEMM_TEST_SHAPE(5, 1, 1)
  LYGEMM_TEST_SHAPE(1, 1, 5)
  LYGEMM_TEST_SHAPE(12, 4, 4)
  LYGEMM_TEST_SHAPE(33, 33, 33)
#undef LYGEMM_TEST_SHAPE
}

/* --------------------------------------------------------------- runtime */
template <class T>
void check_runtime(std::size_t m, std::size_t n, std::size_t k, T alpha, T beta,
                   const std::string& tag) {
  std::mt19937 gen(99u + m * 3 + n * 5 + k * 11);
  guarded<T> A(m * k), B(k * n), C(m * n), C0(m * n), R(m * n);
  A.fill_rand(gen);
  B.fill_rand(gen);
  C0.fill_rand(gen);
  C.fill(T(-7));
  for (std::size_t i = 0; i < m * n; ++i) C.data()[i] = C0.data()[i];
  lygemm::gemm(m, n, k, A.data(), B.data(), C.data(), alpha, beta);
  reference<T>(m, n, k, alpha, A.data(), B.data(), beta, C0.data(), R.data());
  char buf[160];
  std::snprintf(buf, sizeof buf, "%s runtime %zux%zux%zu alpha=%g beta=%g",
                tag.c_str(), m, n, k, double(alpha), double(beta));
  check_equal<T>(C, R, buf);
}

template <class T>
void runtime_suite(const std::string& tag) {
  for (std::size_t m = 1; m <= 17; ++m)
    for (std::size_t n = 1; n <= 17; ++n)
      for (std::size_t k = 1; k <= 17; ++k) check_runtime<T>(m, n, k, T(1), T(0), tag);
  const std::size_t big[] = {19, 20, 24, 31, 32, 33, 40, 48, 64};
  for (std::size_t m : big)
    for (std::size_t n : big)
      for (std::size_t k : big) check_runtime<T>(m, n, k, T(1), T(0), tag);
  check_runtime<T>(37, 5, 41, T(2), T(0.5), tag);
  check_runtime<T>(5, 37, 41, T(2), T(0.5), tag);
  check_runtime<T>(41, 41, 5, T(2), T(0.5), tag);
}

/* ----------------------------------------------------------------- batch */
template <std::size_t M, std::size_t N, std::size_t K, class T>
void check_batch(std::ptrdiff_t pad, std::ptrdiff_t sb_override,
                 const std::string& tag) {
  const std::size_t batch = 4;
  const std::ptrdiff_t sa = static_cast<std::ptrdiff_t>(M * K) + pad;
  const std::ptrdiff_t sb = sb_override >= 0
                                ? sb_override
                                : static_cast<std::ptrdiff_t>(K * N) + pad;
  const std::ptrdiff_t sc = static_cast<std::ptrdiff_t>(M * N) + pad;
  /* a zero stride broadcasts one matrix, so every buffer holds at least one
   * full matrix even when the batch span is smaller */
  const std::size_t span_a = std::max<std::size_t>(batch * sa, M * K);
  const std::size_t span_b = std::max<std::size_t>(batch * sb, K * N);
  const std::size_t span_c = std::max<std::size_t>(batch * sc, M * N);

  std::mt19937 gen(7u + M + N + K);
  guarded<T> A(span_a), B(span_b), C(span_c);
  guarded<T> R(span_c);
  A.fill_rand(gen);
  B.fill_rand(gen);
  C.fill(T(-3));
  R.fill(T(-3));
  for (std::size_t i = 0; i < span_c; ++i) R.data()[i] = C.data()[i];

  lygemm::gemm_batch<M, N, K, T>(batch, A.data(), sa, B.data(), sb, C.data(),
                                 sc, T(1.5), T(0.25));
  /* the reference walks the same layout one element at a time */
  for (std::size_t b = 0; b < batch; ++b) {
    const T* ab = A.data() + b * sa;
    const T* bb = B.data() + b * sb;
    T* cb = R.data() + b * sc;
    for (std::size_t i = 0; i < M; ++i)
      for (std::size_t j = 0; j < N; ++j) {
        T acc = T(0);
        for (std::size_t l = 0; l < K; ++l) acc += ab[i * K + l] * bb[l * N + j];
        cb[i * N + j] = T(1.5) * acc + T(0.25) * cb[i * N + j];
      }
  }
  char buf[160];
  std::snprintf(buf, sizeof buf, "%s batch %zux%zux%zu pad=%td sb=%td",
                tag.c_str(), M, N, K, pad, sb);
  check_equal<T>(C, R, buf);
}

template <class T>
void check_batch_runtime(std::size_t batch, std::size_t m, std::size_t n,
                         std::size_t k, std::ptrdiff_t pad,
                         const std::string& tag) {
  const std::ptrdiff_t sa = static_cast<std::ptrdiff_t>(m * k) + pad;
  const std::ptrdiff_t sb = static_cast<std::ptrdiff_t>(k * n) + pad;
  const std::ptrdiff_t sc = static_cast<std::ptrdiff_t>(m * n) + pad;
  std::mt19937 gen(21u + m + n + k);
  guarded<T> A(batch * sa), B(batch * sb), C(batch * sc), R(batch * sc);
  A.fill_rand(gen);
  B.fill_rand(gen);
  C.fill(T(-3));
  R.fill(T(-3));
  for (std::size_t i = 0; i < batch * sc; ++i) R.data()[i] = C.data()[i];
  lygemm::gemm_batch(batch, m, n, k, A.data(), sa, B.data(), sb, C.data(), sc,
                     T(1.5), T(0.25));
  for (std::size_t b = 0; b < batch; ++b) {
    const T* ab = A.data() + b * sa;
    const T* bb = B.data() + b * sb;
    T* cb = R.data() + b * sc;
    for (std::size_t i = 0; i < m; ++i)
      for (std::size_t j = 0; j < n; ++j) {
        T acc = T(0);
        for (std::size_t l = 0; l < k; ++l) acc += ab[i * k + l] * bb[l * n + j];
        cb[i * n + j] = T(1.5) * acc + T(0.25) * cb[i * n + j];
      }
  }
  char buf[160];
  std::snprintf(buf, sizeof buf, "%s batch(runtime) %zux%zux%zu b=%zu pad=%td",
                tag.c_str(), m, n, k, batch, pad);
  check_equal<T>(C, R, buf);
}

template <std::size_t M, std::size_t N, std::size_t K, class T>
void check_batch_shared(const std::string& tag);

template <class T>
void batch_suite(const std::string& tag) {
  check_batch<4, 4, 4, T>(0, -1, tag);
  check_batch<4, 4, 4, T>(1, -1, tag);
  check_batch<7, 3, 5, T>(3, -1, tag);
  check_batch<1, 8, 16, T>(0, -1, tag);
  check_batch<16, 16, 16, T>(5, -1, tag);
  check_batch<4, 4, 4, T>(2, 0, tag); /* B is broadcast over the batch */
  check_batch<8, 8, 8, T>(1, 0, tag);
  check_batch_runtime<T>(6, 5, 3, 7, 0, tag);
  check_batch_runtime<T>(6, 5, 3, 7, 2, tag);
  check_batch_runtime<T>(3, 12, 12, 12, 1, tag);
  check_batch_runtime<T>(3, 20, 9, 4, 0, tag);
  check_batch_runtime<T>(2, 33, 33, 33, 0, tag);
  check_batch_shared<4, 4, 4, T>(tag);
  check_batch_shared<8, 8, 8, T>(tag);
  check_batch_shared<7, 3, 5, T>(tag);
  check_batch_shared<16, 16, 16, T>(tag);
  check_batch_shared<3, 12, 6, T>(tag);
}

/* ------------------------------------------- batch with a shared operand */
template <std::size_t M, std::size_t N, std::size_t K, class T>
void check_batch_shared(const std::string& tag) {
  constexpr std::size_t batch = 5;
  std::mt19937 gen(11u + M + N + K);

  { /* A is shared by every batch element */
    guarded<T> A(M * K), B(batch * K * N), C(batch * M * N), R(batch * M * N);
    A.fill_rand(gen);
    B.fill_rand(gen);
    C.fill(T(0));
    R.fill(T(0));
    lygemm::gemm_batch_shared_a<M, N, K, T>(
        batch, A.data(), B.data(), static_cast<std::ptrdiff_t>(K * N),
        C.data(), static_cast<std::ptrdiff_t>(M * N), T(1.5), T(0.25));
    for (std::size_t b = 0; b < batch; ++b)
      reference<T>(M, N, K, T(1.5), A.data(), B.data() + b * K * N, T(0.25),
                   R.data() + b * M * N, R.data() + b * M * N);
    check_equal<T>(C, R, tag + " shared_a");
  }
  { /* B is shared by every batch element */
    guarded<T> A(batch * M * K), B(K * N), C(batch * M * N), R(batch * M * N);
    A.fill_rand(gen);
    B.fill_rand(gen);
    C.fill(T(0));
    R.fill(T(0));
    lygemm::gemm_batch_shared_b<M, N, K, T>(
        batch, A.data(), static_cast<std::ptrdiff_t>(M * K), B.data(),
        C.data(), static_cast<std::ptrdiff_t>(M * N), T(1.5), T(0.25));
    for (std::size_t b = 0; b < batch; ++b)
      reference<T>(M, N, K, T(1.5), A.data() + b * M * K, B.data(), T(0.25),
                   R.data() + b * M * N, R.data() + b * M * N);
    check_equal<T>(C, R, tag + " shared_b");
  }
}

/* --------------------------------------------------------- tile API */
/* Composes the tile primitives by hand, the way a custom kernel would. */
template <class T>
void check_tile_api(const std::string& tag) {
  constexpr std::size_t M = 8, N = 12, K = 8, BM = 4, BN = 6;
  std::mt19937 gen(3u);
  guarded<T> A(M * K), B(K * N), C(M * N), R(M * N);
  A.fill_rand(gen);
  B.fill_rand(gen);
  C.fill(T(-7));
  for (std::size_t i = 0; i < M * N; ++i) C.data()[i] = T(0);
  for (std::size_t mi = 0; mi < M / BM; ++mi)
    for (std::size_t ni = 0; ni < N / BN; ++ni) {
      lygemm::tile::fragment<T, BM, BN> frag;
      frag.clear();
      lygemm::tile::mma<K>(A.data() + mi * BM * K, K, B.data() + ni * BN, N,
                           frag);
      frag.store(C.data() + mi * BM * N + ni * BN, N, T(1), T(0));
    }
  reference<T>(M, N, K, T(1), A.data(), B.data(), T(0), R.data(), R.data());
  check_equal<T>(C, R, tag + " tile::fragment + tile::mma");

  lygemm::tile::gemm<M, N, K, T>(A.data(), B.data(), C.data());
  check_equal<T>(C, R, tag + " tile::gemm");
}

/* ------------------------------------------------------------- JIT */
/* The generated kernels are exercised directly: the runtime entry point
 * falls back to the C++ path when the generator declines a shape, so without
 * this a broken encoder could hide behind the fallback.  Each kernel runs on
 * guard page buffers, so an over-read in generated code faults. */
template <class T>
void check_jit(const std::string& tag) {
  using entry = lygemm::detail::jit::entry<T>;
  if (!lygemm::detail::jit::enabled()) return;
  int exercised = 0;
  for (std::size_t m = 1; m <= 10; ++m)
    for (std::size_t n = 1; n <= 10; ++n)
      for (std::size_t k = 1; k <= 10; ++k) {
        const auto fn = entry::get(m, n, k);
        if (fn == nullptr) continue;
        report(entry::get(m, n, k) == fn,
               tag + " jit kernel is cached for one shape");
        bool ok = true;
        for (int trial = 0; trial < 3 && ok; ++trial) {
          std::mt19937 gen(3u * unsigned(m) + 5u * unsigned(n) + unsigned(k) +
                           11u * unsigned(trial));
          guarded<T> A(m * k), B(k * n), C(m * n), R(m * n);
          A.fill_rand(gen);
          B.fill_rand(gen);
          C.fill(T(-7));
          R.fill(T(-7));
          fn(A.data(), B.data(), C.data());
          reference<T>(m, n, k, T(1), A.data(), B.data(), T(0), R.data(),
                       R.data());
          for (std::size_t i = 0; i < m * n; ++i)
            if (!close<T>(C.data()[i], R.data()[i])) {
              ok = false;
              std::printf("FAIL %s jit %zux%zux%zu trial %d element %zu\n",
                          tag.c_str(), m, n, k, trial, i);
              break;
            }
        }
        if (ok) ++exercised;
      }
  if (std::is_same<T, double>::value) /* the generator is double only so far */
    report(exercised > 500,
           tag + " jit generated a kernel for the small shapes (" +
               std::to_string(exercised) + ")");
}

/* ------------------------------------------------------- beta == 0 -> C */
/* beta == 0 must not read C: a NaN filled C has to come out clean. */
template <class T>
void check_no_read_c(const std::string& tag) {
  constexpr std::size_t M = 4, N = 4, K = 4;
  std::mt19937 gen(5u);
  guarded<T> A(M * K), B(K * N), C(M * N);
  A.fill_rand(gen);
  B.fill_rand(gen);
  C.fill(std::numeric_limits<T>::quiet_NaN());
  lygemm::gemm<M, N, K, T>(A.data(), B.data(), C.data());
  bool ok = true;
  for (std::size_t i = 0; i < M * N; ++i)
    if (!std::isfinite(static_cast<double>(C.data()[i]))) ok = false;
  report(ok, tag + " beta=0 does not read C (static 4x4x4)");

  guarded<T> A2(3 * 5), B2(5 * 7), C2(3 * 7);
  A2.fill_rand(gen);
  B2.fill_rand(gen);
  C2.fill(std::numeric_limits<T>::quiet_NaN());
  lygemm::gemm(std::size_t(3), std::size_t(7), std::size_t(5), A2.data(),
               B2.data(), C2.data());
  ok = true;
  for (std::size_t i = 0; i < 3 * 7; ++i)
    if (!std::isfinite(static_cast<double>(C2.data()[i]))) ok = false;
  report(ok, tag + " beta=0 does not read C (runtime 3x7x5)");
}

template <class T>
void run_all(const std::string& tag) {
  static_suite<T>(tag);
  runtime_suite<T>(tag);
  batch_suite<T>(tag);
  check_tile_api<T>(tag);
  check_jit<T>(tag);
  check_no_read_c<T>(tag);
}

} // namespace

int main() {
  std::printf("lygemm %d.%d.%d, LYGEMM_ARCH=%d, %zu byte vectors\n",
              LYGEMM_VERSION_MAJOR, LYGEMM_VERSION_MINOR,
              LYGEMM_VERSION_PATCH, LYGEMM_ARCH,
              sizeof(typename lygemm::detail::default_ops<double>::type::vec));
  run_all<double>("double");
  run_all<float>("float");
  std::printf("%d checks, %d failures\n", g_checks, g_failures);
  return g_failures == 0 ? 0 : 1;
}
