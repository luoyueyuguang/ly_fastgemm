#ifndef LYGEMM_HPP
#define LYGEMM_HPP

#include <immintrin.h>

extern "C" {
    int dgemm_(const char*   transa,
           const char*   transb,
           const int*    m,
           const int*    n,
           const int*    k,
           const double* alpha,
           const double* A,
           const int*    lda,
           const double* B,
           const int*    ldb,
           const double* beta,
           double*       C,
           const int*    ldc);
}

void ly_dgemm(const int& m,
             const int&    n,
             const int&    k,
             const double* __restrict__ A,
             const double* __restrict__ B,
             double* __restrict__ C) noexcept;

#endif // LYGEMM_HPP
