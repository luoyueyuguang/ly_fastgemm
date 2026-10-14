/* Minimal CBLAS declarations, used by the benchmark only.
 *
 * A real installation gets this header from the BLAS development package
 * (cblas.h).  This container has the OpenBLAS library but no BLAS headers, so
 * the small subset the benchmark needs is declared here; the CBLAS ABI below
 * is the standard one and has not changed.
 */
#ifndef LYGEMM_BENCH_CBLAS_H
#define LYGEMM_BENCH_CBLAS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { CblasRowMajor = 101, CblasColMajor = 102 } CBLAS_LAYOUT;
typedef enum {
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113
} CBLAS_TRANSPOSE;

void cblas_dgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa,
                 const CBLAS_TRANSPOSE transb, const int m, const int n,
                 const int k, const double alpha, const double* a,
                 const int lda, const double* b, const int ldb,
                 const double beta, double* c, const int ldc);

void cblas_sgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa,
                 const CBLAS_TRANSPOSE transb, const int m, const int n,
                 const int k, const float alpha, const float* a,
                 const int lda, const float* b, const int ldb,
                 const float beta, float* c, const int ldc);

#ifdef __cplusplus
}
#endif

#endif /* LYGEMM_BENCH_CBLAS_H */
