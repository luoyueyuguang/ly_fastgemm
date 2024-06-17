#include "lygemm.hpp"

void ly_dgemm(const int& m,
             const int&    n,
             const int&    k,
             const double* __restrict__ A,
             const double* __restrict__ B,
             double* __restrict__ C) noexcept
{
    /*
     * This matrix multiply is optimized for m = 4 and m = 2
     * I will try to optimize it more in the future
     *
     * The idea is to use the fact that the matrix is small and we can unroll the loop
     * and use the fact that the matrix is small to use the permute4x64_pd to do the
     * multiplication
     *
     * C = A * B
     *
     * 1. m = 4
     * 2. m = 2
     * 3. other cases is handled by blas
    */
    if(m == 4){
    auto tmp1 = _mm256_setzero_pd();
    auto tmp2 = _mm256_setzero_pd();
    auto tmp3 = _mm256_setzero_pd();
    auto tmp4 = _mm256_setzero_pd();
    for(size_t l = 0; l < k; l += 4)
    {
        const auto atmp1 = _mm256_castpd_si256(_mm256_loadu_pd(A + l));
        const auto atmp2 = _mm256_castpd_si256(_mm256_loadu_pd(A + l + k));
        const auto atmp3 = _mm256_castpd_si256(_mm256_loadu_pd(A + l + 2 * k));
        const auto atmp4 = _mm256_castpd_si256(_mm256_loadu_pd(A + l + 3 * k));

        const auto ba0 = _mm256_loadu_pd(B + l * n);	
        const auto ba1 = _mm256_loadu_pd(B + (l + 1) * n);
        const auto ba2 = _mm256_loadu_pd(B + (l + 2) * n);
        const auto ba3 = _mm256_loadu_pd(B + (l + 3) * n);

        tmp1 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp1, 0x00), ba0, tmp1);
        tmp2 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp2, 0x55), ba1, tmp2);
        tmp3 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp3, 0xaa), ba2, tmp3);
        tmp4 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp4, 0xff), ba3, tmp4);

        tmp2 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp2, 0x00), ba0, tmp2);
        tmp1 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp1, 0x55), ba1, tmp1);
        tmp3 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp3, 0xff), ba3, tmp3);
        tmp4 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp4, 0xaa), ba2, tmp4);
        

        tmp3 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp3, 0x00), ba0, tmp3);
        tmp2 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp2, 0xaa), ba2, tmp2);
        tmp4 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp4, 0x55), ba1, tmp4);
        tmp1 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp1, 0xaa), ba2, tmp1);
        
        
        tmp4 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp4, 0x00), ba0, tmp4);
        tmp3 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp3, 0x55), ba1, tmp3);
        tmp2 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp2, 0xff), ba3, tmp2);
        tmp1 = _mm256_fmadd_pd(_mm256_permute4x64_pd(atmp1, 0xff), ba3, tmp1);
    }
    _mm256_storeu_pd(C, _mm256_add_pd(tmp1, _mm256_loadu_pd(C)));
    _mm256_storeu_pd(C + 1 * n, _mm256_add_pd(tmp2, _mm256_loadu_pd(C + 1 * n)));
    _mm256_storeu_pd(C + 2 * n, _mm256_add_pd(tmp3, _mm256_loadu_pd(C + 2 * n)));
    _mm256_storeu_pd(C + 3 * n, _mm256_add_pd(tmp4, _mm256_loadu_pd(C + 3 * n)));
    }
    else if(m == 2){
            auto ctmp = _mm_loadu_pd(C);
            auto ctmp1 = _mm_loadu_pd(C + 1 * n);
            #pragma unroll
            for(auto l = 0; l < 4; ++l){
                auto atmp1 = _mm_set1_pd(A[0 * k + l]);
                auto atmp2 = _mm_set1_pd(A[1 * k + l]);
                ctmp = _mm_fmadd_pd(atmp1, _mm_loadu_pd(B + l * n), ctmp);
                ctmp1 = _mm_fmadd_pd(atmp2, _mm_loadu_pd(B + l * n), ctmp1);
            }
            _mm_storeu_pd(C, ctmp);
            _mm_storeu_pd(C + 1 * n, ctmp1);
    }
    else{
         auto alpha = 1.0;
        auto beta = 1.0;
        const char transa = 'N', transb = 'N';
        dgemm_(&transa, &transb, &n, &m, &k, &alpha, B, &n, A, &k, &beta, C, &n);
    }
}
