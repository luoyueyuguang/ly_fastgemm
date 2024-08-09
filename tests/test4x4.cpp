#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>

#include "lygemm.hpp"

int main() {
    std::vector<double> A = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };

    std::vector<double> B = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };

    std::vector<double> C(16, 0);
    std::vector<double> C2(16, 0);

    auto alpha = 1.0;
    auto beta = 1.0;
    const char transa = 'N', transb = 'N';
    auto n = 4, m = 4, k = 4;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    ly_dgemm(4, 4, 4, A.data(), B.data(), C.data());
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "\033[1;35mTime difference in 4x4 = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]\033[0m" << std::endl;

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    dgemm_(&transa, &transb, &n, &m, &k, &alpha, B.data(), &n, A.data(), &k, &beta, C2.data(), &n);
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "\033[1;35mTime difference in 4x4 = " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[µs]\033[0m" << std::endl;

    std::cout << "\033[1;36mResult of lygemm\033[0m" << std::endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            std::cout << C[i * 4 + j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\033[1;36mResult of dgemm\033[0m" << std::endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            std::cout << C2[i * 4 + j] << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < 16; i++) {
        if (std::abs(C[i] - C2[i]) > 1e-6) {
            std::cout << "\033[1;31mTest failed\033[0m" << std::endl;
            return 1;
        }
    }

    std::cout << "\033[1;34mTest passed\033[0m" << std::endl;

    return 0;
}
