#include "ref_gemmt.hpp"
#include <stdexcept>

template <typename T>
void gemm(T a, const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B, T b, std::vector<std::vector<T>>& C) {
    size_t m = A.size();
    if (m == 0) throw std::invalid_argument("Matrix A is empty.");
    size_t p = A[0].size();
    for (const auto& row : A) {
        if (row.size() != p) throw std::invalid_argument("A is not a valid matrix.");
    }

    size_t p_B = B.size();
    if (p_B != p) throw std::invalid_argument("A's columns must match B's rows.");
    if (p_B == 0) throw std::invalid_argument("Matrix B is empty.");
    size_t n = B[0].size();
    for (const auto& row : B) {
        if (row.size() != n) throw std::invalid_argument("B is not a valid matrix.");
    }

    if (C.size() != m) throw std::invalid_argument("C's rows must match A's rows.");
    for (const auto& row : C) {
        if (row.size() != n) throw std::invalid_argument("C's columns must match B's columns.");
    }

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            T sum = 0;
            for (size_t k = 0; k < p; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = a * sum + b * C[i][j];
        }
    }
}

template void gemm<float>(float, const std::vector<std::vector<float>>&, const std::vector<std::vector<float>>&, float, std::vector<std::vector<float>>&);
template void gemm<double>(double, const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, double, std::vector<std::vector<double>>&);