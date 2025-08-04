#include "ref_gemvt.hpp"
#include <stdexcept>

template <typename T>
void gemv(T a, const std::vector<std::vector<T>>& A, const std::vector<T>& x, T b, std::vector<T>& y) {
    size_t m = A.size();
    if (m == 0) throw std::invalid_argument("Matrix A is empty.");
    size_t n = A[0].size();
    for (const auto& row : A) {
        if (row.size() != n) throw std::invalid_argument("A is not a valid matrix.");
    }
    if (x.size() != n) throw std::invalid_argument("A's columns must match x's size.");
    if (y.size() != m) throw std::invalid_argument("A's rows must match y's size.");

    for (size_t i = 0; i < m; ++i) {
        T sum = 0;
        for (size_t j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        y[i] = a * sum + b * y[i];
    }
}

template void gemv<float>(float, const std::vector<std::vector<float>>&, const std::vector<float>&, float, std::vector<float>&);
template void gemv<double>(double, const std::vector<std::vector<double>>&, const std::vector<double>&, double, std::vector<double>&);