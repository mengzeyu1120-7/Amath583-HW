#include "ref_dgemv.hpp"
#include <stdexcept>
#include <iostream>

void dgemv(double a, const std::vector<std::vector<double>>& A, const std::vector<double>& x, double b, std::vector<double>& y) {

    if (A.empty() || A[0].empty()) throw std::invalid_argument("Matrix A is empty.");
    size_t m = A.size();
    size_t n = A[0].size();
    if (x.size() != n) throw std::invalid_argument("A.cols != x.size");
    if (y.size() != m) throw std::invalid_argument("A.rows != y.size");


    for (size_t i = 0; i < m; ++i) {
        if (A[i].size() != n) throw std::invalid_argument("A is not a valid matrix (rows have different lengths).");
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        y[i] = a * sum + b * y[i];
    }
}