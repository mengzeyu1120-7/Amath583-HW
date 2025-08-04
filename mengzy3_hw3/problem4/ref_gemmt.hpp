#ifndef REF_GEMMT_HPP
#define REF_GEMMT_HPP

#include <vector>
#include <stdexcept>

template <typename T>
void gemm(T a, const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B, T b, std::vector<std::vector<T>>& C);

extern template void gemm<float>(float, const std::vector<std::vector<float>>&, const std::vector<std::vector<float>>&, float, std::vector<std::vector<float>>&);
extern template void gemm<double>(double, const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, double, std::vector<std::vector<double>>&);

#endif