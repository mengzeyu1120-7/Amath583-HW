#ifndef REF_GEMVT_HPP
#define REF_GEMVT_HPP

#include <vector>
#include <stdexcept>

template <typename T>
void gemv(T a, const std::vector<std::vector<T>>& A, const std::vector<T>& x, T b, std::vector<T>& y);


extern template void gemv<float>(float, const std::vector<std::vector<float>>&, const std::vector<float>&, float, std::vector<float>&);
extern template void gemv<double>(double, const std::vector<std::vector<double>>&, const std::vector<double>&, double, std::vector<double>&);

#endif