#ifndef REF_AXPYT_HPP
#define REF_AXPYT_HPP

#include <vector>
#include <stdexcept>

template <typename T>
void axpy(T a, const std::vector<T>& x, std::vector<T>& y);

extern template void axpy<float>(float a, const std::vector<float>& x, std::vector<float>& y);
extern template void axpy<double>(double a, const std::vector<double>& x, std::vector<double>& y);

#endif