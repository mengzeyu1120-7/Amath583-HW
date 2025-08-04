#ifndef REF_DGEMV_HPP
#define REF_DGEMV_HPP

#include <vector>

void dgemv(double a, const std::vector<std::vector<double>>& A, const std::vector<double>& x, double b, std::vector<double>& y);

#endif