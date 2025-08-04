#include "ref_axpyt.hpp"
#include <stdexcept>

template <typename T>
void axpy(T a, const std::vector<T>& x, std::vector<T>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors x and y must be of the same size.");
    }
    if (x.size() == 0 || y.size() == 0) {
        throw std::invalid_argument("Vectors x and y can't be empty.");
    }
    for (size_t i = 0; i < x.size(); ++i) {
        y[i] += a * x[i];
    }
}

template void axpy<float>(float, const std::vector<float>&, std::vector<float>&);
template void axpy<double>(double, const std::vector<double>&, std::vector<double>&);