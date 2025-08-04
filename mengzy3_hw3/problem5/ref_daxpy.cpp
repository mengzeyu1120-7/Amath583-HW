#include "ref_daxpy.hpp"
#include <stdexcept>
#include <iostream>


void daxpy(double a, const std::vector<double>& x, std::vector<double>& y) {
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