#include <vector>
#include <algorithm>

void swapRows(std::vector<double> &matrix, int nRows, int nCols, int i, int j) {
    for (int c = 0; c < nCols; ++c) {
        std::swap(matrix[c * nRows + i], matrix[c * nRows + j]);
    }
}

void swapCols(std::vector<double> &matrix, int nRows, int nCols, int i, int j) {
    auto col_i_start = matrix.begin() + i * nRows;
    auto col_j_start = matrix.begin() + j * nRows;
    std::swap_ranges(col_i_start, col_i_start + nRows, col_j_start);
}