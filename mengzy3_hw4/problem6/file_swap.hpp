// file_swaps.hpp
#include <fstream>
#include <vector>

void swapRowsInFile(std::fstream &file, int nRows, int nCols, int i, int j) {
    if (i < 0 || i >= nRows || j < 0 || j >= nRows) {
        throw std::out_of_range("Row indices out of bounds");
    }

    double temp_i, temp_j;
    for (int col = 0; col < nCols; ++col) {
    
        std::streampos pos_i = col * nRows * sizeof(double) + i * sizeof(double);
        std::streampos pos_j = col * nRows * sizeof(double) + j * sizeof(double);

        // read row i
        file.seekg(pos_i);
        file.read(reinterpret_cast<char*>(&temp_i), sizeof(double));

        // read row j
        file.seekg(pos_j);
        file.read(reinterpret_cast<char*>(&temp_j), sizeof(double));

        // swap and write
        file.seekp(pos_i);
        file.write(reinterpret_cast<const char*>(&temp_j), sizeof(double));

        file.seekp(pos_j);
        file.write(reinterpret_cast<const char*>(&temp_i), sizeof(double));
    }
    file.flush();
}

void swapColsInFile(std::fstream &file, int nRows, int nCols, int i, int j) {
    if (i < 0 || i >= nCols || j < 0 || j >= nCols) {
        throw std::out_of_range("Column indices out of bounds");
    }

    std::vector<double> buffer_i(nRows);
    std::vector<double> buffer_j(nRows);

    // read column i
    std::streampos start_i = i * nRows * sizeof(double);
    file.seekg(start_i);
    file.read(reinterpret_cast<char*>(buffer_i.data()), nRows * sizeof(double));

    // read column j
    std::streampos start_j = j * nRows * sizeof(double);
    file.seekg(start_j);
    file.read(reinterpret_cast<char*>(buffer_j.data()), nRows * sizeof(double));

    // swap and write
    file.seekp(start_i);
    file.write(reinterpret_cast<const char*>(buffer_j.data()), nRows * sizeof(double));

    file.seekp(start_j);
    file.write(reinterpret_cast<const char*>(buffer_i.data()), nRows * sizeof(double));

    file.flush();
}