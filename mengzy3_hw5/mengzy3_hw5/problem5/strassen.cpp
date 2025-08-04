// @uw.edu
// AMATH 483-583
// strassen.cpp : starter code for Strassen implementation

#include <iostream>
#include <vector>

using namespace std;

template <typename T>
vector<vector<T>> addMatrix(const vector<vector<T>> &A, const vector<vector<T>> &B)
{
    int n = A.size();
    int m = A[0].size();
    vector<vector<T>> C(n, vector<T>(m));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

template <typename T>
vector<vector<T>> subtractMatrix(const vector<vector<T>> &A, const vector<vector<T>> &B)
{
    int n = A.size();
    int m = A[0].size();
    vector<vector<T>> C(n, vector<T>(m));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

template <typename T>
vector<vector<T>> padMatrix(const vector<vector<T>> &matrix, int new_size) {
    int original_size = matrix.size();
    vector<vector<T>> padded(new_size, vector<T>(new_size, T(0)));
    for (int i = 0; i < original_size; ++i)
        for (int j = 0; j < original_size; ++j)
            padded[i][j] = matrix[i][j];
    return padded;
}

template <typename T>
vector<vector<T>> trimMatrix(const vector<vector<T>> &matrix, int original_size) {
    vector<vector<T>> trimmed(original_size, vector<T>(original_size));
    for (int i = 0; i < original_size; ++i)
        for (int j = 0; j < original_size; ++j)
            trimmed[i][j] = matrix[i][j];
    return trimmed;
}

template <typename T>
vector<vector<T>> getSubmatrix(const vector<vector<T>> &matrix, int start_row, int start_col, int rows, int cols) {
    vector<vector<T>> sub(rows, vector<T>(cols, T(0)));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (start_row + i < matrix.size() && start_col + j < matrix[0].size())
                sub[i][j] = matrix[start_row + i][start_col + j];
    return sub;
}

template <typename T>
vector<vector<T>> mergeSubmatrices(const vector<vector<T>> &C11, const vector<vector<T>> &C12,
                                   const vector<vector<T>> &C21, const vector<vector<T>> &C22) {
    int n = C11.size();
    int m = C11[0].size();
    vector<vector<T>> C(2 * n, vector<T>(2 * m));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = C11[i][j];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j + m] = C12[i][j];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i + n][j] = C21[i][j];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i + n][j + m] = C22[i][j];
    return C;
}

template <typename T>
vector<vector<T>> strassenMultiply(const vector<vector<T>> &A, const vector<vector<T>> &B) {
    int original_n = A.size();
    int n = original_n;
    bool is_padded = false;
    
    // Pad matrices to even dimensions if necessary
    if (n % 2 != 0) {
        n += 1;
        is_padded = true;
    }
    
    vector<vector<T>> A_padded = A;
    vector<vector<T>> B_padded = B;
    if (is_padded) {
        A_padded = padMatrix(A, n);
        B_padded = padMatrix(B, n);
    }
    
    // Base case
    if (n == 1) {
        return {{A_padded[0][0] * B_padded[0][0]}};
    } else if (n == 2) {
        T a = A_padded[0][0], b = A_padded[0][1],
          c = A_padded[1][0], d = A_padded[1][1];
        T e = B_padded[0][0], f = B_padded[0][1],
          g = B_padded[1][0], h = B_padded[1][1];
        
        T M1 = (a + d) * (e + h);
        T M2 = (c + d) * e;
        T M3 = a * (f - h);
        T M4 = d * (g - e);
        T M5 = (a + b) * h;
        T M6 = (c - a) * (e + f);
        T M7 = (b - d) * (g + h);
        
        vector<vector<T>> C = {
            {M1 + M4 - M5 + M7, M3 + M5},
            {M2 + M4, M1 - M2 + M3 + M6}
        };
        if (is_padded) C = trimMatrix(C, original_n);
        return C;
    }
    
    // Split matrices into submatrices
    int half = n / 2;
    auto A11 = getSubmatrix(A_padded, 0, 0, half, half);
    auto A12 = getSubmatrix(A_padded, 0, half, half, half);
    auto A21 = getSubmatrix(A_padded, half, 0, half, half);
    auto A22 = getSubmatrix(A_padded, half, half, half, half);
    
    auto B11 = getSubmatrix(B_padded, 0, 0, half, half);
    auto B12 = getSubmatrix(B_padded, 0, half, half, half);
    auto B21 = getSubmatrix(B_padded, half, 0, half, half);
    auto B22 = getSubmatrix(B_padded, half, half, half, half);
    
    // Compute intermediate matrices
    auto S1 = subtractMatrix(B12, B22);
    auto S2 = addMatrix(A11, A12);
    auto S3 = addMatrix(A21, A22);
    auto S4 = subtractMatrix(B21, B11);
    auto S5 = addMatrix(A11, A22);
    auto S6 = addMatrix(B11, B22);
    auto S7 = subtractMatrix(A12, A22);
    auto S8 = addMatrix(B21, B22);
    auto S9 = subtractMatrix(A21, A11);
    auto S10 = addMatrix(B11, B12);
    
    // Recursive multiplication
    auto P1 = strassenMultiply(S5, S6);
    auto P2 = strassenMultiply(S3, B11);
    auto P3 = strassenMultiply(A11, S1);
    auto P4 = strassenMultiply(A22, S4);
    auto P5 = strassenMultiply(S2, B22);
    auto P6 = strassenMultiply(S9, S10);
    auto P7 = strassenMultiply(S7, S8);
    
    // Compute C submatrices
    auto C11 = subtractMatrix(addMatrix(addMatrix(P1, P4), P7), P5);  // M1 + M4 - M5 + M7
    auto C12 = addMatrix(P3, P5);                                     // M3 + M5
    auto C21 = addMatrix(P2, P4);                                     // M2 + M4
    auto C22 = addMatrix(subtractMatrix(P1, P2), addMatrix(P3, P6));  // M1 - M2 + M3 + M6
    
    // Merge and trim
    vector<vector<T>> C_padded = mergeSubmatrices(C11, C12, C21, C22);
    if (is_padded) C_padded = trimMatrix(C_padded, original_n);
    return C_padded;
}

template <typename T>
void printMatrix(const vector<vector<T>> &matrix)
{
    int n = matrix.size();
    int m = matrix[0].size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// int
template vector<vector<int>> addMatrix<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template vector<vector<int>> subtractMatrix<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template vector<vector<int>> strassenMultiply<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template void printMatrix<int>(const vector<vector<int>> &matrix);
// double
template vector<vector<double>> addMatrix<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template vector<vector<double>> subtractMatrix<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template vector<vector<double>> strassenMultiply<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template void printMatrix<double>(const vector<vector<double>> &matrix);
