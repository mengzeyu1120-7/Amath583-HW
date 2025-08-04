#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>
#include <limits>
#include <cblas.h>
#include <lapacke.h>

int main() {
    std::vector<int> sizes = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};

    for (int n : sizes) {
        std::complex<double>* A = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n * n);
        std::complex<double>* b = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n);
        std::complex<double>* z = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n);
        int* ipiv = (int*)malloc(sizeof(int) * n);

        // Initialize A
        srand(0);
        int k = 0;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                A[k] = 0.5 - (double)rand() / RAND_MAX
                     + std::complex<double>(0, 1) * (0.5 - (double)rand() / RAND_MAX);
                if (i == j) A[k] *= static_cast<double>(n);
                k++;
            }
        }

        // Initialize b
        srand(1);
        for (int i = 0; i < n; i++) {
            b[i] = 0.5 - (double)rand() / RAND_MAX
                 + std::complex<double>(0, 1) * (0.5 - (double)rand() / RAND_MAX);
        }

        // Copy b to z (solution vector)
        memcpy(z, b, sizeof(std::complex<double>) * n);

        // Solve A z = b using LAPACK
        int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, 1,
                                 reinterpret_cast<lapack_complex_double*>(A), n,
                                 ipiv,
                                 reinterpret_cast<lapack_complex_double*>(z), 1);
        if (info != 0) {
            std::cerr << "LAPACKE_zgesv failed with info = " << info << " at size " << n << std::endl;
            continue;
        }

        // Compute Az
        std::vector<std::complex<double>> Az(n);
        std::complex<double> alpha(1.0, 0.0);
        std::complex<double> beta(0.0, 0.0);
        cblas_zgemv(CblasRowMajor, CblasNoTrans, n, n,
                    &alpha,
                    A, n,
                    z, 1,
                    &beta,
                    Az.data(), 1);

        // Compute ||b - Az||_2
        double residual = 0.0;
        for (int i = 0; i < n; ++i) {
            std::complex<double> diff = b[i] - Az[i];
            residual += std::norm(diff);
        }
        residual = std::sqrt(residual);

        // Compute ||z||_2
        double norm_z = 0.0;
        for (int i = 0; i < n; ++i) norm_z += std::norm(z[i]);
        norm_z = std::sqrt(norm_z);

        // Compute ||A||_inf
        double norm_inf = 0.0;
        for (int i = 0; i < n; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < n; j++) {
                row_sum += std::abs(A[i * n + j]);
            }
            norm_inf = std::max(norm_inf, row_sum);
        }

        double epsilon = std::numeric_limits<double>::epsilon();
        double normalized_error = residual / (norm_inf * norm_z * epsilon);

        std::cout << n << " "
                  << std::log10(residual) << " "
                  << std::log10(normalized_error) << std::endl;

        free(A);
        free(b);
        free(z);
        free(ipiv);
    }

    return 0;
}
