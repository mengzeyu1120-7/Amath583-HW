#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cblas.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

const int NTRIAL = 3;
const double alpha = 1.0;
const double beta = 0.0;

int main() {
    std::vector<int> n_values;
    for (int n = 2; n <= 16384; n *= 2) {
        n_values.push_back(n);
    }

    std::ofstream outfile("results.txt");
    if (!outfile.is_open()) {
        std::cerr << "Failed to open results.txt!" << std::endl;
        return 1;
    }

    outfile << "Performance Results (GFLOPs)\n";
    outfile << "----------------------------\n";

    for (int n : n_values) {
        size_t size = n * n;
        double *A = new double[size];
        double *B = new double[size];
        double *C = new double[size];

        // Initialize matrices
        for (size_t i = 0; i < size; ++i) {
            A[i] = static_cast<double>(rand()) / RAND_MAX;
            B[i] = static_cast<double>(rand()) / RAND_MAX;
            C[i] = 0.0;
        }

        // Benchmark OpenBLAS
        double openblas_time = 0.0;
        for (int trial = 0; trial < NTRIAL; ++trial) {
            auto start = std::chrono::high_resolution_clock::now();
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        n, n, n, alpha, A, n, B, n, beta, C, n);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            openblas_time += elapsed.count();
        }
        openblas_time /= NTRIAL;
        double openblas_gflops = (2.0 * n * n * n) / (openblas_time * 1e9);

        // Benchmark CUBLAS
        double cublas_time = 0.0;
        cublasHandle_t handle;
        cublasCreate(&handle);

        double *d_A, *d_B, *d_C;
        cudaMalloc(&d_A, size * sizeof(double));
        cudaMalloc(&d_B, size * sizeof(double));
        cudaMalloc(&d_C, size * sizeof(double));

        cublasSetMatrix(n, n, sizeof(double), A, n, d_A, n);
        cublasSetMatrix(n, n, sizeof(double), B, n, d_B, n);
        cublasSetMatrix(n, n, sizeof(double), C, n, d_C, n);

        // Warm-up
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    n, n, n, &alpha, d_A, n, d_B, n, &beta, d_C, n);

        for (int trial = 0; trial < NTRIAL; ++trial) {
            cudaEvent_t start, stop;
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
            cudaEventRecord(start);
            cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                        n, n, n, &alpha, d_A, n, d_B, n, &beta, d_C, n);
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            float milliseconds = 0;
            cudaEventElapsedTime(&milliseconds, start, stop);
            cublas_time += milliseconds / 1e3;
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
        }
        cublas_time /= NTRIAL;
        double cublas_gflops = (2.0 * n * n * n) / (cublas_time * 1e9);

        // Cleanup
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
        cublasDestroy(handle);

        delete[] A;
        delete[] B;
        delete[] C;

        // Write results for this n to file
        outfile << "n = " << n << "\n";
        outfile << "  OpenBLAS: " << openblas_gflops << " GFLOPs\n";
        outfile << "  CUBLAS  : " << cublas_gflops << " GFLOPs\n";
        outfile << "----------------------------\n";
    }

    outfile.close();
    return 0;
}