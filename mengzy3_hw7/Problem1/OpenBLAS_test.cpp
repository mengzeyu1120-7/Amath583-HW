#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cblas.h>

const int NTRIAL = 3;
const int MIN_N = 2;
const int MAX_N = 4096;
const int STEP = 2;

double get_time() {
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main() {
    openblas_set_num_threads(1);
    
    std::ofstream csv("blas_results.csv");
    csv << "n,L1_MFLOPs,L2_MFLOPs,L3_MFLOPs\n";

    for (int n = MIN_N; n <= MAX_N; n *= STEP) {
        // Memory allocation
        std::vector<double> X(n), Y(n);
        std::vector<double> A(n*n), B(n*n), C(n*n);

        // Initialize data
        srand(time(NULL));
        for(auto& v : X) v = rand()/static_cast<double>(RAND_MAX);
        for(auto& v : Y) v = rand()/static_cast<double>(RAND_MAX);
        for(auto& v : A) v = rand()/static_cast<double>(RAND_MAX);
        for(auto& v : B) v = rand()/static_cast<double>(RAND_MAX);
        for(auto& v : C) v = rand()/static_cast<double>(RAND_MAX);

        // --- L1: daxpy ---
        double alpha = 1.0;
        double t_total = 0.0;
        for(int t=0; t<NTRIAL; ++t) {
            std::vector<double> Y_copy = Y;
            double t_start = get_time();
            cblas_daxpy(n, alpha, X.data(), 1, Y_copy.data(), 1);
            t_total += get_time() - t_start;
        }
        double flops_l1 = (2.0 * n * NTRIAL) / t_total;

        // --- L2: dgemv ---
        t_total = 0.0;
        for(int t=0; t<NTRIAL; ++t) {
            std::vector<double> Y_copy(n, 0);
            double t_start = get_time();
            cblas_dgemv(CblasColMajor, CblasNoTrans, 
                       n, n, alpha, A.data(), n, 
                       X.data(), 1, 0.0, Y_copy.data(), 1);
            t_total += get_time() - t_start;
        }
        double flops_l2 = (2.0 * n * n * NTRIAL) / t_total;

        // --- L3: dgemm ---
        t_total = 0.0;
        for(int t=0; t<NTRIAL; ++t) {
            std::vector<double> C(n*n, 0);
            double t_start = get_time();
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                       n, n, n, alpha, 
                       A.data(), n, 
                       B.data(), n, 
                       0.0, C.data(), n);
            t_total += get_time() - t_start;
        }
        double flops_l3 = (2.0 * n * n * n * NTRIAL) / t_total;

        // Write CSV
        csv << n << ","
            << flops_l1/1e6 << ","
            << flops_l2/1e6 << ","
            << flops_l3/1e6 << "\n";
    }
    
    return 0;
}