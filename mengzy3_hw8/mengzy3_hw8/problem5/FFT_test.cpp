#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <fftw3.h>
#include <cufft.h>
#include <cuda_runtime.h>

#define CUDA_CHECK(err) \
    do { \
        cudaError_t err_ = (err); \
        if (err_ != cudaSuccess) { \
            std::cerr << "CUDA error " << err_ << " at " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(err_) << std::endl; \
            std::exit(1); \
        } \
    } while (0)

#define CUFFT_CHECK(err) \
    do { \
        cufftResult err_ = (err); \
        if (err_ != CUFFT_SUCCESS) { \
            std::cerr << "CUFFT error " << err_ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::exit(1); \
        } \
    } while (0)

__global__ void multiply_wave_numbers_kernel(cufftDoubleComplex *fourier,
                                            cufftDoubleComplex *Ax,
                                            cufftDoubleComplex *Ay,
                                            cufftDoubleComplex *Az,
                                            double *kx, double *ky, double *kz,
                                            int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (i >= n || j >= n || k >= n) return;
    
    int idx = i + j * n + k * n * n;
    double re = fourier[idx].x;
    double im = fourier[idx].y;
    
    // Multiply by i*kx (pure imaginary)
    Ax[idx].x = -im * kx[i];
    Ax[idx].y =  re * kx[i];
    
    // Multiply by i*ky
    Ay[idx].x = -im * ky[j];
    Ay[idx].y =  re * ky[j];
    
    // Multiply by i*kz
    Az[idx].x = -im * kz[k];
    Az[idx].y =  re * kz[k];
}

void compute_wave_numbers(int n, std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz) {
    kx.resize(n);
    ky.resize(n);
    kz.resize(n);
    
    for (int i = 0; i < n; ++i) {
        double k = (i <= n/2) ? 2.0 * M_PI * i / n : 2.0 * M_PI * (i - n) / n;
        kx[i] = k;
        ky[i] = k;
        kz[i] = k;
    }
}

double measure_fftw(int n, int ntrial) {
    int total_points = n * n * n;
    std::vector<double> kx, ky, kz;
    compute_wave_numbers(n, kx, ky, kz);
    
    // Allocate memory
    fftw_complex *in = fftw_alloc_complex(total_points);
    fftw_complex *fourier = fftw_alloc_complex(total_points);
    fftw_complex *grad_x = fftw_alloc_complex(total_points);
    fftw_complex *grad_y = fftw_alloc_complex(total_points);
    fftw_complex *grad_z = fftw_alloc_complex(total_points);
    
    // Initialize input (plane wave)
    for (int i = 0; i < total_points; ++i) {
        in[i][0] = 1.0; // Real part
        in[i][1] = 0.0; // Imaginary part
    }
    
    // Create FFTW plans
    fftw_plan plan_forward = fftw_plan_dft_3d(n, n, n, in, fourier, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan plan_backward = fftw_plan_dft_3d(n, n, n, grad_x, grad_x, FFTW_BACKWARD, FFTW_MEASURE);
    
    double total_time = 0.0;
    
    for (int trial = 0; trial < ntrial; ++trial) {
        // Reset input (FFTW overwrites arrays)
        for (int i = 0; i < total_points; ++i) {
            in[i][0] = 1.0;
            in[i][1] = 0.0;
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Forward FFT
        fftw_execute_dft(plan_forward, in, fourier);
        
        // Multiplication in Fourier space
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    int idx = i + j * n + k * n * n;
                    double re = fourier[idx][0];
                    double im = fourier[idx][1];
                    
                    grad_x[idx][0] = -im * kx[i];
                    grad_x[idx][1] =  re * kx[i];
                    
                    grad_y[idx][0] = -im * ky[j];
                    grad_y[idx][1] =  re * ky[j];
                    
                    grad_z[idx][0] = -im * kz[k];
                    grad_z[idx][1] =  re * kz[k];
                }
            }
        }
        
        // Backward FFTs
        fftw_execute_dft(plan_backward, grad_x, grad_x);
        fftw_execute_dft(plan_backward, grad_y, grad_y);
        fftw_execute_dft(plan_backward, grad_z, grad_z);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        total_time += diff.count();
    }
    
    // Cleanup
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(fourier);
    fftw_free(grad_x);
    fftw_free(grad_y);
    fftw_free(grad_z);
    
    return total_time / ntrial;
}

double measure_cufft(int n, int ntrial) {
    int total_points = n * n * n;
    std::vector<double> h_kx, h_ky, h_kz;
    compute_wave_numbers(n, h_kx, h_ky, h_kz);
    
    // Allocate device memory
    cufftDoubleComplex *d_in, *d_fourier, *d_Ax, *d_Ay, *d_Az;
    double *d_kx, *d_ky, *d_kz;
    
    CUDA_CHECK(cudaMalloc(&d_in, sizeof(cufftDoubleComplex) * total_points));
    CUDA_CHECK(cudaMalloc(&d_fourier, sizeof(cufftDoubleComplex) * total_points));
    CUDA_CHECK(cudaMalloc(&d_Ax, sizeof(cufftDoubleComplex) * total_points));
    CUDA_CHECK(cudaMalloc(&d_Ay, sizeof(cufftDoubleComplex) * total_points));
    CUDA_CHECK(cudaMalloc(&d_Az, sizeof(cufftDoubleComplex) * total_points));
    CUDA_CHECK(cudaMalloc(&d_kx, sizeof(double) * n));
    CUDA_CHECK(cudaMalloc(&d_ky, sizeof(double) * n));
    CUDA_CHECK(cudaMalloc(&d_kz, sizeof(double) * n));
    
    // Copy wave numbers to device
    CUDA_CHECK(cudaMemcpy(d_kx, h_kx.data(), sizeof(double) * n, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ky, h_ky.data(), sizeof(double) * n, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_kz, h_kz.data(), sizeof(double) * n, cudaMemcpyHostToDevice));
    
    // Initialize input (plane wave)
    cufftDoubleComplex *h_in = new cufftDoubleComplex[total_points];
    for (int i = 0; i < total_points; ++i) {
        h_in[i].x = 1.0;
        h_in[i].y = 0.0;
    }
    CUDA_CHECK(cudaMemcpy(d_in, h_in, sizeof(cufftDoubleComplex) * total_points, cudaMemcpyHostToDevice));
    delete[] h_in;
    
    // Create CUFFT plans
    cufftHandle plan_forward, plan_inverse;
    CUFFT_CHECK(cufftPlan3d(&plan_forward, n, n, n, CUFFT_Z2Z));
    CUFFT_CHECK(cufftPlan3d(&plan_inverse, n, n, n, CUFFT_Z2Z));
    
    // Configure kernel launch
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (n + blockDim.x - 1) / blockDim.x,
        (n + blockDim.y - 1) / blockDim.y,
        (n + blockDim.z - 1) / blockDim.z
    );
    
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    float total_time = 0.0f;
    
    for (int trial = 0; trial < ntrial; ++trial) {
        CUDA_CHECK(cudaEventRecord(start));
        
        // Forward FFT
        CUFFT_CHECK(cufftExecZ2Z(plan_forward, d_in, d_fourier, CUFFT_FORWARD));
        
        // Multiplication kernel
        multiply_wave_numbers_kernel<<<gridDim, blockDim>>>(
            d_fourier, d_Ax, d_Ay, d_Az, d_kx, d_ky, d_kz, n
        );
        CUDA_CHECK(cudaGetLastError());
        
        // Backward FFTs
        CUFFT_CHECK(cufftExecZ2Z(plan_inverse, d_Ax, d_Ax, CUFFT_INVERSE));
        CUFFT_CHECK(cufftExecZ2Z(plan_inverse, d_Ay, d_Ay, CUFFT_INVERSE));
        CUFFT_CHECK(cufftExecZ2Z(plan_inverse, d_Az, d_Az, CUFFT_INVERSE));
        
        CUDA_CHECK(cudaEventRecord(stop));
        CUDA_CHECK(cudaEventSynchronize(stop));
        
        float milliseconds = 0;
        CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
        total_time += milliseconds / 1000.0f; // Convert to seconds
    }
    
    // Cleanup
    CUFFT_CHECK(cufftDestroy(plan_forward));
    CUFFT_CHECK(cufftDestroy(plan_inverse));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    CUDA_CHECK(cudaFree(d_in));
    CUDA_CHECK(cudaFree(d_fourier));
    CUDA_CHECK(cudaFree(d_Ax));
    CUDA_CHECK(cudaFree(d_Ay));
    CUDA_CHECK(cudaFree(d_Az));
    CUDA_CHECK(cudaFree(d_kx));
    CUDA_CHECK(cudaFree(d_ky));
    CUDA_CHECK(cudaFree(d_kz));
    
    return total_time / ntrial;
}

double compute_flops(int n, double time) {
    double N = n * n * n;
    double log2n = log2(n);
    // Operation count: 60*N*log2(N) + 6*N
    double total_ops = 60 * N * (3 * log2n) + 6 * N;
    return total_ops / time; // FLOPS
}

int main() {
    const int ntrial = 5;
    std::vector<int> n_values = {16, 32, 64, 128, 256};
    std::ofstream outfile("performance.csv");
    
    outfile << "n,FFTW_FLOPS,CUFFT_FLOPS\n";
    
    for (int n : n_values) {
        std::cout << "Running n = " << n << std::endl;
        
        double fftw_time = measure_fftw(n, ntrial);
        double cufft_time = measure_cufft(n, ntrial);
        
        double fftw_flops = compute_flops(n, fftw_time);
        double cufft_flops = compute_flops(n, cufft_time);
        
        outfile << n << "," << fftw_flops << "," << cufft_flops << "\n";
    }
    
    outfile.close();
    std::cout << "Results saved to performance.csv" << std::endl;
    
    return 0;
}