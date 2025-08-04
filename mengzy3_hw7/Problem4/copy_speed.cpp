#include <iostream>
#include <fstream>
#include <iomanip>
#include <cuda_runtime.h>

#define GB (1024ULL * 1024 * 1024)

void checkCudaError(cudaError_t err, const char* msg) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error: " << msg << " - " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main() {
    cudaEvent_t start, stop;
    checkCudaError(cudaEventCreate(&start), "Failed to create start event");
    checkCudaError(cudaEventCreate(&stop), "Failed to create stop event");

    // Open output file
    std::ofstream output_file("copy_speed_results.txt");
    if (!output_file.is_open()) {
        std::cerr << "Failed to open output file!" << std::endl;
        return EXIT_FAILURE;
    }

    // Write header to file
    output_file << "BufferSize(Bytes) H2D_Bandwidth(Bytes/s) D2H_Bandwidth(Bytes/s)" << std::endl;

    for (size_t size = 1; size <= 2 * GB; size *= 2) {
        char *host_ptr = nullptr;
        char *device_ptr = nullptr;

        // Allocate pinned host memory
        cudaError_t err = cudaMallocHost(&host_ptr, size);
        if (err != cudaSuccess) {
            std::cerr << "cudaMallocHost failed for size " << size << ": " << cudaGetErrorString(err) << std::endl;
            continue;
        }

        // Allocate device memory
        err = cudaMalloc(&device_ptr, size);
        if (err != cudaSuccess) {
            std::cerr << "cudaMalloc failed for size " << size << ": " << cudaGetErrorString(err) << std::endl;
            cudaFreeHost(host_ptr);
            continue;
        }

        // Measure H2D transfer
        float h2d_time;
        checkCudaError(cudaEventRecord(start), "H2D start event record failed");
        cudaMemcpy(device_ptr, host_ptr, size, cudaMemcpyHostToDevice);
        checkCudaError(cudaEventRecord(stop), "H2D stop event record failed");
        cudaEventSynchronize(stop);
        checkCudaError(cudaEventElapsedTime(&h2d_time, start, stop), "H2D elapsed time failed");
        float h2d_bandwidth = size / (h2d_time / 1000);

        // Measure D2H transfer
        float d2h_time;
        checkCudaError(cudaEventRecord(start), "D2H start event record failed");
        cudaMemcpy(host_ptr, device_ptr, size, cudaMemcpyDeviceToHost);
        checkCudaError(cudaEventRecord(stop), "D2H stop event record failed");
        cudaEventSynchronize(stop);
        checkCudaError(cudaEventElapsedTime(&d2h_time, start, stop), "D2H elapsed time failed");
        float d2h_bandwidth = size / (d2h_time / 1000);

        // Write formatted results to file
        output_file << size << " " 
                    << std::fixed << std::setprecision(2) << h2d_bandwidth << " " 
                    << d2h_bandwidth << std::endl;

        // Free resources
        cudaFree(device_ptr);
        cudaFreeHost(host_ptr);
    }

    // Cleanup
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    output_file.close();
    return 0;
}