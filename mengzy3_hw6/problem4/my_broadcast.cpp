#include <mpi.h>
#include <iostream>
#include <vector>
#include <chrono>
#include "my_broadcast.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int root = 0;
    std::vector<size_t> message_sizes = {8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456};

    for (size_t size_bytes : message_sizes) {
        int count = size_bytes / sizeof(int);
        std::vector<int> data(count, 0);

        if (rank == root) {
            for (int i = 0; i < count; ++i) data[i] = i;
        }

        // Measure my_broadcast
        MPI_Barrier(MPI_COMM_WORLD);
        auto start = std::chrono::high_resolution_clock::now();
        my_broadcast(data.data(), count, root, MPI_COMM_WORLD);
        auto end = std::chrono::high_resolution_clock::now();
        double time_my = std::chrono::duration<double>(end - start).count();

        // Measure MPI_Bcast
        MPI_Barrier(MPI_COMM_WORLD);
        start = std::chrono::high_resolution_clock::now();
        MPI_Bcast(data.data(), count, MPI_INT, root, MPI_COMM_WORLD);
        end = std::chrono::high_resolution_clock::now();
        double time_mpi = std::chrono::duration<double>(end - start).count();

        if (rank == root) {
            double bandwidth_my = (size_bytes * (size - 1)) / time_my;
            double bandwidth_mpi = size_bytes / time_mpi;
            std::cout << size_bytes << " " << bandwidth_my << " " << bandwidth_mpi << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}