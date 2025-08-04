#ifndef MY_BROADCAST_HPP
#define MY_BROADCAST_HPP

#include <mpi.h>

template <typename T>
void my_broadcast(T* data, int count, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        for (int dest = 0; dest < size; ++dest) {
            if (dest != root) {
                MPI_Send(data, count * sizeof(T), MPI_BYTE, dest, 0, comm);
            }
        }
    } else {
        MPI_Recv(data, count * sizeof(T), MPI_BYTE, root, 0, comm, MPI_STATUS_IGNORE);
    }
}

#endif // MY_BROADCAST_HPP