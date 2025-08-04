#include <iostream>
#include <cmath>
#include <mpi.h>

const double L_true = 35.0/8.0 + std::log(6.0); 

double f_deriv(double x) {
    return 1.0 / x - 0.25 * x;
}

double compute_partial_sum(int start, int end, double dx, double a) {
    double sum = 0.0;
    for (int i = start; i < end; ++i) {
        double x = a + i * dx;
        double dfdx = f_deriv(x);
        sum += std::sqrt(1.0 + dfdx * dfdx) * dx;
    }
    return sum;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc != 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <n>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    int n = std::atoi(argv[1]);
    double a = 1.0;
    double b = 6.0;
    double dx = (b - a) / n;


    int quotient = n / size;
    int remainder = n % size;
    int start, end;
    
    if (rank < remainder) {
        start = rank * (quotient + 1);
        end = start + (quotient + 1);
    } else {
        start = remainder * (quotient + 1) + (rank - remainder) * quotient;
        end = start + quotient;
    }
    
    double start_time = MPI_Wtime();
    double partial_sum = compute_partial_sum(start, end, dx, a);
    double end_time = MPI_Wtime();
    
    double total_sum;
    MPI_Reduce(&partial_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        double error = L_true - total_sum;
        double log_error = std::log10(std::abs(error));
        std::cout << "n=" << n 
                  << ", L_numeric=" << total_sum
                  << ", error=" << error
                  << ", log10(error)=" << log_error
                  << ", Time=" << end_time - start_time << "s" 
                  << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}