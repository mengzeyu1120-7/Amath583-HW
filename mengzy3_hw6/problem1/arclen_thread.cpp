#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
#include <chrono>
#include <string>

using namespace std;

const double a = 1.0;
const double b = 6.0;
const double L_true = 35.0/8.0 + log(6.0); // Analytical solution: 35/8 + ln(6)

// Integrand function for arc length calculation
double integrand(double x) {
    double derivative = 1.0/x - x/4.0;
    return sqrt(1.0 + derivative * derivative);
}

// Thread function to compute partial Riemann sum
void compute_partial(int thread_id, long n, int num_threads, double dx, vector<double>& partials) {
    long chunk = n / num_threads;
    int remainder = n % num_threads;
    long start = thread_id * chunk + (thread_id < remainder ? thread_id : remainder);
    long end = start + chunk + (thread_id < remainder ? 1 : 0);

    double sum = 0.0;
    for (long i = start; i < end; ++i) {
        double x = a + i * dx;
        sum += integrand(x) * dx;
    }
    partials[thread_id] = sum;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <n> <num_threads>\n";
        return 1;
    }

    long n = stol(argv[1]);
    int num_threads = stoi(argv[2]);
    double dx = (b - a) / n;

    vector<double> partials(num_threads, 0.0);
    vector<thread> threads;

    auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(compute_partial, i, n, num_threads, dx, ref(partials));
    }

    for (auto& t : threads) {
        t.join();
    }

    double total = 0.0;
    for (double sum : partials) {
        total += sum;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    double error = abs(total - L_true);

    // Output CSV format: n, num_threads, time, error
    cout << n << "," << num_threads << "," << duration.count() << "," << error << endl;

    return 0;
}