#include <iostream>
#include <vector>
#include <complex>
#include <pthread.h>
#include <cmath>
#include <fstream>

using namespace std;

typedef complex<double> dcomp;

const double PI = 3.141592653589793;
const double EPS0 = 8.854e-12;
const double OMEGA = 2.0 * PI * 14.1e6;
const double K = OMEGA / 299792458.0;

struct GridData {
    int start_row;
    int end_row;
    int grid_size;
    int N_segments;
    double z_min, z_max, y_min, y_max;
    double dz_seg;
    const dcomp* I_vec;
    double* field_results; // Store magnitude for plotting
};

// Core Field Kernel
double calc_field_magnitude(double z_obs, double y_obs, int N_seg, double dz_seg, const dcomp* I_vec) {
    dcomp total_Ez = 0;
    for (int j = 0; j < N_seg; ++j) {
        double z_src = j * dz_seg;
        double dz    = z_obs - z_src;
        double R = sqrt(dz * dz + y_obs * y_obs);
        if (R < 0.001) R = 0.001;

        double cos2 = (dz * dz) / (R * R);
        double sin2 = (y_obs * y_obs) / (R * R);

        dcomp G = exp(dcomp(0, -K * R)) / (4.0 * PI * R);

        // Full near-field Ez kernel: d²G/dz² + k²G
        // = G × [ k²sin²θ + (3cos²θ − 1)(1/R² + jk/R) ]
        dcomp kernel = K * K * sin2
                     + dcomp(3.0 * cos2 - 1.0) * dcomp(1.0 / (R * R), K / R);

        total_Ez += I_vec[j] * G * kernel * dz_seg / dcomp(0, OMEGA * EPS0);
    }
    return abs(total_Ez);
}

void* thread_func(void* arg) {
    GridData* d = (GridData*)arg;
    double dz_grid = (d->z_max - d->z_min) / d->grid_size;
    double dy_grid = (d->y_max - d->y_min) / d->grid_size;

    for (int i = d->start_row; i < d->end_row; ++i) {
        double y_obs = d->y_min + i * dy_grid;
        for (int j = 0; j < d->grid_size; ++j) {
            double z_obs = d->z_min + j * dz_grid;
            d->field_results[i * d->grid_size + j] = calc_field_magnitude(z_obs, y_obs, d->N_segments, d->dz_seg, d->I_vec);
        }
    }
    return NULL;
}

int main() {
    const int GRID = 1000;
    const int NUM_THREADS = 8;

    // Chameleon MPAS 2.0 Physical Specs
    const double TOTAL_LENGTH = 5.55; // Meters (Whip + Extension)
    const int N_SEG = 500;
    const double dz_seg = TOTAL_LENGTH / N_SEG;

    // Pre-calculate Sinusoidal Current Distribution
    vector<dcomp> I_vec(N_SEG);
    for (int j = 0; j < N_SEG; ++j) {
        double z_src = j * dz_seg;
        // Current is max at base (z=0) and 0 at the tip (z=L)
        double current_mag = sin(K * (TOTAL_LENGTH - z_src));
        I_vec[j] = dcomp(current_mag, 0.0);
    }

    double* results = new double[GRID * GRID];
    pthread_t threads[NUM_THREADS];
    GridData td[NUM_THREADS];
    int rows_per_thread = GRID / NUM_THREADS;

    for (int i = 0; i < NUM_THREADS; ++i) {
        td[i].start_row = i * rows_per_thread;
        td[i].end_row = (i == NUM_THREADS - 1) ? GRID : (i + 1) * rows_per_thread;
        td[i].grid_size = GRID;
        td[i].N_segments = N_SEG;

        // Adjust view to match the 5.55m physical height
        td[i].z_min = -2.0; td[i].z_max = 8.0;
        td[i].y_min = 0.1;  td[i].y_max = 10.1;

        td[i].dz_seg = dz_seg;
        td[i].I_vec = I_vec.data();
        td[i].field_results = results;
        pthread_create(&threads[i], NULL, thread_func, &td[i]);
    }

    for (int i = 0; i < NUM_THREADS; ++i) pthread_join(threads[i], NULL);

    // Export to Binary (Fastest for 1M doubles)
    ofstream outfile("field_map.bin", ios::out | ios::binary);
    outfile.write(reinterpret_cast<char*>(results), GRID * GRID * sizeof(double));
    outfile.close();

    delete[] results;
    return 0;
}