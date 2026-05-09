#include <iostream>
#include <vector>
#include <complex>
#include <pthread.h>
#include <cmath>
#include <fstream>

using namespace std;
typedef complex<double> dcomp;

constexpr double PI    = 3.141592653589793;
constexpr double EPS0  = 8.854e-12;
constexpr double OMEGA = 2.0 * PI * 14.1e6;
constexpr double K     = OMEGA / 299792458.0;

// Sandy Hook, CT — ITU-R P.527, New England glacial till
constexpr double SIGMA = 0.002;  // S/m
constexpr double EPS_R = 13.0;   // relative permittivity

// Gauss-Legendre points per sub-interval of the Sommerfeld integral
constexpr int N1     = 128;  // [0,      K    ]  propagating
constexpr int N2     = 256;  // [K,    100·K  ]  near evanescent
constexpr int N3     = 128;  // [100·K, 1000·K]  far evanescent
constexpr int N_QUAD = N1 + N2 + N3;

struct SommerfeldQuad {
    double lam   [N_QUAD];  // quadrature abscissas λ
    dcomp  u0    [N_QUAD];  // √(λ²−K²), branch Im(u₀)≥0
    dcomp  factor[N_QUAD];  // R_TM·λ³/u₀·w/(4π) — geometric kernel per point
    dcomp  S     [N_QUAD];  // dz·Σ Iₙ·exp(−u₀·zₙ) — current array factor
};

struct GridData {
    int    start_row, end_row, grid_size, N_segments;
    double z_min, z_max, y_min, y_max, dz_seg;
    const dcomp*          I_vec;
    double*               field_results;
    const SommerfeldQuad* sq;
};

// Square root with branch Re(u)≥0; tie-break Im(u)≥0 when Re(u)=0
static dcomp branch_sqrt(dcomp z) {
    dcomp u = sqrt(z);
    if (real(u) < 0.0 || (real(u) == 0.0 && imag(u) < 0.0)) u = -u;
    return u;
}

// Fill n Gauss-Legendre nodes x[] and weights w[] mapped to [a,b],
// starting at x[offset] and w[offset].
static void gl_fill(int n, double a, double b,
                    double* x, double* w, int offset) {
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double xi = cos(PI * (i + 0.75) / (n + 0.5));
        for (int iter = 0; iter < 100; ++iter) {
            double p0 = 1.0, p1 = xi;
            for (int j = 2; j <= n; ++j) {
                double p2 = ((2*j - 1)*xi*p1 - (j - 1)*p0) / j;
                p0 = p1; p1 = p2;
            }
            double dp  = n * (p0 - xi*p1) / (1.0 - xi*xi);
            double dxi = p1 / dp;
            xi -= dxi;
            if (fabs(dxi) < 1e-15) break;
        }
        double p0 = 1.0, p1 = xi;
        for (int j = 2; j <= n; ++j) {
            double p2 = ((2*j - 1)*xi*p1 - (j - 1)*p0) / j;
            p0 = p1; p1 = p2;
        }
        double dp = n * (p0 - xi*p1) / (1.0 - xi*xi);
        double wi = 2.0 / ((1.0 - xi*xi) * dp*dp);

        double mid = 0.5*(a + b), half = 0.5*(b - a);
        x[offset + i]         = mid - half*xi;  w[offset + i]         = half*wi;
        x[offset + n - 1 - i] = mid + half*xi;  w[offset + n - 1 - i] = half*wi;
    }
}

void* thread_func(void* arg) {
    auto d = static_cast<GridData*>(arg);
    const SommerfeldQuad* sq = d->sq;
    const double dz_grid = (d->z_max - d->z_min) / d->grid_size;
    const double dy_grid = (d->y_max - d->y_min) / d->grid_size;

    // Precompute running-product seeds for exp(-u0[q]*z_obs):
    //   decay_z0[q]   = exp(-u0[q] * z_min)   — value at first column
    //   decay_step[q] = exp(-u0[q] * dz_grid) — multiplier per column step
    // Avoids one complex exp() per pixel per quadrature point (replace with multiply).
    dcomp decay_z0  [N_QUAD];
    dcomp decay_step[N_QUAD];
    for (int q = 0; q < N_QUAD; ++q) {
        decay_z0  [q] = exp(-sq->u0[q] * d->z_min);
        decay_step[q] = exp(-sq->u0[q] * dz_grid);
    }

    // Per-row Sommerfeld cache: factor[q]*S[q]*J0(lam[q]*y_obs)
    // J0 depends on y_obs (row) but not z_obs (column), so cache per row.
    dcomp row_cache[N_QUAD];
    dcomp decay_cur[N_QUAD];

    for (int i = d->start_row; i < d->end_row; ++i) {
        const double y_obs = d->y_min + i * dy_grid;

        for (int q = 0; q < N_QUAD; ++q) {
            row_cache[q] = sq->factor[q] * sq->S[q] * j0(sq->lam[q] * y_obs);
            decay_cur[q] = decay_z0[q];
        }

        for (int j = 0; j < d->grid_size; ++j) {
            const double z_obs = d->z_min + j * dz_grid;

            if (z_obs <= 0.0) {
                d->field_results[i * d->grid_size + j] = 1e-30;
                // Still step the running products to stay in sync
                for (int q = 0; q < N_QUAD; ++q)
                    decay_cur[q] *= decay_step[q];
                continue;
            }

            // Direct free-space field from each antenna segment
            dcomp total_Ez = 0;
            for (int k = 0; k < d->N_segments; ++k) {
                const double z_src = k * d->dz_seg;
                const double dz    = z_obs - z_src;
                double R           = sqrt(dz*dz + y_obs*y_obs);
                if (R < 0.001) R   = 0.001;

                const double cos2  = dz*dz / (R*R);
                const double sin2  = y_obs*y_obs / (R*R);
                const dcomp  G     = exp(dcomp(0.0, -K*R)) / (4.0*PI*R);
                const dcomp  kern  = K*K*sin2
                                   + dcomp(3.0*cos2 - 1.0)*dcomp(1.0/(R*R), K/R);
                total_Ez += d->I_vec[k] * G * kern * d->dz_seg
                          / dcomp(0.0, OMEGA*EPS0);
            }

            // Sommerfeld ground-reflected field:
            //   Ez_refl = (1/jωε₀) · Σq factor[q]·S[q]·J0(λq·ρ)·exp(−u0q·z)
            dcomp Ez_refl = 0;
            for (int q = 0; q < N_QUAD; ++q)
                Ez_refl += row_cache[q] * decay_cur[q];
            Ez_refl /= dcomp(0.0, OMEGA*EPS0);

            d->field_results[i * d->grid_size + j] = abs(total_Ez + Ez_refl);

            // Advance running products to next column
            for (int q = 0; q < N_QUAD; ++q)
                decay_cur[q] *= decay_step[q];
        }
    }
    return nullptr;
}

int main() {
    constexpr int    GRID         = 1000;
    constexpr int    NUM_THREADS  = 8;
    constexpr double TOTAL_LENGTH = 5.55;
    constexpr int    N_SEG        = 500;
    constexpr double dz_seg       = TOTAL_LENGTH / N_SEG;

    // Sinusoidal current distribution: zero at tip, maximum near base
    vector<dcomp> I_vec(N_SEG);
    for (int j = 0; j < N_SEG; ++j) {
        const double z_src = j * dz_seg;
        I_vec[j] = dcomp(sin(K * (TOTAL_LENGTH - z_src)), 0.0);
    }

    // Ground complex permittivity and wavenumber
    const dcomp eps_c = dcomp(EPS_R, -SIGMA / (OMEGA * EPS0));
    const dcomp K1    = K * branch_sqrt(eps_c);

    // Build Sommerfeld quadrature table
    SommerfeldQuad sq;
    {
        double raw_lam[N_QUAD], raw_w[N_QUAD];
        gl_fill(N1, 0.0,        K,          raw_lam, raw_w, 0);
        gl_fill(N2, K,          100.0*K,    raw_lam, raw_w, N1);
        gl_fill(N3, 100.0*K,    1000.0*K,   raw_lam, raw_w, N1 + N2);

        for (int q = 0; q < N_QUAD; ++q) {
            const double lam = raw_lam[q];
            const dcomp  u0  = branch_sqrt(dcomp(lam*lam - K*K));
            const dcomp  u1  = branch_sqrt(dcomp(lam*lam) - K1*K1);
            const dcomp  Rtm = (eps_c*u0 - u1) / (eps_c*u0 + u1);

            sq.lam   [q] = lam;
            sq.u0    [q] = u0;
            // Absorb quadrature weight and 1/(4π) into factor
            sq.factor[q] = Rtm * (lam*lam*lam / u0) * raw_w[q] / (4.0*PI);

            // Current array factor: S = dz_seg · Σ Iₙ · exp(−u₀·zₙ)
            dcomp S = 0;
            for (int n = 0; n < N_SEG; ++n)
                S += I_vec[n] * exp(-u0 * (n * dz_seg));
            sq.S[q] = S * dz_seg;
        }
    }

    auto results = new double[GRID * GRID];
    pthread_t threads[NUM_THREADS];
    GridData  td[NUM_THREADS];
    const int rows_per_thread = GRID / NUM_THREADS;

    for (int i = 0; i < NUM_THREADS; ++i) {
        td[i].start_row     = i * rows_per_thread;
        td[i].end_row       = (i == NUM_THREADS - 1) ? GRID : (i+1)*rows_per_thread;
        td[i].grid_size     = GRID;
        td[i].N_segments    = N_SEG;
        td[i].z_min = -2.0; td[i].z_max = 8.0;
        td[i].y_min =  0.1; td[i].y_max = 10.1;
        td[i].dz_seg        = dz_seg;
        td[i].I_vec         = I_vec.data();
        td[i].field_results = results;
        td[i].sq            = &sq;
        pthread_create(&threads[i], nullptr, thread_func, &td[i]);
    }

    for (pthread_t thread : threads) pthread_join(thread, nullptr);

    ofstream outfile("field_map.bin", ios::out | ios::binary);
    outfile.write(reinterpret_cast<char*>(results), GRID * GRID * sizeof(double));
    outfile.close();

    delete[] results;
    return 0;
}