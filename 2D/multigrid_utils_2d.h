#ifndef MULTIGRID_UTILS_H
#define MULTIGRID_UTILS_H

#include <cmath>

#include "grid_2d.h"
#include "smoothers_2d.h"


inline std::vector<double> compute_residual(Grid2D& grid) {
    std::vector<double> r((grid.nx+1) * (grid.ny+1), 0.0);
    for (int i = 1; i < grid.nx; i++) {
        for (int j = 1; j < grid.ny; j++) {
            double Au_ij = (-grid.u[grid.idx(i-1,j)]
                            -grid.u[grid.idx(i,j-1)]
                            +4*grid.u[grid.idx(i,j)]
                            -grid.u[grid.idx(i+1,j)]
                            -grid.u[grid.idx(i,j+1)]) / (grid.hx*grid.hx);
        r[grid.idx(i, j)] = grid.f[grid.idx(i, j)] - Au_ij;
        }
    }
    return r;
}

inline double residual_norm(Grid2D& grid) {
    std::vector<double> r = compute_residual(grid);
    double norm = 0.0;
    for (int i = 1; i < grid.nx; i++) {
        for (int j = 1; j < grid.ny; j++) {
            norm += r[grid.idx(i, j)] * r[grid.idx(i, j)];
        }
    }
    return sqrt(norm);
}

inline std::vector<double> restrict(const std::vector<double>& r, int nx, int ny) {
    // numero de intervalos do grid gross em cada direcao
    int nx_c = nx / 2;
    int ny_c = ny / 2;
    
    std::vector<double> r_coarse((nx_c+1) * (ny_c + 1), 0.0);
    
    // Injecao (full-weighting TODO)
    for (int i = 1; i < nx_c; i++) {
        for (int j = 1; j < ny_c; j++) {
            r_coarse[i * (ny_c + 1) + j] = r[2*i*(ny+1) + 2*j];
        }
    }
    return r_coarse;
}

inline std::vector<double> prolongate(const std::vector<double>& e_coarse, int n_coarse) {
    // grid fino tem n_coarse*2 intervalos
    int n_fine = n_coarse * 2;
    std::vector<double> e_fine(n_fine + 1, 0.0);

    for (int j = 0; j < n_coarse; j++) {
        // pontos pares (copia direto)
        e_fine[2*j] = e_coarse[j];
        // pontos impares  (media dos vizinhos; interpolacao linear entre 2 pontos)
        e_fine[2*j+1] = (e_coarse[j] + e_coarse[j+1]) / 2.0;
    }
    // copia o ultimo ponto
    e_fine[n_fine] = e_coarse[n_coarse];

    return e_fine;
}

inline void solve_coarse(Grid& coarse) {
    std::fill(coarse.u.begin(), coarse.u.end(), 0.0);
    for (int k = 0; k < 1000; k++) {
        gauss_seidel(coarse);
    }
}

#endif