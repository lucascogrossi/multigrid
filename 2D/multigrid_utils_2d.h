#ifndef MULTIGRID_UTILS_2D_H
#define MULTIGRID_UTILS_2D_H
#include <cmath>

#include "grid_2d.h"
#include "smoothers_2d.h"


inline std::vector<double> compute_residual(Grid2D& grid) {
    std::vector<double> r((grid.nx+1) * (grid.ny+1), 0.0);
    double hx2 = grid.hx * grid.hx;
    double hy2 = grid.hy * grid.hy;

    for (int i = 1; i < grid.nx; i++) {
        for (int j = 1; j < grid.ny; j++) {
            double Au_ij = (-grid.u[grid.idx(i-1,j)] + 2*grid.u[grid.idx(i,j)] - grid.u[grid.idx(i+1,j)]) / hx2
                         + (-grid.u[grid.idx(i,j-1)] + 2*grid.u[grid.idx(i,j)] - grid.u[grid.idx(i,j+1)]) / hy2;
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

inline std::vector<double> prolongate(const std::vector<double>& e_coarse, int nx_c, int ny_c) {
    // grid fino tem n_coarse*2 intervalos
    int nx_f = nx_c * 2;
    int ny_f = ny_c * 2;
    std::vector<double> e_fine((nx_f+1) * (ny_f+1), 0.0);

    // Interpolacao linear
    for (int i = 0; i <= nx_c; i++) {
        for (int j = 0; j <= ny_c; j++) {
            // caso 1: copia direto
            e_fine[2*i*(ny_f+1) + 2*j] = e_coarse[i*(ny_c+1) + j];

            // caso 2: media horizontal
            if (j < ny_c)
                e_fine[2*i*(ny_f+1) + 2*j+1] = (e_coarse[i*(ny_c+1) + j] + e_coarse[i*(ny_c+1) + j+1]) / 2.0;
            
            // caso 3: media vertical
            if (i < nx_c)
                 e_fine[(2*i+1)*(ny_f+1) + 2*j] = (e_coarse[i*(ny_c+1) + j] + e_coarse[(i+1)*(ny_c+1) + j]) / 2.0;
            
            // caso 4: media dos 4 vizinhos
            if (i < nx_c && j < ny_c)
                e_fine[(2*i+1)*(ny_f+1) + 2*j+1] = (e_coarse[i*(ny_c+1) + j] + e_coarse[i*(ny_c+1) + j+1] +
                                                    e_coarse[(i+1)*(ny_c+1) + j] + e_coarse[(i+1)*(ny_c+1) + j+1]) / 4.0;
        }
    }

    return e_fine;
}

inline void solve_coarse(Grid2D& coarse) {
    std::fill(coarse.u.begin(), coarse.u.end(), 0.0);
    gauss_seidel(coarse);
}

#endif