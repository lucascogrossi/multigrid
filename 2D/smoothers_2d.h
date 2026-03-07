#ifndef SMOOTHERS_H
#define SMOOTHERS_H

#include <vector>
#include "grid_2d.h"

inline void jacobi(Grid2D& grid) {
    std::vector<double> u_new((grid.nx+1) * (grid.ny+1), 0.0);

    for (int i = 1; i < grid.nx; i++) {
        for (int j = 1; j < grid.ny; j++) {
            u_new[grid.idx(i, j)] = (grid.u[grid.idx(i-1, j)] +
                                     grid.u[grid.idx(i+1, j)] +
                                     grid.u[grid.idx(i, j-1)] +
                                     grid.u[grid.idx(i, j+1)]) / 4.0 +
                                     (grid.hx*grid.hx / 4.0) * grid.f[grid.idx(i,j)];
        }
    }
    grid.u = u_new;
}

inline void jacobi_amortecido(Grid2D& grid) {
    std::vector<double> u_new((grid.nx + 1) * (grid.ny+1), 0.0);
    double u_jacobi;
    double omega = 4.0/5.0; // valor otimo para suavizacao

    for (int i = 1; i < grid.nx; i++) {
        for (int j = 1; j < grid.ny; j++) {
            u_jacobi = (grid.u[grid.idx(i-1, j)] +
                        grid.u[grid.idx(i+1, j)] +
                        grid.u[grid.idx(i, j-1)] +
                        grid.u[grid.idx(i, j+1)]) / 4.0 +
                        (grid.hx*grid.hx / 4.0) * grid.f[grid.idx(i,j)];

            u_new[grid.idx(i, j)] = grid.u[grid.idx(i, j)] + omega * (u_jacobi - grid.u[grid.idx(i, j)]);
        }
    }
    grid.u = u_new;
}

inline void gauss_seidel(Grid2D& grid) {
    for (int i = 1; i < grid.nx; i++) {
        for (int j = 1; j < grid.ny; j++) {
            grid.u[grid.idx(i, j)] = (grid.u[grid.idx(i-1, j)] +
                                     grid.u[grid.idx(i+1, j)] +
                                     grid.u[grid.idx(i, j-1)] +
                                     grid.u[grid.idx(i, j+1)]) / 4.0 +
                                     (grid.hx*grid.hx / 4.0) * grid.f[grid.idx(i,j)];
        }
    }
}

// Gauss-Seidel Sobrerelaxado
inline void sor(Grid& grid, double omega) {
    for (int i = 1; i < grid.n; i++) {
        double u_gs = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
        grid.u[i] = grid.u[i] + omega * (u_gs - grid.u[i]);
    }
}

inline void gauss_seidel_rb(Grid& grid) {
    // vermelho: indices impares
    for (int i = 1; i < grid.n; i += 2) {
        grid.u[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
    }
    // preto: indices pares
    for (int i = 2; i < grid.n; i += 2) {
        grid.u[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
    }
}

#endif