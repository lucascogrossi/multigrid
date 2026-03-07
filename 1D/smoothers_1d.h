#ifndef SMOOTHERS_H
#define SMOOTHERS_H

#include <vector>
#include "grid.h"

inline void jacobi(Grid& grid) {
    std::vector<double> u_new(grid.n + 1, 0.0);

    for (int i = 1; i < grid.n; i++) {
        u_new[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
    }
    grid.u = u_new;
}

inline void jacobi_amortecido(Grid& grid) {
    std::vector<double> u_new(grid.n + 1, 0.0);
    double u_jacobi;
    double omega = 2.0/3.0;

    for (int i = 1; i < grid.n; i++) {
        u_jacobi = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
        u_new[i] = grid.u[i] + omega * (u_jacobi - grid.u[i]);
    }
    grid.u = u_new;
}

inline void gauss_seidel(Grid& grid) {
    for (int i = 1; i < grid.n; i++) {
        grid.u[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
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