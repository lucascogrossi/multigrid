#ifndef MULTIGRID_UTILS_H
#define MULTIGRID_UTILS_H

#include <cmath>

#include "grid.h"
#include "smoothers.h"


inline std::vector<double> compute_residual(Grid& grid) {
    std::vector<double> r(grid.n + 1, 0.0);
    for (int i = 1; i < grid.n; i++)
        r[i] = grid.f[i] - (-grid.u[i-1] + 2*grid.u[i] - grid.u[i+1]) / (grid.h*grid.h);
    return r;
}

inline double residual_norm(Grid& grid) {
    std::vector<double> r = compute_residual(grid);
    double norm = 0.0;
    for (int i = 1; i < grid.n; i++)
        norm += r[i] * r[i];
    return sqrt(norm);
}

inline std::vector<double> restrict(const std::vector<double>& r, int  n) {
    // grid grosso tem n/2 intervalos e n/2 - 1 pontos interiores
    int n_coarse = n / 2;
    std::vector<double> r_coarse(n_coarse + 1, 0.0);
    
    // copia pontos pares (injecao)
    for (int j = 1; j < n_coarse; j++) {
        r_coarse[j] = r[2*j];
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