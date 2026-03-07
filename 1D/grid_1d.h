#ifndef GRID_H
#define GRID_H

#include <vector>

struct Grid {
    int n;                 // numero de intervalos
    double h;              // tamanho do intervalo
    double L;              // comprimento do dominio [0, L]
    std::vector<double> u; // solucao
    std::vector<double> f; // termo fonte

    // construtor
    Grid(int n, double L) : n(n), L(L), h(L/n), u(n+1, 0.0), f(n+1, 0.0) {}
};

#endif