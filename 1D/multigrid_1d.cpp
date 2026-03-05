/*
TODO:
- injecao -> full weighting
- thomas algorithm on solve_coarse
*/

#include <iostream>
#include <vector>
#include <cmath>

#include "smoothers.h"
#include "grid.h"
#include "multigrid_utils.h"

void v_cycle(Grid& grid) {

    // Condicao de parada: grid com 1 ponto interior
    if (grid.n == 2) {
        solve_coarse(grid);
        return;
    }

    // 1. pre-suavizacao
    for (int k = 0; k < 2; k++)
        jacobi_amortecido(grid);

    // 2. calcula residuo no grid fino
    std::vector<double> r = compute_residual(grid);

    // 3. restricao: leva residuo para grid grosso
    int n_coarse = grid.n / 2;
    std::vector<double> r_coarse = restrict(r, grid.n);

    // 4. resolve no grid grosso
    Grid coarse_grid(n_coarse, grid.L);
    coarse_grid.f = r_coarse;
    v_cycle(coarse_grid);

    // 5. prolongamento: leva correcao de volta para grid fino
    std::vector<double> e_fine = prolongate(coarse_grid.u, n_coarse);

    // 6. corrige solucao
    for (int i = 1; i < grid.n; i++)
        grid.u[i] += e_fine[i];

    // 7. pos-suavizacao
    for (int k = 0; k < 2; k++)
        jacobi_amortecido(grid);
}

int main(void) {
    // cria grid com 128 intervalos em [0, 1]
    Grid grid(128, 1.0);

    // preenche f
    for (int i = 1; i < grid.n; i++) {
        double x = i * grid.h;
        grid.f[i] = 1.0; // −uxx = 1
                         // a solucao analitica exata eh:
                         // u(x) = (x * (1 - x)) / 2
    }

    for (int k = 0; k < 10; k++) {
        v_cycle(grid);
        std::cout << "residuo " << "k = " << k << " " << residual_norm(grid) << std::endl;
    }
    std::cout << "residuo final: " << residual_norm(grid) << std::endl;

    /*
    std::cout << "\nSolucao aproximada:" << std::endl;
    for (int i = 0; i <= grid.n; i++) {
        double x = i * grid.h;
        std::cout << "u(" << x << ") = " << grid.u[i] << std::endl;
    }
    */

    return 0;
}