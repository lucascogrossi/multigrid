#include <iostream>
#include <vector>
#include <cmath>

struct Grid {
    int n;                 // numero de intervalos
    double h;              // tamanho do intervalo
    std::vector<double> u; // solucao
    std::vector<double> f; // termo fonte

    // construtor
    Grid(int n, double L) : n(n), h(L/n), u(n+1, 0.0), f(n+1, 0.0) {}
};

double residual(Grid& grid) {
    double res = 0.0;
    for (int i = 1; i < grid.n; i++) {
        double Au_i = (-grid.u[i-1] + 2*grid.u[i] - grid.u[i+1]) / (grid.h*grid.h);
        double r_i = grid.f[i] - Au_i;
        res += r_i * r_i;
    }
    return sqrt(res);
}

void jacobi(Grid& grid) {
    std::vector<double> u_new(grid.n + 1, 0.0);

    for (int i = 1; i < grid.n; i++) {
        u_new[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
    }
    grid.u = u_new;
}

void jacobi_amortecido(Grid& grid) {
    std::vector<double> u_new(grid.n + 1, 0.0);
    double u_jacobi;
    double omega = 2.0/3.0;

    for (int i = 1; i < grid.n; i++) {
        u_jacobi = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
        u_new[i] = grid.u[i] + omega * (u_jacobi - grid.u[i]);
    }
    grid.u = u_new;
}

void gauss_seidel(Grid& grid) {
    for (int i = 1; i < grid.n; i++) {
        grid.u[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
    }
}

// Gauss-Seidel Sobrerelaxado
void sor(Grid& grid, double omega) {
    for (int i = 1; i < grid.n; i++) {
        double u_gs = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
        grid.u[i] = grid.u[i] + omega * (u_gs - grid.u[i]);
    }
}

void gauss_seidel_rb(Grid& grid) {
    // vermelho: indices impares
    for (int i = 1; i < grid.n; i += 2) {
        grid.u[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
    }
    // preto: indices pares
    for (int i = 2; i < grid.n; i += 2) {
        grid.u[i] = (grid.u[i-1] + grid.u[i+1]) / 2.0 + (grid.h*grid.h / 2.0) * grid.f[i];
    }
}

std::vector<double> restrict_grid(const std::vector<double>& r, int  n) {
    // grid grosso tem n/2 intervalos e n/2 - 1 pontos interiores
    int n_coarse = n / 2;
    std::vector<double> r_coarse(n_coarse + 1, 0.0);
    
    // copia pontos pares (injecao)
    for (int j = 1; j < n_coarse; j++) {
        r_coarse[j] = r[2*j];
    }
    
    return r_coarse;
}

std::vector<double> prolongate(const std::vector<double>& e_coarse, int n_coarse) {
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

void solve_coarse(Grid& coarse) {
    std::fill(coarse.u.begin(), coarse.u.end(), 0.0);
    for (int k = 0; k < 1000; k++) {
        gauss_seidel(coarse);
    }
}

void bigrid_cycle(Grid& grid) {
    // 1. pre-suavizacao
    for (int k = 0; k < 2; k++)
        jacobi_amortecido(grid);

    // 2. calcula residuo no grid fino
    std::vector<double> r(grid.n + 1, 0.0);
    for (int i = 1; i < grid.n; i++)
        r[i] = grid.f[i] - (-grid.u[i-1] + 2*grid.u[i] - grid.u[i+1]) / (grid.h*grid.h);

    // 3. restricao: leva residuo para grid grosso
    int n_coarse = grid.n / 2;
    std::vector<double> r_coarse = restrict_grid(r, grid.n);

    // 4. resolve no grid grosso
    Grid coarse_grid(n_coarse, 1.0);
    coarse_grid.f = r_coarse;
    solve_coarse(coarse_grid);

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
    // cria grid com 8 intervalos em [0, 1]
    Grid grid(128, 1.0);

    // preenche f
    for (int i = 1; i < grid.n; i++) {
        double x = i * grid.h;
        grid.f[i] = 1.0; // −uxx = 1
                    // a solucao analitica exata eh:
                    // u(x) = (x * (1 - x)) / 2
    }

    // std::cout << "h = " << h << std::endl;
    // std::cout << "pontos interiores: " << n - 1 << std::endl;

    for (int k = 0; k < 10; k++) {
        // jacobi(u, f, n, h);
        // gauss_seidel(u, f, n, h);
        // jacobi_amortecido(u, f, n, h);
        // sor(u, f, n, h, 1.9);
        // gauss_seidel_rb(u, f, n, h);
        // std::cout << "residuo: " << residual(u, f, n, h) << std::endl;
        bigrid_cycle(grid);
        std::cout << "residuo: " << residual(grid) << std::endl;
    }
    std::cout << "residuo final: " << residual(grid) << std::endl;

    std::cout << "\nSolucao aproximada:" << std::endl;
    for (int i = 0; i <= grid.n; i++) {
        double x = i * grid.h;
        std::cout << "u(" << x << ") = " << grid.u[i] << std::endl;
    }

    return 0;


}