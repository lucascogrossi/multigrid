#include <iostream>
#include <vector>
#include <cmath>

double residual(const std::vector<double>& u, const std::vector<double>& f, int n, double h) {
    double res = 0.0;
    for (int i = 1; i < n; i++) {
        double Au_i = (-u[i-1] + 2*u[i] - u[i+1]) / (h*h);
        double r_i = f[i] - Au_i;
        res += r_i * r_i;
    }
    return sqrt(res);
}

void jacobi(std::vector<double>& u, const std::vector<double>& f, int n, double h) {
    std::vector<double> u_new(n + 1, 0.0);

    for (int i = 1; i < n; i++) {
        u_new[i] = (u[i-1] + u[i+1]) / 2.0 + (h*h / 2.0) * f[i];
    }
    u = u_new;
}

void jacobi_amortecido(std::vector<double>& u, const std::vector<double>& f, int n, double h) {
    std::vector<double> u_new(n + 1, 0.0);
    double u_jacobi;
    double omega = 2.0/3.0;

    for (int i = 1; i < n; i++) {
        u_jacobi = (u[i-1] + u[i+1]) / 2.0 + (h*h / 2.0) * f[i];
        u_new[i] = u[i] + omega * (u_jacobi - u[i]);
    }
    u = u_new;
}

void gauss_seidel(std::vector<double>& u, const std::vector<double>& f, int n, double h) {
    for (int i = 1; i < n; i++) {
        u[i] = (u[i-1] + u[i+1]) / 2.0 + (h*h / 2.0) * f[i];
    }
}

// Gauss-Seidel Sobrerelaxado
void sor(std::vector<double>& u, const std::vector<double>& f, int n, double h, double omega) {
    for (int i = 1; i < n; i++) {
        double u_gs = (u[i-1] + u[i+1]) / 2.0 + (h*h / 2.0) * f[i];
        u[i] = u[i] + omega * (u_gs - u[i]);
    }
}

void gauss_seidel_rb(std::vector<double>& u, const std::vector<double>& f, int n, double h) {
    // vermelho: indices impares
    for (int i = 1; i < n; i += 2) {
        u[i] = (u[i-1] + u[i+1]) / 2.0 + (h*h / 2.0) * f[i];
    }
    // preto: indices pares
    for (int i = 2; i < n; i += 2) {
        u[i] = (u[i-1] + u[i+1]) / 2.0 + (h*h / 2.0) * f[i];
    }
}

std::vector<double> restrict_grid(const std::vector<double>& r, int n) {
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

void solve_coarse(std::vector<double>& e, const std::vector<double>& r, int n_coarse, double h_coarse) {
    std::fill(e.begin(), e.end(), 0.0);
    for (int k = 0; k < 1000; k++) {
        gauss_seidel(e, r, n_coarse, h_coarse);
    }
}

void bigrid_cycle(std::vector<double>& u, const std::vector<double>& f, int n, double h) {
    // 1. pre-suavizacao
    for (int k = 0; k < 2; k++)
        jacobi_amortecido(u, f, n, h);

    // 2. calcula residuo no grid fino
    std::vector<double> r(n + 1, 0.0);
    for (int i = 1; i < n; i++)
        r[i] = f[i] - (-u[i-1] + 2*u[i] - u[i+1]) / (h*h);

    // 3. restricao: leva residuo para grid grosso
    int n_coarse = n / 2;
    double h_coarse = 2 * h;
    std::vector<double> r_coarse = restrict_grid(r, n);

    // 4. resolve no grid grosso
    std::vector<double> e_coarse(n_coarse + 1, 0.0);
    solve_coarse(e_coarse, r_coarse, n_coarse, h_coarse);

    // 5. prolongamento: leva correcao de volta para grid fino
    std::vector<double> e_fine = prolongate(e_coarse, n_coarse);

    // 6. corrige solucao
    for (int i = 1; i < n; i++)
        u[i] += e_fine[i];

    // 7. pos-suavizacao
    for (int k = 0; k < 2; k++)
        jacobi_amortecido(u, f, n, h);
}

int main(void) {
    int n = 8;   // numero de intervalos
    double L = 1.0; // dominio [0, L]
    double h = L / n; // distancia entre dois pontos consecutivos (tamanho do intervalo)

    std::vector<double> u(n + 1, 0.0);
    std::vector<double> f(n + 1, 0.0);

    // preenche f
    for (int i = 1; i < n; i++) {
        double x = i * h;
        f[i] = 1.0; // −uxx = 1
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
        bigrid_cycle(u, f, n, h);
        std::cout << "residuo: " << residual(u, f, n, h) << std::endl;
    }
    std::cout << "residuo final: " << residual(u, f, n, h) << std::endl;

    std::cout << "\nSolucao aproximada:" << std::endl;
    for (int i = 0; i <= n; i++) {
        double x = i * h;
        std::cout << "u(" << x << ") = " << u[i] << std::endl;
    }

    return 0;


}