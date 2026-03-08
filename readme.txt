Multigrid 2D

Solves the 2D Poisson equation −∇²u = f on [0,Lx]×[0,Ly] with Dirichlet boundary conditions using the multigrid V-cycle.

Test problem:
  f(x,y) = 2π²sin(πx)sin(πy)
  Exact solution: u(x,y) = sin(πx)sin(πy)

Smoothers:
- Jacobi,
- Damped Jacobi (ω = 4/5),
- Gauss-Seidel,
- Red-Black GS,
- SOR.

Build: g++ -O2 -o multigrid_2d multigrid_2d.cpp
