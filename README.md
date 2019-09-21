# DiffEq
Various differential equation problems and solvers.


Poisson_DST: Poisson equation on the unit square with homogeneous Dirichlet boundary conditions. Solved using DST/FFT. Comparison with builtin Matlab solver. Needs Poisson_DST_square.mat (included here) and PDE Toolkit.

Heat_Wave: solves some heat and wave equations using the Crank-Nicholson scheme. Also examines some stability issues and tests an ill-posed problem (reversing direction of time). 

Trapezoid_Extrapolation: eigenvalue problem, solved with trapezoidal method plus Richardson's extrapolation.

Simple_Trap_Extrap: a simpler application of trapezoidal + Richardon's extrapolation. 
