#!/usr/bin/env julia
#
# Usage Example 4: Nonlinear BVP y'' = sin(y) with homogeneous Dirichlet BCs.
# Demonstrates solve_nonlinear_bvp with analytic derivatives and reports
# Newton convergence diagnostics.

using SpectralTools
using LinearAlgebra: norm
using Printf

g(x, y, yp) = sin(y)
dgdy(x, y, yp) = cos(y)
dg_dyp(x, y, yp) = zero(x)

result = solve_nonlinear_bvp(g;
    dg_dy = dgdy,
    dg_dyp = dg_dyp,
    N = 48,
    basis = :chebyshev,
    domain = (-1.0, 1.0),
    bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)),
    tol = 1e-10,
    maxiter = 15,
)

grid = result.grid
residual = grid.D2 * result.u .- g.(grid.x, result.u, grid.D1 * result.u)
@printf("Nonlinear BVP iterations = %d, residual norm = %.3e (converged = %s)\n",
        result.iterations, norm(residual, Inf), result.converged)
