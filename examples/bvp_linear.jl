#!/usr/bin/env julia
#
# Usage Example 1: Chebyshev linear BVP with variable diffusivity.
#
# Problem: -(1 + x) y''(x) = sin(πx), x ∈ [-1, 1], y(-1) = y(1) = 0.
# Goal: demonstrate solve_linear_bvp with Chebyshev Lobatto grid and
# verify the discrete residual of the manufactured operator.

using SpectralTools
using SpectralTools: Collocation
using LinearAlgebra
using Printf

a(x) = -(1 + x)
b(x) = zero(x)
c(x) = zero(x)
rhs(x) = sinpi(x)  # sin(πx)

params = (N = 48, basis = :chebyshev, domain = (-1.0, 1.0))
solution = solve_linear_bvp(a, b, c, rhs;
    N = params.N,
    basis = params.basis,
    domain = params.domain,
    bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)),
)

grid = solution.grid
A = Diagonal(a.(grid.x)) * grid.D2
residual = A * solution.u .- rhs.(grid.x)
inner_residual = residual[2:end-1]
@printf("Chebyshev linear BVP solved with ‖residual‖∞ = %.3e (N = %d)\n",
        maximum(abs, inner_residual), params.N)
