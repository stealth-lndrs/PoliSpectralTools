#!/usr/bin/env julia
#
# Usage Example 3: Linear BVP solved on Legendre Lobatto grids.
# Problem matches Example 1, allowing direct comparison against a
# Chebyshev reference interpolated via barycentric interpolation.

using SpectralTools
using SpectralTools: Generic
using Printf

a(x) = -(1 + x)
b(x) = zero(x)
c(x) = zero(x)
rhs(x) = sinpi(x)
bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0))

legendre_sol = solve_linear_bvp(a, b, c, rhs;
    N = 40,
    basis = :legendre,
    domain = (-1.0, 1.0),
    bc = bc,
)

cheb_reference = solve_linear_bvp(a, b, c, rhs;
    N = 96,
    basis = :chebyshev,
    domain = (-1.0, 1.0),
    bc = bc,
)

interp_vals, _ = Generic.Bary_Interp(cheb_reference.x, cheb_reference.u, legendre_sol.x)
err = maximum(abs.(legendre_sol.u .- interp_vals))

@printf("Legendre linear BVP: max difference vs. Chebyshev reference = %.3e\n", err)
