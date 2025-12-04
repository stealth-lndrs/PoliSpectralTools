#!/usr/bin/env julia
using SpectralTools

# Manufactured solution: y(x) = x^5 - 2x + 1 on [-1, 1]
exact(x) = x^5 - 2x + 1

coef_a(x) = x^2
coef_b(x) = x
coef_c(x) = -1
forcing(x) = 5x^5 - 1

sol = solve_linear_bvp(coef_a, coef_b, coef_c, forcing;
    basis = :chebyshev,
    domain = (-1.0, 1.0),
    N = 32,
    bc = (left = (:dirichlet, exact(-1.0)), right = (:dirichlet, exact(1.0)))
)

err = maximum(abs.(sol.u .- exact.(sol.x)))
println("Chebyshev linear BVP solved with max error = $(err)")
