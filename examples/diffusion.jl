#!/usr/bin/env julia
using SpectralTools

κ = 1.0
u_exact(x, t) = exp(-π^2 * t / 4) * sin(π * (x + 1) / 2)

u0(x) = u_exact(x, 0.0)

sol = solve_diffusion_1d(u0, (0.0, 0.1);
    diffusivity = κ,
    basis = :chebyshev,
    N = 48,
    dt = 2.5e-4,
    bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))

u_end_exact = u_exact.(sol.x, sol.t[end])
err = maximum(abs.(sol.u[:, end] .- u_end_exact))
println("Diffusion RK4 solver max error = $(err)")
