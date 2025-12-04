using Test
using SpectralTools
using LinearAlgebra: norm
using SpectralTools: Generic
using SpectralTools.BVP
using SpectralTools.PDE
using SpectralTools.Collocation

@testset "Generic Utilities" begin
    @testset "Barycentric interpolation" begin
        x = collect(range(-1.0, 1.0; length = 5))
        y = x .^ 3 .- 2x .+ 1
        vals, P = Generic.Bary_Interp(x, y, [0.0, 0.5])
        @test size(P) == (2, length(x))
        @test vals[1] ≈ 1.0 atol = 1e-12
        @test vals[2] ≈ (0.5^3 - 2*0.5 + 1) atol = 1e-12
    end

    @testset "Differentiation matrix" begin
        grid = build_grid(6; basis = :chebyshev)
        D = Generic.Generalized_Diff_Mat(grid.x)
        ones_vec = ones(length(grid.x))
        @test norm(D * ones_vec, Inf) < 1e-12
    end
end

@testset "BVP Solvers" begin
    exact(x) = x^4
    a(x) = one(x)
    b(x) = zero(x)
    c(x) = zero(x)
    rhs(x) = 12x^2
    bvp = solve_linear_bvp(a, b, c, rhs; N = 28,
        bc = (left = (:dirichlet, exact(-1.0)), right = (:dirichlet, exact(1.0))))
    @test maximum(abs.(bvp.u .- exact.(bvp.x))) < 1e-8

    g(x, y, yp) = exp(y)
    dgdy(x, y, yp) = exp(y)
    sol = solve_nonlinear_bvp(g; dg_dy = dgdy, N = 40,
        bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))
    grid = sol.grid
    residual = grid.D2 * sol.u - exp.(sol.u)
    @test sol.converged
    @test norm(residual, Inf) < 5e-7
end

@testset "PDE Solvers" begin
    u_exact(x, t) = exp(-π^2 * t / 4) * sin(π * (x + 1) / 2)
    u0(x) = u_exact(x, 0.0)
    diff_sol = solve_diffusion_1d(u0, (0.0, 0.02); N = 36, dt = 2e-5,
        bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))
    final_exact = u_exact.(diff_sol.x, diff_sol.t[end])
    diff_err = maximum(abs.(diff_sol.u[:, end] .- final_exact))
    @test diff_err < 5e-5

    wave_exact(x, t) = cos(π * t / 2) * sin(π * (x + 1) / 2)
    wave_sol = solve_wave_1d(x -> wave_exact(x, 0.0), x -> 0.0,
        (0.0, 0.02); N = 36, dt = 3e-3,
        bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))
    wave_err = maximum(abs.(wave_sol.u[:, end] .- wave_exact.(wave_sol.x, wave_sol.t[end])))
    @test wave_err < 1e-3

    u2d(x, y) = sin(π * (x + 1) / 2) * sin(π * (y + 1) / 2)
    forcing(x, y) = - (π^2 / 2) * u2d(x, y)
    poisson = solve_poisson_2d(forcing; Nx = 18, Ny = 18)
    U_exact = [u2d(x, y) for x in poisson.x, y in poisson.y]
    @test maximum(abs.(poisson.u .- U_exact)) < 5e-4
end
