using Test
using SpectralTools
using LinearAlgebra: norm, Diagonal, dot
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

@testset "Usage-driven Tests" begin
    @testset "Chebyshev BVP residual" begin
        a(x) = -(1 + x)
        b(x) = zero(x)
        c(x) = zero(x)
        rhs(x) = sinpi(x)
        sol = solve_linear_bvp(a, b, c, rhs; N = 48,
            basis = :chebyshev,
            domain = (-1.0, 1.0),
            bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))
        grid = sol.grid
        operator = Diagonal(a.(grid.x)) * grid.D2
        residual = operator * sol.u .- rhs.(grid.x)
        @test maximum(abs.(residual[2:end-1])) < 1e-9
        @test abs(sol.u[1]) < 1e-12 && abs(sol.u[end]) < 1e-12
    end

    @testset "Legendre grid properties" begin
        grid = build_grid(18; basis = :legendre)
        sym_err = maximum(abs.(grid.x .+ reverse(grid.x)))
        @test sym_err < 1e-13
        ones_vec = ones(length(grid.x))
        @test norm(grid.D1 * ones_vec, Inf) < 1e-12
    end

    @testset "Nonlinear BVP convergence" begin
        g(x, y, yp) = sin(y)
        dgdy(x, y, yp) = cos(y)
        dg_dyp(x, y, yp) = zero(x)
        sol = solve_nonlinear_bvp(g;
            dg_dy = dgdy,
            dg_dyp = dg_dyp,
            N = 48,
            basis = :chebyshev,
            bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)),
            maxiter = 15,
            tol = 1e-10)
        grid = sol.grid
        residual = grid.D2 * sol.u .- sin.(sol.u)
        @test sol.converged
        @test sol.iterations <= 8
        @test norm(residual, Inf) < 1e-8
    end

    @testset "Diffusion analytic comparison" begin
        u_exact(x, t) = exp(-π^2 * t / 4) * sin(π * (x + 1) / 2)
        u0(x) = u_exact(x, 0.0)
        tspan = (0.0, 0.05)
        sol_coarse = solve_diffusion_1d(u0, tspan;
            diffusivity = 1.0,
            N = 40,
            dt = 1e-5,
            bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))
        exact_final = u_exact.(sol_coarse.x, sol_coarse.t[end])
        err_coarse = maximum(abs.(sol_coarse.u[:, end] .- exact_final))
        @test err_coarse < 5e-5

        sol_fine = solve_diffusion_1d(u0, tspan;
            diffusivity = 1.0,
            N = 40,
            dt = 5e-6,
            bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))
        exact_final_fine = u_exact.(sol_fine.x, sol_fine.t[end])
        err_fine = maximum(abs.(sol_fine.u[:, end] .- exact_final_fine))
        ratio = err_coarse / err_fine
        @test ratio > 1.05
    end

    @testset "Wave energy mixed BCs" begin
        c = 1.0
        u0(x) = sin(π * (x + 1) / 2)
        v0(x) = zero(x)
        flux(t) = cos(5t)
        bc = (left = (:neumann, (x, t) -> flux(t)), right = (:dirichlet, 0.0))
        sol = solve_wave_1d(u0, v0, (0.0, 0.2);
            N = 40, dt = 5e-4, c = c, bc = bc)
        grid = build_grid(40; basis = :chebyshev)
        energy(u_slice, v_half) = begin
            grad = grid.D1 * u_slice
            return 0.5 * (norm(v_half)^2 + c^2 * norm(grad)^2)
        end
        E0 = energy(sol.u[:, 1], sol.v[:, 1])
        Eend = energy(sol.u[:, end], sol.v[:, end])
        drift = abs(Eend - E0) / E0
        @test drift <= 0.2
        left_flux = dot(view(grid.D1, 1, :), sol.u[:, end])
        @test isapprox(left_flux, flux(sol.t[end]); atol = 5e-3)
    end
end
