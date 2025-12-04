module BVP

using LinearAlgebra

using ..Collocation: SpectralGrid, build_grid
using ..BoundaryConditions: normalize_1d_bc, eval_bc_value

export solve_linear_bvp, solve_nonlinear_bvp

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

function coefficient_vector(data, x::AbstractVector)
    if data isa Number
        return fill(float(data), length(x))
    elseif data isa AbstractVector
        length(data) == length(x) || error("Vector coefficient has wrong length.")
        return vec(float.(data))
    elseif data isa Function
        return float.(data.(x))
    else
        error("Unsupported coefficient type $(typeof(data)).")
    end
end

boundary_index(side::Symbol, n::Int) = side === :left ? 1 : n

function derivative_row(grid::SpectralGrid, side::Symbol)
    return side === :left ? grid.D1[1, :] : grid.D1[end, :]
end

function identity_row(T, n::Int, idx::Int)
    row = zeros(T, n)
    row[idx] = one(T)
    return row
end

function apply_linear_bc!(A, rhs, grid::SpectralGrid, bc)
    n = length(grid.x)
    for (side, spec) in pairs(bc)
        idx = boundary_index(side, n)
        xval = grid.x[idx]
        if spec.kind === :dirichlet
            A[idx, :] .= 0
            A[idx, idx] = 1
            rhs[idx] = eval_bc_value(spec, xval, nothing)
        elseif spec.kind === :neumann
            A[idx, :] .= derivative_row(grid, side)
            rhs[idx] = eval_bc_value(spec, xval, nothing)
        elseif spec.kind === :robin
            row = spec.alpha .* identity_row(eltype(A), size(A, 2), idx) .+
                  spec.beta .* derivative_row(grid, side)
            A[idx, :] .= row
            rhs[idx] = eval_bc_value(spec, xval, nothing)
        else
            error("Unknown BC type $(spec.kind)")
        end
    end
end

function impose_bc_on_residual!(J, F, u, grid::SpectralGrid, bc)
    n = length(u)
    for (side, spec) in pairs(bc)
        idx = boundary_index(side, n)
        xval = grid.x[idx]
        if spec.kind === :dirichlet
            J[idx, :] .= 0
            J[idx, idx] = 1
            F[idx] = u[idx] - eval_bc_value(spec, xval, nothing)
        elseif spec.kind === :neumann
            row = derivative_row(grid, side)
            J[idx, :] .= row
            F[idx] = dot(row, u) - eval_bc_value(spec, xval, nothing)
        elseif spec.kind === :robin
            row = spec.alpha .* identity_row(eltype(J), size(J, 2), idx) .+
                  spec.beta .* derivative_row(grid, side)
            J[idx, :] .= row
            F[idx] = spec.alpha * u[idx] + spec.beta * dot(derivative_row(grid, side), u) -
                     eval_bc_value(spec, xval, nothing)
        else
            error("Unsupported BC type $(spec.kind)")
        end
    end
end

# ---------------------------------------------------------------------------
# Linear solver
# ---------------------------------------------------------------------------

"""
    solve_linear_bvp(a, b, c, rhs; kwargs...)

Solve a linear second-order BVP of the form

    a(x) y''(x) + b(x) y'(x) + c(x) y(x) = rhs(x)

using spectral collocation.

Keyword arguments:
- `N`: number of Lobatto nodes (default 32).
- `basis`: `:chebyshev` or `:legendre`.
- `domain`: tuple `(a, b)`.
- `bc`: named tuple describing boundary data, e.g.
  `(left = (:dirichlet, 0.0), right = (:neumann, 1.0))`.
"""
function solve_linear_bvp(a, b, c, rhs; N::Int = 32, basis::Symbol = :chebyshev,
                          domain::Tuple = (-1.0, 1.0),
                          bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)))
    grid = build_grid(N; basis = basis, domain = domain)
    ax = coefficient_vector(a, grid.x)
    bx = coefficient_vector(b, grid.x)
    cx = coefficient_vector(c, grid.x)
    fv = coefficient_vector(rhs, grid.x)

    A = Diagonal(ax) * grid.D2 + Diagonal(bx) * grid.D1 + Diagonal(cx)
    bc_spec = normalize_1d_bc(bc)
    apply_linear_bc!(A, fv, grid, bc_spec)

    sol = A \ fv
    return (x = grid.x, u = sol, grid = grid)
end

# ---------------------------------------------------------------------------
# Nonlinear solver
# ---------------------------------------------------------------------------

function approximate_derivative(fun, x, y, yp, gvals; ε = 1e-7)
    forward = fun.(x, y .+ ε, yp)
    return (forward .- gvals) ./ ε
end

function approximate_derivative_yp(fun, x, y, yp, gvals; ε = 1e-7)
    forward = fun.(x, y, yp .+ ε)
    return (forward .- gvals) ./ ε
end

"""
    solve_nonlinear_bvp(g; kwargs...)

Solve nonlinear problems of the form `y'' = g(x, y, y')` via Newton iteration.

Keyword arguments:
- `dg_dy`, `dg_dyp`: optional analytic derivatives of `g`.
- `basis`, `N`, `domain`, `bc`: same as the linear solver.
- `initial_guess`: vector or function used as the starting iterate.
- `maxiter`, `tol`: Newton parameters.
"""
function solve_nonlinear_bvp(g; dg_dy = nothing, dg_dyp = nothing,
                              N::Int = 32, basis::Symbol = :chebyshev,
                              domain::Tuple = (-1.0, 1.0),
                              bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)),
                              initial_guess = nothing, maxiter::Int = 20, tol = 1e-10)
    grid = build_grid(N; basis = basis, domain = domain)
    u = initial_guess === nothing ? zeros(length(grid.x)) : coefficient_vector(initial_guess, grid.x)
    bc_spec = normalize_1d_bc(bc)

    converged = false
    nit = 0
    for k in 1:maxiter
        nit = k
        Dy = grid.D1 * u
        D2y = grid.D2 * u
        gvals = g.(grid.x, u, Dy)
        dgdy = dg_dy === nothing ? approximate_derivative(g, grid.x, u, Dy, gvals) : dg_dy.(grid.x, u, Dy)
        dgdyp = dg_dyp === nothing ? approximate_derivative_yp(g, grid.x, u, Dy, gvals) : dg_dyp.(grid.x, u, Dy)

        F = D2y .- gvals
        J = copy(grid.D2)
        J .-= Diagonal(dgdy)
        J .-= Diagonal(dgdyp) * grid.D1

        impose_bc_on_residual!(J, F, u, grid, bc_spec)
        δ = J \ (-F)
        u .+= δ
        if norm(F, Inf) < tol && norm(δ, Inf) < tol
            converged = true
            break
        end
    end

    return (x = grid.x, u = u, iterations = nit, converged = converged, grid = grid)
end

end
