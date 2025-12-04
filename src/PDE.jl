module PDE

using LinearAlgebra

using ..Collocation: build_grid
using ..BoundaryConditions: normalize_1d_bc, normalize_2d_bc, eval_bc_value, eval_bc_time_derivative

export solve_diffusion_1d, solve_wave_1d, solve_poisson_2d

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

const SideSpecTuple = Tuple{Symbol, Any}

function call_forcing(force, x, t, u_val)
    if force isa Number
        return float(force)
    elseif force isa Function
        try
            return force(x, t, u_val)
        catch err
            if err isa MethodError
                try
                    return force(x, t)
                catch err2
                    if err2 isa MethodError
                        try
                            return force(x)
                        catch err3
                            if err3 isa MethodError
                                try
                                    return force(t)
                                catch err4
                                    if err4 isa MethodError
                                        return force()
                                    else
                                        rethrow(err4)
                                    end
                                end
                            else
                                rethrow(err3)
                            end
                        end
                    else
                        rethrow(err2)
                    end
                end
            else
                rethrow(err)
            end
        end
    else
        error("Unsupported forcing type $(typeof(force)).")
    end
end

function forcing_vector(force, x, t, u)
    if force === nothing
        return zeros(eltype(u), length(u))
    else
        vals = similar(u)
        for i in eachindex(x)
            vals[i] = call_forcing(force, x[i], t, u[i])
        end
        return vals
    end
end

function initialize_state(data, x)
    if data isa Number
        return fill(float(data), length(x))
    elseif data isa AbstractVector
        length(data) == length(x) || error("State vector has incorrect length.")
        return vec(float.(data))
    elseif data isa Function
        return float.(data.(x))
    else
        error("Unsupported initial condition $(typeof(data)).")
    end
end

boundary_index(side::Symbol, n::Int) = side === :left ? 1 : n

function derivative_row(grid, side::Symbol)
    return side === :left ? grid.D1[1, :] : grid.D1[end, :]
end

function enforce_state!(u, grid, bc_spec, t)
    n = length(u)
    # Dirichlet updates
    if bc_spec.left.kind === :dirichlet
        u[1] = eval_bc_value(bc_spec.left, grid.x[1], t)
    end
    if bc_spec.right.kind === :dirichlet
        u[end] = eval_bc_value(bc_spec.right, grid.x[end], t)
    end
    # Solve Neumann boundaries (can be one or both)
    unknowns = Int[]
    specs = SideSpecTuple[]
    if bc_spec.left.kind === :neumann
        push!(unknowns, 1)
        push!(specs, (:left, bc_spec.left))
    end
    if bc_spec.right.kind === :neumann
        push!(unknowns, n)
        push!(specs, (:right, bc_spec.right))
    end
    isempty(unknowns) && return
    D1 = grid.D1
    known = setdiff(1:n, unknowns)
    M = zeros(eltype(u), length(unknowns), length(unknowns))
    rhs = zeros(eltype(u), length(unknowns))
    for (eq, (side, spec)) in enumerate(specs)
        row = side === :left ? D1[1, :] : D1[end, :]
        rhs_val = eval_bc_value(spec, side === :left ? grid.x[1] : grid.x[end], t)
        rhs[eq] = rhs_val - sum(row[j] * u[j] for j in known)
        for (col, idx) in enumerate(unknowns)
            M[eq, col] = row[idx]
        end
    end
    sol = M \ rhs
    for (idx, val) in zip(unknowns, sol)
        u[idx] = val
    end
end
function apply_derivative_constraints!(du, grid, bc_spec, t)
    n = length(du)
    if bc_spec.left.kind === :dirichlet
        du[1] = eval_bc_time_derivative(bc_spec.left, grid.x[1], t)
    end
    if bc_spec.right.kind === :dirichlet
        du[end] = eval_bc_time_derivative(bc_spec.right, grid.x[end], t)
    end

    unknowns = Int[]
    specs = SideSpecTuple[]
    if bc_spec.left.kind === :neumann
        push!(unknowns, 1)
        push!(specs, (:left, bc_spec.left))
    end
    if bc_spec.right.kind === :neumann
        push!(unknowns, n)
        push!(specs, (:right, bc_spec.right))
    end
    isempty(unknowns) && return
    D1 = grid.D1
    known = setdiff(1:n, unknowns)
    M = zeros(eltype(du), length(unknowns), length(unknowns))
    rhs = zeros(eltype(du), length(unknowns))
    for (eq, (side, spec)) in enumerate(specs)
        row = side === :left ? D1[1, :] : D1[end, :]
        rhs_val = eval_bc_time_derivative(spec, side === :left ? grid.x[1] : grid.x[end], t)
        rhs[eq] = rhs_val - sum(row[j] * du[j] for j in known)
        for (col, idx) in enumerate(unknowns)
            M[eq, col] = row[idx]
        end
    end
    sol = M \ rhs
    for (idx, val) in zip(unknowns, sol)
        du[idx] = val
    end
end

# ---------------------------------------------------------------------------
# Diffusion solver (RK4)
# ---------------------------------------------------------------------------

function solve_diffusion_1d(u0, tspan; basis::Symbol = :chebyshev, domain::Tuple = (-1.0, 1.0),
                            N::Int = 32, diffusivity::Real = 1.0, dt = nothing,
                            bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)),
                            forcing = nothing)
    grid = build_grid(N; basis = basis, domain = domain)
    state = initialize_state(u0, grid.x)
    bc_spec = normalize_1d_bc(bc)
    enforce_state!(state, grid, bc_spec, tspan[1])

    total_time = tspan[2] - tspan[1]
    if dt === nothing
        dt = total_time / max(2000, 10 * N^2)
    end
    times = collect(tspan[1]:dt:tspan[2])
    sol = zeros(eltype(state), length(state), length(times))
    sol[:, 1] .= state

    function rhs(u, t)
        vals = diffusivity .* (grid.D2 * u)
        vals .+= forcing_vector(forcing, grid.x, t, u)
        return vals
    end

    for k in 1:(length(times) - 1)
        t = times[k]
        u = sol[:, k]
        k1 = rhs(u, t)
        apply_derivative_constraints!(k1, grid, bc_spec, t)

        u2 = u .+ 0.5 * dt .* k1
        enforce_state!(u2, grid, bc_spec, t + 0.5 * dt)
        k2 = rhs(u2, t + 0.5 * dt)
        apply_derivative_constraints!(k2, grid, bc_spec, t + 0.5 * dt)

        u3 = u .+ 0.5 * dt .* k2
        enforce_state!(u3, grid, bc_spec, t + 0.5 * dt)
        k3 = rhs(u3, t + 0.5 * dt)
        apply_derivative_constraints!(k3, grid, bc_spec, t + 0.5 * dt)

        u4 = u .+ dt .* k3
        enforce_state!(u4, grid, bc_spec, t + dt)
        k4 = rhs(u4, t + dt)
        apply_derivative_constraints!(k4, grid, bc_spec, t + dt)

        new_state = u .+ (dt / 6.0) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        enforce_state!(new_state, grid, bc_spec, times[k + 1])
        sol[:, k + 1] .= new_state
    end

    return (x = grid.x, t = times, u = sol)
end

# ---------------------------------------------------------------------------
# Wave solver (Leapfrog)
# ---------------------------------------------------------------------------

function solve_wave_1d(u0, v0, tspan; basis::Symbol = :chebyshev, domain::Tuple = (-1.0, 1.0),
                       N::Int = 32, c::Real = 1.0, dt = nothing,
                       bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)),
                       forcing = nothing)
    grid = build_grid(N; basis = basis, domain = domain)
    disp = initialize_state(u0, grid.x)
    vel = initialize_state(v0, grid.x)
    bc_spec = normalize_1d_bc(bc)
    enforce_state!(disp, grid, bc_spec, tspan[1])

    total_time = tspan[2] - tspan[1]
    dt === nothing && (dt = total_time / max(400, 6N))
    times = collect(tspan[1]:dt:tspan[2])
    sol = zeros(eltype(disp), length(disp), length(times))
    vel_hist = zeros(eltype(disp), length(disp), length(times))
    sol[:, 1] .= disp
    vel_hist[:, 1] .= vel

    function acceleration(u, t)
        vals = (c^2) .* (grid.D2 * u)
        vals .+= forcing_vector(forcing, grid.x, t, u)
        return vals
    end

    acc0 = acceleration(disp, times[1])
    apply_derivative_constraints!(acc0, grid, bc_spec, times[1])
    v_half = vel .+ 0.5 * dt .* acc0
    apply_derivative_constraints!(v_half, grid, bc_spec, times[1])

    for k in 1:(length(times) - 1)
        t = times[k]
        u_new = sol[:, k] .+ dt .* v_half
        enforce_state!(u_new, grid, bc_spec, times[k + 1])
        sol[:, k + 1] .= u_new

        acc = acceleration(u_new, times[k + 1])
        apply_derivative_constraints!(acc, grid, bc_spec, times[k + 1])
        v_half .+= dt .* acc
        apply_derivative_constraints!(v_half, grid, bc_spec, times[k + 1])
        vel_hist[:, k + 1] .= v_half
    end

    return (x = grid.x, t = times, u = sol, v = vel_hist)
end

# ---------------------------------------------------------------------------
# 2D Poisson solver
# ---------------------------------------------------------------------------

function call_rhs(f, x, y)
    if f isa Number
        return float(f)
    elseif f isa Function
        return f(x, y)
    elseif f isa AbstractMatrix
        error("Provide forcing as a function or scalar for Poisson solver.")
    else
        error("Unsupported forcing type $(typeof(f)).")
    end
end

function solve_poisson_2d(f; Nx::Int = 32, Ny::Int = 32, basis::Symbol = :chebyshev,
                           domainx::Tuple = (-1.0, 1.0), domainy::Tuple = (-1.0, 1.0),
                           bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0),
                                 bottom = (:dirichlet, 0.0), top = (:dirichlet, 0.0)))
    gridx = build_grid(Nx; basis = basis, domain = domainx)
    gridy = build_grid(Ny; basis = basis, domain = domainy)
    bc_spec = normalize_2d_bc(bc)

    ix = 2:(Nx - 1)
    iy = 2:(Ny - 1)
    length(ix) > 0 && length(iy) > 0 || error("Poisson solver requires at least 3 nodes per dimension.")

    X = gridx.x[ix]
    Y = gridy.x[iy]
    rhs = zeros(length(ix), length(iy))
    for (i, xval) in enumerate(X)
        for (j, yval) in enumerate(Y)
            rhs[i, j] = call_rhs(f, xval, yval)
        end
    end

    Dx2_full = gridx.D2
    Dy2_full = gridy.D2
    Dx2 = Dx2_full[ix, ix]
    Dy2 = Dy2_full[iy, iy]

    left_full = [eval_bc_value(bc_spec.left, (gridx.x[1], gridy.x[j]), nothing) for j in 1:Ny]
    right_full = [eval_bc_value(bc_spec.right, (gridx.x[end], gridy.x[j]), nothing) for j in 1:Ny]
    bottom_full = [eval_bc_value(bc_spec.bottom, (gridx.x[i], gridy.x[1]), nothing) for i in 1:Nx]
    top_full = [eval_bc_value(bc_spec.top, (gridx.x[i], gridy.x[end]), nothing) for i in 1:Nx]

    left_vals = left_full[iy]
    right_vals = right_full[iy]
    bottom_vals = bottom_full[ix]
    top_vals = top_full[ix]

    left_col = Dx2_full[ix, 1]
    right_col = Dx2_full[ix, end]
    bottom_col = Dy2_full[iy, 1]
    top_col = Dy2_full[iy, end]

    for (j, val) in enumerate(left_vals)
        rhs[:, j] .-= left_col .* val
    end
    for (j, val) in enumerate(right_vals)
        rhs[:, j] .-= right_col .* val
    end
    for (i, val) in enumerate(bottom_vals)
        for (j, coeff) in enumerate(bottom_col)
            rhs[i, j] -= coeff * val
        end
    end
    for (i, val) in enumerate(top_vals)
        for (j, coeff) in enumerate(top_col)
            rhs[i, j] -= coeff * val
        end
    end

    Ix = I(length(ix))
    Iy = I(length(iy))
    L = kron(Iy, Dx2) + kron(Dy2, Ix)
    sol_vec = L \ vec(rhs)
    U_inner = reshape(sol_vec, length(ix), length(iy))

    U = zeros(Nx, Ny)
    U[ix, iy] .= U_inner
    for j in 1:Ny
        U[1, j] = left_full[j]
        U[end, j] = right_full[j]
    end
    for i in 1:Nx
        U[i, 1] = bottom_full[i]
        U[i, end] = top_full[i]
    end

    return (x = gridx.x, y = gridy.x, u = U)
end

end
