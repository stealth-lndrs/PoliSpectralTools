# SpectralTools.jl – User Guide

This guide explains how to install SpectralTools, summarizes the major modules,
and documents each solver so that users can quickly evaluate boundary-value and
time-dependent PDE problems via spectral collocation.

## Installation and Quick Start

```julia
julia> import Pkg; Pkg.add(path=\"/path/to/SpectralTools\") # until the package is registered
julia> using SpectralTools
# inspect submodules
julia> SpectralTools.Collocation
```

To develop locally:

```julia
$ git clone <repo-url> SpectralTools
$ cd SpectralTools
$ julia --project=.
julia> using Pkg
julia> Pkg.instantiate()
```

Running the test suite:

```julia
julia> include(\"test/runtests.jl\")
```

## Package Layout

| Module | Purpose |
| --- | --- |
| `SpectralTools.Generic` | Barycentric interpolation, generalized differentiation matrices, polynomial tools |
| `SpectralTools.Chebyshev`, `SpectralTools.Legendre`, `SpectralTools.Fourier` | Basis-specific transforms, quadrature, differentiation |
| `SpectralTools.Collocation` | `SpectralGrid` type and `build_grid` helper for Chebyshev/Legendre Lobatto nodes |
| `SpectralTools.BoundaryConditions` | Normalization of 1D/2D boundary specifications (`Dirichlet`, `Neumann`, `Robin`) |
| `SpectralTools.BVP` | Linear and nonlinear boundary-value problem solvers |
| `SpectralTools.PDE` | Time-dependent solvers: diffusion (MOL+RK4), wave (leapfrog), Poisson (Kronecker) |

## Feature Reference

### Collocation and Boundary Setup

```julia
grid = SpectralTools.Collocation.build_grid(40; basis = :chebyshev, domain = (-1, 1))
bc = SpectralTools.BoundaryConditions.normalize_1d_bc((
    left = (:dirichlet, 1.0),
    right = (:neumann, (x, t) -> cos(t))
))
```

`build_grid(N; basis, domain)` returns nodes in the reference interval and
physical domain, along with differentiation matrices `D1`, `D2` already scaled.

### Linear Boundary-Value Problems

```julia
result = SpectralTools.solve_linear_bvp(a, b, c, rhs;
    N = 48,
    basis = :chebyshev,     # or :legendre
    domain = (-1.0, 1.0),
    bc = (left = (:dirichlet, 0.0), right = (:neumann, 0.0))
)
```

- `a`, `b`, `c`, `rhs`: accept numbers, vectors, or functions evaluated at grid nodes.
- Recommended when the operator is linear and smooth; use Chebyshev for uniform accuracy or Legendre to emphasize interior points.
- Returns `(x, u, grid)` for downstream diagnostics.

### Nonlinear Boundary-Value Problems

```julia
sol = SpectralTools.solve_nonlinear_bvp(g;
    dg_dy = (x, y, yp) -> ∂g∂y,
    dg_dyp = (x, y, yp) -> ∂g∂y′,
    initial_guess = x -> 0.0,
    maxiter = 20,
    tol = 1e-10,
    bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0))
)
```

Use for problems of the form `y'' = g(x, y, y')`. Analytical derivatives speed
up Newton; if omitted, finite differences are used. Inspect `sol.converged` and
`sol.iterations` before trusting the result.

### Diffusion (Method of Lines + RK4)

```julia
diff = SpectralTools.solve_diffusion_1d(u0, (t0, t1);
    diffusivity = κ,
    N = 48,
    dt = :auto,              # optional manual step
    basis = :chebyshev,
    bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0)),
    forcing = (x, t, u) -> f(x, t, u)
)
```

- `u0`: number, vector, or function providing the initial state.
- `forcing`: optional source term that may depend on `(x, t, u)`.
- Dirichlet/Neumann mixed boundaries are supported; the solver enforces BCs at
each Runge–Kutta substep.

### Wave Equation (Leapfrog)

```julia
wave = SpectralTools.solve_wave_1d(u0, v0, (t0, t1);
    c = 1.0,
    N = 60,
    dt = :auto,
    bc = (left = (:neumann, (x, t) -> flux(t)), right = (:dirichlet, 0.0))
)
```

Leapfrog integrates displacement and midstep velocity. Monitor energy
`0.5(‖v‖² + c²‖D₁u‖²)` to check stability; choose `dt` consistent with the CFL
restriction `dt ≤ 2/N` for Chebyshev grids.

### 2D Poisson (Kronecker Solve)

```julia
poisson = SpectralTools.solve_poisson_2d(
    (x, y) -> forcing(x, y);
    Nx = 32,
    Ny = 32,
    basis = :chebyshev,
    domainx = (-1, 1),
    domainy = (-1, 1),
    bc = (left = (:dirichlet, 0.0), right = (:dirichlet, 0.0),
          bottom = (:dirichlet, 0.0), top = (:dirichlet, 0.0))
)
```

Provides `x`, `y`, `u` arrays with BCs imposed and interior solved via
`kron(I, Dy²) + kron(Dx², I)`.

## Example and Test Progress

| # | Scenario | Script/Test | Status |
| - | --- | --- | --- |
| 1 | Chebyshev linear BVP | `examples/bvp_linear.jl`, `Usage-driven Tests` set 1 | ✅ |
| 2 | Diffusion decay (Dirichlet) | `examples/diffusion.jl`, test suite entry | ✅ |
| 3 | Legendre linear BVP | `examples/bvp_legendre.jl`, Legendre grid test | ✅ |
| 4 | Nonlinear BVP Newton | `examples/bvp_nonlinear.jl`, convergence test | ✅ |
| 5 | Wave mixed BCs | `examples/wave_mixed_bc.jl`, energy test | ✅ |
| 6 | Traveling pulse wave | _Pending script/test_ | ⬜ |
| 7 | Forced diffusion | _Pending script/test_ | ⬜ |
| 8 | Poisson square | _Pending script/test_ | ⬜ |
| 9 | Poisson rectangular domain | _Pending script/test_ | ⬜ |
| 10 | Mapping-enhanced solver | _Pending script/test_ | ⬜ |

Legend: ✅ ready/implemented, ⬜ needs implementation.

## Visual Preview

Example diffusion decay curve:

![Diffusion decay](assets/diffusion_decay.png)
![Diffusion decay animation](assets/diffusion_decay.gif)

Wave mode snapshot:

![Wave mode](assets/wave_mode.png)
![Wave animation](assets/wave_mode.gif)

Manufactured Poisson surface:

![Poisson solution](assets/poisson_surface.png)
![Poisson rotation](assets/poisson_surface.gif)

## Next Steps

1. Publish generated docs via GitHub Pages (`docs/web/index.html`).
2. Implement the remaining five scenarios/tests per `EXAMPLES_TEST_PLAN.md`.
3. Run `julia --project --code-coverage=user scripts/run_full_report.jl`
   to collect consolidated pass/fail summaries and coverage percentages.
4. Package registration: once tests/examples complete, tag a release and submit
   to the Julia General registry using `PkgDev.register`.
