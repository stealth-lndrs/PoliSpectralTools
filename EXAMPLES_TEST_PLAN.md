# Usage Examples and Test Plan

This guide enumerates ten showcase examples and ten regression tests that demonstrate the
current SpectralTools capabilities. Each contributor can pick any combination of usage examples
and tests to implement independently. The expectations are simply that every selected item gains a
script under `examples/` and a matching automated test in `test/` (when applicable).

## Basic Usage Examples

### 1. Chebyshev Linear BVP (variable diffusivity)
- **Problem**: Solve `-(1 + x) y''(x) = sin(πx)` on `x ∈ [-1, 1]` with `y(-1) = y(1) = 0`.
- **Tools**: `solve_linear_bvp` with `basis = :chebyshev`, `N ∈ {32, 48}`.
- **Outputs**: Plot of numerical vs. exact solution `y(x) = (sin(πx) + πx cos(πx)) / (π^2 (1 + x))` (derived from integrating twice) and residual norm.
- **Hypothesis**: Spectral collocation yields machine-precision accuracy for smooth data; expect `L∞` error < `1e-9`.

### 2. Diffusion Decay (Dirichlet)
- **Problem**: Heat equation `u_t = u_xx` on `[-1, 1]`, `t ∈ [0, 0.1]`, Dirichlet zero boundaries, initial condition `u(x,0) = sin(π (x + 1)/2)`.
- **Tools**: `solve_diffusion_1d` with `N = 40`, `dt = 2e-4`, `basis = :chebyshev`.
- **Outputs**: compare numerical final state to analytic `u(x, t) = e^{-π^2 t / 4} sin(π (x + 1)/2)`; record error vs. time and confirm stability when reducing/increasing `dt`.
- **Hypothesis**: RK4 time integration preserves the analytic decay, with global error dominated by time step (`O(dt^4)`).

### 3. Legendre Linear BVP
- **Problem**: Same PDE as Usage 1 but solved with `basis = :legendre` to assess node placement effects.
- **Tools**: `solve_linear_bvp` using Legendre Lobatto grid, `N = 36`, plus `build_grid` introspection for node distribution.
- **Outputs**: error curves comparing Chebyshev vs. Legendre for identical `N`; histogram of node spacing illustrating more interior clustering for Chebyshev.
- **Hypothesis**: Legendre nodes yield comparable accuracy but may need slightly larger `N` for endpoint-sensitive problems.

### 4. Nonlinear BVP with Newton Damping
- **Problem**: `y'' = sin(y)` on `[-1, 1]`, with `y(-1) = y(1) = 0`.
- **Tools**: `solve_nonlinear_bvp` with `basis = :chebyshev`, `N = 48`, analytic derivatives `dg_dy = cos(y)`, `dg_dyp = 0`, initial guess `y0(x) = 0`.
- **Outputs**: iteration log (residual norms per Newton step), plot of converged solution, optional comparison with finite-difference solution.
- **Hypothesis**: Newton converges quadratically after first few iterations; expect < 8 iterations to reach tolerance `1e-10`.

### 5. Wave Propagation with Mixed BCs
- **Problem**: `u_tt = c^2 u_xx` on `[-1, 1]`, with left boundary Neumann `u_x(-1, t) = g(t) = cos(5 t)`, right boundary Dirichlet `u(1, t) = 0`, initial displacement `sin(π(x+1)/2)`, initial velocity zero, `t ∈ [0, 2]`, `c = 1`.
- **Tools**: `solve_wave_1d` with `basis = :chebyshev`, `N = 50`, `dt = 1e-3`, forcing term `forcing(x, t, u) = 0`.
- **Outputs**: time series of boundary values vs. imposed data, energy diagnostic `E(t)` verifying ≤ 1% drift, animation-ready snapshots of `u(x, t)`.
- **Hypothesis**: Mixed BC enforcement keeps Neumann flux accurate while satisfying Dirichlet clamp, demonstrating solver flexibility.

### 6. Wave Propagation with Traveling Pulse
- **Problem**: Initialize `u(x, 0) = exp(-50(x + 0.5)^2)`, `u_t(x, 0) = -c * u_x(x, 0)` to create a right-traveling Gaussian pulse on `[-1, 1]` with periodic re-entry suppressed by Dirichlet boundaries (`u(±1, t) = 0`), `c = 1`, `t ∈ [0, 1]`.
- **Tools**: `solve_wave_1d` with `N = 80`, `dt = 5e-4`, optional sponge forcing `forcing(x, t, u) = -0.5 * χ(x) * u` where `χ` damps near boundaries.
- **Outputs**: snapshot grid showing pulse propagation, CFL study varying `dt` to illustrate stability, energy metric with/without sponge.
- **Hypothesis**: Leapfrog remains stable when `dt ≤ 2/N`; adding sponge reduces reflections while slightly dissipating energy.

### 7. Diffusion with Source Forcing
- **Problem**: Solve `u_t = 0.1 u_xx + f(x, t)` with `f(x, t) = e^{-t} cos(πx)`, Dirichlet zero BCs, `u(x, 0) = 0`, `t ∈ [0, 0.5]`.
- **Tools**: `solve_diffusion_1d` with `diffusivity = 0.1`, `N = 30`, `dt = 1e-3`, forcing callback using the optional third argument `(x, t, u)`.
- **Outputs**: compare with numerically generated reference (e.g., high-resolution FD solver) and verify integral of `u` matches expected source injection.
- **Hypothesis**: Source term couples smoothly into RK4; integral of solution should track `∫₀^t ∫ f dx dt` minus diffusion losses.

### 8. Poisson on Square (Manufactured Solution)
- **Problem**: Solve `∇²u = f(x, y)` on `[-1, 1]^2` with Dirichlet BC matching `u_exact(x, y) = sin(π (x + 1)/2) sin(π (y + 1)/2)`. Forcing `f = - (π^2 / 2) u_exact`.
- **Tools**: `solve_poisson_2d` with `Nx = Ny = 24`, `basis = :chebyshev`.
- **Outputs**: heat map of numerical solution, max-error metric, and discussion of how Kronecker operator solves compare to separable solutions.
- **Hypothesis**: Expect spectral convergence; with 24 nodes, error should be < `1e-6`.

### 9. Poisson on Rectangular Domain
- **Problem**: Domain `x ∈ [-1, 2]`, `y ∈ [0, 1]`, Dirichlet BC from `u_exact(x, y) = sin(π (x + 1)/3) sinh(π y / 1)` (requires adjusting to satisfy Laplace). Choose `f = -Δ u_exact`.
- **Tools**: `solve_poisson_2d` with `Nx = 30`, `Ny = 20`, `basis = :chebyshev`, custom `domainx`, `domainy`.
- **Outputs**: demonstrate anisotropic node spacing, error vs. `Nx/Ny` sweeps, confirm domain scaling in differentiation matrices.
- **Hypothesis**: Differing domain lengths require careful scaling; solver should adapt automatically via `build_grid`.

### 10. Mapping Prototype (Optional Extension)
- **Problem**: If `Mapping.jl` is added, illustrate transforming `ξ ∈ [-1, 1]` to `x = sin(π ξ / 2)` for boundary-layer clustering. Apply to linear BVP and compare error vs. standard grid.
- **Tools**: hypothetical `SpectralTools.Mapping.simple_map` applied before `build_grid`.
- **Outputs**: node plots, error improvements near boundaries.
- **Hypothesis**: Conformal mappings reduce error for localized features.

## Julia Package Tests

### 1. Chebyshev BVP Residual Check
- **Objective**: Ensure Usage Example 1 achieves `L∞` error `< 1e-9`.
- **Test Steps**: Run `solve_linear_bvp` with `N = 40`; compute `Ly - f` via stored matrices to verify residual and confirm boundary values match exact data.
- **Acceptance Criteria**: `maximum(abs.(residual)) < 1e-9`, `u[1] ≈ 0`, `u[end] ≈ 0`.

### 2. Legendre Grid Consistency
- **Objective**: Validate `build_grid` with `basis = :legendre` returns symmetric nodes and `D1 * ones ≈ 0`.
- **Test Steps**: Build grid with `N = 18`; check `grid.x[k] ≈ -grid.x[end - k + 1]`, confirm `norm(grid.D1 * ones(N), Inf) < 1e-12`.
- **Acceptance Criteria**: Symmetry holds within `1e-13`, derivative matrix annihilates constants to machine precision.

### 3. Nonlinear BVP Newton Convergence
- **Objective**: Reproduce Usage Example 4 and ensure Newton converges in ≤ 8 iterations with residual < `1e-8`.
- **Test Steps**: Track `sol.iterations` and compute `norm(grid.D2 * sol.u - sin.(sol.u), Inf)`; optionally perturb initial guess to `0.1 * sin(π (x + 1) / 2)` to verify robustness.
- **Acceptance Criteria**: `sol.converged` true, `sol.iterations ≤ 8`, residual satisfies tolerance.

### 4. Diffusion Analytic Comparison
- **Objective**: Use Usage Example 2 parameters to confirm `L∞` error at `t = 0.1` below `5e-5`.
- **Test Steps**: Evaluate analytic solution at final time, compute error, and rerun with halved `dt` to confirm fourth-order convergence (error ratio ≈ 16).
- **Acceptance Criteria**: Baseline error < `5e-5`, refined run error ratio between 14–18.

### 5. Wave Energy Balance with Mixed BCs
- **Objective**: For Usage Example 5, verify discrete energy drift stays below `0.5%` over the simulation.
- **Test Steps**: Compute `E_k = 0.5 * (norm(v_half)^2 + c^2 * norm(D1 * u)^2)` at each stored time; compare start vs. end energies.
- **Acceptance Criteria**: `abs(E_end - E_start) / E_start ≤ 5e-3` and boundary constraints satisfied to `1e-8`.

### 6. Traveling Pulse CFL Study
- **Objective**: Ensure Usage Example 6 respects CFL and detects instability when `dt` exceeds limit.
- **Test Steps**: Run with `dt = 5e-4` (stable) and `dt = 1e-2` (expected unstable). Monitor max displacement growth and flag if stable run deviates < 1% from reference while unstable run diverges.
- **Acceptance Criteria**: Stable run `max|u|` remains bounded (< 2), unstable run exhibits > 10× amplification within 100 steps.

### 7. Diffusion Forcing Consistency
- **Objective**: Confirm solver honors forcing term from Usage Example 7 by comparing numerical integral of `u` against expected energy input.
- **Test Steps**: Compute `∑ u(x_i, t_end) Δx` and compare with trapezoidal integration of source minus diffusion losses estimated via `∑ D2 u`.
- **Acceptance Criteria**: Discrepancy < 2% when `dt = 1e-3`; halving `dt` reduces discrepancy by ≈ 4×.

### 8. Poisson Square Accuracy
- **Objective**: Validate Usage Example 8 achieves spectral accuracy.
- **Test Steps**: Solve with `Nx = Ny = 24`, compute `max|u - u_exact|`, and repeat with `Nx = Ny = 32` to confirm error reduction by at least 10×.
- **Acceptance Criteria**: Base error < `1e-6`, refined error ratio ≥ 8.

### 9. Poisson Rectangular Scaling
- **Objective**: Ensure rescaled domains (Usage Example 9) produce consistent accuracy along both axes.
- **Test Steps**: Compare second-derivative matrices to confirm correct scaling (`D2` factors `(2/(b - a))^2`). Solve for two aspect ratios and check max error below `5e-5`.
- **Acceptance Criteria**: `D2` scaling match analytic factors within `1e-12`; solution error < `5e-5`.

### 10. Mapping Impact (if implemented)
- **Objective**: Quantify benefits of optional mapping utilities.
- **Test Steps**: Run BVP from Usage Example 10 with and without mapping, record boundary-layer error, ensure mapping reduces error by ≥ 3×.
- **Acceptance Criteria**: Documented improvement; fallback to skip if mapping module absent (test should pass trivially with warning).

Each contributor should document assumptions, share diagnostic plots, and ensure their tests run
under `julia --project=. include("test/runtests.jl")` before opening a PR.
