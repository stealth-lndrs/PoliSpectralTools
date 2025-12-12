# PoliSpectralTools – Release Notes

## v0.2.0 (Initial public release)

- Renamed the package from SpectralTools to **PoliSpectralTools**, updated all
  module references, and added compat bounds (`julia = "1.9"`, `FFTW = "1"`,
  `ToeplitzMatrices = "0.8"`).
- Added an MIT license.
- Completed implementations for:
  - `Collocation.build_grid` (Chebyshev and Legendre Lobatto grids).
  - Linear and nonlinear BVP solvers with boundary-condition enforcement.
  - Diffusion, wave, and 2D Poisson PDE solvers.
- Provided five example scripts (`examples/bvp_linear.jl`, `diffusion.jl`,
  `bvp_legendre.jl`, `bvp_nonlinear.jl`, `wave_mixed_bc.jl`) with figures/GIFs.
- Expanded documentation:
  - `docs/USER_GUIDE.md`, `docs/USAGE_EXAMPLES.md`, and static HTML portals with
    MathJax-enabled content.
  - GitHub Pages configuration via `docs/_config.yml`.
- Added `scripts/run_full_report.jl` producing console + markdown summaries with
  optional coverage statistics.
- Test suite covers generic utilities, BVP/PDE solvers, and usage-driven checks
  (11 tests) passing under Julia ≥ 1.9.
