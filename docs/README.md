# Documentation Notes

This folder hosts high-level documentation for SpectralTools.jl.  It is intentionally lightweight so
students can expand it with worked examples or migrate it to Documenter.jl later.

## Suggested Layout

- `overview.md` – outline of modules (Chebyshev, Legendre, Fourier, Generic, Collocation, BVP, PDE).
- `bvp.md` – derivations behind the collocation-based boundary value solvers plus convergence tables.
- `pde.md` – notes on the diffusion, wave, and Poisson solvers, together with stability restrictions.
- `examples/` – short scripts referenced from the README to replicate experiments.

When extending the docs, remember to activate the local environment (`julia --project`) so snippets
run with the same dependencies that power the test suite.
