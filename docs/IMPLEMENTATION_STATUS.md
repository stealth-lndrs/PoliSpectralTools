---
layout: default
title: "Implementation Status"
---

# Implementation Status

This page describes what remains to be implemented for the extended usage
examples and Julia tests defined in `EXAMPLES_TEST_PLAN.md`.  Contributions
should reference the cited slides inside `class_slides/` and follow the guidance
below when preparing scripts, documentation, and automated tests.

## Usage Examples Still Missing

| # | Scenario | Requirements | Slide references |
| --- | --- | --- | --- |
| 6 | Wave propagation with traveling pulse | Implement `examples/wave_pulse.jl` launching a Gaussian pulse `u(x,0) = exp(-50(x+0.5)^2)` with initial velocity `u_t = -c u_x`, Dirichlet BCs on both ends, dt study showing stable vs unstable CFL choices, and optional sponge forcing. Plot snapshots and CFL error table. | `PE_Aula_10_N.pdf` (pp. 2–3) |
| 7 | Diffusion with source forcing | Create `examples/diffusion_forced.jl` solving `u_t = 0.1 u_xx + e^{-t} cos(πx)` on `[-1,1]`, Dirichlet zero walls, initial state zero. Track total “mass” (`∫ u dx`) vs injected energy and compare against a high-resolution finite difference reference. | `PE_Aula_06_N.pdf` (forcing discussion) |
| 8 | Poisson on square | Add `examples/poisson_square.jl` reproducing the manufactured solution `sin(π(x+1)/2) sin(π(y+1)/2)` with Dirichlet BCs. Report convergence when increasing `Nx = Ny`. Include surface/heatmap plots. | `PE_Aula_09_N.pdf` (Programa 16) |
| 9 | Poisson on rectangular domain | Implement `examples/poisson_rectangle.jl` on `x ∈ [-1, 2]`, `y ∈ [0, 1]` with a manufactured solution (e.g., `sin(π(x+1)/3) sinh(π y)`), emphasizing how scaling impacts `Dx²`, `Dy²`. Include anisotropic grid study. | `PE_Aula_09_N.pdf` (pp. 8–13) |
| 10 | Mapping prototype (optional) | Create `examples/mapping_bvp.jl` once `Mapping.jl` exists. Demonstrate mapping `ξ ∈ [-1,1]` to a clustered physical grid (e.g., `x = sin(πξ/2)`), solving a boundary-layer BVP, and comparing errors vs. un-mapped grid. | `PE_Aula_09_N.pdf` (pp. 13–19) |

### Deliverables per example

1. Script under `examples/` (named as suggested above).
2. Corresponding figures/GIFs saved to `docs/assets/`.
3. Documentation snippet added to `docs/USAGE_EXAMPLES.md` and the HTML gallery
   (`docs/web/examples.html`) summarizing the new scenario.

## Julia Tests Still Missing

| # | Test description | Requirements | Linked example |
| --- | --- | --- | --- |
| 6 | Traveling pulse CFL study | Extend `test/runtests.jl` to run the pulse example twice (stable vs unstable dt) and assert the stable run remains bounded while the unstable run diverges rapidly. | Usage Example 6 |
| 7 | Diffusion forcing consistency | Verify the forced diffusion script conserves net integral according to the analytic injection minus diffusion losses. Include dt refinement check. | Usage Example 7 |
| 8 | Poisson square accuracy | Ensure the manufactured square solution achieves `< 1e-6` max error at `Nx=Ny=24` and improves ≥8× when using 32 nodes. | Usage Example 8 |
| 9 | Poisson rectangular scaling | Confirm `Dx²`/`Dy²` matrices reflect domain scaling (`(2/(b-a))^2` factors) and that the rectangular solve hits `< 5e-5` max error. | Usage Example 9 |
| 10 | Mapping impact (optional) | Once mapping exists, add a regression test comparing mapped vs. unmapped errors, expecting ≥3× improvement in the boundary layer region. Skip gracefully if the mapping module is absent. | Usage Example 10 |

### Testing guidance

- Add new `@testset`s beneath the existing “Usage-driven Tests” block in
  `test/runtests.jl`.
- Keep run times under ~5 seconds per test; prefer smaller grids when a
  refinement/ratio check conveys the same message.
- Reuse helper functions (e.g., `build_grid`, `Generic.Bary_Interp`) rather than
  duplicating logic.

## Workflow Checklist

1. **Implementation** – finish the script + figures, update documentation, and
   ensure `examples/<new>.jl` runs.
2. **Testing** – extend `test/runtests.jl`, run
   `julia --project=. -e 'include("test/runtests.jl")'`.
3. **Reporting** – regenerate the full report
   (`julia --project=. scripts/run_full_report.jl`) so the new testsets appear in
   `reports/test_report_<timestamp>.md`.
4. **Docs** – update `docs/USAGE_EXAMPLES.md`, `docs/web/examples.html`, and
   link any new assets. Ensure MathJax renders formulas where needed.

Questions? Post an issue or chat on the project’s discussion board; be sure to
mention which example/test number you are tackling so efforts don’t overlap.
