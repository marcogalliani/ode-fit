# Inverse Solvers Consistency Study (Lotka-Volterra)

This folder contains a simulation study to assess consistency of inverse ODE parameter estimates as the number of observations increases, comparing:

- Tracking approach (`TrackingOdeSolver`)
- Parameter cascading approach (`CascadingInverseSolver`)

The setup mirrors the Lotka-Volterra configuration used in the inverse-solvers test suite:

- True parameters: `alpha=1.20`, `beta=0.45`, `delta=0.12`, `gamma=0.80`
- Fixed initial conditions: `x0=10`, `y0=8`
- Simulation grid: `times_sim = seq(0, 5, by = 0.1)`
- Inner scheme: `gl2`

## Script

- `lv_inverse_solvers_consistency.R`

## What the script does

For each observation count in a grid and each Monte Carlo replicate:

1. Subsamples observation times from `times_sim`.
2. Simulates noisy observations from the true Lotka-Volterra trajectory.
3. Fits both inverse solvers to estimate `alpha, beta, delta, gamma`.
4. Stores parameter estimates and errors.
5. Aggregates errors by method and sample size.
6. Saves CSV summaries and PNG plots.

## Run

From repository root:

```bash
Rscript simulations/inverse-solvers-consistency/lv_inverse_solvers_consistency.R
```

Optional CLI arguments:

1. `n_rep`: number of Monte Carlo replications (default `30`)
2. `n_obs_grid_csv`: comma-separated observation counts (default `6,11,16,21,31,41,51`)
3. `inner_max_iter`: inner solver iterations (default `200`)

Example:

```bash
Rscript simulations/inverse-solvers-consistency/lv_inverse_solvers_consistency.R 20 6,11,21,31,51 200
```

## Outputs

Written to:

- `simulations/inverse-solvers-consistency/results/lv_consistency_raw_results.csv`
- `simulations/inverse-solvers-consistency/results/lv_consistency_summary_l2.csv`
- `simulations/inverse-solvers-consistency/results/lv_consistency_summary_params.csv`
- `simulations/inverse-solvers-consistency/results/lv_consistency_l2_error.png`
- `simulations/inverse-solvers-consistency/results/lv_consistency_param_errors.png`
