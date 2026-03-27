library(ggplot2)
library(gridExtra)
library(reshape2)

source("src/solvers/tracking_ode_solver.R")
source("src/utils/trace_optimisation.R")
source("examples/ode_models.R")
source("examples/plot_utils.R")
source("examples/sensitivity-analysis/sensitivity_utils.R")


# =============================================================================
# run_example -- canonical TrackingOdeSolver workflow
#
# Arguments
# ---------
# cfg          ODE_CONFIGS entry from examples/ode_models.R.
# init_params  Named numeric: initial guess for estimated parameters.
# param_names  Character vector of parameters to estimate.
# lambda       Regularisation weight.
# lower_phys   Named numeric: lower optimisation bounds.
# upper_phys   Named numeric: upper optimisation bounds.
# noise_sd     Noise std dev (NULL = auto: 2% of peak signal).
# max_iter     Max outer iterations.
# seed         Random seed for data generation.
#
# Returns invisibly: list(tracking = <solver>, result = <named numeric>).
# =============================================================================
run_example <- function(cfg,
                        init_params,
                        param_names = names(init_params),
                        lambda = 1e0,
                        lower_phys = NULL,
                        upper_phys = NULL,
                        noise_sd = NULL,
                        max_iter = 200,
                        seed = 123) {
  y0     <- cfg$y0
  t_obs  <- cfg$times_obs
  t_sim  <- cfg$times_sim
  n_vars <- length(y0)
  var_names <- if (!is.null(names(y0))) names(y0) else
    paste0("Var", seq_len(n_vars))

  if (is.null(lower_phys)) {
    lower_phys <- setNames(
      vapply(cfg$uncertain_bounds[param_names],
             function(b) b[1], numeric(1)),
      param_names)
  }
  if (is.null(upper_phys)) {
    upper_phys <- setNames(
      vapply(cfg$uncertain_bounds[param_names],
             function(b) b[2], numeric(1)),
      param_names)
  }

  # 1. Synthetic data
  syn <- generate_synthetic_data(cfg, noise_sd = noise_sd, seed = seed)

  # 2. Tracking solver
  cat("\n=== TrackingOdeSolver  params:", param_names,
      " lambda =", lambda, "===\n")
  tracking <- TrackingOdeSolver$new(
    func_rhs     = cfg$rhs,
    times_sim    = t_sim,
    obs_times    = t_obs,
    obs_values   = syn$obs_data,
    init_state   = function(p) as.numeric(y0),
    fixed_params = cfg$fixed_params,
    lambda       = lambda,
    param_scales = cfg$param_scales
  )
  result <- tracking$optimize_parameters(
    init_theta_physical = init_params,
    param_names         = param_names
  )

  # 3. Parameter recovery report
  true_subset <- unlist(cfg$params[param_names])
  cat("\n--- Parameter recovery ---\n")
  for (nm in param_names)
    cat(sprintf("  %s:  true = %10.4g  |  estimated = %10.4g\n",
                nm, true_subset[nm], result[nm]))

  # 4. Plots
  s <- tracking$last_solver
  plot_state_fit(s$y, syn$obs_data, t_obs, t_sim, var_names)
  plot_residual_forcing(s$u, t_sim, var_names)
  plot_outer_trace(tracking, true_params = as.list(true_subset),
                   cost_field = "cost",
                   cost_label = "J (SSE/n + lambda*||u*||^2)",
                   title = "Parameter Tracking \u2014 Optimisation Progress")

  # 5. Local sensitivity
  cat("\n--- Local sensitivity at true generating parameters ---\n")
  sens_cfg <- make_sens_config(cfg)
  run_local_sensitivity(sens_cfg, verbose = TRUE)

  invisible(list(tracking = tracking, result = result))
}


# =============================================================================
# Example calls
# =============================================================================
# Exponential decay -- estimate k
# run_example(
#   cfg         = ODE_CONFIGS$decay,
#   init_params = c(k = 2.0),
#   lambda      = 1e0
# )

# Lotka-Volterra -- estimate all four parameters
# run_example(
#   cfg         = ODE_CONFIGS$lv,
#   init_params = c(alpha = 0.8, beta = 0.4, delta = 0.4, gamma = 0.6),
#   lambda      = 1e0
# )

# Sestak-Berggren multi-temperature -- estimate E, n, m (A fixed)
run_example(
  cfg         = ODE_CONFIGS$sb,
  init_params = c(E = 50000, n = 3.0, m = 0.5),
  param_names = c("E", "n", "m"),
  lower_phys  = c(E = 10000, n = 0.1, m = 0.1),
  upper_phys  = c(E = 300000, n = 20, m = 20),
  lambda      = 1e1
)
