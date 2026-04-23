library(ggplot2)
library(gridExtra)
library(reshape2)

source("src/solvers/inverse-solvers/load_inverse_solvers.R")
source("src/utils/trace_optimisation.R")
source("examples/ode_models.R")
source("examples/inverse-solvers/inverse_solver_example_helpers.R")
source("examples/plot_utils.R")
source("examples/extra/sensitivity-analysis/sensitivity_utils.R")


# =============================================================================
# run_example -- canonical TrackingInverseSolver workflow
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
  regularizing_ode <- build_inverse_solver_model(cfg)

  # 2. Tracking solver
  cat("\n=== TrackingInverseSolver  params:", param_names,
      " lambda =", lambda, "===\n")
  tracking <- TrackingInverseSolver$new(
    model        = regularizing_ode,
    times_sim    = t_sim,
    obs_times    = t_obs,
    obs_values   = syn$obs_data,
    lambda       = lambda,
    verbose = T
  )
  result <- tracking$optimize_parameters(
    init_theta_physical = init_params,
    param_names         = param_names,
    lower_phys = lower_phys,
    upper_phys = upper_phys)

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


main <- function() {
  run_example(
    cfg = ODE_CONFIGS$lv,
    init_params = c(alpha = 1.0, beta = 0.4, delta = 0.4, gamma = 0.9),
    lambda = 1e1
  )
}

if (sys.nframe() == 0L) {
  main()
}
