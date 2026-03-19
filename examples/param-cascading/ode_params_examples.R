library(ggplot2)
library(gridExtra)
library(reshape2)

source("src/solvers/parameter_cascading.R")
source("src/utils/trace_optimisation.R")
source("examples/ode_models.R")                         # euler_solve, ODE_CONFIGS
source("examples/sensitivity-analysis/sensitivity_utils.R")  # run_local_sensitivity


# =============================================================================
# run_example — canonical CascadingOdeSolver workflow
#
# Arguments
# ---------
# cfg          ODE_CONFIGS entry from examples/ode_models.R.
#              Provides: rhs (solver format), params (true generating
#              parameters), y0, times_obs, times_sim, fixed_params,
#              param_scales, active_params, uncertain_bounds.
# init_params  Named numeric vector: initial guess for the estimated
#              parameters.  Names determine which parameters are estimated.
# param_names  Character vector of parameter names to estimate.
#              Defaults to names(init_params).
# lambda       Regularisation weight for the inner OdeSystemSolver (controls
#              how much of the data mismatch is attributed to u(t) vs.
#              parameter error).
# lower_phys   Named numeric: lower optimisation bounds.
#              Defaults to the lower limits in cfg$uncertain_bounds.
# upper_phys   Named numeric: upper optimisation bounds.
#              Defaults to the upper limits in cfg$uncertain_bounds.
# noise_sd     Std. dev. of additive Gaussian noise added to synthetic data.
#              NULL = auto: 2 % of peak signal amplitude.
# max_iter     Max iterations for the outer parameter optimisation.
# seed         Random seed for data generation.
#
# Workflow
# --------
# 1. Generate synthetic observations from cfg$params (true parameters).
# 2. Run CascadingOdeSolver$optimize_parameters() starting from init_params.
# 3. Report true vs estimated parameter values.
# 4. Plot state fit, residual forcing u(t), and outer parameter trace.
# 5. Run local sensitivity at the true generating parameters:
#    - central-difference Jacobian -> Fisher Information Matrix
#    - FIM rank / condition number (parameter identifiability)
#    - parameter correlation matrix
#    - time-resolved RMS sensitivity plot
#
# Returns invisibly: list(cascading = <solver>, result = <named numeric>).
# =============================================================================
run_example <- function(cfg,
                        init_params,
                        param_names  = names(init_params),
                        lambda       = 1e0,
                        lower_phys   = NULL,
                        upper_phys   = NULL,
                        noise_sd     = NULL,
                        max_iter     = 200,
                        seed         = 123) {
  rhs          <- cfg$rhs
  params       <- cfg$params
  y0           <- cfg$y0
  t_obs        <- cfg$times_obs
  t_sim        <- cfg$times_sim
  fixed_params <- cfg$fixed_params
  param_scales <- cfg$param_scales
  n_vars       <- length(y0)
  var_names    <- if (!is.null(names(y0))) names(y0) else
    paste0("Var", seq_len(n_vars))

  # Default bounds from cfg$uncertain_bounds when available
  if (is.null(lower_phys)) {
    lower_phys <- setNames(
      vapply(cfg$uncertain_bounds[param_names],
             function(b) b[1], numeric(1)),
      param_names
    )
  }
  if (is.null(upper_phys)) {
    upper_phys <- setNames(
      vapply(cfg$uncertain_bounds[param_names],
             function(b) b[2], numeric(1)),
      param_names
    )
  }

  # 1. Synthetic data ---------------------------------------------------------
  y_true  <- euler_solve(y0, t_sim, rhs, params)
  obs_idx <- match(round(t_obs, 10), round(t_sim, 10))
  stopifnot("obs times must be a subset of sim times" = !anyNA(obs_idx))

  sd_use <- if (is.null(noise_sd))
    0.02 * max(abs(y_true[obs_idx, , drop = FALSE]))
  else
    noise_sd
  set.seed(seed)
  obs_data <- y_true[obs_idx, , drop = FALSE] +
    matrix(rnorm(length(obs_idx) * n_vars, 0, sd_use),
           length(obs_idx), n_vars)

  # 2. Cascading solver -------------------------------------------------------
  cat("\n=== CascadingOdeSolver  params:", param_names, " lambda =", lambda,
      "===\n")
  cascading <- CascadingOdeSolver$new(
    func_rhs     = rhs,
    times_sim    = t_sim,
    obs_times    = t_obs,
    obs_values   = obs_data,
    y0           = y0,
    fixed_params = fixed_params,
    lambda       = lambda,
    param_scales = param_scales
  )
  result <- cascading$optimize_parameters(
    init_theta_physical = init_params,
    param_names         = param_names,
    #lower_phys          = lower_phys,
    #upper_phys          = upper_phys
  )

  # 3. Parameter recovery report ----------------------------------------------
  true_subset <- unlist(params[param_names])
  cat("\n--- Parameter recovery ---\n")
  for (nm in param_names)
    cat(sprintf("  %s:  true = %10.4g  |  estimated = %10.4g\n",
                nm, true_subset[nm], result[nm]))

  # 4. Plots ------------------------------------------------------------------
  s           <- cascading$last_solver
  sim_obs_idx <- which(round(t_sim, 10) %in% round(t_obs, 10))

  ps <- lapply(seq_len(n_vars), function(i) {
    ggplot(
      data.frame(t   = t_obs,
                 obs = obs_data[, i],
                 fit = s$y[sim_obs_idx, i]),
      aes(x = t)
    ) +
      geom_point(aes(y = obs), alpha = 0.5, size = 1) +
      geom_line(aes(y = fit), color = "steelblue", linewidth = 1) +
      labs(title = var_names[i], x = "Time", y = "State") +
      theme_minimal()
  })

  pu <- lapply(seq_len(n_vars), function(i) {
    ggplot(
      data.frame(t = t_sim, u = s$u[, i]),
      aes(x = t, y = u)
    ) +
      geom_line(color = "tomato", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(title = paste0(var_names[i], "  u(t)"), x = "Time", y = "u") +
      theme_minimal()
  })

  do.call(grid.arrange,
          c(ps,  list(ncol = min(n_vars, 3), top = "State fit")))
  do.call(grid.arrange,
          c(pu,  list(ncol = min(n_vars, 3), top = "Residual forcing u(t)")))
  plot_param_trace(cascading, true_params = as.list(true_subset))

  # 5. Local sensitivity at true generating parameters -----------------------
  cat("\n--- Local sensitivity at true generating parameters ---\n")
  cat("    (central-difference Jacobian; deSolve integrator)\n")
  sens_cfg <- make_sens_config(cfg)
  run_local_sensitivity(sens_cfg, verbose = TRUE)

  invisible(list(cascading = cascading, result = result))
}


# =============================================================================
# Example calls
# =============================================================================
# Exponential decay — estimate k
run_example(
   cfg         = ODE_CONFIGS$decay,
   init_params = c(k = 2.0),
   lambda      = 1e0
 )

# Lotka-Volterra — estimate all four parameters
run_example(
   cfg         = ODE_CONFIGS$lv,
   init_params = c(alpha = 0.8, beta = 0.4, delta = 0.4, gamma = 0.6),
   lambda      = 1e0
 )

# Sestak-Berggren multi-temperature — estimate E, n, m (A fixed)
# run_example(
#   cfg         = ODE_CONFIGS$sb,
#   init_params = c(E = 50000, n = 3.0, m = 0.5),
#   param_names = c("E", "n", "m"),
#   lower_phys  = c(E = 10000, n = 0.1, m = 0.1),
#   upper_phys  = c(E = 300000, n = 20, m = 20),
#   lambda      = 1e1
# )
