library(ggplot2)
library(gridExtra)
library(reshape2)

source("src/solvers/inverse-solvers/load_inverse_solvers.R")
source("src/utils/trace_optimisation.R")
source("examples/ode_models.R")
source("examples/plot_utils.R")


# =============================================================================
# compare_methods -- run both TrackingOdeSolver and CascadingOdeSolver on the
# same synthetic data and compare parameter recovery, state fit, and convergence
#
# Arguments
# ---------
# cfg          ODE_CONFIGS entry
# init_params  Named numeric: initial guess for estimated parameters
# param_names  Character vector of parameters to estimate
# lambda       Regularisation weight
# lower_phys   Named numeric: lower bounds (NULL = from cfg)
# upper_phys   Named numeric: upper bounds (NULL = from cfg)
# noise_sd     Noise std dev (NULL = auto: 2% of peak signal)
# max_iter     Max outer iterations
# seed         Random seed
#
# Returns invisibly: list(tracking, cascading, result_tracking, result_cascading)
# =============================================================================
compare_methods <- function(cfg,
                            init_params,
                            param_names = names(init_params),
                            lambda      = 1e0,
                            lower_phys  = NULL,
                            upper_phys  = NULL,
                            noise_sd    = NULL,
                            max_iter    = 200,
                            seed        = 123) {
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

  # --- Shared synthetic data ---
  syn <- generate_synthetic_data(cfg, noise_sd = noise_sd, seed = seed)
  true_subset <- unlist(cfg$params[param_names])

  # --- Tracking estimator ---
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
    param_scales = cfg$param_scales,
    inner_max_iter = max_iter
  )
  result_tr <- tracking$optimize_parameters(
    init_theta_physical = init_params,
    param_names         = param_names,
    lower_phys          = lower_phys,
    upper_phys          = upper_phys
  )

  # --- Parameter cascading ---
  cat("\n=== CascadingOdeSolver  params:", param_names,
      " lambda =", lambda, "===\n")
  cascading <- CascadingOdeSolver$new(
    func_rhs     = cfg$rhs,
    times_sim    = t_sim,
    obs_times    = t_obs,
    obs_values   = syn$obs_data,
    init_state   = function(p) as.numeric(y0),
    fixed_params = cfg$fixed_params,
    lambda       = lambda,
    param_scales = cfg$param_scales,
    inner_max_iter = max_iter
  )
  result_pc <- cascading$optimize_parameters(
    init_theta_physical = init_params,
    param_names         = param_names,
    lower_phys          = lower_phys,
    upper_phys          = upper_phys
  )

  # --- Parameter recovery table ---
  cat("\n--- Parameter recovery ---\n")
  cat(sprintf("  %-10s %12s %12s %12s\n",
              "param", "true", "tracking", "cascading"))
  for (nm in param_names)
    cat(sprintf("  %-10s %12.4g %12.4g %12.4g\n",
                nm, true_subset[nm], result_tr[nm], result_pc[nm]))

  # --- Comparison plots ---
  s_tr <- tracking$last_solver
  s_pc <- cascading$last_solver

  # State fit overlay
  ps <- lapply(seq_len(n_vars), function(i) {
    df <- rbind(
      data.frame(t = t_sim, y = s_tr$y[, i], method = "Tracking"),
      data.frame(t = t_sim, y = s_pc$y[, i], method = "Cascading")
    )
    df_obs <- data.frame(t = t_obs, y = syn$obs_data[, i])

    ggplot() +
      geom_point(data = df_obs, aes(x = t, y = y),
                 alpha = 0.4, size = 1) +
      geom_line(data = df, aes(x = t, y = y, color = method),
                linewidth = 0.8) +
      scale_color_manual(values = c(Tracking = "steelblue",
                                    Cascading = "tomato")) +
      labs(title = var_names[i], x = "Time", y = "State") +
      theme_minimal()
  })
  do.call(grid.arrange,
          c(ps, list(ncol = min(n_vars, 3), top = "State fit comparison")))

  # Residual forcing overlay
  pu <- lapply(seq_len(n_vars), function(i) {
    df <- rbind(
      data.frame(t = t_sim, u = s_tr$u[, i], method = "Tracking"),
      data.frame(t = t_sim, u = s_pc$u[, i], method = "Cascading")
    )
    ggplot(df, aes(x = t, y = u, color = method)) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c(Tracking = "steelblue",
                                    Cascading = "tomato")) +
      labs(title = paste0(var_names[i], "  u(t)"), x = "Time", y = "u") +
      theme_minimal()
  })
  do.call(grid.arrange,
          c(pu, list(ncol = min(n_vars, 3),
                     top = "Residual forcing comparison")))

  # Convergence traces side by side
  plot_outer_trace(tracking, true_params = as.list(true_subset),
                   cost_field = "cost",
                   cost_label = "J (SSE/n + lambda*||u*||^2)",
                   title = "Tracking — Optimisation Progress")
  plot_outer_trace(cascading, true_params = as.list(true_subset),
                   cost_field = "sse",
                   cost_label = "SSE (log scale)",
                   title = "Cascading — Optimisation Progress")

  invisible(list(
    tracking         = tracking,
    cascading        = cascading,
    result_tracking  = result_tr,
    result_cascading = result_pc
  ))
}


# =============================================================================
# Example calls
# =============================================================================

# Exponential decay -- estimate k
compare_methods(
  cfg         = ODE_CONFIGS$decay,
  init_params = c(k = 2.0),
  lambda      = 1e0
)

# Lotka-Volterra -- estimate all four parameters
compare_methods(
  cfg         = ODE_CONFIGS$lv,
  init_params = c(alpha = 0.8, beta = 0.4, delta = 0.4, gamma = 0.6),
  lambda      = 1e0
)

# Sestak-Berggren multi-temperature -- estimate E, n, m
compare_methods(
  cfg         = ODE_CONFIGS$sb,
  init_params = c(E = 50000, n = 3.0, m = 0.2),
  param_names = c("E", "n", "m"),
  #lower_phys  = c(E = 10000, n = 0.1, m = 0.1),
  #upper_phys  = c(E = 300000, n = 20, m = 20),
  lambda      = 1e2,
  max_iter = 50
)
