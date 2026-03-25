# =============================================================================
# examples/uncertainty/laplace_ci_examples.R
#
# Laplace / NLS uncertainty quantification for ODE trajectory estimates.
#
# Both CascadingOdeSolver and TrackingOdeSolver expose compute_uncertainty(),
# which implements:
#
#   1. NLS parameter covariance
#        Sigma_theta = sigma^2 * (J^T J)^{-1}
#        sigma^2     = SSE / (n_obs - n_params)
#        J           = sensitivity of y*(t) w.r.t. theta at observed entries
#
#   2. Delta-method confidence band for y*(t)
#        Var(y*_{t,v}) = s_{t,v}^T * Sigma_theta * s_{t,v}
#        CI            = y*_{t,v}  +/-  z_{alpha/2} * sqrt(Var(y*_{t,v}))
#
# The sensitivity J differs between the two methods:
#   Cascading — TOTAL: accounts for how u*(theta) shifts with theta (IFT BVP)
#   Tracking  — PARTIAL: u* held fixed; consistent with the envelope theorem
#
# Run from repo root:
#   source("examples/uncertainty/laplace_ci_examples.R")
# =============================================================================

library(ggplot2)
library(gridExtra)
library(reshape2)

source("src/solvers/parameter_cascading.R")
source("src/solvers/tracking_ode_solver.R")
source("examples/ode_models.R")


# =============================================================================
# plot_ci_band  —  CI ribbon + fit line for one state variable
# =============================================================================
plot_ci_band <- function(t_sim, y_hat, lower, upper,
                         t_obs, obs_col, true_col = NULL,
                         var_name = "y", method_label = "",
                         color = "steelblue") {
  df_sim <- data.frame(t = t_sim, y_hat = y_hat, lower = lower, upper = upper)
  df_obs <- data.frame(t = t_obs, obs = obs_col)

  p <- ggplot(df_sim, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
      fill = color, alpha = 0.25, color = NA
    ) +
    geom_line(aes(y = y_hat), color = color, linewidth = 1) +
    geom_point(
      data = df_obs, aes(y = obs),
      color = "gray30", size = 1.5, alpha = 0.8
    ) +
    labs(
      title = sprintf("%s  [%s]", var_name, method_label),
      x = "Time", y = var_name
    ) +
    theme_minimal(base_size = 11)

  if (!is.null(true_col)) {
    df_true <- data.frame(t = t_sim, y_true = true_col)
    p <- p + geom_line(
      data = df_true, aes(y = y_true),
      color = "gray30", linetype = "dashed", linewidth = 0.7
    )
  }
  p
}


# =============================================================================
# plot_se_time  —  SE(y*(t)) vs time, one curve per state variable
# =============================================================================
plot_se_time <- function(t_sim, se_mat, var_names, method_label = "") {
  colnames(se_mat) <- var_names
  df <- reshape2::melt(
    data.frame(t = t_sim, as.data.frame(se_mat)),
    id.vars = "t", variable.name = "Variable", value.name = "SE"
  )
  ggplot(df, aes(x = t, y = SE, color = Variable)) +
    geom_line(linewidth = 0.9) +
    labs(
      title = sprintf("SE of y*(t)  [%s]", method_label),
      x = "Time", y = "Standard Error", color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = if (length(var_names) > 1) "bottom" else "none")
}


# =============================================================================
# plot_ci_overlay  —  side-by-side Cascading vs Tracking on the same panel
# =============================================================================
plot_ci_overlay <- function(t_sim, ci_c, ci_t,
                            t_obs, obs_col, true_col = NULL,
                            var_idx = 1, var_name = "y", alpha_label = "95") {
  df_all <- rbind(
    data.frame(
      t = t_sim,
      y_hat = ci_c$y_hat[, var_idx],
      lower = ci_c$lower[, var_idx],
      upper = ci_c$upper[, var_idx],
      Method = "Cascading"
    ),
    data.frame(
      t = t_sim,
      y_hat = ci_t$y_hat[, var_idx],
      lower = ci_t$lower[, var_idx],
      upper = ci_t$upper[, var_idx],
      Method = "Tracking"
    )
  )
  df_obs <- data.frame(t = t_obs, obs = obs_col)

  p <- ggplot(df_all, aes(x = t, color = Method, fill = Method)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
      alpha = 0.20, color = NA
    ) +
    geom_line(aes(y = y_hat), linewidth = 1) +
    geom_point(
      data = df_obs, aes(x = t, y = obs),
      inherit.aes = FALSE,
      color = "gray30", size = 1.5, alpha = 0.8
    ) +
    scale_color_manual(values = c(Cascading = "steelblue", Tracking = "tomato")) +
    scale_fill_manual(values = c(Cascading = "steelblue", Tracking = "tomato")) +
    guides(fill = "none") +
    labs(
      title = sprintf(
        "%s  —  %s%% CI: Cascading vs Tracking",
        var_name, alpha_label
      ),
      x = "Time", y = var_name, color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  if (!is.null(true_col)) {
    df_true <- data.frame(t = t_sim, y_true = true_col)
    p <- p + geom_line(
      data = df_true, aes(x = t, y = y_true),
      inherit.aes = FALSE,
      color = "gray30", linetype = "dashed", linewidth = 0.7
    )
  }
  p
}


# =============================================================================
# report_param_ci  —  print parameter estimates with NLS confidence intervals
# =============================================================================
report_param_ci <- function(ci, theta_phys, param_names,
                            true_params = NULL, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  se_t <- sqrt(diag(ci$Sigma_theta))

  cat(sprintf("  sigma^2 = %.4g   dof = %f\n", ci$sigma2, ci$dof))
  cat(sprintf(
    "  %-12s  %10s  %8s  %26s\n",
    "parameter", "estimate", "SE",
    sprintf("%.0f%% CI", 100 * (1 - alpha))
  ))
  cat("  ", strrep("-", 62), "\n", sep = "")
  for (j in seq_along(param_names)) {
    nm <- param_names[j]
    est <- theta_phys[j]
    lo <- est - z * se_t[j]
    hi <- est + z * se_t[j]
    tru <- if (!is.null(true_params) && nm %in% names(true_params)) {
      sprintf("  (true: %g)", true_params[[nm]])
    } else {
      ""
    }
    cat(sprintf(
      "  %-12s  %10.4g  %8.4g  [%10.4g, %10.4g]%s\n",
      nm, est, se_t[j], lo, hi, tru
    ))
  }
  cat("\n")
}


# =============================================================================
# run_uncertainty_example
#
# Fits parameters with BOTH CascadingOdeSolver and TrackingOdeSolver, calls
# compute_uncertainty() on each, and produces:
#
#   Plot 1  —  95% CI band per state variable, Cascading and Tracking side by side
#   Plot 2  —  SE(y*(t)) vs time, both methods
#   Plot 3  —  Overlay: Cascading vs Tracking CI bands on the same panel
#   Console —  Parameter estimates, SE, and 95% CI for both methods
#
# Arguments
# ---------
# cfg         ODE_CONFIGS entry (from examples/ode_models.R)
# init_params Named numeric: initial parameter guess
# param_names Names of parameters to estimate (default: names(init_params))
# lambda      Inner regularisation weight (fixed throughout; not estimated)
# noise_sd    Noise SD for synthetic data (NULL = 2% of peak signal)
# alpha       Significance level for CI (default 0.05 → 95% bands)
# seed        RNG seed for data generation
#
# Returns invisibly: list(cascading, tracking, theta_c, theta_t, ci_c, ci_t)
# =============================================================================
run_uncertainty_example <- function(cfg,
                                    init_params,
                                    param_names = names(init_params),
                                    lambda = 1e0,
                                    noise_sd = NULL,
                                    alpha = 0.05,
                                    seed = 123) {
  y0 <- cfg$y0
  t_obs <- cfg$times_obs
  t_sim <- cfg$times_sim
  n_vars <- length(y0)
  var_names <- if (!is.null(names(y0))) {
    names(y0)
  } else {
    paste0("Var", seq_len(n_vars))
  }
  alpha_pct <- sprintf("%.0f", 100 * (1 - alpha))

  # 1. Synthetic data -----------------------------------------------------------
  syn <- generate_synthetic_data(cfg, noise_sd = noise_sd, seed = seed)
  obs_data <- syn$obs_data
  y_true   <- syn$y_true

  true_subset <- unlist(cfg$params[param_names])

  # 2. Fit: CascadingOdeSolver -------------------------------------------------
  cat(
    "\n=== CascadingOdeSolver: fitting parameters:", param_names,
    " | lambda =", lambda, "===\n"
  )
  cascading <- CascadingOdeSolver$new(
    func_rhs     = cfg$rhs,
    times_sim    = t_sim,
    obs_times    = t_obs,
    obs_values   = obs_data,
    y0           = y0,
    fixed_params = cfg$fixed_params,
    lambda       = lambda,
    param_scales = cfg$param_scales
  )
  theta_c <- cascading$optimize_parameters(init_params, param_names)

  # 3. Fit: TrackingOdeSolver --------------------------------------------------
  cat(
    "\n=== TrackingOdeSolver: fitting parameters:", param_names,
    " | lambda =", lambda, "===\n"
  )
  tracking <- TrackingOdeSolver$new(
    func_rhs     = cfg$rhs,
    times_sim    = t_sim,
    obs_times    = t_obs,
    obs_values   = obs_data,
    y0           = y0,
    fixed_params = cfg$fixed_params,
    lambda       = lambda,
    param_scales = cfg$param_scales
  )
  theta_t <- tracking$optimize_parameters(init_params, param_names)

  # 4. Uncertainty -------------------------------------------------------------
  cat("\n=== Computing uncertainty bands ===\n")
  ci_c <- cascading$compute_uncertainty(param_names = param_names, alpha = alpha)
  ci_t <- tracking$compute_uncertainty(param_names = param_names, alpha = alpha)

  # Console report
  cat("\n--- Cascading: parameter uncertainty ---\n")
  report_param_ci(ci_c, theta_c, param_names, true_subset, alpha)
  cat("--- Tracking:  parameter uncertainty ---\n")
  report_param_ci(ci_t, theta_t, param_names, true_subset, alpha)

  # 5. Plot: per-variable CI bands (Cascading | Tracking) ---------------------
  for (i in seq_len(n_vars)) {
    p_c <- plot_ci_band(
      t_sim, ci_c$y_hat[, i], ci_c$lower[, i], ci_c$upper[, i],
      t_obs, obs_data[, i], y_true[, i],
      var_name = var_names[i], method_label = "Cascading",
      color = "steelblue"
    )
    p_t <- plot_ci_band(
      t_sim, ci_t$y_hat[, i], ci_t$lower[, i], ci_t$upper[, i],
      t_obs, obs_data[, i], y_true[, i],
      var_name = var_names[i], method_label = "Tracking",
      color = "tomato"
    )
    grid.arrange(
      p_c, p_t,
      ncol = 2,
      top = sprintf(
        "%s%%  CI for %s  (lambda = %g, dashed = true)",
        alpha_pct, var_names[i], lambda
      )
    )
  }

  # 6. Plot: SE(y*(t)) vs time -------------------------------------------------
  p_se_c <- plot_se_time(t_sim, ci_c$se, var_names, "Cascading")
  p_se_t <- plot_se_time(t_sim, ci_t$se, var_names, "Tracking")
  grid.arrange(p_se_c, p_se_t,
    ncol = 2,
    top = "SE of y*(t): grows from 0 at t=0 (y0 known)"
  )

  # 7. Plot: overlay per variable ----------------------------------------------
  for (i in seq_len(n_vars)) {
    print(plot_ci_overlay(
      t_sim, ci_c, ci_t,
      t_obs, obs_data[, i], y_true[, i],
      var_idx = i, var_name = var_names[i], alpha_label = alpha_pct
    ))
  }

  invisible(list(
    cascading = cascading, tracking = tracking,
    theta_c   = theta_c,   theta_t  = theta_t,
    ci_c      = ci_c,      ci_t     = ci_t
  ))
}


# =============================================================================
# Example 1: Exponential decay  (1 parameter, 1 variable)
#
# True k = 0.5.  Starting guess k = 0.8.
# The SE band grows monotonically from 0 (y0 known) because uncertainty in k
# accumulates along the trajectory via the ODE dynamics.
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("Example 1: Exponential decay  —  estimate k\n")
cat(strrep("=", 70), "\n")

out_decay <- run_uncertainty_example(
  cfg         = ODE_CONFIGS$decay,
  init_params = c(k = 0.8),
  lambda      = 0.5
)


# =============================================================================
# Example 2: Lotka-Volterra  (4 parameters, 2 variables)
#
# True: alpha=1.1, beta=0.4, delta=0.1, gamma=0.4.
# Starting from the true values so the solver converges quickly.
# The SE bands for x and y will differ in shape: x (prey) and y (predator)
# have different sensitivities to each parameter.
#
# Note: the LV sim grid is seq(0,5,by=0.001) — 5001 steps.  The sensitivity
# BVP is cheap (2x2 matrices), but expect ~10–30 s for each solver.
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("Example 2: Lotka-Volterra  —  estimate alpha, beta, delta, gamma\n")
cat(strrep("=", 70), "\n")

out_lv <- run_uncertainty_example(
  cfg         = ODE_CONFIGS$lv,
  init_params = c(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4),
  lambda      = 1.0
)
