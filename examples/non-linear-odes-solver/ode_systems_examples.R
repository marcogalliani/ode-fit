library(ggplot2)
library(gridExtra)
library(reshape2)
library(plotly)

source("src/solvers/general_ode_system_solver.R")
source("src/utils/trace_optimisation.R")
source("examples/ode_models.R")
source("examples/sensitivity-analysis/sensitivity_utils.R")


# =============================================================================
# run_example — canonical OdeSystemSolver workflow
#
# Arguments
# ---------
# cfg        ODE_CONFIGS entry from examples/ode_models.R.
#            Provides: rhs (solver format), params (true physics), y0,
#            times_obs, times_sim, active_params, uncertain_bounds.
# lambda     Regularisation weight penalising large |u(t)|.
#            The solver estimates u(t) to bridge model-data mismatch;
#            higher lambda forces u(t) ~ 0 (trusts the physics more).
# noise_sd   Std. dev. of additive Gaussian observation noise.
#            NULL = auto: 5 % of the peak signal amplitude.
# max_iter   Max L-BFGS-B iterations for the u(t) optimisation.
# seed       Random seed for reproducible data generation.
#
# Workflow
# --------
# 1. Generate synthetic observations from cfg$params (true parameters).
# 2. Run OdeSystemSolver$optimize() to recover the unknown forcing u(t).
# 3. Plot state fit and residual forcing for every state variable.
# 4. Run local sensitivity at the true generating parameters:
#    - central-difference Jacobian -> Fisher Information Matrix
#    - FIM rank / condition number (parameter identifiability)
#    - parameter correlation matrix
#    - time-resolved RMS sensitivity plot
#
# Returns the fitted OdeSystemSolver object invisibly.
# =============================================================================
run_example <- function(cfg,
                        lambda    = 1e5,
                        noise_sd  = NULL,
                        max_iter  = 200,
                        seed      = 123) {
  rhs    <- cfg$rhs
  params <- cfg$params
  y0     <- cfg$y0
  t_obs  <- cfg$times_obs
  t_sim  <- cfg$times_sim
  n_vars <- length(y0)
  var_names <- if (!is.null(names(y0))) names(y0) else
    paste0("Var", seq_len(n_vars))

  # 1. Synthetic data ---------------------------------------------------------
  y_true  <- euler_solve(y0, t_sim, rhs, params)
  obs_idx <- match(round(t_obs, 10), round(t_sim, 10))
  stopifnot("obs times must be a subset of sim times" = !anyNA(obs_idx))

  sd_use <- if (is.null(noise_sd))
    0.05 * max(abs(y_true[obs_idx, , drop = FALSE]))
  else
    noise_sd
  set.seed(seed)
  obs_data <- y_true[obs_idx, , drop = FALSE] +
    matrix(rnorm(length(obs_idx) * n_vars, 0, sd_use),
           length(obs_idx), n_vars)

  # 2. Solver -----------------------------------------------------------------
  cat("\n=== OdeSystemSolver  lambda =", lambda, "===\n")
  solver <- OdeSystemSolver$new(
    func_rhs   = rhs,
    obs_times  = t_obs,
    times_sim  = t_sim,
    obs_values = obs_data,
    params     = params,
    lambda     = lambda
  )
  solver$optimize(y0 = y0, max_iter = max_iter)

  # 3. Plots ------------------------------------------------------------------
  sim_obs_idx <- which(round(t_sim, 10) %in% round(t_obs, 10))

  ps <- lapply(seq_len(n_vars), function(i) {
    ggplot(
      data.frame(t   = t_obs,
                 obs = obs_data[, i],
                 fit = solver$y[sim_obs_idx, i]),
      aes(x = t)
    ) +
      geom_point(aes(y = obs), alpha = 0.4, size = 1) +
      geom_line(aes(y = fit), color = "steelblue", linewidth = 1) +
      labs(title = var_names[i], x = "Time", y = "State") +
      theme_minimal()
  })

  pu <- lapply(seq_len(n_vars), function(i) {
    ggplot(
      data.frame(t = t_sim, u = solver$u[, i]),
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

  # 4. Local sensitivity at true generating parameters -----------------------
  cat("\n--- Local sensitivity at true generating parameters ---\n")
  cat("    (central-difference Jacobian; deSolve integrator)\n")
  sens_cfg <- make_sens_config(cfg)
  run_local_sensitivity(sens_cfg, verbose = TRUE)

  invisible(solver)
}


# =============================================================================
# Specialised examples (test specific OdeSystemSolver capabilities)
# =============================================================================

# Missing-data example --------------------------------------------------------
# Demonstrates that OdeSystemSolver bridges observation gaps by relying on
# the physics to propagate the state through the missing-data window.
run_missing_data_example <- function(lambda = 1) {
  cat("\n=== OdeSystemSolver: missing-data (Lotka-Volterra) ===\n")

  params    <- ODE_CONFIGS$lv$params
  y0        <- c(10, 10)
  times_obs <- seq(0, 40, by = 0.2)
  times_sim <- sort(unique(c(seq(0, 40, 0.05), times_obs)))

  y_true  <- euler_solve(y0, times_sim, lv_rhs, params)
  obs_idx <- which(times_sim %in% times_obs)
  set.seed(123)
  obs_data <- y_true[obs_idx, ] +
    matrix(rnorm(length(obs_idx) * 2, 0, 1.5), length(obs_idx), 2)

  gap_idx <- which(times_obs >= 10 & times_obs <= 25)
  obs_data[gap_idx, 2] <- NA
  cat(sprintf("  Introduced %d NAs (predator, t in [10, 25]).\n",
              length(gap_idx)))

  solver <- OdeSystemSolver$new(
    func_rhs   = lv_rhs,
    times_sim  = times_sim,
    obs_times  = times_obs,
    obs_values = obs_data,
    params     = params,
    lambda     = lambda
  )
  solver$optimize(y0 = y0, max_iter = 200)

  df <- data.frame(
    time     = times_sim,
    prey_fit = solver$y[, 1],
    pred_fit = solver$y[, 2],
    prey_obs = solver$observations_mapped[, 1],
    pred_obs = solver$observations_mapped[, 2]
  )
  p1 <- ggplot(df, aes(x = time)) +
    geom_point(aes(y = pred_obs), color = "blue",
               size = 2, alpha = 0.5) +
    geom_line(aes(y = pred_fit), color = "blue", linewidth = 1) +
    annotate("rect", xmin = 10, xmax = 25,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "gray") +
    annotate("text", x = 17.5, y = 5, label = "Missing data") +
    labs(title    = "Predator with gap (t = 10 to 25)",
         subtitle = paste("Points = available data",
                          " |  Line = physics-informed fit"),
         y = "Population") +
    theme_minimal()
  print(p1)
  invisible(solver)
}


# Discovery example -----------------------------------------------------------
# The true dynamics contain a hidden external forcing 0.2*sin(t); the solver
# recovers u(t) without any knowledge of the hidden term.
run_discovery_example <- function(lambda = 0.01) {
  cat("\n=== OdeSystemSolver: dynamics-discovery example ===\n")

  k          <- 6e-2
  simple_rhs <- function(y, t, p) -k * y^3   # assumed physics

  times_sim <- seq(0.01, 10, by = 0.01)
  n_steps   <- length(times_sim)

  # Truth: dy/dt = -k*y^(1/3) + 0.2*sin(t)
  y_true    <- numeric(n_steps)
  y_true[1] <- 10
  for (i in seq_len(n_steps - 1))
    y_true[i + 1] <- y_true[i] +
      0.01 * (-k * y_true[i]^(1 / 3) + 0.2 * sin(times_sim[i]))

  set.seed(123)
  obs_data <- matrix(y_true + rnorm(n_steps, 0, 0.3), ncol = 1)

  solver <- OdeSystemSolver$new(
    func_rhs   = simple_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = obs_data,
    params     = list(),
    lambda     = lambda
  )
  solver$optimize(y0 = y_true[1], max_iter = 100)

  df <- data.frame(
    time          = times_sim,
    true_y        = y_true,
    fitted_y      = solver$y[, 1],
    hidden_force  = 0.2 * sin(times_sim),
    recovered_u   = solver$u[, 1]
  )
  p1 <- ggplot(df, aes(x = time)) +
    geom_point(aes(y = true_y),  alpha = 0.3, size = 0.8) +
    geom_line(aes(y = fitted_y), color = "red") +
    labs(title    = "State fit (red) vs truth (dots)",
         subtitle = "Solver uses u(t) to bridge the model-reality gap") +
    theme_minimal()
  p2 <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = hidden_force), linetype = "dashed") +
    geom_line(aes(y = recovered_u),  color = "blue") +
    labs(title    = "Force recovery",
         subtitle = "Dashed: 0.2*sin(t)  |  Blue: recovered u(t)",
         y = "Force") +
    theme_minimal()
  grid.arrange(p1, p2, ncol = 1)
  invisible(solver)
}


# SSE surface example ---------------------------------------------------------
# Visualises the cost landscape in (alpha, beta) space with u = 0, so the
# surface reflects pure parametric misfit.
run_sse_surface_example <- function() {
  cat("\n=== OdeSystemSolver: SSE surface (Lotka-Volterra) ===\n")

  params    <- ODE_CONFIGS$lv$params
  y0        <- c(10, 10)
  times_obs <- seq(0, 40, by = 0.2)
  times_sim <- sort(unique(c(seq(0, 40, 0.01), times_obs)))

  y_true   <- euler_solve(y0, times_sim, lv_rhs, params)
  obs_data <- y_true[times_sim %in% times_obs, ]
  set.seed(123)
  obs_data <- obs_data +
    matrix(rnorm(length(obs_data), 0, 1.5), nrow(obs_data), 2)

  solver <- OdeSystemSolver$new(
    func_rhs   = lv_rhs,
    obs_times  = times_obs,
    times_sim  = times_sim,
    obs_values = obs_data,
    params     = params,
    lambda     = 0
  )

  alpha_seq    <- seq(0.8, 1.4, length.out = 25)
  beta_seq     <- seq(0.2, 0.6, length.out = 25)
  u0_flat      <- rep(0, solver$n_steps * solver$n_vars)
  cost_surface <- matrix(NA, length(alpha_seq), length(beta_seq))

  cat("  Computing surface...")
  for (i in seq_along(alpha_seq))
    for (j in seq_along(beta_seq)) {
      solver$params$alpha <- alpha_seq[i]
      solver$params$beta  <- beta_seq[j]
      cost_surface[i, j]  <- solver$cost_function(u0_flat, y0)
    }
  cat(" done.\n")

  df_surf <- melt(cost_surface)
  colnames(df_surf) <- c("alpha_idx", "beta_idx", "cost")
  df_surf$alpha <- alpha_seq[df_surf$alpha_idx]
  df_surf$beta  <- beta_seq[df_surf$beta_idx]

  p1 <- ggplot(df_surf, aes(x = alpha, y = beta, z = cost)) +
    geom_contour_filled(bins = 15) +
    annotate("point", x = 1.1, y = 0.4,
             color = "red", size = 3, shape = 4) +
    labs(title    = "SSE landscape (u = 0)",
         subtitle = "Red X = true parameters",
         fill = "Cost") +
    theme_minimal()
  print(p1)

  p2 <- plot_ly(x = beta_seq, y = alpha_seq, z = cost_surface) |>
    add_surface() |>
    layout(
      title = "Cost surface (u = 0)",
      scene = list(xaxis = list(title = "Beta"),
                   yaxis = list(title = "Alpha"),
                   zaxis = list(title = "Cost"))
    )
  print(p2)
  invisible(solver)
}


# =============================================================================
# Example calls
# =============================================================================
# run_example(ODE_CONFIGS$lv,    lambda = 1e5)
# run_example(ODE_CONFIGS$sb,    lambda = 0.1)
# run_example(ODE_CONFIGS$decay, lambda = 1e2)
#
# run_missing_data_example()
# run_discovery_example()
# run_sse_surface_example()
