library(ggplot2)
library(gridExtra)
library(reshape2)

source("src/solvers/parameter_cascading.R")
source("src/utils/trace_optimisation.R")

# --- Shared utility: forward Euler integrator (consistent with OdeSystemSolver) ---
euler_solve <- function(y0, times, func, params) {
  n      <- length(times)
  n_vars <- length(y0)
  y      <- matrix(0, n, n_vars)
  y[1, ] <- y0
  dt_vec <- c(diff(times), 0)
  for (i in 1:(n - 1)) {
    dy         <- func(y[i, ], times[i], params)
    y[i + 1, ] <- y[i, ] + dt_vec[i] * dy
  }
  return(y)
}


# EXAMPLE 1: Exponential Decay -- estimate k ----
run_decay_example <- function(init_params = NULL, lambda = NULL) {
  cat("\n=== Example 1: Exponential Decay - Estimating k ===\n")

  decay_rhs <- function(y, t, p) -p$k * y^3

  k_true  <- 0.5
  y0_true <- c(5.0)

  times_obs  <- seq(0, 10, by = 0.1)
  times_fine <- seq(0, 10, by = 0.01)
  times_sim  <- times_fine   # pass clean grid; solver will merge obs_times internally

  # Ground truth + noise
  y_sim_fine <- euler_solve(y0_true, times_fine, decay_rhs, list(k = k_true))
  obs_idx    <- match(round(times_obs, 10), round(times_fine, 10))
  set.seed(42)
  obs_data <- y_sim_fine[obs_idx, , drop = FALSE] +
    matrix(rnorm(length(obs_idx), 0, 0.1), length(obs_idx), 1)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = times_obs,
    obs_values   = obs_data,
    fixed_params = list(),
    lambda       = ifelse(is.null(lambda), 1e0, lambda),
    param_scales = list(k = 1.0)
  )

  init_theta_physical <- if (!is.null(init_params)) init_params else c(k = 2.0)

  result <- cascading$optimize_parameters(
    init_theta_physical = init_theta_physical,
    param_names         = "k",
    lower_phys          = c(k = 0.01),
    upper_phys          = c(k = 10.0)
  )

  cat(sprintf("\n  True k = %.3f | Estimated k = %.3f\n", k_true, result["k"]))

  s       <- cascading$last_solver
  df_fit  <- data.frame(Time = s$times_sim, Fit   = s$y[, 1])
  df_obs  <- data.frame(Time = times_obs,   Obs   = obs_data[, 1])
  df_true <- data.frame(Time = times_fine,  Truth = y_sim_fine[, 1])

  p_state <- ggplot() +
    geom_line(data = df_true, aes(x = Time, y = Truth),
              linetype = "dashed", color = "gray50") +
    geom_point(data = df_obs, aes(x = Time, y = Obs),
               color = "black", alpha = 0.7) +
    geom_line(data = df_fit, aes(x = Time, y = Fit),
              color = "steelblue", linewidth = 1) +
    labs(
      title    = sprintf("Exponential Decay: k_true=%.2f  k_est=%.3f", k_true, result["k"]),
      subtitle = "Dashed: true trajectory  |  Points: noisy obs  |  Blue: cascading fit",
      y = "y(t)"
    ) +
    theme_minimal()

  df_u      <- data.frame(Time = s$times_sim, u = s$u[, 1])
  p_forcing <- ggplot(df_u, aes(x = Time, y = u)) +
    geom_line(color = "tomato") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title    = "Residual Forcing u(t)",
      subtitle = "Small |u(t)| throughout indicates a good parameter estimate"
    ) +
    theme_minimal()

  grid.arrange(p_state, p_forcing, ncol = 1)
  plot_param_trace(cascading, true_params = list(k = k_true))
  return(result)
}


# EXAMPLE 2: Lotka-Volterra -- estimate alpha and gamma ----
run_lotka_volterra_example <- function(init_params = NULL, lambda = NULL) {
  cat("\n=== Example 2: Lotka-Volterra - Estimating alpha and gamma ===\n")

  lotka_volterra <- function(y, t, p) {
    c(p$alpha * y[1] - p$beta  * y[1] * y[2],
      p$delta * y[1] * y[2] - p$gamma * y[2])
  }

  p_true  <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  y0_true <- c(10, 10)

  times_obs  <- seq(0, 30, by = 0.5)
  times_fine <- seq(0, 30, by = 0.1)
  times_sim  <- times_fine   # pass clean grid; solver will merge obs_times internally

  y_sim_fine <- euler_solve(y0_true, times_fine, lotka_volterra, p_true)
  obs_idx    <- match(round(times_obs, 10), round(times_fine, 10))
  set.seed(123)
  obs_data <- y_sim_fine[obs_idx, ] +
    matrix(rnorm(length(obs_idx) * 2, 0, 0.8), length(obs_idx), 2)

  # Fix beta and delta; estimate alpha and gamma
  fixed_params <- list(beta = p_true$beta, delta = p_true$delta)
  param_scales <- list(alpha = 1.0, gamma = 1.0)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = lotka_volterra,
    times_sim    = times_sim,
    obs_times    = times_obs,
    obs_values   = obs_data,
    fixed_params = fixed_params,
    lambda       = ifelse(is.null(lambda), 1e1, lambda),
    param_scales = param_scales
  )

  init_theta_physical <- if (!is.null(init_params)) init_params else c(alpha = 0.8, gamma = 0.6)

  result <- cascading$optimize_parameters(
    init_theta_physical = init_theta_physical,
    param_names         = c("alpha", "gamma"),
    lower_phys          = c(alpha = 0.1, gamma = 0.05),
    upper_phys          = c(alpha = 10, gamma = 10)
  )

  cat(sprintf("\n  True:  alpha=%.2f, gamma=%.2f\n", p_true$alpha, p_true$gamma))
  cat(sprintf("  Est:   alpha=%.3f, gamma=%.3f\n", result["alpha"], result["gamma"]))

  s      <- cascading$last_solver
  df_sim <- data.frame(Time = s$times_sim, Prey = s$y[, 1], Pred = s$y[, 2])
  df_obs <- data.frame(Time = times_obs,
                       Prey_Obs = obs_data[, 1], Pred_Obs = obs_data[, 2])

  p1 <- ggplot() +
    geom_point(data = df_obs, aes(x = Time, y = Prey_Obs),
               color = "forestgreen", alpha = 0.5) +
    geom_line(data = df_sim, aes(x = Time, y = Prey),
              color = "forestgreen", linewidth = 1) +
    labs(
      title = sprintf("Prey  (alpha: %.2f -> %.3f)", p_true$alpha, result["alpha"]),
      y = "Population"
    ) +
    theme_minimal()

  p2 <- ggplot() +
    geom_point(data = df_obs, aes(x = Time, y = Pred_Obs),
               color = "firebrick", alpha = 0.5) +
    geom_line(data = df_sim, aes(x = Time, y = Pred),
              color = "firebrick", linewidth = 1) +
    labs(
      title = sprintf("Predator  (gamma: %.2f -> %.3f)", p_true$gamma, result["gamma"]),
      y = "Population"
    ) +
    theme_minimal()

  df_u <- data.frame(Time = s$times_sim, u_Prey = s$u[, 1], u_Pred = s$u[, 2])
  p3 <- ggplot(df_u, aes(x = Time)) +
    geom_line(aes(y = u_Prey), color = "forestgreen", alpha = 0.8) +
    geom_line(aes(y = u_Pred), color = "firebrick",   alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title    = "Residual Forcing u(t)",
      subtitle = "Green: prey  |  Red: predator"
    ) +
    theme_minimal()

  grid.arrange(p1, p2, p3, ncol = 1)
  plot_param_trace(cascading,
                   true_params = list(alpha = p_true$alpha, gamma = p_true$gamma))
  return(result)
}


# EXAMPLE 3: Sestak-Berggren -- estimate E, n, m ----
run_sb_example <- function(init_params = c(E = 50000, n = 3.0, m = 0.5), lambda = 1e1) {
  cat("\n=== Example 3: Sestak-Berggren - Estimating E, n, m ===\n")

  sb_rhs <- function(x_vec, t, p) {
    eps    <- 1e-6
    x_safe <- pmax(pmin(x_vec, p$C0 - eps), eps)
    k_T    <- p$A * exp(-p$E / (p$R * p$T_vec))
    denom  <- p$C0^(p$m + p$n - 1)
    -(k_T / denom) * ((p$C0 - x_safe)^p$m) * (x_safe^p$n)
  }

  T_vec  <- c(400, 420, 440)
  p_true <- list(A = 1e6, E = 40000, m = 0.5, n = 4.0,
                 C0 = 1.0, R = 8.314, T_vec = T_vec)
  y0_true <- rep(0.99, length(T_vec))

  times_obs  <- seq(0, 40, by = 0.5)
  times_fine <- seq(0, 40, by = 0.1)
  times_sim  <- times_fine   # pass clean grid; solver will merge obs_times internally

  y_sim_fine <- euler_solve(y0_true, times_fine, sb_rhs, p_true)
  obs_idx    <- match(round(times_obs, 10), round(times_fine, 10))
  set.seed(123)
  obs_data <- y_sim_fine[obs_idx, ] +
    matrix(rnorm(length(obs_idx) * length(T_vec), 0, 0.02),
           length(obs_idx), length(T_vec))

  # Fix A, C0, R, T_vec; estimate E, n, m
  fixed_params <- list(A = p_true$A, C0 = p_true$C0,
                       R = p_true$R, T_vec = T_vec)
  param_scales <- list(E = 10000, n = 1, m = 1)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = sb_rhs,
    times_sim    = times_sim,
    obs_times    = times_obs,
    obs_values   = obs_data,
    fixed_params = fixed_params,
    lambda       = lambda,
    param_scales = param_scales
  )

  result <- cascading$optimize_parameters(
    init_theta_physical = init_params,
    param_names         = c("E", "n", "m"),
    lower_phys          = c(E = 10000, n = 0.1, m = 0.1),
    upper_phys          = c(E = 300000, n = 20,  m = 20)
  )

  cat(sprintf("\n  True:  E=%.0f,  n=%.1f,  m=%.2f\n", p_true$E, p_true$n, p_true$m))
  cat(sprintf("  Est:   E=%.0f, n=%.3f, m=%.3f\n", result["E"], result["n"], result["m"]))

  s           <- cascading$last_solver
  temp_labels <- paste0("T=", T_vec, "K")

  df_sim <- data.frame(
    Time = rep(s$times_sim, length(T_vec)),
    Conc = as.vector(s$y),
    Temp = factor(rep(temp_labels, each = s$n_steps), levels = temp_labels)
  )
  df_obs <- data.frame(
    Time = rep(times_obs, length(T_vec)),
    Conc = as.vector(obs_data),
    Temp = factor(rep(temp_labels, each = length(times_obs)), levels = temp_labels)
  )

  colors <- c("firebrick", "steelblue", "forestgreen")

  p_fit <- ggplot() +
    geom_point(data = df_obs, aes(x = Time, y = Conc, color = Temp), alpha = 0.5) +
    geom_line(data = df_sim,  aes(x = Time, y = Conc, color = Temp), linewidth = 1) +
    scale_color_manual(values = colors) +
    labs(
      title = sprintf(
        "Sestak-Berggren: E=%.0f (true %.0f)  n=%.2f (true %.1f)  m=%.2f (true %.2f)",
        result["E"], p_true$E, result["n"], p_true$n, result["m"], p_true$m
      ),
      y = "Concentration", color = "Temperature"
    ) +
    theme_minimal()

  df_u <- data.frame(
    Time = rep(s$times_sim, length(T_vec)),
    u    = as.vector(s$u),
    Temp = factor(rep(temp_labels, each = s$n_steps), levels = temp_labels)
  )

  p_u <- ggplot(df_u, aes(x = Time, y = u, color = Temp)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = colors) +
    labs(title = "Residual Forcing u(t)", y = "u(t)", color = "Temperature") +
    theme_minimal()

  grid.arrange(p_fit, p_u, ncol = 1)
  plot_param_trace(cascading,
                   true_params = list(E = p_true$E, n = p_true$n, m = p_true$m))
  return(result)
}


# --- Run all examples --------------------------------------------------------
# source("examples/param-cascading/ode_params_examples.R")
# run_decay_example()
# run_lotka_volterra_example()
# run_sb_example()
