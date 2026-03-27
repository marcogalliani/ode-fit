# =============================================================================
# examples/plot_utils.R
#
# Shared plotting helpers for solver example scripts.
# =============================================================================

library(ggplot2)
library(gridExtra)
library(reshape2)


# State-fit plot: observed vs fitted for each variable.
plot_state_fit <- function(y_fit, obs_data, t_obs, t_sim, var_names) {
  # check assumption
  stopifnot(nrow(y_fit) == length(t_sim))
  stopifnot(nrow(obs_data) == length(t_obs))
  stopifnot(ncol(y_fit) == length(var_names))
  stopifnot(ncol(obs_data) == length(var_names))
  n_vars <- length(var_names)

  ps <- lapply(seq_len(n_vars), function(i) {
    df_fit <- data.frame(t = t_sim, fit = y_fit[, i])
    df_obs <- data.frame(t = t_obs, obs = obs_data[, i])

    ggplot(mapping = aes(x = t)) +
      geom_point(data = df_obs, aes(y = .data$obs), alpha = 0.5, size = 1) +
      geom_line(data = df_fit, aes(y = .data$fit), color = "steelblue", linewidth = 1) +
      labs(title = var_names[i], x = "Time", y = "State") +
      theme_minimal()
  })

  do.call(grid.arrange,
          c(ps, list(ncol = min(n_vars, 3), top = "State fit")))
  invisible(ps)
}


# Residual forcing u(t) plot for each variable.
plot_residual_forcing <- function(u_mat, t_sim, var_names) {
  n_vars <- length(var_names)

  pu <- lapply(seq_len(n_vars), function(i) {
    ggplot(
      data.frame(t = t_sim, u = u_mat[, i]),
      aes(x = t, y = u)
    ) +
      geom_line(color = "tomato", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(title = paste0(var_names[i], "  u(t)"), x = "Time", y = "u") +
      theme_minimal()
  })

  do.call(grid.arrange,
          c(pu, list(ncol = min(n_vars, 3), top = "Residual forcing u(t)")))
  invisible(pu)
}


# AccelStab fit + validation plot (multi-temperature).
plot_accelstab_fit <- function(y_fit, times_grid, obs_matrix, df_clean,
                               df_val, title = "Sestak-Berggren Model Fit") {
  df_fit <- as.data.frame(y_fit)
  colnames(df_fit) <- colnames(obs_matrix)
  df_fit$time <- times_grid
  df_fit_long <- melt(df_fit, id.vars = "time",
                      variable.name = "Celsius", value.name = "Conc")

  df_obs <- as.data.frame(obs_matrix)
  colnames(df_obs) <- colnames(obs_matrix)
  df_obs$time <- df_clean$time
  df_obs_long <- melt(df_obs, id.vars = "time",
                      variable.name = "Celsius", value.name = "Conc")

  df_val_long <- melt(df_val, id.vars = "time",
                      variable.name = "Celsius", value.name = "Conc")

  p <- ggplot() +
    geom_point(data = df_obs_long,
               aes(x = time, y = Conc, color = Celsius),
               size = 2, alpha = 0.6) +
    geom_line(data = df_fit_long,
              aes(x = time, y = Conc, color = Celsius), size = 1) +
    lims(y = c(0, 100)) +
    labs(title = title, y = "Concentration") +
    geom_point(data = df_val_long,
               aes(x = time, y = Conc, color = Celsius),
               size = 2, pch = 1) +
    theme_minimal()

  print(p)
  invisible(list(plot = p, df_fit_long = df_fit_long))
}
