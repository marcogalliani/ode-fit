library(deSolve)
library(sensitivity)

source("examples/ode_models.R")

# Select configuration ----
# Change this line to switch models.  Available keys: decay, lv, sb,
# sb_aspack, sb_asymptote.  All configurations are defined in
# examples/ode_models.R (ODE_CONFIGS list).
active <- make_sens_config(ODE_CONFIGS$lv)

# Helpers ----
evaluate_trajectory <- function(pars, model, init_state, t_grid,
                                rhs_args = list()) {
  ode_args <- c(
    list(y = init_state, times = t_grid, func = model, parms = pars),
    rhs_args
  )
  as.data.frame(do.call(ode, ode_args))
}

# Collapse a multi-state trajectory to a scalar.
# state_aggregation: how to combine multiple state columns into one series.
# metric:            which summary of that series to return.
extract_metric <- function(trajectory, state_name,
                           metric = c("final", "auc", "max", "mean"),
                           state_aggregation = c("mean", "sum", "first")) {
  metric            <- match.arg(metric)
  state_aggregation <- match.arg(state_aggregation)

  if (!all(state_name %in% names(trajectory)))
    stop("Some requested states are not present in the trajectory output")

  y_mat <- as.matrix(trajectory[, state_name, drop = FALSE])
  y <- switch(state_aggregation,
    mean  = rowMeans(y_mat),
    sum   = rowSums(y_mat),
    first = y_mat[, 1]
  )
  t <- trajectory$time

  switch(metric,
    final = y[length(y)],
    auc   = sum((y[-1] + y[-length(y)]) * diff(t) / 2),
    max   = max(y),
    mean  = mean(y)
  )
}

# Names of the state variables to aggregate (all state names by default).
active_states <- names(active$state)

# Evaluate a scalar metric for every row of a parameter design matrix X.
# Columns of X must correspond to active$uncertain_bounds keys.
evaluate_parameter_sets <- function(X, active,
                                    metric = "final",
                                    state_aggregation = "mean") {
  uncertain_names <- names(active$uncertain_bounds)
  vapply(seq_len(nrow(X)), function(i) {
    pars <- active$params
    pars[uncertain_names] <- as.numeric(X[i, uncertain_names])
    traj <- evaluate_trajectory(pars, active$model, active$state,
                                active$t_grid, active$rhs_args)
    extract_metric(traj, active_states, metric, state_aggregation)
  }, numeric(1))
}

# Evaluate the full trajectory for every row of X.
# Returns a matrix with nrow(X) rows and length(t_grid) columns.
evaluate_trajectory_sets <- function(X, active,
                                     state_aggregation = "mean") {
  uncertain_names <- names(active$uncertain_bounds)
  do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    pars <- active$params
    pars[uncertain_names] <- as.numeric(X[i, uncertain_names])
    traj <- evaluate_trajectory(pars, active$model, active$state,
                                active$t_grid, active$rhs_args)
    y_mat <- as.matrix(traj[, active_states, drop = FALSE])
    switch(state_aggregation,
      mean  = rowMeans(y_mat),
      sum   = rowSums(y_mat),
      first = y_mat[, 1]
    )
  }))
}

sample_uniform <- function(n, bounds) {
  as.data.frame(lapply(bounds, function(b) runif(n, b[1], b[2])))
}

# Scalar Sobol analysis ----
run_global_sobol <- function(active, n_samples = 500,
                             metric = "final", state_aggregation = "mean",
                             nboot = 100, conf = 0.95, seed = 123) {
  set.seed(seed)
  X1 <- sample_uniform(n_samples, active$uncertain_bounds)
  X2 <- sample_uniform(n_samples, active$uncertain_bounds)

  sob <- soboljansen(model = NULL, X1 = X1, X2 = X2, nboot = nboot, conf = conf)
  Y   <- evaluate_parameter_sets(sob$X, active, metric, state_aggregation)
  sob <- tell(sob, Y)

  list(sobol = sob, metric = metric,
       state_aggregation = state_aggregation, n_samples = n_samples)
}

# Time-resolved Sobol analysis ----
# Computes first- and total-order Sobol indices at every time step using a
# single set of ODE evaluations (no extra solves per time point).
# Returns matrices S and T of shape (n_times x n_pars).
run_global_sobol_timeseries <- function(active, n_samples = 300,
                                        state_aggregation = "mean",
                                        nboot = 0, seed = 123) {
  set.seed(seed)
  X1 <- sample_uniform(n_samples, active$uncertain_bounds)
  X2 <- sample_uniform(n_samples, active$uncertain_bounds)

  sob_template <- soboljansen(model = NULL, X1 = X1, X2 = X2, nboot = nboot)

  # Evaluate full trajectories for all design points in one pass
  Y_mat <- evaluate_trajectory_sets(sob_template$X, active, state_aggregation)
  # Y_mat: nrow = N*(k+2),  ncol = length(t_grid)

  n_times   <- ncol(Y_mat)
  par_names <- names(active$uncertain_bounds)
  n_pars    <- length(par_names)

  S_time <- matrix(NA, n_times, n_pars, dimnames = list(NULL, par_names))
  T_time <- matrix(NA, n_times, n_pars, dimnames = list(NULL, par_names))

  for (j in seq_len(n_times)) {
    sob_j <- sob_template
    sob_j <- tell(sob_j, Y_mat[, j])
    S_time[j, ] <- sob_j$S$original
    T_time[j, ] <- sob_j$T$original
  }

  list(S = S_time, T = T_time,
       t_grid = active$t_grid, par_names = par_names,
       n_samples = n_samples)
}

# Run scalar analysis ----
result <- run_global_sobol(
  active            = active,
  n_samples         = 300,
  metric            = "final",
  state_aggregation = "mean",
  nboot             = 100
)

cat("\n--- Global Sensitivity (Sobol-Jansen, scalar metric) ---\n")
cat("Metric       :", result$metric, "|",
    "State aggr.  :", result$state_aggregation, "\n")
cat("Samples      :", result$n_samples, "\n\n")

cat("First-order Sobol indices (S):\n")
print(result$sobol$S)

cat("\nTotal-order Sobol indices (T):\n")
print(result$sobol$T)

interaction_share <- result$sobol$T$original - result$sobol$S$original
names(interaction_share) <- rownames(result$sobol$S)
cat("\nInteraction / nonlinear contribution (T - S):\n")
print(round(interaction_share, 4))

cat("\nInterpretation:\n")
cat("  High S       : large direct (main-effect) contribution to output variance.\n")
cat("  High T       : large total contribution (main + interactions).\n")
cat("  Large (T-S)  : strong interactions or nonlinear parameter coupling.\n")
cat("  Sobol indices are variance fractions, not correlation coefficients.\n")

if (interactive()) plot(result$sobol)

# Run time-resolved analysis ----
ts_result <- run_global_sobol_timeseries(
  active            = active,
  n_samples         = 300,
  state_aggregation = "mean"
)

par_colors <- seq_len(length(ts_result$par_names))

matplot(ts_result$t_grid, ts_result$S,
        type = "l", lty = 1, col = par_colors, ylim = c(0, 1),
        xlab = "Time", ylab = "First-order Sobol index (S)",
        main = "Time-resolved global sensitivity (first-order)")
legend("topright", legend = ts_result$par_names,
       col = par_colors, lty = 1, bty = "n")

matplot(ts_result$t_grid, ts_result$T,
        type = "l", lty = 2, col = par_colors, ylim = c(0, 1),
        xlab = "Time", ylab = "Total-order Sobol index (T)",
        main = "Time-resolved global sensitivity (total-order)")
legend("topright", legend = ts_result$par_names,
       col = par_colors, lty = 2, bty = "n")
