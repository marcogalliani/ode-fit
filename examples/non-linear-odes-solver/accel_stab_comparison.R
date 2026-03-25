source("src/solvers/general_ode_system_solver.R")
source("src/utils/trace_optimisation.R")
source("examples/accel_stab_data.R")

# --- Parameters (pre-fitted starting point for inner solver) ----------------
params <- list(
  k_ref = 8, E = 99489.9, m = 0, n = 8.472189, C0 = 98.05758,
  R = 8.314, T_ref = T_ref, T_vec = T_vec_Kelvin
)
y0 <- rep(params$C0, times = ncol(obs_matrix))
times_grid <- as_times_grid(dt = 0.001)

# --- Inner solver -----------------------------------------------------------
solver <- OdeSystemSolver$new(
  func_rhs  = sb_kref_rhs,
  times_sim = times_grid,
  obs_times = df_clean$time,
  obs_values = obs_matrix,
  params    = params,
  lambda    = 1e2
)
solver$optimize(y0 = y0, max_iter = 100)

# --- Visualization ----------------------------------------------------------
vis <- plot_accelstab_fit(solver$y, times_grid, obs_matrix, df_clean, df_val)

# Residual forcing u(t)
df_diagn <- as.data.frame(solver$u)
colnames(df_diagn) <- colnames(obs_matrix)
df_diagn$time <- times_grid
df_diagn_long <- melt(df_diagn, id.vars = "time",
                      variable.name = "Celsius", value.name = "Conc")
p_diagn <- ggplot() +
  geom_line(data = df_diagn_long,
            aes(x = time, y = Conc, color = Celsius), size = 1) +
  theme_minimal()
print(p_diagn)

# --- AccelStab comparison ---------------------------------------------------
run_accelstab_baseline(vis$df_fit_long)

# --- Optimisation trace -----------------------------------------------------
trace_results <- trace_optimization(solver, y0 = y0, max_iter = 50)
plot(df_clean$time, obs_matrix[, 1], col = "black", pch = 16, cex = 0.5,
     main = "Optimization Progression", xlab = "Time", ylab = "Y",
     ylim = c(80, 110))
history <- trace_results$history
n_snaps <- length(history)
for (i in seq_along(history)) {
  alpha_val <- seq(0.1, 1, length.out = n_snaps)[i]
  lines(times_grid, history[[i]][, 1],
        col = rgb(0, 0, 1, alpha = alpha_val), lwd = 1)
}
lines(times_grid,
      solver$solve_state(
        matrix(trace_results$res$par, solver$n_steps), y0
      )[, 1],
      col = "red", lwd = 2)
