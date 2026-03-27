source("src/solvers/inverse-solvers/load_inverse_solvers.R")
source("examples/accel_stab_data.R")

# --- Parameters and grid ----------------------------------------------------
params       <- as_default_params
params$m     <- 1
y0           <- as_default_y0
times_grid   <- as_times_grid(dt = 0.05)

param_scales <- list(E = 10000, n = 1, C0 = 1, k_ref = 1, m = 1)
fixed_params <- params[!names(params) %in% names(param_scales)]

# --- Tracking solver --------------------------------------------------------
tracking_solver <- TrackingOdeSolver$new(
  func_rhs     = sb_kref_rhs,
  obs_times    = df_clean$time,
  times_sim    = times_grid,
  obs_values   = obs_matrix,
  init_state   = function(p) as.numeric(rep(p$C0,ncol(obs_matrix))),
  fixed_params = fixed_params,
  lambda       = 1e-1,
  param_scales = param_scales
)

init_params <- params[names(param_scales)]
init_params$C0 <- mean(obs_matrix[1,], na.rm=T)

lb_params <- rep(0, length(init_params))
names(lb_params) <- names(init_params)


result <- tracking_solver$optimize_parameters(
  init_theta_physical = unlist(init_params),
  param_names         = names(init_params),
  lower_phys = lb_params
)

# --- Visualization ----------------------------------------------------------
vis <- plot_accelstab_fit(tracking_solver$last_solver$y,
                          times_grid, obs_matrix, df_clean, df_val)

# --- AccelStab comparison ---------------------------------------------------
run_accelstab_baseline(vis$df_fit_long)
