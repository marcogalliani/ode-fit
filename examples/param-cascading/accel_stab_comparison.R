source("src/solvers/inverse-solvers/load_inverse_solvers.R")
source("examples/accel_stab_data.R")

# --- Parameters and grid ----------------------------------------------------
params       <- as_default_params
y0           <- as_default_y0
times_grid   <- as_times_grid(dt = 0.01)

param_scales <- list(E = 10000, n = 1, C0 = 1, k_ref = 1)
fixed_params <- params[!names(params) %in% names(param_scales)]

# --- Cascading solver -------------------------------------------------------
cascading_solver <- CascadingOdeSolver$new(
  func_rhs     = sb_kref_rhs,
  obs_times    = df_clean$time,
  times_sim    = times_grid,
  obs_values   = obs_matrix,
  init_state   = function(p) as.numeric(y0),
  fixed_params = fixed_params,
  lambda       = 1e-3,
  param_scales = param_scales
)

init_params <- params[names(param_scales)]
l_bounds <- c(E = 10000, n = 0.1,
              C0 = max(obs_matrix, na.rm = TRUE) + 0.1, k_ref = 1e-4)
u_bounds <- c(E = 300000, n = 20, C0 = 120, k_ref = 100)

result <- cascading_solver$optimize_parameters(
  init_theta_physical = unlist(init_params),
  param_names         = names(init_params),
  lower_phys          = l_bounds,
  upper_phys          = u_bounds
)

# --- Visualization ----------------------------------------------------------
vis <- plot_accelstab_fit(cascading_solver$last_solver$y,
                          times_grid, obs_matrix, df_clean, df_val)

# --- AccelStab comparison ---------------------------------------------------
run_accelstab_baseline(vis$df_fit_long)
