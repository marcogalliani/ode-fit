# =============================================================================
# examples/ode_models.R
#
# Single source of truth for ODE models shared by:
#   examples/non-linear-odes-solver/
#   examples/param-cascading/
#   examples/sensitivity-analysis/
#
# Tests define their own copies for self-containment and must not source this.
#
# MODEL FORMAT
# ------------
# All RHS functions use the solver-native signature:
#
#   rhs(y, t, p)
#
# where  y  is a numeric vector (state),  t  is a scalar (time), and
# p  is a named list of parameters.  This signature is directly compatible
# with OdeSystemSolver and CascadingOdeSolver.
#
# DESOLVE WRAPPER
# ---------------
# Use  as_desolve(rhs)  to obtain a deSolve / FME-compatible wrapper.
# Extra named arguments passed to the wrapper via  ...  (e.g. T_vec) are
# merged into p before calling rhs, so vector-valued parameters such as
# T_vec can be supplied through the rhs_args mechanism used by deSolve.
#
# CANONICAL CONFIGURATIONS
# ------------------------
# ODE_CONFIGS  is a named list.  Each entry bundles:
#   rhs              solver-format RHS
#   params           full parameter list (may contain vector-valued entries)
#   y0               initial state (named numeric vector)
#   times_obs        default observation grid
#   times_sim        fine integration grid
#   active_params    parameters to vary in sensitivity analysis
#   uncertain_bounds named list of c(lo, hi) for each active parameter
#   fixed_params     parameters held fixed by CascadingOdeSolver
#   param_scales     normalisation scales for CascadingOdeSolver
#
# SENSITIVITY HELPER
# ------------------
# make_sens_config(cfg)  converts an ODE_CONFIGS entry into the dict format
# consumed by local_sensitivity_analysis.R and global_sensitivity_analysis.R.
# Scalar entries of cfg$params become the `params` numeric vector; vector-
# valued entries (e.g. T_vec) are placed in `rhs_args` so they reach rhs
# via the deSolve ... mechanism.
# =============================================================================

# === 1. Exponential decay ===
# dy/dt = -k * y
decay_rhs <- function(y, t, p) -p$k * y


# === 2. Lotka-Volterra (predator-prey) ===
# dx/dt = alpha*x - beta*x*y
# dy/dt = delta*x*y - gamma*y
# Canonical parameters: alpha=1.1, beta=0.4, delta=0.1, gamma=0.4
# Initial state: c(x=10, y=10)
lv_rhs <- function(y, t, p) {
  c(
    p$alpha * y[1] - p$beta * y[1] * y[2],
    p$delta * y[1] * y[2] - p$gamma * y[2]
  )
}


# === 3. Sestak-Berggren (multi-temperature) ===
# One state component per temperature level; T_vec has the same length as y.
# dx_i/dt = -(A * exp(-E/(R*T_i)) / C0^(m+n-1)) * (C0-x_i)^m * x_i^n
# Canonical parameters: A=1e6, E=40000, m=0.5, n=4.0, C0=1.0, R=8.314
# T_vec: c(400, 420, 440),  y0: c(0.99, 0.99, 0.99)
sb_rhs <- function(y, t, p) {
  eps <- 1e-6
  x_eff <- pmax(pmin(y, p$C0 - eps), eps)
  k_T <- p$A * exp(-p$E / (p$R * p$T_vec))
  denom <- p$C0^(p$m + p$n - 1)
  -(k_T / denom) * ((p$C0 - x_eff)^p$m) * (x_eff^p$n)
}


# === 4. Sestak-Berggren AccelStab-package variant ===
# Arrhenius referenced at T_ref, first-order-like reaction of order k3.
# dx_i/dt = exp(k1 - k2*(1/T_i - 1/T_ref)) * (1-x_i)^k3
# Canonical: k1=40, k2=8710.63, k3=3, T_ref=296.65 K
# T_vec: 273.15+c(5,20,32,37),  y0: c(0.01,0.01,0.01,0.01)
sb_aspack_rhs <- function(y, t, p) {
  eps <- 1e-6
  x_eff <- pmin(pmax(as.numeric(y), eps), 1 - eps)
  k_arr <- exp(p$k1 - p$k2 * (1 / p$T_vec - 1 / p$T_ref))
  k_arr * (1 - x_eff)^p$k3
}


# === 5. Sestak-Berggren asymptote variant ===
# Bounded between L0 and C0; Arrhenius referenced at T_ref.
# dx_i/dt = -exp(k1 - k2*(1/T_i-1/T_ref) - (m+n-1)*log(C0-L0)) *
#            (C0-x_i)^m * (x_i-L0)^n
# Canonical: k1=40, k2=8710.63, m=0.1, n=5, C0=99.99, L0=0.01, T_ref=296.65 K
# T_vec: 273.15+c(5,20,32,37),  y0: c(99.99,99.99,99.99,99.99)
sb_asymptote_rhs <- function(y, t, p) {
  eps <- 1e-12
  x_eff <- pmin(pmax(as.numeric(y), p$L0 + eps), p$C0 - eps)
  k_mod <- exp(p$k1 - p$k2 * (1 / p$T_vec - 1 / p$T_ref) -
    (p$m + p$n - 1) * log(p$C0 - p$L0))
  -k_mod * (p$C0 - x_eff)^p$m * (x_eff - p$L0)^p$n
}


# === 6. Sestak-Berggren k_ref variant (AccelStab comparison) ===
# Arrhenius referenced at T_ref with explicit rate constant k_ref.
# dx_i/dt = -(k_ref * exp(-(E/R)*(1/T_i - 1/T_ref)) / C0^(m+n-1)) *
#            (C0 - x_i)^m * x_i^n
# Used by examples/accel_stab_data.R for real-data comparisons.
sb_kref_rhs <- function(y, t, p) {
  eps <- 1e-6
  x_eff <- pmax(pmin(y, p$C0 - eps), eps)
  k_T <- p$k_ref * exp(-(p$E / p$R) * (1 / p$T_vec - 1 / p$T_ref))
  scaling <- p$C0^(p$m + p$n - 1)
  -(k_T / scaling) * (p$C0 - x_eff)^p$m * x_eff^p$n
}


# === Forward Euler integrator ===
# Consistent with OdeSystemSolver / CascadingOdeSolver.  Used to generate
# synthetic ground-truth data in examples.
euler_solve <- function(y0, times, rhs, params) {
  n <- length(times)
  n_vars <- length(y0)
  y <- matrix(0, n, n_vars)
  y[1, ] <- y0
  dt_vec <- c(diff(times), 0)
  for (i in seq_len(n - 1)) {
    y[i + 1, ] <- y[i, ] + dt_vec[i] * rhs(y[i, ], times[i], params)
  }
  y
}


# === deSolve adapter ===
# Returns a deSolve-compatible wrapper around any solver-format rhs(y, t, p).
# Named extra arguments (e.g. T_vec = c(400, 420, 440)) are merged into p.
as_desolve <- function(rhs) {
  function(t, state, parms, ...) {
    p <- c(as.list(parms), list(...))
    list(rhs(as.numeric(state), t, p))
  }
}


# === Sensitivity config builder ===
# Converts an ODE_CONFIGS entry into the dict format consumed by
# local_sensitivity_analysis.R and global_sensitivity_analysis.R.
# Scalar params  -> `params`   (named numeric vector for deSolve parms)
# Vector params  -> `rhs_args` (passed as ... to the deSolve wrapper)
make_sens_config <- function(cfg) {
  all_pars <- cfg$params
  is_scalar <- vapply(all_pars, function(x) length(x) == 1, logical(1))

  list(
    model            = as_desolve(cfg$rhs),
    params           = unlist(all_pars[is_scalar]),
    active_params    = cfg$active_params,
    state            = cfg$y0,
    t_grid           = cfg$times_obs,
    rhs_args         = all_pars[!is_scalar],
    uncertain_bounds = cfg$uncertain_bounds
  )
}


# === Synthetic data generator ===
# Generates noisy observations from an ODE_CONFIGS entry.
# Returns: list(obs_data, y_true, sd_used, obs_idx)
generate_synthetic_data <- function(cfg, noise_sd = NULL, noise_pct = 0.02,
                                    seed = 123) {
  y_true  <- euler_solve(cfg$y0, cfg$times_sim, cfg$rhs, cfg$params)
  obs_idx <- match(round(cfg$times_obs, 10), round(cfg$times_sim, 10))
  stopifnot("obs times must be a subset of sim times" = !anyNA(obs_idx))

  sd_used <- if (is.null(noise_sd))
    noise_pct * max(abs(y_true[obs_idx, , drop = FALSE]))
  else
    noise_sd

  set.seed(seed)
  n_obs  <- length(obs_idx)
  n_vars <- length(cfg$y0)
  obs_data <- y_true[obs_idx, , drop = FALSE] +
    matrix(rnorm(n_obs * n_vars, 0, sd_used), n_obs, n_vars)

  list(obs_data = obs_data, y_true = y_true, sd_used = sd_used,
       obs_idx = obs_idx)
}


# === Canonical configurations ===
ODE_CONFIGS <- list(
  decay = list(
    rhs = decay_rhs,
    params = list(k = 0.5),
    y0 = c(y = 5.0),
    times_obs = seq(0, 10, by = 0.5),
    times_sim = seq(0, 10, by = 0.1),
    active_params = "k",
    uncertain_bounds = list(k = c(0.1, 2.0)),
    fixed_params = list(),
    param_scales = list(k = 1.0)
  ),
  lv = list(
    rhs = lv_rhs,
    params = list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4),
    y0 = c(x = 10, y = 10),
    times_obs = seq(0, 5, by = 0.5),
    times_sim = seq(0, 5, by = 0.001),
    active_params = c("alpha", "beta", "delta", "gamma"),
    uncertain_bounds = list(
      alpha = c(0.9, 1.3),
      beta  = c(0.2, 0.6),
      delta = c(0.05, 0.2),
      gamma = c(0.2, 0.6)
    ),
    fixed_params = list(),
    param_scales = list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)
  ),

  # Sestak-Berggren multi-temperature.
  # Canonical parameters match examples/param-cascading/synthetic_sb_test.R
  # and examples/param-cascading/ode_params_examples.R (run_sb_example).
  # active_params / uncertain_bounds cover the parameters estimated by
  # CascadingOdeSolver (E, m, n); A is treated as known.
  sb = list(
    rhs = sb_rhs,
    params = list(
      A = 1e6, E = 40000, m = 0.5, n = 4.0,
      C0 = 1.0, R = 8.314, T_vec = c(400, 420, 440)
    ),
    y0 = c(x1 = 0.99, x2 = 0.99, x3 = 0.99),
    times_obs = seq(0, 40, by = 0.5),
    times_sim = seq(0, 40, by = 0.1),
    active_params = c("E", "m", "n"),
    uncertain_bounds = list(
      E = c(20000, 80000),
      m = c(0.1, 2.0),
      n = c(1.0, 8.0)
    ),
    fixed_params = list(
      A = 1e6, C0 = 1.0, R = 8.314,
      T_vec = c(400, 420, 440)
    ),
    # Scales cover all params that may be estimated (A included for joint runs)
    param_scales = list(A = 1e5, E = 10000, m = 1.0, n = 1.0)
  ),
  sb_aspack = list(
    rhs = sb_aspack_rhs,
    params = list(
      k1 = 40, k2 = 8710.63, k3 = 3,
      T_ref = 273.15 + 23.5,
      T_vec = 273.15 + c(5, 20, 32, 37)
    ),
    y0 = c(x1 = 0.01, x2 = 0.01, x3 = 0.01, x4 = 0.01),
    times_obs = seq(0, 3, by = 0.1),
    times_sim = seq(0, 3, by = 0.01),
    active_params = c("k1", "k2", "k3"),
    uncertain_bounds = list(
      k1 = c(35, 45),
      k2 = c(7800, 9600),
      k3 = c(2.5, 3.5)
    ),
    fixed_params = list(
      T_ref = 273.15 + 23.5,
      T_vec = 273.15 + c(5, 20, 32, 37)
    ),
    param_scales = list(k1 = 1.0, k2 = 1000.0, k3 = 1.0)
  ),
  sb_asymptote = list(
    rhs = sb_asymptote_rhs,
    params = list(
      k1 = 1, k2 = 1000, m = 0, n = 5,
      C0 = 100 - 0.01, L0 = 50,
      T_ref = 273.15 + 23.5,
      T_vec = 273.15 + c(5, 20, 32, 37)
    ),
    y0 = 100 - c(x1 = 0.01, x2 = 0.01, x3 = 0.01, x4 = 0.01),
    times_obs = seq(0, 10, by = 0.5),
    times_sim = seq(0, 10, by = 0.01),
    active_params = c("k1", "k2", "m", "n", "L0", "C0"),
    uncertain_bounds = list(
      k1 = c(35, 45), k2 = c(7800, 9600),
      m = c(0, 0.5), n = c(3, 8),
      L0 = c(0.001, 0.1), C0 = c(90, 110)
    ),
    fixed_params = list(
      T_ref = 273.15 + 23.5,
      T_vec = 273.15 + c(5, 20, 32, 37)
    ),
    param_scales = list(
      k1 = 1.0, k2 = 1000.0, m = 1.0,
      n = 1.0, L0 = 1.0, C0 = 1.0
    )
  )
)
