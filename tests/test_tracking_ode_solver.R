# =============================================================================
# tests/test_tracking_ode_solver.R
#
# Test suite for TrackingOdeSolver (src/solvers/tracking_ode_solver.R).
#
# HOW TO RUN (from repo root):
#   source("tests/test_tracking_ode_solver.R")
#
# WHAT IS TESTED:
#   TR1. Outer gradient consistency — both gradient methods tested against FD
#        on the full Lotka-Volterra system (4 free parameters).
#        Also checks cross-consistency between the adjoint method
#        (outer_gradient) and the forward-sensitivity method
#        (outer_gradient_sensitivity).
#   TR2. Descent direction — at a non-optimal theta, the negative gradient is
#        a descent direction for the tracking outer objective.
#   TR3. Caching correctness — repeated outer_objective call hits cache.
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/tracking_ode_solver.R")  # also sources general_ode_system_solver.R

# ---------------------------------------------------------------------------
# Shared physics
# ---------------------------------------------------------------------------

decay_rhs <- function(y, t, p) -p$k * y

lv_rhs <- function(y, t, p) c(
  p$alpha * y[1] - p$beta  * y[1] * y[2],
  p$delta * y[1] * y[2] - p$gamma * y[2]
)

euler_solve <- function(rhs, y0, times, params) {
  dt_v <- c(diff(times), 0)
  out  <- matrix(0, length(times), length(y0))
  out[1, ] <- y0
  for (i in seq_len(length(times) - 1))
    out[i + 1, ] <- out[i, ] + dt_v[i] * rhs(out[i, ], times[i], params)
  out
}

# ---------------------------------------------------------------------------
# TR1. Outer gradient consistency — Lotka-Volterra, all 4 parameters.
#
#   Both gradient methods are tested against FD and against each other:
#     • outer_gradient             — adjoint formula: sum_t p[t+1]^T*f_theta*dt
#     • outer_gradient_sensitivity — forward sensitivity sweep
#
#   Because the two methods are mathematically equivalent (adjoint identity),
#   their cross-consistency threshold is set to 0.9999 (near machine precision),
#   versus 0.97 for the cascading solver's two methods (which implement
#   different Gauss-Newton approximations and can diverge).
#
#   Checks (per method vs FD):
#     • Cosine similarity > 0.98 — direction test
#     • Sign agreement per component — tightest magnitude check
# ---------------------------------------------------------------------------
describe("TR1: Outer gradient consistency — Lotka-Volterra (all parameters)", {

  p_true    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  y0_true   <- c(10, 10)
  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.5)

  y_true <- euler_solve(lv_rhs, y0_true, obs_times, p_true)
  set.seed(10)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.3),
                              nrow(y_true), ncol(y_true))

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)
  param_names  <- c("alpha", "beta", "delta", "gamma")

  tracking <- TrackingOdeSolver$new(
    func_rhs     = lv_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    y0           = y0_true,
    fixed_params = list(),
    lambda       = 1.0,
    param_scales = param_scales
  )

  # Test point: each parameter perturbed so all 4 gradient components differ
  theta_test <- c(alpha = 0.8, beta = 0.4, delta = 0.4, gamma = 0.6)

  set.seed(11)
  tracking$outer_objective(theta_test, param_names)
  g_adj  <- tracking$outer_gradient(theta_test, param_names)
  g_sens <- tracking$outer_gradient_sensitivity(theta_test, param_names)

  # Numerical FD reference (central differences)
  eps  <- 1e-3
  g_fd <- numeric(length(theta_test))
  for (j in seq_along(theta_test)) {
    th_p <- theta_test; th_p[j] <- theta_test[j] + eps
    th_m <- theta_test; th_m[j] <- theta_test[j] - eps
    set.seed(11); f_p <- tracking$outer_objective(th_p, param_names)
    set.seed(11); f_m <- tracking$outer_objective(th_m, param_names)
    g_fd[j] <- (f_p - f_m) / (2 * eps)
  }

  cos_sim <- function(a, b)
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-16)
  fmt_vec <- function(x) paste(sprintf("%+.3e", x), collapse = ", ")

  cos_adj   <- cos_sim(g_adj,  g_fd)
  cos_sens  <- cos_sim(g_sens, g_fd)
  cos_cross <- cos_sim(g_adj,  g_sens)

  # --- Adjoint method vs FD ---
  test_that("adjoint gradient direction (cosine similarity) > 0.98", {
    expect_greater_than(
      cos_adj, 0.98,
      sprintf("[adj]  cosine = %.4f\n  adj=(%s)\n   fd=(%s)",
              cos_adj, fmt_vec(g_adj), fmt_vec(g_fd))
    )
  })

  test_that("adjoint gradient component signs agree with FD", {
    expect_true(
      all(sign(g_adj) == sign(g_fd)),
      sprintf("[adj]  sign mismatch\n  adj=(%s)\n   fd=(%s)",
              fmt_vec(g_adj), fmt_vec(g_fd))
    )
  })

  # --- Sensitivity method vs FD ---
  test_that("sensitivity gradient direction (cosine similarity) > 0.98", {
    expect_greater_than(
      cos_sens, 0.98,
      sprintf("[sens] cosine = %.4f\n  sens=(%s)\n    fd=(%s)",
              cos_sens, fmt_vec(g_sens), fmt_vec(g_fd))
    )
  })

  test_that("sensitivity gradient component signs agree with FD", {
    expect_true(
      all(sign(g_sens) == sign(g_fd)),
      sprintf("[sens] sign mismatch\n  sens=(%s)\n    fd=(%s)",
              fmt_vec(g_sens), fmt_vec(g_fd))
    )
  })

  # --- Cross-consistency: adjoint vs forward sensitivity ---
  # The two methods are MATHEMATICALLY IDENTICAL (adjoint identity theorem).
  # Any deviation is purely numerical (different accumulation order), so the
  # threshold is much tighter than for cascading's two methods (0.97 there).
  test_that("adjoint and sensitivity gradients are mutually consistent (cos > 0.9999)", {
    expect_greater_than(
      cos_cross, 0.9999,
      sprintf("cross cosine = %.6f\n  adj =(%s)\n  sens=(%s)",
              cos_cross, fmt_vec(g_adj), fmt_vec(g_sens))
    )
  })
})

# ---------------------------------------------------------------------------
# TR2. Descent direction: at a non-optimal theta, the negative tracking
#      gradient is a descent direction for the tracking outer objective.
# ---------------------------------------------------------------------------
describe("TR2: Descent direction — outer gradient points downhill", {

  true_k    <- 0.4
  times_sim <- seq(0, 4, by = 0.2)
  obs_times <- seq(0, 4, by = 0.5)
  y0_true   <- 3.0

  y_true   <- euler_solve(decay_rhs, y0_true, obs_times, list(k = true_k))
  set.seed(12)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.05), nrow(y_true), 1)

  param_scales <- list(k = 0.5)

  tracking <- TrackingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    y0           = y0_true,
    fixed_params = list(),
    lambda       = 0.3,
    param_scales = param_scales
  )

  theta_test <- c(k = true_k / param_scales$k * 1.4)  # 40% above truth

  set.seed(13)
  f0 <- tracking$outer_objective(theta_test, "k")
  g  <- tracking$outer_gradient(theta_test, "k")

  step_size  <- 0.05
  theta_step <- theta_test - step_size * g
  theta_step <- pmax(theta_step, 1e-3)

  set.seed(13)
  f1 <- tracking$outer_objective(theta_step, "k")

  test_that("taking a step along negative gradient reduces the tracking outer objective", {
    expect_less_than(f1, f0,
                     sprintf("f0 = %.6g, f1 = %.6g", f0, f1))
  })
})

# ---------------------------------------------------------------------------
# TR3. Caching correctness — calling outer_objective twice with the same
#      theta must return the same value without re-running the inner solver.
# ---------------------------------------------------------------------------
describe("TR5: Caching correctness", {

  times_sim <- seq(0, 3, by = 0.2)
  obs_times <- seq(0, 3, by = 0.5)
  y0_true   <- 2.0
  params    <- list(k = 0.5)
  y_true    <- euler_solve(decay_rhs, y0_true, obs_times, params)
  set.seed(18)
  obs_data  <- y_true + matrix(rnorm(length(y_true), 0, 0.05), nrow(y_true), 1)

  tracking <- TrackingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    y0           = y0_true,
    fixed_params = list(),
    lambda       = 0.3,
    param_scales = list(k = 0.5)
  )

  theta <- c(k = 1.0)

  set.seed(19)
  f1 <- tracking$outer_objective(theta, "k")
  f2 <- tracking$outer_objective(theta, "k")   # cache hit — no new inner solve

  test_that("outer_objective returns identical value on repeated call (cache hit)", {
    expect_equal(f1, f2, tol = 0, "values differ: cache may not be working")
  })

  test_that("last_theta is set after the first call", {
    expect_true(!is.null(tracking$last_theta))
  })

  test_that("last_solver is populated after the first call", {
    expect_true(!is.null(tracking$last_solver))
  })
})

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
test_summary()
