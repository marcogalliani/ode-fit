# =============================================================================
# tests/test_general_ode_solver.R
#
# Test suite for OdeSystemSolver (src/solvers/general_ode_system_solver.R).
#
# HOW TO RUN (from repo root):
#   source("tests/test_general_ode_solver.R")
#
# WHAT IS TESTED:
#   T1. Gradient consistency  -- adjoint gradient vs central finite differences
#   T2. Zero-forcing baseline -- cost_function and gradient_function are finite
#       and self-consistent when u = 0.
#   T3. Cost monotonicity     -- optimize() strictly lowers the cost.
#   T4. State recovery        -- exact physics + dense noiseless obs => u* ≈ 0.
#   T5. Forcing recovery      -- known additive forcing is reconstructed by u*.
#   T6. NA handling           -- NAs in observations don't produce NaN/errors.
#   T7. Multi-variable gradient consistency (2D Lotka-Volterra system).
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/general_ode_system_solver.R")

# ---------------------------------------------------------------------------
# Shared physics definitions
# ---------------------------------------------------------------------------

# 1-D exponential decay:  dy/dt = -k * y,  y(t) = y0 * exp(-k*t)
decay_rhs <- function(y, t, p) -p$k * y

# 1-D system with a sinusoidal forcing the model does NOT know about:
#   truth:    dy/dt = -k*y + sin(t)
#   model:    dy/dt = -k*y          (forcing is "hidden")
decay_rhs_no_forcing <- function(y, t, p) -p$k * y

# 2-D Lotka-Volterra
lv_rhs <- function(y, t, p) c(
  p$alpha * y[1] - p$beta  * y[1] * y[2],
  p$delta * y[1] * y[2] - p$gamma * y[2]
)

# Simple forward-Euler helper (for ground-truth generation)
euler_solve <- function(rhs, y0, times, params) {
  dt_v <- c(diff(times), 0)
  out  <- matrix(0, length(times), length(y0))
  out[1, ] <- y0
  for (i in seq_len(length(times) - 1))
    out[i + 1, ] <- out[i, ] + dt_v[i] * rhs(out[i, ], times[i], params)
  out
}

# ---------------------------------------------------------------------------
# T1. Gradient consistency: 1-D exponential decay
# ---------------------------------------------------------------------------
describe("T1: Gradient consistency — 1-D exponential decay", {

  params    <- list(k = 0.5)
  times_sim <- seq(0, 2, by = 0.1)      # 21 steps
  y0        <- 5.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  solver <- OdeSystemSolver$new(
    func_rhs  = decay_rhs,
    times_sim = times_sim,
    obs_times = times_sim,
    obs_values = y_true,
    params    = params,
    lambda    = 0.05
  )

  set.seed(1)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk <- check_gradient(
    fn  = solver$cost_function,
    gr  = solver$gradient_function,
    par = u_test,
    eps = 1e-5,
    y0  = y0
  )

  test_that("max relative error < 1e-3 (adjoint matches finite differences)", {
    expect_less_than(
      chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e, cosine sim = %.6f",
              chk$max_rel_error, chk$cosine_similarity)
    )
  })

  test_that("cosine similarity > 0.9999 (gradient direction is correct)", {
    expect_greater_than(
      chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity)
    )
  })
})

# ---------------------------------------------------------------------------
# T2. Zero-forcing baseline: cost and gradient are finite and consistent
# ---------------------------------------------------------------------------
describe("T2: Zero-forcing baseline — finite values and self-consistency", {

  params    <- list(k = 0.3)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- 2.0
  y_obs     <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  solver <- OdeSystemSolver$new(
    func_rhs   = decay_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = y_obs,
    params     = params,
    lambda     = 0.1
  )

  u_zero <- rep(0, solver$n_steps * solver$n_vars)

  test_that("cost_function(u=0) returns a finite, non-negative scalar", {
    cost <- solver$cost_function(u_zero, y0)
    expect_no_na(cost)
    expect_true(is.finite(cost) && cost >= 0,
                sprintf("cost = %g", cost))
  })

  test_that("gradient_function(u=0) returns a finite vector", {
    g <- solver$gradient_function(u_zero, y0)
    expect_no_na(g)
    expect_true(all(is.finite(g)),
                sprintf("%d non-finite elements", sum(!is.finite(g))))
  })

  test_that("gradient length equals n_steps * n_vars", {
    g <- solver$gradient_function(u_zero, y0)
    expect_equal(length(g), solver$n_steps * solver$n_vars, tol = 0,
                 "gradient length mismatch")
  })
})

# ---------------------------------------------------------------------------
# T3. Cost monotonicity: optimize() must strictly reduce the cost
# ---------------------------------------------------------------------------
describe("T3: Cost monotonicity — optimize() reduces the cost", {

  params    <- list(k = 0.5)
  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.5)
  y0        <- 3.0
  y_true    <- matrix(y0 * exp(-params$k * obs_times), ncol = 1)

  set.seed(2)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.1), nrow(y_true), 1)

  solver <- OdeSystemSolver$new(
    func_rhs   = decay_rhs,
    times_sim  = times_sim,
    obs_times  = obs_times,
    obs_values = obs_data,
    params     = params,
    lambda     = 0.1
  )

  # Cost before optimization (u = 0)
  u_zero       <- rep(0, solver$n_steps * solver$n_vars)
  cost_before  <- solver$cost_function(u_zero, y0)

  set.seed(3)
  solver$optimize(y0 = y0, max_iter = 50)
  cost_after   <- solver$cost_function(as.vector(solver$u), y0)

  test_that("cost after optimize() is strictly lower than with u = 0", {
    expect_less_than(cost_after, cost_before,
                     sprintf("before=%.6g, after=%.6g", cost_before, cost_after))
  })

  test_that("final cost is non-negative and finite", {
    expect_no_na(cost_after)
    expect_true(is.finite(cost_after) && cost_after >= 0)
  })
})

# ---------------------------------------------------------------------------
# T4. State recovery: exact physics + dense noiseless obs => u* ≈ 0
#     When the model is exactly correct and obs cover every grid point,
#     the optimal forcing should be near-zero everywhere.
# ---------------------------------------------------------------------------
describe("T4: State recovery — u* ≈ 0 when physics and obs are perfect", {

  params    <- list(k = 0.4)
  times_sim <- seq(0, 3, by = 0.1)
  y0        <- 4.0

  # Ground truth generated with the SAME Euler integrator the solver uses
  y_true <- euler_solve(decay_rhs, y0, times_sim, params)

  solver <- OdeSystemSolver$new(
    func_rhs   = decay_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = y_true,   # exact, noiseless, dense
    params     = params,
    lambda     = 1.0       # high lambda strongly penalises non-zero u
  )

  set.seed(4)
  solver$optimize(y0 = y0, max_iter = 200)

  rmse_u    <- sqrt(mean(solver$u^2))
  rmse_y    <- sqrt(mean((solver$y - y_true)^2))

  test_that("RMSE of forcing u* is small (< 0.05)", {
    expect_less_than(rmse_u, 0.05,
                     sprintf("RMSE(u) = %.4g", rmse_u))
  })

  test_that("RMSE of fitted state vs ground truth is small (< 0.05)", {
    expect_less_than(rmse_y, 0.05,
                     sprintf("RMSE(y) = %.4g", rmse_y))
  })
})

# ---------------------------------------------------------------------------
# T5. Forcing recovery: known hidden forcing is reconstructed by u*
#     Truth:  dy/dt = -k*y + sin(t)  (sinusoidal forcing)
#     Model:  dy/dt = -k*y            (model ignores forcing)
#     The solver should reconstruct u*(t) ≈ sin(t).
# ---------------------------------------------------------------------------
describe("T5: Forcing recovery — u*(t) correlates with true sin(t) forcing", {

  k         <- 0.3
  params    <- list(k = k)
  times_sim <- seq(0, 6, by = 0.05)   # fine grid, step = 0.05
  y0        <- 2.0

  # Generate ground truth WITH forcing using the correct grid step
  rhs_true <- function(y, t, p) -p$k * y + sin(t)
  y_true   <- euler_solve(rhs_true, y0, times_sim, params)

  set.seed(5)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.05), nrow(y_true), 1)

  solver <- OdeSystemSolver$new(
    func_rhs   = decay_rhs_no_forcing,   # model does NOT know about sin(t)
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = obs_data,
    params     = params,
    lambda     = 1e-3    # small lambda: allow u to be non-zero
  )

  set.seed(6)
  solver$optimize(y0 = y0, max_iter = 200)

  true_force  <- sin(times_sim)
  recovered_u <- solver$u[, 1]

  # Pearson correlation between recovered u and true forcing
  cor_val <- cor(recovered_u, true_force)

  test_that("correlation between u*(t) and sin(t) is > 0.85", {
    expect_greater_than(cor_val, 0.85,
                        sprintf("cor(u*, sin(t)) = %.4f", cor_val))
  })
})

# ---------------------------------------------------------------------------
# T6. NA handling: NAs in obs must not propagate to cost or gradient
# ---------------------------------------------------------------------------
describe("T6: NA handling — observations with missing values", {

  params    <- list(k = 0.5)
  times_sim <- seq(0, 4, by = 0.2)
  obs_times <- times_sim
  y0        <- 3.0
  y_obs     <- matrix(y0 * exp(-params$k * obs_times), ncol = 1)

  # Introduce NAs in the middle stretch
  gap_idx <- which(obs_times >= 1.5 & obs_times <= 2.5)
  y_obs[gap_idx, ] <- NA

  solver <- OdeSystemSolver$new(
    func_rhs   = decay_rhs,
    times_sim  = times_sim,
    obs_times  = obs_times,
    obs_values = y_obs,
    params     = params,
    lambda     = 0.1
  )

  u_test <- rep(0.01, solver$n_steps * solver$n_vars)

  test_that("cost_function does not return NA with missing observations", {
    cost <- solver$cost_function(u_test, y0)
    expect_no_na(cost)
    expect_true(is.finite(cost))
  })

  test_that("gradient_function does not return NA with missing observations", {
    g <- solver$gradient_function(u_test, y0)
    expect_no_na(g)
    expect_true(all(is.finite(g)))
  })

  test_that("optimize() completes without error when NAs are present", {
    set.seed(7)
    res <- solver$optimize(y0 = y0, max_iter = 30)
    expect_no_na(res$value)
    expect_true(is.finite(res$value))
  })
})

# ---------------------------------------------------------------------------
# T7. Gradient consistency: 2-D Lotka-Volterra
#     Checks that the adjoint gradient is correct for a multi-variable system.
#     Designed to FAIL before fixes, PASS afterwards.
# ---------------------------------------------------------------------------
source("src/solvers/general_ode_system_solver.R")
describe("T7: Gradient consistency — 2-D Lotka-Volterra", {

  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)   # compact window for speed
  y0        <- c(5, 3)

  y_true <- euler_solve(lv_rhs, y0, times_sim, params)
  set.seed(8)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.5), nrow(y_true), 2)

  solver <- OdeSystemSolver$new(
    func_rhs   = lv_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = obs_data,
    params     = params,
    lambda     = 0.1,
    method = "cn"
  )

  set.seed(9)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk <- check_gradient(
    fn  = solver$cost_function,
    gr  = solver$gradient_function,
    par = u_test,
    eps = 1e-5,
    y0  = y0
  )

  test_that("max relative error < 1e-3 (2-D adjoint matches finite differences)", {
    expect_less_than(
      chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e, cosine sim = %.6f",
              chk$max_rel_error, chk$cosine_similarity)
    )
  })

  test_that("cosine similarity > 0.9999 (2-D gradient direction correct)", {
    expect_greater_than(
      chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity)
    )
  })
})

# ---------------------------------------------------------------------------
# T8. GL1 gradient consistency — 1-D exponential decay.
#     GL1 (implicit midpoint) is self-adjoint: the exact discrete adjoint
#     uses the same method run backward, so the gradient is exact.
#     Same tight thresholds as T1 (CN).
# ---------------------------------------------------------------------------
describe("T8: GL1 method — gradient consistency with finite differences (1-D decay)", {

  params    <- list(k = 0.5)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 5.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  solver <- OdeSystemSolver$new(
    func_rhs   = decay_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = y_true,
    params     = params,
    lambda     = 0.05,
    method     = "gl1"
  )

  set.seed(10)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.1)

  chk <- check_gradient(
    fn  = solver$cost_function,
    gr  = solver$gradient_function,
    par = u_test,
    eps = 1e-5,
    y0  = y0
  )

  test_that("GL1 max relative FD error < 1e-3 (exact discrete adjoint)", {
    expect_less_than(chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e, cosine sim = %.6f",
              chk$max_rel_error, chk$cosine_similarity))
  })

  test_that("GL1 cosine similarity > 0.9999", {
    expect_greater_than(chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity))
  })
})


# ---------------------------------------------------------------------------
# T9. GL1 gradient consistency — 2-D Lotka-Volterra.
#     Confirms self-adjoint exactness carries over to nonlinear multi-var.
# ---------------------------------------------------------------------------
describe("T9: GL1 method — gradient consistency with finite differences (2-D LV)", {

  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- c(5, 3)

  y_lv <- euler_solve(lv_rhs, y0, times_sim, params)
  set.seed(11)
  obs_lv <- y_lv + matrix(rnorm(length(y_lv), 0, 0.5), nrow(y_lv), 2)

  solver <- OdeSystemSolver$new(
    func_rhs   = lv_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = obs_lv,
    params     = params,
    lambda     = 0.1,
    method     = "gl1"
  )

  set.seed(12)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk_gl1 <- check_gradient(
    fn  = solver$cost_function,
    gr  = solver$gradient_function,
    par = u_test,
    eps = 1e-8,
    y0  = y0
  )

  test_that("GL1 2-D max relative FD error < 1e-3 (exact discrete adjoint)", {
    expect_less_than(chk_gl1$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e, cosine sim = %.6f",
              chk_gl1$max_rel_error, chk_gl1$cosine_similarity))
  })

  test_that("GL1 2-D cosine similarity > 0.9999", {
    expect_greater_than(chk_gl1$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk_gl1$cosine_similarity))
  })

  solver <- OdeSystemSolver$new(
    func_rhs   = lv_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = obs_lv,
    params     = params,
    lambda     = 0.1,
    method     = "gl2"
  )

  set.seed(12)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk_gl2 <- check_gradient(
    fn  = solver$cost_function,
    gr  = solver$gradient_function,
    par = u_test,
    eps = 1e-8,
    y0  = y0
  )

  test_that("GL2 max relative FD error < 1e-3 (exact discrete adjoint)", {
    expect_less_than(chk_gl2$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e, cosine sim = %.6f",
              chk_gl2$max_rel_error, chk_gl2$cosine_similarity))
  })

  test_that("GL2 cosine similarity > 0.9999", {
    expect_greater_than(chk_gl2$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk_gl2$cosine_similarity))
  })
})

# ---------------------------------------------------------------------------
# T11 (GL2 2-D). GL2 gradient consistency — 2-D Lotka-Volterra.
# ---------------------------------------------------------------------------
describe("T11: GL2 method — gradient consistency with finite differences (2-D LV)", {

  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- c(5, 3)

  y_lv <- euler_solve(lv_rhs, y0, times_sim, params)
  set.seed(21)
  obs_lv <- y_lv + matrix(rnorm(length(y_lv), 0, 0.5), nrow(y_lv), 2)

  solver <- OdeSystemSolver$new(
    func_rhs   = lv_rhs,
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = obs_lv,
    params     = params,
    lambda     = 0.1,
    method     = "gl2"
  )

  set.seed(22)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk <- check_gradient(
    fn  = solver$cost_function,
    gr  = solver$gradient_function,
    par = u_test,
    eps = 1e-5,
    y0  = y0
  )

  test_that("GL2 2-D max relative FD error < 1e-3 (exact discrete adjoint)", {
    expect_less_than(chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e, cosine sim = %.6f",
              chk$max_rel_error, chk$cosine_similarity))
  })

  test_that("GL2 2-D cosine similarity > 0.9999", {
    expect_greater_than(chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity))
  })
})


# ---------------------------------------------------------------------------
# T12. optimize() vs optimize_bvp() — 1-D exponential decay (method = "cn").
#      Both methods must reduce cost and produce close trajectories.
# ---------------------------------------------------------------------------
describe("T12: optimize vs optimize_bvp — 1-D decay, cost reduction and trajectory agreement", {

  params    <- list(k = 0.4)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 3.0

  set.seed(13)
  y_obs_10 <- matrix(
    y0 * exp(-params$k * times_sim) + rnorm(length(times_sim), 0, 0.05),
    ncol = 1
  )

  make_solver_10 <- function() {
    OdeSystemSolver$new(
      func_rhs   = decay_rhs,
      times_sim  = times_sim,
      obs_times  = times_sim,
      obs_values = y_obs_10,
      params     = params,
      lambda     = 0.1,
      method     = "cn"
    )
  }

  u_zero <- rep(0, length(times_sim))
  cost0  <- make_solver_10()$cost_function(u_zero, y0)

  solver_bfgs_10 <- make_solver_10()
  set.seed(14)
  solver_bfgs_10$optimize(y0 = y0, max_iter = 100)
  cost_bfgs <- solver_bfgs_10$cost_function(as.vector(solver_bfgs_10$u), y0)

  solver_bvp_10 <- make_solver_10()
  solver_bvp_10$optimize_bvp(y0 = y0)
  cost_bvp <- solver_bvp_10$cost_function(as.vector(solver_bvp_10$u), y0)

  test_that("optimize() strictly reduces cost vs u = 0", {
    expect_less_than(cost_bfgs, cost0,
      sprintf("before = %.6g, after = %.6g", cost0, cost_bfgs))
  })

  test_that("optimize_bvp() strictly reduces cost vs u = 0", {
    expect_less_than(cost_bvp, cost0,
      sprintf("before = %.6g, after = %.6g", cost0, cost_bvp))
  })

  test_that("optimize_bvp() and optimize() costs agree within a factor of 5", {
    expect_less_than(cost_bvp / cost_bfgs, 5,
      sprintf("bvp = %.6g, bfgs = %.6g", cost_bvp, cost_bfgs))
  })

  rmse_traj   <- sqrt(mean((solver_bfgs_10$y - solver_bvp_10$y)^2))
  state_range <- diff(range(y_obs_10))

  test_that("trajectories from optimize and optimize_bvp agree within 20% of state range", {
    expect_less_than(rmse_traj / state_range, 0.2,
      sprintf("RMSE = %.4g, range = %.4g", rmse_traj, state_range))
  })
})


# ---------------------------------------------------------------------------
# T11. optimize() vs optimize_bvp() — 2-D Lotka-Volterra.
#      Both methods must reduce cost and produce close trajectories.
# ---------------------------------------------------------------------------
describe("T13: optimize vs optimize_bvp — 2-D LV, cost reduction and trajectory agreement", {

  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- c(5, 3)

  y_lv_11 <- euler_solve(lv_rhs, y0, times_sim, params)
  set.seed(15)
  obs_lv_11 <- y_lv_11 + matrix(rnorm(length(y_lv_11), 0, 0.5), nrow(y_lv_11), 2)

  make_solver_11 <- function() {
    OdeSystemSolver$new(
      func_rhs   = lv_rhs,
      times_sim  = times_sim,
      obs_times  = times_sim,
      obs_values = obs_lv_11,
      params     = params,
      lambda     = 0.1,
      method     = "cn"
    )
  }

  u_zero <- rep(0, nrow(obs_lv_11) * ncol(obs_lv_11))
  cost0  <- make_solver_11()$cost_function(u_zero, y0)

  solver_bfgs_11 <- make_solver_11()
  set.seed(16)
  solver_bfgs_11$optimize(y0 = y0, max_iter = 100)
  cost_bfgs <- solver_bfgs_11$cost_function(as.vector(solver_bfgs_11$u), y0)

  solver_bvp_11 <- make_solver_11()
  solver_bvp_11$optimize_bvp(y0 = y0)
  cost_bvp <- solver_bvp_11$cost_function(as.vector(solver_bvp_11$u), y0)

  test_that("2-D LV: optimize() strictly reduces cost vs u = 0", {
    expect_less_than(cost_bfgs, cost0,
      sprintf("before = %.6g, after = %.6g", cost0, cost_bfgs))
  })

  test_that("2-D LV: optimize_bvp() strictly reduces cost vs u = 0", {
    expect_less_than(cost_bvp, cost0,
      sprintf("before = %.6g, after = %.6g", cost0, cost_bvp))
  })

  test_that("2-D LV: optimize_bvp() and optimize() costs agree within a factor of 5", {
    expect_less_than(cost_bvp / cost_bfgs, 5,
      sprintf("bvp = %.6g, bfgs = %.6g", cost_bvp, cost_bfgs))
  })

  rmse_traj   <- sqrt(mean((solver_bfgs_11$y - solver_bvp_11$y)^2))
  state_range <- diff(range(obs_lv_11))

  test_that("2-D LV: trajectories agree within 20% of state range", {
    expect_less_than(rmse_traj / state_range, 0.2,
      sprintf("RMSE = %.4g, range = %.4g", rmse_traj, state_range))
  })
})


# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------
test_summary()
