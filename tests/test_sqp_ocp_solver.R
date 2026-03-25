# =============================================================================
# tests/test_sqp_ocp_solver.R
#
# Tests for SqpOcpSolver (src/solvers/sqp_ocp_solver.R).
#
# HOW TO RUN (from repo root):
#   source("tests/test_sqp_ocp_solver.R")
#
# WHAT IS TESTED:
#   SQP1. 1-D decay — convergence and cost reduction
#   SQP2. 1-D decay — trajectory agreement with BFGS (Euler)
#   SQP3. 2-D Lotka-Volterra — convergence and BFGS agreement
#   SQP4. State recovery — u* ≈ 0 when physics is exact
#   SQP5. NA handling — missing observations
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/sqp_ocp_solver.R")
source("src/solvers/general_ode_system_solver.R")

# ---------------------------------------------------------------------------
# Shared physics
# ---------------------------------------------------------------------------
.sqp_decay_rhs <- function(y, t, p) -p$k * y

.sqp_lv_rhs <- function(y, t, p) c(
  p$alpha * y[1] - p$beta  * y[1] * y[2],
  p$delta * y[1] * y[2] - p$gamma * y[2]
)

.sqp_euler <- function(rhs, y0, times, params) {
  dt_v <- c(diff(times), 0)
  out  <- matrix(0, length(times), length(y0))
  out[1, ] <- y0
  for (i in seq_len(length(times) - 1))
    out[i + 1, ] <- out[i, ] + dt_v[i] * rhs(out[i, ], times[i], params)
  out
}

# ---------------------------------------------------------------------------
# SQP1. 1-D decay — convergence and cost reduction
# ---------------------------------------------------------------------------
describe("SQP1: 1-D decay — convergence and cost reduction", {

  params    <- list(k = 0.4)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 3.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  set.seed(30)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.1), nrow(y_true), 1)

  sqp <- SqpOcpSolver$new(
    func_rhs = .sqp_decay_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = obs_data,
    params = params, lambda = 0.1, y0 = y0
  )

  ref <- OdeSystemSolver$new(
    func_rhs = .sqp_decay_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = obs_data,
    params = params, lambda = 0.1, method = "euler"
  )
  cost_zero <- ref$cost_function(rep(0, ref$n_steps * ref$n_vars), y0)

  sol      <- sqp$solve(max_iter = 50, tol = 1e-8)
  cost_sqp <- sqp$cost_function()

  test_that("SQP converges (KKT norm < 1e-6)", {
    expect_less_than(sol$kkt_norm, 1e-6,
      sprintf("kkt_norm = %.3e, iter = %d", sol$kkt_norm, sol$iter))
  })

  test_that("SQP cost < u=0 cost", {
    expect_less_than(cost_sqp, cost_zero,
      sprintf("sqp=%.6g, zero=%.6g", cost_sqp, cost_zero))
  })
})

# ---------------------------------------------------------------------------
# SQP2. 1-D decay — SQP vs BFGS trajectory agreement
# ---------------------------------------------------------------------------
describe("SQP2: 1-D decay — SQP vs BFGS trajectory agreement", {

  params    <- list(k = 0.4)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 3.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  set.seed(31)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.1), nrow(y_true), 1)

  sqp <- SqpOcpSolver$new(
    func_rhs = .sqp_decay_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = obs_data,
    params = params, lambda = 0.1, y0 = y0
  )
  sqp$solve(max_iter = 50, tol = 1e-8)

  ref <- OdeSystemSolver$new(
    func_rhs = .sqp_decay_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = obs_data,
    params = params, lambda = 0.1, method = "euler"
  )
  ref$optimize(y0 = y0, max_iter = 200)

  rmse_y <- sqrt(mean((sqp$y - ref$y)^2))

  test_that("SQP and BFGS state trajectories agree (RMSE < 0.05)", {
    expect_less_than(rmse_y, 0.05,
      sprintf("RMSE(y) = %.4g", rmse_y))
  })

  cost_sqp  <- sqp$cost_function()
  cost_bfgs <- ref$cost_function(as.vector(ref$u), y0)

  test_that("SQP and BFGS costs agree within factor 2", {
    ratio <- max(cost_sqp, cost_bfgs) / min(cost_sqp, cost_bfgs)
    expect_less_than(ratio, 2,
      sprintf("sqp=%.6g, bfgs=%.6g, ratio=%.3f", cost_sqp, cost_bfgs, ratio))
  })
})

# ---------------------------------------------------------------------------
# SQP3. 2-D Lotka-Volterra — convergence and BFGS agreement
# ---------------------------------------------------------------------------
describe("SQP3: 2-D Lotka-Volterra — convergence and BFGS agreement", {

  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- c(5, 3)

  y_lv <- .sqp_euler(.sqp_lv_rhs, y0, times_sim, params)
  set.seed(32)
  obs_data <- y_lv + matrix(rnorm(length(y_lv), 0, 0.5), nrow(y_lv), 2)

  sqp <- SqpOcpSolver$new(
    func_rhs = .sqp_lv_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = obs_data,
    params = params, lambda = 0.1, y0 = y0
  )
  sol <- sqp$solve(max_iter = 50, tol = 1e-6)

  ref <- OdeSystemSolver$new(
    func_rhs = .sqp_lv_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = obs_data,
    params = params, lambda = 0.1, method = "euler"
  )
  cost_zero <- ref$cost_function(rep(0, ref$n_steps * ref$n_vars), y0)
  ref$optimize(y0 = y0, max_iter = 200)

  cost_sqp  <- sqp$cost_function()
  cost_bfgs <- ref$cost_function(as.vector(ref$u), y0)

  test_that("2-D LV: SQP converges (KKT norm < 1e-4)", {
    expect_less_than(sol$kkt_norm, 1e-4,
      sprintf("kkt_norm = %.3e", sol$kkt_norm))
  })

  test_that("2-D LV: SQP cost < u=0 cost", {
    expect_less_than(cost_sqp, cost_zero,
      sprintf("sqp=%.6g, zero=%.6g", cost_sqp, cost_zero))
  })

  rmse_y      <- sqrt(mean((sqp$y - ref$y)^2))
  state_range <- diff(range(obs_data))

  test_that("2-D LV: trajectories agree within 20% of state range", {
    expect_less_than(rmse_y / state_range, 0.2,
      sprintf("RMSE/range = %.4g", rmse_y / state_range))
  })
})

# ---------------------------------------------------------------------------
# SQP4. State recovery — u* ≈ 0 when physics is exact
# ---------------------------------------------------------------------------
describe("SQP4: state recovery — u* ≈ 0 with exact physics", {

  params    <- list(k = 0.4)
  times_sim <- seq(0, 3, by = 0.1)
  y0        <- 4.0
  y_true    <- .sqp_euler(.sqp_decay_rhs, y0, times_sim, params)

  sqp <- SqpOcpSolver$new(
    func_rhs = .sqp_decay_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = y_true,
    params = params, lambda = 1.0, y0 = y0
  )
  sqp$solve(max_iter = 50, tol = 1e-10)

  rmse_u <- sqrt(mean(sqp$u^2))
  rmse_y <- sqrt(mean((sqp$y - y_true)^2))

  test_that("RMSE of u* is small (< 0.05)", {
    expect_less_than(rmse_u, 0.05,
      sprintf("RMSE(u) = %.4g", rmse_u))
  })

  test_that("RMSE of y vs truth is small (< 0.05)", {
    expect_less_than(rmse_y, 0.05,
      sprintf("RMSE(y) = %.4g", rmse_y))
  })
})

# ---------------------------------------------------------------------------
# SQP5. NA handling — missing observations
# ---------------------------------------------------------------------------
describe("SQP5: NA handling — missing observations", {

  params    <- list(k = 0.5)
  times_sim <- seq(0, 4, by = 0.2)
  y0        <- 3.0
  y_obs     <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  gap_idx <- which(times_sim >= 1.5 & times_sim <= 2.5)
  y_obs[gap_idx, ] <- NA

  sqp <- SqpOcpSolver$new(
    func_rhs = .sqp_decay_rhs, times_sim = times_sim,
    obs_times = times_sim, obs_values = y_obs,
    params = params, lambda = 0.1, y0 = y0
  )
  sol <- sqp$solve(max_iter = 50, tol = 1e-8)

  test_that("SQP converges with NA observations", {
    expect_true(sol$converged, sprintf("kkt_norm = %.3e", sol$kkt_norm))
  })

  test_that("cost is finite and non-negative", {
    cost <- sqp$cost_function()
    expect_true(is.finite(cost) && cost >= 0, sprintf("cost = %g", cost))
  })
})

# ---------------------------------------------------------------------------
test_summary()
