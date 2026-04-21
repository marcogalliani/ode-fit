# =============================================================================
# Tests for solve_bvp_colloc and DtOForwardSolver$optimize_bvp.
# =============================================================================

setwd("../../..")
source("src/solvers/forward-solvers/load_forward_solvers.R")


# =============================================================================
# BVC1: solve_bvp_colloc — linear BVP with analytic solution
#
# y'' - y = 0  rewritten as first-order system:
#   z' = [z2; z1]   (z1 = y, z2 = y')
# BCs: z1(0) = 0,  z1(1) = sinh(1)
# Analytic solution: z1(t) = sinh(t), z2(t) = cosh(t)
# =============================================================================
describe("BVC1: solve_bvp_colloc — linear BVP, known solution", {
  ns    <- 201L
  t_g   <- seq(0, 1, length.out = ns)

  F_lin  <- function(t, z) c(z[2], z[1])
  bc_lin <- function(z_l, z_r) c(z_l[1] - 0, z_r[1] - sinh(1))

  z_init       <- matrix(0, ns, 2L)
  z_init[, 1]  <- seq(0, sinh(1), length.out = ns)
  z_init[, 2]  <- 1.0

  sol <- solve_bvp_colloc(F_lin, bc_lin, t_g, z_init, tol = 1e-10)

  test_that("BVC1 Newton converged", {
    expect_true(sol$converged)
  })
  test_that("BVC1 residual < tol", {
    expect_lt(sol$residual_norm, 1e-10)
  })

  err_y  <- max(abs(sol$z[, 1] - sinh(t_g)))
  err_yp <- max(abs(sol$z[, 2] - cosh(t_g)))
  test_that("BVC1 solution z1 accurate (< 5e-5)", {
    expect_lt(err_y, 5e-5)
  })
  test_that("BVC1 solution z2 accurate (< 5e-5)", {
    expect_lt(err_yp, 5e-5)
  })
})


# =============================================================================
# BVC2: solve_bvp_colloc — nonlinear BVP (Bratu problem)
#
# y'' + exp(y) = 0  rewritten as:
#   z' = [z2; -exp(z1)]
# BCs: z1(0) = 0,  z1(1) = 0
# =============================================================================
describe("BVC2: solve_bvp_colloc — nonlinear Bratu BVP", {
  ns  <- 101L
  t_g <- seq(0, 1, length.out = ns)

  F_bratu  <- function(t, z) c(z[2], -exp(z[1]))
  bc_bratu <- function(z_l, z_r) c(z_l[1] - 0, z_r[1] - 0)

  z_init      <- matrix(0, ns, 2L)
  z_init[, 1] <- 0.1 * t_g * (1 - t_g)
  z_init[, 2] <- 0.1 * (1 - 2 * t_g)

  sol <- solve_bvp_colloc(F_bratu, bc_bratu, t_g, z_init, tol = 1e-10)

  test_that("BVC2 Newton converged", {
    expect_true(sol$converged)
  })
  test_that("BVC2 residual < tol", {
    expect_lt(sol$residual_norm, 1e-10)
  })
  test_that("BVC2 left BC satisfied", {
    expect_lt(abs(sol$z[1L, 1L]), 1e-10)
  })
  test_that("BVC2 right BC satisfied", {
    expect_lt(abs(sol$z[ns, 1L]), 1e-10)
  })

  t_mid  <- round(ns / 2)
  dt_mid <- t_g[t_mid + 1L] - t_g[t_mid]
  dz1_fd <- (sol$z[t_mid + 1L, 1L] - sol$z[t_mid - 1L, 1L]) / (2 * dt_mid)
  test_that("BVC2 ODE consistency at midpoint (|dz1/dt - z2| < 1e-3)", {
    expect_lt(abs(dz1_fd - sol$z[t_mid, 2L]), 1e-3)
  })
})


# =============================================================================
# BVC3: optimize_bvp — KKT conditions satisfied at solution
# =============================================================================
describe("BVC3: optimize_bvp — gradient norm = 0 at solution", {
  set.seed(42)
  k      <- 0.5
  lambda <- 0.3
  f_lin  <- function(y, t, p) -p$k * y

  times  <- seq(0, 2, by = 0.05)
  y_true <- exp(-k * times)
  obs_v  <- matrix(y_true + rnorm(length(times), sd = 0.05), ncol = 1L)

  solver <- DtOForwardSolver$new(
    model      = ODEModel$new(rhs = f_lin),
    times_sim  = times,
    obs_times  = times,
    obs_values = obs_v,
    params     = list(k = k),
    lambda     = lambda
  )

  res <- solver$optimize_bvp(y0 = 1.0, tol = 1e-8)

  test_that("BVC3 BVP converged", {
    expect_true(res$converged)
  })

  u_opt  <- as.vector(solver$u)
  g_opt  <- solver$gradient_function(u_opt, y0 = 1.0)
  g_norm <- sqrt(sum(g_opt^2))

  test_that("BVC3 gradient norm at BVP solution < 5e-3", {
    expect_lt(g_norm, 5e-3)
  })
})


# =============================================================================
# BVC4: optimize_bvp vs optimize (BFGS) — trajectories agree
# =============================================================================
describe("BVC4: optimize_bvp vs optimize (BFGS) — trajectory agreement", {
  set.seed(7)
  k      <- 0.8
  lambda <- 0.5
  f_lin  <- function(y, t, p) -p$k * y

  times  <- seq(0, 1, by = 0.02)
  y_true <- exp(-k * times)
  obs_v  <- matrix(y_true + rnorm(length(times), sd = 0.03), ncol = 1L)

  make_solver <- function() {
    DtOForwardSolver$new(
      model      = ODEModel$new(rhs = f_lin),
      times_sim  = times,
      obs_times  = times,
      obs_values = obs_v,
      params     = list(k = k),
      lambda     = lambda
    )
  }

  s_bfgs <- make_solver()
  s_bfgs$optimize(y0 = 1.0, max_iter = 200,
                  reltol = sqrt(.Machine$double.eps))

  s_bvp <- make_solver()
  s_bvp$optimize_bvp(y0 = 1.0, tol = 1e-8)

  state_range <- max(s_bfgs$y) - min(s_bfgs$y)
  max_diff    <- max(abs(s_bvp$y - s_bfgs$y))
  rel_diff    <- max_diff / state_range

  test_that("BVC4 trajectories agree within 5% of state range", {
    expect_lt(rel_diff, 0.05)
  })

  sse_bfgs <- sum((s_bfgs$y - s_bfgs$observations_mapped)^2, na.rm = TRUE)
  sse_bvp  <- sum((s_bvp$y  - s_bvp$observations_mapped)^2,  na.rm = TRUE)
  test_that("BVC4 BVP SSE within 2x of BFGS SSE", {
    expect_lt(sse_bvp / sse_bfgs, 2.0)
  })
})
