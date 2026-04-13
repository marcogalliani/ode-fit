# =============================================================================
# Tests for FEMForwardSolver (src/solvers/fem_fwd_solver.R).
#
# WHAT IS TESTED:
#   F1. Gradient consistency (CN)
#   F2. Gradient consistency (GL2)
#   F3. Gradient consistency (GL1)
#   F4. Load-vector structure
#   F5. Equivalence with DtOForwardSolver
#   F6. Global K matrix (CN)
#   F7. Cost monotonicity
# =============================================================================

setwd("../../..")
source("src/solvers/forward-solvers/load_forward_solvers.R")

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
  out  <- matrix(0, length(times), length(y0)); out[1, ] <- y0
  for (i in seq_len(length(times) - 1))
    out[i + 1, ] <- out[i, ] + dt_v[i] * rhs(out[i, ], times[i], params)
  out
}

# ---------------------------------------------------------------------------
# F1. Gradient consistency — CN, 1-D decay
# ---------------------------------------------------------------------------
describe("F1: FEMForwardSolver CN — gradient consistency (1-D decay)", {
  params    <- list(k = 0.5)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 5.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  solver <- FEMForwardSolver$new(decay_rhs, times_sim, times_sim, y_true,
                                 params, lambda = 0.05, method = "cn")
  set.seed(1)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk <- check_gradient(solver$cost_function, solver$gradient_function,
                        u_test, eps = 1e-5, y0 = y0)

  test_that("CN FEM max relative FD error < 1e-3", {
    expect_lt(chk$max_rel_error, 1e-3)
  })
  test_that("CN FEM cosine similarity > 0.9999", {
    expect_gt(chk$cosine_similarity, 0.9999)
  })
})

# ---------------------------------------------------------------------------
# F2. Gradient consistency — GL2, 2-D Lotka-Volterra
# ---------------------------------------------------------------------------
describe("F2: FEMForwardSolver GL2 — gradient consistency (2-D LV)", {
  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- c(5, 3)
  y_lv      <- euler_solve(lv_rhs, y0, times_sim, params)
  set.seed(2)
  obs_lv    <- y_lv + matrix(rnorm(length(y_lv), 0, 0.5), nrow(y_lv), 2)

  solver <- FEMForwardSolver$new(lv_rhs, times_sim, times_sim, obs_lv,
                                 params, lambda = 0.1, method = "gl2")
  set.seed(3)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk <- check_gradient(solver$cost_function, solver$gradient_function,
                        u_test, eps = 1e-5, y0 = y0)

  test_that("GL2 FEM max relative FD error < 1e-3", {
    expect_lt(chk$max_rel_error, 1e-3)
  })
  test_that("GL2 FEM cosine similarity > 0.9999", {
    expect_gt(chk$cosine_similarity, 0.9999)
  })
})

# ---------------------------------------------------------------------------
# F3. Gradient consistency — GL1, 1-D decay
# ---------------------------------------------------------------------------
describe("F3: FEMForwardSolver GL1 — gradient consistency (1-D decay)", {
  params    <- list(k = 0.5)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 5.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  solver <- FEMForwardSolver$new(decay_rhs, times_sim, times_sim, y_true,
                                 params, lambda = 0.05, method = "gl1")
  set.seed(4)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.1)

  chk <- check_gradient(solver$cost_function, solver$gradient_function,
                        u_test, eps = 1e-5, y0 = y0)

  test_that("GL1 FEM max relative FD error < 1e-3", {
    expect_lt(chk$max_rel_error, 1e-3)
  })
  test_that("GL1 FEM cosine similarity > 0.9999", {
    expect_gt(chk$cosine_similarity, 0.9999)
  })
})

# ---------------------------------------------------------------------------
# F4. Load-vector structure
# ---------------------------------------------------------------------------
describe("F4: Load-vector assembly — exact point evaluation", {
  params    <- list(k = 0.4)
  times_sim <- seq(0, 2, by = 0.2)
  y0        <- 3.0
  ns        <- length(times_sim); nv <- 1
  obs       <- matrix(y0 * exp(-params$k * times_sim) + 0.1, ncol = 1)
  obs[c(3, 5), ] <- NA

  solver <- FEMForwardSolver$new(decay_rhs, times_sim, times_sim, obs,
                                 params, lambda = 0.1, method = "cn")
  u_zero <- rep(0, ns)
  solver$cost_function(u_zero, y0)
  L <- solver$get_load_vector()

  test_that("L has correct dimensions (ns × nv)", {
    expect_equal(nrow(L), ns, tolerance = 0)
    expect_equal(ncol(L), nv, tolerance = 0)
  })
  test_that("L is zero at NA observation nodes", {
    expect_true(all(L[c(3, 5), ] == 0))
  })
  test_that("L equals (2/ns)*(y_fwd - obs) at observed nodes", {
    obs_clean <- solver$observations_mapped; obs_clean[is.na(obs_clean)] <- 0
    obs_idx   <- which(!is.na(solver$observations_mapped[, 1]))
    y_fwd     <- solver$solve_forward(matrix(u_zero, ns, nv), y0)$y
    L2        <- solver$get_load_vector()
    for (j in obs_idx) {
      expected <- (2 / ns) * (y_fwd[j, 1] - obs_clean[j, 1])
      expect_lt(abs(L2[j, 1] - expected), 1e-10)
    }
  })
})

# ---------------------------------------------------------------------------
# F5. Equivalence: FEMForwardSolver and DtOForwardSolver give identical gradients
# ---------------------------------------------------------------------------
describe("F5: FEMForwardSolver == DtOForwardSolver — identical gradients", {
  for (meth in c("cn", "gl1", "gl2")) {
    local({
      m         <- meth
      params    <- list(k = 0.5)
      times_sim <- seq(0, 2, by = 0.1)
      y0        <- 5.0
      y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)
      lam       <- 0.05

      fem <- FEMForwardSolver$new(decay_rhs, times_sim, times_sim, y_true,
                                  params, lam, method = m)
      ode <- DtOForwardSolver$new(decay_rhs, times_sim, times_sim, y_true,
                                  params, lam, method = m)

      set.seed(99)
      u_test <- rnorm(fem$n_steps, sd = 0.05)

      g_fem   <- fem$gradient_function(u_test, y0)
      g_ode   <- ode$gradient_function(u_test, y0)
      cos_sim <- sum(g_fem * g_ode) / (sqrt(sum(g_fem^2)) * sqrt(sum(g_ode^2)))
      max_rel <- max(abs(g_fem - g_ode) / (abs(g_ode) + 1e-12))

      test_that(sprintf("%s: FEM and ODE gradients agree (cos > 0.9999999)", m), {
        expect_gt(cos_sim, 0.9999999)
      })
      test_that(sprintf("%s: FEM and ODE gradients max relative diff < 1e-10", m), {
        expect_lt(max_rel, 1e-10)
      })
    })
  }
})

# ---------------------------------------------------------------------------
# F6. Global K matrix (CN) — K^T lambda = L recovers same adjoint
# ---------------------------------------------------------------------------
describe("F6: Global K (CN) — K^T backward solve agrees with element-wise", {
  params    <- list(k = 0.3)
  times_sim <- seq(0, 1, by = 0.1)
  y0        <- 2.0
  ns        <- length(times_sim); nv <- 1

  set.seed(5)
  obs <- matrix(y0 * exp(-params$k * times_sim) + rnorm(ns, 0, 0.05), ncol = 1)

  solver <- FEMForwardSolver$new(decay_rhs, times_sim, times_sim, obs,
                                 params, lambda = 0.1, method = "cn")
  set.seed(6)
  u_test <- rnorm(ns, sd = 0.05)
  solver$cost_function(u_test, y0)

  adj      <- solver$solve_adjoint_fem()
  lam_elem <- adj$lambda

  Kobj      <- solver$get_global_K()
  K         <- Kobj$K
  L         <- solver$get_load_vector()
  L_free    <- as.vector(L[2:ns, , drop = FALSE])
  lam_global <- matrix(solve(t(K), L_free), ns - 1L, nv)

  max_diff <- max(abs(lam_global - lam_elem[2:ns, , drop = FALSE]))

  test_that("Global K^T solve matches element-wise backward substitution", {
    expect_lt(max_diff, 1e-10)
  })
})

# ---------------------------------------------------------------------------
# F7. Cost monotonicity — optimize() reduces cost
# ---------------------------------------------------------------------------
describe("F7: FEMForwardSolver optimize() — cost strictly decreases", {
  params    <- list(k = 0.5)
  times_sim <- seq(0, 3, by = 0.1)
  obs_times <- seq(0, 3, by = 0.5)
  y0        <- 3.0
  set.seed(7)
  obs_data  <- matrix(y0 * exp(-params$k * obs_times) + rnorm(length(obs_times), 0, 0.1),
                      ncol = 1)

  solver      <- FEMForwardSolver$new(decay_rhs, times_sim, obs_times, obs_data,
                                      params, lambda = 0.1, method = "gl2")
  u_zero      <- rep(0, solver$n_steps * solver$n_vars)
  cost_before <- solver$cost_function(u_zero, y0)

  set.seed(8)
  solver$optimize(y0, max_iter = 50)
  cost_after <- solver$cost_function(as.vector(solver$u), y0)

  test_that("cost after optimize() < cost at u=0", {
    expect_lt(cost_after, cost_before)
  })
  test_that("final cost is finite and non-negative", {
    expect_true(is.finite(cost_after) && cost_after >= 0)
  })
})
