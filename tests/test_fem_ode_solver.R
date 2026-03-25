# =============================================================================
# tests/test_fem_ode_solver.R
#
# Tests for FemOdeSolver (src/solvers/fem_ode_solver.R).
#
# HOW TO RUN (from repo root):
#   source("tests/test_fem_ode_solver.R")
#
# WHAT IS TESTED:
#   F1. Gradient consistency (CN)  — K^T adjoint matches finite differences.
#   F2. Gradient consistency (GL2) — exact stage-level K_e^T matches FD.
#   F3. Gradient consistency (GL1) — midpoint K^T matches FD.
#   F4. Load-vector structure      — L equals the element-wise source term
#       used by OdeSystemSolver, confirming the FEM assembly is correct.
#   F5. Equivalence with OdeSystemSolver — FemOdeSolver and OdeSystemSolver
#       give identical gradients (same maths, different code structure).
#   F6. Global K matrix (CN)       — K^T λ = L is consistent with the
#       element-by-element backward substitution.
#   F7. Cost monotonicity          — optimize() reduces cost.
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/fem_ode_solver.R")
source("src/solvers/general_ode_system_solver.R")

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
describe("F1: FemOdeSolver CN — gradient consistency (1-D decay)", {
  params    <- list(k = 0.5)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 5.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  solver <- FemOdeSolver$new(decay_rhs, times_sim, times_sim, y_true,
                              params, lambda = 0.05, method = "cn")
  set.seed(1)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk <- check_gradient(solver$cost_function, solver$gradient_function,
                         u_test, eps = 1e-5, y0 = y0)

  test_that("CN FEM max relative FD error < 1e-3", {
    expect_less_than(chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e", chk$max_rel_error))
  })
  test_that("CN FEM cosine similarity > 0.9999", {
    expect_greater_than(chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity))
  })
})

# ---------------------------------------------------------------------------
# F2. Gradient consistency — GL2, 2-D Lotka-Volterra
# ---------------------------------------------------------------------------
describe("F2: FemOdeSolver GL2 — gradient consistency (2-D LV)", {
  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- c(5, 3)
  y_lv      <- euler_solve(lv_rhs, y0, times_sim, params)
  set.seed(2)
  obs_lv    <- y_lv + matrix(rnorm(length(y_lv), 0, 0.5), nrow(y_lv), 2)

  solver <- FemOdeSolver$new(lv_rhs, times_sim, times_sim, obs_lv,
                              params, lambda = 0.1, method = "gl2")
  set.seed(3)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)

  chk <- check_gradient(solver$cost_function, solver$gradient_function,
                         u_test, eps = 1e-5, y0 = y0)

  test_that("GL2 FEM max relative FD error < 1e-3", {
    expect_less_than(chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e", chk$max_rel_error))
  })
  test_that("GL2 FEM cosine similarity > 0.9999", {
    expect_greater_than(chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity))
  })
})

# ---------------------------------------------------------------------------
# F3. Gradient consistency — GL1, 1-D decay
# ---------------------------------------------------------------------------
describe("F3: FemOdeSolver GL1 — gradient consistency (1-D decay)", {
  params    <- list(k = 0.5)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 5.0
  y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  solver <- FemOdeSolver$new(decay_rhs, times_sim, times_sim, y_true,
                              params, lambda = 0.05, method = "gl1")
  set.seed(4)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.1)

  chk <- check_gradient(solver$cost_function, solver$gradient_function,
                         u_test, eps = 1e-5, y0 = y0)

  test_that("GL1 FEM max relative FD error < 1e-3", {
    expect_less_than(chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e", chk$max_rel_error))
  })
  test_that("GL1 FEM cosine similarity > 0.9999", {
    expect_greater_than(chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity))
  })
})

# ---------------------------------------------------------------------------
# F4. Load-vector structure
#     L[n,] must equal (2/ns)*(y_n - obs_n) for observed nodes, 0 for NA.
# ---------------------------------------------------------------------------
describe("F4: Load-vector assembly — exact point evaluation", {
  params    <- list(k = 0.4)
  times_sim <- seq(0, 2, by = 0.2)
  y0        <- 3.0
  ns        <- length(times_sim); nv <- 1
  obs       <- matrix(y0 * exp(-params$k * times_sim) + 0.1, ncol = 1)
  obs[c(3, 5), ] <- NA

  solver <- FemOdeSolver$new(decay_rhs, times_sim, times_sim, obs,
                              params, lambda = 0.1, method = "cn")
  u_zero <- rep(0, ns)
  # cost_function runs solve_forward internally, populating fwd_y for get_load_vector
  solver$cost_function(u_zero, y0)
  L <- solver$get_load_vector()

  test_that("L has correct dimensions (ns × nv)", {
    expect_equal(nrow(L), ns, tol = 0, "row count")
    expect_equal(ncol(L), nv, tol = 0, "col count")
  })

  test_that("L is zero at NA observation nodes", {
    expect_true(all(L[c(3, 5), ] == 0),
      sprintf("L at NA rows: %s", paste(L[c(3, 5), ], collapse = ", ")))
  })

  test_that("L equals (2/ns)*(y_fwd - obs) at observed nodes", {
    obs_clean <- solver$observations_mapped; obs_clean[is.na(obs_clean)] <- 0
    obs_idx   <- which(!is.na(solver$observations_mapped[, 1]))
    y_fwd     <- solver$solve_forward(matrix(u_zero, ns, nv), y0)$y
    L2        <- solver$get_load_vector()
    for (j in obs_idx) {
      expected <- (2 / ns) * (y_fwd[j, 1] - obs_clean[j, 1])
      expect_true(abs(L2[j, 1] - expected) < 1e-10,
        sprintf("node %d: L=%.6g expected=%.6g", j, L2[j, 1], expected))
    }
  })
})

# ---------------------------------------------------------------------------
# F5. Equivalence: FemOdeSolver and OdeSystemSolver give identical gradients
# ---------------------------------------------------------------------------
describe("F5: FemOdeSolver ≡ OdeSystemSolver — identical gradients", {
  for (meth in c("cn", "gl1", "gl2")) {
    local({
      m         <- meth
      params    <- list(k = 0.5)
      times_sim <- seq(0, 2, by = 0.1)
      y0        <- 5.0
      y_true    <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)
      lam       <- 0.05

      fem  <- FemOdeSolver$new(decay_rhs, times_sim, times_sim, y_true,
                                params, lam, method = m)
      ode  <- OdeSystemSolver$new(decay_rhs, times_sim, times_sim, y_true,
                                   params, lam, method = m)

      set.seed(99)
      u_test <- rnorm(fem$n_steps, sd = 0.05)

      g_fem <- fem$gradient_function(u_test, y0)
      g_ode <- ode$gradient_function(u_test, y0)

      cos_sim <- sum(g_fem * g_ode) /
                 (sqrt(sum(g_fem^2)) * sqrt(sum(g_ode^2)))
      max_rel <- max(abs(g_fem - g_ode) / (abs(g_ode) + 1e-12))

      test_that(sprintf("%s: FEM and ODE gradients agree (cos > 0.9999999)", m), {
        expect_greater_than(cos_sim, 0.9999999,
          sprintf("cosine = %.9f, max_rel = %.3e", cos_sim, max_rel))
      })
      test_that(sprintf("%s: FEM and ODE gradients max relative diff < 1e-10", m), {
        expect_less_than(max_rel, 1e-10,
          sprintf("max relative diff = %.3e", max_rel))
      })
    })
  }
})

# ---------------------------------------------------------------------------
# F6. Global K matrix (CN) — K^T λ = L recovers same adjoint
# ---------------------------------------------------------------------------
describe("F6: Global K (CN) — K^T backward solve agrees with element-wise", {
  params    <- list(k = 0.3)
  times_sim <- seq(0, 1, by = 0.1)
  y0        <- 2.0
  ns        <- length(times_sim); nv <- 1

  set.seed(5)
  obs <- matrix(y0 * exp(-params$k * times_sim) + rnorm(ns, 0, 0.05), ncol = 1)

  solver <- FemOdeSolver$new(decay_rhs, times_sim, times_sim, obs,
                              params, lambda = 0.1, method = "cn")
  set.seed(6)
  u_test <- rnorm(ns, sd = 0.05)
  solver$cost_function(u_test, y0)   # populate fwd_y

  # Run forward + adjoint element-wise (populates private$fwd_y)
  solver$cost_function(u_test, y0)
  adj      <- solver$solve_adjoint_fem()
  lam_elem <- adj$lambda          # ns × nv; lam_elem[1,] = adjoint at IC node

  # Global K^T: K has shape (ns-1)·nv × (ns-1)·nv (free nodes y[2:ns])
  # Solve K^T λ_free = L_free, then compare λ_free with lam_elem[2:ns, ]
  Kobj   <- solver$get_global_K()
  K      <- Kobj$K
  L      <- solver$get_load_vector()
  L_free <- as.vector(L[2:ns, , drop = FALSE])   # free nodes rows 2:ns

  lam_global <- matrix(solve(t(K), L_free), ns - 1L, nv)  # (ns-1) × nv

  max_diff <- max(abs(lam_global - lam_elem[2:ns, , drop = FALSE]))

  test_that("Global K^T solve matches element-wise backward substitution", {
    expect_less_than(max_diff, 1e-10,
      sprintf("max |lambda_global - lambda_elem[2:ns]| = %.3e", max_diff))
  })
})

# ---------------------------------------------------------------------------
# F7. Cost monotonicity — optimize() reduces cost
# ---------------------------------------------------------------------------
describe("F7: FemOdeSolver optimize() — cost strictly decreases", {
  params    <- list(k = 0.5)
  times_sim <- seq(0, 3, by = 0.1)
  obs_times <- seq(0, 3, by = 0.5)
  y0        <- 3.0
  set.seed(7)
  obs_data  <- matrix(y0 * exp(-params$k * obs_times) + rnorm(length(obs_times), 0, 0.1),
                      ncol = 1)

  solver    <- FemOdeSolver$new(decay_rhs, times_sim, obs_times, obs_data,
                                 params, lambda = 0.1, method = "gl2")
  u_zero    <- rep(0, solver$n_steps * solver$n_vars)
  cost_before <- solver$cost_function(u_zero, y0)

  set.seed(8)
  solver$optimize(y0, max_iter = 50)
  cost_after <- solver$cost_function(as.vector(solver$u), y0)

  test_that("cost after optimize() < cost at u=0", {
    expect_less_than(cost_after, cost_before,
      sprintf("before=%.6g  after=%.6g", cost_before, cost_after))
  })
  test_that("final cost is finite and non-negative", {
    expect_true(is.finite(cost_after) && cost_after >= 0)
  })
})

test_summary()
