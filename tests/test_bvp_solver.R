# =============================================================================
# tests/test_bvp_solver.R
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/bvp_solver.R")
source("src/solvers/general_ode_system_solver.R")

# ─────────────────────────────────────────────────────────────────────────────
# Helper: verify BVP equations are satisfied by (Y, P)
# ─────────────────────────────────────────────────────────────────────────────
bvp_max_residuals <- function(sol, ns, ny, nrhs,
                               A_list, C_list, b_list,
                               D_list, E_list, f_list,
                               E_T, f_T) {
  Y <- sol$Y; P <- sol$P
  max_fwd <- max_bwd <- 0

  for (t in seq_len(ns - 1L)) {
    for (r in seq_len(nrhs)) {
      lhs_fwd <- Y[t + 1L, , r]
      rhs_fwd <- A_list[[t]] %*% Y[t, , r] -
                 C_list[[t]] %*% P[t + 1L, , r] + b_list[[t]][, r]
      max_fwd <- max(max_fwd, max(abs(lhs_fwd - rhs_fwd)))

      lhs_bwd <- P[t, , r]
      rhs_bwd <- D_list[[t]] %*% P[t + 1L, , r] +
                 E_list[[t]] %*% Y[t, , r] + f_list[[t]][, r]
      max_bwd <- max(max_bwd, max(abs(lhs_bwd - rhs_bwd)))
    }
  }

  # Terminal BC residual
  Y_T <- matrix(Y[ns, , ], ny, nrhs)
  P_T <- matrix(P[ns, , ], ny, nrhs)
  res_T <- max(abs(P_T - E_T %*% Y_T - f_T))

  list(forward = max_fwd, backward = max_bwd, terminal = res_T)
}


# ─────────────────────────────────────────────────────────────────────────────
cat("=== BV1: solve_linear_bvp_riccati — scalar, 2 RHS ===\n")
# ─────────────────────────────────────────────────────────────────────────────
describe("BV1", {
  set.seed(42)
  ns <- 20L; ny <- 1L; nrhs <- 2L
  dt <- 0.05; lambda <- 0.5
  a <- 0.3; c_t <- dt^2 / (2 * lambda)

  A_list <- replicate(ns - 1L, matrix(1 + dt * a, ny, ny), simplify = FALSE)
  C_list <- replicate(ns - 1L, matrix(c_t, ny, ny), simplify = FALSE)
  D_list <- A_list
  E_list <- replicate(ns - 1L, matrix(0.1, ny, ny), simplify = FALSE)
  E_T    <- matrix(0.1, ny, ny)

  b_list <- lapply(seq_len(ns - 1L), function(t) matrix(rnorm(ny * nrhs), ny, nrhs))
  f_list <- lapply(seq_len(ns - 1L), function(t) matrix(rnorm(ny * nrhs), ny, nrhs))
  f_T    <- matrix(rnorm(ny * nrhs), ny, nrhs)
  y0_mat <- matrix(c(1, -0.5), ny, nrhs)

  sol <- solve_linear_bvp_riccati(ns, ny, nrhs,
           A_list, C_list, b_list, D_list, E_list, f_list,
           y0_mat, E_T, f_T)
  res <- bvp_max_residuals(sol, ns, ny, nrhs,
           A_list, C_list, b_list, D_list, E_list, f_list, E_T, f_T)

  test_that("forward residual < 1e-12", expect_less_than(res$forward,  1e-12))
  test_that("backward residual < 1e-12", expect_less_than(res$backward, 1e-12))
  test_that("terminal BC residual < 1e-12", expect_less_than(res$terminal, 1e-12))
  test_that("initial condition satisfied",
    expect_less_than(max(abs(sol$Y[1L, , ] - y0_mat)), 1e-12))
})


# ─────────────────────────────────────────────────────────────────────────────
cat("=== BV2: solve_linear_bvp_riccati — ny=2, nrhs=1 ===\n")
# ─────────────────────────────────────────────────────────────────────────────
describe("BV2", {
  set.seed(7)
  ns <- 15L; ny <- 2L; nrhs <- 1L
  dt <- 0.1; lambda <- 1.0; c_t <- dt^2 / (2 * lambda)

  theta <- 0.1
  A0 <- diag(ny) + dt * matrix(c(-0.5, -theta, theta, -0.5), ny, ny)
  A_list <- replicate(ns - 1L, A0, simplify = FALSE)
  C_list <- replicate(ns - 1L, c_t * diag(ny), simplify = FALSE)
  D_list <- replicate(ns - 1L, t(A0), simplify = FALSE)
  E_obs  <- (2 / ns) * diag(ny)
  E_list <- replicate(ns - 1L, E_obs, simplify = FALSE)
  E_T    <- E_obs

  b_list <- lapply(seq_len(ns - 1L), function(t) matrix(sin(t) * c(0.1, -0.1), ny, 1))
  f_list <- lapply(seq_len(ns - 1L), function(t) matrix(cos(t) * c(0.05, 0.05), ny, 1))
  f_T    <- matrix(0, ny, 1)
  y0_mat <- matrix(0, ny, 1)

  sol <- solve_linear_bvp_riccati(ns, ny, nrhs,
           A_list, C_list, b_list, D_list, E_list, f_list,
           y0_mat, E_T, f_T)
  res <- bvp_max_residuals(sol, ns, ny, nrhs,
           A_list, C_list, b_list, D_list, E_list, f_list, E_T, f_T)

  test_that("forward residual < 1e-12", expect_less_than(res$forward,  1e-12))
  test_that("backward residual < 1e-12", expect_less_than(res$backward, 1e-12))
  test_that("terminal BC residual < 1e-12", expect_less_than(res$terminal, 1e-12))
})


# ─────────────────────────────────────────────────────────────────────────────
cat("=== BV3: gradient_function consistent with finite differences of cost_function ===\n")
# ─────────────────────────────────────────────────────────────────────────────
# gradient_function must be the exact gradient of cost_function (CN
# discretisation + exact discrete adjoint).  Verified via central FD.
describe("BV3", {
  k <- 0.5
  f_lin <- function(y, t, p) -p$k * y

  set.seed(1)
  times  <- seq(0, 1, by = 0.1)
  obs_v  <- matrix(exp(-k * times) + rnorm(length(times), sd = 0.05), ncol = 1)
  lambda <- 0.3

  solver <- OdeSystemSolver$new(
    func_rhs = f_lin, times_sim = times,
    obs_times = times, obs_values = obs_v,
    params = list(k = k), lambda = lambda
  )

  set.seed(2)
  u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.1)

  chk <- check_gradient(
    fn  = solver$cost_function,
    gr  = solver$gradient_function,
    par = u_test,
    eps = 1e-5,
    y0  = 1.0
  )

  test_that("max relative FD error < 1e-3 (CN adjoint matches finite differences)",
    expect_less_than(chk$max_rel_error, 1e-3,
      sprintf("max rel error = %.4e, cosine sim = %.6f",
              chk$max_rel_error, chk$cosine_similarity)))

  test_that("cosine similarity > 0.9999",
    expect_greater_than(chk$cosine_similarity, 0.9999,
      sprintf("cosine sim = %.6f", chk$cosine_similarity)))
})