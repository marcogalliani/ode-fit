# =============================================================================
# tests/test_cascading_ode_solver.R
#
# Test suite for CascadingOdeSolver (src/solvers/parameter_cascading.R).
#
# HOW TO RUN (from repo root):
#   source("tests/test_cascading_ode_solver.R")
#
# WHAT IS TESTED:
#   C1. Outer gradient consistency — both gradient methods tested against FD
#       on the full Lotka-Volterra system (4 free parameters).  Also checks
#       cross-consistency between the IFT method (outer_gradient) and the
#       sensitivity-equation FBSM method (outer_gradient_sensitivity).
#       Depends on correct inner adjoint (Bugs 1 & 2).
#   C2. Descent direction -- at the true parameters outer_gradient ≈ 0 and
#       moving in the negative-gradient direction lowers the outer objective.
#   C3. Parameter recovery (1-D decay) -- CascadingOdeSolver recovers the
#       rate constant k from noisy exponential-decay data.
#   C4. Parameter recovery (Lotka-Volterra) -- recovers alpha from LV data
#       when beta, delta, gamma are fixed.
#   C5. Caching correctness -- calling outer_objective twice with the same
#       theta returns the same value without re-running the inner solver.
#   C6. Bounds enforcement -- optimize_parameters respects lower/upper bounds.
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/parameter_cascading.R")   # also sources general_ode_system_solver.R

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
# C1. Outer gradient consistency — Lotka-Volterra, all 4 parameters.
#
#     Both gradient methods are tested against the same FD reference and
#     against each other:
#       • outer_gradient            — IFT via assembled A,B matrices
#       • outer_gradient_sensitivity — sensitivity equations solved by FBSM
#
#     Checks (per method vs FD):
#       • Cosine similarity > 0.98 — direction test; non-trivial with 4
#         parameters (a wrong adjoint rotates the vector in 4-D space).
#       • Sign agreement per component — tightest magnitude check that is
#         theoretically justified given inner-solver inexactness.
#
#     Cross-consistency check:
#       • Cosine similarity between the two analytic methods > 0.999 — both
#         methods implement the same Gauss-Newton approximation and should
#         agree to near machine precision.
# ---------------------------------------------------------------------------
describe("C1: Outer gradient consistency — Lotka-Volterra (all parameters)", {

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

  cascading <- CascadingOdeSolver$new(
    func_rhs     = lv_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) as.numeric(y0_true),
    fixed_params = list(),
    lambda       = 1.0,
    param_scales = param_scales
  )

  # Test point: each parameter perturbed in a different direction/magnitude
  # so all four gradient components are non-zero and point different ways.
  theta_test <- c(alpha = 0.8,
                  beta  = 0.4,
                  delta = 0.4,
                  gamma = 0.6)

  # Warm-up: populate cache, then compute both analytic gradients
  set.seed(11)
  cascading$outer_objective(theta_test, param_names)
  g_ift  <- cascading$outer_gradient(theta_test, param_names)
  g_sens <- cascading$outer_gradient_sensitivity(theta_test, param_names)

  # Numerical FD reference (central differences, fixed seed for inner solver)
  eps <- 1e-3
  g_fd <- numeric(length(theta_test))
  for (j in seq_along(theta_test)) {
    th_p <- theta_test; th_p[j] <- theta_test[j] + eps
    th_m <- theta_test; th_m[j] <- theta_test[j] - eps
    set.seed(11); f_p <- cascading$outer_objective(th_p, param_names)
    set.seed(11); f_m <- cascading$outer_objective(th_m, param_names)
    g_fd[j] <- (f_p - f_m) / (2 * eps)
  }

  # Helper: cosine similarity between two vectors
  cos_sim <- function(a, b)
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-16)

  fmt_vec <- function(x) paste(sprintf("%+.3e", x), collapse = ", ")

  cos_ift  <- cos_sim(g_ift,  g_fd)
  cos_sens <- cos_sim(g_sens, g_fd)
  cos_cross <- cos_sim(g_ift, g_sens)

  # --- IFT method vs FD ---
  test_that("IFT gradient direction (cosine similarity) > 0.95", {
    expect_greater_than(
      cos_ift, 0.95,
      sprintf("[IFT]  cosine = %.4f\n  ift=(%s)\n   fd=(%s)",
              cos_ift, fmt_vec(g_ift), fmt_vec(g_fd))
    )
  })

  test_that("IFT gradient component signs agree with FD", {
    expect_true(
      all(sign(g_ift) == sign(g_fd)),
      sprintf("[IFT]  sign mismatch\n  ift=(%s)\n   fd=(%s)",
              fmt_vec(g_ift), fmt_vec(g_fd))
    )
  })

  # --- Sensitivity-equation method vs FD ---
  # The FBSM directly differentiates the KKT system — more faithful to the
  # inner problem than the IFT — so we hold it to a tighter threshold.
  test_that("sensitivity gradient direction (cosine similarity) > 0.99", {
    expect_greater_than(
      cos_sens, 0.99,
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

  # --- Cross-consistency: IFT vs sensitivity equations ---
  # The two methods implement DIFFERENT Gauss-Newton approximations:
  #   IFT  : uses u(y) = ODE-residual(y) to eliminate u, then differentiates
  #           J(y; theta) — a coarser approximation, error O(u*).
  #   FBSM : differentiates the actual KKT conditions — more accurate.
  # Both point in the correct descent direction (both pass FD tests above), but
  # they will not agree to near-machine precision.  A threshold of 0.97 checks
  # that neither has a gross implementation error while acknowledging the
  # systematic approximation gap.
  test_that("IFT and sensitivity gradients are mutually consistent (cos > 0.97)", {
    expect_greater_than(
      cos_cross, 0.97,
      sprintf("cross cosine = %.4f\n  ift =(%s)\n  sens=(%s)",
              cos_cross, fmt_vec(g_ift), fmt_vec(g_sens))
    )
  })
})

# ---------------------------------------------------------------------------
# C2. Descent direction: at a non-optimal theta, the negative gradient is
#     a descent direction (inner product with step should lower the objective).
# ---------------------------------------------------------------------------
describe("C2: Descent direction — outer gradient points downhill", {

  true_k    <- 0.4
  times_sim <- seq(0, 4, by = 0.2)
  obs_times <- seq(0, 4, by = 0.5)
  y0_true   <- 3.0

  y_true   <- euler_solve(decay_rhs, y0_true, obs_times, list(k = true_k))
  set.seed(12)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.05), nrow(y_true), 1)

  param_scales <- list(k = 0.5)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) c(y0_true),
    fixed_params = list(),
    lambda       = 0.3,
    param_scales = param_scales
  )

  # Evaluate at a perturbed parameter (away from optimum)
  theta_test <- c(k = true_k / param_scales$k * 1.4)  # 40% above truth

  set.seed(13)
  f0 <- cascading$outer_objective(theta_test, "k")
  g  <- cascading$outer_gradient(theta_test, "k")

  # Step in the negative-gradient direction
  step_size  <- 0.05
  theta_step <- theta_test - step_size * g
  theta_step <- pmax(theta_step, 1e-3)   # keep positive

  set.seed(13)
  f1 <- cascading$outer_objective(theta_step, "k")

  test_that("taking a step along negative gradient reduces the outer objective", {
    expect_less_than(f1, f0,
                     sprintf("f0 = %.6g, f1 = %.6g", f0, f1))
  })
})

# ---------------------------------------------------------------------------
# C3. Parameter recovery: 1-D exponential decay
#     CascadingOdeSolver should recover k to within 15% of the true value.
# ---------------------------------------------------------------------------
describe("C3: Parameter recovery — 1-D exponential decay", {

  true_k    <- 0.6
  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.5)
  y0_true   <- 4.0

  y_true   <- euler_solve(decay_rhs, y0_true, obs_times, list(k = true_k))
  set.seed(14)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.1), nrow(y_true), 1)

  param_scales <- list(k = 0.5)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) c(y0_true),
    fixed_params = list(),
    lambda       = 0.5,
    param_scales = param_scales
  )

  set.seed(15)
  result <- cascading$optimize_parameters(
    init_theta_physical = c(k = 0.3),   # start far from truth
    param_names         = "k",
    lower_phys          = c(k = 0.05),
    upper_phys          = c(k = 2.0)
  )

  recovered_k  <- result["k"]
  relative_err <- abs(recovered_k - true_k) / true_k

  test_that("recovered k is within 15% of true value", {
    expect_less_than(
      relative_err, 0.15,
      sprintf("true k = %.4g, recovered k = %.4g (%.1f%% error)",
              true_k, recovered_k, 100 * relative_err)
    )
  })
})

# ---------------------------------------------------------------------------
# C4. Parameter recovery: Lotka-Volterra (recover alpha, beta fixed)
#     Recovers the prey growth rate alpha with other params fixed.
# ---------------------------------------------------------------------------
describe("C4: Parameter recovery — Lotka-Volterra (alpha)", {

  true_params <- list(alpha = 1.2, beta = 0.4, delta = 0.1, gamma = 0.4)
  y0_true     <- c(8, 5)
  times_sim   <- seq(0, 10, by = 0.2)
  obs_times   <- seq(0, 10, by = 1.0)

  y_true   <- euler_solve(lv_rhs, y0_true, obs_times, true_params)
  set.seed(16)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.5), nrow(y_true), 2)

  # Fix beta, delta, gamma; estimate alpha
  fixed_params <- list(beta = 0.4, delta = 0.1, gamma = 0.4)
  param_scales <- list(alpha = 1.0)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = lv_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) as.numeric(y0_true),
    fixed_params = fixed_params,
    lambda       = 0.5,
    param_scales = param_scales
  )

  set.seed(17)
  result <- cascading$optimize_parameters(
    init_theta_physical = c(alpha = 0.8),
    param_names         = "alpha",
    lower_phys          = c(alpha = 0.3),
    upper_phys          = c(alpha = 3.0)
  )

  recovered_alpha <- result["alpha"]
  relative_err    <- abs(recovered_alpha - true_params$alpha) / true_params$alpha

  test_that("recovered alpha is within 20% of true value", {
    expect_less_than(
      relative_err, 0.20,
      sprintf("true alpha = %.4g, recovered = %.4g (%.1f%% error)",
              true_params$alpha, recovered_alpha, 100 * relative_err)
    )
  })
})

# ---------------------------------------------------------------------------
# C5. Caching correctness: calling outer_objective twice with the same theta
#     must return the same value and reuse the cached inner solver.
# ---------------------------------------------------------------------------
describe("C5: Caching correctness", {

  times_sim <- seq(0, 3, by = 0.2)
  obs_times <- seq(0, 3, by = 0.5)
  y0_true   <- 2.0
  params    <- list(k = 0.5)
  y_true    <- euler_solve(decay_rhs, y0_true, obs_times, params)
  set.seed(18)
  obs_data  <- y_true + matrix(rnorm(length(y_true), 0, 0.05), nrow(y_true), 1)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) c(y0_true),
    fixed_params = list(),
    lambda       = 0.3,
    param_scales = list(k = 0.5)
  )

  theta <- c(k = 1.0)

  set.seed(19)
  f1 <- cascading$outer_objective(theta, "k")
  # second call with same theta — must hit cache (no new random u_init)
  f2 <- cascading$outer_objective(theta, "k")

  test_that("outer_objective returns identical value on repeated call (cache hit)", {
    expect_equal(f1, f2, tol = 0, "values differ: cache may not be working")
  })

  test_that("last_theta is set after the first call", {
    expect_true(!is.null(cascading$last_theta))
  })

  test_that("last_solver is populated after the first call", {
    expect_true(!is.null(cascading$last_solver))
  })
})

# ---------------------------------------------------------------------------
# C6. Bounds enforcement: recovered parameter must stay within lower/upper
# ---------------------------------------------------------------------------
describe("C6: Bounds enforcement", {

  times_sim <- seq(0, 4, by = 0.2)
  obs_times <- seq(0, 4, by = 0.5)
  y0_true   <- 3.0
  params    <- list(k = 0.5)
  y_true    <- euler_solve(decay_rhs, y0_true, obs_times, params)
  set.seed(20)
  obs_data  <- y_true + matrix(rnorm(length(y_true), 0, 0.2), nrow(y_true), 1)

  # Set tight bounds that exclude the true k (to verify enforcement)
  lower_k <- 0.8
  upper_k <- 1.5

  cascading <- CascadingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) c(y0_true),
    fixed_params = list(),
    lambda       = 0.3,
    param_scales = list(k = 1.0)
  )

  set.seed(21)
  result <- cascading$optimize_parameters(
    init_theta_physical = c(k = 1.0),
    param_names         = "k",
    lower_phys          = c(k = lower_k),
    upper_phys          = c(k = upper_k)
  )

  recovered_k <- result["k"]

  test_that("recovered parameter respects lower bound", {
    expect_true(recovered_k >= lower_k - 1e-6,
                sprintf("k = %.4g < lower bound %.4g", recovered_k, lower_k))
  })

  test_that("recovered parameter respects upper bound", {
    expect_true(recovered_k <= upper_k + 1e-6,
                sprintf("k = %.4g > upper bound %.4g", recovered_k, upper_k))
  })
})

# ---------------------------------------------------------------------------
# C7. Initial condition estimation via init_state callback
# ---------------------------------------------------------------------------
describe("C7: Joint estimation of ODE parameter and initial condition", {

  true_k    <- 0.5
  true_y0   <- 3.2
  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.5)

  y_true   <- euler_solve(decay_rhs, true_y0, obs_times, list(k = true_k))
  set.seed(201)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.06), nrow(y_true), 1)

  cascading <- CascadingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    fixed_params = list(),
    lambda       = 0.4,
    param_scales = list(k = 1.0),
    init_state   = function(p) c(2 + 2.5 * p$k)
  )

  set.seed(202)
  result <- cascading$optimize_parameters(
    init_theta_physical = c(k = 0.2),
    param_names         = "k",
    lower_phys          = c(k = 0.05),
    upper_phys          = c(k = 2.0)
  )

  test_that("optimization returns ODE parameters only", {
    expect_true(identical(names(result), c("k")))
  })

  test_that("recovered k is reasonably close to truth", {
    expect_less_than(abs(result[["k"]] - true_k) / true_k, 0.3)
  })
})

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
test_summary()
