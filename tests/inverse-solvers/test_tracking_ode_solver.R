# =============================================================================
# tests/test_tracking_ode_solver.R
#
# Test suite for TrackingOdeSolver (src/solvers/tracking_ode_solver.R).
#
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/inverse-solvers/load_inverse_solvers.R")

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
# TR1. Outer gradient consistency — Lotka-Volterra, all 6 parameters.
#
#   Tests both gradient methods against finite differences and each other
#   across N_RAND random points sampled from the LV parameter domain.
#     outer_gradient             — adjoint formula
#     outer_gradient_sensitivity — forward sensitivity sweep
#
#   Checks per point (4 per point × 20 points = 80 total):
#     adj-cos:   adjoint cosine similarity with FD > 0.98
#     adj-sign:  adjoint sign agreement per component
#     sens-cos:  sensitivity cosine similarity with FD > 0.98
#     sens-sign: sensitivity sign agreement per component
# ---------------------------------------------------------------------------
describe("TR1: Outer gradient consistency — random Lotka-Volterra parameters", {

  N_RAND   <- 20
  N_CHECKS <- 4

  p_true    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4, x0 = 10, y0 = 10)
  times_sim <- seq(0, 5, by = 0.1)
  set.seed(10)
  obs_times <- sort(sample(times_sim, 10))
  y_true    <- euler_solve(lv_rhs, unlist(p_true[c("x0","y0")]), times_sim, p_true)
  obs_data  <- y_true[which(times_sim %in% obs_times), ] +
               matrix(rnorm(10 * 2, 0, 0.01), 10, 2)

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0,
                       x0 = 1.0, y0 = 1.0)

  make_tracking <- function() {
    TrackingOdeSolver$new(
      func_rhs     = lv_rhs,
      times_sim    = times_sim,
      obs_times    = obs_times,
      obs_values   = obs_data,
      init_state   = function(p) as.numeric(c(p$x0, p$y0)),
      fixed_params = list(),
      lambda       = 1e2,
      param_scales = param_scales
    )
  }

  cos_sim <- function(a, b)
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-16)
  fmt_vec <- function(x) paste(sprintf("%+.3e", x), collapse = ", ")

  eval_gradients <- function(tracking, theta_test, seed = 11L) {
    set.seed(seed)
    tracking$outer_objective(theta_test, names(theta_test))
    g_adj  <- tracking$outer_gradient(theta_test, names(theta_test))
    g_sens <- tracking$outer_gradient_sensitivity(theta_test, names(theta_test))

    eps    <- 1e-3
    g_fd   <- numeric(length(theta_test))
    for (j in seq_along(theta_test)) {
      th_p <- theta_test; th_p[j] <- theta_test[j] + eps
      th_m <- theta_test; th_m[j] <- theta_test[j] - eps
      set.seed(seed); f_p <- tracking$outer_objective(th_p, names(theta_test))
      set.seed(seed); f_m <- tracking$outer_objective(th_m, names(theta_test))
      g_fd[j] <- (f_p - f_m) / (2 * eps)
    }

    list(g_adj = g_adj, g_sens = g_sens, g_fd = g_fd,
         cos_adj = cos_sim(g_adj, g_fd),
         cos_sens = cos_sim(g_sens, g_fd))
  }

  set.seed(42)
  theta_candidates <- replicate(N_RAND, {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0),
      x0 = runif(1, 2, 20), y0 = runif(1, 2, 20))
  }, simplify = FALSE)

  passes_by_check <- list("adj-cos"=0L, "adj-sign"=0L, "sens-cos"=0L,
                          "sens-sign"=0L)
  failures_by_point <- list()

  for (i in seq_len(N_RAND)) {
    theta_test <- theta_candidates[[i]]
    tracking   <- make_tracking()
    r          <- eval_gradients(tracking, theta_test, seed = 100L + i)

    pt_label <- sprintf("pt%d(α=%.2f β=%.2f δ=%.2f γ=%.2f x0=%.1f y0=%.1f)",
                        i, theta_test["alpha"], theta_test["beta"],
                        theta_test["delta"], theta_test["gamma"],
                        theta_test["x0"], theta_test["y0"])
    failures_at_point <- character(0)

    check <- function(cond, tag, msg) {
      if (isTRUE(cond)) {
        passes_by_check[[tag]] <<- passes_by_check[[tag]] + 1L
      } else {
        failures_at_point <<- c(failures_at_point, sprintf("%s: %s", tag, msg))
      }
    }

    check(r$cos_adj > 0.98, "adj-cos",
          sprintf("cos=%.4f [expected >0.98]", r$cos_adj))
    check(all(sign(r$g_adj) == sign(r$g_fd)), "adj-sign",
          sprintf("mismatch: adj=(%s) vs fd=(%s)", fmt_vec(r$g_adj), fmt_vec(r$g_fd)))
    check(r$cos_sens > 0.98, "sens-cos",
          sprintf("cos=%.4f [expected >0.98]", r$cos_sens))
    check(all(sign(r$g_sens) == sign(r$g_fd)), "sens-sign",
          sprintf("mismatch: sens=(%s) vs fd=(%s)", fmt_vec(r$g_sens), fmt_vec(r$g_fd)))

    if (length(failures_at_point) > 0)
      failures_by_point[[pt_label]] <- failures_at_point
  }

  total_passes <- sum(unlist(passes_by_check))
  total        <- N_RAND * N_CHECKS

  breakdown <- sprintf("%-10s %2d/%2d (%3.0f%%)",
    names(passes_by_check),
    unlist(passes_by_check), N_RAND,
    100L * unlist(passes_by_check) / N_RAND)

  test_name <- sprintf("gradient checks across %d random points: %d/%d (%.0f%%)",
                       N_RAND, total_passes, total, 100L * total_passes / total)

  test_that(test_name, {
    if (length(failures_by_point) > 0) {
      fail_lines <- c("\nFailures by point:")
      for (pt in names(failures_by_point)) {
        fail_lines <- c(fail_lines, sprintf("  %s:", pt),
                        sprintf("    • %s", failures_by_point[[pt]]))
      }
      fail_lines <- c(fail_lines, "\nBreakdown by check type:", sprintf("  %s", breakdown))
      stop(paste(fail_lines, collapse = "\n"))
    }
  })
})

# ---------------------------------------------------------------------------
# TR1a. Outer gradient consistency — Lotka-Volterra with fixed IC.
#
#   Same as TR1 but with fixed initial conditions (x0=10, y0=10).
#   Tests both gradient methods against finite differences across N_RAND
#   random points sampled from the 4-parameter Lotka-Volterra domain.
#     outer_gradient             — adjoint formula
#     outer_gradient_sensitivity — forward sensitivity sweep
#
#   Checks per point (4 per point × 20 points = 80 total):
#     adj-cos:   adjoint cosine similarity with FD > 0.98
#     adj-sign:  adjoint sign agreement per component
#     sens-cos:  sensitivity cosine similarity with FD > 0.98
#     sens-sign: sensitivity sign agreement per component
# ---------------------------------------------------------------------------
describe("TR1a: Outer gradient consistency — fixed initial conditions", {

  N_RAND   <- 20
  N_CHECKS <- 4

  p_true    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4, x0 = 10, y0 = 10)
  times_sim <- seq(0, 5, by = 0.1)
  set.seed(10)
  obs_times <- sort(sample(times_sim, 10))
  y_true    <- euler_solve(lv_rhs, unlist(p_true[c("x0","y0")]), times_sim, p_true)
  obs_data  <- y_true[which(times_sim %in% obs_times), ] +
               matrix(rnorm(10 * 2, 0, 0.01), 10, 2)

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)

  make_tracking <- function() {
    TrackingOdeSolver$new(
      func_rhs     = lv_rhs,
      times_sim    = times_sim,
      obs_times    = obs_times,
      obs_values   = obs_data,
      init_state   = function(p) as.numeric(c(10, 10)),
      fixed_params = list(),
      lambda       = 1e2,
      param_scales = param_scales
    )
  }

  cos_sim <- function(a, b)
    sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-16)
  fmt_vec <- function(x) paste(sprintf("%+.3e", x), collapse = ", ")

  eval_gradients <- function(tracking, theta_test, seed = 11L) {
    set.seed(seed)
    tracking$outer_objective(theta_test, names(theta_test))
    g_adj  <- tracking$outer_gradient(theta_test, names(theta_test))
    g_sens <- tracking$outer_gradient_sensitivity(theta_test, names(theta_test))

    eps    <- 1e-3
    g_fd   <- numeric(length(theta_test))
    for (j in seq_along(theta_test)) {
      th_p <- theta_test; th_p[j] <- theta_test[j] + eps
      th_m <- theta_test; th_m[j] <- theta_test[j] - eps
      set.seed(seed); f_p <- tracking$outer_objective(th_p, names(theta_test))
      set.seed(seed); f_m <- tracking$outer_objective(th_m, names(theta_test))
      g_fd[j] <- (f_p - f_m) / (2 * eps)
    }

    list(g_adj = g_adj, g_sens = g_sens, g_fd = g_fd,
         cos_adj = cos_sim(g_adj, g_fd),
         cos_sens = cos_sim(g_sens, g_fd))
  }

  set.seed(42)
  theta_candidates <- replicate(N_RAND, {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0))
  }, simplify = FALSE)

  passes_by_check <- list("adj-cos"=0L, "adj-sign"=0L, "sens-cos"=0L,
                          "sens-sign"=0L)
  failures_by_point <- list()

  for (i in seq_len(N_RAND)) {
    theta_test <- theta_candidates[[i]]
    tracking   <- make_tracking()
    r          <- eval_gradients(tracking, theta_test, seed = 100L + i)

    pt_label <- sprintf("pt%d(α=%.2f β=%.2f δ=%.2f γ=%.2f)",
                        i, theta_test["alpha"], theta_test["beta"],
                        theta_test["delta"], theta_test["gamma"])
    failures_at_point <- character(0)

    check <- function(cond, tag, msg) {
      if (isTRUE(cond)) {
        passes_by_check[[tag]] <<- passes_by_check[[tag]] + 1L
      } else {
        failures_at_point <<- c(failures_at_point, sprintf("%s: %s", tag, msg))
      }
    }

    check(r$cos_adj > 0.98, "adj-cos",
          sprintf("cos=%.4f [expected >0.98]", r$cos_adj))
    check(all(sign(r$g_adj) == sign(r$g_fd)), "adj-sign",
          sprintf("mismatch: adj=(%s) vs fd=(%s)", fmt_vec(r$g_adj), fmt_vec(r$g_fd)))
    check(r$cos_sens > 0.98, "sens-cos",
          sprintf("cos=%.4f [expected >0.98]", r$cos_sens))
    check(all(sign(r$g_sens) == sign(r$g_fd)), "sens-sign",
          sprintf("mismatch: sens=(%s) vs fd=(%s)", fmt_vec(r$g_sens), fmt_vec(r$g_fd)))

    if (length(failures_at_point) > 0)
      failures_by_point[[pt_label]] <- failures_at_point
  }

  total_passes <- sum(unlist(passes_by_check))
  total        <- N_RAND * N_CHECKS

  breakdown <- sprintf("%-10s %2d/%2d (%3.0f%%)",
    names(passes_by_check),
    unlist(passes_by_check), N_RAND,
    100L * unlist(passes_by_check) / N_RAND)

  test_name <- sprintf("gradient checks (4 params) across %d random points: %d/%d (%.0f%%)",
                       N_RAND, total_passes, total, 100L * total_passes / total)

  test_that(test_name, {
    if (length(failures_by_point) > 0) {
      fail_lines <- c("\nFailures by point:")
      for (pt in names(failures_by_point)) {
        fail_lines <- c(fail_lines, sprintf("  %s:", pt),
                        sprintf("    • %s", failures_by_point[[pt]]))
      }
      fail_lines <- c(fail_lines, "\nBreakdown by check type:", sprintf("  %s", breakdown))
      stop(paste(fail_lines, collapse = "\n"))
    }
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
    init_state   = function(p) c(y0_true),
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

  times_sim <- seq(0, 3, by = 0.1)
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
    init_state   = function(p) c(y0_true),
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
# TR6. Initial condition estimation via init_state callback
# ---------------------------------------------------------------------------
describe("TR6: Joint estimation of ODE parameter and initial condition", {

  true_k    <- 0.55
  true_y0   <- 3.5
  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.5)

  y_true   <- euler_solve(decay_rhs, true_y0, obs_times, list(k = true_k))
  set.seed(101)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.05), nrow(y_true), 1)

  tracking <- TrackingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    fixed_params = list(),
    lambda       = 0.4,
    param_scales = list(k = 1.0),
    init_state   = function(p) c(2 + 3 * p$k)
  )

  set.seed(102)
  result <- tracking$optimize_parameters(
    init_theta_physical = c(k = 0.25),
    param_names         = "k",
    lower_phys          = c(k = 0.05),
    upper_phys          = c(k = 2.0),
    gradient_mode       = "adjoint"
  )

  test_that("optimization returns ODE parameters only", {
    expect_true(identical(names(result), c("k")))
  })

  test_that("recovered k is reasonably close to truth", {
    expect_less_than(abs(result[["k"]] - true_k) / true_k, 0.25)
  })
})

# ---------------------------------------------------------------------------
# TR7. Runtime gradient mode switch with init_state callback
# ---------------------------------------------------------------------------
describe("TR7: Gradient mode switch supports sensitivity mode", {

  true_k    <- 0.4
  true_y0   <- 2.8
  times_sim <- seq(0, 4, by = 0.1)
  obs_times <- seq(0, 4, by = 0.4)

  y_true   <- euler_solve(decay_rhs, true_y0, obs_times, list(k = true_k))
  set.seed(103)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.04), nrow(y_true), 1)

  tracking <- TrackingOdeSolver$new(
    func_rhs     = decay_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    fixed_params = list(),
    lambda       = 1e2,
    param_scales = list(k = 1.0, y0 = 1.0),
    init_state   = function(p) c(p$y0)
  )

  result <- tracking$optimize_parameters(
    init_theta_physical = c(k = 0.7, y0 = 1.5),
    param_names         = c("k", "y0"),
    lower_phys          = c(k = 0.05, y0 = 0.1),
    upper_phys          = c(k = 2.0,  y0 = 10.0),
    gradient_mode       = "sensitivity"
  )

  test_that("sensitivity mode returns finite parameter estimate", {
    expect_true(all(is.finite(result)))
    expect_true(identical(names(result), c("k","y0")))
  })
})

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
test_summary()
