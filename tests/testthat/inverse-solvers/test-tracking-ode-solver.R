# =============================================================================
# Test suite for TrackingOdeSolver (src/solvers/tracking_ode_solver.R).
#
# WHAT IS TESTED:
#   TR1.  Outer gradient consistency — random Lotka-Volterra (6 params)
#   TR1a. Outer gradient consistency — fixed IC (4 params)
#   TR2.  Descent direction
#   TR7.  Gradient mode switch
# =============================================================================

setwd("../../..")
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

fmt_named <- function(x, digits = 3) {
  paste(sprintf("%s=%.*f", names(x), digits, as.numeric(x)), collapse = ", ")
}

check_outer_gradient <- function(tracking, theta_test, seed = 11L) {
  param_names <- names(theta_test)

  safe_chk <- function(grad_fun) {
    tryCatch(
      check_gradient(
        fn  = function(theta, param_names) tracking$outer_objective(theta, param_names),
        gr  = grad_fun,
        par = theta_test,
        eps = 1e-4,
        param_names = param_names
      ),
      error = function(e) {
        list(
          max_rel_error = Inf,
          cosine_similarity = NA_real_,
          is_descent = FALSE,
          error_msg = conditionMessage(e)
        )
      }
    )
  }

  set.seed(seed)
  adj_chk <- safe_chk(
    function(theta, param_names) tracking$outer_gradient(theta, param_names)
  )

  set.seed(seed)
  sens_chk <- safe_chk(
    function(theta, param_names) tracking$outer_gradient_sensitivity(theta, param_names)
  )

  list(
    checks = c(
      adj_rel_error = adj_chk$max_rel_error < 1e-2,
      adj_cos       = adj_chk$cosine_similarity > 0.98,
      sens_rel_error = sens_chk$max_rel_error < 1e-2,
      sens_cos      = sens_chk$cosine_similarity > 0.98
    ),
    metrics = c(
      adj_max_rel_error = adj_chk$max_rel_error,
      adj_cosine        = adj_chk$cosine_similarity,
      sens_max_rel_error = sens_chk$max_rel_error,
      sens_cosine       = sens_chk$cosine_similarity
    ),
    errors = c(
      adj_error = if (!is.null(adj_chk$error_msg)) adj_chk$error_msg else "",
      sens_error = if (!is.null(sens_chk$error_msg)) sens_chk$error_msg else ""
    )
  )
}

# ---------------------------------------------------------------------------
# TR1. Outer gradient consistency — Lotka-Volterra, all 6 parameters
# ---------------------------------------------------------------------------
describe("TR1: Outer gradient consistency — deterministic Lotka-Volterra parameters", {

  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.2)

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0,
                       x0 = 1.0, y0 = 1.0)

  p_true <- c(alpha = 1.20, beta = 0.45, delta = 0.12, gamma = 0.80,
              x0 = 8.0, y0 = 6.0)
  theta_test <- c(alpha = 1.10, beta = 0.50, delta = 0.10, gamma = 0.75,
                  x0 = 7.5, y0 = 6.5)

  y_true <- euler_solve(lv_rhs, unlist(p_true[c("x0", "y0")]), times_sim, as.list(p_true))
  obs_data <- y_true[which(times_sim %in% obs_times), , drop = FALSE]

  tracking <- TrackingOdeSolver$new(
    func_rhs     = lv_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) as.numeric(c(p$x0, p$y0)),
    fixed_params = list(),
    lambda       = 1e2,
    param_scales = param_scales,
    inner_method = "gl1"
  )

  stability <- filter_stable_info(theta_test, function(theta) {
    tracking$outer_objective(theta, names(theta))
  })
  if (!stability$stable) {
    stop(sprintf("Deterministic TR1 set is unstable: %s (value=%s)",
                 stability$reason, as.character(stability$value)))
  }

  chk <- check_outer_gradient(tracking, theta_test, seed = 111L)

  test_that("TR1 deterministic set passes adjoint relative error", {
    expect_true(chk$checks[["adj_rel_error"]],
                info = sprintf("theta=(%s) | adj_rel=%.3e | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["adj_max_rel_error"]],
                               chk$errors[["adj_error"]]))
  })
  test_that("TR1 deterministic set passes adjoint cosine similarity", {
    expect_true(chk$checks[["adj_cos"]],
                info = sprintf("theta=(%s) | adj_cos=%.6f | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["adj_cosine"]],
                               chk$errors[["adj_error"]]))
  })
  test_that("TR1 deterministic set passes sensitivity relative error", {
    expect_true(chk$checks[["sens_rel_error"]],
                info = sprintf("theta=(%s) | sens_rel=%.3e | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["sens_max_rel_error"]],
                               chk$errors[["sens_error"]]))
  })
  test_that("TR1 deterministic set passes sensitivity cosine similarity", {
    expect_true(chk$checks[["sens_cos"]],
                info = sprintf("theta=(%s) | sens_cos=%.6f | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["sens_cosine"]],
                               chk$errors[["sens_error"]]))
  })
})

# ---------------------------------------------------------------------------
# TR1a. Outer gradient consistency — fixed IC, 4 parameters
# ---------------------------------------------------------------------------
describe("TR1a: Outer gradient consistency — deterministic fixed initial conditions", {

  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.2)

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)

  p_true <- c(alpha = 1.20, beta = 0.45, delta = 0.12, gamma = 0.80)
  theta_test <- c(alpha = 1.10, beta = 0.50, delta = 0.10, gamma = 0.75)
  p_data <- c(as.list(p_true), list(x0 = 10, y0 = 10))

  y_true <- euler_solve(lv_rhs, unlist(p_data[c("x0", "y0")]), times_sim, p_data)
  obs_data <- y_true[which(times_sim %in% obs_times), , drop = FALSE]

  tracking <- TrackingOdeSolver$new(
    func_rhs     = lv_rhs,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    init_state   = function(p) as.numeric(c(10, 10)),
    fixed_params = list(),
    lambda       = 1e2,
    param_scales = param_scales,
    inner_method = "euler"
  )

  stability <- filter_stable_info(theta_test, function(theta) {
    tracking$outer_objective(theta, names(theta))
  })
  if (!stability$stable) {
    stop(sprintf("Deterministic TR1a set is unstable: %s (value=%s)",
                 stability$reason, as.character(stability$value)))
  }

  chk <- check_outer_gradient(tracking, theta_test, seed = 211L)

  test_that("TR1a deterministic set passes adjoint relative error", {
    expect_true(chk$checks[["adj_rel_error"]],
                info = sprintf("theta=(%s) | adj_rel=%.3e | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["adj_max_rel_error"]],
                               chk$errors[["adj_error"]]))
  })
  test_that("TR1a deterministic set passes adjoint cosine similarity", {
    expect_true(chk$checks[["adj_cos"]],
                info = sprintf("theta=(%s) | adj_cos=%.6f | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["adj_cosine"]],
                               chk$errors[["adj_error"]]))
  })
  test_that("TR1a deterministic set passes sensitivity relative error", {
    expect_true(chk$checks[["sens_rel_error"]],
                info = sprintf("theta=(%s) | sens_rel=%.3e | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["sens_max_rel_error"]],
                               chk$errors[["sens_error"]]))
  })
  test_that("TR1a deterministic set passes sensitivity cosine similarity", {
    expect_true(chk$checks[["sens_cos"]],
                info = sprintf("theta=(%s) | sens_cos=%.6f | err=%s",
                               fmt_named(theta_test),
                               chk$metrics[["sens_cosine"]],
                               chk$errors[["sens_error"]]))
  })
})
