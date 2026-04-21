# =============================================================================
# Test suite for CascadingOdeSolver (src/solvers/parameter_cascading.R).
#
# WHAT IS TESTED:
#   - Outer gradient consistency
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

check_outer_gradient <- function(cascading, theta_test, param_names, seed = 11L) {
  safe_chk <- function(grad_fun) {
    tryCatch(
      check_gradient(
        fn  = function(theta, param_names) cascading$outer_objective(theta, param_names),
        gr  = grad_fun,
        par = theta_test,
        eps = 1e-3,
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
  ift_chk <- safe_chk(
    function(theta, param_names) cascading$outer_gradient(theta, param_names)
  )

  set.seed(seed)
  sens_chk <- safe_chk(
    function(theta, param_names) cascading$outer_gradient_sensitivity(theta, param_names)
  )

  list(
    ift_chk = ift_chk,
    sens_chk = sens_chk,
    errors = c(
      ift_error = if (!is.null(ift_chk$error_msg)) ift_chk$error_msg else "",
      sens_error = if (!is.null(sens_chk$error_msg)) sens_chk$error_msg else ""
    )
  )
}

# ---------------------------------------------------------------------------
# Outer gradient consistency — Lotka-Volterra, all 4 parameters
# ---------------------------------------------------------------------------
describe("Outer gradient consistency — deterministic Lotka-Volterra (all parameters)", {

  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.5)

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)
  param_names  <- c("alpha", "beta", "delta", "gamma")

  p_true <- c(alpha = 1.20, beta = 0.45, delta = 0.12, gamma = 0.80,
              x0 = 10.0, y0 = 8.0)
  theta_test <- c(alpha = 1.10, beta = 0.50, delta = 0.10, gamma = 0.75)
  y0_true <- as.numeric(c(p_true[["x0"]], p_true[["y0"]]))

  y_true <- euler_solve(lv_rhs, y0_true, obs_times, as.list(p_true))
  obs_data <- y_true

  model <- ODEModel$new(
    rhs = lv_rhs,
    init_state = function(p) y0_true,
    fixed_params = list(x0 = y0_true[1], y0 = y0_true[2]),
    param_scales = param_scales
  )

  cascading <- CascadingOdeSolver$new(
    model        = model,
    times_sim    = times_sim,
    obs_times    = obs_times,
    obs_values   = obs_data,
    lambda       = 1.0,
    inner_method = "gl2"
  )

  stability <- filter_stable_info(theta_test, function(theta) {
    cascading$outer_objective(theta, param_names)
  })
  if (!stability$stable) {
    stop(sprintf("Deterministic cascading set is unstable: %s (value=%s)",
                 stability$reason, as.character(stability$value)))
  }

  chk <- check_outer_gradient(cascading, theta_test, param_names, seed = 311L)

  test_that("Cascading deterministic set passes IFT relative error", {
    expect_true(chk$ift_chk$max_rel_error < 1.0,
                info = sprintf("theta=(%s) | ift_rel=%.3e | err=%s",
                               fmt_named(theta_test),
                               chk$ift_chk$max_rel_error,
                               chk$errors[["ift_error"]]))
  })
  test_that("Cascading deterministic set passes IFT cosine similarity", {
    expect_true(chk$ift_chk$cosine_similarity > 0.95,
                info = sprintf("theta=(%s) | ift_cos=%.6f | err=%s",
                               fmt_named(theta_test),
                               chk$ift_chk$cosine_similarity,
                               chk$errors[["ift_error"]]))
  })
  test_that("Cascading deterministic set passes IFT descent", {
    expect_true(chk$ift_chk$is_descent,
                info = sprintf("theta=(%s) | ift_desc=%s | err=%s",
                               fmt_named(theta_test),
                               as.character(chk$ift_chk$is_descent),
                               chk$errors[["ift_error"]]))
  })
  test_that("Cascading deterministic set passes sensitivity relative error", {
    expect_true(chk$sens_chk$max_rel_error < 1.0,
                info = sprintf("theta=(%s) | sens_rel=%.3e | err=%s",
                               fmt_named(theta_test),
                               chk$sens_chk$max_rel_error,
                               chk$errors[["sens_error"]]))
  })
  test_that("Cascading deterministic set passes sensitivity cosine similarity", {
    expect_true(chk$sens_chk$cosine_similarity > 0.99,
                info = sprintf("theta=(%s) | sens_cos=%.6f | err=%s",
                               fmt_named(theta_test),
                               chk$sens_chk$cosine_similarity,
                               chk$errors[["sens_error"]]))
  })
  test_that("Cascading deterministic set passes sensitivity descent", {
    expect_true(chk$sens_chk$is_descent,
                info = sprintf("theta=(%s) | sens_desc=%s | err=%s",
                               fmt_named(theta_test),
                               as.character(chk$sens_chk$is_descent),
                               chk$errors[["sens_error"]]))
  })
})
