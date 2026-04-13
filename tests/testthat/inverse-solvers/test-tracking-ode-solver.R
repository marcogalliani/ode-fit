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
  adj_chk <- safe_chk(
    function(theta, param_names) tracking$outer_gradient(theta, param_names)
  )

  set.seed(seed)
  sens_chk <- safe_chk(
    function(theta, param_names) tracking$outer_gradient_sensitivity(theta, param_names)
  )

  list(
    checks = c(
      adj_rel_error = adj_chk$max_rel_error < 1e-3,
      adj_cos       = adj_chk$cosine_similarity > 0.98,
      sens_rel_error = sens_chk$max_rel_error < 1e-3,
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
describe("TR1: Outer gradient consistency — random Lotka-Volterra parameters", {

  N_RAND <- 20
  MAX_ATTEMPTS <- 100
  times_sim <- seq(0, 5, by = 0.1)
  n_obs <- 30
  set.seed(10)
  obs_times <- sort(sample(times_sim, n_obs))

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0,
                       x0 = 1.0, y0 = 1.0)

  random_true_candidate <- function() {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0),
      x0 = runif(1, 2, 20), y0 = runif(1, 2, 20))
  }
  random_theta_candidate <- function() {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0),
      x0 = runif(1, 2, 20), y0 = runif(1, 2, 20))
  }

  check_names <- c("adj_rel_error", "adj_cos", "sens_rel_error", "sens_cos")
  test_results <- matrix(FALSE, nrow = length(check_names), ncol = 0,
                         dimnames = list(check_names, NULL))
  trial_details <- list()
  unstable_trials <- list()

  sink_file <- tempfile()
  sink(sink_file)
  on.exit(sink(), add = TRUE)

  attempts <- 0L
  stable_idx <- 0L
  while (stable_idx < N_RAND && attempts < MAX_ATTEMPTS) {
    attempts <- attempts + 1L
    p_true <- as.list(random_true_candidate())
    theta_test <- random_theta_candidate()
    y_true <- euler_solve(lv_rhs, unlist(p_true[c("x0", "y0")]), times_sim, p_true)
    set.seed(1000L + attempts)
    obs_data <- y_true[which(times_sim %in% obs_times), ] # + matrix(rnorm(n_obs * 2, 0, 0.01), n_obs, 2)

    tracking <- TrackingOdeSolver$new(
      func_rhs     = lv_rhs,
      times_sim    = times_sim,
      obs_times    = obs_times,
      obs_values   = obs_data,
      init_state   = function(p) as.numeric(c(p$x0, p$y0)),
      fixed_params = list(),
      lambda       = 1e2,
      param_scales = param_scales
    )

    stability <- filter_stable_info(theta_test, function(theta) {
      tracking$outer_objective(theta, names(theta))
    })
    if (!stability$stable) {
      unstable_trials[[length(unstable_trials) + 1L]] <- list(
        p_true = unlist(p_true),
        theta_test = theta_test,
        reason = stability$reason,
        value = stability$value
      )
      next
    }

    stable_idx <- stable_idx + 1L
    chk <- check_outer_gradient(tracking, theta_test, seed = 100L + attempts)
    test_results <- cbind(test_results, chk$checks)
    trial_details[[stable_idx]] <- list(
      theta_test = theta_test,
      p_true = unlist(p_true),
      metrics = chk$metrics,
      errors = chk$errors
    )
  }

  if (stable_idx < N_RAND) {
    stop(sprintf("Could not collect %d stable trials after %d attempts.", N_RAND, MAX_ATTEMPTS))
  }

  for (chk_cond in rownames(test_results)) {
    failed_idx <- which(!test_results[chk_cond, ])
    test_that(sprintf("%s: all random trials pass", chk_cond), {
      if (length(failed_idx) > 0) {
        fail_lines <- vapply(failed_idx, function(idx) {
          det <- trial_details[[idx]]
          sprintf(
            "trial %d | true=(%s) | theta=(%s) | adj_rel=%.3e | adj_cos=%.6f | sens_rel=%.3e | sens_cos=%.6f",
            idx,
            fmt_named(det$p_true),
            fmt_named(det$theta_test),
            det$metrics[["adj_max_rel_error"]],
            det$metrics[["adj_cosine"]],
            det$metrics[["sens_max_rel_error"]],
            det$metrics[["sens_cosine"]]
          )
        }, character(1L))
        err_lines <- vapply(failed_idx, function(idx) {
          det <- trial_details[[idx]]
          sprintf("trial %d errors | adj=%s | sens=%s", idx,
                  det$errors[["adj_error"]], det$errors[["sens_error"]])
        }, character(1L))
        stop(paste(c("Failing trials:", fail_lines, "Solver errors:", err_lines), collapse = "\n"))
      }
      expect_true(TRUE)
    })
  }
})

# ---------------------------------------------------------------------------
# TR1a. Outer gradient consistency — fixed IC, 4 parameters
# ---------------------------------------------------------------------------
describe("TR1a: Outer gradient consistency — fixed initial conditions", {

  N_RAND <- 20
  MAX_ATTEMPTS <- 100
  times_sim <- seq(0, 5, by = 0.1)
  set.seed(10)
  n_obs <- 20
  obs_times <- sort(sample(times_sim, n_obs))

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)

  random_true_candidate <- function() {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0))
  }
  random_theta_candidate <- function() {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0))
  }

  check_names <- c("adj_rel_error", "adj_cos", "sens_rel_error", "sens_cos")
  test_results <- matrix(FALSE, nrow = length(check_names), ncol = 0,
                         dimnames = list(check_names, NULL))
  trial_details <- list()
  unstable_trials <- list()

  sink_file <- tempfile()
  sink(sink_file)
  on.exit(sink(), add = TRUE)

  attempts <- 0L
  stable_idx <- 0L
  while (stable_idx < N_RAND && attempts < MAX_ATTEMPTS) {
    attempts <- attempts + 1L
    p_true <- as.list(random_true_candidate())
    theta_test <- random_theta_candidate()
    p_data <- c(p_true, list(x0 = 10, y0 = 10))
    y_true <- euler_solve(lv_rhs, unlist(p_data[c("x0", "y0")]), times_sim, p_data)
    set.seed(2000L + attempts)
    obs_data <- y_true[which(times_sim %in% obs_times), ] #+ matrix(rnorm(n_obs * 2, 0, 0.01), n_obs, 2)

    tracking <- TrackingOdeSolver$new(
      func_rhs     = lv_rhs,
      times_sim    = times_sim,
      obs_times    = obs_times,
      obs_values   = obs_data,
      init_state   = function(p) as.numeric(c(10, 10)),
      fixed_params = list(),
      lambda       = 1e2,
      param_scales = param_scales
    )

    stability <- filter_stable_info(theta_test, function(theta) {
      tracking$outer_objective(theta, names(theta))
    })
    if (!stability$stable) {
      unstable_trials[[length(unstable_trials) + 1L]] <- list(
        p_true = unlist(p_true),
        theta_test = theta_test,
        reason = stability$reason,
        value = stability$value
      )
      next
    }

    stable_idx <- stable_idx + 1L
    chk <- check_outer_gradient(tracking, theta_test, seed = 100L + attempts)
    test_results <- cbind(test_results, chk$checks)
    trial_details[[stable_idx]] <- list(
      theta_test = theta_test,
      p_true = unlist(p_true),
      metrics = chk$metrics,
      errors = chk$errors
    )
  }

  if (stable_idx < N_RAND) {
    stop(sprintf("Could not collect %d stable trials after %d attempts.", N_RAND, MAX_ATTEMPTS))
  }

  for (chk_cond in rownames(test_results)) {
    failed_idx <- which(!test_results[chk_cond, ])
    test_that(sprintf("%s: all random trials pass", chk_cond), {
      if (length(failed_idx) > 0) {
        fail_lines <- vapply(failed_idx, function(idx) {
          det <- trial_details[[idx]]
          sprintf(
            "trial %d | true=(%s) | theta=(%s) | adj_rel=%.3e | adj_cos=%.6f | sens_rel=%.3e | sens_cos=%.6f",
            idx,
            fmt_named(det$p_true),
            fmt_named(det$theta_test),
            det$metrics[["adj_max_rel_error"]],
            det$metrics[["adj_cosine"]],
            det$metrics[["sens_max_rel_error"]],
            det$metrics[["sens_cosine"]]
          )
        }, character(1L))
        err_lines <- vapply(failed_idx, function(idx) {
          det <- trial_details[[idx]]
          sprintf("trial %d errors | adj=%s | sens=%s", idx,
                  det$errors[["adj_error"]], det$errors[["sens_error"]])
        }, character(1L))
        stop(paste(c("Failing trials:", fail_lines, "Solver errors:", err_lines), collapse = "\n"))
      }
      expect_true(TRUE)
    })
  }
})
