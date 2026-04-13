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
describe("Outer gradient consistency — Lotka-Volterra (all parameters)", {

  N_RAND <- 20
  MAX_ATTEMPTS <- 100

  times_sim <- seq(0, 5, by = 0.1)
  obs_times <- seq(0, 5, by = 0.5)

  param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)
  param_names  <- c("alpha", "beta", "delta", "gamma")

  random_true_candidate <- function() {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0),
      x0 = runif(1, 2, 20), y0 = runif(1, 2, 20))
  }
  random_theta_candidate <- function() {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0))
  }

  check_names <- c("ift_rel_error", "ift_cos", "ift_descent",
                   "sens_rel_error", "sens_cos", "sens_descent")
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
    y0_true <- as.numeric(c(p_true$x0, p_true$y0))

    y_true <- euler_solve(lv_rhs, y0_true, obs_times, p_true)
    set.seed(100L + attempts)
    obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.3),
                                nrow(y_true), ncol(y_true))

    cascading <- CascadingOdeSolver$new(
      func_rhs     = lv_rhs,
      times_sim    = times_sim,
      obs_times    = obs_times,
      obs_values   = obs_data,
      init_state   = function(p) y0_true,
      fixed_params = list(),
      lambda       = 1.0,
      param_scales = param_scales
    )

    stability <- filter_stable_info(theta_test, function(theta) {
      cascading$outer_objective(theta, param_names)
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
    chk <- check_outer_gradient(cascading, theta_test, param_names, seed = 200L + attempts)
    test_results <- cbind(test_results, c(
      ift_rel_error = chk$ift_chk$max_rel_error < 1e-3,
      ift_cos       = chk$ift_chk$cosine_similarity > 0.95,
      ift_descent   = chk$ift_chk$is_descent,
      sens_rel_error = chk$sens_chk$max_rel_error < 1e-3,
      sens_cos      = chk$sens_chk$cosine_similarity > 0.99,
      sens_descent  = chk$sens_chk$is_descent
    ))
    trial_details[[stable_idx]] <- list(
      p_true = unlist(p_true),
      theta_test = theta_test,
      ift_max_rel_error = chk$ift_chk$max_rel_error,
      ift_cosine = chk$ift_chk$cosine_similarity,
      ift_descent = chk$ift_chk$is_descent,
      sens_max_rel_error = chk$sens_chk$max_rel_error,
      sens_cosine = chk$sens_chk$cosine_similarity,
      sens_descent = chk$sens_chk$is_descent,
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
            "trial %d | true=(%s) | theta=(%s) | ift_rel=%.3e | ift_cos=%.6f | ift_desc=%s | sens_rel=%.3e | sens_cos=%.6f | sens_desc=%s",
            idx,
            fmt_named(det$p_true),
            fmt_named(det$theta_test),
            det$ift_max_rel_error,
            det$ift_cosine,
            as.character(det$ift_descent),
            det$sens_max_rel_error,
            det$sens_cosine,
            as.character(det$sens_descent)
          )
        }, character(1L))
        err_lines <- vapply(failed_idx, function(idx) {
          det <- trial_details[[idx]]
          sprintf("trial %d errors | ift=%s | sens=%s", idx,
                  det$errors[["ift_error"]], det$errors[["sens_error"]])
        }, character(1L))
        stop(paste(c("Failing trials:", fail_lines, "Solver errors:", err_lines), collapse = "\n"))
      }
      expect_true(TRUE)
    })
  }
})
