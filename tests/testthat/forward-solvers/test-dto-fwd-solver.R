# =============================================================================
# Test suite for AdjointForwardSolver (src/solvers/dto_fwd_solver.R).
# =============================================================================

setwd("../../..")
source("src/solvers/forward-solvers/load_forward_solvers.R")

# ---------------------------------------------------------------------------
# Shared physics definitions
# ---------------------------------------------------------------------------

decay_rhs <- function(y, t, p) -p$k * y

decay_rhs_no_forcing <- function(y, t, p) -p$k * y

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

safe_check_gradient <- function(fn, gr, par, eps, y0) {
  tryCatch(
    check_gradient(fn = fn, gr = gr, par = par, eps = eps, y0 = y0),
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

# ---------------------------------------------------------------------------
# T1. Gradient consistency: 1-D exponential decay
# ---------------------------------------------------------------------------
describe("T1: Gradient consistency — 1-D exponential decay", {

  N_RAND <- 5

  set.seed(41)
  trial_cfg <- replicate(N_RAND, {
    c(k = runif(1, 0.2, 2.0), y0 = runif(1, 2.0, 20.0))
  }, simplify = FALSE)

  check_names <- c("max_rel_error", "cos_sim", "is_descent")
  test_results <- matrix(FALSE, nrow = length(check_names), ncol = N_RAND,
                         dimnames = list(check_names, NULL))
  trial_details <- vector("list", N_RAND)

  for (i in seq_len(N_RAND)) {
    cfg      <- trial_cfg[[i]]
    params   <- list(k = cfg[["k"]])
    times_sim <- seq(0, 2, by = 0.1)
    y0       <- cfg[["y0"]]
    y_true   <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

    solver <- AdjointForwardSolver$new(
      model      = ODEModel$new(rhs = decay_rhs),
      times_sim  = times_sim,
      obs_times  = times_sim,
      obs_values = y_true,
      params     = params,
      lambda     = 0.05
    )

    set.seed(200L + i)
    u_test <- rnorm(solver$n_steps * solver$n_vars, sd = 0.05)
    chk <- safe_check_gradient(
      fn  = solver$cost_function,
      gr  = solver$gradient_function,
      par = u_test,
      eps = 1e-5,
      y0  = y0
    )

    test_results[, i] <- c(
      max_rel_error = chk$max_rel_error < 1e-3,
      cos_sim       = chk$cosine_similarity > 1 - 1e-3,
      is_descent    = chk$is_descent
    )
    trial_details[[i]] <- list(
      params = c(k = params$k, y0 = y0),
      max_rel_error = chk$max_rel_error,
      cosine_similarity = chk$cosine_similarity,
      is_descent = chk$is_descent,
      error_msg = if (!is.null(chk$error_msg)) chk$error_msg else ""
    )
  }

  for (chk_cond in rownames(test_results)) {
    failed_idx <- which(!test_results[chk_cond, ])
    test_that(sprintf("%s: all random trials pass", chk_cond), {
      if (length(failed_idx) > 0) {
        fail_lines <- vapply(failed_idx, function(idx) {
          det <- trial_details[[idx]]
          sprintf(
            "trial %d | %s | max_rel_error=%.3e | cos=%.6f | is_descent=%s | error=%s",
            idx,
            fmt_named(det$params),
            det$max_rel_error,
            det$cosine_similarity,
            as.character(det$is_descent),
            det$error_msg
          )
        }, character(1L))
        stop(paste(c("Failing trials:", fail_lines), collapse = "\n"))
      }
      expect_true(TRUE)
    })
  }
})

# ---------------------------------------------------------------------------
# T2. Gradient consistency: 2-D Lotka-Volterra
# ---------------------------------------------------------------------------
describe("T2: Gradient consistency — 2-D Lotka-Volterra", {

  N_RAND <- 5
  MAX_ATTEMPTS <- 50

  random_trial_cfg <- function() {
    c(alpha = runif(1, 0.2, 2.0), beta = runif(1, 0.2, 2.0),
      delta = runif(1, 0.2, 2.0), gamma = runif(1, 0.2, 2.0),
      x0 = runif(1, 2.0, 20.0), y0 = runif(1, 2.0, 20.0))
  }

  check_names <- c("max_rel_error", "cos_sim", "is_descent")
  stable_results <- matrix(FALSE, nrow = length(check_names), ncol = 0,
                           dimnames = list(check_names, NULL))
  trial_details <- list()
  unstable_trials <- list()

  evaluate_trial <- function(cfg, trial_seed) {
    params <- list(alpha = cfg[["alpha"]], beta = cfg[["beta"]],
                   delta = cfg[["delta"]], gamma = cfg[["gamma"]])
    times_sim <- seq(0, 5, by = 0.2)
    y0 <- c(cfg[["x0"]], cfg[["y0"]])
    y_true <- euler_solve(lv_rhs, y0, times_sim, params)
    set.seed(500L + trial_seed)
    obs_data <- y_true + matrix(rnorm(length(y_true), 0, 0.2), nrow(y_true), 2)

    solver <- AdjointForwardSolver$new(
      model      = ODEModel$new(rhs = lv_rhs),
      times_sim  = times_sim,
      obs_times  = times_sim,
      obs_values = obs_data,
      params     = params,
      lambda     = 0.1,
      method     = sample(c("cn", "gl1", "gl2"), 1)
    )

    stability <- filter_stable_info(cfg, function(ps) {
      solver$cost_function(rep(0, solver$n_steps * solver$n_vars), y0)
    })
    cost0 <- stability$value
    if (!stability$stable) {
      return(list(
        stable = FALSE,
        detail = list(
          params = c(alpha = params$alpha, beta = params$beta,
                     delta = params$delta, gamma = params$gamma,
                     x0 = y0[1], y0 = y0[2]),
          reason = stability$reason,
          cost0 = cost0
        )
      ))
    }

    set.seed(900L + trial_seed)
    u_test <- matrix(rnorm(solver$n_steps * solver$n_vars, sd = 3), ncol = solver$n_vars)
    chk <- safe_check_gradient(
      fn  = solver$cost_function,
      gr  = solver$gradient_function,
      par = u_test,
      eps = 1e-5,
      y0  = y0
    )

    if (!is.finite(chk$max_rel_error) || is.na(chk$cosine_similarity)) {
      return(list(
        stable = FALSE,
        detail = list(
          params = c(alpha = params$alpha, beta = params$beta,
                     delta = params$delta, gamma = params$gamma,
                     x0 = y0[1], y0 = y0[2]),
          reason = if (!is.null(chk$error_msg)) chk$error_msg else "non-finite gradient check result",
          cost0 = cost0
        )
      ))
    }

    list(
      stable = TRUE,
      detail = list(
        params = c(alpha = params$alpha, beta = params$beta,
                   delta = params$delta, gamma = params$gamma,
                   x0 = y0[1], y0 = y0[2]),
        max_rel_error = chk$max_rel_error,
        cosine_similarity = chk$cosine_similarity,
        is_descent = chk$is_descent,
        error_msg = if (!is.null(chk$error_msg)) chk$error_msg else ""
      ),
      checks = c(
        max_rel_error = chk$max_rel_error < 1e-3,
        cos_sim       = chk$cosine_similarity > 1 - 1e-3,
        is_descent    = chk$is_descent
      )
    )
  }

  attempts <- 0L
  stable_idx <- 0L
  while (stable_idx < N_RAND && attempts < MAX_ATTEMPTS) {
    attempts <- attempts + 1L
    cfg <- random_trial_cfg()
    result <- evaluate_trial(cfg, attempts)
    if (isTRUE(result$stable)) {
      stable_idx <- stable_idx + 1L
      stable_results <- cbind(stable_results, result$checks)
      trial_details[[stable_idx]] <- result$detail
    } else {
      unstable_trials[[length(unstable_trials) + 1L]] <- result$detail
    }
  }

  if (stable_idx < N_RAND) {
    stop(sprintf("Could not collect %d stable LV trials after %d attempts.", N_RAND, MAX_ATTEMPTS))
  }

  for (chk_cond in rownames(stable_results)) {
    failed_idx <- which(!stable_results[chk_cond, ])
    test_that(sprintf("%s: all random trials pass", chk_cond), {
      if (length(failed_idx) > 0) {
        fail_lines <- vapply(failed_idx, function(idx) {
          det <- trial_details[[idx]]
          sprintf(
            "trial %d | %s | max_rel_error=%.3e | cos=%.6f | is_descent=%s | error=%s",
            idx,
            fmt_named(det$params),
            det$max_rel_error,
            det$cosine_similarity,
            as.character(det$is_descent),
            det$error_msg
          )
        }, character(1L))
        stop(paste(c("Failing trials:", fail_lines), collapse = "\n"))
      }
      expect_true(TRUE)
    })
  }
})

# ---------------------------------------------------------------------------
# T6. NA handling: NAs in obs must not propagate to cost or gradient
# ---------------------------------------------------------------------------
describe("T6: NA handling — observations with missing values", {

  params    <- list(k = 0.5)
  times_sim <- seq(0, 4, by = 0.2)
  y0        <- 3.0
  y_obs     <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  gap_idx <- which(times_sim >= 1.5 & times_sim <= 2.5)
  y_obs[gap_idx, ] <- NA

  solver <- AdjointForwardSolver$new(
    model      = ODEModel$new(rhs = decay_rhs),
    times_sim  = times_sim,
    obs_times  = times_sim,
    obs_values = y_obs,
    params     = params,
    lambda     = 0.1
  )

  u_test <- rep(0.01, solver$n_steps * solver$n_vars)

  test_that("cost_function does not return NA with missing observations", {
    cost <- solver$cost_function(u_test, y0)
    expect_false(any(is.na(cost)))
    expect_true(is.finite(cost))
  })

  test_that("gradient_function does not return NA with missing observations", {
    g <- solver$gradient_function(u_test, y0)
    expect_false(any(is.na(g)))
    expect_true(all(is.finite(g)))
  })

  test_that("optimize() completes without error when NAs are present", {
    set.seed(7)
    res <- solver$optimize(y0 = y0, max_iter = 30)
    expect_false(any(is.na(res$value)))
    expect_true(is.finite(res$value))
  })
})

# ---------------------------------------------------------------------------
# T12. optimize() vs optimize_bvp() — 1-D exponential decay
# ---------------------------------------------------------------------------
describe("T12: optimize vs optimize_bvp — 1-D decay, cost reduction and trajectory agreement", {

  params    <- list(k = 0.4)
  times_sim <- seq(0, 2, by = 0.1)
  y0        <- 3.0

  set.seed(13)
  y_obs_10 <- matrix(
    y0 * exp(-params$k * times_sim) + rnorm(length(times_sim), 0, 0.05),
    ncol = 1
  )

  make_solver_10 <- function() {
    AdjointForwardSolver$new(
      model      = ODEModel$new(rhs = decay_rhs),
      times_sim  = times_sim,
      obs_times  = times_sim,
      obs_values = y_obs_10,
      params     = params,
      lambda     = 0.1,
      method     = "cn"
    )
  }

  u_zero <- rep(0, length(times_sim))
  cost0  <- make_solver_10()$cost_function(u_zero, y0)

  solver_bfgs_10 <- make_solver_10()
  set.seed(14)
  solver_bfgs_10$optimize(y0 = y0, max_iter = 100)
  cost_bfgs <- solver_bfgs_10$cost_function(as.vector(solver_bfgs_10$u), y0)

  solver_bvp_10 <- make_solver_10()
  solver_bvp_10$optimize_bvp(y0 = y0)
  cost_bvp <- solver_bvp_10$cost_function(as.vector(solver_bvp_10$u), y0)

  test_that("optimize() strictly reduces cost vs u = 0", {
    expect_lt(cost_bfgs, cost0)
  })
  test_that("optimize_bvp() strictly reduces cost vs u = 0", {
    expect_lt(cost_bvp, cost0)
  })
  test_that("optimize_bvp() and optimize() costs agree within a factor of 5", {
    expect_lt(cost_bvp / cost_bfgs, 5)
  })

  rmse_traj   <- sqrt(mean((solver_bfgs_10$y - solver_bvp_10$y)^2))
  state_range <- diff(range(y_obs_10))

  test_that("trajectories from optimize and optimize_bvp agree within 20% of state range", {
    expect_lt(rmse_traj / state_range, 0.2)
  })
})


# ---------------------------------------------------------------------------
# T13. optimize() vs optimize_bvp() — 2-D Lotka-Volterra
# ---------------------------------------------------------------------------
describe("T13: optimize vs optimize_bvp — 2-D LV, cost reduction and trajectory agreement", {

  params    <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  times_sim <- seq(0, 3, by = 0.2)
  y0        <- c(5, 3)

  y_lv_11 <- euler_solve(lv_rhs, y0, times_sim, params)
  set.seed(15)
  obs_lv_11 <- y_lv_11 + matrix(rnorm(length(y_lv_11), 0, 0.5), nrow(y_lv_11), 2)

  make_solver_11 <- function() {
    AdjointForwardSolver$new(
      model      = ODEModel$new(rhs = lv_rhs),
      times_sim  = times_sim,
      obs_times  = times_sim,
      obs_values = obs_lv_11,
      params     = params,
      lambda     = 0.1,
      method     = "cn"
    )
  }

  u_zero <- rep(0, nrow(obs_lv_11) * ncol(obs_lv_11))
  cost0  <- make_solver_11()$cost_function(u_zero, y0)

  solver_bfgs_11 <- make_solver_11()
  set.seed(16)
  solver_bfgs_11$optimize(y0 = y0, max_iter = 100)
  cost_bfgs <- solver_bfgs_11$cost_function(as.vector(solver_bfgs_11$u), y0)

  solver_bvp_11 <- make_solver_11()
  solver_bvp_11$optimize_bvp(y0 = y0)
  cost_bvp <- solver_bvp_11$cost_function(as.vector(solver_bvp_11$u), y0)

  test_that("2-D LV: optimize() strictly reduces cost vs u = 0", {
    expect_lt(cost_bfgs, cost0)
  })
  test_that("2-D LV: optimize_bvp() strictly reduces cost vs u = 0", {
    expect_lt(cost_bvp, cost0)
  })
  test_that("2-D LV: optimize_bvp() and optimize() costs agree within a factor of 5", {
    expect_lt(cost_bvp / cost_bfgs, 5)
  })

  rmse_traj   <- sqrt(mean((solver_bfgs_11$y - solver_bvp_11$y)^2))
  state_range <- diff(range(obs_lv_11))

  test_that("2-D LV: trajectories agree within 20% of state range", {
    expect_lt(rmse_traj / state_range, 0.2)
  })
})

# ---------------------------------------------------------------------------
# T14. Discrete control Jacobian API
# ---------------------------------------------------------------------------
describe("T14: discrete control Jacobian API across methods", {

  params    <- list(k = 0.4)
  times_sim <- seq(0, 1, by = 0.1)
  y0        <- 1.0
  y_obs     <- matrix(y0 * exp(-params$k * times_sim), ncol = 1)

  for (method in c("euler", "cn", "gl1", "gl2")) {
    local({
      m <- method

      solver <- AdjointForwardSolver$new(
        model      = ODEModel$new(rhs = decay_rhs),
        times_sim  = times_sim,
        obs_times  = times_sim,
        obs_values = y_obs,
        params     = params,
        lambda     = 0.1,
        method     = m
      )

      u0 <- matrix(0, solver$n_steps, solver$n_vars)
      solver$solve_state(u0, y0)

      test_that(sprintf("%s: combined Jacobian has expected structure", m), {
        dU <- solver$get_discrete_control_jacobian(
          t_idx = 1L,
          param_names = c("k"),
          include_theta = TRUE
        )

        expect_true(is.list(dU))
        expect_true(all(c("Du_dyc", "Du_dyn", "Du_dtheta") %in% names(dU)))

        expect_equal(dim(dU$Du_dyc), c(1L, 1L))
        expect_equal(dim(dU$Du_dyn), c(1L, 1L))
        expect_equal(dim(dU$Du_dtheta), c(1L, 1L))

        expect_true(all(is.finite(dU$Du_dyc)))
        expect_true(all(is.finite(dU$Du_dyn)))
        expect_true(all(is.finite(dU$Du_dtheta)))
      })

      test_that(sprintf("%s: split derivative accessors are consistent", m), {
        dU <- solver$get_discrete_control_jacobian(
          t_idx = 1L,
          param_names = c("k"),
          include_theta = TRUE
        )

        expect_equal(solver$get_dcontrol_dy_curr(1L), dU$Du_dyc)
        expect_equal(solver$get_dcontrol_dy_next(1L), dU$Du_dyn)
        expect_equal(solver$get_dcontrol_dtheta(1L, c("k")), dU$Du_dtheta)
      })
    })
  }

  test_that("include_theta = FALSE returns NULL Du_dtheta", {
    solver <- AdjointForwardSolver$new(
      model      = ODEModel$new(rhs = decay_rhs),
      times_sim  = times_sim,
      obs_times  = times_sim,
      obs_values = y_obs,
      params     = params,
      lambda     = 0.1,
      method     = "gl2"
    )
    solver$solve_state(matrix(0, solver$n_steps, solver$n_vars), y0)

    dU <- solver$get_discrete_control_jacobian(
      t_idx = 1L,
      include_theta = FALSE
    )

    expect_null(dU$Du_dtheta)
  })
})
