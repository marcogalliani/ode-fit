source("src/solvers/forward-solvers/load_forward_solvers.R")
AdjointForwardSolver <- get("AdjointForwardSolver", envir = .GlobalEnv)

lv_rhs <- function(y, t, p) c(
  p$alpha * y[1] - p$beta * y[1] * y[2],
  p$delta * y[1] * y[2] - p$gamma * y[2]
)

euler_solve <- function(rhs, y0, times, params) {
  dt_v <- c(diff(times), 0)
  out <- matrix(0, length(times), length(y0))
  out[1, ] <- y0
  for (i in seq_len(length(times) - 1)) {
    out[i + 1, ] <- out[i, ] + dt_v[i] * rhs(out[i, ], times[i], params)
  }
  out
}

fmt_named <- function(x, digits = 4) {
  paste(sprintf("%s=%.*f", names(x), digits, as.numeric(x)), collapse = ", ")
}

main <- function() {
  params <- list(alpha = 2.969, beta = 1.314, delta = 0.524, gamma = 0.3953)
  y0 <- c(6.387, 16.26)
  times_ref <- seq(0, 40, by = 0.01)
  times_sim <- seq(0, 40, by = 0.2)

  y_ref <- euler_solve(lv_rhs, y0, times_ref, params)
  obs_idx <- seq(1L, length(times_ref), by = 20L)
  set.seed(914L)
  obs_data <- y_ref[obs_idx, , drop = FALSE] +
    matrix(rnorm(length(times_sim) * 2, 0, 0.2), length(times_sim), 2)

  cat("Deterministic DtO overflow example\n")
  cat(sprintf("params: %s\n", fmt_named(c(alpha = params$alpha, beta = params$beta,
                                           delta = params$delta, gamma = params$gamma,
                                           x0 = y0[1], y0 = y0[2]))))
  cat(sprintf("Reference trajectory finite: %s\n", as.character(all(is.finite(y_ref)))))
  cat(sprintf("Reference max |y|: %.6f\n\n", max(abs(y_ref))))

  for (method in c("cn", "gl1", "gl2")) {
    result <- tryCatch({
      solver <- AdjointForwardSolver$new(
        func_rhs   = lv_rhs,
        times_sim  = times_sim,
        obs_times  = times_sim,
        obs_values = obs_data,
        params     = params,
        lambda     = 1e2,
        method     = method
      )

      u0 <- rep(0, solver$n_steps * solver$n_vars)
      sol <- solver$solve_state_adjoint(matrix(u0, solver$n_steps, solver$n_vars), y0)
      cost0 <- solver$cost_function(u0, y0)

      sse0 <- sum((solver$observations_mapped - sol$y)^2, na.rm = TRUE)
      reg0 <- solver$lambda * sum(matrix(solver$dt_vec, nrow = solver$n_steps, ncol = solver$n_vars) * matrix(u0^2, ncol = solver$n_vars))

      sprintf("method=%s | state_finite=%s | max|state|=%.6f | sse=%s | reg=%s | cost(u=0)=%s",
              method,
              as.character(all(is.finite(sol$y))),
              max(abs(sol$y), na.rm = TRUE),
              as.character(sse0),
              as.character(reg0),
              as.character(cost0))
    }, error = function(e) {
      sprintf("method=%s | ERROR: %s", method, conditionMessage(e))
    })

    cat(result, "\n")
  }

  cat("\nInterpretation:\n")
  cat("- The dense reference Lotka-Volterra rollout stays finite.\n")
  cat("- The DtO discretization can still hit a singular linear solve or overflow in the objective path on the same draw.\n")
  cat("- This is a numerical issue in the discrete solver/objective, not an ill-defined continuous ODE.\n")
}

if (sys.nframe() == 0L) {
  main()
}
