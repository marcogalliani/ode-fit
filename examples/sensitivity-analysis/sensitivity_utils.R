# =============================================================================
# examples/sensitivity-analysis/sensitivity_utils.R
#
# Shared utilities for local sensitivity analysis.  Sourced by:
#   examples/sensitivity-analysis/local_sensitivity_analysis.R
#   examples/non-linear-odes-solver/ode_systems_examples.R
#   examples/param-cascading/ode_params_examples.R
#
# Requires: deSolve, MASS, examples/ode_models.R (for make_sens_config).
# =============================================================================

library(deSolve)
library(MASS)


# ODE solver wrapper ---------------------------------------------------------
# Solves the ODE for a (possibly perturbed) parameter set, substituting
# `pars` into `base_pars` so that inactive parameters remain unchanged.
# Returns a data frame with columns  time, <state_vars...>  as returned by
# deSolve::ode().
ode_solver <- function(pars, model, init_state, t_grid,
                       rhs_args = list(), base_pars = NULL) {
  full_pars <- if (is.null(base_pars)) {
    pars
  } else {
    bp <- base_pars
    bp[names(pars)] <- pars
    bp
  }
  ode_args <- c(
    list(y = init_state, times = t_grid, func = model, parms = full_pars),
    rhs_args
  )
  as.data.frame(do.call(ode, ode_args))
}


# Central-difference Jacobian ------------------------------------------------
# Returns a matrix J of shape  (n_times * n_states)  x  n_active_params.
# Column j of J is  (traj(p + h/2*e_j) - traj(p - h/2*e_j)) / h  where
# traj() stacks all state columns into a single vector (row-major over time).
# Perturbation step: h = |p_j| * eps + eps  (relative + absolute).
compute_jacobian <- function(model, pars, init_state, t_grid,
                             active_params = names(pars),
                             rhs_args = list(), eps = 1e-5) {
  build_args <- function(p)
    c(list(y = init_state, times = t_grid, func = model, parms = p), rhs_args)

  n_obs  <- length(as.matrix(do.call(ode, build_args(pars)))[, -1])
  n_pars <- length(active_params)

  J <- matrix(0, nrow = n_obs, ncol = n_pars,
               dimnames = list(NULL, active_params))

  for (i in seq_along(active_params)) {
    nm <- active_params[i]
    h  <- abs(pars[nm]) * eps + eps

    p_p <- pars; p_p[nm] <- pars[nm] + h / 2
    p_m <- pars; p_m[nm] <- pars[nm] - h / 2

    v_p <- as.matrix(do.call(ode, build_args(p_p)))[, -1, drop = FALSE]
    v_m <- as.matrix(do.call(ode, build_args(p_m)))[, -1, drop = FALSE]
    J[, i] <- (v_p - v_m) / h
  }
  J
}


# Full local sensitivity analysis --------------------------------------------
# cfg    — output of make_sens_config() (from examples/ode_models.R)
# verbose — print FIM diagnostics and correlation matrix
#
# What is assessed:
#   J        — Jacobian of the trajectory w.r.t. active_params at cfg$params
#   FIM      — J^T J  (Fisher Information Matrix, unit noise assumed)
#   rank/cond — identifiability: rank-deficient or ill-conditioned FIM means
#               the data cannot uniquely determine all active parameters
#   COR      — parameter correlation matrix derived from (J^T J)^{-1}
#              High off-diagonal entries indicate co-varying parameters
#   sens_rms — RMS sensitivity across states at each time step; shows when
#              each parameter most influences the trajectory
#
# Returns invisibly: list(J, FIM, COR, sens_rms)
run_local_sensitivity <- function(cfg, verbose = TRUE) {
  J <- compute_jacobian(
    model         = cfg$model,
    pars          = cfg$params,
    init_state    = cfg$state,
    t_grid        = cfg$t_grid,
    active_params = cfg$active_params,
    rhs_args      = cfg$rhs_args
  )

  FIM <- t(J) %*% J
  sv  <- svd(FIM)$d
  cond_fim <- max(sv) / max(min(sv), .Machine$double.eps)
  rank_fim <- qr(FIM)$rank

  COV <- tryCatch(
    solve(FIM),
    error = function(e) {
      message("  [sensitivity] FIM is singular; using Moore-Penrose pseudo-inverse.")
      ginv(FIM + 1e-15 * diag(nrow(FIM)))
    }
  )
  COR <- cov2cor(COV)

  if (verbose) {
    cat(sprintf(
      "  FIM rank: %d / %d   |   condition number: %s\n",
      rank_fim, ncol(FIM), format(cond_fim, scientific = TRUE, digits = 3)
    ))
    cat("  Singular values:", paste(signif(sv, 3), collapse = ", "), "\n")
    cat("  Parameter correlation at true parameters:\n")
    print(round(COR, 4))
  }

  # Time-resolved RMS sensitivity
  n_times  <- length(cfg$t_grid)
  n_states <- length(cfg$state)
  n_pars   <- length(cfg$active_params)

  sens_rms <- matrix(0, n_times, n_pars,
                     dimnames = list(NULL, cfg$active_params))
  for (i in seq_len(n_pars)) {
    J_par        <- matrix(J[, i], nrow = n_times, ncol = n_states)
    sens_rms[, i] <- sqrt(rowMeans(J_par^2))
  }

  par_colors <- seq_len(n_pars)
  matplot(cfg$t_grid, sens_rms,
          type = "l", lty = 1, col = par_colors, lwd = 1.5,
          xlab = "Time", ylab = "RMS sensitivity",
          main = "Local sensitivity at true parameters")
  legend("topright", legend = cfg$active_params,
         col = par_colors, lty = 1, bty = "n")

  invisible(list(J = J, FIM = FIM, COR = COR, sens_rms = sens_rms))
}
