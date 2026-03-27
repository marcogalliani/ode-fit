library(R6)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)



# =============================================================================
# OdeSystemSolver
#
# Physics-informed smoother estimating unknown additive forcing u(t).
#
# Inner problem:
#   min_u  (1/ns)*SSE(y)  +  lambda * integral(u^2) dt
#   s.t.   dy/dt = f(y, t, theta) + u(t),   y(t_1) = y_0
#
# `method` selects the DtO scheme: "euler", "cn", "gl1", or "gl2".
# The corresponding DtOSolver (from ode_solvers.R) wraps a DtO scheme that
# provides both the forward integrator and its consistent discrete adjoint.
# =============================================================================
OdeSystemSolver <- R6Class("OdeSystemSolver",

  private = list(
    dto_solver         = NULL,
    cache_u            = NULL,
    cache_y            = NULL,
    cache_p            = NULL,
    cache_grad_contrib = NULL,

    # source_fn(t_idx) = (2/ns) * (y[t] - obs[t]),  NA obs -> 0
    make_source_fn = function(y_curr) {
      ns  <- self$n_steps
      obs <- self$observations_mapped
      function(t_idx) {
        r <- y_curr[t_idx, ] - obs[t_idx, ]
        r[is.na(r)] <- 0
        (2 / ns) * r
      }
    }
  ),

  public = list(
    func_rhs = NULL, params = NULL, lambda = NULL,
    times_sim = NULL, n_steps = NULL, dt_vec = NULL, n_vars = NULL,
    method = NULL,

    observations_mapped = NULL,
    y = NULL, u = NULL, p = NULL,

    initialize = function(func_rhs, times_sim, obs_times, obs_values, params, lambda, method = "gl2") {
      self$func_rhs <- func_rhs
      obs_times  <- round(obs_times,  digits = 10)
      times_sim  <- sort(unique(round(c(times_sim, obs_times), digits = 10)))
      self$times_sim <- times_sim
      self$params    <- params
      self$lambda    <- lambda
      self$n_steps   <- length(times_sim)
      self$n_vars    <- ncol(obs_values)
      self$dt_vec    <- c(diff(times_sim), 0)
      self$method    <- method
      private$dto_solver <- make_dto_solver(method)

      self$observations_mapped <- matrix(NA, nrow = self$n_steps, ncol = self$n_vars)
      self$observations_mapped[times_sim %in% obs_times, ] <- obs_values
      self$y <- matrix(0, self$n_steps, self$n_vars)
      self$u <- matrix(0, self$n_steps, self$n_vars)
      self$p <- matrix(0, self$n_steps, self$n_vars)
    },

    get_jacobian = function(y_vec, t_val) {
      n   <- length(y_vec)
      J   <- matrix(0, n, n)
      eps <- 1e-7
      for (j in seq_len(n)) {
        y_p <- y_vec; y_p[j] <- y_p[j] + eps
        y_m <- y_vec; y_m[j] <- y_m[j] - eps
        J[, j] <- (self$func_rhs(y_p, t_val, self$params) -
                   self$func_rhs(y_m, t_val, self$params)) / (2 * eps)
      }
      J
    },

    solve_state = function(u_mat, y0) {
      u_fn   <- function(t) u_mat[pmax(1L, pmin(findInterval(t, self$times_sim), self$n_steps)), ]
      jac_fn <- function(y, t) self$get_jacobian(y, t)
      rhs    <- function(y, t) self$func_rhs(y, t, self$params) + u_fn(t)
      private$dto_solver$solve_state(rhs, y0, self$times_sim, jac_fn)
    },

    solve_adjoint = function(y_fwd) {
      rhs_cont  <- function(y, t) self$func_rhs(y, t, self$params)
      jac_fn    <- function(y, t) self$get_jacobian(y, t)
      source_fn <- private$make_source_fn(y_fwd)
      pT        <- rep(0, self$n_vars)
      private$dto_solver$solve_adjoint(rhs_cont, pT, jac_fn, source_fn)
    },

    solve_state_adjoint = function(u_mat, y0) {
      fwd <- self$solve_state(u_mat, y0)
      if (!all(is.finite(fwd$y)))
        return(list(y = fwd$y, p = NULL, grad_contrib = NULL, converged = FALSE))
      bwd <- self$solve_adjoint(fwd$y)
      list(y = fwd$y, p = bwd$p, grad_contrib = bwd$grad_contrib, converged = TRUE)
    },

    cost_function = function(u_flat, y0) {
      ns    <- self$n_steps; nv <- self$n_vars
      u_mat <- matrix(u_flat, ns, nv)

      sol <- self$solve_state_adjoint(u_mat, y0)
      if (!sol$converged || !all(is.finite(sol$y))) return(1e20)

      private$cache_u            <- u_flat
      private$cache_y            <- sol$y
      private$cache_p            <- sol$p
      private$cache_grad_contrib <- sol$grad_contrib

      sse    <- sum((self$observations_mapped - sol$y)^2, na.rm = TRUE)
      dt_l   <- c(0, self$dt_vec[seq_len(ns - 1L)])
      w_trap <- (dt_l + self$dt_vec) / 2
      reg    <- self$lambda * sum(w_trap * rowSums(u_mat^2))
      (1 / ns) * sse + reg
    },

    gradient_function = function(u_flat, y0) {
      ns <- self$n_steps; nv <- self$n_vars

      cache_valid <- !is.null(private$cache_u) &&
                     length(private$cache_u) == length(u_flat) &&
                     isTRUE(all.equal(private$cache_u, u_flat, tolerance = 0))

      if (!cache_valid) {
        u_mat <- matrix(u_flat, ns, nv)
        sol   <- self$solve_state_adjoint(u_mat, y0)
        private$cache_u            <- u_flat
        private$cache_y            <- sol$y
        private$cache_p            <- sol$p
        private$cache_grad_contrib <- sol$grad_contrib
      }

      u_mat  <- matrix(u_flat, ns, nv)
      dt_l   <- c(0, self$dt_vec[seq_len(ns - 1L)])
      w_trap <- matrix((dt_l + self$dt_vec) / 2, ns, nv)
      as.vector(2 * self$lambda * w_trap * u_mat + private$cache_grad_contrib)
    },

    optimize_bvp = function(y0, z_init = NULL,
                            max_iter = 50L, tol = 1e-8, verbose = FALSE) {
      ns <- self$n_steps; ny <- self$n_vars
      jac_fn <- function(y, t) self$get_jacobian(y, t)

      if (length(y0) == 1L && is.na(y0)) {
        first_row <- which(!is.na(self$observations_mapped[, 1L]))[1L]
        y0 <- self$observations_mapped[first_row, ]; y0[is.na(y0)] <- 0
      } else if (any(is.na(y0))) {
        for (v in seq_len(ny)) {
          if (is.na(y0[v])) {
            first_v <- which(!is.na(self$observations_mapped[, v]))[1L]
            y0[v] <- if (!is.na(first_v)) self$observations_mapped[first_v, v] else 0
          }
        }
      }

      obs_clean   <- self$observations_mapped; obs_clean[is.na(obs_clean)] <- 0
      obs_mask    <- matrix(as.numeric(!is.na(self$observations_mapped)), ns, ny)
      two_lam     <- 2 * self$lambda
      two_over_ns <- 2 / ns

      t_keys  <- as.character(round(self$times_sim, 10))
      t_lut   <- setNames(seq_len(ns), t_keys)
      get_idx <- function(t_val) t_lut[[as.character(round(t_val, 10))]]

      F_rhs_inner <- function(t_val, z) {
        y_t <- z[seq_len(ny)]; p_t <- z[ny + seq_len(ny)]
        ti  <- get_idx(t_val)
        fy  <- self$func_rhs(y_t, t_val, self$params)
        Jfy <- self$get_jacobian(y_t, t_val)
        c(fy - p_t / two_lam,
          as.vector(-(t(Jfy) %*% p_t)) -
            two_over_ns * obs_mask[ti, ] * (y_t - obs_clean[ti, ]))
      }

      obs_T  <- obs_clean[ns, ]; mask_T <- obs_mask[ns, ]
      bc_inner <- function(z_l, z_r) {
        c(z_l[seq_len(ny)] - y0,
          z_r[ny + seq_len(ny)] - two_over_ns * mask_T * (z_r[seq_len(ny)] - obs_T))
      }

      if (is.null(z_init)) {
        euler_dto <- make_dto_solver("euler")
        rhs_euler <- function(y, t) self$func_rhs(y, t, self$params)
        y_guess   <- euler_dto$solve_state(rhs_euler, y0, self$times_sim)$y
        source_fn <- private$make_source_fn(y_guess)
        pT        <- rep(0, ny)
        p_guess   <- euler_dto$solve_adjoint(rhs_euler, pT, jac_fn, source_fn)$p
        z_init    <- cbind(y_guess, p_guess)
      }

      sol <- solve_bvp_colloc(
        F_rhs       = F_rhs_inner,
        bc_residual = bc_inner,
        t_grid      = self$times_sim,
        z_init      = z_init,
        max_iter    = max_iter,
        tol         = tol,
        verbose     = verbose
      )

      y_sol <- sol$z[, seq_len(ny), drop = FALSE]
      p_sol <- sol$z[, ny + seq_len(ny), drop = FALSE]
      self$y <- y_sol; self$p <- p_sol; self$u <- -p_sol / two_lam
      sol
    },

    optimize = function(y0, max_iter = 100, u_init = NULL,
                        reltol = sqrt(.Machine$double.eps)) {
      if (is.null(u_init)) u_init <- rep(0, self$n_steps * self$n_vars)

      if (length(y0) == 1 && is.na(y0)) {
        first_row <- which(!is.na(self$observations_mapped[, 1]))[1]
        y0 <- self$observations_mapped[first_row, ]; y0[is.na(y0)] <- 0
      } else if (any(is.na(y0))) {
        for (v in seq_len(self$n_vars)) {
          if (is.na(y0[v])) {
            first_v <- which(!is.na(self$observations_mapped[, v]))[1]
            y0[v] <- if (!is.na(first_v)) self$observations_mapped[first_v, v] else 0
          }
        }
      }

      res <- optim(par = u_init, fn = self$cost_function, gr = self$gradient_function,
                   y0 = y0, method = "BFGS",
                   control = list(maxit = max_iter, reltol = reltol, trace = 1))

      self$u <- matrix(res$par, self$n_steps, self$n_vars)
      if (!is.null(private$cache_y)) {
        self$y <- private$cache_y; self$p <- private$cache_p
      } else {
        sol <- self$solve_state_adjoint(self$u, y0)
        self$y <- sol$y; self$p <- sol$p
      }
      res
    }
  )
)
