library(R6)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)

# =============================================================================
# DtOForwardSolver
#
# Physics-informed smoother estimating unknown additive forcing u(t).
#
# Inner problem:
#   min_u  (1/ns)*SSE(y)  +  lambda * integral(u^2) dt
#   s.t.   dy/dt = f(y, t, theta) + u(t),   y(t_1) = y_0
#
# `method` selects the DtO scheme: "euler", "cn", "gl1", or "gl2".
# The corresponding OCPDtOSolver (from ode_solvers.R) wraps a DtO scheme that
# provides both the forward integrator and its consistent discrete adjoint.
#
# =============================================================================
DtOForwardSolver <- R6Class("DtOForwardSolver",

  private = list(
    dto_solver         = NULL,
    cache_u            = NULL,
    cache_y            = NULL,
    cache_p            = NULL,
    cache_grad_contrib = NULL,

    log_debug = function(fmt, ...) {
      if (!isTRUE(self$verbose)) return(invisible(NULL))
      message(sprintf(paste0("[DtOForwardSolver] ", fmt), ...))
      invisible(NULL)
    },

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
    model = NULL, params = NULL, lambda = NULL,
    times_sim = NULL, n_steps = NULL, dt_vec = NULL, n_vars = NULL,
    method = NULL,
    verbose = FALSE,

    observations_mapped = NULL,
    y = NULL, u = NULL, p = NULL,

    #' Initialize DtOForwardSolver runtime state.
    initialize = function(model, times_sim, obs_times, obs_values, params, lambda,
                          method = "gl2", verbose = FALSE) {
      if (is.null(model) || !inherits(model, "ODEModel")) {
        stop("model must be an ODEModel instance")
      }

      self$model <- model
      obs_times  <- round(obs_times,  digits = 10)
      times_sim  <- sort(unique(round(c(times_sim, obs_times), digits = 10)))
      self$times_sim <- times_sim
      self$params    <- params
      self$lambda    <- lambda
      self$n_steps   <- length(times_sim)
      self$n_vars    <- ncol(obs_values)
      self$dt_vec    <- c(diff(times_sim), 0)
      self$method    <- method
      self$verbose   <- isTRUE(verbose)
      private$dto_solver <- make_dto_solver(method)

      self$observations_mapped <- matrix(NA, nrow = self$n_steps, ncol = self$n_vars)
      self$observations_mapped[times_sim %in% obs_times, ] <- obs_values
      self$y <- matrix(0, self$n_steps, self$n_vars)
      self$u <- matrix(0, self$n_steps, self$n_vars)
      self$p <- matrix(0, self$n_steps, self$n_vars)
    },
    

    #' Compute loss function gradient w.r.t. ODE params from discrete sensitivities.
    #'
    #' @description
    #' Computes `dH/dtheta` for the outer parameter problem by propagating
    #' state sensitivities through differentiated discrete control maps and
    #' contracting them with residuals.
    #'
    #' @details Builds sensitivity tensor `S[time, state, param]` via one-step
    #'   recursion using discrete Jacobian blocks from
    #'   `get_discrete_control_jacobian()`.
    #' 
    #' @note If `self$y` or `self$p` contains non-finite values, returns a zero
    #'   gradient vector instead of stopping.
    compute_parameter_gradient_adjoint = function(param_names,
                                                  init_state_jacobian = NULL,
                                                  return_normalized = FALSE,
                                                  scales = NULL) {
      if (is.null(self$y) || is.null(self$p)) {
        stop("State/adjoint trajectories are not available. Run optimize() first.")
      }

      if (!all(is.finite(self$y)) || !all(is.finite(self$p))) {
        grad <- rep(0, length(param_names))
        names(grad) <- param_names
        return(grad)
      }

      ns <- self$n_steps
      np <- length(param_names)
      nv <- self$n_vars

      S <- array(0, c(ns, nv, np))
      if (!is.null(init_state_jacobian)) {
        S[1L, , ] <- init_state_jacobian
      }

      for (t in seq_len(ns - 1L)) {
        dU <- self$get_discrete_control_jacobian(
          t_idx = t,
          param_names = param_names,
          include_theta = TRUE
        )

        rhs <- dU$Du_dyc %*% S[t, , ] + dU$Du_dtheta
        S[t + 1L, , ] <- tryCatch(
          -solve(dU$Du_dyn, rhs),
          error = function(e) -qr.solve(dU$Du_dyn, rhs)
        )
      }

      resid <- ifelse(is.na(self$observations_mapped), 0,
                      self$y - self$observations_mapped)
      grad_phys <- vapply(
        seq_len(np),
        function(j) sum((2 / ns) * resid * S[, , j]),
        numeric(1L)
      )

      names(grad_phys) <- param_names
      if (!return_normalized) {
        return(grad_phys)
      }

      if (is.null(scales) || length(scales) != np) {
        stop("scales must be provided with one entry per parameter when return_normalized = TRUE")
      }
      grad_norm <- grad_phys * as.numeric(scales)
      names(grad_norm) <- param_names
      grad_norm
    },

    #' Differentiate one-step discrete control map.
    #'
    #' Returns Jacobian blocks of the local discrete control equation with
    #' respect to `y_t`, `y_{t+1}`, and optionally `theta`.
    #'
    #' @details Delegates to
    #'   `private$dto_solver$differentiate_discrete_control_map()` using
    #'   internally constructed Jacobian callbacks.
    get_discrete_control_jacobian = function(t_idx, param_names = NULL,
                                             eps = 1e-7, include_theta = TRUE) {
      if (is.null(self$y)) {
        stop("State trajectory is not available. Run optimize() first.")
      }
      if (is.null(param_names)) {
        param_names <- names(self$params)
      }

      jac_fn <- function(y, t) self$model$jacobian_state(y, t, self$params)
      param_jac_fn <- NULL
      if (include_theta) {
        param_jac_fn <- function(y, t) {
          self$model$jacobian_param(
            y = y, t = t, params = self$params,
            param_names = param_names,
            eps = eps
          )
        }
      }
      private$dto_solver$differentiate_discrete_control_map(t_idx, jac_fn, param_jac_fn)
    },

    # Convenience accessor for `dU/dy_curr`.
    get_dcontrol_dy_curr = function(t_idx, eps = 1e-7) {
      self$get_discrete_control_jacobian(
        t_idx = t_idx,
        eps = eps,
        include_theta = FALSE
      )$Du_dyc
    },

    # Convenience accessor for `dU/dy_next`.
    get_dcontrol_dy_next = function(t_idx, eps = 1e-7) {
      self$get_discrete_control_jacobian(
        t_idx = t_idx,
        eps = eps,
        include_theta = FALSE
      )$Du_dyn
    },

    # Convenience accessor for `dU/dtheta`.
    get_dcontrol_dtheta = function(t_idx, param_names = NULL, eps = 1e-7) {
      self$get_discrete_control_jacobian(
        t_idx = t_idx,
        param_names = param_names,
        eps = eps,
        include_theta = TRUE
      )$Du_dtheta
    },

    # Solve forward state trajectory for fixed controls.
    solve_state = function(u_mat, y0) {
      u_fn   <- function(t) u_mat[pmax(1L, pmin(findInterval(t, self$times_sim), self$n_steps)), ]
      jac_fn <- function(y, t) self$model$jacobian_state(y, t, self$params)
      rhs    <- function(y, t) self$model$rhs(y, t, self$params) + u_fn(t)
      private$dto_solver$solve_state(rhs, y0, self$times_sim, jac_fn)
    },

    # Solve discrete adjoint trajectory for a forward state path.
    solve_adjoint = function(y_fwd) {
      rhs_cont  <- function(y, t) self$model$rhs(y, t, self$params)
      jac_fn    <- function(y, t) self$model$jacobian_state(y, t, self$params)
      source_fn <- private$make_source_fn(y_fwd)
      pT        <- rep(0, self$n_vars)
      private$dto_solver$solve_adjoint(rhs_cont, pT, jac_fn, source_fn)
    },

    # Solve state and adjoint in one pipeline call.
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
      if (!sol$converged || !all(is.finite(sol$y))) {
        private$cache_u            <- u_flat
        private$cache_y            <- sol$y
        private$cache_p            <- sol$p
        private$cache_grad_contrib <- rep(0, length(u_flat))
        return(1e20)
      }

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
      grad_len <- length(u_flat)

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

      if (is.null(private$cache_grad_contrib) || length(private$cache_grad_contrib) != grad_len) {
        private$cache_grad_contrib <- rep(0, grad_len)
      }

      u_mat  <- matrix(u_flat, ns, nv)
      dt_l   <- c(0, self$dt_vec[seq_len(ns - 1L)])
      w_trap <- matrix((dt_l + self$dt_vec) / 2, ns, nv)
      as.vector(2 * self$lambda * w_trap * u_mat + private$cache_grad_contrib)
    },

    optimize_bvp = function(y0, z_init = NULL,
                            max_iter = 50L, tol = 1e-8, verbose = NULL) {
      ns <- self$n_steps; ny <- self$n_vars
      if (is.null(verbose)) verbose <- self$verbose
      jac_fn <- function(y, t) self$model$jacobian_state(y, t, self$params)

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
        fy  <- self$model$rhs(y_t, t_val, self$params)
        Jfy <- self$model$jacobian_state(y_t, t_val, self$params)
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
        rhs_euler <- function(y, t) self$model$rhs(y, t, self$params)
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
                        reltol = sqrt(.Machine$double.eps), verbose = NULL) {
      if (is.null(verbose)) verbose <- self$verbose
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

      private$log_debug("Starting optimization: method=%s max_iter=%d reltol=%.3e",
            self$method, max_iter, reltol)
      res <- optim(par = u_init, fn = self$cost_function, gr = self$gradient_function,
                   y0 = y0, method = "BFGS",
           control = list(maxit = max_iter,
                  reltol = reltol,
                  trace = if (isTRUE(verbose)) 1 else 0))

      self$u <- matrix(res$par, self$n_steps, self$n_vars)
      if (!is.null(private$cache_y)) {
        self$y <- private$cache_y; self$p <- private$cache_p
      } else {
        sol <- self$solve_state_adjoint(self$u, y0)
        self$y <- sol$y; self$p <- sol$p
      }
      private$log_debug("Optimization complete: converged_code=%d final_value=%.6e",
                        res$convergence, res$value)
      res
    }
  )
)
