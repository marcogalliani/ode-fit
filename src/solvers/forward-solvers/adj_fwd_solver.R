library(R6)

# =============================================================================
# AdjointForwardSolver
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
AdjointForwardSolver <- R6Class("AdjointForwardSolver",
  inherit = ForwardSolverBase,

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
    #' Initialize AdjointForwardSolver runtime state.
    initialize = function(model, times_sim, obs_times, obs_values, params, lambda,
                          method = "gl2", verbose = FALSE) {
      self$initialize_forward_solver(
        model = model,
        times_sim = times_sim,
        obs_times = obs_times,
        obs_values = obs_values,
        params = params,
        lambda = lambda,
        method = method,
        verbose = verbose
      )
      private$dto_solver <- make_dto_solver(method)
    },

    #' Partial sensitivity matrix S[ns, nv, np]:
    #'    S[t, v, j] = dy_{t,v} / d theta_j_physical  (u* held fixed)
    compute_sensitivity_matrix = function(params, param_names) {
      ns     <- self$n_steps
      nv     <- self$n_vars
      np     <- length(param_names)

      S <- array(0, c(ns, nv, np))
      J0 <- self$model$init_state_jacobian_fd(params, param_names)
      S[1L, , ] <- J0

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
      return(S)   # [ns, nv, np]
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
    compute_parameter_gradient_sens = function(param_names,
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

      S <- self$compute_sensitivity_matrix(self$params, param_names)

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

    #' Compute loss function gradient w.r.t. the ODE params through adjoint method
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

      grad_phys <- rep(0, np)
      names(grad_phys) <- param_names

      if (!is.null(init_state_jacobian)) {
        grad_phys <- grad_phys + as.numeric(crossprod(self$p[1L, ], init_state_jacobian))
      }

      # Adjoint Sweep (t = 1 to ns - 1)
      for (t in seq_len(ns - 1L)) {
        # Get the partial derivative of the discrete step equation w.r.t parameters
        dU <- self$get_discrete_control_jacobian(
          t_idx = t,
          param_names = param_names,
          include_theta = TRUE
        )
        # Envelope theorem term under the solver convention (+p^T U):
        # dH/dtheta = - sum_t dt[t] * p[t+1]^T (dU/dtheta)_t.
        grad_phys <- grad_phys - self$dt_vec[t] *
          as.numeric(crossprod(self$p[t + 1L, ], dU$Du_dtheta))
      }

      if (!return_normalized) {
        return(grad_phys)
      }

      if (is.null(scales) || length(scales) != np) {
        stop("scales must be provided with one entry per parameter when return_normalized = TRUE")
      }
      
      grad_norm <- grad_phys * as.numeric(scales)
      names(grad_norm) <- param_names
      return(grad_norm)
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
      bvp <- BvpForwardSolver$new(
        model = self$model,
        times_sim = self$times_sim,
        obs_times = self$times_sim,
        obs_values = self$observations_mapped,
        params = self$params,
        lambda = self$lambda,
        verbose = if (is.null(verbose)) self$verbose else verbose
      )

      sol <- bvp$optimize_bvp(y0 = y0, z_init = z_init, max_iter = max_iter, tol = tol)
      self$y <- bvp$y
      self$p <- bvp$p
      self$u <- bvp$u
      sol
    },

    optimize = function(y0, max_iter = 100, u_init = NULL,
                        reltol = sqrt(.Machine$double.eps), verbose = NULL) {
      if (is.null(verbose)) verbose <- self$verbose
      if (is.null(u_init)) u_init <- rep(0, self$n_steps * self$n_vars)

      y0 <- self$sanitize_initial_state(y0)

      self$log_debug("Starting optimization: method=%s max_iter=%d reltol=%.3e",
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
      self$log_debug("Optimization complete: converged_code=%d final_value=%.6e",
                        res$convergence, res$value)
      res
    }
  )
)
