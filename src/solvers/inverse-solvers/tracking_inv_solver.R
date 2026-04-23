library(R6)

# =============================================================================
# TrackingOdeSolver
# =============================================================================
TrackingOdeSolver <- R6Class("TrackingOdeSolver",
  inherit = InverseSolverBase,
  public = list(
    initialize = function(model, times_sim, obs_times, obs_values,
                          lambda,
                          inner_max_iter = 200,
                          inner_reltol   = sqrt(.Machine$double.eps),
              inner_method   = "gl2",
              verbose        = FALSE) {
      self$initialize_estimator(model, times_sim, obs_times, obs_values,
                                lambda, inner_max_iter, inner_reltol,
                inner_method, verbose)
    },

    # =========================================================================
    # 1. Outer Objective: H(theta) = J(y*(theta), u*(theta); theta)
    #    = (1/ns)*SSE + lambda*||u*||^2   — the FULL inner cost at the optimum
    # =========================================================================
    outer_objective = function(theta_norm, param_names) {
      # --- Cache hit ---
      if (!is.null(self$last_theta) && all(theta_norm == self$last_theta)) {
        solver <- self$last_solver
      } else {
        # --- Run inner optimisation ---
        p_phys <- self$unpack_physical(theta_norm, param_names)
        y0_phys <- self$eval_init_state(p_phys)
        solver <- self$inner_solver_class$new(
          model      = self$model,
          times_sim  = self$times_sim,
          obs_times  = self$obs_times,
          obs_values = self$obs_values,
          params     = p_phys,
          lambda     = self$lambda,
          method     = self$inner_method,
          verbose    = FALSE
        )
        solver$optimize(y0       = y0_phys,
                        u_init   = NULL,
                        max_iter = self$inner_max_iter,
                        reltol   = self$inner_reltol,
                        verbose  = FALSE)

        self$last_u      <- as.vector(solver$u)
        self$last_theta  <- theta_norm
        self$last_solver <- solver

        # Record outer iteration (use actual y0 after NA processing)
        y0_eff <- solver$y[1L, ]
        j_val  <- solver$cost_function(as.vector(solver$u), y0_eff)
        p_vals <- theta_norm * self$get_scales_vector(param_names)
        self$history[[length(self$history) + 1L]] <- list(
          iter   = length(self$history) + 1L,
          params = setNames(p_vals, param_names),
          cost   = j_val
        )
      }

      if (!all(is.finite(solver$y))) return(1e20)

      # Full inner cost J at the optimum (SSE/ns + lambda*||u*||^2)
      y0_eff <- solver$y[1L, ]
      j_val  <- solver$cost_function(as.vector(solver$u), y0_eff)

      p_vals <- theta_norm * self$get_scales_vector(param_names)
      self$log_iter(length(self$history), "J", j_val,
                    setNames(p_vals, param_names))
      return(j_val)
    },

    # =========================================================================
    # 2. Outer Gradient (Adjoint Method)
    #
    # By the envelope theorem, the outer gradient is equivalent to (***):
    #
    #   dH/dtheta_j = sum_{t=1}^{T-1} p[t+1]^T * f_theta_j(y*[t], t) * dt[t]
    #
    # where p is the inner adjoint already stored in solver$p after optimize().
    # This avoids an explicit forward sensitivity sweep.
    # =========================================================================
    outer_gradient = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s      <- self$last_solver
      if (is.null(s) || is.null(s$y) || is.null(s$p)) {
        return(rep(0, length(param_names)))
      }

      if (!all(is.finite(s$y)) || !all(is.finite(s$p))) {
        return(rep(0, length(param_names)))
      }

      p_phys <- self$unpack_physical(theta_norm, param_names)
      J0 <- self$model$init_state_jacobian_fd(p_phys, param_names)
      scales <- self$get_scales_vector(param_names)
      s$compute_parameter_gradient_adjoint(
        param_names = param_names,
        init_state_jacobian = J0,
        return_normalized = TRUE,
        scales = scales
      )
    },

    # Outer Gradient (Forward Sensitivity Method) — delegates to the inner solver
    #   dH/dtheta_j = (2/ns) * sum_t r[t]^T * S[t, , j]
    # Note that we just need the state sensitivity and not the control sens.
    # as dH/du = 0, as the inner solver already converged
    outer_gradient_sensitivity = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s      <- self$last_solver
      if (is.null(s) || is.null(s$y) || is.null(s$p)) {
        return(rep(0, length(param_names)))
      }

      if (!all(is.finite(s$y)) || !all(is.finite(s$p))) {
        return(rep(0, length(param_names)))
      }

      p_phys <- self$unpack_physical(theta_norm, param_names)
      J0 <- self$model$init_state_jacobian_fd(p_phys, param_names)
      scales <- self$get_scales_vector(param_names)
      s$compute_parameter_gradient_sens(
        param_names = param_names,
        init_state_jacobian = J0,
        return_normalized = TRUE,
        scales = scales
      )
    },
    outer_gradient_dispatch = function(theta_norm, param_names, gradient_mode = "adjoint") {
      mode <- match.arg(gradient_mode, c("adjoint", "sensitivity"))
      if (mode == "adjoint") {
        return(self$outer_gradient(theta_norm, param_names))
      }
      self$outer_gradient_sensitivity(theta_norm, param_names)
    },

    # =========================================================================
    # 4. Outer Parameter Optimisation  (L-BFGS-B)
    # =========================================================================
    optimize_parameters = function(init_theta_physical, param_names,
                                   lower_phys = NULL, upper_phys = NULL,
                                   gradient_mode = c("adjoint", "sensitivity")) {
      gradient_mode <- match.arg(gradient_mode)
      prep <- self$prepare_theta_normalized(param_names,
                                            init_theta_physical,
                                            lower_phys,
                                            upper_phys)
      scales <- prep$scales
      init_theta_norm <- prep$init
      lower_norm <- prep$lower
      upper_norm <- prep$upper

      self$history <- list()

      self$log_info("Starting constrained parameter tracking (L-BFGS-B)")
      if (isTRUE(self$verbose)) {
        self$log_info("Initial guess (normalized): %s",
                      paste(sprintf("%.6g", init_theta_norm), collapse = ", "))
      }

      old_warn <- getOption("warn")
      options(warn = 0)
      on.exit(options(warn = old_warn), add = TRUE)

      res <- optim(
        par         = init_theta_norm,
        fn          = function(par, param_names) self$outer_objective(par, param_names),
        gr          = function(par, param_names) {
          self$outer_gradient_dispatch(par, param_names, gradient_mode = gradient_mode)
        },
        param_names = param_names,
        method      = "L-BFGS-B",
        lower       = lower_norm,
        upper       = upper_norm,
        control     = list(maxit = 50,
                           trace = if (isTRUE(self$verbose)) 1 else 0)
      )

      final_params <- res$par * scales
      names(final_params) <- param_names
      self$log_info("Optimization complete: convergence=%d value=%.6e",
                    res$convergence, res$value)
      self$log_info("Estimated parameters: %s",
                    paste(sprintf("%s=%.6g", names(final_params), final_params), collapse = ", "))
      return(final_params)
    }
  )
)
