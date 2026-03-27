library(R6)
library(ggplot2)
library(gridExtra)

source("src/solvers/general_ode_system_solver.R")
source("src/solvers/parameter_estimator_base.R")

# =============================================================================
# TrackingOdeSolver
# =============================================================================
TrackingOdeSolver <- R6Class("TrackingOdeSolver",
  inherit = ParameterEstimatorBase,
  public = list(
    initialize = function(func_rhs, times_sim, obs_times, obs_values,
                fixed_params, lambda, param_scales,
                init_state,
                          inner_max_iter = 200,
                          inner_reltol   = sqrt(.Machine$double.eps)) {
      self$initialize_estimator(func_rhs, times_sim, obs_times, obs_values,
                                fixed_params, lambda, param_scales,
                                init_state, inner_max_iter, inner_reltol)
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
          func_rhs   = self$func_rhs,
          times_sim  = self$times_sim,
          obs_times  = self$obs_times,
          obs_values = self$obs_values,
          params     = p_phys,
          lambda     = self$lambda
        )
        solver$optimize(y0       = y0_phys,
                        u_init   = NULL,
                        max_iter = self$inner_max_iter,
                        reltol   = self$inner_reltol)

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
      cat(sprintf("Iter | Params: %s | J: %.4f\n",
                  paste(round(p_vals, 2), collapse = ","), j_val))
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
      ns     <- s$n_steps
      np     <- length(param_names)
      grad_phys <- numeric(np)
      p_phys <- self$unpack_physical(theta_norm, param_names)
      J0 <- self$init_state_jacobian_fd(p_phys, param_names)
      grad_phys <- grad_phys + as.vector(t(J0) %*% s$p[1L, ])

      for (t in seq_len(ns - 1L)) {
        dt    <- s$dt_vec[t]
        Fth   <- self$get_param_jacobian(s$y[t, ], s$times_sim[t],
                                         p_phys, param_names)
        p_next <- s$p[t + 1L, ]
        grad_phys <- grad_phys + as.vector(t(Fth) %*% p_next) * dt
      }

      scales <- self$get_scales_vector(param_names)
      return(grad_phys * scales)
    },

    outer_gradient_dispatch = function(theta_norm, param_names, gradient_mode = "adjoint") {
      mode <- match.arg(gradient_mode, c("adjoint", "sensitivity"))
      if (mode == "adjoint") {
        return(self$outer_gradient(theta_norm, param_names))
      }
      self$outer_gradient_sensitivity(theta_norm, param_names)
    },

    # =========================================================================
    # 3. Partial sensitivity matrix S[ns, nv, np]:
    #    S[t, v, j] = dy_{t,v} / d theta_j_physical  (u* held fixed)
    # =========================================================================
    compute_sensitivity_matrix = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s      <- self$last_solver
      p_phys <- self$unpack_physical(theta_norm, param_names)
      ns     <- s$n_steps
      nv     <- s$n_vars
      np     <- length(param_names)

      S <- array(0, c(ns, nv, np))
      J0 <- self$init_state_jacobian_fd(p_phys, param_names)
      S[1L, , ] <- J0

      for (t in seq_len(ns - 1L)) {
        dt    <- s$dt_vec[t]
        Fth   <- self$get_param_jacobian(s$y[t, ], s$times_sim[t],
                                         p_phys, param_names)
        Fy  <- s$get_jacobian(s$y[t, ], s$times_sim[t])
        At  <- diag(nv) + dt * Fy
        for (j in seq_len(np)) {
          S[t + 1L, , j] <- At %*% S[t, , j] + dt * Fth[, j]
        }
      }
      return(S)   # [ns, nv, np]
    },

    # Outer Gradient (Forward Sensitivity Method) — delegates to
    # compute_sensitivity_matrix and applies the chain rule.
    #   dH/dtheta_j = (2/ns) * sum_t r[t]^T * S[t, , j]
    outer_gradient_sensitivity = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s     <- self$last_solver
      ns    <- s$n_steps
      np    <- length(param_names)
      resid <- ifelse(is.na(s$observations_mapped), 0,
                      s$y - s$observations_mapped)

      S <- self$compute_sensitivity_matrix(theta_norm, param_names)

      grad_phys <- vapply(seq_len(np),
                          function(j) sum((2 / ns) * resid * S[, , j]),
                          numeric(1L))

      scales <- self$get_scales_vector(param_names)
      return(grad_phys * scales)
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

      cat("=== Starting Constrained Parameter Tracking (L-BFGS-B) ===\n")
      cat("Initial Guess (Norm):", round(init_theta_norm, 4), "\n")

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
        control     = list(maxit = 30, trace = 1)
      )

      final_params <- res$par * scales
      names(final_params) <- param_names
      cat("\n=== Optimization Complete ===\n")
      print(final_params)
      return(final_params)
    }
  )
)
