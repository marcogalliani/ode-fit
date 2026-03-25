library(R6)
library(ggplot2)
library(gridExtra)

source("src/solvers/general_ode_system_solver.R")

# =============================================================================
# TrackingOdeSolver
# =============================================================================
TrackingOdeSolver <- R6Class("TrackingOdeSolver",
  public = list(
    inner_solver_class = NULL,
    times_sim = NULL, obs_times = NULL, obs_values = NULL, y0 = NULL,
    func_rhs = NULL, fixed_params = NULL, lambda = NULL,
    param_scales = NULL,
    inner_max_iter = NULL, inner_reltol = NULL,

    # State cache
    last_theta  = NULL,
    last_solver = NULL,
    last_u      = NULL,   # warm-start: optimal u from previous outer iteration

    # Optimisation history: list of list(iter, params, cost)
    history = NULL,

    initialize = function(func_rhs, times_sim, obs_times, obs_values, y0,
                          fixed_params, lambda, param_scales,
                          inner_max_iter = 200,
                          inner_reltol   = sqrt(.Machine$double.eps)) {
      self$func_rhs        <- func_rhs
      self$times_sim       <- times_sim
      self$obs_times       <- obs_times
      self$obs_values      <- obs_values
      self$y0              <- y0
      self$fixed_params    <- fixed_params
      self$lambda          <- lambda
      self$inner_solver_class <- OdeSystemSolver
      self$param_scales    <- param_scales
      self$inner_max_iter  <- inner_max_iter
      self$inner_reltol    <- inner_reltol
      self$history         <- list()
    },

    # --- Helper: map normalised theta vector to physical params list -----------
    get_physical_params = function(theta_norm, param_names) {
      curr <- self$fixed_params
      for (i in seq_along(param_names)) {
        name       <- param_names[i]
        curr[[name]] <- theta_norm[i] * self$param_scales[[name]]
      }
      return(curr)
    },

    # --- Numerical Jacobian of f w.r.t. physical parameters (nv x np) --------
    # Central-difference scheme; same convention as CascadingOdeSolver.
    get_param_jacobian = function(y, t, p_phys, param_names) {
      nv  <- length(y)
      np  <- length(param_names)
      eps <- 1e-7
      J   <- matrix(0, nv, np)
      for (j in seq_len(np)) {
        dth <- eps * max(abs(p_phys[[param_names[j]]]), 1)
        p_p <- p_phys; p_p[[param_names[j]]] <- p_phys[[param_names[j]]] + dth
        p_m <- p_phys; p_m[[param_names[j]]] <- p_phys[[param_names[j]]] - dth
        J[, j] <- (self$func_rhs(y, t, p_p) - self$func_rhs(y, t, p_m)) / (2 * dth)
      }
      return(J)
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
        p_phys <- self$get_physical_params(theta_norm, param_names)
        solver <- self$inner_solver_class$new(
          func_rhs   = self$func_rhs,
          times_sim  = self$times_sim,
          obs_times  = self$obs_times,
          obs_values = self$obs_values,
          params     = p_phys,
          lambda     = self$lambda
        )
        solver$optimize(y0       = self$y0,
                        u_init   = NULL,
                        max_iter = self$inner_max_iter,
                        reltol   = self$inner_reltol)

        self$last_u      <- as.vector(solver$u)
        self$last_theta  <- theta_norm
        self$last_solver <- solver

        # Record outer iteration (use actual y0 after NA processing)
        y0_eff <- solver$y[1L, ]
        j_val  <- solver$cost_function(as.vector(solver$u), y0_eff)
        p_vals <- theta_norm * unlist(self$param_scales[param_names])
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

      p_vals <- theta_norm * unlist(self$param_scales[param_names])
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
      p_phys <- self$get_physical_params(theta_norm, param_names)
      ns     <- s$n_steps
      np     <- length(param_names)

      grad_phys <- numeric(np)
      for (t in seq_len(ns - 1L)) {
        dt    <- s$dt_vec[t]
        Fth   <- self$get_param_jacobian(s$y[t, ], s$times_sim[t],
                                         p_phys, param_names)  # nv x np
        p_next <- s$p[t + 1L, ]                                # nv
        # Contribution: p[t+1]^T * f_theta * dt  (np vector)
        grad_phys <- grad_phys + as.vector(t(Fth) %*% p_next) * dt
      }

      scales <- unlist(self$param_scales[param_names])
      return(grad_phys * scales)
    },

    # =========================================================================
    # 3. Partial sensitivity matrix S[ns, nv, np]:
    #    S[t, v, j] = dy_{t,v} / d theta_j_physical  (u* held fixed)
    #
    # Pure forward sensitivity sweep:
    #   S[t+1] = (I + dt*f_y[t]) * S[t]  +  dt * f_theta(y[t], t; theta)
    #   S[0]   = 0
    #
    # 'Partial' because the envelope theorem treats u*(theta) as constant
    # when differentiating H_T w.r.t. theta.  This is consistent with the
    # outer gradient (outer_gradient / outer_gradient_sensitivity) but gives
    # a lower bound on the true uncertainty in y* relative to the total
    # sensitivity used by CascadingOdeSolver.
    # =========================================================================
    compute_sensitivity_matrix = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s      <- self$last_solver
      p_phys <- self$get_physical_params(theta_norm, param_names)
      ns     <- s$n_steps
      nv     <- s$n_vars
      np     <- length(param_names)

      S <- array(0, c(ns, nv, np))
      for (t in seq_len(ns - 1L)) {
        dt  <- s$dt_vec[t]
        Fy  <- s$get_jacobian(s$y[t, ], s$times_sim[t])
        Fth <- self$get_param_jacobian(s$y[t, ], s$times_sim[t],
                                       p_phys, param_names)
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

      scales <- unlist(self$param_scales[param_names])
      return(grad_phys * scales)
    },

    # =========================================================================
    # 4. Outer Parameter Optimisation  (L-BFGS-B)
    # =========================================================================
    optimize_parameters = function(init_theta_physical, param_names,
                                   lower_phys = NULL, upper_phys = NULL) {
      init_theta_norm <- init_theta_physical /
        unlist(self$param_scales[param_names])

      if (is.null(lower_phys)) {
        lower_norm <- rep(1e-5, length(param_names))
      } else {
        lower_norm <- lower_phys / unlist(self$param_scales[param_names])
      }

      if (is.null(upper_phys)) {
        upper_norm <- rep(Inf, length(param_names))
      } else {
        upper_norm <- upper_phys / unlist(self$param_scales[param_names])
      }

      self$history <- list()

      cat("=== Starting Constrained Parameter Tracking (L-BFGS-B) ===\n")
      cat("Initial Guess (Norm):", round(init_theta_norm, 4), "\n")

      res <- optim(
        par         = init_theta_norm,
        fn          = self$outer_objective,
        gr          = self$outer_gradient,
        param_names = param_names,
        method      = "L-BFGS-B",
        lower       = lower_norm,
        control     = list(maxit = 50, trace = 1)
      )

      final_params <- res$par * unlist(self$param_scales[param_names])
      cat("\n=== Optimization Complete ===\n")
      print(final_params)
      return(final_params)
    }
  )
)
