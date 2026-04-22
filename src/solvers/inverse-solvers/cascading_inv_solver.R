library(R6)
  
# Parameter cascading solver
CascadingInverseSolver <- R6Class("CascadingInverseSolver",
  inherit = InverseSolverBase,
  public = list(
    initialize = function(model, times_sim, obs_times, obs_values,
                          lambda,
                          inner_max_iter = 200,
                          inner_reltol = sqrt(.Machine$double.eps),
              inner_method = "gl2",
              verbose = FALSE) {
      self$initialize_estimator(model, times_sim, obs_times, obs_values,
                                lambda, inner_max_iter, inner_reltol,
                inner_method, verbose)
    },
    
    # 1. Outer Objective
    outer_objective = function(theta_norm, param_names) {
      # Caching logic
      if(!is.null(self$last_theta) && all(theta_norm == self$last_theta)) {
        solver <- self$last_solver
      } else {
        p_phys <- self$unpack_physical(theta_norm, param_names)
        y0_phys <- self$eval_init_state(p_phys)
        solver <- self$inner_solver_class$new(
          model = self$model, times_sim = self$times_sim,
          obs_times = self$obs_times, obs_values = self$obs_values,
          params = p_phys, lambda = self$lambda, method = self$inner_method,
          verbose = FALSE
        )
        
        solver$optimize(y0 = y0_phys, u_init = NULL, #self$last_u,
                        max_iter = self$inner_max_iter,
                        reltol = self$inner_reltol,
                        verbose = FALSE)

        self$last_u     <- as.vector(solver$u)   # save for next warm start
        self$last_theta <- theta_norm
        self$last_solver <- solver

        # Record this outer iteration
        p_vals <- theta_norm * self$get_scales_vector(param_names)
        resid_rec <- solver$y - solver$observations_mapped
        sse_rec   <- sum(resid_rec^2, na.rm = TRUE)
        self$history[[length(self$history) + 1]] <- list(
          iter   = length(self$history) + 1L,
          params = setNames(p_vals, param_names),
          sse    = sse_rec
        )
      }

      # Data Misfit (SSE); guard against NaN state (e.g. ODE blew up at extreme params)
      if (!all(is.finite(solver$y))) return(1e20)
      resid <- solver$y - solver$observations_mapped
      sse <- sum(resid^2, na.rm = TRUE)

      # Log output
      p_vals <- theta_norm * self$get_scales_vector(param_names)
      self$log_iter(length(self$history), "SSE", sse,
                    setNames(p_vals, param_names))

      return(sse)
    },
    # Computation of
    outer_gradient = function(theta_norm, param_names) {
      self$outer_gradient_sensitivity(theta_norm, param_names)
    },
    # Total sensitivity matrix S[ns, nv, np] computed from the
    # linearized optimality system A * S = -B using scheme-level
    # defect linearization.
    compute_sensitivity_matrix = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s      <- self$last_solver
      ns     <- s$n_steps
      nv     <- s$n_vars
      np     <- length(param_names)
      NT <- ns * nv

      p_phys <- self$unpack_physical(theta_norm, param_names)
      obs_flat <- as.vector(s$observations_mapped)
      dt_l <- c(0, s$dt_vec[seq_len(ns - 1L)])
      w_trap <- (dt_l + s$dt_vec) / 2

      A <- matrix(0, NT, NT)
      # Hessian of inner data term wrt states: (1/ns) * SSE(y).
      diag(A) <- (2 / ns) * as.numeric(!is.na(obs_flat))

      for (t in seq_len(ns - 1L)) {
        idx_c <- seq(t, NT, by = ns)
        idx_n <- seq(t + 1L, NT, by = ns)

        dU <- s$get_discrete_control_jacobian(t, include_theta = FALSE)
        Du_dyc <- dU$Du_dyc
        Du_dyn <- dU$Du_dyn
        wt <- w_trap[t]

        A[idx_c, idx_c] <- A[idx_c, idx_c] + 2 * self$lambda * wt * (t(Du_dyc) %*% Du_dyc)
        A[idx_n, idx_n] <- A[idx_n, idx_n] + 2 * self$lambda * wt * (t(Du_dyn) %*% Du_dyn)
        A[idx_c, idx_n] <- A[idx_c, idx_n] + 2 * self$lambda * wt * (t(Du_dyc) %*% Du_dyn)
        A[idx_n, idx_c] <- A[idx_n, idx_c] + 2 * self$lambda * wt * (t(Du_dyn) %*% Du_dyc)
      }
      diag(A) <- diag(A) + 1e-9 * max(1, max(abs(diag(A))))

      B <- matrix(0, NT, np)
      base_params <- s$params
      for (t in seq_len(ns - 1L)) {
        idx_c <- seq(t, NT, by = ns)
        idx_n <- seq(t + 1L, NT, by = ns)

        dU <- s$get_discrete_control_jacobian(t, param_names, include_theta = TRUE)
        Du_dyc <- dU$Du_dyc
        Du_dyn <- dU$Du_dyn
        Du_dtheta <- dU$Du_dtheta
        wt <- w_trap[t]

        for (j in seq_len(np)) {
          du_dtheta_j <- Du_dtheta[, j]
          B[idx_c, j] <- B[idx_c, j] + 2 * self$lambda * wt * as.vector(t(Du_dyc) %*% du_dtheta_j)
          B[idx_n, j] <- B[idx_n, j] + 2 * self$lambda * wt * as.vector(t(Du_dyn) %*% du_dtheta_j)

          # Additional cross term: d(Du_dy)/dtheta^T * u
          nm <- param_names[j]
          th0 <- base_params[[nm]]
          dth <- 1e-6 * max(abs(th0), 1)

          s$params[[nm]] <- th0 + dth
          dU_p <- s$get_discrete_control_jacobian(t, include_theta = FALSE)

          s$params[[nm]] <- th0 - dth
          dU_m <- s$get_discrete_control_jacobian(t, include_theta = FALSE)

          s$params[[nm]] <- th0

          dDu_dyc_dth <- (dU_p$Du_dyc - dU_m$Du_dyc) / (2 * dth)
          dDu_dyn_dth <- (dU_p$Du_dyn - dU_m$Du_dyn) / (2 * dth)
          u_t <- s$u[t, ]

          B[idx_c, j] <- B[idx_c, j] + 2 * self$lambda * wt * as.vector(t(dDu_dyc_dth) %*% u_t)
          B[idx_n, j] <- B[idx_n, j] + 2 * self$lambda * wt * as.vector(t(dDu_dyn_dth) %*% u_t)
        }
      }
      s$params <- base_params

      J0 <- self$init_state_jacobian_fd(p_phys, param_names)
      idx_t0 <- seq(1L, NT, by = ns)

      # y(t0) is fixed by init_state(theta): enforce this as a hard boundary
      # condition by solving only on free state blocks (t >= 2).
      idx_free <- setdiff(seq_len(NT), idx_t0)
      A_ff <- A[idx_free, idx_free, drop = FALSE]
      A_f0 <- A[idx_free, idx_t0, drop = FALSE]

      S_flat <- matrix(0, NT, np)
      for (j in seq_len(np)) {
        rhs_j <- -B[idx_free, j] - as.vector(A_f0 %*% J0[, j])
        s_free_j <- tryCatch(
          solve(A_ff, rhs_j),
          error = function(e) qr.solve(A_ff, rhs_j)
        )
        S_flat[idx_free, j] <- s_free_j
        S_flat[idx_t0, j] <- J0[, j]
      }

      S_arr <- array(0, c(ns, nv, np))
      for (j in seq_len(np)) {
        for (t in seq_len(ns)) {
          idx_t <- seq(t, NT, by = ns)
          S_arr[t, , j] <- S_flat[idx_t, j]
        }
      }

      S_arr
    },

    # Outer gradient via sensitivity equations — delegates to
    # compute_sensitivity_matrix and applies the chain rule.
    outer_gradient_sensitivity = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s <- self$last_solver
      if (is.null(s) || is.null(s$y) || !all(is.finite(s$y))) {
        return(rep(0, length(param_names)))
      }

      np <- length(param_names)
      resid <- ifelse(is.na(s$observations_mapped), 0,
                      s$y - s$observations_mapped)

      S_arr <- self$compute_sensitivity_matrix(theta_norm, param_names)

      grad_phys <- vapply(
        seq_len(np),
        function(j) sum(2 * resid * S_arr[, , j]),
        numeric(1L)
      )

      scales <- self$get_scales_vector(param_names)
      grad_phys * scales
    },

    # Optimization Routine
    optimize_parameters = function(init_theta_physical, param_names,
                                   lower_phys = NULL, upper_phys = NULL) {
      prep <- self$prepare_theta_normalized(param_names,
                                            init_theta_physical,
                                            lower_phys,
                                            upper_phys)
      scales <- prep$scales
      init_theta_norm <- prep$init
      lower_norm <- prep$lower
      upper_norm <- prep$upper
      
      self$history <- list()   # reset trace for this run

      self$log_info("Starting constrained parameter cascading (L-BFGS-B)")
      if (isTRUE(self$verbose)) {
        self$log_info("Initial guess (normalized): %s",
                      paste(sprintf("%.6g", init_theta_norm), collapse = ", "))
      }

      old_warn <- getOption("warn")
      options(warn = 0)
      on.exit(options(warn = old_warn), add = TRUE)
      
      # 3. Run L-BFGS-B
      res <- optim(
        par = init_theta_norm,
        fn = self$outer_objective,
        gr = self$outer_gradient_sensitivity,
        param_names = param_names,
        method = "L-BFGS-B",
        lower = lower_norm,
        upper = upper_norm,
        control = list(
          maxit = 50,
          trace = if (isTRUE(self$verbose)) 1 else 0
        )
      )
      
      # 4. De-normalize and return
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