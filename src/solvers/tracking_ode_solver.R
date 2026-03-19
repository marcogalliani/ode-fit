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
    # Total sensitivity matrix S[ns, nv, np]:
    #   S[t, v, j] = dy*_{t,v} / d theta_j   (TOTAL derivative)
    #
    # Unlike the partial forward sweep (compute_sensitivity_matrix), this
    # accounts for the indirect path  theta -> u*(theta) -> y*  using the
    # same Riccati BVP as CascadingOdeSolver$compute_sensitivity_matrix.
    # The BVP couples the forward state sensitivity with a backward adjoint
    # that captures how the optimal forcing u* shifts when theta changes.
    # Used exclusively in compute_uncertainty; gradient computation continues
    # to use the partial sensitivity (envelope theorem).
    # =========================================================================
    compute_total_sensitivity_matrix = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s      <- self$last_solver
      p_phys <- self$get_physical_params(theta_norm, param_names)
      ns     <- s$n_steps
      nv     <- s$n_vars
      np     <- length(param_names)
      lambda <- self$lambda

      obs_mask <- matrix(as.numeric(!is.na(s$observations_mapped)), ns, nv)

      A_list  <- vector("list", ns - 1L)
      D_list  <- vector("list", ns - 1L)
      B_list  <- vector("list", ns - 1L)
      E_list  <- vector("list", ns - 1L)
      F_list  <- vector("list", ns - 1L)
      c_vec   <- numeric(ns - 1L)
      zero_np <- matrix(0, nv, np)

      for (t in seq_len(ns - 1L)) {
        dt  <- s$dt_vec[t]
        Fy  <- s$get_jacobian(s$y[t, ], s$times_sim[t])
        Fth <- self$get_param_jacobian(s$y[t, ], s$times_sim[t],
                                       p_phys, param_names)
        A_list[[t]] <- diag(nv) + dt * Fy
        D_list[[t]] <- diag(nv) + dt * t(Fy)
        B_list[[t]] <- dt * Fth
        c_vec[t]    <- dt^2 / (2 * lambda)

        E_obs <- matrix(0, nv, nv)
        diag(E_obs) <- (2 / ns) * obs_mask[t, ]
        E_list[[t]] <- E_obs
        F_list[[t]] <- zero_np
      }

      C_list <- lapply(c_vec, function(ct) ct * diag(nv))
      E_T    <- matrix(0, nv, nv)
      diag(E_T) <- (2 / ns) * obs_mask[ns, ]

      bvp <- solve_linear_bvp_riccati(
        ns     = ns, ny = nv, nrhs = np,
        A_list = A_list, C_list = C_list, b_list = B_list,
        D_list = D_list, E_list = E_list, f_list = F_list,
        y0_mat = matrix(0, nv, np),
        E_T    = E_T, f_T = NULL
      )
      return(bvp$Y)   # [ns, nv, np]
    },

    # =========================================================================
    # Inner criterion Hessian w.r.t. y — model uncertainty and corrected DOF.
    #
    # The inner criterion  J(u) = (1/n_s)||y(u)-y_obs||^2 + lambda*||u||^2
    # viewed as a function of the trajectory y, has Hessian
    #
    #   H_yy = (2/n_s)*I_obs + 2*lambda * L^T L
    #
    # where L is the linearised ODE difference operator:
    #   L_{t,t}   = -(I/dt_t + J_t),   L_{t,t+1} = I/dt_t
    #
    # H_yy is block-tridiagonal and positive definite (for lambda > 0).
    # y_0 is treated as known (large diagonal pin on the first block).
    #
    # Solving  H_yy Z = (2/n_s) I_obs  (one column per observed entry) gives:
    #
    #   trace_hat  = sum_k Z[obs_k, k]   hat matrix trace = effective inner DOF
    #   var_inner  = rowSums(Z^2)         inner variance factor (mult. by sigma2)
    #
    # The corrected residual DOF for sigma2:
    #   dof = n_obs - n_params - trace_hat
    #
    # Returns: list(trace_hat, var_inner [ns x nv])
    # =========================================================================
    compute_inner_variance = function() {
      s      <- self$last_solver
      ns     <- s$n_steps
      nv     <- s$n_vars
      lambda <- self$lambda
      N      <- ns * nv
      obs_mask <- !is.na(s$observations_mapped)

      # Build H_yy [N x N], flat index = (t-1)*nv + v  (row-major)
      H_yy <- matrix(0, N, N)

      # Data term: (2/n_s) * I_obs
      for (t in seq_len(ns)) {
        rows <- ((t - 1L) * nv + 1L):(t * nv)
        H_yy[rows, rows] <- H_yy[rows, rows] + diag(2 / ns * obs_mask[t, ], nv)
      }

      # ODE regularization: 2*lambda * L^T L
      for (t in seq_len(ns - 1L)) {
        dt     <- s$dt_vec[t]
        Jt     <- s$get_jacobian(s$y[t, ], s$times_sim[t])
        L_tt   <- -(diag(nv) / dt + Jt)   # L block at (t, t)
        L_ttp1 <- diag(nv) / dt             # L block at (t, t+1)
        rows_t   <- ((t - 1L) * nv + 1L):(t * nv)
        rows_tp1 <- (t * nv + 1L):((t + 1L) * nv)

        H_yy[rows_t,   rows_t  ] <- H_yy[rows_t,   rows_t  ] + 2*lambda * crossprod(L_tt)
        H_yy[rows_tp1, rows_tp1] <- H_yy[rows_tp1, rows_tp1] + 2*lambda * crossprod(L_ttp1)
        off <- 2*lambda * crossprod(L_tt, L_ttp1)
        H_yy[rows_t,   rows_tp1] <- H_yy[rows_t,   rows_tp1] + off
        H_yy[rows_tp1, rows_t  ] <- H_yy[rows_tp1, rows_t  ] + t(off)
      }

      # Pin y_0 (t=1): initial condition is known, not estimated
      diag(H_yy)[1:nv] <- diag(H_yy)[1:nv] + 1e10

      # Observed flat indices (row-major: t slow, v fast)
      obs_flat <- integer(0)
      for (t in seq_len(ns)) for (v in seq_len(nv))
        if (obs_mask[t, v]) obs_flat <- c(obs_flat, (t - 1L) * nv + v)
      n_f <- length(obs_flat)

      # Build RHS: N x n_f,  column k has (2/n_s) at obs_flat[k]
      RHS <- matrix(0, N, n_f)
      for (k in seq_len(n_f)) RHS[obs_flat[k], k] <- 2 / ns

      # Solve H_yy Z = RHS  (block-tridiagonal: positive definite)
      Z <- tryCatch(solve(H_yy, RHS),
                    error = function(e) solve(H_yy + 1e-10 * diag(N), RHS))

      # Hat matrix trace = effective DOF consumed by inner optimisation
      trace_hat <- sum(Z[cbind(obs_flat, seq_len(n_f))])

      # Inner variance factor per state entry: var_inner * sigma2 = Var(y*)_inner
      var_inner <- matrix(rowSums(Z * Z), ns, nv, byrow = TRUE)

      list(trace_hat = trace_hat, var_inner = var_inner)
    },

    # =========================================================================
    # Laplace / NLS confidence bands for y*(t; theta_hat).
    #
    # Two sources of uncertainty are combined:
    #   1. Parameter uncertainty  — delta method with TOTAL sensitivity
    #      (Riccati BVP; accounts for theta -> u*(theta) -> y* path)
    #   2. Inner optimisation uncertainty — sandwich formula via H_yy^{-1}
    #      (captures how y* varies with data noise at fixed theta)
    #
    # The DOF for sigma^2 is corrected for the effective parameters used
    # by the inner optimisation:
    #   dof = n_obs - n_params - trace_hat
    # where trace_hat = tr(H_yy^{-1} * (2/n_s) * I_obs) is the hat matrix
    # trace of the inner problem.
    #
    # Returns: list with fields
    #   y_hat        — fitted trajectory [ns x nv]
    #   se_param     — SE from parameter uncertainty only
    #   se_inner     — SE from inner optimisation uncertainty
    #   se           — total SE = sqrt(se_param^2 + se_inner^2)
    #   lower/upper  — total CI (use for plotting)
    #   lower_param/upper_param — parameter-only CI (for comparison)
    #   Sigma_theta  — NLS parameter covariance [np x np]
    #   sigma2       — corrected observation noise  SSE / dof
    #   trace_hat    — effective inner DOF (hat matrix trace)
    #   dof          — residual degrees of freedom
    # =========================================================================
    compute_uncertainty = function(theta_norm = NULL, param_names, alpha = 0.05) {
      if (is.null(theta_norm)) theta_norm <- self$last_theta
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s  <- self$last_solver
      ns <- s$n_steps
      nv <- s$n_vars
      np <- length(param_names)

      # --- Inner Hessian (no sigma2 needed) ---
      inner     <- self$compute_inner_variance()
      trace_hat <- inner$trace_hat

      # --- Corrected sigma2 ---
      obs_mask <- !is.na(s$observations_mapped)
      sse      <- sum((s$y - s$observations_mapped)^2, na.rm = TRUE)
      n_obs    <- sum(obs_mask)
      dof      <- max(n_obs - np - trace_hat, 1)
      sigma2   <- sse / dof

      # --- Parameter uncertainty (total sensitivity via Riccati BVP) ---
      S_phys <- self$compute_total_sensitivity_matrix(theta_norm, param_names)
      S_mat  <- matrix(S_phys, ns * nv, np)
      J_obs  <- S_mat[as.vector(obs_mask), , drop = FALSE]
      JtJ    <- t(J_obs) %*% J_obs
      ridge  <- 1e-10 * max(abs(diag(JtJ)), 1)
      Sigma_theta <- sigma2 * solve(JtJ + ridge * diag(np))

      SL        <- S_mat %*% Sigma_theta
      var_param <- matrix(rowSums(SL * S_mat), ns, nv)

      # --- Inner uncertainty ---
      var_inner <- sigma2 * inner$var_inner

      # --- Combine ---
      se_param <- sqrt(pmax(var_param, 0))
      se_inner <- sqrt(pmax(var_inner, 0))
      se_total <- sqrt(pmax(var_param + var_inner, 0))

      z <- qnorm(1 - alpha / 2)
      list(
        y_hat       = s$y,
        se_param    = se_param,
        se_inner    = se_inner,
        se          = se_total,
        lower       = s$y - z * se_total,
        upper       = s$y + z * se_total,
        lower_param = s$y - z * se_param,
        upper_param = s$y + z * se_param,
        Sigma_theta = Sigma_theta,
        sigma2      = sigma2,
        trace_hat   = trace_hat,
        dof         = dof
      )
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
