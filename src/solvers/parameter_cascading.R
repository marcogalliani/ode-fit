library(R6)
library(ggplot2)
library(gridExtra)

source("src/solvers/general_ode_system_solver.R")
  
# Parameter cascading solver
CascadingOdeSolver <- R6Class("CascadingOdeSolver",
  public = list(
    inner_solver_class = NULL,
    times_sim = NULL, obs_times = NULL, obs_values = NULL, y0 = NULL,
    func_rhs = NULL, fixed_params = NULL, lambda = NULL,
    param_scales = NULL,
    inner_max_iter = NULL, inner_reltol = NULL,

    # State Cache
    last_theta = NULL,
    last_solver = NULL,
    last_u = NULL,   # warm-start: optimal u from previous outer iteration

    # Optimisation history: list of list(iter, params, sse)
    history = NULL,

    initialize = function(func_rhs, times_sim, obs_times, obs_values, y0,
                          fixed_params, lambda, param_scales,
                          inner_max_iter = 200, inner_reltol = sqrt(.Machine$double.eps)) {
      self$func_rhs <- func_rhs
      self$times_sim <- times_sim
      self$obs_times <- obs_times
      self$obs_values <- obs_values
      self$y0 <- y0
      self$fixed_params <- fixed_params
      self$lambda <- lambda
      self$inner_solver_class <- OdeSystemSolver
      self$param_scales <- param_scales
      self$inner_max_iter <- inner_max_iter
      self$inner_reltol   <- inner_reltol
      self$history <- list()
    },
    
    # --- Helper: Map vector to params ---
    get_physical_params = function(theta_norm, param_names) {
      curr <- self$fixed_params
      for(i in seq_along(param_names)) {
        name <- param_names[i]
        curr[[name]] <- theta_norm[i] * self$param_scales[[name]]
      }
      return(curr)
    },
    
    # 1. Outer Objective
    outer_objective = function(theta_norm, param_names) {
      # Caching logic
      if(!is.null(self$last_theta) && all(theta_norm == self$last_theta)) {
        solver <- self$last_solver
      } else {
        p_phys <- self$get_physical_params(theta_norm, param_names)
        solver <- self$inner_solver_class$new(
          func_rhs = self$func_rhs, times_sim = self$times_sim,
          obs_times = self$obs_times, obs_values = self$obs_values,
          params = p_phys, lambda = self$lambda
        )
        
        solver$optimize(y0 = self$y0, u_init = NULL, #self$last_u,
                        max_iter = self$inner_max_iter, reltol = self$inner_reltol)

        self$last_u     <- as.vector(solver$u)   # save for next warm start
        self$last_theta <- theta_norm
        self$last_solver <- solver

        # Record this outer iteration
        p_vals <- theta_norm * unlist(self$param_scales[param_names])
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
      p_vals <- theta_norm * unlist(self$param_scales[param_names])
      cat(sprintf("Iter | Params: %s | SSE: %.4f\n", paste(round(p_vals, 2), collapse=","), sse))

      return(sse)
    },
    # Computation of
    outer_gradient = function(theta_norm, param_names) {
      # 1. Sync the inner solver
      if(is.null(self$last_theta) || !all(theta_norm == self$last_theta)) {
        self$outer_objective(theta_norm, param_names)
      }
      s <- self$last_solver
      p_phys <- self$get_physical_params(theta_norm, param_names)
      
      nv <- s$n_vars
      ns <- s$n_steps
      NT <- ns * nv
      np <- length(theta_norm)
      eps <- 1e-7
      
      # --- Step A: Construct the d^2 J/dy^2 
      # It consists of the Data Misfit part and the Physics part.
      A <- matrix(0, NT, NT)
      
      # 1. Data Misfit Part: (2/ns) * I for observed indices
      obs_flat <- as.vector(s$observations_mapped)
      diag(A) <- (2 / ns) * as.numeric(!is.na(obs_flat))
      
      # 2. Physics Part: Linearizing the Smoothing Penalty
      for(t in 1:(ns - 1)) {
        dt <- s$dt_vec[t]
        idx_c <- seq(t, NT, by = ns)      # Indices for state at time t
        idx_n <- seq(t + 1, NT, by = ns)  # Indices for state at time t+1
        
        # Linearize the ODE dynamics: f_y = df/dy
        f_y <- s$get_jacobian(s$y[t, ], s$times_sim[t])
        
        # Sensitivities of the physics error u = (y_next - y_curr)/dt - f
        Du_dyc <- -diag(1/dt, nv) - f_y
        Du_dyn <-  diag(1/dt, nv)
        
        # Add the contribution to the system matrix (2 * lambda * Du^T * Du)
        A[idx_c, idx_c] <- A[idx_c, idx_c] + 2 * self$lambda * (t(Du_dyc) %*% Du_dyc)
        A[idx_n, idx_n] <- A[idx_n, idx_n] + 2 * self$lambda * (t(Du_dyn) %*% Du_dyn)
        A[idx_c, idx_n] <- A[idx_c, idx_n] + 2 * self$lambda * (t(Du_dyc) %*% Du_dyn)
        A[idx_n, idx_c] <- A[idx_n, idx_c] + 2 * self$lambda * (t(Du_dyn) %*% Du_dyc)
      }      
      diag(A) <- diag(A) + 1e-9 * max(abs(diag(A))) # Tikhonov regularization for stability

      # --- Step B: Construct d^2/dtheta dy  J_inner
      B <- matrix(0, NT, np)
      
      for(j in 1:np) {
        p_pert <- p_phys
        p_pert[[param_names[j]]] <- p_phys[[param_names[j]]] + eps
        
        for(t in 1:(ns - 1)) {
          dt <- s$dt_vec[t]
          idx_c <- seq(t, NT, by = ns)
          idx_n <- seq(t + 1, NT, by = ns)
          
          # Change in physics error w.r.t parameter: du/dtheta = -df/dtheta
          f_p <- self$func_rhs(s$y[t, ], s$times_sim[t], p_pert)
          f_0 <- self$func_rhs(s$y[t, ], s$times_sim[t], p_phys)
          du_dtheta <- -(f_p - f_0) / eps
          
          f_y <- s$get_jacobian(s$y[t, ], s$times_sim[t])
          Du_dyc <- -diag(1/dt, nv) - f_y
          Du_dyn <-  diag(1/dt, nv)
          
          # Cross-derivative contribution: 2 * lambda * (Du/dy)^T * (du/dtheta)
          B[idx_c, j] <- B[idx_c, j] + 2 * self$lambda * as.vector(t(Du_dyc) %*% du_dtheta)
          B[idx_n, j] <- B[idx_n, j] + 2 * self$lambda * as.vector(t(Du_dyn) %*% du_dtheta)
        }
      }
      
      # --- Step C: Solve for the Sensitivity S ---
      # Solve optimality system differentiated w.r.t. theta: A * S = -B
      S <- solve(A, -B)

      # y0 is fixed (does not depend on theta); enforce dY_0/dtheta = 0.
      # Without this, the IFT would assign nonzero sensitivity to the initial
      # condition, which is wrong.  Indices of t=1 in the flat state vector:
      idx_t0 <- seq(1L, NT, by = ns)
      S[idx_t0, ] <- 0
      
      # --- Step D: Total Outer Gradient ---
      # Objective: sum((y - obs)^2)
      # Derivative w.r.t y: 2 * (y - obs)
      y_flat <- as.vector(s$y)
      resid <- y_flat - obs_flat
      resid[is.na(resid)] <- 0
      grad_y_outer <- 2 * resid
      
      # Chain rule: dJ/dtheta = dJ/dy * dy/dtheta
      grad_theta_phys <- as.vector(t(grad_y_outer) %*% S)
      
      # Normalize back for the optimizer (BFGS expects gradient w.r.t theta_norm)
      final_grad <- numeric(np)
      for(j in seq_along(param_names)) {
        final_grad[j] <- grad_theta_phys[j] * self$param_scales[[param_names[j]]]
      }
      
      return(final_grad)
    },
    # Jacobian of func_rhs w.r.t. y at arbitrary params (nv x nv).
    # Needed for the second-order correction in outer_gradient_sensitivity.
    compute_jac_y_ext = function(y_vec, t_val, p_phys) {
      nv <- length(y_vec)
      J  <- matrix(0, nv, nv)
      h  <- 1e-5
      for (j in seq_len(nv)) {
        y_p <- y_vec; y_p[j] <- y_vec[j] + h
        y_m <- y_vec; y_m[j] <- y_vec[j] - h
        J[, j] <- (self$func_rhs(y_p, t_val, p_phys) -
                   self$func_rhs(y_m, t_val, p_phys)) / (2 * h)
      }
      return(J)
    },

    # Numerical Jacobian of f w.r.t. physical parameters theta_phys (nv x np).
    # Uses the same central-difference scheme as GeneralOdeSolver$get_jacobian.
    get_param_jacobian = function(y, t, p_phys, param_names) {
      nv  <- length(y)
      np  <- length(param_names)
      f0  <- self$func_rhs(y, t, p_phys)
      eps <- 1e-7
      J   <- matrix(0, nv, np)
      for (j in seq_len(np)) {
        dth <- eps * max(abs(p_phys[[param_names[j]]]), 1)
        p_p <- p_phys
        p_p[[param_names[j]]] <- p_phys[[param_names[j]]] + dth
        p_m <- p_phys
        p_m[[param_names[j]]] <- p_phys[[param_names[j]]] - dth
        J[, j] <- (self$func_rhs(y, t, p_p) - self$func_rhs(y, t, p_m)) / (2 * dth)
      }
      return(J)
    },

    # Outer gradient via sensitivity equations — Riccati-based direct BVP solve.
    #
    # The KKT differentiation yields a coupled two-point BVP in
    # (Y_t = dY_t/dtheta, P_t = dP_t/dtheta):
    #
    #   Forward:   Y[t+1] = A_t*Y[t] + B_t - c_t*P[t+1]       Y[1]=0
    #   Backward:  P[t]   = D_t*P[t+1] + E_t*Y[t]              P[T]=E_T*Y[T]
    #
    # where:
    #   A_t = I + dt*f_y,   D_t = A_t^T
    #   B_t = dt*f_theta             (nv x np, state sensitivity forcing)
    #   c_t = dt^2/(2*lambda)        (scalar, from KKT: S_u_t = -dt*P[t+1]/(2λ))
    #   E_t = (2/ns)*diag(mask_t)    (Gauss-Newton: SOC1 G_t and SOC2 F_t dropped)
    #   E_T = (2/ns)*diag(mask_T)
    #
    # Why SOC terms are dropped (Gauss-Newton approximation):
    #   SOC1 (G_t) and SOC2 (F_t) both scale as O(p_inner).  When parameters are
    #   far from optimal, p_inner is large.  The Riccati matrix R_t grows because
    #   spectral_radius(A_t*D_t) > 1 for marginally-stable ODEs (e.g. LV).
    #   R_t*F_t and E_t*R_t blow up exponentially, corrupting s_t and hence Y.
    #   Dropping both terms keeps E_t = O(1/ns) and gives a self-correcting
    #   mechanism: R_t*q_t ≈ -s_t at late steps, bounding the computation.
    #   The resulting gradient uses a PSD Gauss-Newton Hessian approximation and
    #   is always a descent direction.
    outer_gradient_sensitivity = function(theta_norm, param_names) {
      if (is.null(self$last_theta) || !all(theta_norm == self$last_theta))
        self$outer_objective(theta_norm, param_names)

      s      <- self$last_solver
      ns     <- s$n_steps
      nv     <- s$n_vars
      np     <- length(param_names)
      lambda <- self$lambda

      p_phys   <- self$get_physical_params(theta_norm, param_names)
      obs_mask <- matrix(as.numeric(!is.na(s$observations_mapped)), ns, nv)
      resid    <- ifelse(is.na(s$observations_mapped), 0,
                         s$y - s$observations_mapped)

      # --- Precompute per-step matrices (Gauss-Newton: F_t = 0) ---
      A_list <- vector("list", ns - 1L)   # nv x nv
      D_list <- vector("list", ns - 1L)   # nv x nv
      B_list <- vector("list", ns - 1L)   # nv x np
      E_list <- vector("list", ns - 1L)   # nv x nv
      F_list <- vector("list", ns - 1L)   # nv x np  (all zeros)
      c_vec  <- numeric(ns - 1L)

      zero_np <- matrix(0, nv, np)

      for (t in seq_len(ns - 1L)) {
        dt  <- s$dt_vec[t]
        Fy  <- s$get_jacobian(s$y[t, ], s$times_sim[t])
        Fth <- self$get_param_jacobian(s$y[t, ], s$times_sim[t], p_phys, param_names)

        A_list[[t]] <- diag(nv) + dt * Fy
        D_list[[t]] <- diag(nv) + dt * t(Fy)
        B_list[[t]] <- dt * Fth
        c_vec[t]    <- dt^2 / (2 * lambda)

        E_obs <- matrix(0, nv, nv); diag(E_obs) <- (2 / ns) * obs_mask[t, ]
        E_list[[t]] <- E_obs
        F_list[[t]] <- zero_np   # SOC2 dropped for stability
      }

      # --- Solve linear BVP via Riccati sweep ---
      # C[t] = c_t * I (KKT coupling: S_u_t = -dt*P[t+1]/(2λ))
      C_list <- lapply(c_vec, function(ct) ct * diag(nv))
      E_T    <- matrix(0, nv, nv); diag(E_T) <- (2 / ns) * obs_mask[ns, ]

      bvp <-solve_linear_bvp_riccati(
          ns     = ns, ny = nv, nrhs = np,
          A_list = A_list, C_list = C_list, b_list = B_list,
          D_list = D_list, E_list = E_list, f_list = F_list,
          y0_mat = matrix(0, nv, np),
          E_T    = E_T, f_T = NULL
        )

      Y_arr <- bvp$Y
      P_arr <- bvp$P

      # --- Gradient ---
      grad_phys <- vapply(seq_len(np),
                          function(j) sum(2 * resid * Y_arr[, , j]),
                          numeric(1L))

      scales <- unlist(self$param_scales[param_names])
      return(grad_phys * scales)
    },

    # Optimization Routine
    optimize_parameters = function(init_theta_physical, param_names, lower_phys = NULL, upper_phys = NULL) {
      
      # 1. Normalize initial guess
      init_theta_norm <- init_theta_physical / unlist(self$param_scales[param_names])
      
      # 2. Handle Bounds
      # If no bounds are provided, default to a very small positive number for stability
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
      
      self$history <- list()   # reset trace for this run

      cat("=== Starting Constrained Parameter Cascading (L-BFGS-B) ===\n")
      cat("Initial Guess (Norm):", round(init_theta_norm, 4), "\n")
      
      # 3. Run L-BFGS-B
      res <- optim(
        par = init_theta_norm,
        fn = self$outer_objective,
        gr = self$outer_gradient_sensitivity,
        param_names = param_names,
        method = "L-BFGS-B",
        lower = lower_norm,
        #upper = upper_norm,
        control = list(
          maxit = 50, 
          trace = 1     # Set to 1 to monitor convergence in the console
        )
      )
      
      # 4. De-normalize and return
      final_params <- res$par * unlist(self$param_scales[param_names])
      
      cat("\n=== Optimization Complete ===\n")
      print(final_params)
      return(final_params)
    }
  )
)