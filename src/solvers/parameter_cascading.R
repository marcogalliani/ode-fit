library(R6)
library(ggplot2)
library(gridExtra)

source("src/solvers/general_ode_system_solver.R")
  
# Parameter cascading solver
CascadingOdeSolver <- R6Class("CascadingOdeSolver",
  public = list(
    inner_solver_class = NULL,
    times_sim = NULL, obs_times = NULL, obs_values = NULL,
    func_rhs = NULL, fixed_params = NULL, lambda = NULL,
    param_scales = NULL,
    
    # State Cache
    last_theta = NULL,
    last_solver = NULL,
    
    initialize = function(func_rhs, times_sim, obs_times, obs_values, fixed_params, lambda, param_scales) {
      self$func_rhs <- func_rhs
      self$times_sim <- times_sim
      self$obs_times <- obs_times
      self$obs_values <- obs_values
      self$fixed_params <- fixed_params
      self$lambda <- lambda
      self$inner_solver_class <- OdeSystemSolver
      self$param_scales <- param_scales
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
        solver$optimize(y0 = NA, max_iter = 200) 
        self$last_theta <- theta_norm
        self$last_solver <- solver
      }
      
      # Data Misfit (SSE)
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
      
      # --- Step A: Construct the Linearized Optimality Matrix (A) ---
      # This matrix represents d/dy (Grad_y J_inner)
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
      diag(A) + .Machine$double.eps * max(abs(diag(A))) # thikonov regularization for stability

      # --- Step B: Construct the Parametric Forcing (B) ---
      # This represents d/dtheta (Grad_y J_inner)
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
      # Solve the linearized optimality system: A * S = -B
      S <- solve(A, -B)
      
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
    # Optimization Routine
    optimize_parameters = function(init_theta_physical, param_names, 
                                   lower_phys = NULL, upper_phys = NULL) {
      
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
      
      cat("=== Starting Constrained Parameter Cascading (L-BFGS-B) ===\n")
      cat("Initial Guess (Norm):", round(init_theta_norm, 4), "\n")
      
      # 3. Run L-BFGS-B
      # This method handles box constraints and is robust for large-scale problems
      res <- optim(
        par = init_theta_norm,
        fn = self$outer_objective,
        gr = self$outer_gradient,
        param_names = param_names,
        method = "L-BFGS-B",
        lower = lower_norm,
        upper = upper_norm,
        control = list(
          maxit = 100, 
          factr = 1e7,   # Roughly 1e-8 relative tolerance
          pgtol = 1e-5,  # Projected gradient tolerance
          trace = 1      # Set to 1 to monitor convergence in the console
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