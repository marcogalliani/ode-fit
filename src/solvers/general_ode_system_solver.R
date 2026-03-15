library(R6)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)

source("src/solvers/bvp_solver.R")


## General Physics-Informed Smoother----
OdeSystemSolver <- R6Class("OdeSystemSolver",
  public = list(
    func_rhs = NULL, params = NULL, lambda = NULL,
    times_sim = NULL,      # Simulation Time Grid
    n_steps = NULL, 
    dt_vec = NULL,         # Vector of time steps (can be variable)
    n_vars = NULL,
    
    observations_mapped = NULL, # Matrix matching times_sim (mostly NAs)
    
    y = NULL, u = NULL, p = NULL,
    
    # Initialize: Maps sparse data to the fine grid
    initialize = function(func_rhs, times_sim, obs_times, obs_values, params, lambda) {
      self$func_rhs <- func_rhs
      
      obs_times  <- round(obs_times,  digits = 10)
      times_sim  <- sort(unique(round(c(times_sim, obs_times), digits = 10)))
      self$times_sim <- times_sim
      self$params <- params
      self$lambda <- lambda
      self$n_steps <- length(times_sim)
      self$n_vars <- ncol(obs_values)

      self$dt_vec <- c(diff(times_sim), 0)

      self$observations_mapped <- matrix(NA, nrow = self$n_steps, ncol = self$n_vars)
      self$observations_mapped[times_sim %in% obs_times, ] <- obs_values
      # Initialize State Matrices
      self$y <- matrix(0, self$n_steps, self$n_vars)
      self$u <- matrix(0, self$n_steps, self$n_vars)
      self$p <- matrix(0, self$n_steps, self$n_vars)
    },
    
    # Numerical Jacobian (Central Difference)
    get_jacobian = function(y_vec, t_val) {
      n <- length(y_vec)
      J <- matrix(0, n, n)
      eps <- 1e-7
      f0 <- self$func_rhs(y_vec, t_val, self$params)
      for (j in 1:n) {
        y_p <- y_vec; y_p[j] <- y_p[j] + eps
        y_m <- y_vec; y_m[j] <- y_m[j] - eps
        f_p <- self$func_rhs(y_p, t_val, self$params)
        f_m <- self$func_rhs(y_m, t_val, self$params)
        J[, j] <- (f_p - f_m) / (2 * eps)
      }
      return(J)
    },
    
    # Forward Solver (Explicit Euler with Variable dt)
    # -> dy/dt = f(t,y)  with non-linear f
    solve_state = function(u_mat, y0) {
      y_new <- matrix(0, self$n_steps, self$n_vars)
      y_new[1, ] <- y0
      
      for (t in 1:(self$n_steps - 1)) {
        y_prev <- y_new[t, ]
        dt <- self$dt_vec[t]
        
        dy <- self$func_rhs(y_prev, self$times_sim[t], self$params)
        y_new[t+1, ] <- y_prev + dt * (dy + u_mat[t, ])
      }
      return(y_new)
    },
    
    # Backward Solver (Adjoint) - Handles NAs and Variable dt
    # -dp/dt = J^Tp + res
    solve_adjoint = function(y_curr) {
      p_new <- matrix(0, self$n_steps, self$n_vars)

      # terminal condition is p[T] = (2/n_steps) * r[T]
      # if r[T] = 0, then the terminal conditional is 0
      resid_T <- y_curr[self$n_steps, ] - self$observations_mapped[self$n_steps, ]
      resid_T[is.na(resid_T)] <- 0
      p_new[self$n_steps, ] <- (2 / self$n_steps) * resid_T

      # Backward loop from T-1 down to 1
      # Adjoint update: p[t] = (I + dt*J^T) * p[t+1]  +  (2/n_steps) * r[t]
      for (t in (self$n_steps - 1):1) {
        dt    <- self$dt_vec[t]
        y_now <- y_curr[t, ]

        # 1. Residual at current time t
        resid <- y_now - self$observations_mapped[t, ]
        resid[is.na(resid)] <- 0

        forcing <- (2 / self$n_steps) * resid

        # 2. Jacobian at current state
        J <- self$get_jacobian(y_now, self$times_sim[t])

        # 3. Adjoint update
        grad_prop  <- t(J) %*% p_new[t + 1, ]
        p_new[t, ] <- p_new[t + 1, ] + dt * as.vector(grad_prop) + forcing
      }
      return(p_new)
    },
    
    # Cost & Gradient
    cost_function = function(u_flat, y0) {
      u_mat <- matrix(u_flat, self$n_steps, self$n_vars)
      y_sim <- self$solve_state(u_mat, y0)

      # Guard: forward Euler can blow up at extreme parameters
      if (!all(is.finite(y_sim))) return(1e20)

      # SSE only on non-NA slots
      sse <- sum((self$observations_mapped - y_sim)^2, na.rm = TRUE)
      reg <- self$lambda * sum(u_mat^2)
      return((1/self$n_steps) * sse + reg)
    },
    
    gradient_function = function(u_flat, y0) {
      u_mat <- matrix(u_flat, self$n_steps, self$n_vars)

      # Forward pass: get the state trajectory for the current u.
      y_sim <- self$solve_state(u_mat, y0)

      # Guard: ODE blew up — return regularisation gradient (pushes u → 0)
      if (!all(is.finite(y_sim))) return(2 * self$lambda * u_flat)

      # Build the linear BVP that jointly represents
      #   Forward:  y[t+1] = A[t]*y[t] + b[t]              (C = 0: u is fixed)
      #   Backward: p[t]   = D[t]*p[t+1] + E[t]*y[t]
      # and solve it with solve_linear_bvp_riccati to obtain (Y = y, P = adjoint).
      # With C = 0 the Riccati matrix R[t] = 0 throughout, so the forward sweep
      # reproduces y_sim exactly and the backward substitution gives the adjoint.
      ns <- self$n_steps; nv <- self$n_vars

      obs_mask  <- matrix(as.numeric(!is.na(self$observations_mapped)), ns, nv)
      obs_clean <- self$observations_mapped; obs_clean[is.na(obs_clean)] <- 0

      A_list <- vector("list", ns - 1L)
      C_list <- vector("list", ns - 1L)
      D_list <- vector("list", ns - 1L)
      E_list <- vector("list", ns - 1L)
      b_list <- vector("list", ns - 1L)
      f_list <- vector("list", ns - 1L)

      for (t in seq_len(ns - 1L)) {
        dt <- self$dt_vec[t]
        Jt <- self$get_jacobian(y_sim[t, ], self$times_sim[t])
        At <- diag(nv) + dt * Jt

        A_list[[t]] <- At
        C_list[[t]] <- matrix(0, nv, nv)                        # no p→y coupling
        D_list[[t]] <- diag(nv) + dt * t(Jt)
        E_obs <- matrix(0, nv, nv)
        diag(E_obs) <- (2 / ns) * obs_mask[t, ]
        E_list[[t]] <- E_obs

        # b[t] encodes the (nonlinear) forward step via residual:
        #   b[t] = y[t+1] - A[t]*y[t]  so that A[t]*y[t] + b[t] = y[t+1]
        b_list[[t]] <- matrix(y_sim[t + 1L, ] - At %*% y_sim[t, ], nv, 1)
        # Backward forcing: -(2/ns)*mask[t]*obs_clean[t] so that
        # E[t]*y[t] + f[t] = (2/ns)*mask[t]*(y[t] - obs_clean[t])
        f_list[[t]] <- matrix(-(2 / ns) * obs_mask[t, ] * obs_clean[t, ], nv, 1)
      }

      E_T <- matrix(0, nv, nv)
      diag(E_T) <- (2 / ns) * obs_mask[ns, ]
      f_T <- matrix(-(2 / ns) * obs_mask[ns, ] * obs_clean[ns, ], nv, 1)

      bvp <- solve_linear_bvp_riccati(
        ns = ns, ny = nv, nrhs = 1L,
        A_list = A_list, C_list = C_list, b_list = b_list,
        D_list = D_list, E_list = E_list, f_list = f_list,
        y0_mat = matrix(y0, nv, 1),
        E_T    = E_T, f_T = f_T
      )

      p_sim <- matrix(bvp$P[, , 1L], ns, nv)

      # dJ/du[t] = 2*lambda*u[t] + dt[t]*p[t+1]
      p_shifted <- rbind(p_sim[-1, , drop = FALSE], matrix(0, 1, nv))
      dt_col    <- matrix(self$dt_vec, nrow = ns, ncol = nv)
      grad      <- 2 * self$lambda * u_mat + p_shifted * dt_col

      return(as.vector(grad))
    },
    
    optimize = function(y0, max_iter = 100, u_init = NULL, reltol = sqrt(.Machine$double.eps)) {
      # u_init: optional warm-start forcing (flat vector, length n_steps*n_vars).
      # Passing the previous optimal u avoids restarting from zero at each outer
      # iteration, which keeps the inner solver on the same local-minimum branch
      # and makes y*(theta) vary smoothly with theta.
      # When NULL (first call or cold start), fall back to u=0.
      if (is.null(u_init)) {
        u_init <- rep(0, self$n_steps * self$n_vars)
      }
      
      # If y0 is not provided (scalar NA), initialise from the first observation row
      if (length(y0) == 1 && is.na(y0)) {
        first_row <- which(!is.na(self$observations_mapped[, 1]))[1]
        y0 <- self$observations_mapped[first_row, ]
        y0[is.na(y0)] <- 0   # fallback for variables without an obs at t=0
      } else if (any(is.na(y0))) {
        # Vector y0 with some NA entries: fill per-variable from first available obs
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
      self$y <- self$solve_state(self$u, y0)
      self$p <- self$solve_adjoint(self$y)
      return(res)
    }
  )
)
