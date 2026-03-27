# =============================================================================
# src/solvers/sqp_ocp_solver.R
#
# SQP (Sequential Quadratic Programming) solver for the ODE-constrained
# optimal control problem using Euler discretization.
#
# Problem:
#   min_{u}  (1/ns)*SSE(y) + lambda * integral(u^2) dt
#   s.t.     y[t+1] = y[t] + h[t]*(f(y[t],t) + u[t]),  y[1] = y0
#
# Works in the full (y, u, p) KKT space. Each iteration:
#   1. Evaluate KKT residuals
#   2. Assemble the linearised KKT matrix (Gauss-Newton Hessian)
#   3. Solve the linear system for (delta_y, delta_u, delta_p)
#   4. Line-search on KKT norm, update (y, u, p)
#
# Lagrangian convention:  L = J - sum p[t+1]^T * c_t
# where c_t = y[t+1] - y[t] - h[t]*(f(y[t],t) + u[t])
#
# KKT conditions (stationarity of L):
#   State:    y[t+1] - y[t] - h[t]*(f(y[t],t) + u[t]) = 0
#   Adjoint:  p[t] = p[t+1] + h[t]*J[t]^T*p[t+1] + (2/ns)*mask[t]*(y[t]-obs[t])
#   Terminal: p[ns] = (2/ns)*mask[ns]*(y[ns]-obs[ns])
#   Control:  2*lambda*w[t]*u[t] + h[t]*p[t+1] = 0
# =============================================================================

library(R6)

SqpOcpSolver <- R6Class("SqpOcpSolver",

  private = list(
    obs_clean = NULL,
    obs_mask  = NULL,
    w_trap    = NULL,

    get_jacobian = function(y_vec, t_val) {
      n   <- length(y_vec)
      J   <- matrix(0, n, n)
      eps <- 1e-7
      for (j in seq_len(n)) {
        y_p <- y_vec; y_p[j] <- y_p[j] + eps
        y_m <- y_vec; y_m[j] <- y_m[j] - eps
        J[, j] <- (self$func_rhs(y_p, t_val, self$params) -
                   self$func_rhs(y_m, t_val, self$params)) / (2 * eps)
      }
      J
    },

    kkt_residuals = function(y, u, p, J_cache = NULL, f_cache = NULL) {
      ns <- self$n_steps; ny <- self$n_vars
      h  <- self$dt_vec[seq_len(ns - 1L)]
      two_over_ns <- 2 / ns

      R_state <- matrix(0, ns - 1L, ny)
      R_adj   <- matrix(0, ns, ny)
      R_ctrl  <- matrix(0, ns - 1L, ny)

      for (t in seq_len(ns - 1L)) {
        f_t <- if (!is.null(f_cache)) f_cache[[t]]
               else self$func_rhs(y[t, ], self$times_sim[t], self$params)
        J_t <- if (!is.null(J_cache)) J_cache[[t]]
               else private$get_jacobian(y[t, ], self$times_sim[t])

        R_state[t, ] <- y[t + 1L, ] - y[t, ] - h[t] * (f_t + u[t, ])
        R_adj[t, ]   <- two_over_ns * private$obs_mask[t, ] *
                           (y[t, ] - private$obs_clean[t, ]) -
                         p[t, ] + p[t + 1L, ] +
                         h[t] * as.vector(t(J_t) %*% p[t + 1L, ])
        R_ctrl[t, ]  <- 2 * self$lambda * private$w_trap[t] * u[t, ] +
                         h[t] * p[t + 1L, ]
      }

      R_adj[ns, ] <- two_over_ns * private$obs_mask[ns, ] *
                       (y[ns, ] - private$obs_clean[ns, ]) - p[ns, ]

      list(R_state = R_state, R_adj = R_adj, R_ctrl = R_ctrl,
           R_ic = y[1L, ] - self$y0)
    },

    kkt_norm = function(res) {
      max(abs(res$R_state), abs(res$R_adj), abs(res$R_ctrl), abs(res$R_ic))
    }
  ),

  public = list(
    func_rhs  = NULL, params = NULL, lambda = NULL,
    times_sim = NULL, n_steps = NULL, dt_vec = NULL, n_vars = NULL,
    y0 = NULL,

    observations_mapped = NULL,
    y = NULL, u = NULL, p = NULL,

    initialize = function(func_rhs, times_sim, obs_times, obs_values,
                          params, lambda, y0) {
      self$func_rhs <- func_rhs
      obs_times  <- round(obs_times, digits = 10)
      times_sim  <- sort(unique(round(c(times_sim, obs_times), digits = 10)))
      self$times_sim <- times_sim
      self$params    <- params
      self$lambda    <- lambda
      self$n_steps   <- length(times_sim)
      self$n_vars    <- ncol(obs_values)
      self$dt_vec    <- c(diff(times_sim), 0)
      self$y0        <- y0

      self$observations_mapped <- matrix(NA, self$n_steps, self$n_vars)
      self$observations_mapped[times_sim %in% obs_times, ] <- obs_values

      private$obs_clean <- self$observations_mapped
      private$obs_clean[is.na(private$obs_clean)] <- 0
      private$obs_mask  <- matrix(as.numeric(!is.na(self$observations_mapped)),
                                  self$n_steps, self$n_vars)

      dt_l <- c(0, self$dt_vec[seq_len(self$n_steps - 1L)])
      private$w_trap <- (dt_l + self$dt_vec) / 2

      self$y <- matrix(0, self$n_steps, self$n_vars)
      self$u <- matrix(0, self$n_steps, self$n_vars)
      self$p <- matrix(0, self$n_steps, self$n_vars)
    },

    solve = function(max_iter = 50L, tol = 1e-8, verbose = FALSE) {
      ns  <- self$n_steps; ny <- self$n_vars
      h   <- self$dt_vec[seq_len(ns - 1L)]
      lam <- self$lambda
      two_over_ns <- 2 / ns
      y0  <- self$y0
      I_ny <- diag(ny)

      # unknowns: dy[2..ns], du[1..ns-1], dp[1..ns]
      N <- (3L * ns - 2L) * ny

      # column offsets
      off_y <- 0L
      off_u <- (ns - 1L) * ny
      off_p <- 2L * (ns - 1L) * ny

      cy <- function(t) off_y + (t - 2L) * ny + seq_len(ny)  # dy[t], t=2..ns
      cu <- function(t) off_u + (t - 1L) * ny + seq_len(ny)  # du[t], t=1..ns-1
      cp <- function(t) off_p + (t - 1L) * ny + seq_len(ny)  # dp[t], t=1..ns

      # row offsets (same block ordering: state, adjoint, control)
      off_sr <- 0L
      off_ar <- (ns - 1L) * ny
      off_cr <- (ns - 1L) * ny + ns * ny

      sr <- function(t) off_sr + (t - 1L) * ny + seq_len(ny)  # state, t=1..ns-1
      ar <- function(t) off_ar + (t - 1L) * ny + seq_len(ny)  # adjoint, t=1..ns
      cr <- function(t) off_cr + (t - 1L) * ny + seq_len(ny)  # control, t=1..ns-1

      # --- initial guess: Euler forward (u=0) + consistent adjoint ---
      dto     <- make_dto_solver("euler")
      rhs_fwd <- function(y, t) self$func_rhs(y, t, self$params)
      jac_fn  <- function(y, t) private$get_jacobian(y, t)

      y_curr <- dto$solve_state(rhs_fwd, y0, self$times_sim, jac_fn)$y
      u_curr <- matrix(0, ns, ny)

      source_fn <- function(t_idx) {
        r <- y_curr[t_idx, ] - private$obs_clean[t_idx, ]
        r[!as.logical(private$obs_mask[t_idx, ])] <- 0
        two_over_ns * r
      }
      p_curr <- dto$solve_adjoint(rhs_fwd, rep(0, ny), jac_fn, source_fn)$p

      converged <- FALSE
      kkt_n <- Inf

      for (iter in seq_len(max_iter)) {
        J_cache <- vector("list", ns - 1L)
        f_cache <- vector("list", ns - 1L)
        for (t in seq_len(ns - 1L)) {
          J_cache[[t]] <- private$get_jacobian(y_curr[t, ], self$times_sim[t])
          f_cache[[t]] <- self$func_rhs(y_curr[t, ], self$times_sim[t], self$params)
        }

        res   <- private$kkt_residuals(y_curr, u_curr, p_curr, J_cache, f_cache)
        kkt_n <- private$kkt_norm(res)

        if (verbose)
          cat(sprintf("  [SQP] iter %2d  |KKT|_inf = %.3e\n", iter, kkt_n))
        if (kkt_n < tol) { converged <- TRUE; break }

        # --- assemble linearised KKT system ---
        K   <- matrix(0, N, N)
        rhs <- numeric(N)

        for (t in seq_len(ns - 1L)) {
          J_t  <- J_cache[[t]]
          JT_t <- t(J_t)
          w_t  <- private$w_trap[t]

          # state:  dy[t+1] - (I+h*J)*dy[t] - h*du[t] = -R_state[t]
          K[sr(t), cy(t + 1L)] <- I_ny
          if (t >= 2L) K[sr(t), cy(t)] <- -(I_ny + h[t] * J_t)
          K[sr(t), cu(t)] <- -h[t] * I_ny
          rhs[sr(t)] <- -res$R_state[t, ]

          # adjoint at t: (2/ns)*diag(mask)*dy[t] - dp[t] + (I+h*J^T)*dp[t+1] = -R_adj[t]
          if (t >= 2L)
            K[ar(t), cy(t)] <- two_over_ns * diag(private$obs_mask[t, ], nrow = ny)
          K[ar(t), cp(t)]      <- -I_ny
          K[ar(t), cp(t + 1L)] <- I_ny + h[t] * JT_t
          rhs[ar(t)] <- -res$R_adj[t, ]

          # control: 2*lam*w*du[t] + h*dp[t+1] = -R_ctrl[t]
          K[cr(t), cu(t)]      <- 2 * lam * w_t * I_ny
          K[cr(t), cp(t + 1L)] <- h[t] * I_ny
          rhs[cr(t)] <- -res$R_ctrl[t, ]
        }

        # terminal adjoint: (2/ns)*diag(mask)*dy[ns] - dp[ns] = -R_adj[ns]
        K[ar(ns), cy(ns)] <- two_over_ns * diag(private$obs_mask[ns, ], nrow = ny)
        K[ar(ns), cp(ns)] <- -I_ny
        rhs[ar(ns)] <- -res$R_adj[ns, ]

        # --- solve linear system ---
        delta <- solve(K, rhs)

        dy <- matrix(0, ns, ny)
        du <- matrix(0, ns, ny)
        dp <- matrix(0, ns, ny)
        for (t in 2L:ns)         dy[t, ] <- delta[cy(t)]
        for (t in seq_len(ns-1)) du[t, ] <- delta[cu(t)]
        for (t in seq_len(ns))   dp[t, ] <- delta[cp(t)]

        # --- Armijo backtracking on KKT norm ---
        step <- 1.0
        for (k in seq_len(15L)) {
          kkt_try <- private$kkt_norm(
            private$kkt_residuals(y_curr + step * dy,
                                  u_curr + step * du,
                                  p_curr + step * dp))
          if (kkt_try < kkt_n) break
          step <- step * 0.5
        }

        y_curr <- y_curr + step * dy
        u_curr <- u_curr + step * du
        p_curr <- p_curr + step * dp
      }

      self$y <- y_curr; self$u <- u_curr; self$p <- p_curr

      list(converged = converged, iter = iter, kkt_norm = kkt_n)
    },

    cost_function = function() {
      ns  <- self$n_steps
      sse <- sum((self$observations_mapped - self$y)^2, na.rm = TRUE)
      reg <- self$lambda * sum(private$w_trap * rowSums(self$u^2))
      (1 / ns) * sse + reg
    }
  )
)
