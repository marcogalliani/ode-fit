library(R6)
library(deSolve)
library(Matrix)
library(abind)
library(fda) # Required for B-spline smoothing

#### 1. Data Smoother (Pre-smoothing) ----
# Converts discrete, noisy observations into a continuous smooth curve y_hat(t).
DataSmoother <- R6Class("DataSmoother",
  public = list(
    fd_obj = NULL, # Functional Data Object
    times = NULL,
    initialize = function(times, observations) {
      self$times <- times

      # Create B-spline basis
      # Rule of thumb: n_basis approx equal to number of observations + order
      n_obs <- length(times)
      basis <- create.bspline.basis(rangeval = range(times), nbasis = max(4, n_obs + 2), norder = 4)

      # Smooth the data
      # lambda_smooth can be tuned, here set small to trust data
      fdPar <- fdPar(basis, Lfdobj = 2, lambda = 1e-4)

      # Smooth each dimension
      smooth_res <- smooth.basis(times, observations, fdPar)
      self$fd_obj <- smooth_res$fd
    },

    # Evaluate smoothed curve at specific times
    eval = function(t_eval) {
      eval.fd(t_eval, self$fd_obj)
    },

    # Evaluate derivative of smoothed curve (optional, useful for gradient matching init)
    eval_deriv = function(t_eval) {
      eval.fd(t_eval, self$fd_obj, Lfdobj = 1)
    }
  )
)

#### 2. Numerical Derivatives Helper ----
NumDeriv <- R6Class("NumDeriv",
  public = list(
    func = NULL,
    initialize = function(func) {
      self$func <- func
    },
    grad_x = function(y, t, p) {
      n <- length(y)
      J <- matrix(0, n, n)
      h <- 1e-7
      f0 <- self$func(y, t, p)
      for (i in 1:n) {
        y_p <- y
        y_p[i] <- y_p[i] + h
        J[, i] <- (self$func(y_p, t, p) - f0) / h
      }
      return(J)
    },
    grad_theta = function(y, t, p) {
      np <- length(p)
      ny <- length(y)
      J <- matrix(0, ny, np)
      h <- 1e-7
      f0 <- self$func(y, t, p)
      for (i in 1:np) {
        p_p <- p
        p_p[i] <- p_p[i] + h
        J[, i] <- (self$func(y, t, p_p) - f0) / h
      }
      return(J)
    },
    hess_xx = function(y, t, p) {
      n <- length(y)
      H <- array(0, dim = c(n, n, n))
      h <- 1e-5
      J0 <- self$grad_x(y, t, p)
      for (i in 1:n) {
        y_p <- y
        y_p[i] <- y_p[i] + h
        J_p <- self$grad_x(y_p, t, p)
        H[, , i] <- (J_p - J0) / h
      }
      return(H)
    },
    hess_xtheta = function(y, t, p) {
      n <- length(y)
      np <- length(p)
      H <- array(0, dim = c(n, np, n))
      h <- 1e-5
      J0 <- self$grad_theta(y, t, p)
      for (i in 1:n) {
        y_p <- y
        y_p[i] <- y_p[i] + h
        J_p <- self$grad_theta(y_p, t, p)
        H[, , i] <- (J_p - J0) / h
      }
      return(H)
    }
  )
)

#### 3. Inner Solver (Integral Form with Smoothed Data) ----
OdeSystemSolver <- R6Class("OdeSystemSolver",
  public = list(
    func_rhs = NULL, params = NULL, lambda = NULL,
    times_sim = NULL, n_steps = NULL, dt_vec = NULL, n_vars = NULL,
    smoother = NULL, # Instance of DataSmoother
    y_hat_mapped = NULL, # Pre-evaluated smooth curve on grid
    y = NULL, u = NULL, p = NULL,
    initialize = function(func_rhs = NULL, times_sim = NULL, smoother_obj = NULL) {
      if (!is.null(func_rhs)) self$setup(func_rhs, times_sim, smoother_obj)
    },
    setup = function(func_rhs, times_sim, smoother_obj) {
      self$func_rhs <- func_rhs
      self$smoother <- smoother_obj

      # Simulation grid is just the fine grid now (no need to insert obs times)
      self$times_sim <- times_sim
      self$n_steps <- length(times_sim)
      self$dt_vec <- c(diff(times_sim), 0)

      # Pre-evaluate the smoothed target curve on the fine grid
      # Result is N_steps x N_vars
      self$y_hat_mapped <- self$smoother$eval(self$times_sim)
      self$n_vars <- ncol(self$y_hat_mapped)

      self$y <- matrix(0, self$n_steps, self$n_vars)
      self$u <- matrix(0, self$n_steps, self$n_vars)
      self$p <- matrix(0, self$n_steps, self$n_vars)
    },
    # Numerical Jacobian of f with respect to y (df/dy)
    get_jacobian = function(y_vec, t_val) {
      n <- length(y_vec)
      J <- matrix(0, n, n)
      eps <- 1e-7
      if (any(!is.finite(y_vec))) {
        return(J)
      }
      f0 <- self$func_rhs(y_vec, t_val, self$params)
      for (j in 1:n) {
        y_p <- y_vec
        y_p[j] <- y_p[j] + eps
        f_p <- self$func_rhs(y_p, t_val, self$params)
        J[, j] <- (f_p - f0) / eps
      }
      return(J)
    },
    # State equation: dy/dt = f(y, t, params) + u(t)
    solve_state = function(u_mat, y0) {
      y_new <- matrix(0, self$n_steps, self$n_vars)
      y_new[1, ] <- y0
      for (t in 1:(self$n_steps - 1)) {
        y_prev <- y_new[t, ]
        if (any(abs(y_prev) > 1e10)) {
          y_new[(t + 1):self$n_steps, ] <- 1e10
          return(y_new)
        }
        dt <- self$dt_vec[t]
        dy <- self$func_rhs(y_prev, self$times_sim[t], self$params)
        y_new[t + 1, ] <- y_prev + dt * (dy + u_mat[t, ])
      }
      return(y_new)
    },
    # Adjoint equation: dp/dt = -df/dy'p - 2(x - y_hat)
    solve_adjoint = function(y_curr) {
      p_new <- matrix(0, self$n_steps, self$n_vars)
      if (any(y_curr > 1e9)) {
        return(p_new)
      }

      for (t in (self$n_steps - 1):1) {
        dt <- self$dt_vec[t]
        y_now <- y_curr[t, ]

        # Integral Form Forcing: 2 * (x(t) - y_hat(t))
        # y_hat_mapped contains the smoothed curve values
        resid <- y_now - self$y_hat_mapped[t, ]
        forcing <- 2.0 * resid

        J <- self$get_jacobian(y_now, self$times_sim[t])
        grad_prop <- t(J) %*% p_new[t + 1, ]

        p_new[t, ] <- p_new[t + 1, ] + dt * (as.vector(grad_prop) + forcing)
      }
      return(p_new)
    },
    cost_function = function(u_flat, y0) {
      u_mat <- matrix(u_flat, self$n_steps, self$n_vars)
      y_sim <- self$solve_state(u_mat, y0)
      if (any(!is.finite(y_sim)) || any(abs(y_sim) > 1e9)) {
        return(1e20)
      }

      # Integral Cost: Integral( ||x - y_hat||^2 )
      sq_err <- (self$y_hat_mapped - y_sim)^2
      dt_col <- matrix(self$dt_vec, nrow = self$n_steps, ncol = self$n_vars, byrow = FALSE)

      data_term <- sum(sq_err * dt_col)
      reg_term <- sum(u_mat^2 * dt_col)

      return(data_term + self$lambda * reg_term)
    },
    gradient_function = function(u_flat, y0) {
      u_mat <- matrix(u_flat, self$n_steps, self$n_vars)
      y_sim <- self$solve_state(u_mat, y0)
      if (any(y_sim > 1e9)) {
        return(as.vector(2 * self$lambda * u_mat))
      }

      p_sim <- self$solve_adjoint(y_sim)
      dt_col <- matrix(self$dt_vec, nrow = self$n_steps, ncol = self$n_vars, byrow = FALSE)

      grad <- (2 * self$lambda * u_mat + p_sim) * dt_col
      return(as.vector(grad))
    },
    optimize = function(y0, u_init = NULL, max_iter = 100) {
      start_par <- if (!is.null(u_init)) as.vector(u_init) else if (all(self$u == 0)) rnorm(self$n_steps * self$n_vars, sd = 0.01) else as.vector(self$u)
      res <- tryCatch(
        {
          optim(
            par = start_par, fn = self$cost_function, gr = self$gradient_function, y0 = y0,
            method = "BFGS", control = list(maxit = max_iter, trace = 0)
          )
        },
        error = function(e) list(par = start_par, value = 1e20)
      )
      self$u <- matrix(res$par, self$n_steps, self$n_vars)
      self$y <- self$solve_state(self$u, y0)
      self$p <- self$solve_adjoint(self$y)
      return(res)
    }
  )
)

#### 4. Outer Solver (Exact Gradient with Smoothed Data) ----
ExactGradientEstimator <- R6Class("ExactGradientEstimator",
  public = list(
    inner_solver = NULL,
    n_params = NULL,
    deriv_helper = NULL,
    outer_cost_history = NULL,
    initialize = function(inner_solver, n_params) {
      self$inner_solver <- inner_solver
      self$n_params <- n_params
      self$deriv_helper <- NumDeriv$new(inner_solver$func_rhs)
      self$outer_cost_history <- c()
    },
    solve_sensitivity_adjoint = function(y_traj, p_traj) {
      N <- self$inner_solver$n_steps
      D <- self$inner_solver$n_vars
      dt_vec <- self$inner_solver$dt_vec
      times <- self$inner_solver$times_sim
      theta <- self$inner_solver$params

      psi <- matrix(0, nrow = N, ncol = 2 * D)

      mat_weight <- 1 / (2 * self$inner_solver$lambda) * diag(D)

      C_TC <- diag(D)

      for (t in (N - 1):1) {
        dt <- dt_vec[t]
        y_val <- y_traj[t, ]
        p_val <- p_traj[t, ]

        df_x <- self$deriv_helper$grad_x(y_val, times[t], theta)
        df_xx <- self$deriv_helper$hess_xx(y_val, times[t], theta)

        term_xx <- matrix(0, D, D)
        for (d in 1:D) term_xx <- term_xx + df_xx[, , d] * p_val[d]

        lower_left <- -term_xx + C_TC

        M <- rbind(
          cbind(df_x, mat_weight),
          cbind(lower_left, -t(df_x))
        )

        resid <- y_val - self$inner_solver$y_hat_mapped[t, ]
        forcing_x <- 2.0 * resid
        forcing_p <- mat_weight %*% p_val
        forcing <- c(forcing_x, forcing_p)

        d_psi <- forcing - t(M) %*% psi[t + 1, ]
        psi[t, ] <- psi[t + 1, ] - dt * as.vector(d_psi)
      }
      return(psi)
    },
    outer_grad_fn = function(theta, y0) {
      if (!all(self$inner_solver$params == theta)) self$outer_cost_fn(theta, y0)

      psi <- self$solve_sensitivity_adjoint(self$inner_solver$y, self$inner_solver$p)

      grad_accum <- numeric(self$n_params)
      N <- self$inner_solver$n_steps
      D <- self$inner_solver$n_vars

      for (t in 1:(N - 1)) {
        dt <- self$inner_solver$dt_vec[t]
        y_val <- self$inner_solver$y[t, ]
        p_val <- self$inner_solver$p[t, ]
        psi_val <- psi[t, ]

        df_theta <- self$deriv_helper$grad_theta(y_val, self$inner_solver$times_sim[t], theta)

        # --- FIX: CORRECT TENSOR CONTRACTION ---
        # Hessian H[k, m, i] corresponds to d/dx_i ( df_k / dtheta_m )
        # We need sum_k ( p_k * H[k, m, i] )
        df_xtheta <- self$deriv_helper$hess_xtheta(y_val, self$inner_solver$times_sim[t], theta)
        term_xtheta <- matrix(0, nrow = D, ncol = self$n_params)

        for (i in 1:D) { # Loop over rows of the result (x component)
          slice <- df_xtheta[, , i] # Shape [f_comp, param]
          # Contract p with f_comp (1st dim)
          term_xtheta[i, ] <- t(p_val) %*% slice
        }

        bottom_part <- -term_xtheta
        Der_F_param <- rbind(df_theta, bottom_part)

        term <- -1 * (t(psi_val) %*% Der_F_param)
        grad_accum <- grad_accum + as.vector(term) * dt
      }
      return(grad_accum)
    },
    outer_cost_fn = function(theta, y0) {
      self$inner_solver$params <- theta
      self$inner_solver$optimize(y0, max_iter = 500)
      return(self$inner_solver$cost_function(as.vector(self$inner_solver$u), y0))
    },
    optimize = function(init_params, y0, max_iter = 50) {
      cat("Starting Exact Gradient Estimator (Corrected)...\n")
      fn <- function(p) {
        v <- self$outer_cost_fn(p, y0)
        cat(sprintf("Cost: %.4f | P: %s\n", v, paste(round(p, 3), collapse = " ")))
        return(v)
      }
      gr <- function(p) self$outer_grad_fn(p, y0)

      res <- optim(
        par = init_params, fn = fn, gr = gr, method = "Nelder-Mead",
        # lower = rep(1e-3, self$n_params), upper = rep(Inf, self$n_params),
        control = list(maxit = max_iter, trace = 1)
      )
      return(res)
    }
  )
)

#### 5. Test Run ----
run_exact_tracking <- function() {
  lotka_volterra <- function(y, t, p) c(p[1] * y[1] - p[2] * y[1] * y[2], p[3] * y[1] * y[2] - p[4] * y[2])
  true_params <- c(1.5, 0.4, 0.1, 0.4)

  # 1. Generate Raw Data
  times_sim <- seq(0, 5, by = 0.05)
  out <- ode(c(10, 10), times_sim, function(t, y, p) list(lotka_volterra(y, t, p)), true_params)

  # Add noise
  obs_data <- out[, -1] + matrix(rnorm(2 * length(times_sim), sd = 0.3), ncol = 2)
  obs_times <- times_sim

  # 2. Pre-smooth Data
  smoother <- DataSmoother$new(obs_times, obs_data)

  # 3. Setup Inner Solver with Smoother
  inner <- OdeSystemSolver$new(lotka_volterra, times_sim, smoother)

  # Note: Lambda must be tuned to balance the Integral cost
  inner$lambda <- 5e-2
  inner$params <- c(1, 1, 1, 1)

  # 4. Run Estimator
  estimator <- ExactGradientEstimator$new(inner, n_params = 4)
  res <- estimator$optimize(c(1, 0.5, 0.3, 0.5), c(10, 10), max_iter = 100)

  print("True:")
  print(true_params)
  print("Est:")
  print(res$par)

  # 5. Plot
  y_final <- inner$y
  y_target <- inner$y_hat_mapped

  par(mfrow = c(1, 2))
  plot(times_sim, y_target[, 1], type = "l", lty = 2, ylim = range(y_target), main = "Prey (Smoothed vs Fitted)")
  lines(times_sim, y_final[, 1], col = "blue", lwd = 2)
  points(obs_times, obs_data[, 1], col = "grey", pch = 16, cex = 0.5)

  plot(times_sim, y_target[, 2], type = "l", lty = 2, ylim = range(y_target), main = "Predator (Smoothed vs Fitted)")
  lines(times_sim, y_final[, 2], col = "red", lwd = 2)
  points(obs_times, obs_data[, 2], col = "grey", pch = 16, cex = 0.5)
}

run_exact_tracking()
ExactGradientEstimator <- R6Class("ExactGradientEstimator",
  public = list(
    inner_solver = NULL, n_params = NULL, deriv_helper = NULL, outer_cost_history = NULL,
    initialize = function(inner_solver, n_params) {
      self$inner_solver <- inner_solver
      self$n_params <- n_params
      self$deriv_helper <- NumDeriv$new(inner_solver$func_rhs)
      self$outer_cost_history <- c()
    },
    solve_sensitivity_adjoint = function(y_traj, p_traj) {
      N <- self$inner_solver$n_steps
      D <- self$inner_solver$n_vars
      dt_vec <- self$inner_solver$dt_vec
      times <- self$inner_solver$times_sim
      theta <- self$inner_solver$params

      psi <- matrix(0, nrow = N, ncol = 2 * D)

      mat_weight_neg <- -1 / (2 * self$inner_solver$lambda) * diag(D)
      mat_weight_pos <- 1 / (2 * self$inner_solver$lambda) * diag(D)


      C_TC <- diag(D)

      for (t in (N - 1):1) {
        dt <- dt_vec[t]
        y_val <- y_traj[t, ]
        p_val <- p_traj[t, ]

        df_x <- self$deriv_helper$grad_x(y_val, times[t], theta)
        df_xx <- self$deriv_helper$hess_xx(y_val, times[t], theta)

        term_xx <- matrix(0, D, D)
        for (d in 1:D) term_xx <- term_xx + df_xx[d, , ] * p_val[d]

        lower_left <- -term_xx + C_TC

        M <- rbind(
          cbind(df_x, mat_weight_neg),
          cbind(lower_left, -t(df_x))
        )

        resid <- y_val - self$inner_solver$y_hat_mapped[t, ]
        resid[is.na(resid)] <- 0

        forcing_x <- 2.0 * resid
        forcing_p <- mat_weight_pos %*% p_val
        forcing <- c(forcing_x, forcing_p)

        d_psi <- forcing - t(M) %*% psi[t + 1, ]
        psi[t, ] <- psi[t + 1, ] + dt * as.vector(d_psi)
      }
      return(psi)
    },
    outer_grad_fn = function(theta, y0) {
      if (!all(self$inner_solver$params == theta)) self$outer_cost_fn(theta, y0)

      psi <- self$solve_sensitivity_adjoint(self$inner_solver$y, self$inner_solver$p)

      grad_accum <- numeric(self$n_params)
      N <- self$inner_solver$n_steps
      D <- self$inner_solver$n_vars

      for (t in 1:(N - 1)) {
        dt <- self$inner_solver$dt_vec[t]
        y_val <- self$inner_solver$y[t, ]
        p_val <- self$inner_solver$p[t, ]
        psi_val <- psi[t, ]

        df_theta <- self$deriv_helper$grad_theta(y_val, self$inner_solver$times_sim[t], theta)
        df_xtheta <- self$deriv_helper$hess_xtheta(y_val, self$inner_solver$times_sim[t], theta)

        term_xtheta <- matrix(0, nrow = D, ncol = self$n_params)
        for (d in 1:D) term_xtheta <- term_xtheta + df_xtheta[, , d] * p_val[d]
        bottom_part <- -term_xtheta

        Der_F_param <- rbind(df_theta, bottom_part)

        term <- -1 * (t(psi_val) %*% Der_F_param)
        grad_accum <- grad_accum + as.vector(term) * dt
      }
      return(grad_accum)
    },
    outer_cost_fn = function(theta, y0) {
      self$inner_solver$params <- theta
      self$inner_solver$optimize(y0, max_iter = 500)
      return(self$inner_solver$cost_function(as.vector(self$inner_solver$u), y0))
    },
    optimize = function(init_params, y0, max_iter = 50) {
      cat("Starting Exact Gradient Estimator (Corrected Signs)...\n")

      fn <- function(p) {
        v <- self$outer_cost_fn(p, y0)
        cat(sprintf("Cost: %.4f | P: %s\n", v, paste(round(p, 3), collapse = " ")))
        return(v)
      }
      gr <- function(p) self$outer_grad_fn(p, y0)

      res <- optim(
        par = init_params, fn = fn, gr = gr, method = "L-BFGS-B",
        lower = rep(1e-3, self$n_params), upper = rep(Inf, self$n_params),
        control = list(maxit = max_iter, trace = 1)
      )
      return(res)
    }
  )
)

check_gradient_correctness <- function() {
  cat("\n--- RUNNING GRADIENT CHECK ---\n")

  # 1. Setup minimal problem
  lotka_volterra <- function(y, t, p) c(p[1] * y[1] - p[2] * y[1] * y[2], p[3] * y[1] * y[2] - p[4] * y[2])
  times <- seq(0, 2, by = 0.05) # Short horizon for speed
  obs_data <- matrix(10 + rnorm(2 * length(times)), ncol = 2)
  smoother <- DataSmoother$new(times, obs_data)
  inner <- OdeSystemSolver$new(lotka_volterra, times, smoother)
  inner$lambda <- 1.0
  inner$params <- c(1, 1, 1, 1)

  estimator <- ExactGradientEstimator$new(inner, n_params = 4)
  y0 <- c(10, 10)
  test_params <- c(1.0, 0.5, 0.2, 0.5)

  # 2. Compute Analytic Gradient (Your Code)
  cat("Computing Analytic Gradient...\n")
  grad_analytic <- estimator$outer_grad_fn(test_params, y0)

  # 3. Compute Numerical Gradient (Finite Difference)
  cat("Computing Finite Difference Gradient...\n")
  grad_num <- numeric(4)
  eps <- 1e-8

  base_cost <- estimator$outer_cost_fn(test_params, y0)

  for (i in 1:4) {
    p_up <- test_params
    p_up[i] <- p_up[i] + eps
    cost_up <- estimator$outer_cost_fn(p_up, y0)

    # Forward difference: (f(x+h) - f(x))/h
    grad_num[i] <- (cost_up - base_cost) / eps
  }

  # 4. Compare
  cat("\n--- RESULTS ---\n")
  df <- data.frame(
    Param = 1:4,
    Analytic = grad_analytic,
    Numeric = grad_num,
    Diff = abs(grad_analytic - grad_num),
    Rel_Error = abs(grad_analytic - grad_num) / (abs(grad_num) + 1e-8)
  )
  print(df)

  if (mean(df$Rel_Error) < 1e-3) {
    cat("\nSUCCESS: Analytic gradient matches numerical gradient.\n")
  } else {
    cat("\nFAILURE: Gradient mismatch detected.\n")
  }
}

library(ggplot2)
library(dplyr)
library(viridis)

# 1. Setup the Problem (Same as before)
setup_problem <- function() {
  lotka_volterra <- function(y, t, p) c(p[1] * y[1] - p[2] * y[1] * y[2], p[3] * y[1] * y[2] - p[4] * y[2])
  true_params <- c(1.1, 0.4, 0.1, 0.4)
  times <- seq(0, 5, by = 0.05)

  # Generate Data
  out <- ode(c(10, 10), times, function(t, y, p) list(lotka_volterra(y, t, p)), true_params)
  obs_data <- out[, -1] + matrix(rnorm(2 * length(times), sd = 0.2), ncol = 2)

  # Smoother
  smoother <- DataSmoother$new(times, obs_data)

  # Solver
  inner <- OdeSystemSolver$new(lotka_volterra, times, smoother)
  inner$lambda <- 1e1

  return(list(inner = inner, true_params = true_params))
}

# 2. The Scanning Function
scan_landscape <- function(p1_range, p2_range, p_fixed_indices, p_fixed_values) {
  prob <- setup_problem()
  inner <- prob$inner

  # Create Grid
  grid <- expand.grid(p1 = p1_range, p2 = p2_range)
  grid$cost <- NA

  cat(sprintf("Scanning %d points... \n", nrow(grid)))

  for (i in 1:nrow(grid)) {
    # Construct full parameter vector
    params <- numeric(4)
    params[p_fixed_indices] <- p_fixed_values

    # We assume we are scanning the first two free parameters (or whatever is passed)
    # Here we map grid columns to the remaining slots.
    # For simplicity in this script, let's assume we vary Param 1 and Param 2.
    params[1] <- grid$p1[i]
    params[2] <- grid$p2[i]

    # Update Solver
    inner$params <- params

    # Reset Inner Solver to zero (Cold Start) to ensure fairness
    # (Warm start might introduce hysteresis/noise in the plot)
    inner$u <- matrix(0, inner$n_steps, inner$n_vars)

    # Run Inner Optimization
    # CRITICAL: Use high max_iter to reduce "numerical noise"
    res <- inner$optimize(c(10, 10), max_iter = 500)

    # Store Cost
    # Note: If optimization failed significantly, cost might be massive.
    # We log10 the cost for better visualization.
    grid$cost[i] <- res$value

    if (i %% 10 == 0) cat(".")
  }
  cat("\nDone.\n")
  return(grid)
}

# 3. Run and Plot
run_surface_test <- function() {
  # We vary Theta 1 (Growth) and Theta 2 (Predation)
  # True values: 1.1 and 0.4
  p1_seq <- seq(0.8, 1.4, length.out = 15)
  p2_seq <- seq(0.2, 0.6, length.out = 15)

  # Fixed params: Theta 3 and 4 (0.1, 0.4)
  results <- scan_landscape(p1_seq, p2_seq, c(3, 4), c(0.1, 0.4))

  # Plotting
  # We use log(cost) because bad parameters produce EXPLOSIVE costs
  ggplot(results, aes(x = p1, y = p2, z = log10(cost))) +
    geom_contour_filled(bins = 20) +
    geom_point(aes(x = 1.1, y = 0.4), color = "red", shape = 4, size = 5, stroke = 2) + # True Minimum
    scale_fill_viridis_d(option = "magma") +
    labs(
      title = "Optimization Landscape (Log Cost)",
      subtitle = "Red X marks the true parameters",
      x = "Param 1 (Prey Growth)",
      y = "Param 2 (Predation Rate)",
      fill = "Log10 Cost"
    ) +
    theme_minimal()
}
