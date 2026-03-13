source("src/solvers/general_ode_system_solver.R")
source("src/utils/trace_optimisation.R")

library(plotly)

# Utilities----

## enable debugging for solver methods----
# OdeSystemSolver$debug("gradient_function")

# EXAMPLE 1: Lotka-Volterra (Predator-Prey)----
run_user_example <- function() {
  cat("\n=== Running Lotka-Volterra Example ===\n")
  # Physics
  n_vars <- 2
  lotka_volterra <- function(y, t, p) {
    dx <- p$alpha * y[1] - p$beta * y[1] * y[2]
    dy <- p$delta * y[1] * y[2] - p$gamma * y[2]
    return(c(dx, dy))
  }
  # Parameters
  params <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  # 2. Setup Data
  times_obs <- seq(0, 40, by = 0.2)
  times_sim <- seq(0, 40, 0.01)
  times_sim <- sort(unique((c(times_sim, times_obs))))

  dt_vec <- c(diff(times_sim), 0)

  # Generate synthetic truth with mystery forcing
  solve_ode <- function(y0, times_sim, func_rhs, params, n_vars) {
    n_steps <- length(times_sim)

    y_new <- matrix(0, n_steps, n_vars)
    y_new[1, ] <- y0

    for (t in 1:(n_steps - 1)) {
      y_prev <- y_new[t, ]
      dt <- dt_vec[t]

      dy <- func_rhs(y_prev, times_sim[t], params)
      y_new[t + 1, ] <- y_prev + dt * dy
    }
    return(y_new)
  }


  y0_true <- c(10, 10)
  y_true <- solve_ode(y0_true, times_sim, lotka_volterra, params, n_vars)

  y_true <- y_true[times_sim %in% times_obs, ]

  # Add noise
  set.seed(123)
  obs_data <- y_true + matrix(rnorm(length(y_true), 0, 1.5), nrow(y_true), 2)

  # 3. Run Solver
  solver <- OdeSystemSolver$new(
    func_rhs = lotka_volterra,
    obs_times = times_obs,
    times_sim = times_sim,
    obs_values = obs_data,
    params = params,
    lambda = 1e5
  )
  solver$optimize(y0 = y0_true, max_iter = 200)

  # 4. Plot
  df <- data.frame(
    Time = times_obs,
    Prey_Obs = obs_data[, 1],
    Pred_Obs = obs_data[, 2],
    Prey_Fit = solver$y[times_sim %in% times_obs, 1],
    Pred_Fit = solver$y[times_sim %in% times_obs, 2],
    Prey_Misfit = solver$u[times_sim %in% times_obs, 1],
    Pred_Misfit = solver$u[times_sim %in% times_obs, 2]
  )

  p1 <- ggplot(df, aes(x = Time)) +
    geom_point(aes(y = Prey_Obs), color = "red", alpha = 0.3) +
    geom_line(aes(y = Prey_Fit), color = "red", size = 1) +
    labs(title = "Prey", y = "Pop") +
    theme_minimal()

  p2 <- ggplot(df, aes(x = Time)) +
    geom_point(aes(y = Pred_Obs), color = "blue", alpha = 0.3) +
    geom_line(aes(y = Pred_Fit), color = "blue", size = 1) +
    labs(title = "Predator", y = "Pop") +
    theme_minimal()

  grid.arrange(p1, p2, ncol = 1)

  # Diagnostic: model misfit
  u1 <- ggplot(df, aes(x = Time)) +
    geom_line(aes(y = Prey_Misfit), color = "red", size = 1) +
    labs(title = "Prey", y = "Pop") +
    theme_minimal()

  u2 <- ggplot(df, aes(x = Time)) +
    geom_line(aes(y = Pred_Misfit), color = "blue", size = 1) +
    labs(title = "Predator", y = "Pop") +
    theme_minimal()

  grid.arrange(u1, u2, ncol = 1)

  # Trace optimisation plots
  trace_results <- trace_optimization(solver, y0 = y0_true, max_iter = 50)
  # -> points
  plot(times_obs, obs_data[, 1],
    col = "black", pch = 16, cex = 0.5,
    main = "Optimization Progression", xlab = "Time", ylab = "Y"
  )

  # -> history
  history <- trace_results$history
  n_snaps <- length(history)
  for (i in seq_along(history)) {
    # Calculate a fade from light blue to dark blue
    alpha_val <- seq(0.1, 1, length.out = n_snaps)[i]
    lines(times_sim, history[[i]][, 1],
      col = rgb(0, 0, 1, alpha = alpha_val), lwd = 1
    )
  }
  # -> final fit
  lines(times_sim, solver$solve_state(matrix(trace_results$res$par, solver$n_steps), y0_true)[, 1],
    col = "red", lwd = 2
  )
  legend("topleft",
    legend = c("Data", "Intermediate Steps", "Final Fit"),
    col = c("black", "blue", "red"), lty = c(NA, 1, 1), pch = c(16, NA, NA)
  )
}

# EXAMPLE 2: Sestak-Berggren model----
run_sb_parallel_test <- function() {
  cat("\n=== Running Sestak-Berggren Parallel Test ===\n")

  # physics
  sb_parallel_rhs <- function(x_vec, t, p) {
    # Parallel independent reactions
    eps <- 1e-6
    x_safe <- pmax(pmin(x_vec, p$C0 - eps), eps)

    # Arrhenius equation
    k_T <- p$A * exp(-p$E / (p$R * p$T_vec))
    denom <- p$C0^(p$m + p$n - 1)

    # Rate: -k * (1-alpha)^m * alpha^n
    dx_dt <- -(k_T / denom) * ((p$C0 - x_safe)^p$m) * (x_safe^p$n)

    return(dx_dt)
  }

  T_vec <- c(400, 420, 440)
  L <- length(T_vec)

  params <- list(
    A = 1e6, E = 60000, m = 0.5, n = 1.0,
    C0 = 1.0, R = 8.314, T_vec = T_vec
  )
  # Data
  # Coarse observations (every 0.5s)
  times_obs <- seq(0, 50, by = 0.5)
  # Fine simulation grid (every 0.05s) for stability
  times_sim <- seq(0, 50, by = 0.05)
  times_sim <- sort(unique(c(times_sim, times_obs))) # Ensure obs times are in sim grid

  n_vars <- L
  y0 <- rep(0.99 * params$C0, L)

  solve_ode <- function(y0, times, func, par) {
    dt_v <- c(diff(times), 0)
    n <- length(times)
    out <- matrix(0, n, length(y0))
    out[1, ] <- y0
    for (i in 1:(n - 1)) {
      dy <- func(out[i, ], times[i], par)
      out[i + 1, ] <- out[i, ] + dt_v[i] * dy
    }
    return(out)
  }

  # Generate Truth
  y_true_sim <- solve_ode(y0, times_sim, sb_parallel_rhs, params)

  # Extract observations at specific times
  obs_indices <- which(times_sim %in% times_obs)
  y_true_obs <- y_true_sim[obs_indices, ]

  # Add Noise
  set.seed(42)
  obs_data <- y_true_obs + matrix(rnorm(length(y_true_obs), 0, 0.2), nrow(y_true_obs), L)

  # Optimisation
  solver <- OdeSystemSolver$new(
    func_rhs = sb_parallel_rhs,
    times_sim = times_sim, # The solver runs on this grid
    obs_times = times_obs, # Data exists at these times
    obs_values = obs_data, # The data matrix
    params = params,
    lambda = 0.1 # Regularization
  )

  solver$optimize(y0 = y0, max_iter = 200)

  # plotting
  melt_mat <- function(mat, time_vec, type_name) {
    df <- melt(mat)
    colnames(df) <- c("TimeIdx", "Reactor", "Value")
    df$Time <- time_vec[df$TimeIdx]
    df$Type <- type_name
    return(df)
  }

  # Use solver$observations_mapped for plotting (aligned to times_sim with NAs)
  df_obs <- melt_mat(solver$observations_mapped, times_sim, "Observed")
  df_fit <- melt_mat(solver$y, times_sim, "Fitted")
  df_u <- melt_mat(solver$u, times_sim, "Reconstructed Forcing")

  # Remove NAs from obs for cleaner legend/plotting
  df_obs <- na.omit(df_obs)

  p1 <- ggplot() +
    geom_point(data = df_obs, aes(x = Time, y = Value, color = factor(Reactor)), alpha = 0.4) +
    geom_line(data = df_fit, aes(x = Time, y = Value, color = factor(Reactor), group = Reactor), size = 1) +
    labs(title = "Parallel Sestak-Berggren: Concentration", y = "Conc", color = "Reactor") +
    theme_minimal()

  p2 <- ggplot(data = df_u, aes(x = Time, y = Value, color = factor(Reactor))) +
    geom_line(size = 0.8) +
    labs(title = "Reconstructed Dynamics Misfit u(t)", y = "u(t)") +
    theme_minimal()

  grid.arrange(p1, p2, ncol = 1)
}

# EXAMPLE 3: Missing data test----
run_missing_data_test <- function() {
  cat("\n=== Running Missing Data (NA) Test ===\n")

  # 1. Physics
  lotka_volterra <- function(y, t, p) {
    dx <- p$alpha * y[1] - p$beta * y[1] * y[2]
    dy <- p$delta * y[1] * y[2] - p$gamma * y[2]
    return(c(dx, dy))
  }
  params <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)

  # 2. Setup Data
  times_obs <- seq(0, 40, by = 0.2)
  times_sim <- seq(0, 40, by = 0.05) # Finer grid for smoother derivatives
  times_sim <- sort(unique(c(times_sim, times_obs)))

  n_vars <- 2
  y0_true <- c(10, 10)

  # Synthetic Truth Generator
  solve_ode <- function(y0, times, func, par) {
    dt_v <- c(diff(times), 0)
    out <- matrix(0, length(times), 2)
    out[1, ] <- y0
    for (i in 1:(length(times) - 1)) {
      dy <- func(out[i, ], times[i], par)
      out[i + 1, ] <- out[i, ] + dt_v[i] * dy
    }
    return(out)
  }

  y_true_sim <- solve_ode(y0_true, times_sim, lotka_volterra, params)

  # Extract Obs
  obs_indices <- which(times_sim %in% times_obs)
  obs_data <- y_true_sim[obs_indices, ]

  # Add Noise
  set.seed(123)
  obs_data <- obs_data + matrix(rnorm(length(obs_data), 0, 1.5), nrow(obs_data), 2)

  # 3. Introduce NAs (The Test)
  # Delete data between t=10 and t=25 for the Predator (Column 2)
  # Map time to indices in the OBS grid
  idx_gap <- which(times_obs >= 10 & times_obs <= 25)
  obs_data[idx_gap, 2] <- NA

  cat("Introduced", sum(is.na(obs_data)), "NAs into the dataset.\n")

  # 4. Run Optimization
  solver <- OdeSystemSolver$new(
    func_rhs = lotka_volterra,
    times_sim = times_sim,
    obs_times = times_obs,
    obs_values = obs_data, # Contains the NAs
    params = params,
    lambda = 1
  )

  # Note: The solver's initialize() maps obs_data to the fine grid.
  # The NAs in obs_data are preserved in solver$observations_mapped.

  solver$optimize(y0 = y0_true, max_iter = 200)

  # 5. Plotting
  # We construct a DF on the SIMULATION grid
  df <- data.frame(
    Time = times_sim,
    Prey_Fit = solver$y[, 1],
    Pred_Fit = solver$y[, 2],
    # Map observations to the simulation grid for plotting (aligned)
    Prey_Obs = solver$observations_mapped[, 1],
    Pred_Obs = solver$observations_mapped[, 2]
  )

  p1 <- ggplot(df, aes(x = Time)) +
    # Plot observations (NAs will automatically be skipped by geom_point)
    geom_point(aes(y = Pred_Obs), color = "blue", size = 2, alpha = 0.5) +
    # Plot fitted line (should exist everywhere, bridging the gap)
    geom_line(aes(y = Pred_Fit), color = "blue", size = 1) +
    annotate("rect", xmin = 10, xmax = 25, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "gray") +
    annotate("text", x = 17.5, y = 5, label = "Missing Data Region") +
    labs(
      title = "Predator Reconstruction with Gap (t=10 to 25)",
      subtitle = "Blue Points = Available Data, Line = Physics-Informed Fit",
      y = "Population"
    ) +
    theme_minimal()

  print(p1)
}

# EXAMPLE 4: non-zero forcing----
run_discovery_test <- function(lambda = 0.01) {
  cat("\n=== Running Discovery Test (Missing Physics) ===\n")

  # 1. Physics: Simple Linear System
  k <- 6e-2
  # Solver thinks the physics is just dy/dt = -0.5 * y^2
  simple_rhs <- function(y, t, p) {
    return(-k * y^3)
  }

  times_sim <- seq(0.01, 10, by = 0.01)
  n_steps <- length(times_sim)

  # 2. Generate Truth with an "Unknown" External Force
  # Truth: dy/dt = -0.5*y^2 + sin(t)  <-- Solver doesn't know about sin(t)
  y_true <- numeric(n_steps)
  y_true[1] <- 10

  for (i in 1:(n_steps - 1)) {
    dt <- 0.1
    dy <- -k * y_true[i]^(1/3) + 0.2*sin(times_sim[i]) # The hidden force
    y_true[i + 1] <- y_true[i] + dt * dy
  }

  obs_data <- matrix(y_true + rnorm(n_steps, 0, 0.3), ncol = 1)

  # 3. Solver
  solver <- OdeSystemSolver$new(
    func_rhs = simple_rhs,
    times_sim = times_sim,
    obs_times = times_sim,
    obs_values = obs_data,
    params = list(),
    lambda = lambda
  )

  # 4. Optimization
  solver$optimize(y0 = y_true[1], max_iter = 100)

  # 5. Compare u(t) with the actual hidden sin(t)
  plot_df <- data.frame(
    Time = times_sim,
    Hidden_Force = sin(times_sim),
    Recovered_U = solver$u[, 1],
    Fitted_Y = solver$y[, 1],
    True_Y = y_true
  )

  p1 <- ggplot(plot_df, aes(x = Time)) +
    geom_point(aes(y = True_Y), alpha = 0.3) +
    geom_line(aes(y = Fitted_Y), color = "red") +
    labs(title = "State Fit (Red) vs Truth (Dots)", subtitle = "Optimizer uses u(t) to bridge the gap") +
    theme_minimal()

  p2 <- ggplot(plot_df, aes(x = Time)) +
    geom_line(aes(y = Hidden_Force), linetype = "dashed") +
    geom_line(aes(y = Recovered_U), color = "blue") +
    labs(title = "Force Recovery", y = "Force", subtitle = "Dashed: actual sin(t), Blue: Recovered u(t)") +
    theme_minimal()

  grid.arrange(p1, p2, ncol = 1)


  # Trace optimisation plots
  trace_results <- trace_optimization(solver, y0 = y_true[1], max_iter = 100)
  # -> points
  plot(times_sim, obs_data[, 1],
    col = "black", pch = 16, cex = 0.5,
    main = "Optimization Progression", xlab = "Time", ylab = "Y"
  )

  # -> history
  history <- trace_results$history
  n_snaps <- length(history)
  for (i in seq_along(history)) {
    # Calculate a fade from light blue to dark blue
    alpha_val <- seq(0.1, 1, length.out = n_snaps)[i]
    lines(times_sim, history[[i]][, 1],
      col = rgb(0, 0, 1, alpha = alpha_val), lwd = 1
    )
  }
  # -> final fit
  lines(times_sim, solver$solve_state(matrix(trace_results$res$par, solver$n_steps), 5)[, 1],
    col = "red", lwd = 2
  )
  legend("topleft",
    legend = c("Data", "Intermediate Steps", "Final Fit"),
    col = c("black", "blue", "red"),
    lty = c(NA, 1, 1), pch = c(16, NA, NA)
  )

  generate_optimization_gif(solver, trace_results, y0 = y_true[1])
}

# EXAMPLE 5: opt surface----
run_sse_surface_test <- function() {
  cat("\n=== Generating SSE Surface using OdeSystemSolver ===\n")

  # 1. Setup Common Physics & Data (Same as your example)
  lotka_volterra <- function(y, t, p) {
    dx <- p$alpha * y[1] - p$beta * y[1] * y[2]
    dy <- p$delta * y[1] * y[2] - p$gamma * y[2]
    return(c(dx, dy))
  }

  # Truth Parameters
  p_true <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  y0_true <- c(10, 10)

  # Time and Data Setup
  times_obs <- seq(0, 40, by = 0.2)
  times_sim <- sort(unique(c(seq(0, 40, 0.01), times_obs)))

  # Generate Truth (reuse your logic)
  dt_vec <- c(diff(times_sim), 0)
  y_new <- matrix(0, length(times_sim), 2)
  y_new[1, ] <- y0_true
  for (t in 1:(length(times_sim) - 1)) {
    y_prev <- y_new[t, ]
    dy <- lotka_volterra(y_prev, times_sim[t], p_true)
    y_new[t + 1, ] <- y_prev + dt_vec[t] * dy
  }

  # Extract Obs and add Noise
  obs_data <- y_new[times_sim %in% times_obs, ]
  set.seed(123)
  obs_data <- obs_data + matrix(rnorm(length(obs_data), 0, 1.5), nrow(obs_data), 2)

  # 2. Initialize Solver Once
  solver <- OdeSystemSolver$new(
    func_rhs = lotka_volterra,
    obs_times = times_obs,
    times_sim = times_sim,
    obs_values = obs_data,
    params = p_true,
    lambda = 0 # Lambda doesn't matter when u=0, but good practice
  )

  # 3. Define Grid for Alpha and Beta
  alpha_seq <- seq(0.8, 1.4, length.out = 25)
  beta_seq <- seq(0.2, 0.6, length.out = 25)

  # Matrix to store costs
  cost_surface <- matrix(NA, nrow = length(alpha_seq), ncol = length(beta_seq))

  # 4. Loop using the Solver's Cost Function
  # We use u=0 to see the pure parametric misfit
  u0_flat <- rep(0, solver$n_steps * solver$n_vars)

  cat("Calculating surface points...")
  for (i in seq_along(alpha_seq)) {
    for (j in seq_along(beta_seq)) {
      # Update solver parameters in place
      solver$params$alpha <- alpha_seq[i]
      solver$params$beta <- beta_seq[j]

      # Calculate cost using your exact solver logic
      # This runs solve_state(u=0) -> calculates SSE -> returns cost
      cost_surface[i, j] <- solver$cost_function(u0_flat, y0_true)
    }
  }
  cat(" Done.\n")

  # 5. Visualisation
  # A. Contour Plot with ggplot2
  df_surf <- melt(cost_surface)
  colnames(df_surf) <- c("Alpha_Idx", "Beta_Idx", "Cost")
  df_surf$Alpha <- alpha_seq[df_surf$Alpha_Idx]
  df_surf$Beta <- beta_seq[df_surf$Beta_Idx]

  p1 <- ggplot(df_surf, aes(x = Alpha, y = Beta, z = Cost)) +
    geom_contour_filled(bins = 15) +
    annotate("point", x = 1.1, y = 0.4, color = "red", size = 3, shape = 4) + # Truth
    labs(
      title = "SSE Landscape (Using OdeSystemSolver)",
      subtitle = "Red X = True Parameters",
      fill = "Cost"
    ) +
    theme_minimal()

  print(p1)

  # B. Interactive 3D Surface with Plotly
  p2 <- plot_ly(x = beta_seq, y = alpha_seq, z = cost_surface) %>%
    add_surface() %>%
    layout(
      title = "Cost Function Surface (u=0)",
      scene = list(
        xaxis = list(title = "Beta"),
        yaxis = list(title = "Alpha"),
        zaxis = list(title = "Cost (SSE)")
      )
    )

  print(p2)
}