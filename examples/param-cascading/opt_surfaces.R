library(ggplot2)
library(gridExtra)
library(reshape2)
library(plotly)

source("src/solvers/parameter_cascading.R")


run_surface_comparison <- function() {
  cat("\n=== Running Surface Comparison (Naive vs Outer) ===\n")

  # 1. Setup ----------------------------------------------------------------
  lotka_volterra <- function(y, t, p) {
    c(p$alpha * y[1] - p$beta  * y[1] * y[2],
      p$delta * y[1] * y[2] - p$gamma * y[2])
  }

  p_true  <- list(alpha = 1.1, beta = 0.4, delta = 0.1, gamma = 0.4)
  y0_true <- c(10, 10)

  times_obs <- seq(0, 40, by = 0.8)
  times_sim <- sort(unique(c(seq(0, 40, 0.05), times_obs)))

  # Generate truth via forward Euler (same integrator as the solver)
  dt_vec  <- c(diff(times_sim), 0)
  y_new   <- matrix(0, length(times_sim), 2)
  y_new[1, ] <- y0_true
  for (t in 1:(length(times_sim) - 1)) {
    y_new[t + 1, ] <- y_new[t, ] + dt_vec[t] *
      lotka_volterra(y_new[t, ], times_sim[t], p_true)
  }

  obs_idx  <- which(times_sim %in% times_obs)
  obs_data <- y_new[obs_idx, ]
  set.seed(123)
  obs_data <- obs_data + matrix(rnorm(length(obs_data), 0, 1),
                                nrow(obs_data), 2)

  # 2. Initialize one GeneralOdeSolver (params updated in-place per grid point)
  # Lambda is the same for naive and outer; it controls how much u(t) can deviate
  # from zero in the outer surface.
  solver <- GeneralOdeSolver$new(
    func_rhs  = lotka_volterra,
    obs_times = times_obs,
    times_sim = times_sim,
    obs_values = obs_data,
    params    = p_true,
    lambda    = 1e3
  )

  # 3. Grid -----------------------------------------------------------------
  alpha_seq <- seq(0.9, 1.3, length.out = 12)
  beta_seq  <- seq(0.3, 0.5, length.out = 12)

  mat_naive <- matrix(NA, length(alpha_seq), length(beta_seq))
  mat_outer <- matrix(NA, length(alpha_seq), length(beta_seq))

  u0_mat <- matrix(0, solver$n_steps, solver$n_vars)
  total  <- length(alpha_seq) * length(beta_seq)
  count  <- 0

  cat(sprintf("Computing surfaces (%d points)...\n", total))

  for (i in seq_along(alpha_seq)) {
    for (j in seq_along(beta_seq)) {
      count <- count + 1
      if (count %% 15 == 0)
        cat(sprintf("  Progress: %d / %d\n", count, total))

      solver$params$alpha <- alpha_seq[i]
      solver$params$beta  <- beta_seq[j]

      # A. Naive surface: SSE when the forcing is fixed to zero.
      #    Bug 5 fix: must call solve_state() with the current parameters
      #    to get a fresh u=0 trajectory — reading solver$y here would return
      #    the state from the previous grid point's optimisation.
      y_naive          <- solver$solve_state(u0_mat, y0_true)
      residual_naive   <- y_naive - solver$observations_mapped
      mat_naive[i, j]  <- sum(residual_naive^2, na.rm = TRUE)

      # B. Outer surface: SSE after the inner optimisation has absorbed
      #    as much misfit as possible via u(t).
      solver$optimize(y0 = y0_true, max_iter = 100)
      residual_outer  <- solver$y - solver$observations_mapped
      mat_outer[i, j] <- sum(residual_outer^2, na.rm = TRUE)
    }
  }
  cat("Done.\n")

  # 4. ggplot2 contour plots ------------------------------------------------
  make_contour_df <- function(mat) {
    df <- melt(mat)
    colnames(df) <- c("Alpha_Idx", "Beta_Idx", "Cost")
    df$Alpha <- alpha_seq[df$Alpha_Idx]
    df$Beta  <- beta_seq[df$Beta_Idx]
    df$logCost <- log1p(df$Cost)
    df
  }

  plot_contour <- function(mat, title) {
    df <- make_contour_df(mat)
    ggplot(df, aes(x = Alpha, y = Beta, z = logCost)) +
      geom_contour_filled(bins = 20) +
      annotate("point", x = 1.1, y = 0.4,
               color = "white", size = 4, shape = 4, stroke = 2) +
      labs(title = title, x = expression(alpha), y = expression(beta),
           fill = "log(1 + SSE)") +
      theme_minimal()
  }

  p1 <- plot_contour(mat_naive,
    "Naive surface  (u = 0)\nSharp minimum — every parameter error is visible")
  p2 <- plot_contour(mat_outer,
    "Outer surface  (u optimised)\nFlatter valley — u(t) absorbs ODE residuals")

  grid.arrange(p1, p2, ncol = 2)

  # 5. Interactive 3-D surfaces (plotly) ------------------------------------
  p3 <- plot_ly() %>%
    add_surface(x = ~beta_seq, y = ~alpha_seq, z = ~log1p(mat_naive),
                colorscale = "Viridis", name = "Naive") %>%
    layout(
      title = "Naive surface  (u = 0, log SSE)",
      scene = list(
        xaxis = list(title = "Beta"),
        yaxis = list(title = "Alpha"),
        zaxis = list(title = "log(1 + SSE)")
      )
    )
  print(p3)

  p4 <- plot_ly() %>%
    add_surface(x = ~beta_seq, y = ~alpha_seq, z = ~log1p(mat_outer),
                colorscale = "Viridis", name = "Outer") %>%
    layout(
      title = "Outer surface  (u optimised, log SSE)",
      scene = list(
        xaxis = list(title = "Beta"),
        yaxis = list(title = "Alpha"),
        zaxis = list(title = "log(1 + SSE)")
      )
    )
  print(p4)

  list(naive = mat_naive, outer = mat_outer,
       alpha = alpha_seq,  beta  = beta_seq)
}
