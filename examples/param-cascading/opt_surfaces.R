library(ggplot2)
library(gridExtra)
library(reshape2)
library(plotly)

source("src/solvers/inverse-solvers/load_inverse_solvers.R")
source("examples/ode_models.R")


run_surface_comparison <- function() {
  cat("\n=== Running Surface Comparison (Naive vs Outer) ===\n")

  cfg     <- ODE_CONFIGS$lv
  p_true  <- cfg$params
  y0_true <- unname(cfg$y0)

  times_obs <- seq(0, 40, by = 0.8)
  times_sim <- sort(unique(c(seq(0, 40, 0.05), times_obs)))

  y_new    <- euler_solve(y0_true, times_sim, lv_rhs, p_true)
  obs_idx  <- which(times_sim %in% times_obs)
  obs_data <- y_new[obs_idx, ]
  set.seed(123)
  obs_data <- obs_data + matrix(rnorm(length(obs_data), 0, 1),
                                nrow(obs_data), 2)

  solver <- OdeSystemSolver$new(
    func_rhs  = lv_rhs,
    obs_times = times_obs,
    times_sim = times_sim,
    obs_values = obs_data,
    params    = p_true,
    lambda    = 1e3
  )

  # Grid
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

      y_naive          <- solver$solve_state(u0_mat, y0_true)
      residual_naive   <- y_naive - solver$observations_mapped
      mat_naive[i, j]  <- sum(residual_naive^2, na.rm = TRUE)

      solver$optimize(y0 = y0_true, max_iter = 100)
      residual_outer  <- solver$y - solver$observations_mapped
      mat_outer[i, j] <- sum(residual_outer^2, na.rm = TRUE)
    }
  }
  cat("Done.\n")

  # ggplot2 contour plots
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
    "Naive surface  (u = 0)\nSharp minimum \u2014 every parameter error is visible")
  p2 <- plot_contour(mat_outer,
    "Outer surface  (u optimised)\nFlatter valley \u2014 u(t) absorbs ODE residuals")

  grid.arrange(p1, p2, ncol = 2)

  # Interactive 3-D surfaces (plotly)
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
