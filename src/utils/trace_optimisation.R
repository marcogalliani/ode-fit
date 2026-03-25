library(ggplot2)
library(gridExtra)

# Plot the outer optimisation trace from a solver that has been run.
#
# Arguments:
#   solver       - a CascadingOdeSolver or TrackingOdeSolver after optimize_parameters()
#   true_params  - named list of ground-truth values, e.g. list(k=0.5, alpha=1.1)
#   cost_field   - name of the cost entry in history (default "sse")
#   cost_label   - y-axis label for the cost plot
#   title        - top-level title for the arranged plot
#
# Returns: invisibly, the list of ggplot objects.
plot_outer_trace <- function(solver, true_params = NULL,
                             cost_field = "sse",
                             cost_label = "SSE (log scale)",
                             title = "Optimisation Progress") {
  history <- solver$history
  if (length(history) == 0) {
    message("No optimisation history found. Run optimize_parameters() first.")
    return(invisible(NULL))
  }

  iters      <- sapply(history, `[[`, "iter")
  costs      <- sapply(history, `[[`, cost_field)
  params_mat <- do.call(rbind, lapply(history, `[[`, "params"))
  param_names <- colnames(params_mat)

  df_cost <- data.frame(iter = iters, cost = costs)
  p_cost <- ggplot(df_cost, aes(x = iter, y = cost)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 1.5) +
    scale_y_log10() +
    labs(title = "Outer Objective", x = "Outer Iteration", y = cost_label) +
    theme_minimal()

  param_plots <- lapply(param_names, function(nm) {
    df_p <- data.frame(iter = iters, value = params_mat[, nm])
    p <- ggplot(df_p, aes(x = iter, y = value)) +
      geom_line(color = "tomato", linewidth = 0.8) +
      geom_point(color = "tomato", size = 1.5) +
      labs(title = nm, x = "Outer Iteration", y = nm) +
      theme_minimal()
    if (!is.null(true_params) && nm %in% names(true_params)) {
      p <- p + geom_hline(
        yintercept = true_params[[nm]], linetype = "dashed", color = "gray40"
      ) +
        annotate("text",
          x = min(iters), y = true_params[[nm]],
          label = paste0("true=", true_params[[nm]]),
          hjust = 0, vjust = -0.5, size = 3, color = "gray40"
        )
    }
    p
  })

  n_params <- length(param_names)
  param_row <- if (n_params == 1) {
    param_plots[[1]]
  } else {
    arrangeGrob(grobs = param_plots, nrow = 1)
  }

  grid.arrange(p_cost, param_row, nrow = 2, top = title)

  invisible(c(list(p_cost), param_plots))
}

plot_param_trace <- function(cascading, true_params = NULL) {
  plot_outer_trace(cascading, true_params,
                   cost_field = "sse", cost_label = "SSE (log scale)",
                   title = "Parameter Cascading \u2014 Optimisation Progress")
}


trace_optimization <- function(solver, y0, max_iter = 100) {
  history_y <- list()
  u_init <- rep(0, solver$n_steps * solver$n_vars)
  iter <- 0
  
  # 1. Gradient Wrapper: Captures the state ONLY on parameter updates
  wrapped_grad <- function(u_flat, y0) {
    iter <<- iter + 1
    
    # Reconstruct matrix and solve state for the snapshot
    u_mat <- matrix(u_flat, solver$n_steps, solver$n_vars)
    current_y <- solver$solve_state(u_mat, y0)
    
    # Store the trajectory
    history_y[[length(history_y) + 1]] <<- current_y
    
    cat("Iteration:", iter, " - Snapshot captured.\n")
    
    # Call the actual gradient logic
    return(solver$gradient_function(u_flat, y0))
  }
  
  # 2. Run the optimization
  # Note: we use the standard cost function but the wrapped gradient
  res <- optim(
    par = u_init, 
    fn = solver$cost_function, 
    gr = wrapped_grad, 
    y0 = y0, 
    method = "BFGS", 
    control = list(maxit = max_iter)
  )
  
  return(list(res = res, history = history_y))
}


library(magick)

# Generates an animated GIF of the inner-solver optimisation trajectory.
#
# Efficiency: O(n) renders via cumulative magick compositing — avoids the
# O(n²) cost of redrawing all ghost lines from scratch for every frame.
#
# Arguments:
#   solver        - a GeneralOdeSolver instance after optimize()
#   trace_results - output of trace_optimization()
#   y0            - initial state vector
#   filename      - output path for the GIF
#   fps           - frames per second
#   n_var         - which state variable to plot (default: 1)
generate_optimization_gif <- function(solver, trace_results, y0,
                                      filename = "optimization.gif",
                                      fps = 5, n_var = 1) {
  history   <- trace_results$history
  n_snaps   <- length(history)
  times_sim <- solver$times_sim
  obs_data  <- solver$observations_mapped
  W <- 800L; H <- 600L

  # Axis limits computed once
  traj_vals <- unlist(lapply(history, function(h) h[, n_var]))
  ylim <- range(c(obs_data[, n_var], traj_vals), na.rm = TRUE)
  xlim <- range(times_sim)

  # Blank plot with matching coordinate system for overlay layers.
  # Must use identical xlim/ylim and default margins so pixels align.
  overlay_plot <- function() {
    plot(xlim, ylim, type = "n", xlab = "", ylab = "",
         axes = FALSE, ann = FALSE, xaxs = "r", yaxs = "r")
  }

  # 1. Base image: axes + observation points (rendered once)
  base_img <- image_graph(width = W, height = H, res = 96)
  plot(times_sim, obs_data[, n_var],
       col = "black", pch = 16, cex = 0.5,
       xlim = xlim, ylim = ylim,
       xlab = "Time", ylab = "Y", xaxs = "r", yaxs = "r")
  dev.off()

  # 2. Build frames via cumulative compositing.
  #    ghost_img accumulates faint past trajectories; only ONE new line is
  #    rendered per iteration instead of redrawing all previous lines.
  ghost_img  <- base_img
  frame_list <- vector("list", n_snaps + 10L)

  for (i in seq_len(n_snaps)) {
    # 2a. After first frame, bake the previous step into the ghost layer
    if (i > 1L) {
      ghost_layer <- image_graph(width = W, height = H, res = 96,
                                 bg = "transparent")
      overlay_plot()
      lines(times_sim, history[[i - 1L]][, n_var],
            col = rgb(0, 0, 1, 0.15), lwd = 1)
      dev.off()
      ghost_img <- image_composite(ghost_img, ghost_layer, operator = "over")
    }

    # 2b. Current step rendered as a separate layer and composited on top
    cur_layer <- image_graph(width = W, height = H, res = 96,
                             bg = "transparent")
    overlay_plot()
    lines(times_sim, history[[i]][, n_var], col = "blue", lwd = 2)
    dev.off()

    frame <- image_composite(ghost_img, cur_layer, operator = "over")
    frame <- image_annotate(frame,
                            paste0("Iteration: ", i),
                            size = 18, location = "+10+10", color = "black")
    frame_list[[i]] <- frame
  }

  # 3. Final convergence frames (red), reusing the last ghost_img
  final_y <- solver$solve_state(
    matrix(trace_results$res$par, solver$n_steps, solver$n_vars), y0
  )[, n_var]

  final_layer <- image_graph(width = W, height = H, res = 96, bg = "transparent")
  overlay_plot()
  lines(times_sim, final_y, col = "red", lwd = 3)
  dev.off()
  final_frame <- image_composite(ghost_img, final_layer, operator = "over")
  final_frame <- image_annotate(final_frame, "Converged",
                                size = 18, location = "+10+10", color = "darkred")
  for (k in seq_len(10L)) frame_list[[n_snaps + k]] <- final_frame

  # 4. Animate and write
  cat("Animating...\n")
  animation <- image_animate(image_join(frame_list), fps = fps)
  image_write(animation, filename)
  cat("GIF saved as:", filename, "\n")
}

