#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

.data <- NULL
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("n_obs", "l2_error", "method", "abs_err", "param"))
}

# Ensure relative file paths work regardless of invocation directory.
argv <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", argv, value = TRUE)
if (length(file_arg) == 1L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  setwd(normalizePath(file.path(script_dir, "../..")))
}

plot_consistency_results <- function(
  out_dir = "simulations/inverse-solvers-consistency/results"
) {
  raw_path <- file.path(out_dir, "lv_consistency_raw_results.csv")

  if (!file.exists(raw_path)) {
    stop(sprintf("Missing raw results file: %s", normalizePath(raw_path, mustWork = FALSE)))
  }

  res_raw <- read.csv(raw_path, stringsAsFactors = FALSE)
  res_ok <- res_raw[(res_raw$ok %in% TRUE) & !is.na(res_raw$l2_error), , drop = FALSE]

  p_l2 <- ggplot(
    res_ok,
    aes(x = factor(.data$n_obs), y = .data$l2_error, fill = .data$method)
  ) +
    geom_boxplot(outlier.alpha = 0.35, position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      title = "Lotka-Volterra consistency: total parameter error (replica distribution)",
      x = "Number of observations",
      y = "L2 error across parameters",
      fill = "Method"
    ) +
    theme_minimal(base_size = 12)

  melt_err <- rbind(
    data.frame(method = res_ok$method, n_obs = res_ok$n_obs, param = "alpha", abs_err = abs(res_ok$err_alpha)),
    data.frame(method = res_ok$method, n_obs = res_ok$n_obs, param = "beta", abs_err = abs(res_ok$err_beta)),
    data.frame(method = res_ok$method, n_obs = res_ok$n_obs, param = "delta", abs_err = abs(res_ok$err_delta)),
    data.frame(method = res_ok$method, n_obs = res_ok$n_obs, param = "gamma", abs_err = abs(res_ok$err_gamma))
  )

  p_param <- ggplot(
    melt_err,
    aes(x = factor(.data$n_obs), y = .data$abs_err, fill = .data$method)
  ) +
    geom_boxplot(outlier.alpha = 0.35, position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ param, scales = "free_y") +
    labs(
      title = "Lotka-Volterra consistency: parameter-wise absolute error (replica distribution)",
      x = "Number of observations",
      y = "Absolute error",
      fill = "Method"
    ) +
    theme_minimal(base_size = 12)

  l2_plot_path <- file.path(out_dir, "lv_consistency_l2_error.png")
  param_plot_path <- file.path(out_dir, "lv_consistency_param_errors.png")

  ggsave(filename = l2_plot_path, plot = p_l2, width = 8, height = 5, dpi = 140)
  ggsave(filename = param_plot_path, plot = p_param, width = 9, height = 6, dpi = 140)

  list(
    l2_plot = normalizePath(l2_plot_path),
    param_plot = normalizePath(param_plot_path)
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- "simulations/inverse-solvers-consistency/results"

  if (length(args) >= 1L && nzchar(args[[1L]])) {
    out_dir <- args[[1L]]
  }

  out <- plot_consistency_results(out_dir = out_dir)

  cat("\nPlots generated successfully.\n")
  cat(sprintf("L2 plot: %s\n", out$l2_plot))
  cat(sprintf("Parameter-wise plot: %s\n", out$param_plot))
}

if (sys.nframe() == 0L) {
  main()
}
