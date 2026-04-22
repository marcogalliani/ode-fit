#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# Ensure relative source() calls work regardless of invocation directory.
argv <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", argv, value = TRUE)
if (length(file_arg) == 1L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  setwd(normalizePath(file.path(script_dir, "../..")))
}

source("src/solvers/inverse-solvers/load_inverse_solvers.R")

lv_rhs <- function(y, t, p) c(
  p$alpha * y[1] - p$beta * y[1] * y[2],
  p$delta * y[1] * y[2] - p$gamma * y[2]
)

euler_solve <- function(rhs, y0, times, params) {
  dt_v <- c(diff(times), 0)
  out <- matrix(0, nrow = length(times), ncol = length(y0))
  out[1, ] <- y0
  for (i in seq_len(length(times) - 1L)) {
    out[i + 1L, ] <- out[i, ] + dt_v[i] * rhs(out[i, ], times[i], params)
  }
  out
}

make_obs_index <- function(n_times, n_obs) {
  idx <- round(seq(1, n_times, length.out = n_obs))
  sort(unique(as.integer(idx)))
}

safe_fit_solver <- function(method,
                            model,
                            times_sim,
                            obs_times,
                            obs_values,
                            init_theta,
                            param_names,
                            lower_phys,
                            upper_phys,
                            lambda,
                            inner_method,
                            inner_max_iter) {
  t0 <- proc.time()[["elapsed"]]

  solver <- tryCatch(
    {
      if (identical(method, "tracking")) {
        TrackingOdeSolver$new(
          model = model,
          times_sim = times_sim,
          obs_times = obs_times,
          obs_values = obs_values,
          lambda = lambda,
          inner_method = inner_method,
          inner_max_iter = inner_max_iter,
          verbose = FALSE
        )
      } else {
        CascadingInverseSolver$new(
          model = model,
          times_sim = times_sim,
          obs_times = obs_times,
          obs_values = obs_values,
          lambda = lambda,
          inner_method = inner_method,
          inner_max_iter = inner_max_iter,
          verbose = FALSE
        )
      }
    },
    error = function(e) e
  )

  if (inherits(solver, "error")) {
    return(list(
      ok = FALSE,
      est = rep(NA_real_, length(param_names)),
      value = NA_real_,
      elapsed = proc.time()[["elapsed"]] - t0,
      error = conditionMessage(solver)
    ))
  }

  fit <- tryCatch(
    {
      if (identical(method, "tracking")) {
        solver$optimize_parameters(
          init_theta_physical = init_theta,
          param_names = param_names,
          lower_phys = lower_phys,
          upper_phys = upper_phys,
          gradient_mode = "adjoint"
        )
      } else {
        solver$optimize_parameters(
          init_theta_physical = init_theta,
          param_names = param_names,
          lower_phys = lower_phys,
          upper_phys = upper_phys
        )
      }
    },
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      ok = FALSE,
      est = rep(NA_real_, length(param_names)),
      value = NA_real_,
      elapsed = proc.time()[["elapsed"]] - t0,
      error = conditionMessage(fit)
    ))
  }

  names(fit) <- param_names
  last_hist <- if (length(solver$history) > 0L) solver$history[[length(solver$history)]] else NULL
  final_value <- NA_real_
  if (!is.null(last_hist)) {
    if (!is.null(last_hist$cost)) {
      final_value <- as.numeric(last_hist$cost)
    } else if (!is.null(last_hist$sse)) {
      final_value <- as.numeric(last_hist$sse)
    }
  }

  list(
    ok = TRUE,
    est = as.numeric(fit[param_names]),
    value = final_value,
    elapsed = proc.time()[["elapsed"]] - t0,
    error = ""
  )
}

run_consistency_study <- function(n_obs_grid = c(6, 11, 16, 21, 31, 41, 51),
                                  n_rep = 30,
                                  seed = 123,
                                  noise_frac = 0.05,
                                  out_dir = "simulations/inverse-solvers-consistency/results",
                                  inner_method = "cn",
                                  inner_max_iter = 200,
                                  show_progress = TRUE) {
  set.seed(seed)

  times_sim <- seq(0, 3, by = 0.01)
  n_times <- length(times_sim)
  n_obs_grid <- sort(unique(as.integer(n_obs_grid[n_obs_grid >= 2 & n_obs_grid <= n_times])))
  
  if (length(n_obs_grid) == 0L) {
    stop("n_obs_grid must contain values in [2, length(times_sim)]")
  }

  p_true <- c(alpha = 1.20, beta = 0.45, delta = 0.12, gamma = 0.80, x0 = 10.0, y0 = 8.0)
  theta_true <- p_true[c("alpha", "beta", "delta", "gamma")]
  theta_init <- c(alpha = 0.5, beta = 0.8, delta = 0.4, gamma = 0.55)
  param_names <- names(theta_true)

  lower_phys <- c(alpha = 0, beta = 0, delta = 0, gamma = 0)
  upper_phys <- c(alpha = Inf, beta = Inf, delta = Inf, gamma = Inf)

  y0_true <- as.numeric(c(p_true[["x0"]], p_true[["y0"]]))
  y_true <- euler_solve(lv_rhs, y0_true, times_sim, as.list(p_true))
  noise_sd <- pmax(noise_frac * apply(abs(y_true), 2L, max), 1e-6)

  model <- ODEModel$new(
    rhs = lv_rhs,
    init_state = function(p) as.numeric(c(p$x0, p$y0)),
    fixed_params = list(x0 = p_true[["x0"]], y0 = p_true[["y0"]]),
    param_scales = list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0)
  )

  methods <- c("tracking", "cascading")
  lambda_map <- c(tracking = 1e2, cascading = 1e2)

  rows <- vector("list", length(methods) * length(n_obs_grid) * n_rep)
  row_id <- 1L
  total_reps <- length(n_obs_grid) * n_rep
  rep_counter <- 0L
  pb <- NULL

  if (isTRUE(show_progress) && total_reps > 0L) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = total_reps,
      style = 3,
      initial = 0,
      char = "="
    )
    on.exit(close(pb), add = TRUE)
  }

  for (n_obs in n_obs_grid) {
    idx <- make_obs_index(n_times, n_obs)
    obs_times <- times_sim[idx]
    obs_true <- y_true[idx, , drop = FALSE]

    for (rep_id in seq_len(n_rep)) {
      noise_mat <- matrix(
        rnorm(length(obs_true), mean = 0, sd = rep(noise_sd, each = nrow(obs_true))),
        nrow = nrow(obs_true),
        ncol = ncol(obs_true)
      )
      obs_values <- obs_true + noise_mat

      for (method in methods) {
        fit <- safe_fit_solver(
          method = method,
          model = model,
          times_sim = times_sim,
          obs_times = obs_times,
          obs_values = obs_values,
          init_theta = theta_init,
          param_names = param_names,
          lower_phys = lower_phys,
          upper_phys = upper_phys,
          lambda = unname(lambda_map[method]),
          inner_method = inner_method,
          inner_max_iter = inner_max_iter
        )

        est <- fit$est
        names(est) <- param_names
        err <- est - theta_true

        rows[[row_id]] <- data.frame(
          method = method,
          n_obs = n_obs,
          rep = rep_id,
          ok = fit$ok,
          elapsed_sec = fit$elapsed,
          objective = fit$value,
          alpha = est[["alpha"]],
          beta = est[["beta"]],
          delta = est[["delta"]],
          gamma = est[["gamma"]],
          err_alpha = err[["alpha"]],
          err_beta = err[["beta"]],
          err_delta = err[["delta"]],
          err_gamma = err[["gamma"]],
          l2_error = sqrt(sum(err^2)),
          error_msg = fit$error,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }

      rep_counter <- rep_counter + 1L
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, rep_counter)
      }
    }
  }

  res_raw <- do.call(rbind, rows)
  res_raw <- res_raw[seq_len(row_id - 1L), , drop = FALSE]

  agg_fun <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))

  sum_l2 <- aggregate(l2_error ~ method + n_obs, data = res_raw, FUN = agg_fun)
  sum_l2 <- cbind(sum_l2[c("method", "n_obs")], as.data.frame(sum_l2$l2_error))
  names(sum_l2)[3:4] <- c("l2_mean", "l2_sd")

  melt_err <- rbind(
    data.frame(method = res_raw$method, n_obs = res_raw$n_obs, param = "alpha", abs_err = abs(res_raw$err_alpha)),
    data.frame(method = res_raw$method, n_obs = res_raw$n_obs, param = "beta", abs_err = abs(res_raw$err_beta)),
    data.frame(method = res_raw$method, n_obs = res_raw$n_obs, param = "delta", abs_err = abs(res_raw$err_delta)),
    data.frame(method = res_raw$method, n_obs = res_raw$n_obs, param = "gamma", abs_err = abs(res_raw$err_gamma))
  )

  sum_param <- aggregate(abs_err ~ method + n_obs + param, data = melt_err, FUN = agg_fun)
  sum_param <- cbind(sum_param[c("method", "n_obs", "param")], as.data.frame(sum_param$abs_err))
  names(sum_param)[4:5] <- c("abs_err_mean", "abs_err_sd")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(res_raw, file.path(out_dir, "lv_consistency_raw_results.csv"), row.names = FALSE)
  write.csv(sum_l2, file.path(out_dir, "lv_consistency_summary_l2.csv"), row.names = FALSE)
  write.csv(sum_param, file.path(out_dir, "lv_consistency_summary_params.csv"), row.names = FALSE)

  p_l2 <- ggplot(sum_l2, aes(x = n_obs, y = l2_mean, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = pmax(l2_mean - l2_sd, 0), ymax = l2_mean + l2_sd), width = 0.6) +
    labs(
      title = "Lotka-Volterra consistency: total parameter error",
      x = "Number of observations",
      y = "Mean L2 error across parameters",
      color = "Method"
    ) +
    theme_minimal(base_size = 12)

  p_param <- ggplot(sum_param, aes(x = n_obs, y = abs_err_mean, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = pmax(abs_err_mean - abs_err_sd, 0), ymax = abs_err_mean + abs_err_sd), width = 0.6) +
    facet_wrap(~ param, scales = "free_y") +
    labs(
      title = "Lotka-Volterra consistency: parameter-wise absolute error",
      x = "Number of observations",
      y = "Mean absolute error",
      color = "Method"
    ) +
    theme_minimal(base_size = 12)

  ggsave(filename = file.path(out_dir, "lv_consistency_l2_error.png"), plot = p_l2,
         width = 8, height = 5, dpi = 140)
  ggsave(filename = file.path(out_dir, "lv_consistency_param_errors.png"), plot = p_param,
         width = 9, height = 6, dpi = 140)

  fail_tab <- aggregate(ok ~ method + n_obs, data = res_raw, FUN = function(x) sum(!x))
  names(fail_tab)[3] <- "n_fail"

  list(
    raw = res_raw,
    summary_l2 = sum_l2,
    summary_param = sum_param,
    failures = fail_tab,
    out_dir = normalizePath(out_dir)
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  n_rep <- 3
  n_obs_grid <- c(51,101,256)
  inner_max_iter <- 100
  if (length(args) >= 1L) {
    maybe_rep <- suppressWarnings(as.integer(args[[1L]]))
    if (!is.na(maybe_rep) && maybe_rep > 0L) {
      n_rep <- maybe_rep
    }
  }

  if (length(args) >= 2L) {
    grid_tokens <- strsplit(args[[2L]], ",", fixed = TRUE)[[1L]]
    grid_vals <- suppressWarnings(as.integer(trimws(grid_tokens)))
    grid_vals <- grid_vals[!is.na(grid_vals)]
    if (length(grid_vals) > 0L) {
      n_obs_grid <- grid_vals
    }
  }

  if (length(args) >= 3L) {
    maybe_inner <- suppressWarnings(as.integer(args[[3L]]))
    if (!is.na(maybe_inner) && maybe_inner > 0L) {
      inner_max_iter <- maybe_inner
    }
  }

  out <- run_consistency_study(
    n_rep = n_rep,
    n_obs_grid = n_obs_grid,
    inner_max_iter = inner_max_iter
  )

  cat("\nConsistency study completed.\n")
  cat(sprintf("Output directory: %s\n", out$out_dir))
  cat("\nFailures by method and sample size:\n")
  print(out$failures)

  cat("\nMean L2 error summary:\n")
  print(out$summary_l2)
}

if (sys.nframe() == 0L) {
  main()
}
