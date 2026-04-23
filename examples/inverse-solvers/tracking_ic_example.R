source("src/solvers/inverse-solvers/load_inverse_solvers.R")

lv_rhs <- function(y, t, p) c(
  p$alpha * y[1] - p$beta  * y[1] * y[2],
  p$delta * y[1] * y[2] - p$gamma * y[2]
)

euler_solve <- function(rhs, y0, times, params) {
  dt_v <- c(diff(times), 0)
  out  <- matrix(0, length(times), length(y0))
  out[1, ] <- y0
  for (i in seq_len(length(times) - 1))
    out[i + 1, ] <- out[i, ] + dt_v[i] * rhs(out[i, ], times[i], params)
  out
}

times_sim <- round(seq(0, 5, by = 0.01),2)
obs_times <- round(seq(0, 5, by = 0.1),2)

param_scales <- list(alpha = 1.0, beta = 1.0, delta = 1.0, gamma = 1.0,
                       x0 = 1.0, y0 = 1.0)

p_true <- c(alpha = 1.20, beta = 0.45, delta = 0.12, gamma = 0.80,
              x0 = 8.0, y0 = 6.0)

theta_test <- c(alpha = 0.8, beta = 0.60, delta = 0.20, gamma = 0.55,
                  x0 = 6.5, y0 = 4.5)

y_true <- euler_solve(lv_rhs, unlist(p_true[c("x0", "y0")]), times_sim, as.list(p_true))
obs_data <- y_true[which(times_sim %in% obs_times), , drop = FALSE]

model <- ODEModel$new(
    rhs = lv_rhs,
    init_state = function(p) as.numeric(c(p$x0, p$y0)),
    fixed_params = list(),
    param_scales = param_scales
)

tracking <- TrackingInverseSolver$new(
model        = model,
times_sim    = times_sim,
obs_times    = obs_times,
obs_values   = obs_data,
lambda       = 1e1,
inner_method = "gl2",
verbose = T
)


lower_phys <- rep(0,6)
names(lower_phys) <- names(theta_test)

upper_phys <- rep(Inf,6)
names(upper_phys) <- names(theta_test)

result <- tracking$optimize_parameters(
    init_theta_physical = theta_test,
    param_names         = names(theta_test),
    lower_phys = lower_phys, upper_phys = upper_phys)


result
p_true


for(i in 1:ncol(obs_data)){
  plot(obs_times, obs_data[,i])
  lines(times_sim, tracking$last_solver$y[,i])

  plot(times_sim, tracking$last_solver$u[,i], type="l")
}

