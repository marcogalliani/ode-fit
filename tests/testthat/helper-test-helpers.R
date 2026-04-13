check_gradient <- function(fn, gr, par, eps = 1e-5, ...) {
  safe_eval <- function(x) {
    out <- tryCatch(fn(x, ...), error = function(e) NA_real_)
    if (length(out) != 1 || !is.finite(out)) return(NA_real_)
    as.numeric(out)
  }

  g_analytic  <- gr(par, ...)
  n           <- length(par)
  g_numerical <- numeric(n)

  for (i in seq_len(n)) {
    p_plus  <- par; p_plus[i]  <- par[i] + eps
    p_minus <- par; p_minus[i] <- par[i] - eps
    f_plus  <- safe_eval(p_plus)
    f_minus <- safe_eval(p_minus)
    g_numerical[i] <- (f_plus - f_minus) / (2 * eps)
  }

  rel_err <- abs(g_analytic - g_numerical) / (abs(g_numerical) + 1e-10)
  norm_a  <- sqrt(sum(g_analytic^2))
  norm_n  <- sqrt(sum(g_numerical^2))
  cos_sim <- if (norm_a == 0 && norm_n == 0) {
    1
  } else if (norm_a > 0 && norm_n > 0) {
    sum(g_analytic * g_numerical) / (norm_a * norm_n)
  } else {
    0
  }

  alpha <- 1.0
  rho   <- 0.5
  c1    <- 1e-4
  f0 <- safe_eval(par)
  f_trial <- safe_eval(par - alpha * g_analytic)
  while (is.finite(f0) && is.finite(f_trial) && f_trial > f0 - c1 * alpha * sum(g_analytic^2)) {
    alpha <- alpha * rho
    if (alpha < 1e-12) break
    f_trial <- safe_eval(par - alpha * g_analytic)
  }
  is_descent <- is.finite(f0) && is.finite(f_trial) && (f_trial < f0)

  list(
    max_rel_error     = max(rel_err),
    cosine_similarity = cos_sim,
    g_analytic        = g_analytic,
    g_numerical       = g_numerical,
    is_descent        = is_descent
  )
}

filter_stable_info <- function(param_set, evaluator) {
  tryCatch({
    val <- evaluator(param_set)
    stable <- is.numeric(val) && length(val) == 1L && is.finite(val)
    list(
      stable = stable,
      value = if (is.numeric(val) && length(val) > 0L) as.numeric(val[[1L]]) else NA_real_,
      reason = if (stable) "" else "non-finite stability value"
    )
  }, error = function(e) {
    list(stable = FALSE, value = NA_real_, reason = conditionMessage(e))
  })
}

filter_stable <- function(param_set, evaluator) {
  isTRUE(filter_stable_info(param_set, evaluator)$stable)
}
