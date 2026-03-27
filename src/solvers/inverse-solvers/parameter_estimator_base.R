library(R6)

ParameterEstimatorBase <- R6Class("ParameterEstimatorBase",
  public = list(
    inner_solver_class = NULL,
    times_sim = NULL, obs_times = NULL, obs_values = NULL,
    func_rhs = NULL, fixed_params = NULL, lambda = NULL,
    param_scales = NULL,
    init_state = NULL,
    n_vars = NULL,
    inner_max_iter = NULL, inner_reltol = NULL,

    last_theta = NULL,
    last_solver = NULL,
    last_u = NULL,
    history = NULL,

    initialize_estimator = function(func_rhs, times_sim, obs_times, obs_values,
                                    fixed_params, lambda, param_scales,
                                    init_state,
                                    inner_max_iter = 200,
                                    inner_reltol = sqrt(.Machine$double.eps)) {
      self$func_rhs <- func_rhs
      self$times_sim <- times_sim
      self$obs_times <- obs_times
      self$obs_values <- obs_values
      self$n_vars <- self$infer_n_vars(obs_values)
      if (is.null(init_state) || !is.function(init_state)) {
        stop("init_state must be provided as a function(p)")
      }
      self$fixed_params <- fixed_params
      self$lambda <- lambda
      self$inner_solver_class <- OdeSystemSolver
      self$param_scales <- param_scales
      self$init_state <- init_state
      self$inner_max_iter <- inner_max_iter
      self$inner_reltol <- inner_reltol
      self$history <- list()
    },

    infer_n_vars = function(obs_values) {
      if (is.null(dim(obs_values))) return(1L)
      as.integer(ncol(obs_values))
    },

    get_scales_vector = function(param_names) {
      scales <- unlist(self$param_scales[param_names], use.names = TRUE)
      if (length(scales) != length(param_names) || any(is.na(scales))) {
        missing <- param_names[is.na(scales)]
        stop(sprintf("Missing scale for parameter(s): %s", paste(missing, collapse = ", ")))
      }
      as.numeric(scales)
    },

    prepare_theta_normalized = function(param_names,
                                        init_theta_physical,
                                        lower_phys = NULL,
                                        upper_phys = NULL) {
      theta_names <- names(init_theta_physical)
      if (is.null(theta_names)) {
        stop("init_theta_physical must be a named vector")
      }
      missing <- setdiff(param_names, theta_names)
      if (length(missing) > 0L) {
        stop(sprintf("Missing initial value for parameter(s): %s", paste(missing, collapse = ", ")))
      }

      scales <- self$get_scales_vector(param_names)
      lower <- rep(-Inf, length(param_names))
      upper <- rep(Inf, length(param_names))

      if (!is.null(lower_phys)) {
        if (is.null(names(lower_phys))) stop("lower_phys must be a named vector")
        idx <- match(param_names, names(lower_phys))
        sel <- !is.na(idx)
        lower[sel] <- lower_phys[idx[sel]]
      }
      if (!is.null(upper_phys)) {
        if (is.null(names(upper_phys))) stop("upper_phys must be a named vector")
        idx <- match(param_names, names(upper_phys))
        sel <- !is.na(idx)
        upper[sel] <- upper_phys[idx[sel]]
      }

      init_norm <- as.numeric(init_theta_physical[param_names] / scales)
      names(init_norm) <- param_names
      lower_norm <- as.numeric(lower / scales)
      names(lower_norm) <- param_names
      upper_norm <- as.numeric(upper / scales)
      names(upper_norm) <- param_names

      list(scales = scales, init = init_norm, lower = lower_norm, upper = upper_norm)
    },

    unpack_physical = function(theta_norm, param_names) {
      curr <- self$fixed_params
      for (i in seq_along(param_names)) {
        name <- param_names[i]
        scale <- self$param_scales[[name]]
        if (is.null(scale)) stop(sprintf("Missing scale for parameter '%s'", name))
        curr[[name]] <- theta_norm[i] * scale
      }
      curr
    },

    eval_init_state = function(params_phys) {
      y0 <- as.numeric(self$init_state(params_phys))
      if (length(y0) != self$n_vars) {
        stop(sprintf("init_state must return a vector of length %d", self$n_vars))
      }
      if (!all(is.finite(y0))) {
        stop("init_state returned non-finite values")
      }
      y0
    },

    init_state_jacobian_fd = function(params_phys, param_names, eps = 1e-7) {
      np <- length(param_names)
      J <- matrix(0, self$n_vars, np)

      for (j in seq_len(np)) {
        nm <- param_names[j]
        dth <- eps * max(abs(params_phys[[nm]]), 1)
        p_p <- params_phys
        p_m <- params_phys
        p_p[[nm]] <- params_phys[[nm]] + dth
        p_m[[nm]] <- params_phys[[nm]] - dth
        y0_p <- self$eval_init_state(p_p)
        y0_m <- self$eval_init_state(p_m)
        J[, j] <- (y0_p - y0_m) / (2 * dth)
      }

      J
    },

    get_param_jacobian = function(y, t, p_phys, param_names, eps = 1e-7) {
      nv  <- length(y)
      np  <- length(param_names)
      J   <- matrix(0, nv, np)
      for (j in seq_len(np)) {
        dth <- eps * max(abs(p_phys[[param_names[j]]]), 1)
        p_p <- p_phys
        p_m <- p_phys
        p_p[[param_names[j]]] <- p_phys[[param_names[j]]] + dth
        p_m[[param_names[j]]] <- p_phys[[param_names[j]]] - dth
        J[, j] <- (self$func_rhs(y, t, p_p) - self$func_rhs(y, t, p_m)) / (2 * dth)
      }
      J
    }
  )
)