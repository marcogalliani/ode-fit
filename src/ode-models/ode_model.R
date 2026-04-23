library(R6)

.fd_state_jacobian <- function(rhs_fn, y_vec, t_val, params, eps = 1e-7) {
  n <- length(y_vec)
  J <- matrix(0, n, n)

  for (j in seq_len(n)) {
    y_p <- y_vec
    y_m <- y_vec
    y_p[j] <- y_p[j] + eps
    y_m[j] <- y_m[j] - eps

    J[, j] <- (rhs_fn(y_p, t_val, params) - rhs_fn(y_m, t_val, params)) / (2 * eps)
  }

  J
}

.fd_param_jacobian <- function(rhs_fn, y_vec, t_val, params, param_names = NULL, eps = 1e-7) {
  if (is.null(param_names)) {
    param_names <- names(params)
  }
  if (is.null(param_names) || length(param_names) == 0L) {
    stop("param_names must be non-empty for parameter Jacobian computation")
  }

  nv <- length(y_vec)
  np <- length(param_names)
  J <- matrix(0, nv, np)

  for (j in seq_len(np)) {
    nm <- param_names[j]
    if (is.null(params[[nm]])) {
      stop(sprintf("Parameter '%s' not found in params", nm))
    }

    dth <- eps * max(abs(params[[nm]]), 1)
    p_p <- params
    p_m <- params
    p_p[[nm]] <- params[[nm]] + dth
    p_m[[nm]] <- params[[nm]] - dth

    J[, j] <- (rhs_fn(y_vec, t_val, p_p) - rhs_fn(y_vec, t_val, p_m)) / (2 * dth)
  }

  colnames(J) <- param_names
  J
}

ODEModel <- R6Class("ODEModel",
  public = list(
    rhs_fun = NULL,
    jacobian_state_fun = NULL,
    jacobian_param_fun = NULL,
    init_state_fun = NULL,
    fixed_params = NULL,
    param_scales = NULL,

    initialize = function(rhs,
                          init_state = NULL,
                          fixed_params = list(),
                          param_scales = list(),
                          jacobian_state = NULL,
                          jacobian_param = NULL) {
      if (!is.function(rhs)) {
        stop("rhs must be a function with signature rhs(y, t, params)")
      }
      if (!is.null(init_state) && !is.function(init_state)) {
        stop("init_state must be NULL or a function with signature init_state(params)")
      }
      if (!is.null(jacobian_state) && !is.function(jacobian_state)) {
        stop("jacobian_state must be NULL or a function")
      }
      if (!is.null(jacobian_param) && !is.function(jacobian_param)) {
        stop("jacobian_param must be NULL or a function")
      }

      self$rhs_fun <- rhs
      self$jacobian_state_fun <- jacobian_state
      self$jacobian_param_fun <- jacobian_param
      self$init_state_fun <- init_state
      self$fixed_params <- as.list(fixed_params)
      self$param_scales <- as.list(param_scales)
    },

    merge_params = function(params = list()) {
      if (is.null(params)) {
        params <- list()
      }
      params <- as.list(params)
      merged <- self$fixed_params

      if (length(params) == 0L) {
        return(merged)
      }

      nms <- names(params)
      if (is.null(nms) || any(!nzchar(nms))) {
        stop("params must be a named list/vector")
      }

      for (nm in nms) {
        merged[[nm]] <- params[[nm]]
      }

      merged
    },

    rhs = function(y, t, params = list()) {
      self$rhs_fun(y, t, self$merge_params(params))
    },

    jacobian_state = function(y, t, params = list(), eps = 1e-7) {
      merged <- self$merge_params(params)

      if (!is.null(self$jacobian_state_fun)) {
        return(self$jacobian_state_fun(y, t, merged))
      }

      .fd_state_jacobian(self$rhs_fun, y, t, merged, eps = eps)
    },

    jacobian_param = function(y, t, params = list(), param_names = NULL, eps = 1e-7) {
      merged <- self$merge_params(params)

      if (!is.null(self$jacobian_param_fun)) {
        return(self$jacobian_param_fun(y, t, merged, param_names))
      }

      .fd_param_jacobian(
        rhs_fn = self$rhs_fun,
        y_vec = y,
        t_val = t,
        params = merged,
        param_names = param_names,
        eps = eps
      )
    },

    init_state = function(params = list()) {
      if (is.null(self$init_state_fun)) {
        stop("init_state_fun is not defined for this model")
      }
      self$init_state_fun(self$merge_params(params))
    },

    init_state_jacobian_fd = function(params_phys, param_names, eps = 1e-7) {
      if (is.null(param_names) || length(param_names) == 0L) {
        stop("param_names must be non-empty for initial-state Jacobian computation")
      }

      y0 <- as.numeric(self$init_state(params_phys))
      nv <- length(y0)
      np <- length(param_names)
      J <- matrix(0, nv, np)

      for (j in seq_len(np)) {
        nm <- param_names[j]
        dth <- eps * max(abs(params_phys[[nm]]), 1)
        p_p <- params_phys
        p_m <- params_phys
        p_p[[nm]] <- params_phys[[nm]] + dth
        p_m[[nm]] <- params_phys[[nm]] - dth
        y0_p <- as.numeric(self$init_state(p_p))
        y0_m <- as.numeric(self$init_state(p_m))
        J[, j] <- (y0_p - y0_m) / (2 * dth)
      }

      J
    },

    get_scales = function(param_names) {
      scales <- unlist(self$param_scales[param_names], use.names = TRUE)
      if (length(scales) != length(param_names) || any(is.na(scales))) {
        missing <- param_names[is.na(scales)]
        stop(sprintf("Missing scale for parameter(s): %s", paste(missing, collapse = ", ")))
      }
      as.numeric(scales)
    }
  )
)
