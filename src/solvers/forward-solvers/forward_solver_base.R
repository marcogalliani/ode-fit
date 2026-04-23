library(R6)

ForwardSolverBase <- R6Class("ForwardSolverBase",
  public = list(
    model = NULL,
    params = NULL,
    lambda = NULL,
    times_sim = NULL,
    n_steps = NULL,
    dt_vec = NULL,
    n_vars = NULL,
    method = NULL,
    verbose = FALSE,
    observations_mapped = NULL,
    y = NULL,
    u = NULL,
    p = NULL,

    log_debug = function(fmt, ...) {
      if (!isTRUE(self$verbose)) return(invisible(NULL))
      message(sprintf(paste0("[", class(self)[1], "] ", fmt), ...))
      invisible(NULL)
    },

    infer_n_vars = function(obs_values) {
      if (is.null(dim(obs_values))) return(1L)
      as.integer(ncol(obs_values))
    },

    normalize_time_grid = function(times_sim, obs_times) {
      obs_times <- round(obs_times, digits = 10)
      sort(unique(round(c(times_sim, obs_times), digits = 10)))
    },

    build_observations_map = function(times_sim, obs_times, obs_values, n_steps, n_vars) {
      obs_mat <- matrix(NA, nrow = n_steps, ncol = n_vars)
      obs_mat[times_sim %in% obs_times, ] <- obs_values
      obs_mat
    },

    sanitize_initial_state = function(y0) {
      if (length(y0) == 1L && is.na(y0)) {
        first_row <- which(!is.na(self$observations_mapped[, 1L]))[1L]
        y0 <- self$observations_mapped[first_row, ]
        y0[is.na(y0)] <- 0
        return(y0)
      }

      if (!any(is.na(y0))) {
        return(y0)
      }

      y0_out <- y0
      for (v in seq_len(self$n_vars)) {
        if (is.na(y0_out[v])) {
          first_v <- which(!is.na(self$observations_mapped[, v]))[1L]
          y0_out[v] <- if (!is.na(first_v)) self$observations_mapped[first_v, v] else 0
        }
      }
      y0_out
    },

    trapezoidal_weights = function() {
      dt_l <- c(0, self$dt_vec[seq_len(self$n_steps - 1L)])
      (dt_l + self$dt_vec) / 2
    },

    initialize_forward_solver = function(model, times_sim, obs_times, obs_values,
                                         params, lambda, method = NULL,
                                         verbose = FALSE) {
      if (is.null(model) || !inherits(model, "ODEModel")) {
        stop("model must be an ODEModel instance")
      }

      self$model <- model
      self$times_sim <- self$normalize_time_grid(times_sim, obs_times)
      self$params <- params
      self$lambda <- lambda
      self$n_steps <- length(self$times_sim)
      self$n_vars <- self$infer_n_vars(obs_values)
      self$dt_vec <- c(diff(self$times_sim), 0)
      self$method <- method
      self$verbose <- isTRUE(verbose)

      self$observations_mapped <- self$build_observations_map(
        times_sim = self$times_sim,
        obs_times = round(obs_times, digits = 10),
        obs_values = obs_values,
        n_steps = self$n_steps,
        n_vars = self$n_vars
      )

      self$y <- matrix(0, self$n_steps, self$n_vars)
      self$u <- matrix(0, self$n_steps, self$n_vars)
      self$p <- matrix(0, self$n_steps, self$n_vars)
    }

  )
)