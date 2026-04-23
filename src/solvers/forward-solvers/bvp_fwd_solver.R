library(R6)

# =============================================================================
# BvpForwardSolver
#
# Solves the inner OCP via the continuous first-order optimality BVP:
#   y' = f(y,t,theta) - p/(2*lambda)
#   p' = -f_y(y,t,theta)^T p - (2/ns) * mask * (y - obs)
# with boundary conditions on y(t0) and p(tN).
# =============================================================================
BvpForwardSolver <- R6Class("BvpForwardSolver",
  inherit = ForwardSolverBase,

  public = list(
    initialize = function(model, times_sim, obs_times, obs_values,
                          params, lambda, verbose = FALSE) {
      self$initialize_forward_solver(
        model = model,
        times_sim = times_sim,
        obs_times = obs_times,
        obs_values = obs_values,
        params = params,
        lambda = lambda,
        method = "bvp",
        verbose = verbose
      )
    },

    optimize_bvp = function(y0, z_init = NULL,
                            max_iter = 50L, tol = 1e-8, verbose = NULL) {
      ns <- self$n_steps
      ny <- self$n_vars
      if (is.null(verbose)) verbose <- self$verbose

      y0 <- self$sanitize_initial_state(y0)

      obs_clean <- self$observations_mapped
      obs_clean[is.na(obs_clean)] <- 0
      obs_mask <- matrix(as.numeric(!is.na(self$observations_mapped)), ns, ny)
      two_lam <- 2 * self$lambda
      two_over_ns <- 2 / ns

      t_keys <- as.character(round(self$times_sim, 10))
      t_lut <- setNames(seq_len(ns), t_keys)
      get_idx <- function(t_val) t_lut[[as.character(round(t_val, 10))]]

      F_rhs_inner <- function(t_val, z) {
        y_t <- z[seq_len(ny)]
        p_t <- z[ny + seq_len(ny)]
        ti <- get_idx(t_val)
        fy <- self$model$rhs(y_t, t_val, self$params)
        Jfy <- self$model$jacobian_state(y_t, t_val, self$params)
        c(
          fy - p_t / two_lam,
          as.vector(-(t(Jfy) %*% p_t)) - two_over_ns * obs_mask[ti, ] * (y_t - obs_clean[ti, ])
        )
      }

      obs_T <- obs_clean[ns, ]
      mask_T <- obs_mask[ns, ]
      bc_inner <- function(z_l, z_r) {
        c(
          z_l[seq_len(ny)] - y0,
          z_r[ny + seq_len(ny)] - two_over_ns * mask_T * (z_r[seq_len(ny)] - obs_T)
        )
      }

      if (is.null(z_init)) {
        euler_dto <- make_dto_solver("euler")
        rhs_euler <- function(y, t) self$model$rhs(y, t, self$params)
        jac_fn <- function(y, t) self$model$jacobian_state(y, t, self$params)

        y_guess <- euler_dto$solve_state(rhs_euler, y0, self$times_sim)$y
        source_fn <- function(t_idx) {
          r <- y_guess[t_idx, ] - self$observations_mapped[t_idx, ]
          r[is.na(r)] <- 0
          (2 / ns) * r
        }
        p_guess <- euler_dto$solve_adjoint(rhs_euler, rep(0, ny), jac_fn, source_fn)$p
        z_init <- cbind(y_guess, p_guess)
      }

      sol <- solve_bvp_colloc(
        F_rhs = F_rhs_inner,
        bc_residual = bc_inner,
        t_grid = self$times_sim,
        z_init = z_init,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose
      )

      y_sol <- sol$z[, seq_len(ny), drop = FALSE]
      p_sol <- sol$z[, ny + seq_len(ny), drop = FALSE]
      self$y <- y_sol
      self$p <- p_sol
      self$u <- -p_sol / two_lam
      sol
    }
  )
)
