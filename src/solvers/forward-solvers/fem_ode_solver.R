# =============================================================================
# src/solvers/fem_ode_solver.R  —  EXPERIMENTAL
#
# Space-time Galerkin (DG-in-time) ODE solver.
#
# The forward ODE  dy/dt = f(y,t,θ) + u  is discretized element-by-element,
# producing a block-bidiagonal system  K y = b(u).  The adjoint is recast
# as  K^T λ = L,  where L is assembled as an ns×nv load matrix before the
# backward pass (exact point evaluation, no quadrature smearing).
#
# Joint state–adjoint structure
# ─────────────────────────────
# Each private forward solver runs the same Newton iteration as the
# corresponding solve_ode integrator and caches per-element data so that
# solve_adjoint_fem can back-substitute without recomputing any Jacobian.
#
#   Forward  : assemble K element-by-element  →  store fwd_elem_mats
#   Backward : transpose each K_e, back-substitute  →  K^T λ = L
#
# What is cached depends on the Galerkin family:
#
#   cG methods (CN, GL1): element DOFs are the two boundary nodes only.
#     K_e = (M_n, C_n)   nv×nv diagonal and sub-diagonal blocks.
#     Adjoint: M_n^T λ_n = L[n] − C_n^T λ_{n+1}  (scalar back-sub).
#
#   dG method (GL2 = dG(1)): each element has TWO interior stage values
#     Y1, Y2 (Gauss points) that are NOT shared with adjacent elements.
#     K_e is a 2nv×2nv system coupling these internal DOFs.
#     Cached: (Y1, Y2, Jf1, Jf2)  →  K_e^T is the 2nv×2nv stage-adjoint
#     system solved once per element during the backward pass.
#
# For a linear ODE K is the same for both passes; for a nonlinear ODE the
# cached linearised K (from the converged forward Newton iterate) is reused.
#
# Method ↔ DG-in-time equivalence:
#   "cn"  ↔ dG(0)/cG(1) trapezoidal      O(h²)
#   "gl1" ↔ cG(1) 1-point Gauss          O(h²)
#   "gl2" ↔ dG(1) 2-point Gauss-Legendre  O(h⁴)
#
# Global K structure (CN, free nodes 2 … ns):
#
#   K = [ M_1          0    0  … ]     M_n = I − h_n/2 · J_{n+1}  (right node)
#       [ C_2    M_2   0    0  … ]     C_n = −(I + h_n/2 · J_n)   (left  node)
#       [  0     C_3   M_3  0  … ]
#       [  …                      ]
#
# K^T is block-upper-bidiagonal; solving K^T λ = L backward from n=ns to 1
# is the exact discrete adjoint of the forward integrator.
# =============================================================================

library(R6)


# -----------------------------------------------------------------------------
# Assemble the global K matrix for the linearised CN step (dense, for testing).
# Returns list(K, b_ic) where K is (ns-1)nv × (ns-1)nv.
# -----------------------------------------------------------------------------
assemble_global_K_cn <- function(y_fwd, times, dt_vec, jac_fn) {
  ns <- nrow(y_fwd); nv <- ncol(y_fwd)
  N  <- ns - 1L
  K  <- matrix(0, N * nv, N * nv)
  b_ic <- rep(0, N * nv)

  for (n in seq_len(N)) {
    h    <- dt_vec[n]
    J_l  <- jac_fn(y_fwd[n,     ], times[n    ])
    J_r  <- jac_fn(y_fwd[n + 1L,], times[n + 1L])
    M_n  <- diag(nv) - (h / 2) * J_r
    C_n  <- -(diag(nv) + (h / 2) * J_l)
    idx  <- (n - 1L) * nv + seq_len(nv)

    K[idx, idx] <- K[idx, idx] + M_n
    if (n > 1L) {
      idx_prev <- (n - 2L) * nv + seq_len(nv)
      K[idx, idx_prev] <- K[idx, idx_prev] + C_n
    } else {
      b_ic[idx] <- b_ic[idx] - C_n %*% y_fwd[1L, ]
    }
  }
  list(K = K, b_ic = b_ic)
}


# =============================================================================
# FemOdeSolver
# =============================================================================
FemOdeSolver <- R6Class("FemOdeSolver",

  private = list(
    fwd_y         = NULL,   # ns×nv trajectory
    fwd_elem_mats = NULL,   # list[ns-1]: element data cached during forward pass

    cache_u            = NULL,
    cache_y            = NULL,
    cache_lambda       = NULL,
    cache_grad_contrib = NULL,

    # ------------------------------------------------------------------
    # Load vector L = E^T r  (exact point evaluation, no quadrature).
    # L[n,] = (2/ns)(y_n − obs_n); NA observations contribute 0.
    # ------------------------------------------------------------------
    assemble_load = function(y_fwd) {
      ns  <- self$n_steps; nv <- self$n_vars
      obs <- self$observations_mapped
      L   <- matrix(0, ns, nv)
      for (j in seq_len(ns)) {
        r        <- y_fwd[j, ] - obs[j, ]
        r[is.na(r)] <- 0
        L[j, ]   <- (2 / ns) * r
      }
      L
    },

    # ------------------------------------------------------------------
    # Forward CN — Newton per element, caches (M_n, C_n).
    #   M_n = I − h/2 J_right  (right node, converged)
    #   C_n = −(I + h/2 J_left) (left  node)
    # ------------------------------------------------------------------
    forward_cn = function(u_mat, y0) {
      ns <- self$n_steps; nv <- self$n_vars
      times <- self$times_sim; dt <- self$dt_vec
      y    <- matrix(0, ns, nv); y[1L, ] <- y0
      elem <- vector("list", ns - 1L)
      for (t in seq_len(ns - 1L)) {
        h  <- dt[t]
        ft <- self$func_rhs(y[t, ], times[t], self$params) + u_mat[t, ]
        yn <- y[t, ] + h * ft
        for (k in seq_len(10L)) {
          fn  <- self$func_rhs(yn, times[t + 1L], self$params) + u_mat[t + 1L, ]
          res <- yn - y[t, ] - (h / 2) * (ft + fn)
          if (max(abs(res)) < 1e-12) break
          J_r <- self$get_jacobian(yn, times[t + 1L])
          yn  <- yn - as.vector(solve(diag(nv) - (h / 2) * J_r, res))
        }
        J_l    <- self$get_jacobian(y[t, ], times[t])
        J_r    <- self$get_jacobian(yn,     times[t + 1L])
        elem[[t]] <- list(
          M = diag(nv) - (h / 2) * J_r,
          C = -(diag(nv) + (h / 2) * J_l),
          h = h
        )
        y[t + 1L, ] <- yn
      }
      list(y = y, elem_mats = elem)
    },

    # ------------------------------------------------------------------
    # Forward GL1 — implicit midpoint, caches midpoint Jacobian Jm.
    #   M_n = I − h/2 Jm,  C_n = −(I + h/2 Jm)  (same midpoint Jm)
    #   ym = (y[t] + y[t+1]) / 2  (holds by construction)
    # ------------------------------------------------------------------
    forward_gl1 = function(u_mat, y0) {
      ns <- self$n_steps; nv <- self$n_vars
      times <- self$times_sim; dt <- self$dt_vec
      y    <- matrix(0, ns, nv); y[1L, ] <- y0
      elem <- vector("list", ns - 1L)
      for (t in seq_len(ns - 1L)) {
        h     <- dt[t]
        t_mid <- (times[t] + times[t + 1L]) / 2
        u_mid <- u_mat[t, ]   # floor interpolation: consistent with u_fn(t_mid) in solve_gl1
        ym    <- y[t, ] + (h / 2) * (self$func_rhs(y[t, ], t_mid, self$params) + u_mid)
        for (k in seq_len(10L)) {
          fn  <- self$func_rhs(ym, t_mid, self$params) + u_mid
          res <- ym - y[t, ] - (h / 2) * fn
          if (max(abs(res)) < 1e-12) break
          Jm  <- self$get_jacobian(ym, t_mid)
          ym  <- ym - as.vector(solve(diag(nv) - (h / 2) * Jm, res))
        }
        Jm        <- self$get_jacobian(ym, t_mid)
        elem[[t]] <- list(
          M  = diag(nv) - (h / 2) * Jm,
          C  = -(diag(nv) + (h / 2) * Jm),
          h  = h
        )
        y[t + 1L, ] <- 2 * ym - y[t, ]
      }
      list(y = y, elem_mats = elem)
    },

    # ------------------------------------------------------------------
    # Forward GL2 (dG(1)) — 2-stage Gauss-Legendre.
    # Y1, Y2 are the TWO interior DOFs of each dG element; they are NOT
    # nodal values and are not shared across elements.  Caching (Y1, Y2,
    # Jf1, Jf2) gives the full 2nv×2nv element stiffness K_e needed for
    # the exact discrete adjoint.
    # Stage forcing: both stages use u_mat[t,] (floor interpolation,
    # consistent with solve_ode's u_fn).
    # ------------------------------------------------------------------
    forward_gl2 = function(u_mat, y0) {
      ns <- self$n_steps; nv <- self$n_vars
      times <- self$times_sim; dt <- self$dt_vec
      y    <- matrix(0, ns, nv); y[1L, ] <- y0
      elem <- vector("list", ns - 1L)
      for (t in seq_len(ns - 1L)) {
        h   <- dt[t]; yt <- y[t, ]
        t1  <- times[t] + .gl2_c1 * h; t2 <- times[t] + .gl2_c2 * h
        us  <- u_mat[t, ]
        K1  <- self$func_rhs(yt, times[t], self$params) + us; K2 <- K1
        for (k in seq_len(15L)) {
          Y1  <- yt + h * (.gl2_A11 * K1 + .gl2_A12 * K2)
          Y2  <- yt + h * (.gl2_A21 * K1 + .gl2_A22 * K2)
          f1  <- self$func_rhs(Y1, t1, self$params) + us
          f2  <- self$func_rhs(Y2, t2, self$params) + us
          G1  <- K1 - f1; G2 <- K2 - f2; res <- c(G1, G2)
          if (max(abs(res)) < 1e-12) break
          Jf1    <- self$get_jacobian(Y1, t1); Jf2 <- self$get_jacobian(Y2, t2)
          Jblock <- rbind(
            cbind(diag(nv) - h * .gl2_A11 * Jf1, -h * .gl2_A12 * Jf1),
            cbind(-h * .gl2_A21 * Jf2, diag(nv) - h * .gl2_A22 * Jf2)
          )
          dK <- as.vector(solve(Jblock, -res))
          K1 <- K1 + dK[seq_len(nv)]; K2 <- K2 + dK[nv + seq_len(nv)]
        }
        Y1  <- yt + h * (.gl2_A11 * K1 + .gl2_A12 * K2)
        Y2  <- yt + h * (.gl2_A21 * K1 + .gl2_A22 * K2)
        Jf1 <- self$get_jacobian(Y1, t1); Jf2 <- self$get_jacobian(Y2, t2)
        elem[[t]] <- list(Y1 = Y1, Y2 = Y2, Jf1 = Jf1, Jf2 = Jf2, h = h)
        y[t + 1L, ] <- yt + h * 0.5 * (K1 + K2)
      }
      list(y = y, elem_mats = elem)
    }
  ),

  public = list(
    func_rhs = NULL, params = NULL, lambda_reg = NULL,
    times_sim = NULL, n_steps = NULL, dt_vec = NULL, n_vars = NULL,
    method = NULL,
    observations_mapped = NULL,
    y = NULL, u = NULL, p = NULL,

    initialize = function(func_rhs, times_sim, obs_times, obs_values,
                          params, lambda, method = "gl2") {
      self$func_rhs  <- func_rhs
      obs_times  <- round(obs_times,  digits = 10)
      times_sim  <- sort(unique(round(c(times_sim, obs_times), digits = 10)))
      self$times_sim  <- times_sim
      self$params     <- params
      self$lambda_reg <- lambda
      self$n_steps    <- length(times_sim)
      self$n_vars     <- ncol(obs_values)
      self$dt_vec     <- c(diff(times_sim), 0)
      self$method     <- method

      self$observations_mapped <- matrix(NA, self$n_steps, self$n_vars)
      self$observations_mapped[times_sim %in% obs_times, ] <- obs_values
      self$y <- matrix(0, self$n_steps, self$n_vars)
      self$u <- matrix(0, self$n_steps, self$n_vars)
      self$p <- matrix(0, self$n_steps, self$n_vars)
    },

    get_jacobian = function(y_vec, t_val) {
      n <- length(y_vec); J <- matrix(0, n, n); eps <- 1e-7
      for (j in seq_len(n)) {
        yp <- y_vec; yp[j] <- yp[j] + eps
        ym <- y_vec; ym[j] <- ym[j] - eps
        J[, j] <- (self$func_rhs(yp, t_val, self$params) -
                   self$func_rhs(ym, t_val, self$params)) / (2 * eps)
      }
      J
    },

    # ------------------------------------------------------------------
    # Joint forward pass: run the per-method Newton integrator, store y
    # and cache element stiffness matrices for the adjoint.
    # ------------------------------------------------------------------
    solve_forward = function(u_mat, y0) {
      res <- switch(self$method,
        cn  = private$forward_cn( u_mat, y0),
        gl1 = private$forward_gl1(u_mat, y0),
        gl2 = private$forward_gl2(u_mat, y0),
        stop("Unknown method: ", self$method)
      )
      private$fwd_y         <- res$y
      private$fwd_elem_mats <- res$elem_mats
      list(y = res$y, aux = NULL)
    },

    # ------------------------------------------------------------------
    # Adjoint solve: K^T λ = L
    #
    # Step 1 — assemble L = E^T r (exact point evaluation, O(ns)).
    # Step 2 — backward substitution through K^T using CACHED element
    #          matrices — no Jacobian recomputed here.
    #
    # CN / GL1 element (step n → n+1):
    #   M_n^T λ_{n+1}  (diagonal)
    #   C_n^T λ_n contribution to λ_{n+1}  (off-diagonal)
    # K^T backward at node n:
    #   M_{n-1}^T λ_n = L[n,] − C_n^T λ_{n+1}
    #
    # GL2: exact 2nv×2nv stage-adjoint system per element.
    # ------------------------------------------------------------------
    solve_adjoint_fem = function() {
      ns     <- self$n_steps; nv <- self$n_vars
      elem   <- private$fwd_elem_mats
      L      <- private$assemble_load(private$fwd_y)
      lambda <- matrix(0, ns, nv)

      if (self$method %in% c("cn", "gl1")) {
        # Terminal: M_{ns-1}^T λ[ns] = L[ns]
        lambda[ns, ] <- as.vector(solve(t(elem[[ns - 1L]]$M), L[ns, ]))

        # Backward: n = ns-1 down to 1
        for (n in seq.int(ns - 1L, 1L)) {
          # Off-diagonal: -C_n^T λ[n+1]  where C_n = elem[[n]]$C
          off   <- as.vector(-t(elem[[n]]$C) %*% lambda[n + 1L, ])
          rhs_n <- L[n, ] + off
          if (n > 1L) {
            # Diagonal: M_{n-1}^T λ[n] = rhs_n
            lambda[n, ] <- as.vector(solve(t(elem[[n - 1L]]$M), rhs_n))
          } else {
            lambda[1L, ] <- rhs_n   # no left element at IC node
          }
        }

        dt_vec <- self$dt_vec
        if (self$method == "cn") {
          dt_l <- c(0, dt_vec[seq_len(ns - 1L)])
          grad_contrib <- matrix(dt_vec / 2, ns, nv) *
                            rbind(lambda[-1L, , drop = FALSE], matrix(0, 1, nv)) +
                          matrix(dt_l  / 2, ns, nv) * lambda
        } else {   # gl1
          grad_contrib <- matrix(dt_vec, ns, nv) *
                            rbind(lambda[-1L, , drop = FALSE], matrix(0, 1, nv))
        }

      } else if (self$method == "gl2") {
        # dG(1): each element has 2nv interior stage DOFs.
        # K_e^T maps the nodal adjoint λ_{n+1} to stage adjoints (μ1, μ2),
        # which then accumulate into λ_n and grad_contrib.
        lam_sum      <- matrix(0, ns, nv)
        lambda[ns, ] <- L[ns, ]
        for (n in seq.int(ns - 1L, 1L)) {
          e   <- elem[[n]]
          h   <- e$h; Jf1 <- e$Jf1; Jf2 <- e$Jf2
          pt1 <- lambda[n + 1L, ]
          Ke_T <- rbind(
            cbind(diag(nv) - h * .gl2_A11 * t(Jf1), -h * .gl2_A21 * t(Jf2)),
            cbind(-h * .gl2_A12 * t(Jf1), diag(nv) - h * .gl2_A22 * t(Jf2))
          )
          lam  <- as.vector(solve(Ke_T, c(h * 0.5 * pt1, h * 0.5 * pt1)))
          lam1 <- lam[seq_len(nv)]; lam2 <- lam[nv + seq_len(nv)]
          lambda[n, ]  <- pt1 + as.vector(t(Jf1) %*% lam1) +
                                 as.vector(t(Jf2) %*% lam2) + L[n, ]
          lam_sum[n, ] <- lam1 + lam2
        }
        grad_contrib <- lam_sum

      } else {
        stop("Unknown method: ", self$method)
      }

      list(lambda = lambda, grad_contrib = grad_contrib)
    },

    solve_forward_adjoint = function(u_mat, y0) {
      fwd <- self$solve_forward(u_mat, y0)
      if (!all(is.finite(fwd$y)))
        return(list(y = fwd$y, lambda = NULL, grad_contrib = NULL, converged = FALSE))
      adj <- self$solve_adjoint_fem()
      list(y = fwd$y, lambda = adj$lambda, grad_contrib = adj$grad_contrib,
           converged = TRUE)
    },

    cost_function = function(u_flat, y0) {
      ns <- self$n_steps; nv <- self$n_vars
      u_mat <- matrix(u_flat, ns, nv)
      sol   <- self$solve_forward_adjoint(u_mat, y0)
      if (!sol$converged || !all(is.finite(sol$y))) return(1e20)

      private$cache_u            <- u_flat
      private$cache_y            <- sol$y
      private$cache_lambda       <- sol$lambda
      private$cache_grad_contrib <- sol$grad_contrib

      sse    <- sum((self$observations_mapped - sol$y)^2, na.rm = TRUE)
      dt_l   <- c(0, self$dt_vec[seq_len(ns - 1L)])
      w_trap <- (dt_l + self$dt_vec) / 2
      reg    <- self$lambda_reg * sum(w_trap * rowSums(u_mat^2))
      (1 / ns) * sse + reg
    },

    gradient_function = function(u_flat, y0) {
      ns <- self$n_steps; nv <- self$n_vars

      cache_valid <- !is.null(private$cache_u) &&
                     length(private$cache_u) == length(u_flat) &&
                     isTRUE(all.equal(private$cache_u, u_flat, tolerance = 0))

      if (!cache_valid) {
        u_mat <- matrix(u_flat, ns, nv)
        sol   <- self$solve_forward_adjoint(u_mat, y0)
        private$cache_u            <- u_flat
        private$cache_y            <- sol$y
        private$cache_lambda       <- sol$lambda
        private$cache_grad_contrib <- sol$grad_contrib
      }

      u_mat  <- matrix(u_flat, ns, nv)
      dt_l   <- c(0, self$dt_vec[seq_len(ns - 1L)])
      w_trap <- matrix((dt_l + self$dt_vec) / 2, ns, nv)
      as.vector(2 * self$lambda_reg * w_trap * u_mat + private$cache_grad_contrib)
    },

    get_load_vector = function() {
      if (is.null(private$fwd_y)) stop("Run solve_forward first.")
      private$assemble_load(private$fwd_y)
    },

    # Expose the global K matrix (CN only, linearised around current y_fwd).
    get_global_K = function() {
      if (self$method != "cn")
        stop("assemble_global_K_cn only available for method='cn'.")
      if (is.null(private$fwd_y)) stop("Run solve_forward first.")
      jac_fn <- function(y, t) self$get_jacobian(y, t)
      assemble_global_K_cn(private$fwd_y, self$times_sim, self$dt_vec, jac_fn)
    },

    optimize = function(y0, max_iter = 100, u_init = NULL,
                        reltol = sqrt(.Machine$double.eps)) {
      if (is.null(u_init)) u_init <- rep(0, self$n_steps * self$n_vars)
      if (length(y0) == 1 && is.na(y0)) {
        first_row <- which(!is.na(self$observations_mapped[, 1]))[1]
        y0 <- self$observations_mapped[first_row, ]; y0[is.na(y0)] <- 0
      } else if (any(is.na(y0))) {
        for (v in seq_len(self$n_vars)) {
          if (is.na(y0[v])) {
            fv <- which(!is.na(self$observations_mapped[, v]))[1]
            y0[v] <- if (!is.na(fv)) self$observations_mapped[fv, v] else 0
          }
        }
      }
      res <- optim(par = u_init,
                   fn  = self$cost_function,
                   gr  = self$gradient_function,
                   y0  = y0, method = "BFGS",
                   control = list(maxit = max_iter, reltol = reltol, trace = 1))
      self$u <- matrix(res$par, self$n_steps, self$n_vars)
      if (!is.null(private$cache_y)) {
        self$y <- private$cache_y; self$p <- private$cache_lambda
      } else {
        sol    <- self$solve_forward_adjoint(self$u, y0)
        self$y <- sol$y; self$p <- sol$lambda
      }
      res
    }
  )
)
