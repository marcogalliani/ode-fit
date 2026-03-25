# =============================================================================
# src/solvers/ode_solvers.R
#
# Discretize-then-Optimize (DtO) schemes for ODE-constrained optimal control.
#
# Each scheme implements a forward integrator and its consistent discrete
# adjoint derived from the discrete forward equations. The adjoint is NOT
# obtained by discretizing the continuous adjoint — it is the exact adjoint
# of the discrete forward map.
#
# Architecture:
#   DtOScheme (abstract)       — forward/adjoint interface
#   ├── EulerScheme            — explicit Euler, order 1
#   ├── CrankNicolsonScheme    — implicit trapezoidal, order 2 (self-adjoint)
#   ├── GL1Scheme              — implicit midpoint / cG(1), order 2 (self-adjoint)
#   └── GL2Scheme              — 2-stage Gauss-Legendre, order 4 (self-adjoint)
#
#   DtOSolver  — unified OCP solver wrapping any DtOScheme
#
# Scheme interface:
#   forward(rhs, y0, times, jac_fn, source_fn) -> list(y = ns*nv, aux)
#   adjoint(rhs, pT, jac_fn, source_fn)        -> list(p = ns*nv, grad_contrib = ns*nv)
#
# grad_contrib provides the ODE-constraint contribution to dJ/du, so the
# full control gradient is:  2*lambda*w*u + grad_contrib
#
# Dispatcher: make_dto_solver(method)  method in {"euler","cn","gl1","gl2"}
# Backward-compat: solve_ode(method, rhs, y0, times, jac_fn)
# =============================================================================

library(R6)

.num_jac <- function(rhs, y, t, eps = 1e-7) {
  n <- length(y); J <- matrix(0, n, n)
  for (j in seq_len(n)) {
    yp <- y; yp[j] <- y[j] + eps
    ym <- y; ym[j] <- y[j] - eps
    J[, j] <- (rhs(yp, t) - rhs(ym, t)) / (2 * eps)
  }
  J
}

# =============================================================================
# DtOScheme — abstract base class
#
# forward() must store times and y_fwd in private so that adjoint() can
# use the forward trajectory for Jacobian evaluation.
# =============================================================================
DtOScheme <- R6Class("DtOScheme",
  private = list(times = NULL, y_fwd = NULL, aux = NULL),
  public = list(
    forward = function(rhs, y0, times, jac_fn = NULL, source_fn = NULL)
      stop("Abstract: override in subclass"),
    adjoint = function(rhs, pT, jac_fn = NULL, source_fn = NULL)
      stop("Abstract: override in subclass")
  )
)

# =============================================================================
# Euler (explicit, order 1) — not self-adjoint
#
# Forward:  y[t+1] = y[t] + dt * rhs(y[t], t[t])
# Adjoint:  p[t]   = p[t+1] + dt * J(t)^T * p[t+1]
# Gradient: grad_contrib[t] = dt[t] * p[t+1]
# =============================================================================
EulerScheme <- R6Class("EulerScheme", inherit = DtOScheme,
  public = list(
    forward = function(rhs, y0, times, jac_fn = NULL, source_fn = NULL) {
      ns <- length(times); nv <- length(y0)
      dt <- diff(times)
      y  <- matrix(0, ns, nv); y[1L, ] <- y0
      if (!is.null(source_fn)) y[1L, ] <- y[1L, ] + source_fn(1L)
      for (t in seq_len(ns - 1L)) {
        y[t + 1L, ] <- y[t, ] + dt[t] * rhs(y[t, ], times[t])
        if (!is.null(source_fn)) y[t + 1L, ] <- y[t + 1L, ] + source_fn(t + 1L)
      }
      private$times <- times; private$y_fwd <- y
      list(y = y, aux = NULL)
    },

    adjoint = function(rhs, pT, jac_fn = NULL, source_fn = NULL) {
      times <- private$times; y_fwd <- private$y_fwd
      if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
      ns <- length(times); nv <- length(pT)
      dt <- diff(times)
      p  <- matrix(0, ns, nv)
      p[ns, ] <- pT
      if (!is.null(source_fn)) p[ns, ] <- p[ns, ] + source_fn(ns)
      for (t in seq.int(ns - 1L, 1L)) {
        J      <- jac_fn(y_fwd[t, ], times[t])
        p[t, ] <- p[t + 1L, ] + dt[t] * as.vector(t(J) %*% p[t + 1L, ])
        if (!is.null(source_fn)) p[t, ] <- p[t, ] + source_fn(t)
      }
      dt_r <- c(dt, 0)
      grad_contrib <- matrix(dt_r, ns, nv) *
                        rbind(p[-1L, , drop = FALSE], matrix(0, 1, nv))
      list(p = p, grad_contrib = grad_contrib)
    }
  )
)

# =============================================================================
# Crank-Nicolson (implicit trapezoidal, order 2) — self-adjoint
#
# Forward:  y[t+1] = y[t] + (h/2)*(rhs(y[t],t) + rhs(y[t+1],t+1))  [Newton]
# Adjoint:  (I - h_l/2*J^T)*p[t] = (I + h_r/2*J^T)*p[t+1] + source
# Gradient: grad_contrib[t] = (dt_r/2)*p[t+1] + (dt_l/2)*p[t]
# =============================================================================
CrankNicolsonScheme <- R6Class("CrankNicolsonScheme", inherit = DtOScheme,
  public = list(
    forward = function(rhs, y0, times, jac_fn = NULL, source_fn = NULL) {
      ns <- length(times); nv <- length(y0)
      if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
      dt <- diff(times)
      y  <- matrix(0, ns, nv); y[1L, ] <- y0
      for (t in seq_len(ns - 1L)) {
        h  <- dt[t]
        ft <- rhs(y[t, ], times[t])
        yn <- y[t, ] + h * ft
        if (!is.null(source_fn)) yn <- yn + source_fn(t + 1L)
        for (k in seq_len(10L)) {
          fn  <- rhs(yn, times[t + 1L])
          res <- yn - y[t, ] - (h / 2) * (ft + fn)
          if (!is.null(source_fn)) res <- res - source_fn(t + 1L)
          if (max(abs(res)) < 1e-12) break
          yn  <- yn - as.vector(solve(diag(nv) - (h / 2) * jac_fn(yn, times[t + 1L]), res))
        }
        y[t + 1L, ] <- yn
      }
      private$times <- times; private$y_fwd <- y
      list(y = y, aux = NULL)
    },

    adjoint = function(rhs, pT, jac_fn = NULL, source_fn = NULL) {
      times  <- private$times; y_fwd <- private$y_fwd
      if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
      ns     <- length(times); nv <- length(pT)
      dt_vec <- diff(times)
      p      <- matrix(0, ns, nv)
      Jf_T   <- jac_fn(y_fwd[ns, ], times[ns])
      rhs_ns <- pT
      if (!is.null(source_fn)) rhs_ns <- rhs_ns + source_fn(ns)
      p[ns, ] <- as.vector(solve(diag(nv) - (dt_vec[ns - 1L] / 2) * t(Jf_T), rhs_ns))
      for (t in seq.int(ns - 1L, 1L)) {
        Jf_t  <- jac_fn(y_fwd[t, ], times[t])
        dt_l  <- if (t > 1L) dt_vec[t - 1L] else 0
        rhs_t <- as.vector((diag(nv) + (dt_vec[t] / 2) * t(Jf_t)) %*% p[t + 1L, ])
        if (!is.null(source_fn)) rhs_t <- rhs_t + source_fn(t)
        p[t, ] <- as.vector(solve(diag(nv) - (dt_l / 2) * t(Jf_t), rhs_t))
      }
      dt_r         <- c(dt_vec, 0)
      dt_l_vec     <- c(0, dt_vec[seq_len(ns - 1L)])
      grad_contrib <- matrix(dt_r / 2, ns, nv) *
                        rbind(p[-1L, , drop = FALSE], matrix(0, 1, nv)) +
                      matrix(dt_l_vec / 2, ns, nv) * p
      list(p = p, grad_contrib = grad_contrib)
    }
  )
)

# =============================================================================
# GL1 — 1-stage Gauss-Legendre (implicit midpoint, order 2) — self-adjoint
#
# Forward:  ym = solve(ym - y[t] - h/2*rhs(ym, t_mid) = 0); y[t+1] = 2*ym - y[t]
# Adjoint:  uses midpoint Jacobians for consistent DtO adjoint
# Gradient: grad_contrib[t] = dt[t] * p[t+1]
# =============================================================================
GL1Scheme <- R6Class("GL1Scheme", inherit = DtOScheme,
  public = list(
    forward = function(rhs, y0, times, jac_fn = NULL, source_fn = NULL) {
      ns <- length(times); nv <- length(y0)
      if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
      dt <- diff(times)
      y  <- matrix(0, ns, nv); y[1L, ] <- y0
      if (!is.null(source_fn)) y[1L, ] <- y[1L, ] + source_fn(1L)
      for (t in seq_len(ns - 1L)) {
        h     <- dt[t]
        t_mid <- (times[t] + times[t + 1L]) / 2
        ym    <- y[t, ] + (h / 2) * rhs(y[t, ], t_mid)
        for (k in seq_len(10L)) {
          res <- ym - y[t, ] - (h / 2) * rhs(ym, t_mid)
          if (max(abs(res)) < 1e-12) break
          ym  <- ym - as.vector(solve(diag(nv) - (h / 2) * jac_fn(ym, t_mid), res))
        }
        y[t + 1L, ] <- 2 * ym - y[t, ]
        if (!is.null(source_fn)) y[t + 1L, ] <- y[t + 1L, ] + source_fn(t + 1L)
      }
      private$times <- times; private$y_fwd <- y
      list(y = y, aux = NULL)
    },

    adjoint = function(rhs, pT, jac_fn = NULL, source_fn = NULL) {
      times  <- private$times; y_fwd <- private$y_fwd
      if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
      ns     <- length(times); nv <- length(pT)
      dt_vec <- diff(times)
      jac_mid <- function(tl, tr)
        jac_fn((y_fwd[tl, ] + y_fwd[tr, ]) / 2, (times[tl] + times[tr]) / 2)
      p <- matrix(0, ns, nv)
      rhs_ns <- pT
      if (!is.null(source_fn)) rhs_ns <- rhs_ns + source_fn(ns)
      p[ns, ] <- as.vector(
        solve(diag(nv) - (dt_vec[ns - 1L] / 2) * t(jac_mid(ns - 1L, ns)), rhs_ns)
      )
      for (t in seq.int(ns - 1L, 1L)) {
        Jf_r  <- jac_mid(t, t + 1L)
        Jf_l  <- if (t > 1L) jac_mid(t - 1L, t) else matrix(0, nv, nv)
        dt_l  <- if (t > 1L) dt_vec[t - 1L] else 0
        rhs_t <- as.vector((diag(nv) + (dt_vec[t] / 2) * t(Jf_r)) %*% p[t + 1L, ])
        if (!is.null(source_fn)) rhs_t <- rhs_t + source_fn(t)
        p[t, ] <- as.vector(solve(diag(nv) - (dt_l / 2) * t(Jf_l), rhs_t))
      }
      dt_r         <- c(dt_vec, 0)
      grad_contrib <- matrix(dt_r, ns, nv) *
                        rbind(p[-1L, , drop = FALSE], matrix(0, 1, nv))
      list(p = p, grad_contrib = grad_contrib)
    }
  )
)

# =============================================================================
# GL2 — 2-stage Gauss-Legendre (order 4) — self-adjoint
#
# Butcher: c=(1/2 +/- sqrt(3)/6), b=(1/2,1/2),
#          A11=A22=1/4, A12=1/4-sqrt(3)/6, A21=1/4+sqrt(3)/6
# Stores stage values in aux (required for exact discrete adjoint).
# =============================================================================
.gl2_s3  <- sqrt(3) / 6
.gl2_A11 <- 0.25;             .gl2_A12 <- 0.25 - .gl2_s3
.gl2_A21 <- 0.25 + .gl2_s3;  .gl2_A22 <- 0.25
.gl2_c1  <- 0.5 - .gl2_s3;   .gl2_c2  <- 0.5 + .gl2_s3

GL2Scheme <- R6Class("GL2Scheme", inherit = DtOScheme,
  public = list(
    forward = function(rhs, y0, times, jac_fn = NULL, source_fn = NULL) {
      ns <- length(times); nv <- length(y0)
      if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
      dt     <- diff(times)
      y      <- matrix(0, ns, nv); y[1L, ] <- y0
      stages <- vector("list", ns - 1L)
      if (!is.null(source_fn)) y[1L, ] <- y[1L, ] + source_fn(1L)
      for (t in seq_len(ns - 1L)) {
        h   <- dt[t]; yt <- y[t, ]
        t1  <- times[t] + .gl2_c1 * h; t2 <- times[t] + .gl2_c2 * h
        K1  <- rhs(yt, times[t]); K2 <- K1
        for (k in seq_len(15L)) {
          Y1  <- yt + h * (.gl2_A11 * K1 + .gl2_A12 * K2)
          Y2  <- yt + h * (.gl2_A21 * K1 + .gl2_A22 * K2)
          G1  <- K1 - rhs(Y1, t1)
          G2  <- K2 - rhs(Y2, t2)
          res <- c(G1, G2)
          if (max(abs(res)) < 1e-12) break
          Jf1    <- jac_fn(Y1, t1); Jf2 <- jac_fn(Y2, t2)
          Jblock <- rbind(
            cbind(diag(nv) - h * .gl2_A11 * Jf1,  -h * .gl2_A12 * Jf1),
            cbind(-h * .gl2_A21 * Jf2,  diag(nv) - h * .gl2_A22 * Jf2)
          )
          dK <- as.vector(solve(Jblock, -res))
          K1 <- K1 + dK[seq_len(nv)]; K2 <- K2 + dK[nv + seq_len(nv)]
        }
        Y1 <- yt + h * (.gl2_A11 * K1 + .gl2_A12 * K2)
        Y2 <- yt + h * (.gl2_A21 * K1 + .gl2_A22 * K2)
        stages[[t]]  <- rbind(Y1, Y2)
        y[t + 1L, ] <- yt + h * 0.5 * (K1 + K2)
        if (!is.null(source_fn)) y[t + 1L, ] <- y[t + 1L, ] + source_fn(t + 1L)
      }
      private$times <- times; private$y_fwd <- y; private$aux <- stages
      list(y = y, aux = stages)
    },

    adjoint = function(rhs, pT, jac_fn = NULL, source_fn = NULL) {
      times  <- private$times; aux <- private$aux
      if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
      ns      <- length(times); nv <- length(pT)
      dt_vec  <- diff(times)
      lam_sum <- matrix(0, ns, nv)
      p       <- matrix(0, ns, nv)
      p[ns, ] <- pT
      if (!is.null(source_fn)) p[ns, ] <- p[ns, ] + source_fn(ns)
      for (t in seq.int(ns - 1L, 1L)) {
        h   <- dt_vec[t]
        Y1  <- aux[[t]][1L, ]; Y2 <- aux[[t]][2L, ]
        t1  <- times[t] + .gl2_c1 * h; t2 <- times[t] + .gl2_c2 * h
        Jf1 <- jac_fn(Y1, t1); Jf2 <- jac_fn(Y2, t2)
        pt1 <- p[t + 1L, ]
        Mlam <- rbind(
          cbind(diag(nv) - h * .gl2_A11 * t(Jf1),  -h * .gl2_A21 * t(Jf2)),
          cbind(-h * .gl2_A12 * t(Jf1),  diag(nv) - h * .gl2_A22 * t(Jf2))
        )
        lam  <- as.vector(solve(Mlam, c(h * 0.5 * pt1, h * 0.5 * pt1)))
        lam1 <- lam[seq_len(nv)]; lam2 <- lam[nv + seq_len(nv)]
        p[t, ] <- pt1 + as.vector(t(Jf1) %*% lam1) + as.vector(t(Jf2) %*% lam2)
        if (!is.null(source_fn)) p[t, ] <- p[t, ] + source_fn(t)
        lam_sum[t, ] <- lam1 + lam2
      }
      list(p = p, grad_contrib = lam_sum)
    }
  )
)

# =============================================================================
# DtOSolver — unified OCP solver wrapping any DtOScheme
#
# Provides a uniform interface for forward and adjoint solves regardless
# of the underlying discretization method.
# =============================================================================
DtOSolver <- R6Class("DtOSolver",
  private = list(scheme = NULL),
  public = list(
    initialize = function(scheme) {
      private$scheme <- scheme
    },
    solve_state = function(rhs, y0, times, jac_fn = NULL, source_fn = NULL) {
      private$scheme$forward(rhs, y0, times, jac_fn, source_fn)
    },
    solve_adjoint = function(rhs, pT, jac_fn = NULL, source_fn = NULL) {
      private$scheme$adjoint(rhs, pT, jac_fn, source_fn)
    }
  )
)

# =============================================================================
# Dispatcher and backward-compat wrappers
# =============================================================================
make_dto_scheme <- function(method) {
  switch(method,
    euler = EulerScheme$new(),
    cn    = CrankNicolsonScheme$new(),
    gl1   = GL1Scheme$new(),
    gl2   = GL2Scheme$new(),
    stop("Unknown method: ", method)
  )
}

make_dto_solver <- function(method) {
  DtOSolver$new(make_dto_scheme(method))
}

make_ode_solver <- make_dto_solver

solve_ode <- function(method, rhs, y0, times, jac_fn = NULL)
  make_dto_solver(method)$solve_state(rhs, y0, times, jac_fn)
