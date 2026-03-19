# =============================================================================
# src/solvers/ode_solvers.R
#
# Generic ODE integrators.  Solves:
#   dy/dt = rhs(y, t),   y(times[1]) = y0
#
# Interface (uniform across all methods):
#   solve_X(rhs, y0, times, jac_fn = NULL, mode = "forward")
#   solve_ode(method, rhs, y0, times, jac_fn = NULL, mode = "forward")
#
# Arguments:
#   rhs      function(y, t) -> numeric[nv]
#   y0       initial condition at times[1] (forward) or times[ns] (backward)
#   times    numeric[ns] ascending time grid
#   jac_fn   function(y, t) -> matrix[nv x nv]  (default: central FD from rhs)
#   mode     "forward" : integrate times[1] -> times[ns]
#            "backward": integrate times[ns] -> times[1],  y0 at times[ns]
#              internally reverses times and output; same code runs both ways
#
# Returns: list(y = ns×nv trajectory, aux = method-specific or NULL)
#   gl2 forward: aux is a list[ns-1] of 2×nv stage matrices per interval.
#   backward / all other methods: aux = NULL.
#
# Supported methods: "euler", "cn", "gl1", "gl2"
# =============================================================================


.num_jac <- function(rhs, y, t, eps = 1e-7) {
  n <- length(y); J <- matrix(0, n, n)
  for (j in seq_len(n)) {
    yp <- y; yp[j] <- y[j] + eps
    ym <- y; ym[j] <- y[j] - eps
    J[, j] <- (rhs(yp, t) - rhs(ym, t)) / (2 * eps)
  }
  J
}

# Euler (explicit, order 1)
solve_euler <- function(rhs, y0, times, mode = "forward") {
  if (mode == "backward") {
    res <- solve_euler(rhs, y0, rev(times))
    n   <- nrow(res$y)
    return(list(y = res$y[n:1L, , drop = FALSE], aux = NULL))
  }
  ns <- length(times); nv <- length(y0)
  dt <- diff(times)
  y  <- matrix(0, ns, nv); y[1L, ] <- y0
  for (t in seq_len(ns - 1L))
    y[t + 1L, ] <- y[t, ] + dt[t] * rhs(y[t, ], times[t])
  list(y = y, aux = NULL)
}

# Crank-Nicolson (implicit trapezoidal, order 2)
solve_cn <- function(rhs, y0, times, jac_fn = NULL, mode = "forward") {
  if (mode == "backward") {
    res <- solve_cn(rhs, y0, rev(times), jac_fn)
    n   <- nrow(res$y)
    return(list(y = res$y[n:1L, , drop = FALSE], aux = NULL))
  }
  ns <- length(times); nv <- length(y0)
  if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
  dt <- diff(times)
  y  <- matrix(0, ns, nv); y[1L, ] <- y0

  for (t in seq_len(ns - 1L)) {
    h  <- dt[t]
    ft <- rhs(y[t, ], times[t])
    yn <- y[t, ] + h * ft
    for (k in seq_len(10L)) {
      fn  <- rhs(yn, times[t + 1L])
      res <- yn - y[t, ] - (h / 2) * (ft + fn)
      if (max(abs(res)) < 1e-12) break
      yn  <- yn - as.vector(solve(diag(nv) - (h / 2) * jac_fn(yn, times[t + 1L]), res))
    }
    y[t + 1L, ] <- yn
  }
  list(y = y, aux = NULL)
}

# GL1 — 1-stage Gauss-Legendre (implicit midpoint, order 2)
solve_gl1 <- function(rhs, y0, times, jac_fn = NULL, mode = "forward") {
  if (mode == "backward") {
    res <- solve_gl1(rhs, y0, rev(times), jac_fn)
    n   <- nrow(res$y)
    return(list(y = res$y[n:1L, , drop = FALSE], aux = NULL))
  }
  ns <- length(times); nv <- length(y0)
  if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
  dt <- diff(times)
  y  <- matrix(0, ns, nv); y[1L, ] <- y0

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
  }
  list(y = y, aux = NULL)
}

# GL2 — 2-stage Gauss-Legendre (order 4)
# Butcher tableau: c1=1/2-√3/6, c2=1/2+√3/6, b1=b2=1/2,
#   A11=A22=1/4, A12=1/4-√3/6, A21=1/4+√3/6.
# For forward mode, stores stage values in aux (needed for exact discrete adjoint).
.gl2_s3  <- sqrt(3) / 6
.gl2_A11 <- 0.25;             .gl2_A12 <- 0.25 - .gl2_s3
.gl2_A21 <- 0.25 + .gl2_s3;  .gl2_A22 <- 0.25
.gl2_c1  <- 0.5 - .gl2_s3;   .gl2_c2  <- 0.5 + .gl2_s3

solve_gl2 <- function(rhs, y0, times, jac_fn = NULL, mode = "forward") {
  if (mode == "backward") {
    res <- solve_gl2(rhs, y0, rev(times), jac_fn)
    n   <- nrow(res$y)
    return(list(y = res$y[n:1L, , drop = FALSE], aux = rev(res$aux)))
  }
  ns <- length(times); nv <- length(y0)
  if (is.null(jac_fn)) jac_fn <- function(y, t) .num_jac(rhs, y, t)
  dt     <- diff(times)
  y      <- matrix(0, ns, nv); y[1L, ] <- y0
  stages <- vector("list", ns - 1L)

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
  }
  list(y = y, aux = stages)
}

# Dispatcher
solve_ode <- function(method, rhs, y0, times, jac_fn = NULL, mode = "forward") {
  switch(method,
    euler = solve_euler(rhs, y0, times,         mode = mode),
    cn    = solve_cn(   rhs, y0, times, jac_fn, mode = mode),
    gl1   = solve_gl1(  rhs, y0, times, jac_fn, mode = mode),
    gl2   = solve_gl2(  rhs, y0, times, jac_fn, mode = mode),
    stop("Unknown method: ", method)
  )
}


# =============================================================================
# Adjoint integrators
#
# Exact discrete adjoints of the corresponding forward integrators.
# Generic in the sense that physics enters only through source_fn and jac_fn.
#
# Interface (uniform across all methods):
#   solve_adjoint_X(y_fwd, aux, times, dt_vec, jac_fn, source_fn)
#   solve_adjoint_ode(method, y_fwd, aux, times, dt_vec, jac_fn, source_fn)
#
# Arguments:
#   y_fwd      ns×nv forward trajectory (used for Jacobian evaluation points)
#   aux        method-specific auxiliary data: NULL for cn/gl1/euler, stage
#              list[ns-1] for gl2
#   times      ns ascending time grid
#   dt_vec     c(diff(times), 0)
#   jac_fn     function(y, t) -> nv×nv Jacobian of the ODE rhs w.r.t. y
#   source_fn  function(t_idx) -> numeric[nv]  additive source at grid index t
#              (e.g. observation residual scaled by 2/ns)
#
# Returns: list(p = ns×nv adjoint, grad_contrib = ns×nv)
#   grad_contrib encodes the physics term in the gradient of u:
#     grad[t] = reg_term[t] + grad_contrib[t]
#   euler: grad_contrib = NULL (euler is typically used for warm-starts only)
# =============================================================================

solve_adjoint_euler <- function(y_fwd, aux, times, dt_vec, jac_fn, source_fn) {
  ns <- nrow(y_fwd); nv <- ncol(y_fwd)
  p  <- matrix(0, ns, nv)
  p[ns, ] <- source_fn(ns)
  for (t in seq.int(ns - 1L, 1L)) {
    J      <- jac_fn(y_fwd[t, ], times[t])
    p[t, ] <- p[t + 1L, ] + dt_vec[t] * as.vector(t(J) %*% p[t + 1L, ]) + source_fn(t)
  }
  list(p = p, grad_contrib = NULL)
}

solve_adjoint_cn <- function(y_fwd, aux, times, dt_vec, jac_fn, source_fn) {
  ns <- nrow(y_fwd); nv <- ncol(y_fwd)
  p  <- matrix(0, ns, nv)

  Jf_T    <- jac_fn(y_fwd[ns, ], times[ns])
  p[ns, ] <- as.vector(solve(t(diag(nv) - (dt_vec[ns - 1L] / 2) * Jf_T), source_fn(ns)))

  for (t in seq.int(ns - 1L, 1L)) {
    Jf_t   <- jac_fn(y_fwd[t, ], times[t])
    dt_l   <- if (t > 1L) dt_vec[t - 1L] else 0
    rhs    <- t(diag(nv) + (dt_vec[t] / 2) * Jf_t) %*% p[t + 1L, ] + source_fn(t)
    p[t, ] <- as.vector(solve(t(diag(nv) - (dt_l / 2) * Jf_t), rhs))
  }

  dt_l <- c(0, dt_vec[seq_len(ns - 1L)])
  grad_contrib <- matrix(dt_vec / 2, ns, nv) * rbind(p[-1L, , drop = FALSE], matrix(0, 1, nv)) +
                  matrix(dt_l  / 2, ns, nv) * p
  list(p = p, grad_contrib = grad_contrib)
}

solve_adjoint_gl1 <- function(y_fwd, aux, times, dt_vec, jac_fn, source_fn) {
  ns <- nrow(y_fwd); nv <- ncol(y_fwd)
  p  <- matrix(0, ns, nv)

  jac_mid <- function(tl, tr)
    jac_fn((y_fwd[tl, ] + y_fwd[tr, ]) / 2, (times[tl] + times[tr]) / 2)

  p[ns, ] <- as.vector(
    solve(t(diag(nv) - (dt_vec[ns - 1L] / 2) * jac_mid(ns - 1L, ns)), source_fn(ns))
  )
  for (t in seq.int(ns - 1L, 1L)) {
    Jf_r   <- jac_mid(t, t + 1L)
    Jf_l   <- if (t > 1L) jac_mid(t - 1L, t) else matrix(0, nv, nv)
    dt_l   <- if (t > 1L) dt_vec[t - 1L] else 0
    rhs    <- t(diag(nv) + (dt_vec[t] / 2) * Jf_r) %*% p[t + 1L, ] + source_fn(t)
    p[t, ] <- as.vector(solve(t(diag(nv) - (dt_l / 2) * Jf_l), rhs))
  }

  grad_contrib <- matrix(dt_vec, ns, nv) * rbind(p[-1L, , drop = FALSE], matrix(0, 1, nv))
  list(p = p, grad_contrib = grad_contrib)
}

solve_adjoint_gl2 <- function(y_fwd, aux, times, dt_vec, jac_fn, source_fn) {
  ns  <- nrow(y_fwd); nv <- ncol(y_fwd)
  p   <- matrix(0, ns, nv)
  lam_sum <- matrix(0, ns, nv)

  p[ns, ] <- source_fn(ns)
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
    p[t, ]       <- pt1 + as.vector(t(Jf1) %*% lam1) + as.vector(t(Jf2) %*% lam2) +
                    source_fn(t)
    lam_sum[t, ] <- lam1 + lam2
  }
  list(p = p, grad_contrib = lam_sum)
}

# Adjoint dispatcher
solve_adjoint_ode <- function(method, y_fwd, aux, times, dt_vec, jac_fn, source_fn) {
  switch(method,
    euler = solve_adjoint_euler(y_fwd, aux, times, dt_vec, jac_fn, source_fn),
    cn    = solve_adjoint_cn(   y_fwd, aux, times, dt_vec, jac_fn, source_fn),
    gl1   = solve_adjoint_gl1(  y_fwd, aux, times, dt_vec, jac_fn, source_fn),
    gl2   = solve_adjoint_gl2(  y_fwd, aux, times, dt_vec, jac_fn, source_fn),
    stop("Unknown method: ", method)
  )
}
