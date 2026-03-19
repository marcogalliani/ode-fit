# =============================================================================
# src/solvers/bvp_solver.R
#
# Discrete two-point BVP solver (Riccati-based).
# Provides:
#   solve_linear_bvp_riccati  — exact one-pass solver for linear BVPs
# =============================================================================


# -----------------------------------------------------------------------------
# solve_linear_bvp_riccati
#
# Solves the coupled discrete linear two-point BVP:
#   Forward:   y[t+1] = A[t]*y[t] - C[t]*p[t+1] + b[t]    (y[1] = y0_mat)
#   Backward:  p[t]   = D[t]*p[t+1] + E[t]*y[t]  + f[t]   (p[T] = E_T*y[T] + f_T)
#
# The system is solved by the Riccati forward sweep (ansatz y[t] = R[t]*p[t] + s[t])
# followed by one backward substitution. Unconditionally stable.
#
# Arguments:
#   ns        integer        number of time steps
#   ny        integer        state / adjoint dimension
#   nrhs      integer        simultaneous RHS (np for sensitivity, 1 for inner BVP)
#   A_list    list[ns-1]     A[t]: ny x ny propagation (forward)
#   C_list    list[ns-1]     C[t]: ny x ny backward coupling (forward eq)
#   b_list    list[ns-1]     b[t]: ny x nrhs forcing (forward)
#   D_list    list[ns-1]     D[t]: ny x ny propagation (backward)
#   E_list    list[ns-1]     E[t]: ny x ny forward coupling (backward eq)
#   f_list    list[ns-1]     f[t]: ny x nrhs forcing (backward)
#   y0_mat    ny x nrhs      initial condition for y (zeros for sensitivity BVP)
#   E_T       ny x ny        terminal coupling:  p[T] = E_T*y[T] + f_T
#   f_T       ny x nrhs      terminal forcing    (default zeros)
#
# Returns: list(Y, P)  where Y and P are ns x ny x nrhs arrays.
# -----------------------------------------------------------------------------
solve_linear_bvp_riccati <- function(ns, ny, nrhs,
                                      A_list, C_list, b_list,
                                      D_list, E_list, f_list,
                                      y0_mat, E_T, f_T = NULL) {
  if (is.null(f_T)) f_T <- matrix(0, ny, nrhs)

  R_list <- vector("list", ns)        # R[t]: ny x ny  (same for all RHS)
  s_list <- vector("list", ns)        # s[t]: ny x nrhs
  M_list <- vector("list", ns - 1L)   # M[t] = (I - E[t]*R[t])^{-1}

  R_list[[1L]] <- matrix(0, ny, ny)
  s_list[[1L]] <- y0_mat

  # --- Riccati forward sweep ---
  for (t in seq_len(ns - 1L)) {
    Rt <- R_list[[t]];  st <- s_list[[t]]
    At <- A_list[[t]];  Ct <- C_list[[t]]
    Dt <- D_list[[t]];  Et <- E_list[[t]]
    bt <- b_list[[t]];  ft <- f_list[[t]]

    Mt <- solve(diag(ny) - Et %*% Rt)
    M_list[[t]] <- Mt

    qt <- Mt %*% (Et %*% st + ft)
    R_list[[t + 1L]] <- At %*% Rt %*% Mt %*% Dt - Ct
    s_list[[t + 1L]] <- At %*% (Rt %*% qt + st) + bt
  }

  # --- Terminal condition:  (I - R[T]*E_T)*y[T] = s[T] + R[T]*f_T ---
  R_T <- R_list[[ns]];  s_T <- s_list[[ns]]
  Y_T <- solve(diag(ny) - R_T %*% E_T, s_T + R_T %*% f_T)
  P_T <- E_T %*% Y_T + f_T

  Y_arr <- array(0, c(ns, ny, nrhs))
  P_arr <- array(0, c(ns, ny, nrhs))
  Y_arr[ns, , ] <- Y_T
  P_arr[ns, , ] <- P_T

  # --- Backward substitution ---
  for (t in seq.int(ns - 1L, 1L)) {
    Mt    <- M_list[[t]]
    Dt    <- D_list[[t]];  Et <- E_list[[t]]
    st    <- s_list[[t]];  ft <- f_list[[t]]
    Rt    <- R_list[[t]]
    Pnext <- matrix(P_arr[t + 1L, , ], ny, nrhs)

    Pt <- Mt %*% (Dt %*% Pnext + Et %*% st + ft)
    Yt <- Rt %*% Pt + st

    P_arr[t, , ] <- Pt
    Y_arr[t, , ] <- Yt
  }

  list(Y = Y_arr, P = P_arr)
}


# Internal helper: numerical Jacobian of func_rhs w.r.t. y.
.bvp_jac_y <- function(func_rhs, y, t, params, eps = 1e-7) {
  ny <- length(y)
  J  <- matrix(0, ny, ny)
  for (j in seq_len(ny)) {
    yp <- y; yp[j] <- y[j] + eps
    ym <- y; ym[j] <- y[j] - eps
    J[, j] <- (func_rhs(yp, t, params) - func_rhs(ym, t, params)) / (2 * eps)
  }
  J
}


# =============================================================================
# solve_bvp_colloc
#
# General nonlinear two-point BVP solver using trapezoidal collocation and
# Newton iteration — analogous to MATLAB bvp4c/bvp5c but without adaptive
# mesh refinement.  Does NOT require a pre-linearised system; it handles
# nonlinear F directly through Newton's method.
#
# Solves:
#   dz/dt = F(t, z),   bc(z(t[1]), z(t[ns])) = 0
#
# via the trapezoidal collocation conditions on each interval [t_k, t_{k+1}]:
#   z[t+1] - z[t] - (dt/2)*(F(t, z[t]) + F(t+1, z[t+1])) = 0    (nz eqs)
# together with the nz boundary conditions.
#
# The full system has N = ns*nz unknowns and N equations:
#   rows 1:nz           — bc(z[1], z[ns]) = 0
#   rows nz+1:N         — collocation equations, nz per interval
#
# The Jacobian of the collocation blocks is formed analytically from jac_F;
# the BC rows use numerical central differences.
#
# Arguments:
#   F_rhs         function(t, z) -> numeric[nz]  combined RHS
#   bc_residual   function(z_left, z_right) -> numeric[nz]
#                   first (nz/2) entries  = left  BCs
#                   last  (nz/2) entries  = right BCs
#   t_grid        numeric[ns]   ascending time grid
#   z_init        matrix[ns x nz]  initial guess
#   jac_F         optional function(t, z) -> matrix[nz x nz]
#                 Jacobian of F_rhs w.r.t. z; if NULL computed by central FD
#   max_iter      Newton iteration limit             (default 50)
#   tol           convergence: max|R|_inf < tol      (default 1e-8)
#   verbose       print Newton progress              (default FALSE)
#
# Returns: list(z, converged, iter, residual_norm)
#   z             matrix[ns x nz]  solution on t_grid
#   converged     logical
#   iter          integer — Newton iterations used
#   residual_norm numeric — final max|R|_inf
# =============================================================================
solve_bvp_colloc <- function(F_rhs, bc_residual, t_grid, z_init,
                              jac_F = NULL, max_iter = 50L, tol = 1e-8,
                              verbose = FALSE) {

  ns  <- length(t_grid)
  nz  <- ncol(z_init)
  N   <- ns * nz
  FD_eps <- 1e-7

  # Default Jacobian: central finite differences on F_rhs
  if (is.null(jac_F)) {
    jac_F <- function(t, z) {
      J <- matrix(0, nz, nz)
      for (j in seq_len(nz)) {
        zp <- z; zp[j] <- z[j] + FD_eps
        zm <- z; zm[j] <- z[j] - FD_eps
        J[, j] <- (F_rhs(t, zp) - F_rhs(t, zm)) / (2 * FD_eps)
      }
      J
    }
  }

  # Flat index range for time step t (row-major: z[t,v] at (t-1)*nz+v)
  idx <- function(t) (t - 1L) * nz + seq_len(nz)

  # Flatten initial guess (row-major)
  z_flat <- as.vector(t(z_init))

  # ---- Residual R(z_flat) ∈ R^N ----
  build_R <- function(zf) {
    R <- numeric(N)
    # BC block
    R[seq_len(nz)] <- bc_residual(zf[idx(1L)], zf[idx(ns)])
    # Trapezoidal collocation blocks
    for (t in seq_len(ns - 1L)) {
      dt  <- t_grid[t + 1L] - t_grid[t]
      zt  <- zf[idx(t)];       zt1 <- zf[idx(t + 1L)]
      Ft  <- F_rhs(t_grid[t],        zt)
      Ft1 <- F_rhs(t_grid[t + 1L],  zt1)
      R[nz + (t - 1L) * nz + seq_len(nz)] <- zt1 - zt - (dt / 2) * (Ft + Ft1)
    }
    R
  }

  # ---- Jacobian dR/dz (exploits block-bidiagonal sparsity) ----
  #
  # BC rows  (1:nz): depend on z[1] and z[ns]  → numerical FD over those cols
  # Colloc row-block t: depends on z[t] and z[t+1]
  #   d(colloc)/d(z_t)   = -I - (dt/2)*jac_F(t,   z_t)
  #   d(colloc)/d(z_t+1) =  I - (dt/2)*jac_F(t+1, z_t+1)
  build_J <- function(zf) {
    J <- matrix(0, N, N)
    h <- FD_eps

    # BC rows: perturb only z[1] and z[ns] columns
    for (j in c(idx(1L), idx(ns))) {
      zp <- zf; zp[j] <- zf[j] + h
      zm <- zf; zm[j] <- zf[j] - h
      J[seq_len(nz), j] <-
        (bc_residual(zp[idx(1L)], zp[idx(ns)]) -
         bc_residual(zm[idx(1L)], zm[idx(ns)])) / (2 * h)
    }

    # Collocation rows: analytical Jacobian from jac_F
    for (t in seq_len(ns - 1L)) {
      dt      <- t_grid[t + 1L] - t_grid[t]
      row_idx <- nz + (t - 1L) * nz + seq_len(nz)
      zt      <- zf[idx(t)];  zt1 <- zf[idx(t + 1L)]
      Jt      <- jac_F(t_grid[t],        zt)
      Jt1     <- jac_F(t_grid[t + 1L],  zt1)
      J[row_idx, idx(t)]      <- -(diag(nz) + (dt / 2) * Jt)
      J[row_idx, idx(t + 1L)] <-   diag(nz) - (dt / 2) * Jt1
    }
    J
  }

  converged  <- FALSE
  iter_final <- max_iter

  for (iter in seq_len(max_iter)) {
    R0       <- build_R(z_flat)
    res_norm <- max(abs(R0))
    if (verbose)
      cat(sprintf("  [BVP colloc] iter %2d  |R|_inf = %.3e\n", iter, res_norm))
    if (res_norm < tol) { converged <- TRUE; iter_final <- iter; break }

    J_mat <- build_J(z_flat)
    dz    <- tryCatch(
      solve(J_mat, -R0),
      error = function(e) solve(J_mat + 1e-10 * diag(N), -R0)
    )

    # Armijo backtracking line search
    step <- 1.0
    for (k in seq_len(12L)) {
      if (max(abs(build_R(z_flat + step * dz))) < res_norm) break
      step <- step * 0.5
    }
    z_flat <- z_flat + step * dz
  }

  list(
    z             = matrix(z_flat, ns, nz, byrow = TRUE),
    converged     = converged,
    iter          = iter_final,
    residual_norm = max(abs(build_R(z_flat)))
  )
}
