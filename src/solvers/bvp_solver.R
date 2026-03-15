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
