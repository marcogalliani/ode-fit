make.PinballLik <- function()
{
  # 1. THE VALUE FUNCTION ($fn)
  fn <- function(data, times, devals, pars, more) {
    # Extract Tau and Epsilon (smoothing factor) from 'more'
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    difs[is.na(difs)] <- 0
    weights <- checkweights(more$weights, more$whichobs, mat(difs))
    weights <- mat(weights)
    
    # Smoothed Pinball Loss Formula: 0.5 * sqrt(u^2 + eps) + (tau - 0.5) * u
    # Note: We sum this loss, rather than squaring differences
    smooth_abs <- sqrt(difs^2 + eps)
    loss_val <- 0.5 * smooth_abs + (tau - 0.5) * difs
    
    f <- apply(weights * loss_val, 1, sum)
    return(f)
  }
  
  # 2. GRADIENT W.R.T STATES ($dfdx)
  dfdx <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    difs[is.na(difs)] <- 0
    weights <- checkweights(more$weights, more$whichobs, mat(difs))
    weights <- mat(weights)
    
    # First Derivative of Loss w.r.t residual (u)
    # rho'(u) = 0.5 * u / sqrt(u^2 + eps) + (tau - 0.5)
    denom <- sqrt(difs^2 + eps)
    rho_prime <- 0.5 * (difs / denom) + (tau - 0.5)
    
    # Chain rule: dL/dx = dL/du * du/dx
    # du/dx = -dfdx (since u = y - f(x))
    # So we return: -1 * weights * rho_prime * dfdx
    
    dfdx_eval <- more$dfdx(times, devals, pars, more$more)
    
    g <- c()
    # Apply weights and rho_prime
    # Note: equivalent to SSE's 'difs' but using 'rho_prime'
    weighted_diffs <- weights * rho_prime
    
    for (i in 1:dim(dfdx_eval)[3]) {
      g <- cbind(g, apply(weighted_diffs * dfdx_eval[, , i], 1, sum))
    }
    return(-1 * g)
  }
  
  # 3. GRADIENT W.R.T Y ($dfdy)
  dfdy <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    difs[is.na(difs)] <- 0
    weights <- checkweights(more$weights, more$whichobs, mat(difs))
    weights <- mat(weights)
    
    denom <- sqrt(difs^2 + eps)
    rho_prime <- 0.5 * (difs / denom) + (tau - 0.5)
    
    return(weights * rho_prime)
  }
  
  # 4. GRADIENT W.R.T PARAMETERS ($dfdp)
  dfdp <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    difs[is.na(difs)] <- 0
    weights <- checkweights(more$weights, more$whichobs, mat(difs))
    weights <- mat(weights)
    
    denom <- sqrt(difs^2 + eps)
    rho_prime <- 0.5 * (difs / denom) + (tau - 0.5)
    
    dfdp_eval <- more$dfdp(times, devals, pars, more$more)
    
    g <- c()
    weighted_diffs <- weights * rho_prime
    
    for (i in 1:dim(dfdp_eval)[3]) {
      g <- cbind(g, apply(weighted_diffs * dfdp_eval[, , i], 1, sum))
    }
    return(-1 * g) 
  }
  
  # 5. HESSIAN W.R.T STATES ($d2fdx2)
  d2fdx2 <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    difs[is.na(difs)] <- 0
    
    dfdx_eval <- more$dfdx(times, devals, pars, more$more)
    d2fdx2_eval <- more$d2fdx2(times, devals, pars, more$more)
    
    weights <- checkweights(more$weights, more$whichobs, mat(difs))
    weights <- mat(weights)
    
    # Precompute Rho derivatives
    denom <- sqrt(difs^2 + eps)
    rho_prime <- 0.5 * (difs / denom) + (tau - 0.5)
    
    # Second derivative w.r.t residual
    # rho''(u) = 0.5 * eps / (u^2 + eps)^(3/2)
    rho_dbl_prime <- 0.5 * eps / (denom^3)
    
    H <- array(0, c(dim(devals), dim(devals)[2]))
    
    # Formula: d2L/dx2 = rho''(u)*(df/dx)^2 - rho'(u)*(d2f/dx2)
    # The minus sign on the second term is because u = y - f(x), so second deriv of u wrt x is negative
    
    for (i in 1:dim(d2fdx2_eval)[3]) {
      for (j in 1:dim(d2fdx2_eval)[4]) {
        # Term 1: Curvature from the loss function (rho_dbl_prime) * Jacobian squared
        term1 <- weights * rho_dbl_prime * dfdx_eval[, , j] * dfdx_eval[, , i]
        
        # Term 2: Gradient from the loss function (rho_prime) * Hessian of the model
        term2 <- -1 * weights * rho_prime * d2fdx2_eval[, , i, j]
        
        H[, i, j] <- apply(term1 + term2, 1, sum)
      }
    }
    return(H)
  }
  
  # 6. MIXED PARTIALS ($d2fdxdy)
  d2fdxdy <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    denom <- sqrt(difs^2 + eps)
    rho_dbl_prime <- 0.5 * eps / (denom^3)
    
    dfdx_eval <- more$dfdx(times, devals, pars, more$more)
    weights <- checkweights(more$weights, more$whichobs, mat(dfdx_eval[,,1]))
    weights <- mat(weights)
    weights[is.na(data)] <- 0
    
    # Mixed partial is: rho''(u) * (-dfdx)
    for (i in 1:dim(dfdx_eval)[3]) {
      dfdx_eval[, , i] <- -1 * weights * rho_dbl_prime * dfdx_eval[, , i]
    }
    return(aperm(dfdx_eval, c(1, 3, 2)))
  }
  
  # 7. HESSIAN W.R.T Y ($d2fdy2)
  d2fdy2 <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    denom <- sqrt(difs^2 + eps)
    rho_dbl_prime <- 0.5 * eps / (denom^3)
    
    r <- array(0, c(dim(data), dim(data)[2]))
    ind <- cbind(rep(1:dim(data)[1], dim(data)[2]), 
                 kronecker(cbind(1:dim(data)[2], 1:dim(data)[2]), rep(1, dim(data)[1])))
    
    weights <- checkweights(more$weights, more$whichobs, mat(data))
    weights <- mat(weights)
    
    # In SSE this is constant 2. Here it is rho''(u).
    r[ind] <- weights * rho_dbl_prime
    return(r)
  }
  
  # 8. MIXED X and P ($d2fdxdp)
  d2fdxdp <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    difs[is.na(difs)] <- 0
    weights <- checkweights(more$weights, more$whichobs, mat(difs))
    weights <- mat(weights)
    
    denom <- sqrt(difs^2 + eps)
    rho_prime <- 0.5 * (difs / denom) + (tau - 0.5)
    rho_dbl_prime <- 0.5 * eps / (denom^3)
    
    dfdx_eval <- more$dfdx(times, devals, pars, more$more)
    dfdp_eval <- more$dfdp(times, devals, pars, more$more)
    d2fdxdp_eval <- more$d2fdxdp(times, devals, pars, more$more)
    
    H <- array(0, c(dim(devals), length(pars)))
    
    for (i in 1:dim(d2fdxdp_eval)[3]) {
      for (j in 1:dim(d2fdxdp_eval)[4]) {
        # Similar logic to d2fdx2:
        # Term 1: rho''(u) * (df/dx) * (df/dp)
        term1 <- weights * rho_dbl_prime * dfdx_eval[, , i] * dfdp_eval[, , j]
        # Term 2: -rho'(u) * (d2f/dxdp)
        term2 <- -1 * weights * rho_prime * d2fdxdp_eval[, , i, j]
        
        H[, i, j] <- apply(term1 + term2, 1, sum)
      }
    }
    return(H)
  }
  
  # 9. MIXED Y and P ($d2fdydp)
  d2fdydp <- function(data, times, devals, pars, more) {
    tau <- if(!is.null(more$tau)) more$tau else 0.5
    eps <- if(!is.null(more$eps)) more$eps else 1e-6
    
    fdevals <- more$fn(times, devals, pars, more$more)
    difs <- data - fdevals
    denom <- sqrt(difs^2 + eps)
    rho_dbl_prime <- 0.5 * eps / (denom^3)
    
    dfdp_eval <- more$dfdp(times, devals, pars, more$more)
    weights <- checkweights(more$weights, more$whichobs, mat(dfdp_eval[,,1]))
    weights <- mat(weights)
    weights[is.na(data)] <- 0
    
    for (i in 1:dim(dfdp_eval)[3]) {
      dfdp_eval[, , i] <- -1 * weights * rho_dbl_prime * dfdp_eval[, , i]
    }
    return(dfdp_eval)
  }
  
  return(list(fn = fn, dfdx = dfdx, dfdy = dfdy, dfdp = dfdp, 
              d2fdx2 = d2fdx2, d2fdxdy = d2fdxdy, d2fdy2 = d2fdy2, 
              d2fdxdp = d2fdxdp, d2fdydp = d2fdydp))
}