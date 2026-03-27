# =============================================================================
# tests/test_ode_accuracy.R
#
# Convergence-rate tests for ODE integrators in src/solvers/ode_solvers.R.
#
# HOW TO RUN (from repo root):
#   source("tests/test_ode_accuracy.R")
#
# WHAT IS TESTED (scalar ODE  dy/dt = -y,  y(0) = 1,  y(t) = exp(-t)):
#   A1. Euler  — global error O(h^1)
#   A2. CN     — global error O(h^2)
#   A3. GL1    — global error O(h^2)
#   A4. GL2    — global error O(h^4)
#
# Test strategy:
#   Compute max |y_num - y_exact| on [0, 1] at four step sizes
#   h = 0.1, 0.05, 0.025, 0.0125.  Estimate the convergence rate from
#   log-log regression of error vs h.  Pass if estimated rate >= expected - 0.3.
# =============================================================================

source("tests/test_helpers.R")
source("src/solvers/numerical-solvers/load_numerical_solvers.R")

# ---------------------------------------------------------------------------
# Shared setup: scalar exponential decay, exact solution exp(-t)
# ---------------------------------------------------------------------------
.acc_rhs <- function(y, t) -y
.acc_jac <- function(y, t) matrix(-1, 1, 1)
.acc_y0  <- c(1.0)
.acc_T   <- 1.0
.acc_hs  <- c(0.1, 0.05, 0.025, 0.0125)

.convergence_rate <- function(method) {
  errs <- sapply(.acc_hs, function(h) {
    times <- seq(0, .acc_T, by = h)
    y_num <- solve_ode(method, .acc_rhs, .acc_y0, times, .acc_jac)$y[, 1]
    y_ex  <- exp(-times)
    max(abs(y_num - y_ex))
  })
  # Least-squares log-log slope
  log_h <- log(.acc_hs); log_e <- log(errs)
  cov(log_h, log_e) / var(log_h)
}

# ---------------------------------------------------------------------------
# A1. Euler — order 1
# ---------------------------------------------------------------------------
describe("A1: Euler integrator — convergence rate O(h^1)", {
  rate <- .convergence_rate("euler")
  test_that("Euler convergence rate >= 0.7 (expected ~1)", {
    expect_greater_than(rate, 0.7,
      sprintf("estimated rate = %.3f", rate))
  })
  test_that("Euler convergence rate <= 1.5 (not super-convergent)", {
    expect_less_than(rate, 1.5,
      sprintf("estimated rate = %.3f", rate))
  })
})

# ---------------------------------------------------------------------------
# A2. CN — order 2
# ---------------------------------------------------------------------------
describe("A2: CN integrator — convergence rate O(h^2)", {
  rate <- .convergence_rate("cn")
  test_that("CN convergence rate >= 1.7 (expected ~2)", {
    expect_greater_than(rate, 1.7,
      sprintf("estimated rate = %.3f", rate))
  })
  test_that("CN convergence rate <= 2.5", {
    expect_less_than(rate, 2.5,
      sprintf("estimated rate = %.3f", rate))
  })
})

# ---------------------------------------------------------------------------
# A3. GL1 — order 2
# ---------------------------------------------------------------------------
describe("A3: GL1 integrator — convergence rate O(h^2)", {
  rate <- .convergence_rate("gl1")
  test_that("GL1 convergence rate >= 1.7 (expected ~2)", {
    expect_greater_than(rate, 1.7,
      sprintf("estimated rate = %.3f", rate))
  })
  test_that("GL1 convergence rate <= 2.5", {
    expect_less_than(rate, 2.5,
      sprintf("estimated rate = %.3f", rate))
  })
})

# ---------------------------------------------------------------------------
# A4. GL2 — order 4
# ---------------------------------------------------------------------------
describe("A4: GL2 integrator — convergence rate O(h^4)", {
  rate <- .convergence_rate("gl2")
  test_that("GL2 convergence rate >= 3.7 (expected ~4)", {
    expect_greater_than(rate, 3.7,
      sprintf("estimated rate = %.3f", rate))
  })
  test_that("GL2 convergence rate <= 4.5", {
    expect_less_than(rate, 4.5,
      sprintf("estimated rate = %.3f", rate))
  })
})

# ---------------------------------------------------------------------------
# A5. Absolute accuracy sanity: at h=0.01 each method should be close
# ---------------------------------------------------------------------------
describe("A5: Absolute accuracy at h=0.01 (sanity)", {
  times_fine <- seq(0, 1, by = 0.01)
  y_ex       <- exp(-times_fine)

  thresholds <- list(euler = 5e-3, cn = 1e-4, gl1 = 1e-4, gl2 = 1e-8)
  for (method in c("euler", "cn", "gl1", "gl2")) {
    local({
      m   <- method
      thr <- thresholds[[m]]
      y_num <- solve_ode(m, .acc_rhs, .acc_y0, times_fine, .acc_jac)$y[, 1]
      err   <- max(abs(y_num - y_ex))
      test_that(sprintf("%s: max error at h=0.01 < %.0e", m, thr), {
        expect_less_than(err, thr,
          sprintf("max error = %.3e", err))
      })
    })
  }
})

test_summary()
