# =============================================================================
# tests/run_all_tests.R
#
# Entry point — runs all test suites and prints a combined summary.
#
# HOW TO RUN (from repo root):
#   source("tests/run_all_tests.R")
# =============================================================================

cat("================================================================\n")
cat(" ode-fit test suite\n")
cat("================================================================\n")

# Ensure the working directory is the repo root regardless of how this
# script is invoked (Rscript, source() from RStudio, etc.).
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path) |> dirname())
} else {
  setwd(dirname(sys.frame(1)$ofile) |> dirname())
}

# Reset the shared result store before running all suites
source("tests/test_helpers.R")
.test_env$results <- list()

# ---- Suite 1: OdeSystemSolver ----
cat("\n### OdeSystemSolver ###\n")
local({
  source("tests/forward-solvers/test_general_ode_solver.R", local = TRUE)
})

# ---- Suite 2: CascadingOdeSolver ----
cat("\n### CascadingOdeSolver ###\n")
local({
  source("tests/inverse-solvers/test_cascading_ode_solver.R", local = TRUE)
})

# ---- Suite 3: TrackingOdeSolver ----
cat("\n### TrackingOdeSolver ###\n")
local({
  source("tests/inverse-solvers/test_tracking_ode_solver.R", local = TRUE)
})

# ---- Suite 4: ODE integrator accuracy ----
cat("\n### ODE integrator accuracy ###\n")
local({
  source("tests/numerical-solvers/test_ode_accuracy.R", local = TRUE)
})

# ---- Suite 5: BVP Solver (linear Riccati) ----
cat("\n### BVP Solver (Riccati) ###\n")
local({
  source("tests/numerical-solvers/test_bvp_solver.R", local = TRUE)
})

# ---- Suite 6: BVP Collocation solver ----
cat("\n### BVP Collocation + optimize_bvp ###\n")
local({
  source("tests/numerical-solvers/test_bvp_colloc.R", local = TRUE)
})

# ---- Suite 7: SQP OCP Solver ----
cat("\n### SQP OCP Solver ###\n")
local({
  source("tests/forward-solvers/test_sqp_ocp_solver.R", local = TRUE)
})

# ---- Combined summary ----
cat("\n================================================================\n")
test_summary()
cat("================================================================\n")
