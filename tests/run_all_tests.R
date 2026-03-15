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
  source("tests/test_general_ode_solver.R", local = TRUE)
})

# ---- Suite 2: CascadingOdeSolver ----
cat("\n### CascadingOdeSolver ###\n")
local({
  source("tests/test_cascading_ode_solver.R", local = TRUE)
})

# ---- Suite 3: BVP Solver ----
cat("\n### BVP Solver ###\n")
local({
  source("tests/test_bvp_solver.R", local = TRUE)
})

# ---- Combined summary ----
cat("\n================================================================\n")
test_summary()
cat("================================================================\n")
