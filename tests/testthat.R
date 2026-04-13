library(testthat)

# Numerical ODE solver
# test_dir("tests/testthat/numerical-solvers", reporter = "summary")

# Forward solvers
test_file("tests/testthat/forward-solvers/test-dto-fwd-solver.R")

# Inverse solvers
test_file("tests/testthat/inverse-solvers/test-cascading-ode-solver.R")
test_file("tests/testthat/inverse-solvers/test-tracking-ode-solver.R")
