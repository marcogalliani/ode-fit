source("examples/param-cascading/ode_params_examples.R")

# Reproduce the original synthetic SB test using the canonical run_example.
# True parameters match ODE_CONFIGS$sb; we estimate E, n, m with A fixed.
run_example(
  cfg         = ODE_CONFIGS$sb,
  init_params = c(E = 50000, n = 3.0, m = 0.5),
  param_names = c("E", "n", "m"),
  lower_phys  = c(E = 10000, n = 0.1, m = 0.1),
  upper_phys  = c(E = 300000, n = 20,  m = 20),
  lambda      = 1e0
)
