build_inverse_solver_model <- function(cfg) {
  ODEModel$new(
    rhs = cfg$rhs,
    init_state = function(params) as.numeric(cfg$y0),
    fixed_params = cfg$fixed_params,
    param_scales = cfg$param_scales
  )
}
