library(FME)

source("examples/ode_models.R")
source("examples/sensitivity-analysis/sensitivity_utils.R")

# Select configuration -------------------------------------------------------
# Change this line to switch models.  Available keys: decay, lv, sb,
# sb_aspack, sb_asymptote.  All configurations are defined in
# examples/ode_models.R (ODE_CONFIGS list).
active <- make_sens_config(ODE_CONFIGS$lv)

# ---- 1. Local sensitivity via central-difference Jacobian ------------------
cat("\n=== Local sensitivity analysis ===\n")
sens <- run_local_sensitivity(active, verbose = TRUE)

# ---- 2. Comparison: FME::sensFun (forward differences, normalised) ---------
# sensFun perturbs each active parameter by a small relative fraction and
# computes forward-difference sensitivities.  The resulting SENS data frame
# has columns: x (time), var (state name), <param_1>, ..., <param_k>.
# Correlation of columns reveals parameter co-variation along the trajectory.
cat("\n--- FME::sensFun parameter correlation ---\n")
sens_fme <- sensFun(
  f          = ode_solver,
  parms      = active$params[active$active_params],
  model      = active$model,
  init_state = active$state,
  t_grid     = active$t_grid,
  rhs_args   = active$rhs_args,
  base_pars  = active$params
)
print(round(cor(sens_fme[, -(1:2)]), 4))
plot(sens_fme,
     main = paste("sensFun —", paste(active$active_params, collapse = ", ")))

# ---- 3. Identifiability summary --------------------------------------------
cat("\n--- Identifiability summary ---\n")
sv   <- svd(sens$FIM)$d
rank <- qr(sens$FIM)$rank
cond <- max(sv) / max(min(sv), .Machine$double.eps)

if (rank < length(active$active_params)) {
  cat("  WARNING: FIM is rank-deficient — parameters are not jointly identifiable.\n")
} else if (cond > 1e10) {
  cat("  WARNING: FIM is ill-conditioned (cond =", format(cond, scientific = TRUE),
      ") — near-collinear parameters.\n")
} else {
  cat("  All active parameters appear identifiable at these parameter values.\n")
}

# ---- 4. Notes on the FME vs central-difference comparison ------------------
# FME::sensFun uses forward (one-sided) differences and normalises
# sensitivities by p / y, yielding dimensionless elasticities.
# compute_jacobian() uses central (two-sided) differences (O(h^2) accuracy)
# without normalisation, yielding absolute dy/dp values.
# Both should agree on the sign and rank ordering of parameter influence;
# differences in magnitude reflect the normalisation convention only.
