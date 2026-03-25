# =============================================================================
# examples/accel_stab_data.R
#
# Shared data loading and preprocessing for AccelStab antigenicity comparisons.
# Sourced by:
#   examples/non-linear-odes-solver/accel_stab_comparison.R
#   examples/param-cascading/accel_stab_comparison.R
#   examples/tracking-estimator/accel_stab_comparison.R
# =============================================================================

library(AccelStab)
library(magrittr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

source("examples/ode_models.R")
source("examples/plot_utils.R")

# --- Data loading and preprocessing -----------------------------------------
data("antigenicity")

df_clean <- antigenicity %>%
  subset(validA == 0) %>%
  group_by(time, Celsius) %>%
  summarise(conc = mean(conc, na.rm = TRUE), .groups = "drop") %>%
  arrange(time) %>%
  pivot_wider(names_from = Celsius, values_from = conc) %>%
  arrange(time)

df_val <- antigenicity %>%
  subset(validA == 1) %>%
  group_by(time, Celsius) %>%
  summarise(conc = mean(conc, na.rm = TRUE), .groups = "drop") %>%
  arrange(time) %>%
  pivot_wider(names_from = Celsius, values_from = conc) %>%
  arrange(time)

obs_matrix   <- df_clean %>% dplyr::select(-time) %>% as.matrix()
temp_cols     <- as.numeric(colnames(obs_matrix))
T_vec_Kelvin  <- temp_cols + 273.15
T_ref         <- mean(T_vec_Kelvin)

# --- Shared default parameters and grids -----------------------------------
as_default_params <- list(
  k_ref = 1, E = 50000, m = 0, n = 5, C0 = 100,
  R = 8.314, T_ref = T_ref, T_vec = T_vec_Kelvin
)

as_default_y0 <- rep(as_default_params$C0, times = ncol(obs_matrix))

as_times_grid <- function(dt = 0.001) {
  sort(unique(c(seq(0, 3, dt), df_clean$time)))
}

# --- AccelStab baseline comparison -----------------------------------------
run_accelstab_baseline <- function(df_fit_long) {
  antigenicity$Validate <- as.factor(
    ifelse(antigenicity$time <= 0.5, 0, 1)
  )
  res <- step1_down(
    data = antigenicity, y = "conc",
    .time = "time", C = "Celsius",
    validation = "Validate",
    parms = list(k1 = 50, k2 = 10000, k3 = 3, c0 = 100),
    max_time_pred = 3
  )
  p <- step1_plot_pred(res, yname = "Antigenicity") +
    geom_line(data = df_fit_long,
              aes(x = time, y = Conc, color = Celsius),
              linewidth = 1, linetype = "dotted") +
    lims(y = c(0, 100)) +
    labs(title = "AccelStab vs solver fit", y = "Concentration")
  print(p)
  invisible(res)
}
