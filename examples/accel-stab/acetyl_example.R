library(AccelStab)
library(ggplot2)

library(magrittr)
library(dplyr)
library(tidyr)

source("src/solvers/inverse-solvers/load_inverse_solvers.R")
source("examples/ode_models.R")
 
# data
dat = read.table("data/acetyl.txt", sep = "\t", header = T)

# fit
fit = step1_down(
    data = dat, y = "y", .time = "days", C = "Temp",
    temp_pred_C = 5, max_time_pred = 20) 
summary(fit$fit)

# plots
step1_plot_T(step1_down_object = fit,
  focus_T = 5,
  xname = "Days",
  yname = "Acetyl",
  ribbon = TRUE
) + theme_minimal()
 
step1_plot_pred(step1_down_object = fit,
  xname = "Days",
  yname = "Acetyl",
)

# fit with the cascading solver

## data manipulation
df_acetyl <- dat %>%
    group_by(days, Temp) %>%
    summarise(y = mean(y, na.rm=T), .groups = "drop") %>%
    arrange(days) %>%
    pivot_wider(names_from = Temp, values_from = y) %>%
    arrange(days)
df_acetyl[['5']] <- NA
df_acetyl <- df_acetyl[,c("days","5","15","25","37","50")]

obs_data <- df_acetyl %>% select(-days) %>% as.matrix()
temps <- as.numeric(colnames(obs_data))
temps_K <- 273.15 + temps
T_ref <- mean(temps_K)

times_grid <- seq(0,20,0.5)
y0 <- rep(mean(obs_data[1,],na.rm=T), length(temps))

## ode parameters
params <- list(
  k_ref = 1, E = 50000, m = 0, n = 5, C0 = 100,
  R = 8.314, T_ref = T_ref, T_vec = temps_K
)
param_scales <- list(E = 10000, n = 1, C0 = 1, k_ref = 1)
fixed_params <- params[!names(params) %in% names(param_scales)]

## solver
cascading_solver <- CascadingOdeSolver$new(
  func_rhs     = sb_kref_rhs,
  obs_times    = df_acetyl$days,
  times_sim    = times_grid,
  obs_values   = obs_data,
  init_state   = function(p) as.numeric(y0),
  fixed_params = fixed_params,
  lambda       = 1e2,
  param_scales = param_scales
)

init_params <- params[names(param_scales)]

result <- cascading_solver$optimize_parameters(
  init_theta_physical = unlist(init_params),
  param_names         = names(init_params)
)


df_fit_plot <- data.frame(cascading_solver$last_solver$y)
names(df_fit_plot) <- colnames(obs_data)
df_fit_plot$Days <- times_grid

df_fit_plot <- melt(df_fit_plot, id.vars = "Days",
                      variable.name = "Celsius", value.name = "Acetyl")


step1_plot_pred(step1_down_object = fit,
  xname = "Days",
  yname = "Acetyl",
) + 
    geom_line(data = df_fit_plot, aes(x = Days, y = Acetyl, color = Celsius)) +
    theme_minimal()

