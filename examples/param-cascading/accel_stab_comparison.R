source("src/solvers/general_ode_system_solver.R")
source("src/solvers/parameter_cascading.R")

library(AccelStab)
library(magrittr)
library(dplyr)
library(tidyr)
library(reshape2)

# Data----

data("antigenicity")

# Processing
df_clean <- antigenicity %>%
  subset(validA == 0) %>%
  group_by(time, Celsius) %>%
  summarise(conc = mean(conc, na.rm=TRUE), .groups = "drop") %>%
  arrange(time) %>%
  pivot_wider(names_from = Celsius, values_from = conc) %>%
  arrange(time)

df_val <- antigenicity %>%
  subset(validA == 1) %>%
  group_by(time, Celsius) %>%
  summarise(conc = mean(conc, na.rm=TRUE), .groups = "drop") %>%
  arrange(time) %>%
  pivot_wider(names_from = Celsius, values_from = conc) %>%
  arrange(time)

# Prepare Solver Inputs
obs_matrix <- df_clean %>% dplyr::select(-time) %>% as.matrix()

times_grid <- seq(0,3,0.01)
times_grid <- sort(unique((c(times_grid, df_clean$time))))

temp_cols <- as.numeric(colnames(obs_matrix))
T_vec_Kelvin <- temp_cols + 273.15


# reference temperetarue of the study
T_ref <- mean(T_vec_Kelvin)

# SB model----
sestak_berggren_rhs <- function(x_vec, t, p) {
  eps <- 1e-6
  x_safe <- pmax(pmin(x_vec, p$C0 - eps), eps)
  
  # k(T) = k_ref * exp( -E/R * (1/T - 1/T_ref) )
  k_T <- p$k_ref * exp(-(p$E / p$R) * (1/p$T_vec - 1/p$T_ref))
  
  scaling <- p$C0^(p$m + p$n - 1)
  term_prod <- (p$C0 - x_safe)^p$m
  term_react <- x_safe^p$n
  
  dx_dt <- - (k_T / scaling) * term_prod * term_react
  return(dx_dt)
}

## Parameters----
params <-list(
  k_ref = 1, 
  E = 50000,
  m = 0,
  n = 5,
  C0 = 100,
  R = 8.314,
  T_ref = T_ref,      
  T_vec = T_vec_Kelvin
)

# Initial Condition
y0 <-rep(params$C0,times=ncol(obs_matrix)) # Fill missing t=0 with C0

# Parameter cascading----
# Setup
param_scales <- list(E = 10000, n = 1, C0 = 1, k_ref = 1)
fixed_params <- params[!names(params) %in% names(param_scales)]


cascading_solver <- CascadingOdeSolver$new(
  func_rhs = sestak_berggren_rhs,
  obs_times = df_clean$time,
  times_sim = times_grid,
  obs_values = obs_matrix,
  y0 = y0,
  fixed_params = fixed_params,
  lambda = 1e-3,
  param_scales = param_scales
)

# Run
init_params <- params[names(param_scales)]

l_bounds <- c(E = 10000, n = 0.1, C0 = max(obs_matrix, na.rm=TRUE) + 0.1, k_ref = 1e-4)
u_bounds <- c(E = 300000, n = 20, C0 = 120, k_ref = 100)

# Run the solver
result <- cascading_solver$optimize_parameters(
  init_theta_physical = unlist(init_params),
  param_names = names(init_params),
  lower_phys = l_bounds,
  upper_phys = u_bounds
)

# --- Visualization ---
df_fit <- as.data.frame(cascading_solver$last_solver$y)
colnames(df_fit) <- colnames(obs_matrix)
df_fit$time <- times_grid
df_fit_long <- melt(df_fit, id.vars="time", variable.name="Celsius", value.name="Conc")

df_obs <- as.data.frame(obs_matrix)
colnames(df_obs) <- colnames(obs_matrix)
df_obs$time <- df_clean$time
df_obs_long <- melt(df_obs, id.vars="time", variable.name="Celsius", value.name="Conc")

df_val_long <- melt(df_val, id.vars="time", variable.name="Celsius", value.name="Conc")

p1 <- ggplot() +
  geom_point(data=df_obs_long, aes(x=time, y=Conc, color=Celsius), size=2, alpha=0.6) +
  geom_line(data=df_fit_long, aes(x=time, y=Conc, color=Celsius), size=1) +
  lims(y=c(0,100)) +
  labs(title = "Sestak-Berggren Model Fit", y = "Concentration") + 
  geom_point(data=df_val_long, aes(x=time, y=Conc, color=Celsius), size=2, pch=1) +
  theme_minimal()

print(p1)

# AccelStab comparison----
antigenicity$Validate = as.factor(ifelse(antigenicity$time <= 0.5, 0, 1))

res = step1_down(
  data = antigenicity, y = "conc", 
  .time = "time", C = "Celsius",
  validation = "Validate",
  parms = list(k1 = 50, k2 = 10000, k3 = 3, c0 = 100),
  max_time_pred = 3)

step1_plot_pred(res, yname = "Antigenicity") +
  geom_line(data=df_fit_long, aes(x=time, y=Conc, color=Celsius), linewidth=1,linetype = "dotted") +
  lims(y=c(0,100)) +
  labs(title = "Sestak-Berggren Model Fit", y = "Concentration")

