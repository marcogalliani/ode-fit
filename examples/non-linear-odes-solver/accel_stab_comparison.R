source("src/solvers/general_ode_system_solver.R")
source("src/utils/trace_optimisation.R")


library(AccelStab)
library(magrittr)
library(dplyr)
library(tidyr)
library(reshape2)

# Debug----
# OdeSystemSolver$debug("optimize")

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

times_grid <- seq(0,3,0.001)
times_grid <- sort(unique((c(times_grid, df_clean$time))))

temp_cols <- as.numeric(colnames(obs_matrix))
T_vec_Kelvin <- temp_cols + 273.15


# reference temperetarue of the study
T_ref <- mean(T_vec_Kelvin)
R <- 8.314

# SB model----
sestak_berggren_rhs <- function(x_vec, t, p) {
  eps <- 1e-6
  x_safe <- pmax(pmin(x_vec, p$C0 - eps), eps)
  
  # k(T) = k_ref * exp( -E/R * (1/T - 1/T_ref) )
  k_T <- p$k_ref * exp(-(p$E / R) * (1/T_vec_Kelvin - 1/T_ref))
  
  scaling <- p$C0^(p$m + p$n - 1)
  term_prod <- (p$C0 - x_safe)^p$m
  term_react <- x_safe^p$n
  
  dx_dt <- - (k_T / scaling) * term_prod * term_react
  return(dx_dt)
}

## Parameters----
params <-list(
  k_ref = 8, 
  E = 99489.9,
  m = 0,
  n = 8.472189,
  C0 = 98.05758
)

# Initial Condition
y0 <-rep(params$C0,times=ncol(obs_matrix)) # Fill missing t=0 with C0

# 3. Solver
solver <- OdeSystemSolver$new(
  func_rhs = sestak_berggren_rhs,
  times_sim = times_grid,
  obs_times = df_clean$time,
  obs_values = obs_matrix,
  params = params,
  lambda = 1e2
)


# 4. Optimization
solver$optimize(y0 = rep(params$C0,times=ncol(obs_matrix)), max_iter = 100)


# Visualization----
## fit----
df_fit <- as.data.frame(solver$y)
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

## diagnostics----
df_diagn <- as.data.frame(solver$u)
colnames(df_diagn) <- colnames(obs_matrix)
df_diagn$time <- times_grid
df_diagn_long <- melt(df_diagn, id.vars="time", variable.name="Celsius", value.name="Conc")

p_diagn <- ggplot() +
  geom_line(data=df_diagn_long, aes(x=time, y=Conc, color=Celsius), size=1) +
  theme_minimal()

print(p_diagn)

# AccelStab comparison----
antigenicity$Validate = as.factor(ifelse(antigenicity$time <= 0.5, 0, 1))

res = step1_down(
  data = antigenicity, y = "conc", 
  .time = "time", C = "Celsius",
  validation = "Validate",
  parms = list(k1 = 50, k2 = 10000, k3 = 3, c0 = 100),
  max_time_pred = 3)

pdf("nl_smoother.pdf", width = 8, height = 5)
step1_plot_pred(res, yname = "Antigenicity") +
  geom_line(data=df_fit_long, aes(x=time, y=Conc, color=Celsius), linewidth=1,linetype = "twodash") +
  lims(y=c(30,100)) +
  labs(title = "Sestak-Berggren Model Fit", y = "Concentration")
dev.off()


# Trace the optimisation----
trace_results <- trace_optimization(solver, y0 = y0, max_iter = 50)
# -> points  
plot(df_clean$time, obs_matrix[,1], col="black", pch=16, cex=0.5, 
     main="Optimization Progression", xlab="Time", ylab="Y",ylim=c(80,110))

# -> history 
history <- trace_results$history
n_snaps <- length(history)
for(i in seq_along(history)) {
  # Calculate a fade from light blue to dark blue
  alpha_val <- seq(0.1, 1, length.out = n_snaps)[i]
  lines(times_grid, history[[i]][,1], 
        col = rgb(0, 0, 1, alpha = alpha_val), lwd = 1)
}
# -> final fit
lines(times_grid, solver$solve_state(matrix(trace_results$res$par, solver$n_steps), y0)[,1], 
      col="red", lwd=2)

