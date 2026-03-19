library(R6)
library(ggplot2)
library(gridExtra)

# ==============================================================================
# R6 Class: ODE Sestak-Berggren Optimization
# Fits a physics-informed forcing term u(t) to the SB kinetics model.
# Equation: d(alpha)/dt = k(T) * alpha^m * (1-alpha)^n + u(t)
# ==============================================================================
OdeSestakBerggren <- R6Class("OdeSestakBerggren",
                             public = list(
                               # Physics Parameters
                               A = NULL,       # Pre-exponential factor
                               E = NULL,       # Activation Energy
                               m = NULL,       # Reaction order 1
                               n = NULL,       # Reaction order 2
                               R_gas = 8.314,  # Gas constant
                               
                               # Regularization
                               lambda = NULL,
                               
                               # System State
                               n_steps = NULL,
                               dt = NULL,
                               times = NULL,
                               Temp = NULL,    # Temperature profile T(t) vector
                               
                               # State Vectors
                               alpha = NULL,   # Extent of reaction (0 to 1)
                               u = NULL,       # Unknown control/forcing term
                               p = NULL,       # Adjoint state
                               
                               # Observations
                               observations = NULL,
                               
                               # --------------------------------------------------------------------------
                               # Initialization
                               # --------------------------------------------------------------------------
                               initialize = function(params, lambda, times, Temp, observations) {
                                 self$A <- params$A
                                 self$E <- params$E
                                 self$m <- params$m
                                 self$n <- params$n
                                 
                                 self$lambda <- lambda
                                 self$times <- times
                                 self$n_steps <- length(times)
                                 self$dt <- times[2] - times[1]
                                 self$Temp <- Temp
                                 self$observations <- observations
                                 
                                 # Init vectors
                                 self$alpha <- numeric(self$n_steps)
                                 self$u <- numeric(self$n_steps)
                                 self$p <- numeric(self$n_steps)
                               },
                               
                               # --------------------------------------------------------------------------
                               # Helper: Sestak-Berggren Rate Calculation
                               # f(alpha, t) = A * exp(-E/RT) * alpha^m * (1-alpha)^n
                               # --------------------------------------------------------------------------
                               sb_rate = function(a_val, t_idx) {
                                 # Clamp alpha to (0, 1) to avoid NaN in powers
                                 eps <- 1e-6
                                 a_safe <- max(min(a_val, 1 - eps), eps)
                                 
                                 k <- self$A * exp(-self$E / (self$R_gas * self$Temp[t_idx]))
                                 rate <- k * (a_safe ^ self$m) * ((1 - a_safe) ^ self$n)
                                 return(rate)
                               },
                               
                               # --------------------------------------------------------------------------
                               # Helper: Derivative of SB Rate w.r.t alpha
                               # df/dalpha = Rate * [ m/alpha - n/(1-alpha) ]
                               # --------------------------------------------------------------------------
                               sb_deriv = function(a_val, t_idx) {
                                 eps <- 1e-6
                                 a_safe <- max(min(a_val, 1 - eps), eps)
                                 
                                 rate <- self$sb_rate(a_safe, t_idx)
                                 
                                 # Analytical derivative chain rule
                                 term1 <- self$m / a_safe
                                 term2 <- self$n / (1 - a_safe)
                                 
                                 df_da <- rate * (term1 - term2)
                                 return(df_da)
                               },
                               
                               # --------------------------------------------------------------------------
                               # Forward Solver: Explicit Euler
                               # alpha[t] = alpha[t-1] + dt * ( SB_Rate + u[t] )
                               # --------------------------------------------------------------------------
                               solve_state = function(u_curr, alpha0) {
                                 alpha_new <- numeric(self$n_steps)
                                 alpha_new[1] <- alpha0
                                 
                                 for (t in 2:self$n_steps) {
                                   # Calculate kinetics based on previous state
                                   kinetics <- self$sb_rate(alpha_new[t-1], t-1)
                                   
                                   # Update
                                   alpha_new[t] <- alpha_new[t-1] + self$dt * (kinetics + u_curr[t])
                                   
                                   # Hard clamp to physics limits [0,1] isn't applied strictly inside the 
                                   # loop to allow gradient flow, but soft clamping happens in rate calc.
                                 }
                                 return(alpha_new)
                               },
                               
                               # --------------------------------------------------------------------------
                               # Backward Solver: Adjoint Equation
                               # --------------------------------------------------------------------------
                               solve_adjoint = function(alpha_curr) {
                                 p_new <- numeric(self$n_steps)
                                 p_new[self$n_steps] <- 0 # Terminal condition
                                 
                                 # Backward loop
                                 for (t in (self$n_steps):2) {
                                   
                                   # 1. Forcing term from data mismatch (derivative of SSE)
                                   # Scale by 1/N to match typical mean loss scaling
                                   forcing <- (2 / self$n_steps) * (alpha_curr[t-1] - self$observations[t-1])
                                   
                                   # 2. Jacobian of the dynamics (df/dalpha)
                                   df_da <- self$sb_deriv(alpha_curr[t-1], t-1)
                                   
                                   # 3. Discrete Adjoint Update (Explicit Euler Adjoint)
                                   # p[t-1] = p[t] * (1 + dt * df_da) + forcing_term (simplified)
                                   # Note: The exact form depends on the Lagrangian definition.
                                   # For a_t = a_{t-1} + dt*f(...), sens is (1 + dt*df/da).
                                   
                                   p_new[t-1] <- p_new[t] * (1 + self$dt * df_da) - self$observations[t-1] * 0 # (Correction below)
                                   
                                   # Correct Adjoint accumulation for Least Squares:
                                   # p[t-1] = p[t] + dt * (df_da * p[t]) - forcing
                                   # Note: Depending on sign convention of Lagrangian, 'forcing' might be + or -
                                   # Here we use: p_prev = p_next + dt * p_next * J_f + gradient_of_loss_wrt_state
                                   
                                   # Loss derivative w.r.t state (t-1)
                                   dL_da <- (2 / self$n_steps) * (alpha_curr[t-1] - self$observations[t-1])
                                   
                                   p_new[t-1] <- p_new[t] + self$dt * (df_da * p_new[t]) + dL_da
                                 }
                                 return(p_new)
                               },
                               
                               # --------------------------------------------------------------------------
                               # Cost Function
                               # --------------------------------------------------------------------------
                               cost_function = function(u_flat, alpha0) {
                                 alpha_temp <- self$solve_state(u_flat, alpha0)
                                 
                                 sse <- sum((self$observations - alpha_temp)^2)
                                 reg <- self$lambda * sum(u_flat^2)
                                 
                                 return((1/self$n_steps) * sse + reg)
                               },
                               
                               # --------------------------------------------------------------------------
                               # Gradient Function
                               # Returns dJ/du
                               # --------------------------------------------------------------------------
                               gradient_function = function(u_flat, alpha0) {
                                 # 1. Forward
                                 alpha_temp <- self$solve_state(u_flat, alpha0)
                                 
                                 # 2. Backward
                                 p_temp <- self$solve_adjoint(alpha_temp)
                                 
                                 # 3. Compute Gradient w.r.t u
                                 # For Explicit Euler: a_t = ... + dt * u_t
                                 # The sensitivity da_t/du_t = dt
                                 # Gradient component from adjoint = p[t] * dt (approx, aligned indices)
                                 
                                 # Shifting: p[t] is the Lagrange multiplier for the step t-1 -> t.
                                 # u[t] drives that step.
                                 
                                 grad_u <- 2 * self$lambda * u_flat
                                 
                                 # Add adjoint contribution
                                 # We align indices carefully. u[t] affects alpha[t].
                                 # In the discrete Lagrangian: sum p[t] * (alpha[t] - alpha[t-1] - dt*u[t] ...)
                                 # dL/du[t] = 2*lam*u[t] - p[t] * dt
                                 
                                 grad_u <- grad_u - p_temp * self$dt
                                 
                                 return(grad_u)
                               },
                               
                               # --------------------------------------------------------------------------
                               # Optimization Loop
                               # --------------------------------------------------------------------------
                               optimize = function(alpha0, u_init_guess, max_iter = 200) {
                                 res <- optim(
                                   par = u_init_guess,
                                   fn = self$cost_function,
                                   gr = self$gradient_function,
                                   alpha0 = alpha0,
                                   method = "BFGS",
                                   control = list(maxit = max_iter, trace = 1, fnscale = 1) 
                                 )
                                 
                                 self$u <- res$par
                                 self$alpha <- self$solve_state(self$u, alpha0)
                                 self$p <- self$solve_adjoint(self$alpha)
                                 
                                 return(res)
                               }
                             )
)

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

main <- function() {
  cat("Running Sestak-Berggren Kinetics Optimization...\n")
  
  # 1. Setup Time and Temperature
  T_start <- 273  # Kelvin
  beta <- 0   # Heating rate (10 K/min converted to K/sec)
  T_end <- 800
  
  # Time grid based on heating rate
  duration <- 300000#(T_end - T_start) / beta
  n_steps <- 10
  times <- seq(0, duration, length.out = n_steps)
  dt <- times[2] - times[1]
  
  # Linear Temperature Ramp
  Temp_profile <- T_start + beta * times
  
  # 2. Physics Parameters (Synthetic Truth)
  # Typical values for a polymer degradation or crystallization
  params_true <- list(A = 1e4, E = 50000, m = 0, n = 10) 
  
  # 3. Generate Synthetic "Truth" Data
  # We assume the true process has a slight deviation from the SB model
  # Deviation: A burst of acceleration in reaction mid-way
  u_true <- rep(0, n_steps)
  
  # Run simulation to get ground truth alpha
  model_gen <- OdeSestakBerggren$new(params_true, 0, times, Temp_profile, rep(0, n_steps))
  alpha_true <- model_gen$solve_state(u_true, alpha0 = 0.01)
  
  # Add noise to observation
  set.seed(123)
  obs_noise <- alpha_true + rnorm(n_steps, 0, 0.1)
  # Clamp noisy observations to realistic range [0,1] for plotting
  obs_noise <- pmax(pmin(obs_noise, 1), 0)
  
  # 4. Setup Estimation Model
  # We give the model the correct physics parameters, but it doesn't know about u_true
  # It must discover u(t) to fit the data.
  model_fit <- OdeSestakBerggren$new(params_true, lambda = 1e3, times, Temp_profile, obs_noise)

  # 5. Optimize
  u_guess <- rep(0, n_steps)
  cat("Starting Optimization...\n")
  result <- model_fit$optimize(alpha0 = 0, u_init_guess = u_guess)
  cat("Convergence:", result$convergence, "\n")
  
  # 6. Visualization
  df_res <- data.frame(
    Time = times,
    Temperature = Temp_profile,
    Obs = obs_noise,
    Alpha_Fitted = model_fit$alpha,
    Alpha_True = alpha_true,
    Control_Est = model_fit$u,
    Control_True = u_true
  )
  
  # Plot 1: Reaction Extent (Alpha)
  p1 <- ggplot(df_res, aes(x=Time)) +
    geom_point(aes(y=Obs), color="grey", alpha=0.5) +
    geom_line(aes(y=Alpha_True), color="black", linetype="dashed") +
    geom_line(aes(y=Alpha_Fitted), color="blue", size=1) +
    labs(title = "Sestak-Berggren Fit", y = "Alpha") + theme_minimal()
  
  # Plot 2: Recovered Physics Deviation (Control u)
  p2 <- ggplot(df_res, aes(x=Time)) +
    geom_line(aes(y=Control_True), color="black", linetype="dashed", size=0.8) +
    geom_line(aes(y=Control_Est), color="red", size=1) +
    labs(title = "Recovered Dynamics Deviation (u)", y = "u(t)") + theme_minimal()
  
  # Plot 3: Temperature Profile (Context)
  p3 <- ggplot(df_res, aes(x=Time, y=Temperature)) +
    geom_line(color="orange") +
    labs(title = "Temperature Profile", y = "T (K)") + theme_minimal()
  
  grid.arrange(p1, p2, p3, ncol=1)
}

main()