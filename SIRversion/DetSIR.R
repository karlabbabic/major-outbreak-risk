# Install and load the deSolve package
# install.packages("deSolve")
library(deSolve)
library(ggplot2)

# Define the SIR model equations
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -(beta * S * I)
    dI <- (beta * S * I ) - mu * I
    dR <- mu * I
    return(list(c(dS, dI, dR)))
  })
}

# Initial conditions
N <- 1000        # Total population
I <- 10          # Initial number of infected
S <- N - I0     # Initial number of susceptible
R <- 0          # Initial number of recovered
initial_state <- c(S = S, I = I, R = R)

# Parameters
R0 <- 2
beta <- (mu * R0) / N     # Infection rate
mu <- 1/2     # Recovery rate
parameters <- c(beta = beta, mu = mu, N = N)

# Time points
times <- seq(0, 100, by = 1)

# Solve the differential equations
sir_output <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)

# Convert output to a data frame
sir_output_df <- as.data.frame(sir_output)

# Plot the results
ggplot(sir_output_df, aes(x = time)) +
  geom_line(aes(y = I, color = "Infected")) +
  labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
  theme_minimal()


ggplot() +
  # Stochastic mean trajectory
  geom_line(data = mean_trajectory, aes(x = time_bin, y = mean_I, color = "Stochastic Mean"), size = 1.2) +
  
  # Deterministic trajectory
  geom_line(data = sir_output_df, aes(x = time, y = I, color = "Deterministic"), size = 1.2, linetype = "dashed") +
  
  labs(x = "Time", y = "Number of Infected Individuals", color = "Trajectory Type") +
  theme_minimal()

