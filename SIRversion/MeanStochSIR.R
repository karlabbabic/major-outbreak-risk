library(ggplot2)
library(dplyr)

set.seed(11)

# Parameters
N <- 1000
R0 <- 2
mu <- 1/2
beta <- (mu * R0) / N
num_simulations <- 1000

# Create an empty dataframe to store results
results1 <- data.frame(t = numeric(0), S = numeric(0), I = numeric(0), R = numeric(0), Simulation = integer(0))

for (sim in 1:num_simulations) {
  
  # Initial conditions for each simulation
  S <- N - 1
  I <- 10
  R <- 0
  t <- 0
  
  # Store initial conditions
  sim_results <- data.frame(t = t, S = S, I = I, R = R, Simulation = sim)
  
  # Simulation loop
  while (I > 0) {
    
    # Generate random numbers
    r1 <- runif(1)
    r2 <- runif(1)
    
    # Time to next event
    rate <- beta * I * S + mu * I
    if (rate == 0) break
    
    t <- t - log(r1) / rate
    
    # Determine type of event
    if (r2 < beta * I * S / rate) {
      # Infection event
      S <- S - 1
      I <- I + 1
    } else {
      # Recovery event
      I <- I - 1
      R <- R + 1
    }
    
    # Record results
    sim_results <- rbind(sim_results, data.frame(t = t, S = S, I = I, R = R, Simulation = sim))
  }
  
  # Append to overall results
  results1 <- rbind(results1, sim_results)
}

# Plot all simulations
ggplot(results1, aes(x = t, y = I, group = Simulation)) +
  geom_line(alpha = 0.3) +  # Set transparency for better visibility
  labs(x = "Time", y = "Number of Infected Individuals", color = "Simulation") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide legend for clarity

bin_width <- 1
results_binned <- results1 %>%
  mutate(time_bin = floor(t / bin_width) * bin_width) %>%
  group_by(time_bin, Simulation) %>%
  summarise(I = mean(I), .groups = "drop")

# Compute mean trajectory
mean_trajectory <- results_binned %>%
  group_by(time_bin) %>%
  summarise(mean_I = mean(I), .groups = "drop")

ggplot() +
  geom_line(data = results_binned, aes(x = time_bin, y = I, group = Simulation), alpha = 0.2, color = "grey") +
  geom_line(data = mean_trajectory, aes(x = time_bin, y = mean_I), color = "blue", size = 1.2) +
  labs(x = "Time", y = "Number of Infected Individuals") +
  theme_minimal()
