## Stochastic SIR

library(ggplot2)

set.seed(11)

# Initial conditions
N <- 1000
S <- N-1
I <- 1
R <- 0
t <- 0
R0 <- 2
mu <- 1/2
beta <- (mu*R0)/N

results <- data.frame(t = t, S = S, I = I, R = R)

# Simulation loop
while (I > 0){
  
  # Generate random nos. 0-1
  r1 = runif(1)
  r2 = runif(1)
  
  # Time to next event
  t <- t - 1/(beta*I*S + mu*I)*log(r1)
  
  # Determine type of event
  if (r2 < beta*I*S / (beta*I*S + mu*I)) {
    # Infection event
    S = S - 1
    I = I + 1
  }
  else {
    # Recovery event
    I = I - 1
    R = R + 1
  }
  
  # Record results
  results <- rbind(results, data.frame(t = t, S = S, I = I, R = R))
}

results

# Plot
ggplot(results, aes(x = t)) +
  geom_line(aes(y = I, color = "Infected")) +
  labs(x = "time", y = "Number of Individuals", color = "Compartment") +
  theme_minimal()

ggplot(results, aes(x = t)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = I, color = "Infected")) +
  geom_line(aes(y = R, color = "Recovered")) +
  labs(x = "time", y = "Number of Individuals", color = "Compartment") +
  theme_minimal()


# Add event type column
results$event_type <- NA  # Start with NA

for (i in 2:nrow(results)) {
  if (results$S[i] < results$S[i - 1]) {
    results$event_type[i] <- "inf"  # Infection happened
  } else if (results$R[i] > results$R[i - 1]) {
    results$event_type[i] <- "rec"  # Recovery happened
  }
}

print(results)
