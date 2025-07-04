# Renewal model one sim initial test

set.seed(42)

# Parameters

R <- 1.5  # Reproduction number
timesteps <- 200  # Total time steps (weeks)
w <- c(0.2, 0.3, 0.3, 0.1, 0.1)  # Serial interval weights (weeks 1-5)

I <- numeric(timesteps)
I[1] <- 1 # Initial number of cases

# Simulate from week 2 onwards
for (t in 2:timesteps) {
  lambda <- 0
  for (s in 1:min(t - 1, 5)) {
    lambda <- lambda + w[s] * I[t - s]  # Prob of a secondary case at t weeks (2, 3, ...)
  }
  I[t] <- rpois(1, R * lambda)          # *R to get mean/expected number of these cases at t
  # Draw from Poisson to get random number of cases at t with that mean
}
I

# Plot
plot(I, type = "o", xlab = "Week", ylab = "Number of cases", main = "Renewal process epidemic simulation")
grid()

