#RMsim, RMsims, RMsimcheck, RMsimearlyabort, RMTMUnif all for old method of matching 
# - probabilities not comparable to ML method as exact same data in same order not seen
# so one could just get 'lucky' early on based on data seen first

set.seed(42)

# Parameters
R <- 1.5  # Reproduction number
timesteps <- 200  # Total time steps (weeks)
w <- c(0.2, 0.3, 0.3, 0.1, 0.1)  # Serial interval weights (weeks 1–5)
target <- c(1, 2, 0)  # Target pattern to match

# Function to simulate one trajectory
simulate_renewal <- function(R, timesteps, w) {
  I <- numeric(timesteps)
  I[1] <- 1  # Initial number of cases
  
  for (t in 2:timesteps) {
    lambda <- sum(w[1:min(t - 1, length(w))] * rev(I[(t - min(t - 1, length(w))):(t - 1)]))
    I[t] <- rpois(1, R * lambda)
  }
  
  return(I)
}

# Run simulations until 1000 matches are found
get_matching_simulations <- function(n_matches, R, timesteps, w, target) {
  matches <- list()
  count <- 0
  total_runs <- 0
  
  while (count < n_matches) {
    sim <- simulate_renewal(R, timesteps, w)
    total_runs <- total_runs + 1
    
    if (all(sim[1:3] == target)) {
      count <- count + 1
      matches[[count]] <- sim
    }
    
    if (total_runs %% 100 == 0) cat("Checked", total_runs, "simulations, found", count, "matches\n")
  }
  
  cat("Finished. Total simulations run:", total_runs, "\n")
  return(matches)
}

# Run and store 1000 matching simulations
matching_sims <- get_matching_simulations(n_matches = 100, R = R, timesteps = timesteps, w = w, target = c(1, 2, 0))
