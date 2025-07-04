#RMsim, RMsims, RMsimcheck, RMsimearlyabort, RMTMUnif all for old method of matching 
# - probabilities not comparable to ML method as exact same data in same order not seen
# so one could just get 'lucky' early on based on data seen first

set.seed(42)

# Parameters
R <- 1.5
timesteps <- 200
w <- c(0.2, 0.3, 0.3, 0.1, 0.1)
target <- c(1, 2, 0)
threshold <- 200  # Major outbreak threshold

# Function to simulate one trajectory
simulate_renewal <- function(R, timesteps, w) {
  I <- numeric(timesteps)
  I[1] <- 1
  
  for (t in 2:timesteps) {
    lambda <- sum(w[1:min(t - 1, length(w))] * rev(I[(t - min(t - 1, length(w))):(t - 1)]))
    I[t] <- rpois(1, R * lambda)
  }
  
  return(I)
}

# Function to run and collect simulations with target match and label outbreaks
get_matching_simulations <- function(n_matches, R, timesteps, w, target, threshold) {
  matches <- list()
  outbreak_labels <- logical(n_matches)
  count <- 0
  total_runs <- 0
  
  while (count < n_matches) {
    sim <- simulate_renewal(R, timesteps, w)
    total_runs <- total_runs + 1
    
    if (all(sim[1:3] == target)) {
      count <- count + 1
      matches[[count]] <- sim
      outbreak_labels[count] <- sum(sim) > threshold
    }
    
    if (total_runs %% 5000 == 0) cat("Checked", total_runs, "simulations, found", count, "matches\n")
  }
  
  cat("Finished. Total simulations run:", total_runs, "\n")
  
  return(list(trajectories = matches, major_outbreak = outbreak_labels))
}

# Run the simulation matching
result <- get_matching_simulations(
  n_matches = 100,
  R = R,
  timesteps = timesteps,
  w = w,
  target = target,
  threshold = threshold
)

# Example: check proportion of major outbreaks
mean(result$major_outbreak)  # Proportion of matches leading to major outbreak
