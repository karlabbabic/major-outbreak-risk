#RMsim, RMsims, RMsimcheck RMsimearlyabort, RMTMUnif all for old method of matching 
# - probabilities not comparable to ML method as exact same data in same order not seen
# so one could just get 'lucky' early on based on data seen first

set.seed(42)

# Parameters

timesteps <- 200
w <- c(0.2, 0.3, 0.3, 0.1, 0.1)
target <- c(1, 3, 2)
threshold <- 200  # Major outbreak threshold

# Simulate only first 3 weeks
simulate_first3 <- function(R, w) {
  I <- numeric(3)
  I[1] <- 1
  
  for (t in 2:3) {
    lambda <- sum(w[1:min(t - 1, length(w))] * rev(I[(t - min(t - 1, length(w))):(t - 1)]))
    I[t] <- rpois(1, R * lambda)
  }
  
  return(I)
}

# Full simulation from t = 4 onwards
simulate_remaining <- function(R, timesteps, w, I) {
  for (t in 4:timesteps) {
    lambda <- sum(w[1:min(t - 1, length(w))] * rev(I[(t - min(t - 1, length(w))):(t - 1)]))
    I[t] <- rpois(1, R * lambda)
  }
  return(I)
}

# Main wrapper with random R per attempt
get_matching_simulations_randomR <- function(n_matches, timesteps, w, target, threshold) {
  matches <- list()
  outbreak_labels <- logical(n_matches)
  used_Rs <- numeric(n_matches)
  attempt_indices <- numeric(n_matches)  # <- NEW: tracks how many sims had been run at each match
  count <- 0
  total_runs <- 0
  
  while (count < n_matches) {
    R <- runif(1, 1, 2)  # Draw R from Uniform(1, 2)
    I3 <- simulate_first3(R, w)
    total_runs <- total_runs + 1
    
    if (all(I3 == target)) {
      I_full <- numeric(timesteps)
      I_full[1:3] <- I3
      I_full <- simulate_remaining(R, timesteps, w, I_full)
      
      count <- count + 1
      matches[[count]] <- I_full
      outbreak_labels[count] <- sum(I_full) > threshold
      used_Rs[count] <- R
      attempt_indices[count] <- total_runs  # <- LOG the attempt count at time of match
    }
    
    if (total_runs %% 1000 == 0) cat("Checked", total_runs, "simulations, found", count, "matches\n")
  }
  
  cat("Finished. Total simulations run:", total_runs, "\n")
  
  return(list(
    trajectories = matches,
    major_outbreak = outbreak_labels,
    R_values = used_Rs,
    attempt_indices = attempt_indices  # <- Include this in the returned list
  ))
}

# mean(result$major_outbreak)  # Proportion of major outbreaks

target <- c(1, 2, 0)
set.seed(42)
result <- get_matching_simulations_randomR(
  n_matches = 5000,
  timesteps = 200,
  w = w,
  target = target,
  threshold = 200
)

#tm_attempts <- c(50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000, 40000, 50000)
tm_labels <- result$major_outbreak
tm_indices <- result$attempt_indices  # how many sims had been attempted at each match

tm_pmos <- sapply(tm_indices, function(n) {
  matched_up_to_n <- which(tm_indices <= n)
  if (length(matched_up_to_n) == 0) return(NA)
  mean(tm_labels[matched_up_to_n])
})
tm_pmos

