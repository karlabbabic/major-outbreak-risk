#RMsim, RMsims, RMsimcheck, RMsimearlyabort, RMTMUnif all for old method of matching 
# - probabilities not comparable to ML method as exact same data in same order not seen
# so one could just get 'lucky' early on based on data seen first

set.seed(42)

# Parameters
R <- 1.5
timesteps <- 200
w <- c(0.2, 0.3, 0.3, 0.1, 0.1)
target <- c(1, 2, 0)
threshold <- 200

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

# Full simulation assuming first 3 already matched
simulate_remaining <- function(R, timesteps, w, I) {
  for (t in 4:timesteps) {
    lambda <- sum(w[1:min(t - 1, length(w))] * rev(I[(t - min(t - 1, length(w))):(t - 1)]))
    I[t] <- rpois(1, R * lambda)
  }
  return(I)
}

# Wrapper with early abort logic
get_matching_simulations_fast <- function(n_matches, R, timesteps, w, target, threshold) {
  matches <- list()
  outbreak_labels <- logical(n_matches)
  count <- 0
  total_runs <- 0
  
  while (count < n_matches) {
    I3 <- simulate_first3(R, w)
    total_runs <- total_runs + 1
    
    if (all(I3 == target)) {
      I_full <- numeric(timesteps)
      I_full[1:3] <- I3
      I_full <- simulate_remaining(R, timesteps, w, I_full)
      
      count <- count + 1
      matches[[count]] <- I_full
      outbreak_labels[count] <- sum(I_full) > threshold
    }
    
    if (total_runs %% 1000 == 0) cat("Checked", total_runs, "simulations, found", count, "matches\n")
  }
  
  cat("Finished. Total simulations run:", total_runs, "\n")
  
  return(list(trajectories = matches, major_outbreak = outbreak_labels))
}

target <- c(1, 1, 0)

set.seed(42)

result <- get_matching_simulations_fast(
  n_matches = 12,
  R = R,
  timesteps = timesteps,
  w = w,
  target = target,
  threshold = threshold
)

mean(result$major_outbreak)
