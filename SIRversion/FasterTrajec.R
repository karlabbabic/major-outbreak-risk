set.seed(11)

# Parameters
N <- 1000
R0 <- 2
mu <- 0.2
beta <- (mu * R0) / N
num_simulations <- 1000    
observed_counts <- c(1, 2, 0)
match_days <- c(1, 2, 3)

matches <- numeric(num_simulations)  # Store only final I + R
match_count <- 0
sim_count <- 0

start_time <- Sys.time()

while (match_count < num_simulations) {
  S <- N - 1
  I <- 1
  R <- 0
  t <- 0
  
  t_vec <- numeric()
  S_vec <- numeric()
  I_vec <- numeric()
  R_vec <- numeric()
  
  t_vec[1] <- t
  S_vec[1] <- S
  I_vec[1] <- I
  R_vec[1] <- R
  
  while (I > 0 && t < max_time) {
    rate <- beta * I * S + mu * I
    if (rate == 0) break
    t <- t - log(runif(1)) / rate
    
    if (runif(1) < beta * I * S / rate) {
      S <- S - 1
      I <- I + 1
    } else {
      I <- I - 1
      R <- R + 1
    }
    
    idx <- length(t_vec) + 1
    t_vec[idx] <- t
    S_vec[idx] <- S
    I_vec[idx] <- I
    R_vec[idx] <- R
  }
  
  # Detect infection events
  infection_indices <- which(c(TRUE, diff(S_vec) < 0))
  infection_days <- ceiling(t_vec[infection_indices])
  infection_days[1] <- 1  # ensure initial infection counts on day 1
  
  # count infections per day
  counts <- tabulate(infection_days[infection_days %in% match_days], nbins = max(match_days))
  if (length(counts) < length(match_days)) {
    counts <- c(counts, rep(0, length(match_days) - length(counts)))
  }
  
  #If match store final outbreak size
  if (all(counts == observed_counts)) {
    match_count <- match_count + 1
    matches[match_count] <- I + R  # store only final value
    
    if (match_count %% 100 == 0) {
      cat("Matches found:", match_count, "| Simulations tried:", sim_count, "\n")
    }
  }
  
  sim_count <- sim_count + 1
}

end_time <- Sys.time()
cat("Total matches found:", match_count, "\n")
cat("Simulations tried:", sim_count, "\n")
cat("Time taken:", end_time - start_time, "\n")

# Final outbreak prob estimate
prob_major_outbreak <- mean(matches > 50)
cat("Estimated probability of major outbreak:", prob_major_outbreak, "\n")
