library(parallel)

set.seed(11)

# Parameters
N <- 1000
R0 <- 2
mu <- 0.2
beta <- (mu * R0) / N
num_matches_needed <- 100000
max_time <- 50
observed_counts <- c(1, 2, 1)
match_days <- c(1, 2, 3)

# Define the simulation function
simulate_until_match <- function(dummy = NULL) {
  repeat {
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
    
    while (I > 0) {
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
    
    infection_indices <- which(c(TRUE, diff(S_vec) < 0))
    infection_days <- ceiling(t_vec[infection_indices])
    infection_days[1] <- 1
    
    counts <- tabulate(infection_days[infection_days %in% match_days], nbins = max(match_days))
    if (length(counts) < length(match_days)) {
      counts <- c(counts, rep(0, length(match_days) - length(counts)))
    }
    
    if (all(counts == observed_counts)) {
      return(I + R)
    }
  }
}

# Detect available cores and set up the cluster
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)

# Export necessary variables and function to each worker
clusterExport(cl, varlist = c("simulate_until_match", "N", "R0", "mu", "beta", 
                              "observed_counts", "match_days"), envir = environment())

# Run simulations in parallel
start_time <- Sys.time()

results <- parLapply(cl, 1:num_matches_needed, simulate_until_match)

stopCluster(cl)
end_time <- Sys.time()

matches <- unlist(results)
prob_major_outbreak <- mean(matches > 50)

cat("Total matches found:", length(matches), "\n")
cat("Time taken:", end_time - start_time, "\n")
cat("Estimated probability of major outbreak:", prob_major_outbreak, "\n")
