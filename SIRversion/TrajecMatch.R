library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(11)

# Parameters
N <- 1000
R0 <- 2
mu <- 0.2
beta <- (mu * R0) / N
num_simulations <- 100000
S <- N - 1
I <- 1
R <- 0
t <- 0

# Observed daily infections (example data)
observed_counts <- c(1, 2, 0)
match_days <- c(1, 2, 3)

# Function to simulate and check first three days
simulate_sir_with_early_abort <- function(S, I, R, beta, mu, observed_counts, match_days) {
  
  # Record simulation results
  results <- data.frame(t = t, S = S, I = I, R = R)
  
  # Track cumulative infections at the end of each day
  daily_infections <- rep(0, length(match_days))
  
  while (I > 0) {
    rate <- beta * I * S + mu * I
    if (rate == 0) break
    
    # Time to next event
    t <- t - log(runif(1)) / rate
    
    # Decide event type
    if (runif(1) < beta * I * S / rate) {
      S <- S - 1
      I <- I + 1
    } else {
      I <- I - 1
      R <- R + 1
    }
    
    results <- rbind(results, data.frame(t = t, S = S, I = I, R = R))
  }
  
  ## For early abort:
  
  # Mark infection events where S drops
  results <- results %>%
    mutate(
      infection_event = c(TRUE, diff(S) < 0),  # mark first row as infection
      day = ceiling(t)
    )
  
  # Ensure initial infection is counted on day 1
  results$day[1] <- 1
  
  # Filter just the infection events
  infection_events <- results %>%
    filter(infection_event)

  # Count number of infection events per day
  daily_summaries <- infection_events %>%
    filter(day %in% match_days) %>%
    group_by(day) %>%
    summarise(total_infections = n(), .groups = "drop")
  
  # Fill in missing days with 0s
  full_summary <- data.frame(day = match_days) %>%
    left_join(daily_summaries, by = "day") %>%
    mutate(total_infections = replace_na(total_infections, 0))

  # Early rejection if the trajectory doesn't match the observed data
  if (!all(full_summary$total_infections == observed_counts)) {
    return(NULL)
  }
  
  return(results)
}
  
# Run and collect matched simulations
matches <- list()
count <- 0
  
start_time <- Sys.time()
  
# Run simulation
while (length(matches) < num_simulations) {
  sim <- simulate_sir_with_early_abort(S, I, R,
                                       beta, mu,
                                       observed_counts, match_days)
  # Saves only the matches
  if (!is.null(sim)) {
    matches[[length(matches) + 1]] <- sim
  }
  
  # Print progress
  count <- count + 1
  if (count %% 2 == 0) cat("Simulations tried:", count, " | Matches found:", length(matches), "\n")
}

end_time <- Sys.time()
cat("Total matches found:", length(matches), "\n")
cat("Time taken for test:", end_time - start_time, "\n")

# Calculate proportion of major outbreaks (>= 50 infections)
check_major_outbreak <- function(sim) {
  return(tail(sim$I + sim$R, 1) > 50)
}

major_outbreaks <- sapply(matches, check_major_outbreak)
probability_major_outbreak <- mean(major_outbreaks)
cat("Estimated Probability of Major Outbreak:", probability_major_outbreak, "\n")
