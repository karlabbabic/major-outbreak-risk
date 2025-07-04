set.seed(123)  

# Parameters
N <- 1000       
R0 <- 2         
mu <- 0.2       # Recovery rate
beta <- (mu * R0) / N
n_sims <- 10000 # Number of simulations
threshold <- 50 # Outbreak definition: 50 total infections

# Function to simulate one SIR trajectory and return 3-day incidence + outbreak flag
simulate_sir <- function(N, beta, mu, threshold = 50) {
  S <- N - 1
  I <- 1
  R <- 0
  t <- 0
  cum_inf <- 1
  results <- data.frame(t = 0, S = S, I = I, R = R)
  
  while (I > 0 && cum_inf < threshold) {
    rate <- beta * I * S + mu * I
    if (rate == 0) break
    
    t <- t - log(runif(1)) / rate
    
    if (runif(1) < beta * I * S / rate) {
      S <- S - 1
      I <- I + 1
      cum_inf <- cum_inf + 1
    } else {
      I <- I - 1
      R <- R + 1
    }
    
    results <- rbind(results, data.frame(t = t, S = S, I = I, R = R))
  }
  
  # Compute new infections at each step, counting the initial one at t = 0
  results$new_inf <- c(1, diff(N - results$S))
  results$day <- floor(results$t)
  daily_new_inf <- tapply(results$new_inf, results$day, sum)
  
  # Extract first 3 days' new infections
  inf_0 <- ifelse("0" %in% names(daily_new_inf), daily_new_inf["0"], 0)
  inf_1 <- ifelse("1" %in% names(daily_new_inf), daily_new_inf["1"], 0)
  inf_2 <- ifelse("2" %in% names(daily_new_inf), daily_new_inf["2"], 0)
  
  outbreak <- as.integer(cum_inf >= threshold)
  
  return(c(day0 = inf_0, day1 = inf_1, day2 = inf_2, outbreak = outbreak))
}

# Run all sims
sim_data <- replicate(n_sims, simulate_sir(N, beta, mu, threshold), simplify = TRUE)
sim_df <- as.data.frame(t(sim_data))

write.csv(sim_df, "data_matrix_incidence_KB.csv", row.names = FALSE)

# Preview
head(sim_df)

