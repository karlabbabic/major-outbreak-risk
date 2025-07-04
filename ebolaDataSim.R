# Simulations for datasets with Ebola ws 3 days + 5 days

set.seed(42)

# Parameters

T <- 200  # Total time steps (weeks)
#w <- c(0.2, 0.3, 0.3, 0.1, 0.1) # Serial interval weights (weeks 1-5)
w  # Serial interval weights for ebola (weeks)
S <- length(w)
n_sim <- 1000000                     # Number of simulations - increase and compare with proportion
threshold <- 200         # Outbreak threshold

# --------------------------
# Step 2: Run simulations
# --------------------------

# Store 3 days + outbreak indicator
summary_mat_ebola_small <- matrix(0, nrow = n_sim, ncol = 4)
colnames(summary_mat_ebola_small) <- c("t1", "t2", "t3", "outbreak_label")

# For each simulation, record the number of cases in the first three timesteps, 
# and whether each outbreak had more or fewer than 200 infections in total (summed over all timesteps)
# Simulation loop

# Simulation loop
for (sim in 1:n_sim) {
  R <- runif(1, min = 1, max = 2)
  I <- numeric(T)
  I[1] <- 1   # Initial case
  
  for (t in 2:T) {
    lambda <- 0
    for (s in 1:min(t - 1, S)) {
      lambda <- lambda + w[s] * I[t - s]
    }
    I[t] <- rpois(1, R * lambda)
  }
  
  #for (t in 2:T) {
  #  lambda <- sum(w[1:min(t - 1, S)] * I[(t - 1):(t - min(t - 1, S))])
  
  #  I[t] <- rpois(1, R * lambda)
  #}
  
  # Store first 3 days and outbreak flag
  summary_mat_ebola_small[sim, 1:3] <- I[1:3]
  summary_mat_ebola_small[sim, 4] <- as.integer(sum(I) > threshold)
}

# summary_mat_prior
# summary_mat_prior[7,]
# summary_mat_ebola
# 
# ->>> #write.csv(sim_df_prior, "data_matrix_incidence_rm.csv", row.names = FALSE)
# count(subset(sim_df_prior, t1 == 1 & t2 ==0 & t3 ==0 & outbreak_label == 1))/count(subset(sim_df_prior, t1 == 1 & t2 ==0 & t3 ==0))
# 
# head(sim_df_prior)

View(summary_mat_ebola_small)

set.seed(42)

# Parameters

T <- 200  # Total time steps (weeks)
#w <- c(0.2, 0.3, 0.3, 0.1, 0.1) # Serial interval weights (weeks 1-5)
w  # Serial interval weights for ebola (weeks 1)
S <- length(w)
n_sim <- 1000000                     # Number of simulations - increase and compare with proportion
threshold <- 200         # Outbreak threshold

# --------------------------
# Step 2: Run simulations
# --------------------------

# Store 5 days + outbreak indicator
summary_mat_ebola_moredays_test <- matrix(0, nrow = n_sim, ncol = 6)
colnames(summary_mat_ebola_moredays_test) <- c("t1", "t2", "t3", "t4", "t5", "outbreak_label")

# For each simulation, record the number of cases in the first three timesteps, 
# and whether each outbreak had more or fewer than 200 infections in total (summed over all timesteps)
# Simulation loop

# Simulation loop
for (sim in 1:n_sim) {
  R <- runif(1, min = 1, max = 2)
  I <- numeric(T)
  I[1] <- 1   # Initial case
  
  for (t in 2:T) {
    lambda <- 0
    for (s in 1:min(t - 1, S)) {
      lambda <- lambda + w[s] * I[t - s]
    }
    I[t] <- rpois(1, R * lambda)
  }
  
  #for (t in 2:T) {
  #  lambda <- sum(w[1:min(t - 1, S)] * I[(t - 1):(t - min(t - 1, S))])
  
  #  I[t] <- rpois(1, R * lambda)
  #}
  
  # Store first 3 days and outbreak flag
  summary_mat_ebola_moredays_test[sim, 1:5] <- I[1:5]
  summary_mat_ebola_moredays_test[sim, 6] <- as.integer(sum(I) > threshold)
}
