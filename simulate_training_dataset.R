library(data.table)

# Settings

max_n <- 30000000
store_until <- 3
threshold <- 100
max_time <- 1000

# parameter for poisson prior on I_1
lambda_I1 <- 1

seed <- 123

output_dir <- "data"
output_file <- file.path(output_dir, "training_data_pois_1_30mil.rds")

w <- readRDS("data/serial_interval_weights.rds")

# Renewal–Poisson simulation function
sim_check <- function(R, w, I1, threshold = 100, max_time = 1000, store_until = 3) {
  I <- integer(max_time)
  I[1] <- I1
  
  I_store <- integer(store_until)
  I_store[1] <- I1
  
  zero_run <- 0
  serial_len <- length(w)
  cum_cases <- I1
  
  if (I1 >= threshold) {
    return(list(I_early = I_store, outcome = "major_outbreak"))
  }
  
  for (t in 2:max_time) {
    lambda_t <- 0
    max_s <- min(t - 1, serial_len)
    
    for (s in 1:max_s) {
      lambda_t <- lambda_t + w[s] * I[t - s]
    }
    
    I[t] <- rpois(1, R * lambda_t)
    cum_cases <- cum_cases + I[t]
    
    if (t <= store_until) {
      I_store[t] <- I[t]
    }
    
    if (I[t] == 0L) {
      zero_run <- zero_run + 1
      if (zero_run >= serial_len) {
        return(list(I_early = I_store, outcome = "extinction"))
      }
    } else {
      zero_run <- 0
    }
    
    if (cum_cases >= threshold) {
      return(list(I_early = I_store, outcome = "major_outbreak"))
    }
  }
  # safety
  list(I_early = I_store, outcome = "undetermined")
}

# Poisson draw
pois_prior <- function(lambda) {
  x <- 0
  while (x == 0) {
    x <- rpois(1, lambda)
  }
  x
}

# R_0 uniform dist [0,10]
R_vals <- runif(max_n, 0, 10)

# storage
I_mat <- matrix(0L, nrow = max_n, ncol = store_until)
outcomes <- character(max_n)
I1_vals <- integer(max_n)

cat("Simulating dataset of size:", format(max_n, big.mark = ","), "\n")
cat("Initial cases drawn from Poisson with lambda =", lambda_I1, "\n")

# Main simulation loop
for (i in seq_len(max_n)) {
  set.seed(seed + i)
  
  # For a point-mass prior at I1 = 1, use:
  # I1_i <- 1
  # For a Poisson prior on I1, use:
  I1_i <- pois_prior(lambda_I1)
  
  sim <- sim_check(
    R = R_vals[i],
    w = w,
    I1 = I1_i,
    threshold = threshold,
    max_time = max_time,
    store_until = store_until
  )
  
  I1_vals[i] <- I1_i
  I_mat[i, ] <- sim$I_early
  outcomes[i] <- sim$outcome
  
  if (i %% 100000 == 0) {
    cat("Simulated", format(i, big.mark = ","), "/", format(max_n, big.mark = ","), "\n")
    flush.console()
  }
}

# Build dataset
dt <- as.data.table(I_mat)
setnames(dt, paste0("I", 1:store_until))

dt[, R := R_vals]
dt[, outcome := outcomes]
dt[, outcome_bin := factor(ifelse(outcome == "major_outbreak", 1, 0), levels = c(0, 1))]

cat("Outcome counts:\n")
print(table(dt$outcome))

# Save
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(dt, output_file)

cat("Saved to:", output_file, "\n")
