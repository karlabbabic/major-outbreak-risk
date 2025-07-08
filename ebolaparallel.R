# Load libraries
library(parallel)

# Serial interval weights (Ebola)
w <- c(0.353243306197408, 0.306732467550487, 0.185973712389417, 0.090314213909739,
       0.038916037429299, 0.015562681290287, 0.005914097619519, 0.002165705491243,
       0.000771063663235, 0.000268526120092, 0.00009186497592, 0.000030970511807,
       0.000010313741967, 0.000003399029052, 0.000001110190441, 0.000000359793426,
       0.000000115807735, 0.00000003705049, 0.000000011789884, 0.000000003733592,
       0.000000001177206, 3.69713e-10, 1.15695e-10, 3.6086e-11, 1.1222e-11,
       3.48e-12, 1.076e-12, 3.32e-13, 1.02e-13, 3.1e-14, 1e-14, 3e-15, 1e-15)

# Simulation parameters
set.seed(42)
T <- 200                 # Total time steps
S <- length(w)           # Serial interval length
threshold <- 200         # Outbreak threshold
n_sim <- 1000000         # Total simulations
chunk_size <- 100000     # Simulations per batch
n_chunks <- ceiling(n_sim / chunk_size)
cores <- detectCores() - 1  # Leave 1 core free

# Simulation function for a single run
simulate_outbreak <- function(dummy) {
  R <- runif(1, 1, 2)
  I <- numeric(T)
  I[1] <- 1
  for (t in 2:T) {
    lambda <- 0
    for (s in 1:min(t - 1, S)) {
      lambda <- lambda + w[s] * I[t - s]
    }
    I[t] <- rpois(1, R * lambda)
  }
  c(I[1:5], as.integer(sum(I) > threshold))
}

# Loop through chunks
for (chunk in 1:n_chunks) {
  message("Running chunk ", chunk, " of ", n_chunks, " ...")
  chunk_result <- mclapply(1:chunk_size, simulate_outbreak, mc.cores = cores)
  mat <- do.call(rbind, chunk_result)
  colnames(mat) <- c("t1", "t2", "t3", "t4", "t5", "outbreak_label")
  
  # Save results
  saveRDS(mat, paste0("eboladataset_chunk_", chunk, ".RDS"))
  write.csv(mat, paste0("eboladataset_chunk_", chunk, ".csv"), row.names = FALSE)
}
