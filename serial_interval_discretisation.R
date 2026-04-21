## Serial interval discretisation (based on Ogi-Gittins et al)

# Parameters of the gamma dist (WHO Ebola Response Team)
mean_si <- 15.3/7
sd_si <- 9.3/7

# Convert mean and sd to shape and rate of gamma
shape_si <- (mean_si / sd_si)^2
rate_si <- mean_si / (sd_si^2)

# Gamma density function
g <- function(u) {
  dgamma(u, shape = shape_si, rate = rate_si)
}

# Compute w_k for k = 2, 3, ..., k_max
compute_wk <- function(k) {
  integrand <- function(u) {
    (1 - abs(u - k)) * g(u)
  }
  integrate(integrand, lower = k - 1, upper = k + 1)$value
}

k_max <- 10
k_vals <- 2:k_max
wk_vals <- sapply(k_vals, compute_wk)

# Compute w_1 such that sum(w_k) = 1
w1 <- 1 - sum(wk_vals)

# Combine
w <- c(w1, wk_vals)
w

dir.create("data", showWarnings = FALSE, recursive = TRUE)
saveRDS(w, file = "data/serial_interval_weights.rds")
