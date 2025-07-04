## Serial interval discretisation (based on Ogi-Gittins et al)

# Parameters of the gamma dist (not rounded)
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

k_max <- 50
k_vals <- 2:k_max
wk_vals <- sapply(k_vals, compute_wk)

# threshold <- 1e-14
# wk_vals_trimmed <- ifelse(wk_vals < threshold, 0, wk_vals)
# w1 <- 1 - sum(wk_vals_trimmed)

# Compute w_1 such that sum(w_k) = 1
w1 <- 1 - sum(wk_vals)

# Combine
w <- c(w1, wk_vals)
w
