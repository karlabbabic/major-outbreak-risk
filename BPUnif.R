# Analytic PMO calculation for candidate trajectories

library(data.table)

# Load serial interval
w <- readRDS("data/serial_interval_weights.rds")

# Major outbreak probability given R
compute_major_outbreak_prob <- function(R, I, w, tol = 1e-10, max_iter = 1000) {
  I <- as.numeric(I)
  T <- length(I)
  
  if (T > 1) {
    w <- w[seq_len(T - 1)]
  } else {
    w <- numeric(0)
  }
  
  q_old <- 0.5
  q_new <- q_old
  
  for (i in seq_len(max_iter)) {
    q_new <- exp(R * (q_old - 1))
    if (abs(q_new - q_old) < tol) break
    q_old <- q_new
  }
  
  q <- q_new
  
  prob_no_outbreak <- 1
  
  for (t in seq_len(T)) {
    if (t == T) {
      term <- q^I[t]
    } else {
      w_slice <- w[1:(T - t)]
      term <- (exp(R * (1 - sum(w_slice)) * (q - 1)))^I[t]
    }
    prob_no_outbreak <- prob_no_outbreak * term
  }
  
  1 - prob_no_outbreak
}

# Likelihood function
compute_likelihood <- function(R, I, w) {
  I <- as.numeric(I)
  T <- length(I)
  
  if (T > 1) {
    w <- w[seq_len(T - 1)]
  } else {
    w <- numeric(0)
  }
  
  likelihood <- 1
  
  if (T >= 2) {
    for (t in 2:T) {
      lambda_t <- 0
      for (s in 1:(t - 1)) {
        lambda_t <- lambda_t + R * I[t - s] * w[s]
      }
      likelihood <- likelihood * dpois(I[t], lambda = lambda_t)
    }
  }
  
  likelihood
}

# Prior on R
prior_R <- function(L = 0, U = 10) {
  list(
    pdf = function(R) ifelse(R >= L & R <= U, 1 / (U - L), 0),
    support = c(L, U)
  )
}

# Unnormalised posterior
posterior_R <- function(R, I, w, prior) {
  compute_likelihood(R, I, w) * prior$pdf(R)
}

# Full analytic major outbreak probability
compute_analytic_outbreak_prob <- function(I, w, prior) {
  I <- as.numeric(I)
  stopifnot(length(w) >= length(I) - 1)
  
  if (length(I) > 1) {
    w_use <- w[seq_len(length(I) - 1)]
  } else {
    w_use <- numeric(0)
  }
  
  a <- prior$support[1]
  b <- prior$support[2]
  
  integrand <- Vectorize(function(R) {
    posterior_R(R, I, w_use, prior) *
      compute_major_outbreak_prob(R, I, w_use)
  })
  
  normaliser <- Vectorize(function(R) {
    posterior_R(R, I, w_use, prior)
  })
  
  norm_constant <- integrate(
    normaliser,
    lower = a,
    upper = b,
    rel.tol = 1e-6
  )$value
  
  result <- integrate(
    integrand,
    lower = a,
    upper = b,
    rel.tol = 1e-6
  )$value / norm_constant
  
  result
}

# example use
# 
# compute_analytic_outbreak_prob(I = c(1,2,0),
#                                w = w,
#                                prior = prior_R(0,10))

# ----------------------------------------

# For downstream underreporting analysis
# Compute the analytic PMO over the full grid of candidate latent incidence sequences

generate_candidate_sequences <- function(Y, buffer, I1 = NULL, fix_I1 = FALSE) {
  T <- length(Y)
  I_ranges <- vector("list", T)
  
  for (t in seq_len(T)) {
    I_ranges[[t]] <- Y[t]:(Y[t] + buffer)
  }
  
  # e.g. when prior on I_1 is a point mass with P(I_1=1)=1
  if (fix_I1) {
    if (is.null(I1)) stop("If fix_I1 = TRUE, provide I1.")
    I_ranges[[1]] <- I1
  }
  
  grid <- expand.grid(I_ranges)
  colnames(grid) <- paste0("I", seq_len(T))
  grid
}

# Settings
Y_obs <- c(1, 0, 0)
buffer <- 25

prior_lower <- 0
prior_upper <- 10

out_dir <- "pmo_outputs"
out_file_csv <- file.path(out_dir, "analytic_pmo_3_weeks.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Build candidate grid
grid <- as.data.table(
  generate_candidate_sequences(
    Y = Y_obs,
    buffer = buffer,
    I1 = NULL, 
    fix_I1 = FALSE
  )
)

cat("Candidate sequences:", nrow(grid), "\n")

prior <- prior_R(prior_lower, prior_upper)

grid[, PMO_analytic := apply(.SD, 1, function(Iseq) {
  compute_analytic_outbreak_prob(
    I = as.numeric(Iseq),
    w = w,
    prior = prior
  )
}), .SDcols = c("I1", "I2", "I3")]


# Save
fwrite(grid, out_file_csv)

cat("Saved analytic PMO grid to:\n")
cat(" -", out_file_csv, "\n")
