# ---- Major outbreak probability given R ----
compute_major_outbreak_prob <- function(R, I1, I2, I3, I4, I5, w1, w2, w3, w4, tol = 1e-10, max_iter = 1000) {
  # Step 1: Solve q = exp(R * (q - 1)) numerically
  q_old <- 0.5  # initial guess
  for (i in 1:max_iter) {
    q_new <- exp(R * (q_old - 1))
    if (abs(q_new - q_old) < tol) break
    q_old <- q_new
  }
  q <- q_new
  
  # Step 2: Compute the probability of no major outbreak
  term1 <- (exp(R * (1 - w1 - w2 - w3 - w4) * (q - 1)))^I1
  term2 <- (exp(R * (1 - w1 - w2 - w3) * (q - 1)))^I2
  term3 <- (exp(R * (1 - w1 - w2) * (q - 1)))^I3
  term4 <- (exp(R * (1 - w1) * (q - 1)))^I4
  term5 <- q^I5
  
  prob_no_outbreak <- term1 * term2 * term3 * term4 *term5
  
  # Step 3: Compute the probability of a major outbreak
  return(1 - prob_no_outbreak)
}

# ---- Likelihood function (Equation 2) ----
compute_likelihood <- function(R, I1, I2, I3, I4, I5, w) {
  I <- c(I1, I2, I3, I4, I5)
  likelihood <- 1
  
  for (t in 2:5) {
    lambda_t <- 0
    for (s in 1:(t - 1)) {
      lambda_t <- lambda_t + R * I[t - s] * w[s]
    }
    likelihood <- likelihood * (lambda_t^I[t] * exp(-lambda_t) / factorial(I[t]))
    #likelihood <- likelihood * (lambda_t^I[t] * exp(-lambda_t) / exp(lgamma(I[t] + 1))) #(Poisson-based) 
  }
  
  return(likelihood)
}

# ---- Uniform prior over R in [1, 2] ----
prior_R <- function(R) {
  ifelse(R >= 1 & R <= 2, 1, 0)
}

# ---- Unnormalised posterior: likelihood × prior ----
posterior_R <- function(R, I1, I2, I3, I4, I5, w) {
  compute_likelihood(R, I1, I2, I3, I4, I5, w) * prior_R(R)
}

# ---- Bayesian major outbreak probability (Equation 1) ----
compute_bayesian_outbreak_prob <- function(I1, I2, I3, I4, I5, w1, w2, w3, w4) {
  w <- c(w1, w2, w3, w4)

  integrand <- Vectorize(function(R) {
    posterior_R(R, I1, I2, I3, I4, I5, w) *
      compute_major_outbreak_prob(R, I1, I2, I3, I4, I5, w1, w2, w3, w4)
  })
  
  normaliser <- Vectorize(function(R) {
    posterior_R(R, I1, I2, I3, I4, I5, w)
  })
  
  norm_constant <- integrate(normaliser, lower = 1, upper = 2)$value
  result <- integrate(integrand, lower = 1, upper = 2)$value / norm_constant
  # norm_constant <- integrate(normaliser, lower = 1, upper = 2, subdivisions = 1000)$value
  # result <- integrate(integrand, lower = 1, upper = 2, subdivisions = 1000)$value / norm_constant
  
  return(result)
}

# Observed incidence over 5 days
I1 <- 1
I2 <- 3
I3 <- 0
I4 <- 2
I5 <- 1

# Infectiousness weights (must sum to 1)
w1 <- 0.2
w2 <- 0.3
w3 <- 0.3
w4 <- 0.1  # Note: w5 is not used because only 4 weights are relevant to 5 time points

# Ebola weights
# w1 <- w[1]
# w2 <- w[2]

# Optional: check weight sum
# if (abs(w1 + w2 + w3 + w4 - 1) > 1e-8) {
#   warning("Infectiousness weights do not sum to 1.")
# }

# Compute Bayesian probability of a major outbreak
bayes_prob <- compute_bayesian_outbreak_prob(I1, I2, I3, I4, I5, w1, w2, w3, w4)
cat("Bayesian probability of a major outbreak:", round(bayes_prob, 4), "\n")
true_prob <- bayes_prob
