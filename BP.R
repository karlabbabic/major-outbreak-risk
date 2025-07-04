compute_major_outbreak_prob <- function(R, I1, I2, I3, w1, w2, tol = 1e-10, max_iter = 1000) {
  # Step 1: Solve q = exp(R * (q - 1)) numerically
  q_old <- 0.5  # initial guess
  for (i in 1:max_iter) {
    q_new <- exp(R * (q_old - 1))
    if (abs(q_new - q_old) < tol) break
    q_old <- q_new
  }
  q <- q_new
  
  # Step 2: Compute the probability of no major outbreak
  term1 <- (exp(R * (1 - w1 - w2) * (q - 1)))^I1
  term2 <- (exp(R * (1 - w1) * (q - 1)))^I2
  term3 <- q^I3
  
  prob_no_outbreak <- term1 * term2 * term3
  
  # Step 3: Compute the probability of a major outbreak
  prob_major_outbreak <- 1 - prob_no_outbreak
  
  return(prob_major_outbreak)
}

# Example parameters
R <- 1.5
I1 <- 1
I2 <- 0
I3 <- 0
w1 <- 0.2
w2 <- 0.3

prob <- compute_major_outbreak_prob(R, I1, I2, I3, w1, w2)
cat("Probability of a major outbreak:", prob, "\n")
