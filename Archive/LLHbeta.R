# Log likelihood of beta

log_likelihood_beta <- function(beta, mu, S, I, N, t, event_type) {
  log_lik <- 0
  
  for (i in 2:length(t)) {
    delta_t <- results$t[i] - results$t[i - 1]
    rate_inf <- (beta * results$S[i] * results$I[i])
    rate_rec <- mu * results$I[i]
    
    total_rate <- rate_inf + rate_rec
    
    # Contribution from waiting time
    log_lik <- log_lik - total_rate * delta_t
    
    # Contribution from the event itself
    if (!is.na(event_type[i])) {
      if (event_type[i] == "I") {
        log_lik <- log_lik + log(rate_inf)
      } else if (event_type[i] == "R") {
        log_lik <- log_lik + log(rate_rec)
      }
    }
    
  }
  
  return(log_lik)
}

event_type <- c(NA, NA, "I", "R")

log_likelihood_beta(beta = 0.5, mu, S, I, N, t, event_type)