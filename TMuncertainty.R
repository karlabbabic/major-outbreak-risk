# CI for TM method

library(HDInterval)
library(ggplot2)
library(scales)

# Observed data
k <- 11 
N <- 12 

# Posterior is Beta(k + 1, N - k + 1)
a <- k + 1
b <- N - k + 1

# Equal-tailed 95% credible interval
equal_tailed <- qbeta(c(0.025, 0.975), shape1 = a, shape2 = b)

# Highest Density interval
posterior_samples <- rbeta(1e6, a, b)
hdi_interval <- hdi(posterior_samples, credMass = 0.95)

cat("Equal-tailed 95% CI: [", round(equal_tailed[1], 3), ",", round(equal_tailed[2], 3), "]\n")
cat("Highest Density 95% CI: [", round(hdi_interval[1], 3), ",", round(hdi_interval[2], 3), "]\n")

# Calc TM PMO and 95% equal-tailed credible intervals
tm_labels <- as.numeric(as.character(sim_df_prior$outbreak_label[match_rows]))
tm_cumsum <- cumsum(tm_labels)
tm_counts <- seq_along(tm_labels)

# Compute mean PMO and intervals at each point
tm_pmo <- tm_cumsum / tm_counts

# Equal-tailed 95% credible intervals using Beta posterior
lower_ci <- qbeta(0.025, shape1 = tm_cumsum + 1, shape2 = tm_counts - tm_cumsum + 1)
upper_ci <- qbeta(0.975, shape1 = tm_cumsum + 1, shape2 = tm_counts - tm_cumsum + 1)

# Data frame with CIs
tm_df <- data.frame(
  sims = attempt_indices,
  pmo = tm_pmo,
  lower = lower_ci,
  upper = upper_ci,
  method = "Trajectory matching"
)

ggplot(tm_df, aes(x = sims)) +
  # 95% CI ribbon
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  
  # PMO estimate line
  geom_line(aes(y = pmo), colour = "steelblue", size = 1) +
  
  # CI bounds
  geom_line(aes(y = lower), colour = "steelblue", linetype = "dashed") +
  geom_line(aes(y = upper), colour = "steelblue", linetype = "dashed") +
  geom_hline(yintercept = true_prob, linetype = "dashed", colour = "red", size = 0.7) +
  annotate("rect",
           xmin = min(tm_df$sims), xmax = max(all_df$sims),
           ymin = lower_bound, ymax = upper_bound,
           alpha = 0.2, fill = "red")+
  
  scale_x_log10(labels = comma) +
  scale_y_continuous(limits = c(0, 1.05)) +
  
  labs(
    title = "Trajectory Matching Uncertainty",
    subtitle = "Estimated PMO with 95% credible intervals for [1, 1, 2]",
    x = "Number of matching simulations (log scale)",
    y = "Estimated PMO",
    #caption = "Credible intervals: Beta(k + 1, n − k + 1)"
  ) +
  
  theme_minimal(base_size = 14)
