## Alt classifer for CIs

library(brms)

sim_df_prior

sim_df_prior$outbreak_label <- as.factor(sim_df_prior$outbreak_label)

set.seed(123)
shuffled_indices <- sample(nrow(sim_df_prior))

sub_indices_test <- shuffled_indices[1:100]  # Cumulative: always the first n rows
sampletest <- sim_df_prior[sub_indices_test, ]

# Fit logistic regression model
fit_bayes <- brm(
  outbreak_label ~ t1 + t2 + t3,
  data = sampletest,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,  # Use parallel cores if available
  seed = 42
)

#saveRDS(fit_bayes, file = 'bayesianNN.rds')

new_data <- data.frame(t1 = 1, t2 = 0, t3 = 0)

# Posterior predictive draws (probability scale)
pred_draws <- posterior_epred(fit_bayes, newdata = new_data)

# Summary statistics
mean_prob <- mean(pred_draws)
mean_prob
lower_CI <- quantile(pred_draws, 0.025)
upper_CI <- quantile(pred_draws, 0.975)

cat("Estimated PMO:", round(mean_prob, 3),
    "\n95% CI:", round(lower_CI, 3), "-", round(upper_CI, 3))
