# --- Libraries ---
library(ranger)
library(ggplot2)
library(dplyr)
library(scales)

# --- Define Functions ---

# Major outbreak probability given R
compute_major_outbreak_prob <- function(R, I1, I2, I3, w1, w2, tol = 1e-10, max_iter = 1000) {
  q_old <- 0.5
  for (i in 1:max_iter) {
    q_new <- exp(R * (q_old - 1))
    if (abs(q_new - q_old) < tol) break
    q_old <- q_new
  }
  q <- q_new
  term1 <- (exp(R * (1 - w1 - w2) * (q - 1)))^I1
  term2 <- (exp(R * (1 - w1) * (q - 1)))^I2
  term3 <- q^I3
  prob_no_outbreak <- term1 * term2 * term3
  return(1 - prob_no_outbreak)
}

# Likelihood function
compute_likelihood <- function(R, I1, I2, I3, w) {
  I <- c(I1, I2, I3)
  likelihood <- 1
  for (t in 2:3) {
    lambda_t <- 0
    for (s in 1:(t - 1)) {
      lambda_t <- lambda_t + R * I[t - s] * w[s]
    }
    likelihood <- likelihood * (lambda_t^I[t] * exp(-lambda_t) / factorial(I[t]))
  }
  return(likelihood)
}

# Prior and posterior
prior_R <- function(R) ifelse(R >= 1 & R <= 2, 1, 0)
posterior_R <- function(R, I1, I2, I3, w) {
  compute_likelihood(R, I1, I2, I3, w) * prior_R(R)
}

# Bayesian major outbreak probability
compute_bayesian_outbreak_prob <- function(I1, I2, I3, w1, w2) {
  w <- c(w1, w2)
  integrand <- function(R_vals) {
    sapply(R_vals, function(R) {
      posterior_R(R, I1, I2, I3, w) * compute_major_outbreak_prob(R, I1, I2, I3, w1, w2)
    })
  }
  normaliser <- function(R_vals) sapply(R_vals, function(R) posterior_R(R, I1, I2, I3, w))
  norm_constant <- integrate(normaliser, lower = 1, upper = 2)$value
  result <- integrate(integrand, lower = 1, upper = 2)$value / norm_constant
  return(result)
}

# Train ML models once at various subset sizes
train_ml_models <- function(sim_df_prior, train_sizes) {
  models <- list()
  set.seed(123)
  for (i in seq_along(train_sizes)) {
    n <- train_sizes[i]
    sub_train <- sim_df_prior[sample(nrow(sim_df_prior), n), ]
    sub_train$outbreak_label <- as.factor(sub_train$outbreak_label)
    rf_model <- ranger(outbreak_label ~ t1 + t2 + t3, data = sub_train, probability = TRUE)
    models[[i]] <- list(size = n, model = rf_model)
  }
  return(models)
}

# Predict with stored models for a specific input [I1, I2, I3]
predict_with_models <- function(models, I1, I2, I3) {
  target_input <- data.frame(t1 = I1, t2 = I2, t3 = I3)
  preds <- sapply(models, function(entry) {
    predict(entry$model, data = target_input)$predictions[, "1"]
  })
  data.frame(sims = sapply(models, function(x) x$size), pmo = preds, method = "ML (Random Forest")
}

# Plot convergence for given input and results
plot_pmo_convergence <- function(I1, I2, I3, w1 = 0.2, w2 = 0.3, sim_df_prior, ml_results) {
  target_input <- data.frame(t1 = I1, t2 = I2, t3 = I3)
  bayes_prob <- compute_bayesian_outbreak_prob(I1, I2, I3, w1, w2)
  
  match_rows <- which(sim_df_prior$t1 == I1 & sim_df_prior$t2 == I2 & sim_df_prior$t3 == I3)
  tm_labels <- sim_df_prior$outbreak_label[match_rows]
  tm_attempts <- seq_along(match_rows)
  tm_pmos <- cumsum(as.numeric(as.character(tm_labels))) / tm_attempts
  tm_df <- data.frame(sims = tm_attempts, pmo = tm_pmos, method = "Trajectory matching")
  
  all_df <- bind_rows(ml_results, tm_df)
  lower_bound <- bayes_prob * 0.95
  upper_bound <- bayes_prob * 1.05
  
  ggplot(all_df, aes(x = sims, y = pmo, colour = method, shape = method)) +
    annotate("rect", xmin = min(all_df$sims), xmax = max(all_df$sims),
             ymin = lower_bound, ymax = upper_bound, alpha = 0.2, fill = "red") +
    geom_hline(yintercept = bayes_prob, linetype = "dashed", colour = "red", linewidth = 1) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_x_log10(
      breaks = c(50, 100, 200, 500, 1000, 2000, 5000, 1e4, 2e4, 5e4),
      labels = label_comma()
    ) +
    scale_y_continuous(limits = c(0, 1.05)) +
    scale_color_manual(values = c("ML (Random Forest)" = "forestgreen", "Trajectory matching" = "blue")) +
    scale_shape_manual(values = c("ML (Random Forest)" = 16, "Trajectory matching" = 17)) +
    labs(
      title = "Convergence of PMO estimates",
      subtitle = paste0("Estimated probability of major outbreak for [", I1, ", ", I2, ", ", I3, "]"),
      x = "Number of simulations (log scale)",
      y = "Estimated PMO"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Define training sizes and train models ONCE
train_sizes <- c(50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000)
ml_models <- train_ml_models(sim_df_prior, train_sizes)

# For a specific early pattern
ml_results <- predict_with_models(ml_models, I1 = 1, I2 = 2, I3 = 0)

# Define target early pattern
target_input <- data.frame(t1 = 1, t2 = 2, t3 = 0)

# Define training sizes to test
n_vals <- c(50, 100, 200, 350, 400, 450, 500, 700, 750, 800, 850, 900, 1000, 1500, 2000, 3000, 5000, 5500, 6000, 6500, 7000, 7500, 10000, 20000, 30000, 40000, 50000)
ml_preds <- numeric(length(n_vals))

set.seed(123)  # For reproducibility

# Loop through training sizes
for (i in seq_along(n_vals)) {
  set.seed(123)
  n <- n_vals[i]
  
  # Sample n rows from the full dataset
  sub_train <- sim_df_prior[sample(nrow(sim_df_prior), n), ]
  sub_train$outbreak_label <- as.factor(sub_train$outbreak_label)
  
  # Train random forest model
  rf_model <- ranger(
    outbreak_label ~ t1 + t2 + t3,
    data = sub_train,
    probability = TRUE,
    mtry = 3
  )
  
  # Predict probability of major outbreak for target pattern
  ml_preds[i] <- predict(rf_model, data = target_input)$predictions[, "1"]
  
  cat(sprintf("Trained on %d sims → Predicted PMO: %.4f\n", n, ml_preds[i]))
}


ml_df <- data.frame(
  sims = n_vals,
  pmo = ml_preds,
  method = "ML (Random Forest)"
)

# Generate convergence plot
plot_pmo_convergence(1, 2, 0, sim_df_prior = sim_df_prior, ml_results = ml_results)

