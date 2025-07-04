### Script for pred outputs for x3 methods -> plot (could turn into a singular plot fn)
                                                    # (see 'Plots.R' for inspo code)

library(scales)

#### BP ######
# Observed incidence over 3 days
I1 <- 1
I2 <- 2
I3 <- 0

bayes_prob <- compute_bayesian_outbreak_prob(I1, I2, I3, w1, w2)
cat("Bayesian probability of a major outbreak:", round(bayes_prob, 4), "\n")
true_prob <- bayes_prob

#### ML preds ####
set.seed(42)
# Define a new input pattern
target_input <- data.frame(t1 = 1, t2 = 0, t3 = 0)
ml_preds <- predict_pmo(target_input, rf_models_n2, n_vals)
print(ml_preds)

### TM preds ########

target <- c(1, 2, 2)
match_rows <- c()
attempt_indices <- c()

for (i in seq_along(shuffled_indices)) {
  idx <- shuffled_indices[i]
  sim <- sim_df_prior[idx, ]
  
  if (sim$t1 == target[1] && sim$t2 == target[2] && sim$t3 == target[3]) {
    match_rows <- c(match_rows, idx)
    attempt_indices <- c(attempt_indices, i)  # Simulation index in shuffled order
  }
}

tm_labels <- as.numeric(as.character(sim_df_prior$outbreak_label[match_rows]))
tm_pmos <- cumsum(tm_labels) / seq_along(tm_labels)
tm_df <- data.frame(sims = attempt_indices, pmo = tm_pmos, method = "Trajectory matching")

# plot
ml_df <- data.frame(
  sims = n_vals,
  pmo = ml_preds,
  method = "ML (Random Forest)"
)

all_df <- rbind(ml_df, tm_df)

# Define tolerance band
lower_bound <- true_prob * 0.95
upper_bound <- true_prob * 1.05

# Logical vector: is each prediction within the band?
in_band <- ml_preds >= lower_bound & ml_preds <= upper_bound

# Find the first index where all *following* values remain within the band
first_consistent_index <- NA
for (i in seq_along(in_band)) {
  if (all(in_band[i:length(in_band)])) {
    first_consistent_index <- i
    break
  }
}
# Retrieve corresponding training size
critical_n_ml <- if (!is.na(first_consistent_index)) n_vals[first_consistent_index] else NA

# Find first point where TM enters the band
### For TM predictions
# Logical vector for RM predictions within the band
in_band_tm <- tm_pmos >= lower_bound & tm_pmos <= upper_bound

# Find first index where all following TM values are within band
first_consistent_index_tm <- NA
for (i in seq_along(in_band_tm)) {
  if (all(in_band_tm[i:length(in_band_tm)])) {
    first_consistent_index_tm <- i
    break
  }
}
# Retrieve corresponding simulation size
critical_n_tm <- if (!is.na(first_consistent_index_tm)) attempt_indices[first_consistent_index_tm] else NA

ggplot(all_df, aes(x = sims, y = pmo, colour = method, shape = method)) +
  # ±5% band
  annotate("rect",
           xmin = min(all_df$sims), xmax = max(all_df$sims),
           ymin = lower_bound, ymax = upper_bound,
           alpha = 0.2, fill = "red") +
  
  # True PMO line
  geom_hline(yintercept = true_prob, linetype = "dashed", colour = "red", size = 1) +
  
  # Vertical lines for critical convergence (optional)
  geom_vline(xintercept = critical_n_ml, linetype = "dotted", colour = "forestgreen", size = 0.8) +
  annotate("text", x = critical_n_ml, y = 0.3, label = paste0("n = ", critical_n_ml),
           angle = 0, vjust = -0.5, hjust = 1, size = 4, colour = "forestgreen") +
  geom_vline(xintercept = critical_n_tm, linetype = "dotted", colour = "blue", size = .8) +
  annotate("text", x = critical_n_tm, y = .2, label = paste0("n = ", critical_n_tm),
           angle = 0, vjust = -0.5, hjust = 1, size = 4, colour = "blue") +

  # # Add ribbon for trajectory matching CI
  # geom_ribbon(data = tm_df, aes(ymin = lower, ymax = upper, fill = method),
  #             alpha = 0.2, colour = NA) +

  # Main data curves
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  scale_x_log10(labels = comma) +
  #scale_x_log10(labels = comma, limits = c(50, max(all_df$sims)))
  scale_y_continuous(limits = c(0, 1.05)) +
  
  labs(
    title = "Convergence of PMO estimates",
    subtitle = "Estimated probability of major outbreak for [1, 2, 0]",
    x = "Number of simulations (log scale)",
    y = "Estimated PMO",
    colour = "Estimation method",
    shape = "Estimation method"
  ) +
  
  scale_color_manual(values = c(
    "ML (Random Forest)" = "forestgreen",
    "Trajectory matching" = "blue"
  )) +
  scale_shape_manual(values = c(
    "ML (Random Forest)" = 16,      # solid circle
    "Trajectory matching" = 17      # solid triangle
  )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
