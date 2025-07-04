library(ggplot2)
library(scales)
library(dplyr)

# Combine all methods into one data frame
all_df <- rbind(ml_df, tm_df, xgb_df)

# Define ±5% tolerance band around true PMO
lower_bound <- true_prob * 0.95
upper_bound <- true_prob * 1.05

# --- Convergence points ---

# Random Forest
in_band_rf <- ml_preds >= lower_bound & ml_preds <= upper_bound
first_consistent_index_rf <- NA
for (i in seq_along(in_band_rf)) {
  if (all(in_band_rf[i:length(in_band_rf)])) {
    first_consistent_index_rf <- i
    break
  }
}
critical_n_rf <- if (!is.na(first_consistent_index_rf)) n_vals[first_consistent_index_rf] else NA

# XGBoost
in_band_xgb <- xgb_preds >= lower_bound & xgb_preds <= upper_bound
first_consistent_index_xgb <- NA
for (i in seq_along(in_band_xgb)) {
  if (all(in_band_xgb[i:length(in_band_xgb)])) {
    first_consistent_index_xgb <- i
    break
  }
}
critical_n_xgb <- if (!is.na(first_consistent_index_xgb)) n_vals[first_consistent_index_xgb] else NA

# Trajectory Matching
in_band_tm <- tm_pmos >= lower_bound & tm_pmos <= upper_bound
first_consistent_index_tm <- NA
for (i in seq_along(in_band_tm)) {
  if (all(in_band_tm[i:length(in_band_tm)])) {
    first_consistent_index_tm <- i
    break
  }
}
critical_n_tm <- if (!is.na(first_consistent_index_tm)) attempt_indices[first_consistent_index_tm] else NA

# --- Plot ---
ggplot(all_df, aes(x = sims, y = pmo, colour = method, shape = method)) +
  
  # Tolerance band
  annotate("rect",
           xmin = min(all_df$sims), xmax = max(all_df$sims),
           ymin = lower_bound, ymax = upper_bound,
           alpha = 0.2, fill = "red") +
  
  # True Bayesian PMO line
  geom_hline(yintercept = true_prob, linetype = "dashed", colour = "red", size = 1) +
  
  # Convergence lines and annotations
  {if (!is.na(critical_n_rf)) geom_vline(xintercept = critical_n_rf, linetype = "dotted", colour = "forestgreen", size = 0.8)} +
  {if (!is.na(critical_n_rf)) annotate("text", x = critical_n_rf, y = 0.25, label = paste0("n = ", critical_n_rf), angle = 0, vjust = -0.5, hjust = 1, size = 4, colour = "forestgreen")} +
  
  {if (!is.na(critical_n_xgb)) geom_vline(xintercept = critical_n_xgb, linetype = "dotted", colour = "darkorange", size = 0.8)} +
  {if (!is.na(critical_n_xgb)) annotate("text", x = critical_n_xgb, y = 0.35, label = paste0("n = ", critical_n_xgb), angle = 0, vjust = -0.5, hjust = 1, size = 4, colour = "darkorange")} +
  
  {if (!is.na(critical_n_tm)) geom_vline(xintercept = critical_n_tm, linetype = "dotted", colour = "blue", size = 0.8)} +
  {if (!is.na(critical_n_tm)) annotate("text", x = critical_n_tm, y = 0.15, label = paste0("n = ", critical_n_tm), angle = 0, vjust = -0.5, hjust = 1, size = 4, colour = "blue")} +
  
  # Method curves
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  # Axes
  scale_x_log10(labels = comma) +
  scale_y_continuous(limits = c(0, 1.05)) +
  
  # Aesthetic mapping
  scale_color_manual(values = c(
    "ML (Random Forest)" = "forestgreen",
    "ML (XGBoost)" = "darkorange",
    "Trajectory matching" = "blue"
  )) +
  scale_shape_manual(values = c(
    "ML (Random Forest)" = 16,
    "ML (XGBoost)" = 15,
    "Trajectory matching" = 17
  )) +
  
  labs(
    title = "Convergence of PMO estimates",
    subtitle = "Estimated probability of major outbreak for [1, 1, 2]",
    x = "Number of simulations (log scale)",
    y = "Estimated PMO",
    colour = "Estimation method",
    shape = "Estimation method"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
