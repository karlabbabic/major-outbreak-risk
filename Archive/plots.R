## i) Extra early n_vals for adding extra subsets for training
# ii) old plot style 
# iii) and signular plot function trial

library(ranger)
library(ggplot2)
library(dplyr)

###########
# Densify in convergence region
extra_n <- seq(0, 50, by = 10)

# Combine and deduplicate
n_vals_updated <- sort(unique(c(n_vals, extra_n)))

# Extract existing names from rf_models (assumes naming like "n_5000")
existing_n <- as.integer(sub("n_", "", names(rf_models)))
rf_models_n2 <- rf_models
# Identify new training sizes
new_n <- setdiff(n_vals_updated, existing_n)
for (n in new_n) {
  set.seed(123)  # Same seed for consistency
  sub_train <- sim_df_prior[sample(nrow(sim_df_prior), n), ]
  sub_train$outbreak_label <- as.factor(sub_train$outbreak_label)
  
  rf_model <- ranger(
    outbreak_label ~ t1 + t2 + t3,
    data = sub_train,
    probability = TRUE,
    mtry = 3
  )
  
  rf_models_n2[[paste0("n_", n)]] <- rf_model
  cat(sprintf("Trained model for n = %d\n", n))
}
# Reorder rf_models to match the new n_vals_updated
rf_models <- rf_models[order(as.integer(sub("n_", "", names(rf_models))))]
####

early_n <- seq(5, 45, by = 5)
n_vals_updated <- sort(unique(c(n_vals, early_n)))

# Extract already trained sizes from rf_models
existing_n <- as.integer(sub("n_", "", names(rf_models)))

# Find new ones that need training
new_early_n <- setdiff(early_n, existing_n)

for (n in new_early_n) {
  set.seed(123)  # for consistency
  sub_train <- sim_df_prior[sample(nrow(sim_df_prior), n), ]
  sub_train$outbreak_label <- as.factor(sub_train$outbreak_label)
  
  rf_model <- ranger(
    outbreak_label ~ t1 + t2 + t3,
    data = sub_train,
    probability = TRUE,
    mtry = 3
  )
  
  rf_models_n2[[paste0("n_", n)]] <- rf_model
  cat(sprintf("Trained small-n model for n = %d\n", n))
}

rf_models_n2 <- rf_models_n2[order(as.integer(sub("n_", "", names(rf_models_n2))))]

ml_preds <- predict_pmo(target_input, rf_models, as.integer(sub("n_", "", names(rf_models))))

plot(n_vals, ml_preds, type = "b", col = "forestgreen", pch = 16, log = "x",
     ylim = c(0, 1),
     xlab = "Number of training simulations",
     ylab = "Estimated PMO for [1, 2, 0]",
     main = "ML convergence to true outbreak probability")

abline(h = 0.7702, col = "red", lty = 2, lwd = 2)
lines(tm_attempts, tm_pmos, type = "b", col = "blue", pch = 17)

legend("bottomright",
       legend = c("ML estimate", "True (branching process)"),
       col = c("forestgreen", "red"),
       lty = c(1, 2),
       pch = c(16, NA),
       lwd = c(2, 2))

# Your ML and TM convergence plot
plot(n_vals, ml_preds, type = "b", col = "forestgreen", pch = 16, log = "x",
     ylim = c(0, 1),
     xlab = "Number of simulations",
     ylab = "Estimated PMO for [1, 1, 0]",
     main = "Convergence of PMO estimates")

true_prob <- 0.7702  # from BP estimate
true_prob <- 0.6052 
true_prob <- 0.9920523
true_prob <- 0.9997
lower_bound <- true_prob * 0.95
upper_bound <- true_prob * 1.05

# Shaded band for ±5% around the true value
polygon(
  x = c(min(n_vals), max(n_vals), max(n_vals), min(n_vals)),
  y = c(lower_bound, lower_bound, upper_bound, upper_bound),
  col = adjustcolor("red", alpha.f = 0.15),
  border = NA
)

# Add the true PMO line
abline(h = true_prob, col = "red", lty = 2, lwd = 2)

# TM estimates (if added)
lines(tm_indices, tm_indices, type = "b", col = "blue", pch = 17)

# ML again (on top)
lines(n_vals, ml_preds, type = "b", col = "forestgreen", pch = 16)

# Legend
legend("bottomright",
       legend = c("ML estimate", "Trajectory matching", "True (branching process)", "±5% band"),
       col = c("forestgreen", "blue", "red", adjustcolor("red", alpha.f = 0.15)),
       lty = c(1, 1, 2, NA), pch = c(16, 17, NA, 15), pt.cex = c(1, 1, NA, 2),
       border = NA, bty = "n")



plot(n_vals, ml_preds, type = "b", col = "forestgreen", pch = 16, log = "x",
     ylim = c(0, 1),
     xlab = "Number of simulations",
     ylab = "Estimated PMO for [1, 2, 5]",
     main = "Convergence of PMO estimates")

# Add ±5% tolerance band
polygon(
  x = c(min(tm_indices), max(tm_indices), max(tm_indices), min(tm_indices)),
  y = c(lower_bound, lower_bound, upper_bound, upper_bound),
  col = adjustcolor("red", alpha.f = 0.15),
  border = NA
)

# Add true value line
abline(h = true_prob, col = "red", lty = 2, lwd = 2)

# Add smoothed TM curve
lines(tm_indices, tm_pmos, col = "blue", lwd = 2)
points(tm_indices, tm_pmos, col = "blue", pch = 17)

# Redraw ML curve on top
lines(n_vals, ml_preds, type = "b", col = "forestgreen", pch = 16)

# Optional: mark first point ML/TM enters ±5% band
within_band_ml <- which(ml_preds >= lower_bound & ml_preds <= upper_bound)
if (length(within_band_ml)) abline(v = n_vals[min(within_band_ml)], col = "forestgreen", lty = 3)

within_band_tm <- which(tm_pmos >= lower_bound & tm_pmos <= upper_bound)
if (length(within_band_tm)) abline(v = tm_indices[min(within_band_tm)], col = "blue", lty = 3)

# Add legend
legend("bottomright", inset = -.1,
       legend = c("ML estimate", "Trajectory matching", "True (branching process)", "±5% band"),
       col = c("forestgreen", "blue", "red", adjustcolor("red", alpha.f = 0.15)),
       lty = c(1, 1, 2, NA), pch = c(16, 17, NA, 15), pt.cex = c(1, 1, NA, 2),
       border = NA, bty = "n")

text(x = critical_n_ml, y = 0.05, labels = paste("ML:", critical_n_ml), col = "forestgreen", pos = 4)
text(x = critical_n_tm, y = 0.1, labels = paste("TM:", critical_n_tm), col = "blue", pos = 4)

################!!##################

ml_df <- data.frame(
  sims = n_vals,
  pmo = ml_preds,
  method = "ML (Random Forest)"
)

tm_df <- data.frame(
  sims = tm_indices,
  pmo = tm_pmos,
  method = "Trajectory matching"
)

all_df <- rbind(ml_df, tm_df)

library(scales)

true_prob <- 0.6052 

# Define tolerance band
lower_bound <- true_prob * 0.95
upper_bound <- true_prob * 1.05

# Find first point where ML enters the band
within_band_ml <- which(ml_preds >= lower_bound & ml_preds <= upper_bound)
critical_n_ml <- if (length(within_band_ml)) n_vals[min(within_band_ml)] else NA

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
within_band_tm <- which(tm_pmos >= lower_bound & tm_pmos <= upper_bound)
critical_n_tm <- if (length(within_band_tm)) tm_indices[min(within_band_tm)] else NA


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
  geom_vline(xintercept = critical_n_tm, linetype = "dotted", colour = "blue", size = .8) +
  
  # Main data curves
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  scale_x_log10(labels = comma) +
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

########################################

library(ggplot2)
library(dplyr)
library(scales)
library(grid)

make_pmo_plot <- function(ml_df, tm_df, true_prob, pattern_label = "[I1, I2, I3]") {
  
  # Prepare ±5% tolerance band
  lower_bound <- true_prob * 0.95
  upper_bound <- true_prob * 1.05
  
  # Label methods
  ml_df$method <- "ML (Random Forest)"
  tm_df$method <- "Trajectory matching"
  
  # Combine data
  all_df <- bind_rows(ml_df, tm_df)
  
  # Critical points (first time within 5% band)
  critical_n_ml <- ml_df$sims[which.min(abs(ml_df$pmo - true_prob) <= 0.05 * true_prob)]
  critical_n_tm <- tm_df$sims[which.min(abs(tm_df$pmo - true_prob) <= 0.05 * true_prob)]
  
  # Plot
  p <- ggplot(all_df, aes(x = sims, y = pmo, colour = method, shape = method)) +
    
    # ±5% band
    annotate("rect",
             xmin = min(all_df$sims), xmax = max(all_df$sims),
             ymin = lower_bound, ymax = upper_bound,
             alpha = 0.2, fill = "red") +
    
    # Horizontal true PMO line
    geom_hline(yintercept = true_prob, colour = "red", linetype = "dashed", linewidth = 1) +
    
    # Convergence points (vertical lines)
    geom_vline(xintercept = critical_n_ml, colour = "forestgreen", linetype = "dotted") +
    geom_vline(xintercept = critical_n_tm, colour = "blue", linetype = "dotted") +
    
    # Data points and lines
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    
    # Log scale for x-axis
    scale_x_log10(
      breaks = c(50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000),
      labels = c("50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K", "50K")
    ) +
    
    scale_y_continuous(limits = c(0, 1.05)) +
    
    scale_color_manual(values = c("ML (Random Forest)" = "forestgreen", "Trajectory matching" = "blue")) +
    scale_shape_manual(values = c("ML (Random Forest)" = 16, "Trajectory matching" = 17)) +
    
    labs(
      title = "Convergence of PMO estimates",
      subtitle = paste("Estimated probability of major outbreak for", pattern_label),
      x = "Number of simulations (log scale)",
      y = "Estimated PMO"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey80"),
      panel.grid.major.y = element_line(colour = "grey85"),
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}

ml_df <- data.frame(sims = n_vals, pmo = ml_preds)
tm_df <- data.frame(sims = tm_indices, pmo = tm_pmos)
true_prob <- 0.987
pattern_label <- "[1, 2, 5]"

make_pmo_plot(ml_df, tm_df, true_prob, pattern_label)
