for (m in c(1, 2, 3)) {
  model <- ranger(outbreak_label ~ ., data = train, probability = TRUE, mtry = m)
  preds <- predict(model, data = test)$predictions[, "1"]
  auc <- pROC::roc(test$outbreak_label, preds)$auc
  cat("mtry =", m, "→ AUC:", round(auc, 4), "\n")
}

rf_m2 <- ranger(outbreak_label ~ ., data = train, mtry = 2, probability = TRUE)
rf_m3 <- ranger(outbreak_label ~ ., data = train, mtry = 3, probability = TRUE)

X_test <- test[, 1:3]
y_test <- as.numeric(as.character(test$outbreak_label))

# Get probabilities from both models
pred_m2 <- predict(rf_m2, data = X_test)$predictions[, "1"]
pred_m3 <- predict(rf_m3, data = X_test)$predictions[, "1"]

brier_m2 <- mean((pred_m2 - y_test)^2)
brier_m3 <- mean((pred_m3 - y_test)^2)

cat("Brier (mtry=2):", round(brier_m2, 4), "\n")
cat("Brier (mtry=3):", round(brier_m3, 4), "\n")

library(ggplot2)
library(dplyr)

calib_plot_data <- function(probs, truth, bins = 10) {
  df <- data.frame(prob = probs, truth = truth)
  df$bin <- cut(df$prob, breaks = bins)
  
  df %>%
    group_by(bin) %>%
    summarise(
      mean_prob = mean(prob),
      observed_freq = mean(truth),
      n = n()
    ) %>%
    filter(!is.na(bin))
}

# Create data for both models
calib2 <- calib_plot_data(pred_m2, y_test)
calib3 <- calib_plot_data(pred_m3, y_test)

# Plot
ggplot() +
  geom_line(data = calib2, aes(x = mean_prob, y = observed_freq), colour = "blue", size = 1.2) +
  geom_line(data = calib3, aes(x = mean_prob, y = observed_freq), colour = "red", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Predicted probability", y = "Observed frequency",
       title = "Calibration plot: RF mtry = 2 (blue) vs mtry = 3 (red)") +
  theme_minimal()