# Ranger (rf but faster?)

library(ranger)

set.seed(123)

# Train the ranger model
ranger_model <- ranger(
  formula = as.factor(outbreak_label) ~ .,  # Target variable vs all predictors
  data = train[,-5],
  num.trees = 500,         # Number of trees
  mtry = 3,                # Number of features to use per split
  probability = TRUE,       # Enables probability predictions
  importance = "impurity",  # Enables feature importance calculation
  classification = TRUE,    # Binary classification
  seed = 123                # Ensures reproducibility
)

y_pred_probs <- predict(ranger_model, data = test[,1:4])$predictions[, 2]

roc_curve <- roc(test$outbreak_label, y_pred_probs)
print(auc(roc_curve))  # Print AUC value
plot(roc_curve)

cor(y_pred_probs, test$mechProb)

residuals_mech <- abs(test$outbreak_label - test$mechProb)
residuals_ml <- abs(test$outbreak_label - y_pred_probs)

cor(residuals_ml, residuals_mech)


importance <- ranger_model$variable.importance
print(importance)

barplot(sort(importance, decreasing = TRUE), main = "Feature Importance", col = "blue")

cor_matrix <- cor(data_matrix[, c("t0", "t1", "t3", "t4")])
print(cor_matrix)

ranger_model <- ranger(target ~ ., data = dataset, importance = "permutation")
ranger_model <- ranger(
  formula = as.factor(outbreak_label) ~ .,  # Target variable vs all predictors
  data = train[,-5],
  num.trees = 500,         # Number of trees
  mtry = 3,                # Number of features to use per split
  probability = TRUE,       # Enables probability predictions
  importance = "permutation",  # Enables feature importance calculation
  classification = TRUE,    # Binary classification
  seed = 123                # Ensures reproducibility
)

print(ranger_model$variable.importance)

y_pred_probs <- predict(ranger_model, data = test[,1:4])$predictions[, 2]

roc_curve <- roc(test$outbreak_label, y_pred_probs)
print(auc(roc_curve))  # Print AUC value
plot(roc_curve)

