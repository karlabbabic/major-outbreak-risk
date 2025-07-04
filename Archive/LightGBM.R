library(lightgbm)
library(readxl)
library(pROC)

data_matrix <- read_excel("data_matrix.xlsx", col_names = FALSE)
head(data_matrix)

colnames(data_matrix) <- c("t0", "t1", "t3", "t4", "mechProb", "outbreak_label")
head(data_matrix)


# Split into training and test sets
set.seed(123)
n <- nrow(data_matrix)
trainIndex <- sample(1:n, size = 0.7 * n, replace = FALSE)
train <- data_matrix[trainIndex, ]
test <- data_matrix[-trainIndex, ]

# Prepare matrices for LightGBM
X_train <- as.matrix(train[, 1:4])  # Predictors
y_train <- as.numeric(train$outbreak_label) # Target

X_test <- as.matrix(test[, 1:4])
y_test <- as.numeric(test$outbreak_label)
dtrain <- lgb.Dataset(data = X_train, label = y_train)
dtest <- lgb.Dataset(data = X_test, label = y_test)
params <- list(
  objective = "binary",       # Binary classification
  metric = "auc",             # Evaluation metric (AUC)
  boosting = "gbdt",          # Gradient boosting decision trees
  num_leaves = 31,            # Maximum leaves per tree
  learning_rate = 0.1,        # Learning rate
  feature_fraction = 0.8,     # Fraction of features to use per iteration
  bagging_fraction = 0.8,     # Fraction of data to use per iteration
  bagging_freq = 1            # Perform bagging every iteration
)
# Train the LightGBM model
model <- lgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,              # Number of boosting rounds
  valids = list(test = dtest), # Validation data for monitoring
  early_stopping_rounds = 10, # Stop if validation doesn't improve for 10 rounds
  verbose = 1                 # Show progress
)
# Predict probabilities
lgb_pred_probs <- predict(model, X_test)
roc_curve <- roc(y_test, lgb_pred_probs)
auc_value <- auc(roc_curve)
print(auc_value)  # AUC value
plot(roc_curve)   # ROC curve

roc_curve <- roc(test$mechProb, lgb_pred_probs)
auc_value <- auc(roc_curve)
print(auc_value)  # AUC value
plot(roc_curve)   # ROC curve

cor(lgb_pred_probs, test$mechProb) 

residuals_mech <- abs(test$outbreak_label - test$mechProb)
residuals_ml <- abs(test$outbreak_label - lgb_pred_probs)

cor(residuals_ml, residuals_mech)
