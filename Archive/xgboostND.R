library(xgboost)
library(pROC)

set.seed(123)

data_matrix <- read_excel("data_matrix.xlsx", col_names = FALSE)
head(data_matrix)

colnames(data_matrix) <- c("t0", "t1", "t3", "t4", "mechProb", "outbreak_label")
head(data_matrix)

n=nrow(data_matrix)
n

trainIndex=sample(1:n, size=.7*n, replace=FALSE)
train=data_matrix[trainIndex,]
test=data_matrix[-trainIndex,]

dim(train)
head(train)

X <- as.matrix(train[, 1:4]) # Predictor variables

y <- as.numeric(as.character(train$outbreak_label)) # Target variable (last col)

# Create a DMatrix

dtrain <- xgb.DMatrix(data = X, label = y)

# Set parameters for XGBoost

params <- list(
  
  booster = "gbtree", # Use tree-based models
  
  objective = "binary:logistic", # For binary classification
  
  eval_metric = "auc", # Evaluation metric
  
  eta = 0.1, # Learning rate
  
  max_depth = 6, # Maximum tree depth
  
  subsample = 0.8, # Fraction of data to use for training each tree
  
  colsample_bytree = 0.8 # Fraction of features to use for training each tree
  
)

# Train the model

xgb_model <- xgb.train(
  
  params = params,
  
  data = dtrain,
  
  nrounds = 500, # Number of boosting iterations (trees)
  
  watchlist = list(train = dtrain), # Monitor training
  
  print_every_n = 50 # Print progress every 50 rounds
  
)

# Prepare the test dataset
X_test <- as.matrix(test[, 1:4])
y_test <- as.numeric(as.character(test$outbreak_label))
dtest <- xgb.DMatrix(data = X_test, label = y_test)

# Predict probabilities
y_pred <- predict(xgb_model, newdata = dtest)

# Evaluate performance
xgb_roc_curve <- roc(y_test, y_pred)
auc(xgb_roc_curve)  # Calculate AUC
plot(xgb_roc_curve)


mechProb <- test$mechProb

mech_roc_curve <- roc(test$outbreak_label, mechProb)
mech_auc <- auc(mech_roc_curve)
print(mech_auc)
plot(mech_roc_curve)

mechProb <- data_matrix$mechProb
mech_roc_curve <- roc(data_matrix$outbreak_label, mechProb)
mech_auc <- auc(mech_roc_curve)
print(mech_auc)
plot(mech_roc_curve)


cor(y_pred_probs, test$mechProb)  # Correlation between ML and mechanistic probabilities



# Binarise predictions using a threshold (e.g., 0.5)
y_pred_labels <- ifelse(y_pred > 0.5, 1, 0)

# Evaluate performance
confusion_matrix <- table(Predicted = y_pred_labels, Actual = test$outbreak_label)
print(confusion_matrix)


mechProb <- test$mechProb

mechPred_labels <- ifelse(mechProb > 0.5, 1, 0)

# Evaluate performance
mech_confusion_matrix <- table(Predicted = mechPred_labels, Actual = test$outbreak_label)
print(mech_confusion_matrix)

mechProb <- data_matrix$mechProb

mechPred_labels <- ifelse(mechProb > 0.5, 1, 0)

# Evaluate performance
mech_confusion_matrix <- table(Predicted = mechPred_labels, Actual = data_train$outbreak_label)
print(mech_confusion_matrix)

cv <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,                  # 5-fold cross-validation
  early_stopping_rounds = 10, # Stop if performance doesn't improve
  maximize = TRUE,
  verbose = TRUE
)

# Use the optimal number of rounds from cross-validation
best_nrounds <- cv$best_iteration

# Retrain the model with the best number of boosting rounds
final_model <- xgb.train(
  params = params,             # Use the same parameters
  data = dtrain,               # Training data
  nrounds = best_nrounds,      # Optimal number of rounds from CV
  watchlist = list(train = dtrain),  # Optional: monitor training performance
  print_every_n = 50           # Optional: print progress every 50 iterations
)

y_pred <- predict(final_model, newdata = dtest)

calibration_xgb <- calibration(test$outbreak_label, y_pred, classProbs = TRUE)
calibration_mech <- calibration(test$outbreak_label, mechProb, classProbs = TRUE)

plot(calibration_xgb, main = "XGBoost Calibration")
plot(calibration_mech, main = "mechProb Calibration")

# Evaluate performance
final_roc_curve <- roc(y_test, y_pred)
auc(final_roc_curve)  # Calculate AUC
auc(xgb_roc_curve)

plot(final_roc_curve, col = "blue", main = "ROC Curve Comparison", lwd = 0.3)
plot(mech_roc_curve, col = "red", add = TRUE, lwd = 0.3)
legend("bottomright", legend = c("Classifier (XGBoost)", "Mechanistic Model"), col = c("blue", "red"), lwd = 2)
