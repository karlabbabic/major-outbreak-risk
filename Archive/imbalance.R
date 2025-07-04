imbalance_ratio <- 0.05  # 5% outbreaks, 95% no outbreaks

set.seed(123)

data_matrix <- read_excel("data_matrix.xlsx", col_names = FALSE)
head(data_matrix)

colnames(data_matrix) <- c("t0", "t1", "t3", "t4", "mechProb", "outbreak_label")
head(data_matrix)

n <- nrow(data_matrix)  # Total number of samples

# Generate target labels
outbreak_label <- sample(c(0, 1), size = n, replace = TRUE, prob = c(1 - imbalance_ratio, imbalance_ratio))
table(outbreak_label)

prop.table(table(data_matrix$outbreak_label))

# Total desired sample size
desired_sample_size <- 251956

# Calculate the number of samples for each class
n_class_0 <- round(0.95 * desired_sample_size)  # 95% of the sample
n_class_1 <- desired_sample_size - n_class_0   # 5% of the sample

print(c(n_class_0, n_class_1))

set.seed(123)  # For reproducibility

# Separate the data by class
class_0 <- data_matrix[data_matrix$outbreak_label == 0, ]
class_1 <- data_matrix[data_matrix$outbreak_label == 1, ]

# Sample from each class
sampled_class_0 <- class_0[sample(1:nrow(class_0), n_class_0, replace = FALSE), ]
sampled_class_1 <- class_1[sample(1:nrow(class_1), n_class_1, replace = FALSE), ]

# Combine the sampled data
sampled_data <- rbind(sampled_class_0, sampled_class_1)

# Shuffle the combined dataset (optional)
sampled_data <- sampled_data[sample(1:nrow(sampled_data)), ]

# Check proportions in the new sample
table(sampled_data$outbreak_label)
prop.table(table(sampled_data$outbreak_label))

nrow(sampled_data)

# Random forest model

n=nrow(sampled_data)
n

trainIndex=sample(1:n, size=0.7*n, replace=FALSE)
train=sampled_data[trainIndex,]
test=sampled_data[-trainIndex,]

#trainIndex <- createDataPartition(data_matrix$outbreak_label, p = 0.7, list = FALSE)
#train <- data_matrix[trainIndex,]
#test <- data_matrix[-trainIndex,]

dim(train)
dim(test)
head(train)

rf_model <- randomForest(as.factor(outbreak_label) ~ ., data = train[,-5], 
                         ntree = 500, mtry = 3, importance = TRUE, probability = TRUE)

# Check feature importance
importance(rf_model)

# Predict probabilities
pred_probs <- predict(rf_model, newdata = test[,1:4], type = "prob")[, 2]

# Binarise predictions using a threshold (e.g., 0.5)
pred_labels <- ifelse(pred_probs > 0.5, 1, 0)

# Evaluate performance
confusion_matrix <- table(Predicted = pred_labels, Actual = test$outbreak_label)
print(confusion_matrix)

# Calculate ROC-AUC
roc_curve <- roc(test$outbreak_label, pred_probs)
auc <- auc(roc_curve)
print(auc)
plot(roc_curve)

test[,5]

mechProb <- test$mechProb

mechPred_labels <- ifelse(mechProb > 0.5, 1, 0)

# Evaluate performance
mech_confusion_matrix <- table(Predicted = mechPred_labels, Actual = test$outbreak_label)
print(mech_confusion_matrix)

mech_roc_curve <- roc(test$outbreak_label, mechProb)
mech_auc <- auc(mech_roc_curve)
print(mech_auc)
plot(mech_roc_curve)

#######

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
