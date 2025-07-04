library(lightgbm)
library(pROC)

#RN data (SIR)
data_matrix_incidence <- read.csv("data_matrix_incidence.csv")
colnames(data_matrix_incidence) <- c("t0", "t1", "t2", "outbreak_label")
head(data_matrix_incidence)

#KB data (SIR)
data_matrix_incidence <- read.csv("data_matrix_incidence_KB.csv")
colnames(data_matrix_incidence) <- c("t0", "t1", "t2", "outbreak_label")
head(data_matrix_incidence)

#KB data (Renewal model)
data_matrix_incidence <- read.csv("data_matrix_incidence_KB_m2.csv")
colnames(data_matrix_incidence) <- c("t0", "t1", "t2", "outbreak_label")
head(data_matrix_incidence)
data_matrix_incidence[] <- lapply(data_matrix_incidence, as.numeric)


# Split into training and test sets
set.seed(123)
n <- nrow(data_matrix_incidence)
trainIndex <- sample(1:n, size = 0.7 * n, replace = FALSE)
train <- data_matrix_incidence[trainIndex, ]
test <- data_matrix_incidence[-trainIndex, ]

dim(test)
dim(train)

# Prepare matrices for LightGBM
X_train <- as.matrix(train[, -4])  # Predictors
y_train <- as.numeric(train$outbreak_label) # Target

X_test <- as.matrix(test[, -4])
y_test <- as.numeric(test$outbreak_label)
dtrain <- lgb.Dataset(data = X_train, label = y_train)
dtest <- lgb.Dataset(data = X_test, label = y_test)
params <- list(
  objective = "binary",       # Binary classification
  metric = "auc",             # Evaluation metric (AUC)
  boosting = "gbdt",          # Gradient boosting decision trees
  num_leaves = 31,            # Maximum leaves per tree
  learning_rate = 0.1,        
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

###############

## resim. with same params as trajec.
###

X_query <- matrix(data = c(1, 2, 0), nrow = 1, ncol = 3)

# Predict prob
ml_pred <- predict(model, X_query)
ml_pred

subset(sim_df, day3 > 4 & outbreak == 1)
subset(data_matrix_incidence, t2 > 4 & outbreak_label == 1)
subset(data_matrix_incidence, t0==1& t1==1& t2==6)
subset(sim_df, day1==1& day2==1& day3==6)

log1p(sim_df$day1)  # log(1 + x)
log1p(sim_df$day2)
log1p(sim_df$day3)


# t0 need to be 1?

#importance?
