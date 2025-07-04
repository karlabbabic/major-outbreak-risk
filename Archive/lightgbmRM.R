# --------------------------
# Step 1: Set up parameters
# --------------------------
set.seed(42)

R <- 1.5                  # Reproduction number
T <- 200                 # Total time steps
w <- c(0.2, 0.3, 0.3, 0.1, 0.1)  # Serial interval distribution
S <- length(w)           # Serial interval length
n_sim <- 10000          # Number of simulations
threshold <- 200         # Outbreak threshold

# --------------------------
# Step 2: Run simulations
# --------------------------
# Store: day1, day2, day3, outbreak_label
summary_mat <- matrix(0, nrow = n_sim, ncol = 4)
colnames(summary_mat) <- c("day1", "day2", "day3", "outbreak")

for (sim in 1:n_sim) {
  I <- numeric(T)
  I[1] <- 1  # Initial condition: always start with 1 case
  
  for (t in 2:T) {
    lambda <- sum(w[1:min(t - 1, S)] * I[(t - 1):(t - min(t - 1, S))])
    I[t] <- rpois(1, R * lambda)
  }
  
  summary_mat[sim, 1:3] <- I[1:3]
  summary_mat[sim, 4] <- as.integer(sum(I) > threshold)
  
  if (sim %% 10000 == 0) cat("Completed", sim, "simulations\n")
}

# Convert to data frame
sim_df <- as.data.frame(summary_mat)
colnames(sim_df) <- c("t1", "t2", "t3", "outbreak_label")

# --------------------------
# Step 3: Train LightGBM
# --------------------------
library(lightgbm)

# Split into train/test
set.seed(123)
n <- nrow(sim_df)
trainIndex <- sample(1:n, size = 0.7 * n, replace = FALSE)
train <- sim_df[trainIndex, ]
test  <- sim_df[-trainIndex, ]

# Prepare matrices
X_train <- as.matrix(train[, 1:3])
y_train <- train$outbreak_label
X_test  <- as.matrix(test[, 1:3])
y_test  <- test$outbreak_label

dtrain <- lgb.Dataset(data = X_train, label = y_train)
dtest  <- lgb.Dataset(data = X_test, label = y_test)

# Set parameters
params <- list(
  objective = "binary",
  metric = "auc",
  boosting = "gbdt",
  num_leaves = 31,
  learning_rate = 0.1,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 1
)

# Train model
model <- lgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  valids = list(test = dtest),
  early_stopping_rounds = 10,
  verbose = 0
)

# --------------------------
# Step 4: Predict for [1, 2, 0]
# --------------------------
X_query <- matrix(c(1, 2, 0), nrow = 1)
ml_pred <- predict(model, X_query)
cat("Predicted probability of major outbreak for [1, 2, 0]:", round(ml_pred, 4), "\n")
