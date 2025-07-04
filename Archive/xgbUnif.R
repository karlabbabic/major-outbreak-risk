
xgb_models <- vector("list", length(n_vals))
names(xgb_models) <- paste0("n_", n_vals)

# Shuffle once
set.seed(123)
shuffled_indices <- sample(nrow(sim_df_prior))

# Train models
for (i in seq_along(n_vals)) {
  n <- n_vals[i]
  sub_indices <- shuffled_indices[1:n]
  sub_train <- sim_df_prior[sub_indices, ]
  
  # Prepare training data
  x <- as.matrix(sub_train[, c("t1", "t2", "t3")])
  y <- as.numeric(as.character(sub_train$outbreak_label))  # Ensure 0/1
  
  dtrain <- xgb.DMatrix(data = x, label = y)
  
  xgb_models[[i]] <- xgboost(
    data = dtrain,
    objective = "binary:logistic",
    eval_metric = "logloss",
    nrounds = 50,
    verbose = 0
  )
  
  cat(sprintf("Trained XGBoost model for n = %d\n", n))
}

# Define your input pattern
target <- c(t1 = 1, t2 = 0, t3 = 0)
target_input <- as.matrix(data.frame(t1 = target[1], t2 = target[2], t3 = target[3]))

# Predict using all stored models
xgb_preds <- sapply(xgb_models, function(model) predict(model, target_input))

# Wrap into a dataframe for plotting

xgb_df <- data.frame(sims = n_vals, pmo = xgb_preds, method = "ML (XGBoost)")
