# i) Function to train ranger (random forest) classifer on increasing sizes of data
# ii) Predict pmo function for ML model

library(ranger)
library(ggplot2)
library(dplyr)

####
# Define training sizes
n_vals <- c(seq(5, 100, by = 5), seq(200, 1000, by = 100),
            seq(2000, 10000, by = 1000), seq(20000, 100000, by = 10000),
            seq(200000, 1000000, by = 100000))

# Read and prepare data (should be a df)
sim_df_prior <- read.csv("data_matrix_incidence_rm.csv")
sim_df_prior$outbreak_label <- as.factor(sim_df_prior$outbreak_label)

train_ranger_models_over_sizes <- function(data, label_col, n_vals,
                                           seed = 123, ntree = 500, num.threads = 1) {
  set.seed(seed)
  
  # Shuffle indices once
  shuffled_indices <- sample(nrow(data))
  
  # Train models on cumulative subsets (prepare output list)
  rf_models <- vector("list", length(n_vals))
  names(rf_models) <- paste0("n_", n_vals)
  
  for (i in seq_along(n_vals)) {
    n <- n_vals[i]
    sub_indices <- shuffled_indices[1:n]  # Cumulative: always the first n rows
    sub_train <- data[sub_indices, ]
    
    predictors <- setdiff(names(sub_train), label_col)
    formula <- as.formula(paste(label_col, "~", paste(predictors, collapse = " + ")))
    
    rf_models[[i]] <- ranger(
      formula = formula,
      data = sub_train,
      probability = TRUE,
      mtry = length(predictors),
      num.trees = ntree,
      num.threads = num.threads
    )
    
    cat(sprintf("Trained model for n = %d\n", n))
  }
  
  return(rf_models)
}

rf_models <- train_ranger_models_over_sizes(
  data = sim_df_prior,
  label_col = "outbreak_label",
  n_vals = n_vals,
  ntree = 1000,          # Increase number of trees
  num.threads = 4        # Use 4 CPU cores
)

# Save!
# saveRDS(rf_models, file = "[...].rds")

# model <- readRDS("rf_model_3days.rds")

######

# Predict pmo function for ML model
predict_pmo <- function(input_pattern, rf_models, n_vals) {
  preds <- numeric(length(rf_models))
  for (i in seq_along(rf_models)) {
    preds[i] <- predict(rf_models[[i]], data = input_pattern)$predictions[, "1"]
  }
  names(preds) <- n_vals
  return(preds)
}

## Test:
# Define a new input pattern
target_input <- data.frame(t1 = 1, t2 = 2, t3 = 0)

# rf_models <- readRDS("[...].rds")
# Get predictions for all training sizes
ml_preds <- predict_pmo(target_input, rf_models, n_vals)

print(ml_preds)