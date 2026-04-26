# 1. Trains classifier on a given simulated dataset
# Rerun this script separately for each dataset, updating input_file and out_dir.
# 2. Generates PMO estimates using Methods 2 and 3 (trajectory matching & ML)

library(data.table)
library(ranger)

# load training data
master_file <- "data/training_data_pois_1_30mil.rds"
if (!file.exists(master_file)) stop("Master dataset not found: ", master_file)
master_dt <- readRDS(master_file)

rf_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
# Select number of weeks of interest (e.g. train on 3 weeks of data) and outcome (minor/major)
dt <- as.data.table(master_dt)[, .(I1, I2, I3, outcome_bin)]

dt[, outcome_bin := factor(as.character(outcome_bin), levels = c("0", "1"))]

n <- nrow(dt)

cat("Training on full dataset\n")
cat("Rows =", format(n, big.mark = ","), "\n")

# Train classifier
model <- ranger(
  outcome_bin ~ I1 + I2 + I3,
  data = dt,
  probability = TRUE,
  mtry = 3,
  num.trees = 1000,
  num.threads = rf_threads
)

# TM pmos = proportion of each {I1,I2,I3} trajectory in training set 
TM <- dt[
  ,
  .(pmo_tm = mean(as.integer(as.character(outcome_bin)))),
  by = .(I1, I2, I3)
][]

# for ML PMO estimation:
# precompute PMOs for all possible sequences (for efficiency for underreporting analysis)

# Candidate sequence generator
generate_candidate_sequences <- function(Y, buffer, I1 = NULL, fix_I1 = FALSE) {
  T <- length(Y)
  I_ranges <- vector("list", T)
  
  for (t in seq_len(T)) {
    I_ranges[[t]] <- Y[t]:(Y[t] + buffer)
  }
  
  # e.g. when prior on I_1 is a point mass with P(I_1=1)=1
  if (fix_I1) {
    if (is.null(I1)) stop("If fix_I1 = TRUE, provide I1.")
    I_ranges[[1]] <- I1
  }
  
  grid <- expand.grid(I_ranges)
  colnames(grid) <- paste0("I", seq_len(T))
  grid
}

# ML PMO function
predict_pmo_grid <- function(grid_df, rf_model, positive_class = "1") {
  preds <- predict(rf_model, data = grid_df)$predictions
  as.numeric(preds[, positive_class])
}

# limit generated sequences
buffer <- 20
Y_obs <- c(1, 0, 0)

grid <- as.data.table(
  generate_candidate_sequences(
    Y = Y_obs,
    buffer = buffer,
    I1 = NULL, 
    fix_I1 = FALSE
    )
)

# apply predict fn on sequences
grid[, pmo_ml := predict_pmo_grid(grid[, .(I1, I2, I3)], model)]

# Save outputs
out_dir <- "pmo_outputs/prior_1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# save TM PMOs
fwrite(TM, file.path(out_dir, "tm_pmo_3_weeks.csv"))
# save ML PMOs
fwrite(grid, file.path(out_dir, "ml_pmo_3_weeks.csv"))
# save classifier
saveRDS(model, file.path(out_dir, "rf_3_weeks.rds"))

cat("Finished training on full dataset\n")
