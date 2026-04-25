# Convergence analysis for ML PMO estimates.
# Train 2-week and 3-week ML models on increasing training-set sizes,
# predict PMO over grid of incidence sequences and save merged outputs
# for downstream convergence plotting.

library(data.table)
library(ranger)

# Settings
rf_num_trees <- 1000
rf_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

master_file <- "training_data/training_data_pointmass_I5_1mil.rds"
base_dir <- "training_convergence"

if (!file.exists(master_file)) {
  stop("Master dataset not found: ", master_file)
}

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading master dataset...\n")
master_dt <- readRDS(master_file)
cat("Master dataset loaded. Rows:", nrow(master_dt), "\n")

# find earliest n where a positive outcome exists
pos_idx <- which(master_dt$outcome_bin == "1")[1]

if (is.na(pos_idx)) {
  stop("No positive (outcome_bin == 1) found in the master dataset. Cannot train RF probability model.")
}

cat("First positive outcome occurs at row n =", pos_idx, "\n")

# Define training sizes

sizes_small  <- 1:20000
sizes_mid    <- seq(30000, 100000, by = 10000)
sizes_large <- seq(200000, 1000000, by = 100000)

sizes <- c(sizes_small, sizes_mid, sizes_large)
# Start sizes only once we have seen at least one positive
sizes <- sizes[sizes >= pos_idx]

if (length(sizes) == 0) {
  stop("No training sizes remain after filtering.")
}

cat("Training will start at n =", min(sizes), "and run", length(sizes), "sizes total.\n")

# Helper functions: generate grid of inputs -> compute PMOs

generate_candidate_sequences <- function(Y, buffer, I1 = 1, fix_I1 = FALSE) {
  T <- length(Y)
  I_ranges <- vector("list", T)
  for (t in 1:T) I_ranges[[t]] <- Y[t]:(Y[t] + buffer)
  if (fix_I1) I_ranges[[1]] <- I1
  grid <- expand.grid(I_ranges)
  colnames(grid) <- paste0("I", 1:T)
  grid
}

# PMO predict function
predict_pmo_grid <- function(grid_df, rf_model, positive_class = "1") {
  preds <- predict(rf_model, data = grid_df)$predictions
  as.numeric(preds[, positive_class])
}

# Fixed 3-week grid

I1 <- 1       # Optionally fix initial case
buffer <- 5   # Cap I_t
Y_obs <- c(I1, 0, 0)

grid_base <- generate_candidate_sequences(
  Y = Y_obs,
  buffer = buffer,
  I1 = I1,
  fix_I1 = TRUE
)
grid_base <- as.data.table(grid_base)

# MAIN LOOP OVER TRAINING SIZES

for (n in sizes) {
  cat("\n==============================\n")
  cat("Training for n =", n, "\n")
  cat("==============================\n")

  dt <- master_dt[1:n]

  # Model A: I1, I2 (model trained on 2 weeks of data)
  
  model_A <- ranger(
    outcome_bin ~ I1 + I2,
    data = dt[, .(I1, I2, outcome_bin)],
    probability = TRUE,
    mtry = 2,
    num.trees = rf_num_trees,
    num.threads = rf_threads
  )

  # Model B: I1, I2, I3 (model trained on 3 weeks of data)
  
  model_B <- ranger(
    outcome_bin ~ I1 + I2 + I3,
    data = dt[, .(I1, I2, I3, outcome_bin)],
    probability = TRUE,
    mtry = 3,
    num.trees = rf_num_trees,
    num.threads = rf_threads
  )

  # PMO grid
  gridA <- unique(grid_base[, .(I1, I2)])
  gridA$pmo_ml <- predict_pmo_grid(gridA[, .(I1, I2)], model_A)

  gridB <- copy(grid_base)
  gridB$pmo_ml <- predict_pmo_grid(gridB[, .(I1, I2, I3)], model_B)

  # Save outputs
  out_dir <- file.path("training_convergence", paste0("n_", n))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save PMO estimates with corresponding incidence data
  fwrite(gridA, file.path(out_dir, "ml_pmo_2_weeks.csv"))
  fwrite(gridB, file.path(out_dir, "ml_pmo_3_weeks.csv"))
  
  # Save classfiers/models
  saveRDS(model_A, file.path(out_dir, "rf_2_weeks.rds"))
  saveRDS(model_B, file.path(out_dir, "rf_3_weeks.rds"))

  cat("Finished n =", n, "\n")
}

cat("\nAll training sizes completed successfully.\n")

# Merge PMO outputs across n

dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("^n_", basename(dirs))]
dirs <- dirs[order(as.numeric(sub("n_", "", basename(dirs))))]

get_n <- function(path) {
  as.integer(sub("n_", "", basename(path)))
}

# Helper: merge a given filename across all n
merge_files <- function(filename) {
  
  out_list <- list()
  
  for (d in dirs) {
    
    n_val <- get_n(d)
    f <- file.path(d, filename)
    
    if (!file.exists(f)) next
    
    dt <- fread(f)
    dt[, n_train := n_val]
    
    out_list[[length(out_list) + 1]] <- dt
  }
  
  if (length(out_list) == 0) {
    stop("No files found for ", filename)
  }
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

cat("Merging PMO 2-week outputs...\n")
pmo_ml_2_weeks <- merge_files("ml_pmo_2_weeks.csv")
setorder(pmo_ml_2_weeks, I1, I2, n_train)
fwrite(pmo_ml_2_weeks, file.path(base_dir, "merged_pmo_2_weeks.csv"))

cat("Merging PMO 3-week outputs...\n")
pmo_ml_3_weeks <- merge_files("ml_pmo_3_weeks.csv")
setorder(pmo_ml_3_weeks, I1, I2, I3, n_train)
fwrite(pmo_ml_3_weeks, file.path(base_dir, "merged_pmo_3_weeks.csv"))

cat("All merging complete.\n")