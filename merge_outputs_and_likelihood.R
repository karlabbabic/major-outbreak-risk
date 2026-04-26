# Merging PMO output files (by method and training data)
# and calculating likelihood 
# file output: for each sequence: PMO (calculated using each of the methods)
#                                 likelihood, prior on I_1, likelihood*prior P(I)

library(data.table)

source("compute_analytic_pmo.R")

# Paths
analytic_file   <- "pmo_outputs/analytic_pmo_3_weeks.csv"

ml_file_prior_1 <- "pmo_outputs/prior_1/ml_pmo_3_weeks.csv"
tm_file_prior_1 <- "pmo_outputs/prior_1/tm_pmo_3_weeks.csv"

ml_file_prior_2 <- "pmo_outputs/prior_2/ml_pmo_3_weeks.csv"
tm_file_prior_2 <- "pmo_outputs/prior_2/tm_pmo_3_weeks.csv"

ml_file_prior_3 <- "pmo_outputs/prior_3/ml_pmo_3_weeks.csv"
tm_file_prior_3 <- "pmo_outputs/prior_3/tm_pmo_3_weeks.csv"

out_dir  <- "combined_outputs"
out_file <- file.path(out_dir, "pmo_and_likelihood_3_weeks.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read files

read_ml <- function(file, suffix) {
  dt <- fread(file)
  
  if (!all(c("I1", "I2", "I3", "pmo_ml") %in% names(dt))) {
    stop(file, " is missing one of: I1, I2, I3, pmo_ml")
  }
  
  out <- dt[, .(I1, I2, I3, pmo_ml)]
  setnames(out, "pmo_ml", paste0("ML_PMO_pois", suffix))
  unique(out) # Remove accidental duplicates if any
}

read_tm <- function(file, suffix) {
  dt <- fread(file)
  
  if (!all(c("I1", "I2", "I3", "pmo_tm") %in% names(dt))) {
    stop(file, " is missing one of: I1, I2, I3, pmo_tm")
  }
  
  out <- dt[, .(I1, I2, I3, pmo_tm)]
  setnames(out, "pmo_tm", paste0("TM_PMO_pois", suffix))
  unique(out)
}

analytic_dt <- fread(analytic_file)

# Keep / rename needed columns
if (!all(c("I1", "I2", "I3", "PMO_analytic") %in% names(analytic_dt))) {
  stop("analytic file is missing one of: I1, I2, I3, PMO_analytic")
}
analytic_dt <- unique(analytic_dt[, .(I1, I2, I3, PMO_analytic)])

ml1 <- read_ml(ml_file_prior_1, 1)
tm1 <- read_tm(tm_file_prior_1, 1)

ml2 <- read_ml(ml_file_prior_2, 2)
tm2 <- read_tm(tm_file_prior_2, 2)

ml3 <- read_ml(ml_file_prior_3, 3)
tm3 <- read_tm(tm_file_prior_3, 3)

combined_dt <- Reduce(
  function(x, y) merge(x, y, by = c("I1", "I2", "I3"), all = FALSE),
  list(analytic_dt, ml1, tm1, ml2, tm2, ml3, tm3)
)

setorder(combined_dt, I1, I2, I3)

if (nrow(combined_dt) == 0) {
  stop("No common I1/I2/I3 sequences found across analytic, ML and TM files.")
}

cat("Rows after common-sequence merge:", nrow(combined_dt), "\n")

# Calculate likelihood for each sequence

w <- readRDS("data/serial_interval_weights.rds")

# Attach likelihood column
combined_dt[, likelihood := apply(
  .SD,
  1,
  function(Iseq) {
    integrated_p_I(as.numeric(Iseq), w, R_max = 10)
  }
), .SDcols = c("I1", "I2", "I3")]

# Add I1 priors and posterior columns

prior_pois <- function(i1, lambda) {
  vals <- dpois(i1, lambda = lambda)
  vals / (1 - dpois(0, lambda = lambda))
}

add_I1_prior_column <- function(dt, lambda) {
  prior_col <- paste0("Prior_pois_", lambda)
  post_col  <- paste0("Posterior_pois_", lambda)
  
  pi1 <- prior_pois(dt$I1, lambda = lambda)
  unnorm <- dt$likelihood * pi1
  s <- sum(unnorm, na.rm = TRUE)
  
  if (!is.finite(s) || s <= 0) {
    stop("All posterior weights are zero for: ", post_col)
  }
  
  dt[[prior_col]] <- pi1
  dt[[post_col]] <- unnorm / s
  
  dt
}

combined_dt <- add_I1_prior_column(combined_dt, 1)
combined_dt <- add_I1_prior_column(combined_dt, 2)
combined_dt <- add_I1_prior_column(combined_dt, 3)

# should sum to 1
cat("Posterior checks:\n")
cat("  Posterior_pois_1:", sum(combined_dt$Posterior_pois_1), "\n")
cat("  Posterior_pois_2:", sum(combined_dt$Posterior_pois_2), "\n")
cat("  Posterior_pois_3:", sum(combined_dt$Posterior_pois_3), "\n")

# Save
fwrite(combined_dt, out_file)

cat("Saved combined file to:", out_file, "\n")
cat("Rows:", nrow(combined_dt), "\n")
cat("Columns:", paste(names(combined_dt), collapse = ", "), "\n")
print(head(combined_dt))