# Calculate PMO under an underreporting model for a chosen observed sequence Y.
# Starting from posterior weights over latent incidence sequences, update to
# P(I | Y, p) across a grid of reporting probabilities, then save both
# posterior-weighted PMO summaries and row-level posterior outputs for plotting.

library(data.table)

# Observed sequence
Y <- c(1, 1, 2)
# Max Yt+It
buffer <- 18

in_file <- "combined_outputs/pmo_and_likelihood_3_weeks.csv"
out_dir <- "combined_outputs/underreporting"
out_file <- file.path(
  out_dir,
  sprintf("summary_over_p_Y_%s.csv", paste(Y, collapse = "_"))
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Binomial reporting model
p_Y_given_I <- function(Y, I, p) {
  if (any(I < Y)) return(0)
  prod(dbinom(Y, size = I, prob = p))
}

# Renormalise posterior after multiplying by binomial term
update_posterior_with_binomial <- function(base_post, binom) {
  unnorm <- base_post * binom
  s <- sum(unnorm, na.rm = TRUE)
  
  if (!is.finite(s) || s <= 0) {
    return(rep(NA_real_, length(base_post)))
  }
  
  unnorm / s
}

# Proper weighted mean on available rows
weighted_mean_available <- function(x, w) {
  keep <- !is.na(x) & !is.na(w)
  s <- sum(w[keep], na.rm = TRUE)
  
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  
  sum(x[keep] * w[keep], na.rm = TRUE) / s
}

# Read data
dt <- fread(in_file)

required_cols <- c(
  "I1", "I2", "I3",
  "Posterior_pois_1", "Posterior_pois_2", "Posterior_pois_3",
  "PMO_analytic",
  "ML_PMO_pois1", "ML_PMO_pois2", "ML_PMO_pois3",
  "TM_PMO_pois1", "TM_PMO_pois2", "TM_PMO_pois3"
)

missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Restrict to common sequence support: I_t <= Y_t + buffer
# and I_t >= Y_t
I_cols <- c("I1", "I2", "I3")
if (length(Y) != length(I_cols)) {
  stop("Length of Y must match number of I columns.")
}

upper <- Y + buffer

dt <- dt[
  I1 >= Y[1] & I1 <= upper[1] &
    I2 >= Y[2] & I2 <= upper[2] &
    I3 >= Y[3] & I3 <= upper[3]
]

if (nrow(dt) == 0) {
  stop("No rows remain after applying common buffer filter.")
}

cat("Using common capped support with buffer =", buffer, "\n")
cat("Y =", paste(Y, collapse = ","), "\n")
cat("Max allowed I =", paste(upper, collapse = ","), "\n")
cat("Rows retained =", nrow(dt), "\n")

# Loop over reporting prob, p and collect summaries
p_grid <- seq(0.3, 1, by = 0.01)
res_list <- vector("list", length(p_grid))
posterior_list <- vector("list", length(p_grid))

for (idx in seq_along(p_grid)) {
  p <- p_grid[idx]
  message("Processing p = ", sprintf("%.2f", p))
  
  binom_only <- apply(dt[, .(I1, I2, I3)], 1, function(Iseq) {
    p_Y_given_I(Y = Y, I = as.numeric(Iseq), p = p)
  })
  
  post1 <- update_posterior_with_binomial(dt$Posterior_pois_1, binom_only)
  post2 <- update_posterior_with_binomial(dt$Posterior_pois_2, binom_only)
  post3 <- update_posterior_with_binomial(dt$Posterior_pois_3, binom_only)
  
  posterior_list[[idx]] <- data.table(
    Y1 = Y[1],
    Y2 = Y[2],
    Y3 = Y[3],
    buffer = buffer,
    p = p,
    
    I1 = dt$I1,
    I2 = dt$I2,
    I3 = dt$I3,
    I_sequence = paste(dt$I1, dt$I2, dt$I3, sep = "-"),
    
    binom_likelihood = binom_only,
    
    P_I_prior_pois1 = dt$Posterior_pois_1,
    P_I_prior_pois2 = dt$Posterior_pois_2,
    P_I_prior_pois3 = dt$Posterior_pois_3,
    
    P_I_given_Y_pois1 = post1,
    P_I_given_Y_pois2 = post2,
    P_I_given_Y_pois3 = post3,
    
    PMO_analytic = dt$PMO_analytic,
    
    ML_PMO_pois1 = dt$ML_PMO_pois1,
    ML_PMO_pois2 = dt$ML_PMO_pois2,
    ML_PMO_pois3 = dt$ML_PMO_pois3,
    
    TM_PMO_pois1 = dt$TM_PMO_pois1,
    TM_PMO_pois2 = dt$TM_PMO_pois2,
    TM_PMO_pois3 = dt$TM_PMO_pois3
  )
  
  res_list[[idx]] <- data.table(
    Y1 = Y[1],
    Y2 = Y[2],
    Y3 = Y[3],
    buffer = buffer,
    max_I1 = upper[1],
    max_I2 = upper[2],
    max_I3 = upper[3],
    n_sequences = nrow(dt),
    p = p,
    
    posterior_sum_pois1 = sum(post1, na.rm = TRUE),
    posterior_sum_pois2 = sum(post2, na.rm = TRUE),
    posterior_sum_pois3 = sum(post3, na.rm = TRUE),
    
    analytic_pmo_pois1 = weighted_mean_available(dt$PMO_analytic, post1),
    analytic_pmo_pois2 = weighted_mean_available(dt$PMO_analytic, post2),
    analytic_pmo_pois3 = weighted_mean_available(dt$PMO_analytic, post3),
    
    ml_pmo_pois1 = weighted_mean_available(dt$ML_PMO_pois1, post1),
    ml_pmo_pois2 = weighted_mean_available(dt$ML_PMO_pois2, post2),
    ml_pmo_pois3 = weighted_mean_available(dt$ML_PMO_pois3, post3),
    
    tm_pmo_pois1 = weighted_mean_available(dt$TM_PMO_pois1, post1),
    tm_pmo_pois2 = weighted_mean_available(dt$TM_PMO_pois2, post2),
    tm_pmo_pois3 = weighted_mean_available(dt$TM_PMO_pois3, post3),
    
    ml_support_mass_pois1 = sum(post1[!is.na(dt$ML_PMO_pois1)], na.rm = TRUE),
    ml_support_mass_pois2 = sum(post2[!is.na(dt$ML_PMO_pois2)], na.rm = TRUE),
    ml_support_mass_pois3 = sum(post3[!is.na(dt$ML_PMO_pois3)], na.rm = TRUE),
    
    tm_support_mass_pois1 = sum(post1[!is.na(dt$TM_PMO_pois1)], na.rm = TRUE),
    tm_support_mass_pois2 = sum(post2[!is.na(dt$TM_PMO_pois2)], na.rm = TRUE),
    tm_support_mass_pois3 = sum(post3[!is.na(dt$TM_PMO_pois3)], na.rm = TRUE)
  )
}

summary_dt <- rbindlist(res_list)
posterior_dt <- rbindlist(posterior_list)
posterior_out_file <- file.path(
  out_dir,
  sprintf("posterior_rows_over_p_Y_%s.csv", paste(Y, collapse = "_"))
)

fwrite(posterior_dt, posterior_out_file)
fwrite(summary_dt, out_file)

cat("Saved summary to:", out_file, "\n")
print(summary_dt)
