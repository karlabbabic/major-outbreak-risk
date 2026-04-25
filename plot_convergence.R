# Plot convergence of PMO estimates for a chosen incidence sequence.
# For the target sequence I, compare ML and TM estimates as the
# number of simulations increases, against the analytic ('true') PMO.

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

source("compute_analytic_pmo.R")

get_running_outbreak_prop <- function(dt, sequence_list) {
  
  dt <- as.data.table(dt)
  
  match_condition <- rep(TRUE, nrow(dt))
  
  for (var in names(sequence_list)) {
    match_condition <- match_condition & (dt[[var]] == sequence_list[[var]])
  }
  
  matched_rows <- which(match_condition)
  
  if (length(matched_rows) == 0) {
    return(data.table(
      master_row = integer(),
      match_index = integer(),
      prop_major_outbreak = numeric()
    ))
  }
  
  outcomes <- dt$outcome_bin[matched_rows]
  
  running_prop <- cumsum(outcomes == 1) / seq_along(outcomes)
  
  data.table(
    master_row = matched_rows,               # total simulations processed
    match_index = seq_along(matched_rows),   # number of matches seen
    prop_major_outbreak = running_prop
  )
}

# Main plotting function
plot_convergence <- function(target,
                             master_dt,
                             width = 0.1,
                             panel_label = NULL,
                             merged_dir = ".") {
  
  len <- length(target)
  col_names <- paste0("I", seq_len(len))
  
  # ML PMOs: load merged file. otherwise load classifier and use predict() on target

  file_suffix <- paste0(len, "_weeks")
  
  ml_df <- fread(file.path(
    merged_dir,
    paste0("merged_pmo_", file_suffix, ".csv")
  ))
  
  for (j in seq_along(target)) {
    ml_df <- ml_df[get(col_names[j]) == target[j]]
  }
  
  ml_df <- ml_df %>%
    select(n_train, pmo_ml) %>%
    rename(sims = n_train, pmo = pmo_ml) %>%
    mutate(method = "ML")
  
  # TM PMOs: compute from master simulations
  sequence_list <- setNames(as.list(target), col_names)
  tm_running <- get_running_outbreak_prop(master_dt, sequence_list)
  
  if (nrow(tm_running) == 0) {
    stop("Sequence not found in master dataset.")
  }
  
  tm_df <- tm_running %>%
    transmute(
      sims = master_row,               # total simulations
      pmo  = prop_major_outbreak,
      match_index = match_index,
      method = "TM"
    )
  
  w <- readRDS("data/serial_interval_weights.rds")
  true_prob <- compute_analytic_outbreak_prob(
    I = target,
    w = w,
    prior = prior_R(0,10)
  )
  
  lower_bound <- max(0, true_prob - width/2)
  upper_bound <- min(1, true_prob + width/2)
  
  # Convergence calculation
  get_convergence <- function(df) {
    df <- df %>% arrange(sims)
    in_band <- df$pmo >= lower_bound & df$pmo <= upper_bound
    cum_consistent <- rev(cumall(rev(in_band)))
    idx <- which(cum_consistent)[1]
    if (!is.na(idx)) df$sims[idx] else NA
  }
  
  critical_n_ml <- get_convergence(ml_df)
  critical_n_tm <- get_convergence(tm_df)
  
  matches_required_tm <- tm_df %>%
    filter(sims == critical_n_tm) %>%
    pull(match_index)
  
  # Combine for plotting
  all_df <- bind_rows(
    ml_df %>% select(sims, pmo, method),
    tm_df %>% select(sims, pmo, method)
  )
  
  # Plot
  p <- ggplot(all_df, aes(x = sims, y = pmo, colour = method)) +
    
    annotate("rect",
             xmin = min(all_df$sims),
             xmax = max(all_df$sims),
             ymin = lower_bound,
             ymax = upper_bound,
             fill = "#ec4d53",
             alpha = 0.25) +
    
    geom_hline(yintercept = true_prob,
               linetype = "dashed",
               colour = "#ec4d53",
               linewidth = 1) +
    
    geom_vline(xintercept = critical_n_ml,
               linetype = "dotted",
               colour = "#4daf4a",
               linewidth = 0.9) +
    
    geom_vline(xintercept = critical_n_tm,
               linetype = "dotted",
               colour = "#377eb8",
               linewidth = 0.9) +
    
    geom_step(linewidth = 1.2) +
    
    scale_x_log10(
      breaks = 10^(0:6),
      labels = scales::comma
    ) +
    
    scale_y_continuous(limits = c(0, 1.05),
                       labels = function(x)
                         formatC(x, format = "fg", digits = 2)) +
    
    scale_colour_manual(values = c(
      "ML" = "#4daf4a",
      "TM" = "#377eb8"
    )) +
    
    labs(
      x = "Number of simulations (log scale)",
      y = "PMO"
    ) +
    
    theme_minimal(base_size = 7) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6)
    )
  
  # Add annotation text
  if (!is.na(critical_n_ml)) {
    p <- p + annotate("text",
                      x = critical_n_ml,
                      y = 0.45,
                      label = paste0(format(critical_n_ml, big.mark=","), " simulations"),
                      colour = "#4daf4a",
                      hjust = -0.05,
                      size = 2)
  }
  
  if (!is.na(critical_n_tm)) {
    p <- p + annotate("text",
                      x = critical_n_tm,
                      y = 0.35,
                      label = paste0(format(critical_n_tm, big.mark=","), 
                                     " simulations,\n",
                                     format(matches_required_tm, big.mark=","), " matches"),
                      colour = "#377eb8",
                      hjust = -0.05,
                      size = 2)
  }
  
  # Print convergence info
  cat("\nSequence:", paste(target, collapse = ", "), "\n")
  cat("ML converged at:", critical_n_ml, "simulations\n")
  cat("TM converged at:", critical_n_tm, "simulations\n")
  cat("TM matches required:", matches_required_tm, "\n\n")
  
  return(p)
}

# example

master_dt <- readRDS("training_data/training_data_pointmass_I5_1mil.rds")
ml_df_external <- read.csv("merged_ml_pmo_3_weeks.csv")

p1 <- plot_convergence(
  c(1,2,0),
  master_dt = master_dt,
  width = 0.1,
  merged_dir = "training_convergence"
)
p1
