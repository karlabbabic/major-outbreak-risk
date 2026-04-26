# Plot ML, TM and analytic PMO estimates for selected 3-week (or t-week) 
# incidence sequences at a fixed training size n, using merged ML convergence 
# outputs and TM estimates computed directly from the master simulation dataset.

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

source("compute_analytic_pmo.R")

base_dir <- "training_convergence"

pmo_ml <- fread(file.path("training_convergence", "merged_pmo_3_weeks.csv"))

plot_pmo_selected_sequences <- function(results_selected, title = NULL) {
  df_long <- results_selected %>%
    select(Sequence, Analytic, ML, TM) %>%
    pivot_longer(cols = c("Analytic", "ML", "TM"),
                 names_to = "Method", values_to = "PMO")
  
  ggplot(df_long, aes(x = Sequence, y = PMO)) +
    geom_col(data = subset(df_long, Method == "Analytic"),
             fill = "#ec4d53", width = 0.6, alpha = 0.5) +
    geom_point(data = subset(df_long, Method == "ML"),
               colour = "#4daf4a", size = 2.7, position = position_nudge(x = +0.15)) +
    geom_point(data = subset(df_long, Method == "TM"),
               colour = "blue", size = 2.7, position = position_nudge(x = -0.15)) +
    scale_y_continuous(labels = function(x) formatC(x, format = "fg", digits = 2)) +
    labs(
      #title = title,
      x = "Observed incidence",
      y = "PMO"
      #subtitle = "Analytic (bar), ML (green), TM (blue)"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 7)
    )
}

# Build results using ML + TM (from master simulations)

build_results_for_n <- function(Ys, n_value, w, master_dt) {
  
  master_dt <- as.data.table(master_dt)
  
  # Turn Ys into table
  sel <- bind_rows(lapply(Ys, function(y)
    tibble(I1 = y[1], I2 = y[2], I3 = y[3])
  )) %>%
    mutate(
      Sequence = paste0("{", I1, ",", I2, ",", I3, "}"),
      Sequence = factor(Sequence, levels = Sequence)
    )
  
  # ML PMO at this n
  ml_n <- pmo_ml %>%
    filter(n_train == n_value) %>%
    select(I1, I2, I3, ML = pmo_ml)
  
  # TM at this n (computed from master simulations)
  tm_n <- sel %>%
    rowwise() %>%
    mutate(
      TM = {
        seq_list <- list(I1 = I1, I2 = I2, I3 = I3)
        running <- get_running_outbreak_prop(master_dt, seq_list)
        
        if (nrow(running) == 0) {
          NA_real_
        } else {
          running %>%
            filter(master_row <= n_value) %>%
            slice_tail(n = 1) %>%
            pull(prop_major_outbreak)
        }
      }
    ) %>%
    ungroup() %>%
    select(I1, I2, I3, TM)
  
  # Join ML + TM
  results <- sel %>%
    left_join(ml_n, by = c("I1","I2","I3")) %>%
    left_join(tm_n, by = c("I1","I2","I3"))
  
  # Analytic probability
  results <- results %>%
    rowwise() %>%
    mutate(
      Analytic = compute_analytic_outbreak_prob(
        I = c(I1, I2, I3),
        w = w,
        prior = prior_R(0,10)
      )
    ) %>%
    ungroup()
  
  return(results)
}

# Example usage

n_value <- 300000
w <- readRDS("data/serial_interval_weights.rds")

# Example selected 3-week sequences
Ys <- list(
  c(1, 0, 0),
  c(1, 1, 0),
  c(1, 2, 0),
  c(1, 1, 2)
)

results_selected_n <- build_results_for_n(
  Ys,
  n_value = n_value,
  w = w,
  master_dt = training_data_pointmass_I5_1mil
)

p<-plot_pmo_selected_sequences(
  results_selected_n,
  title = paste("PMO at n =", n_value)
)
p
