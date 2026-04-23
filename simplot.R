library(dplyr)
library(ggplot2)

# Parameters
set.seed(42)
T <- 200  # Total time steps
w <- readRDS("data/serial_interval_weights.rds")
S <- length(w)
n_plot <- 450    
threshold <- 100    # Outbreak threshold

# Storage
sim_mat <- matrix(0, nrow = n_plot, ncol = T)
outbreak_flag <- numeric(n_plot)
R_vec <- rep(1.5, n_plot)

# Simulations
for (sim in 1:n_plot) {
  R_sim <- R_vec[sim]
  I <- integer(T)
  I[1] <- 1
  for (t in 2:T) {
    lambda <- 0
    for (s in 1:min(t - 1, S)) {
      lambda <- lambda + w[s] * I[t - s]
    }
    I[t] <- rpois(1, R_sim*lambda)
  }
  sim_mat[sim, ] <- I
  outbreak_flag[sim] <- as.integer(sum(I) >= threshold)
}

# Convert to long format for ggplot
df <- as.data.frame(sim_mat)
df$sim_id <- 1:nrow(sim_mat)
df$outbreak <- factor(outbreak_flag)
df$outbreak <- factor(outbreak_flag, levels = c(1, 0), labels = c("Yes", "No"))
df_long <- reshape2::melt(df, id.vars = c("sim_id", "outbreak"))
df_long$time <- as.integer(gsub("V", "", df_long$variable))

df_long <- df_long %>%
  group_by(sim_id) %>%
  arrange(time) %>%
  mutate(tot_cases = cumsum(value)) %>%
  ungroup()

# Plot cumulative epidemic trajectories
ggplot(df_long, aes(x = time, y = tot_cases, group = sim_id, colour = outbreak)) +
  geom_line(alpha = 0.5, linewidth = 0.4) +
  coord_cartesian(ylim = c(0, 200), xlim = c(0, 100)) +
  labs(
    #title = "Simulated Epidemic Trajectories (Cumulative Cases)",
    x = "Time (weeks)",
    y = "Cumulative cases",
    colour = "Major outbreak?"
  ) +
  scale_colour_manual(values = c("No" = "#00a88a", "Yes" = "#f36f5b")) +
  theme_minimal()+
  theme(axis.title  = element_text(size = 8),
        axis.text   = element_text(size = 6),
        legend.position = "none")

# Calculate total infections per simulation
total_infections <- rowSums(sim_mat)

bins <- cut(
  total_infections,
  breaks = c(seq(0, 100, by = 10), Inf),
  right = FALSE,
  labels = c(
    "0-9", "10-19", "20-29", "30-39", "40-49",
    "50-59", "60-69", "70-79", "80-89", "90-99",
    "100+"
  )
)

# Binned plot - final outbreak sizes
df_bins <- data.frame(Bin = bins)
summary_df <- df_bins %>%
  group_by(Bin) %>%
  summarise(Count = n()
  )

# Plot
ggplot(summary_df, aes(x = Bin, y = Count)) +
  geom_bar(stat = "identity", fill = "#00a88a", width = 1) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Total number of cases",
    y = "Number of outbreaks"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 6))
