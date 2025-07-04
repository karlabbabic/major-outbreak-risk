# Plot to visualise minor vs major outbreaks 

set.seed(42)

# Parameters
R_range <- c(1, 2)
T <- 200  # Total time steps
w <- c(0.2, 0.3, 0.3, 0.1, 0.1)
S <- length(w)
n_plot <- 100       # Fewer simulations for plotting
threshold <- 200    # Outbreak threshold

# Storage
sim_mat <- matrix(0, nrow = n_plot, ncol = T)
outbreak_flag <- numeric(n_plot)

# Simulations
for (sim in 1:n_plot) {
  R_sim <- runif(1, min = R_range[1], max = R_range[2])
  I <- numeric(T)
  I[1] <- 1
  
  for (t in 2:T) {
    lambda <- sum(w[1:min(t - 1, S)] * I[(t - 1):(t - min(t - 1, S) + 1)])
    I[t] <- rpois(1, R_sim * lambda)
  }
  
  sim_mat[sim, ] <- I
  outbreak_flag[sim] <- as.integer(sum(I) > threshold)
}

library(ggplot2)
library(reshape2)

# Convert to long format for ggplot
df <- as.data.frame(sim_mat)
df$sim_id <- 1:n_plot
df$outbreak <- factor(outbreak_flag)
df$outbreak <- factor(outbreak_flag, levels = c(1, 0), labels = c("Yes", "No"))

df_long <- melt(df, id.vars = c("sim_id", "outbreak"))
df_long$time <- as.integer(gsub("V", "", df_long$variable))

# Plot
ggplot(df_long, aes(x = time, y = value, group = sim_id, colour = outbreak)) +
  geom_line(alpha = 0.5, size = 0.9) +
  coord_cartesian(ylim = c(0, 75), xlim = c(0,75)) +
  labs(
    title = "Simulated Epidemic Trajectories",
    x = "Time (days)",
    y = "Incident cases",
    colour = "Outbreak"
  ) +
  scale_colour_manual(values = c("No" = "#fc8383", "Yes" = "steelblue"))+
  theme_minimal()
