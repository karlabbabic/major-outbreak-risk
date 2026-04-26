# Plot underreporting results for a given prior on I1
# (i) bar plots of the posterior distribution over latent incidence
# sequences P(I | Y, p) for a chosen reporting probability, 
# (ii) PMO comparison panels across reporting probabilities, observed sequences,
# and estimation methods.

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)

# i. Plot posterior probs for sequences for a given observed sequence and reporting prob
Y <- c(1, 2, 0) # Observed sequence
p_pick <- 0.5 # reporting prob

folder <- "combined_outputs/underreporting"
dt <- fread(file.path(
  folder, 
  sprintf("posterior_rows_over_p_Y_%s.csv", paste(Y, collapse = "_"))
))

# corresponding to middle panels in underreporting figure
plot_dt <- dt %>%
  filter(p == p_pick) %>%
  arrange(desc(P_I_given_Y_pois1)) %>%
  slice_head(n = 10) %>%
  mutate(
    I_label = paste0("{", I1, ",", I2, ",", I3, "}"),
    I_label = factor(I_label, levels = I_label)
  )

ggplot(plot_dt, aes(x = I_label, y = P_I_given_Y_pois1)) +
  geom_col(fill = "#169bd5", width = 0.8) +
  labs(
    x = "True incidence data",
    y = "Posterior probability"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12))
  )

# ii. PMO outputs over a range of reporting probabilities and inputs

# load summary files
files <- file.path(
  folder,
  c(
    "summary_over_p_Y_1_0_0.csv",
    "summary_over_p_Y_1_1_0.csv",
    "summary_over_p_Y_1_2_0.csv",
    "summary_over_p_Y_1_1_2.csv"
  )
)

read_one_summary <- function(f) {
  df <- read.csv(f, check.names = FALSE)
  
  has_Y_cols <- all(c("Y1", "Y2", "Y3") %in% names(df))
  
  if (has_Y_cols) {
    df %>%
      mutate(
        Y = paste(Y1, Y2, Y3, sep = "_"),
        Y_lab = paste0("{", Y1, ",", Y2, ",", Y3, "}")
      )
  } else {
    y_string <- sub("^summary_over_p_Y_(.*)\\.csv$", "\\1", basename(f))
    
    df %>%
      mutate(
        Y = y_string,
        Y_lab = paste0("{", gsub("_", ",", y_string), "}")
      )
  }
}

summary_df <- bind_rows(lapply(files, read_one_summary))

plot_df <- summary_df %>%
  pivot_longer(
    cols = matches("^(analytic|ml|tm)_pmo_pois[123]$"),
    names_to = "name",
    values_to = "PMO"
  ) %>%
  mutate(
    pmo_source = case_when(
      grepl("^analytic_", name) ~ "Analytic",
      grepl("^ml_", name) ~ "ML",
      grepl("^tm_", name) ~ "TM"
    ),
    method = case_when(
      grepl("pois1$", name) ~ "pois_1",
      grepl("pois2$", name) ~ "pois_2",
      grepl("pois3$", name) ~ "pois_3"
    )
  )
Y_colours <- c(
  "{1,1,2}" = "#c77cff",
  "{1,2,0}" = "#1f97d4",
  "{1,1,0}" = "#6a9f2e",
  "{1,0,0}" = "#f26a61"
)

make_pmo_panel <- function(plot_df,
                           method_pick = c("pois_1", "pois_2", "pois_3"),
                           circle_p = 0.5,
                           line_width = 2,
                           base_size = 11) {
  
  method_pick <- match.arg(method_pick)
  
  df_sub <- plot_df %>%
    filter(method == method_pick)
  
  # Analytic line for each Y
  df_line <- df_sub %>%
    filter(pmo_source == "Analytic")
  
  # ML/TM points only at plotted grid values
  df_point <- df_sub %>%
    filter(pmo_source %in% c("ML", "TM")) %>%
    filter(abs(p * 20 - round(p * 20)) < 1e-8) %>%   # keep multiples of 0.05
    filter(p >= 0.3 - 1e-8, p <= 1 + 1e-8)
  
  # Right-end labels
  label_df <- df_line %>%
    group_by(Y_lab) %>%
    slice_max(order_by = p, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Optional highlighted circle p = 0.5 on [1,2,0]
  circle_df <- df_point %>%
    filter(Y == "1_2_0", pmo_source == "TM") %>%
    filter(abs(p - circle_p) == min(abs(p - circle_p))) %>%
    slice(1)
  
  ggplot() +
    geom_line(
      data = df_line,
      aes(x = p, y = PMO, colour = Y_lab, group = Y_lab),
      linewidth = line_width
    ) +
    geom_text(
      data = label_df,
      aes(x = p, y = PMO, label = Y_lab, colour = Y_lab),
      hjust = 0,
      nudge_x = 0.02
    ) +
    scale_colour_manual(values = Y_colours, guide = "none") +
    ggnewscale::new_scale_colour() +
    geom_point(
      data = circle_df,
      aes(x = p, y = PMO),
      shape = 21,
      size = 5,
      stroke = 5,
      colour = "goldenrod2",
      fill = NA
    ) +
    # TM underneath
    geom_point(
      data = df_point %>% filter(pmo_source == "TM"),
      aes(x = p, y = PMO),
      colour = "#004d7b",
      alpha = 0.75,
      size = 8,
      shape = 16
    ) +
    
    # ML on top, with white outline
    geom_point(
      data = df_point %>% filter(pmo_source == "ML"),
      aes(x = p, y = PMO),
      fill = "#00a448",
      colour = "white",
      alpha = 0.75,
      size = 4,
      stroke = 0.8,
      shape = 21
    ) +
    
    scale_x_continuous(
      breaks = seq(0.3, 1, by = 0.1),
      limits = c(0.3, 1.03),
      labels = function(x) formatC(x, format = "fg", digits = 2)
    ) +
    scale_y_continuous(
      limits = c(0, 1.05),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) formatC(x, format = "fg", digits = 2)
    ) +
    labs(
      x = expression(paste("Reporting probability (", italic(p), ")")),
      y = "PMO"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      #panel.grid.minor = element_line(linewidth = 2),
      # axis.title = element_text(size = 16),
      # axis.text = element_text(size = 14)
      # ,
      plot.margin = margin(5.5, 30, 5.5, 5.5)
    ) +
    coord_cartesian(clip = "off")
}

p_c <- make_pmo_panel(plot_df, method_pick = "pois_1")
p_f <- make_pmo_panel(plot_df, method_pick = "pois_2")
p_i <- make_pmo_panel(plot_df, method_pick = "pois_3")

p_c
p_f
p_i
