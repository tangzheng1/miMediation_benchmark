load("global_test_long.RData")
load("global_test_summary.RData")
library(tidyverse)
library(ggplot2)
library(patchwork)

my_colors <- c(
  "CAMRA"      = "red",
  "LDM-med"    = "blue",  
  "microHIMA"  = "grey60",  
  "MarZIC"     = "green",  
  "multimedia" = "yellow", 
  "CRAmed"     = "orange",  
  "CMM"        = "darkgreen",
  "PERMANOVA-med" = "purple",
  "MODIMA"     = "pink",
  "MedTest"    = "skyblue"
)

target_levels <- c(
  "Complete Null",
  "Exposure-only",
  "Outcome-only",
  "Disjoint (Balanced)",
  "Disjoint (Imbalanced)"
)

######## Figure ########
###### Fig3 Quantile-quantile plots of p-values from global mediation tests when p = 200 ######

plot3_data <- function(long_data, target_p = 200) {
  df <- long_data 
  clean_df <- df %>%
    filter(p == target_p, num2==0) %>%
    mutate(
      scenario = case_when(
        num1A == 0 & num1B == 0 ~ "Complete Null",
        num1A > 0  & num1B == 0 ~ "Exposure-only",
        num1A == 0 & num1B > 0  ~ "Outcome-only",
        num1A > 0  & num1B > 0 & d == 0.5 ~ "Disjoint (Balanced +/-)",
        num1A > 0  & num1B > 0 & d != 0.5 ~ "Disjoint (Dominant +)",
        TRUE ~ "Other"
      )
    ) %>%
    filter(scenario != "Other") %>%
    mutate(
      scenario = factor(scenario, levels = target_levels),
      n_lab = factor(paste0("n = ", n), levels = c("n = 200", "n = 400", "n = 800"))
    )
  
  qq_df <- clean_df %>%
    group_by(scenario, n_lab, method) %>%
    arrange(global_p) %>%
    mutate(
      expected = -log10(ppoints(n())),
      observed = -log10(pmax(global_p, 1e-5)),  # truncate at 10^-5
      upper = -log10(qbeta(0.025, seq_along(global_p), n() - seq_along(global_p) + 1)),
      lower = -log10(qbeta(0.975, seq_along(global_p), n() - seq_along(global_p) + 1))
    ) %>%
    ungroup()
  
  return(qq_df)
}

draw_plot3 <- function(qq_data, p_val) {
  ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed") +
    geom_point(aes(color = method), size = 1.2, alpha = 1) +
    facet_grid(scenario ~ n_lab, scales = "free") +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 4))) +
    scale_x_continuous(limits = c(0, NA), expand = c(0.05, 0)) +
    scale_y_continuous(limits = c(0, 5), expand = c(0.05, 0)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 14),     
      legend.title = element_text(size = 16)     
    ) +
    labs(
      x = expression(Expected ~ -log[10](italic(P))),
      y = expression(Observed ~ -log[10](italic(P))),
      color = "Method"
    )
}

plot3_qq_data <- plot3_data(global_test_long, target_p = 200)
final_fig3 <- draw_plot3(plot3_qq_data, p_val = 200)


###### FigS3 Quantile-quantile plots of p-values from global mediation tests when p = 400 ######

plotS3_qq_data <- plot3_data(global_test_long, target_p = 400)
final_S3 <- draw_plot3(plot3_qq_data, p_val = 400)

######## Table ########
###### Table S2 Empirical Type I error rates for global mediation tests across mediation-null simulation settings ######

global_test_summary <- global_res$summary
table_s2_data <- global_test_summary %>%
  filter(num2 == 0) %>% # mediation null
  filter(method %in% global_benchmark_methods) %>%
  mutate(
    Scenario = case_when(
      num1A == 0  & num1B == 0  ~ "Complete Null",
      num1A == 10 & num1B == 0  ~ "Exposure-only",
      num1A == 0  & num1B == 10 ~ "Outcome-only",
      num1A == 10 & num1B == 10 & d == 0.5 ~ "Disjoint (Balanced)",
      num1A == 10 & num1B == 10 & d == 0.9 ~ "Disjoint (Imbalanced)"
    )) %>%
  mutate(
    Scenario = factor(Scenario, levels = target_levels)
  ) %>%
  dplyr::select(Scenario, p, n, method, rejection_rate)

table_s2_final <- table_s2_data %>%
  pivot_wider(
    names_from = method,
    values_from = rejection_rate
  ) %>%
  arrange(Scenario, p, n) %>%
  dplyr::select(Scenario, p, n, everything())

###### Table S3 Computation time across various dataset dimensions ######
tbl_runtime_raw <- global_res$long_data %>%
  group_by(method, n, p) %>%
  summarise(
    med = median(runtime_sec, na.rm = TRUE),
    q1  = quantile(runtime_sec, 0.25, na.rm = TRUE),
    q3  = quantile(runtime_sec, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    val_str = sprintf("%.2f (%.2f, %.2f)", med, q1, q3),
    setting = sprintf("n=%d, p=%d", n, p)
  )

table_s3 <- tbl_runtime_raw %>%
  select(method, setting, val_str) %>%
  pivot_wider(names_from = setting, values_from = val_str)
