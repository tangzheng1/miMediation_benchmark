load("taxon_level_summary.RData")
load("taxon_level_long.RData")
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
) # define the colors

benchmark_methods_taxon <- c("HIMA", "LDM", "MarZIC", "multimedia", "CRAmed")

target_levels <- c(
  "Complete Null",
  "Exposure-only",
  "Outcome-only",
  "Disjoint (Balanced +/-)",
  "Disjoint (Dominant +)"
)

###### Table ######

###### Table S3 Computation time across various dataset dimensions ######

tbl_runtime_raw <- taxon_level_long %>%
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


###### Figure ######

###### Fig 2 Empirical FDR of taxon-level mediation tests ######

plot_fig2_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>%       
  filter(num2 > 0) %>%            # mediation_signal
  filter(method %in% taxon_benchmark_methods) %>%
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

draw_fdr_barplot <- function(data, d_val) {
  sub_data <- data %>% filter(d == d_val)
  
  ggplot(sub_data, aes(x = num2, y = fdr_mean, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", size = 0.2, width = 0.7) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey", size = 0.8) +
    facet_grid(p_lab ~ n_lab) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 0.8)) +
    scale_fill_manual(values = my_colors) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),       
      legend.position = "bottom",         
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0, face = "bold", size = 16), 
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 14),     
      legend.title = element_text(size = 16)     
    ) +
    labs(
      title = paste0("Balanced (d = ", d_val, ")"),
      x = "Number of True Mediators",
      y = "Empirical FDR",
      fill = "Method"
    )
}

p1 <- draw_fdr_barplot(plot_fig2_data, 0.5) + 
  labs(title = "a. Balanced +/- ")

p2 <- draw_fdr_barplot(plot_fig2_data, 0.9) + 
  labs(title = "b. Dominant +")

final_fig2 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")

###### Fig4 FDR–power tradeoff for taxon-level mediator discovery ######

TARGET_FDR <- 0.05

fig4_data <- taxon_level_summary %>%
  filter(num2 > 0) %>% 
  group_by(template, n, p, d, num1A, num1B, num2, method) %>%
  summarise(
    valid_configs = sum(fdr_mean <= TARGET_FDR, na.rm = TRUE),
    max_power = if (valid_configs > 0) {
      max(power_mean[fdr_mean <= TARGET_FDR], na.rm = TRUE)
    } else {
      NA_real_ 
    },
    .groups = "drop"
  )

plot_fig4_data <- fig4_data %>%
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

draw_fig4_aligned <- function(data, d_val) {
  sub_data <- data %>% filter(d == d_val)
  ggplot(sub_data, aes(x = num2, y = max_power, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", size = 0.2, width = 0.7) +
    facet_grid(p_lab ~ n_lab) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1)) + 
    scale_fill_manual(values = my_colors) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),       
      legend.position = "bottom",         
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0, face = "bold", size = 16),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 14),     
      legend.title = element_text(size = 16)    
    ) +
    
    labs(
      title = paste0("Max Power (FDR <= 0.05) | d = ", d_val),
      y = "Power@FDR<=0.05 ",
      x = "Number of True Mediators", 
      fill = "Method"
    )
}

p1 <- draw_fig4_aligned(plot_fig4_data, d_val = 0.5) + 
  labs(title = "a. Balanced +/-")

p2 <- draw_fig4_aligned(plot_fig4_data, d_val = 0.9) + 
  labs(title = "b. Dominant +")

final_fig4 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")

###### Fig. S2 power of taxon level tests for benchmark methods ######

plot_S2_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>%       
  filter(num2 > 0) %>%            # mediation_signal
  filter(method %in% taxon_benchmark_methods) %>%
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

draw_power_barplot <- function(data, d_val) {
  sub_data <- data %>% filter(d == d_val)
  
  ggplot(sub_data, aes(x = num2, y = power_mean, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", size = 0.2, width = 0.7) +
    facet_grid(p_lab ~ n_lab) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1)) +
    scale_fill_manual(values = my_colors) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),       
      legend.position = "bottom",         
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0, face = "bold", size = 16),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 14),     
      legend.title = element_text(size = 16)     
    ) +
    labs(
      title = paste0("Balanced (d = ", d_val, ")"),
      x = "Number of True Mediators",
      y = "Power",
      fill = "Method"
    )
}

p1 <- draw_power_barplot(plot_S2_data, 0.5) + 
  labs(title = "a. Balanced +/-")

p2 <- draw_power_barplot(plot_S2_data, 0.9) + 
  labs(title = "b. Dominant +")

final_S2 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")


###### Fig S4 Empirical FDR of taxon-level mediation tests including CAMRA ######

plot_S4_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>%       
  filter(num2 > 0) %>%            # mediation_signal
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

p1 <- draw_fdr_barplot(plot_S4_data, 0.5) + 
  labs(title = "a. Balanced +/-")

p2 <- draw_fdr_barplot(plot_S4_data, 0.9) + 
  labs(title = "b. Dominant +")

final_S4 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")

###### Fig S5 Power of taxon-level mediation tests including CAMRA ######

plot_S5_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>%       
  filter(num2 > 0) %>%            # mediation_signal
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

p1 <- draw_power_barplot(plot_S5_data, 0.5) + 
  labs(title = "a. Balanced +/-")

p2 <- draw_power_barplot(plot_S5_data, 0.9) + 
  labs(title = "b. Dominant +")

final_S5 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")
