###### Fig. S2 power of taxon level tests for benchmark methods ######

plot_S2_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>%       
  filter(num2 > 0) %>%            # mediation_signal
  filter(method %in% benchmark_methods) %>%
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
