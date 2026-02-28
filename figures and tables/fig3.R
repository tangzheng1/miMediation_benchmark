###### Fig3 Quantile-quantile plots of p-values from global mediation tests when p = 200 ######
load("global_test_long.RData")
load("global_test_summary.RData")

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
      scenario = factor(scenario, levels = c(
        "Complete Null", "Exposure-only", "Outcome-only",
        "Disjoint (Balanced +/-)", "Disjoint (Dominant +)"
      )),
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
