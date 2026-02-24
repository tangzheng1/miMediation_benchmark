library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)

summarize_results_folder <- function(
    results_dir = "rds_results", # change to the name of your untarred folder
    alpha_grid  = c(0.01, 0.05, 0.10, 0.20),
    pattern     = "\\.rds$",
    recursive   = FALSE
) {
  
  if (!dir.exists(results_dir)) stop("Directory not found: ", results_dir)
  
  files <- list.files(results_dir, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (length(files) == 0) stop("No .rds files found in: ", results_dir)
  
  cat("Found", length(files), "files. Parsing metadata...\n")

  parse_meta <- function(f) {
    bn <- basename(f)
    re <- "^template_(.+?)_n_(\\d+)_p_?(\\d+)_d_?([0-9.]+)_num1A_(\\d+)_num1B_(\\d+)_num2_(\\d+)_seed_(\\d+)\\.rds$"
    m <- str_match(bn, re)
    
    if (any(is.na(m))) {
      return(tibble(file = f, parse_ok = FALSE))
    }
    
    tibble(
      file = f,
      template = m[2],
      n = as.integer(m[3]),
      p = as.integer(m[4]),
      d = as.numeric(m[5]),
      num1A = as.integer(m[6]),
      num1B = as.integer(m[7]),
      num2  = as.integer(m[8]),
      seed  = as.integer(m[9]),
      parse_ok = TRUE
    )
  }
  
  meta <- map_dfr(files, parse_meta)
  
  if (sum(!meta$parse_ok) > 0) {
    warning("Skipping ", sum(!meta$parse_ok), " files with invalid filenames.")
    meta <- meta %>% filter(parse_ok)
  }
  
  metrics_one_file <- function(f, meta_row) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    
    if (!is.list(obj) || is.null(obj$p_mat) || is.null(obj$idx_true)) {
      return(NULL)
    }
    
    p_mat <- obj$p_mat
    idx_true <- sort(unique(as.integer(obj$idx_true)))
    p_use <- ncol(p_mat)
    idx_true <- idx_true[idx_true >= 1 & idx_true <= p_use]
    n_true <- length(idx_true)
    
    methods <- rownames(p_mat)
    if (is.null(methods)) methods <- paste0("method_", seq_len(nrow(p_mat)))
    
    out <- map_dfr(methods, function(m) {
      pv <- as.numeric(p_mat[m, ])
      pv[!is.finite(pv)] <- 1 
      
      map_dfr(alpha_grid, function(a) {
        disc <- which(pv <= a)
        tp <- sum(disc %in% idx_true)
        fp <- length(disc) - tp
        
        tibble(
          method = m,
          alpha  = a,
          n_true = n_true,
          n_disc = length(disc),
          TP = tp,
          FP = fp,
          
          power = if (n_true > 0) tp / n_true else NA_real_,
          FDR   = if (length(disc) > 0) fp / length(disc) else 0,
          any_fp = if (n_true == 0) (length(disc) > 0) else (fp > 0) 
        )
      })
    })
    
    if (!is.null(obj$runtime)) {
      rt <- obj$runtime
      out <- out %>% mutate(runtime_sec = as.numeric(rt[match(method, names(rt))]))
    } else {
      out$runtime_sec <- NA_real_
    }
    
    bind_cols(meta_row, out)
  }
  
  cat("Processing files to calculate Power & FDR...\n")
  
  long_list <- vector("list", nrow(meta))
  pb <- txtProgressBar(min = 0, max = nrow(meta), style = 3)
  
  for (i in seq_len(nrow(meta))) {
    long_list[[i]] <- metrics_one_file(meta$file[i], meta[i, ])
    if (i %% 100 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  
  long_all <- bind_rows(long_list)
  
  if (nrow(long_all) == 0) stop("No valid data extracted.")
 
  summary_setting <- long_all %>%
    group_by(template, n, p, d, num1A, num1B, num2, method, alpha) %>%
    summarise(
    n_rep = n(),
    power_mean = mean(power, na.rm = TRUE),
    fdr_mean   = mean(FDR, na.rm = TRUE),
    runtime_median = median(runtime_sec, na.rm = TRUE),
    runtime_q1     = quantile(runtime_sec, 0.25, na.rm = TRUE), 
    runtime_q3     = quantile(runtime_sec, 0.75, na.rm = TRUE), 
    runtime_disp   = sprintf("%.2f (%.2f, %.2f)", 
                             median(runtime_sec, na.rm = TRUE), 
                             quantile(runtime_sec, 0.25, na.rm = TRUE),
                             quantile(runtime_sec, 0.75, na.rm = TRUE)),
    type1_error = mean(any_fp, na.rm = TRUE),
    
    .groups = "drop")
  
  return(list(long = long_all, summary = summary_setting))
}

my_colors <- c(
  "CAMRA"      = "red",
  "LDM-med"  = "blue",  
  "microHIMA"    = "grey60",  
  "MarZIC"     = "green",  
  "multimedia" = "yellow", 
  "CRAmed"     = "orange",  
  "CMM"        = "darkgreen",
  "PERMANOVA-med" = "purple",
  "MODIMA"     = "pink",
  "MedTest"    = "skyblue"
)

benchmark_methods_taxon <- c("HIMA", "LDM", "MarZIC", "multimedia", "CRAmed")

target_levels <- c(
  "Complete Null",
  "Exposure-only",
  "Outcome-only",
  "Disjoint (Balanced +/-)",
  "Disjoint (Dominant +)"
)

load("taxon_level_summary.RData")
load("taxon_level_long.RData")

###### Table ######

###### Table S3 Any-discoverate for benchmark methods under 0.05 nominal level ######

table_s3_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>% # fixed nominal level at 0.05
  filter(num2 == 0) %>% # mediation null
  filter(method %in% taxon_benchmark_methods) %>%
  mutate(
    Scenario = case_when(
      num1A == 0  & num1B == 0  ~ "Complete Null",
      num1A == 10 & num1B == 0  ~ "Exposure-only",
      num1A == 0  & num1B == 10 ~ "Outcome-only",
      num1A == 10 & num1B == 10 & d == 0.5 ~ "Disjoint (Balanced +/-)",
      num1A == 10 & num1B == 10 & d == 0.9 ~ "Disjoint (Dominant +)"
    )) %>%
  mutate(
    Scenario = factor(Scenario, levels = target_levels)
  ) %>%
  dplyr::select(Scenario, p, n, method, type1_error)

table_s3_final <- table_s3_data %>%
  pivot_wider(
    names_from = method,
    values_from = type1_error
  ) %>%
  arrange(Scenario, p, n) %>%
  dplyr::select(Scenario, p, n, everything())


###### Table S4 Any-false-discovery rate of CAMRA ######

table_s4_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>% # fixed nominal level at 0.05
  filter(num2 == 0) %>% # mediation null
  filter(method == "CAMRA") %>%
  mutate(
    Scenario = case_when(
      num1A == 0  & num1B == 0  ~ "Complete Null",
      num1A == 10 & num1B == 0  ~ "Exposure-only",
      num1A == 0  & num1B == 10 ~ "Outcome-only",
      num1A == 10 & num1B == 10 & d == 0.5 ~ "Disjoint (Balanced +/-)",
      num1A == 10 & num1B == 10 & d == 0.9 ~ "Disjoint (Dominant +)"
    )) %>%
  mutate(
    Scenario = factor(Scenario, levels = target_levels)
  ) %>%
  dplyr::select(Scenario, p, n, method, type1_error)

table_s4_final <- table_s4_data %>%
  pivot_wider(
    names_from = method,
    values_from = type1_error
  ) %>%
  arrange(Scenario, p, n) %>%
  dplyr::select(Scenario, p, n, everything())

###### Table S5 Computation time across various dataset dimensions ######

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

table_s5 <- tbl_runtime_raw %>%
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

