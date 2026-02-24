library(tidyverse)
library(stringr)

summarize_global_results <- function(
    results_dir = "results_Global",  
    alpha_grid  = 0.05,              
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
  
  extract_one_file <- function(f, meta_row) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.list(obj) || is.null(obj$summary_mat)) {
      return(NULL)
    }
    
    sm <- obj$summary_mat
                    
    if (!is.matrix(sm)) sm <- as.matrix(sm)
    if (!all(c("global_p", "runtime_sec") %in% rownames(sm))) {
      return(NULL)
    }
    
    methods <- colnames(sm)
    if (is.null(methods)) return(NULL)
    
    p_vals <- as.numeric(sm["global_p", ])
    rt_sec <- as.numeric(sm["runtime_sec", ])

    is_null_scen <- (meta_row$num2 == 0)
    
    tibble(
      method      = methods,
      global_p    = p_vals,
      runtime_sec = rt_sec
    ) %>%
      mutate(
        is_rejected = if_else(is.finite(global_p), global_p <= alpha_grid, FALSE),
        is_null     = is_null_scen
      )
  }

  long_list <- vector("list", nrow(meta))
  pb <- txtProgressBar(min = 0, max = nrow(meta), style = 3)
  
  for (i in seq_len(nrow(meta))) {
    res <- extract_one_file(meta$file[i], meta[i, ])
    if (!is.null(res)) {
      long_list[[i]] <- bind_cols(meta[i, ], res)
    }
    if (i %% 100 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  
  long_all <- bind_rows(long_list)
  
  if (nrow(long_all) == 0) stop("No valid data extracted.")
  
  summary_table <- long_all %>%
    group_by(template, n, p, d, num1A, num1B, num2, is_null, method) %>%
    summarise(
      n_rep = n(),
      n_na  = sum(is.na(global_p)),
      rejection_rate = mean(is_rejected, na.rm = TRUE),
      rt_mean = mean(runtime_sec, na.rm = TRUE),
      rt_sd   = sd(runtime_sec, na.rm = TRUE),
      rt_q1   = quantile(runtime_sec, 0.25, na.rm = TRUE), # 25%
      rt_med  = median(runtime_sec, na.rm = TRUE),         # 50%
      rt_q3   = quantile(runtime_sec, 0.75, na.rm = TRUE), # 75%
      .groups = "drop"
    ) %>%
    mutate(metric_type = if_else(is_null, "Type_I_Error", "Power"))
  return(list(
    long_data = long_all,      
    summary   = summary_table  
  ))
}
