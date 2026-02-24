library(tidyverse)
library(stringr)

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
