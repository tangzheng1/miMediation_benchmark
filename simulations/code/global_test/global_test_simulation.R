##############################################################################
# load packages, source, and data
##############################################################################

pkgs <- c(
  "glmnet", "pscl", "plyr", "hdi", "compositions", "parallel",
  "Iso", "cp4p", "HDMT", "tidyverse", "LDM",
  "harmonicmeanp", "precrec", "multimedia",
  "ccmm", "MarZIC", "HIMA", "PALM", "gtools",
  "MultiMed", "permute", "vegan", "matrixStats", "energy"
)
for (pkg in pkgs) library(pkg, character.only = TRUE)

source("LDM_fun.R")
source("MODIMA.R")
source("MedTest.R")
load("GALAXYMicrobLiver_study.RData")

##############################################################################
# SECTION 1: SparCC covariance estimation
##############################################################################
SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  p <- ncol(x); n <- nrow(x)
  x <- x + 1
  
  cov.w <- cor.w <- matrix(0, p, p)
  indLow <- lower.tri(cov.w, diag = TRUE)
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax)
  
  for (i in 1:imax) {
    y <- t(apply(x, 1, function(row) gtools::rdirichlet(n = 1, alpha = row)))
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin)
    covs[, i] <- cov_cor$cov.w[indLow]
    cors[, i] <- cov_cor$cor.w[indLow]
  }
  
  cov.w[indLow] <- apply(covs, 1, median)
  cor.w[indLow] <- apply(cors, 1, median)
  cov.w <- cov.w + t(cov.w); diag(cov.w) <- diag(cov.w) / 2
  cor.w <- cor.w + t(cor.w); diag(cor.w) <- 1
  
  list(cov.w = cov.w, cor.w = cor.w)
}

SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  x <- log(x); p <- ncol(x)
  TT <- stats::var(x)
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT
  
  rowT0 <- rowSums(T0)
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2)) / (p - 2)
  var.w[var.w < Vmin] <- Vmin
  
  Is <- sqrt(1 / var.w)
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5
  cor.w[cor.w <= -1] <- -1; cor.w[cor.w >= 1] <- 1
  
  Lmat <- diag(rep(p - 2, p)) + 1
  rp <- NULL; cp <- rep(TRUE, p); k <- 0
  
  while (k < kmax && sum(cp) > 3) {
    T02 <- T0; curr_cor.w <- cor.w; diag(curr_cor.w) <- 0
    if (!is.null(rp)) curr_cor.w[rp] <- 0
    
    n_rp <- which.max(abs(curr_cor.w))
    if (abs(curr_cor.w[n_rp]) >= alpha) {
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)))
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp)
      rp <- c(rp, n_rp); T02[rp] <- 0; cp <- (diag(Lmat) > 0)
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]))
      var.w[var.w <= Vmin] <- Vmin
      Is <- sqrt(1 / var.w)
      cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5
      cor.w[cor.w <= -1] <- -1; cor.w[cor.w >= 1] <- 1
    } else { break }
    k <- k + 1
  }
  
  Is <- sqrt(var.w)
  cov.w <- cor.w * Is * rep(Is, each = p)
  list(cov.w = cov.w, cor.w = cor.w)
}

##############################################################################
# SECTION 2: Benchmarking methods
##############################################################################

# ── LDM-med ─────────────────────────────────────────────────────────────────
ldm_sim <- function(count_m, treat_cov, y) {
  rand_id <- paste0(sample(letters, 12), collapse = "")
  mat_name <- paste0("M_mat_", rand_id); meta_name <- paste0("meta_", rand_id)
  assign(mat_name, as.matrix(count_m), envir = .GlobalEnv)
  assign(meta_name, data.frame(trt = treat_cov, Y = y), envir = .GlobalEnv)
  on.exit({
    if (exists(mat_name, envir = .GlobalEnv)) rm(list = mat_name, envir = .GlobalEnv)
    if (exists(meta_name, envir = .GlobalEnv)) rm(list = meta_name, envir = .GlobalEnv)
  }, add = TRUE)
  
  fmla_call <- parse(text = paste0(mat_name, " ~ trt + Y"))[[1]]
  tm <- system.time({
    res <- ldm_new(formula = fmla_call, data = get(meta_name, envir = .GlobalEnv),
                   seed = 67817, fdr.nominal = 0.05, test.mediation = TRUE)
  })
  runtime_sec <- unname(tm["elapsed"])
  p_global <- res$med.p.global.omni
  
  P <- as.matrix(res$p.otu.omni)
  p_joint <- MultiMed::medTest.SBMH(P[1, ], P[2, ], MCP.type = "FDR",
                                    t1 = 0.05 / 2, t2 = 0.05 / 2)
  
  rn <- rownames(P)
  if (!is.null(rn) && all(c("trt", "Y") %in% rn)) {
    p_EM <- as.numeric(P["trt", ]); p_MY <- as.numeric(P["Y", ])
  } else {
    p_EM <- as.numeric(P[1, ]); p_MY <- as.numeric(P[2, ])
  }
  
  taxa_names <- colnames(P)
  if (is.null(taxa_names)) taxa_names <- paste0("taxon_", seq_len(ncol(P)))
  
  det <- res$med.detected.otu.omni
  discoveries <- integer(ncol(P)); names(discoveries) <- taxa_names
  if (is.logical(det)) { discoveries <- as.integer(det) }
  else if (is.numeric(det)) {
    valid <- as.integer(det); valid <- valid[valid >= 1 & valid <= ncol(P)]
    if (length(valid) > 0) discoveries[valid] <- 1L
  } else if (is.character(det)) { discoveries[taxa_names %in% det] <- 1L }
  
  list(discoveries = which(discoveries == 1), p_med = p_joint,
       runtime_sec = runtime_sec, global_p = p_global)
}


# ── CMM ──────────────────────────────────────────────────────────────────────
ccmm_sim <- function(count1, treat1, y1, sudo_count = 0.5,
                     method = c("normal", "boot")) {
  method <- match.arg(method)
  t0 <- proc.time()[["elapsed"]]
  treat1_vec <- as.vector(treat1)
  M <- (count1 + sudo_count) / rowSums(count1 + sudo_count)
  
  if (method == "boot") {
    res_ccmm <- ccmm::ccmm(y = as.numeric(y1), M = M, tr = treat1_vec, n.boot = 500)
    ci <- res_ccmm$IDE.CIs; IDEs <- res_ccmm$IDEs
    z_alpha2 <- qnorm(0.975)
    se <- (ci[2, ] - ci[1, ]) / (2 * z_alpha2)
    bad_se <- !is.finite(se) | se <= 0
    z_val <- IDEs / se; z_val[bad_se] <- NA_real_
    p_val <- 2 * pnorm(-abs(z_val))
    p_adj_cmm <- p.adjust(p_val, method = "fdr")
    idx_cmm <- which(p_adj_cmm < 0.05)
    global_p <- ifelse(res_ccmm$TIDE.CI[1] > 0 | res_ccmm$TIDE.CI[2] < 0, 1e-6, 1)
  } else {
    res_ccmm <- ccmm::ccmm(y1, M, treat1_vec, method.est.cov = "normal")
    se <- sqrt(res_ccmm$Var.IDEs)
    z_val <- res_ccmm$IDEs / se
    p_raw <- 2 * pnorm(-abs(z_val))
    p_adj_cmm <- p.adjust(p_raw, method = "BH")
    idx_cmm <- which(p_adj_cmm < 0.05)
    se_tide <- sqrt(res_ccmm$Var.TIDE)
    global_p <- 2 * pnorm(-abs(res_ccmm$TIDE / se_tide))
  }
  
  runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
  list(discoveries = idx_cmm, p_med = p_adj_cmm,
       runtime_sec = runtime_sec, global_p = global_p)
}


# ── permanovaFL ──────────────────────────────────────────────────────────────

permanovaFL_sim <- function(count_m, treat_cov, y,
                            seed = 67817, n.perm.max = 2000,
                            square.dist = TRUE, center.dist = TRUE) {
  t0_all <- proc.time()
  M <- as.matrix(count_m); storage.mode(M) <- "numeric"; n <- nrow(M)
  stopifnot(length(y) == n)
  
  if (is.data.frame(treat_cov) || is.matrix(treat_cov)) {
    tc <- as.data.frame(treat_cov); stopifnot(nrow(tc) == n, ncol(tc) >= 1)
    trt <- tc[[1]]; Z <- if (ncol(tc) >= 2) tc[, -1, drop = FALSE] else NULL
  } else {
    stopifnot(length(treat_cov) == n); trt <- treat_cov; Z <- NULL
  }
  
  meta_df <- data.frame(trt = as.numeric(trt), Y = as.numeric(y))
  if (!is.null(Z)) {
    Zdf <- as.data.frame(Z); colnames(Zdf) <- make.names(colnames(Zdf), unique = TRUE)
    meta_df <- cbind(meta_df, Zdf)
  }
  
  # Assign to global env (required by permanovaFL formula interface)
  rand_id <- paste0(sample(letters, 12), collapse = "")
  mat_name <- paste0("M_mat_", rand_id); meta_name <- paste0("meta_", rand_id)
  assign(mat_name, M, envir = .GlobalEnv); assign(meta_name, meta_df, envir = .GlobalEnv)
  on.exit({
    if (exists(mat_name, envir = .GlobalEnv)) rm(list = mat_name, envir = .GlobalEnv)
    if (exists(meta_name, envir = .GlobalEnv)) rm(list = meta_name, envir = .GlobalEnv)
  }, add = TRUE)
  
  dist.list <- list(
    bray    = as.matrix(vegan::vegdist(M, method = "bray")),
    jaccard = as.matrix(vegan::vegdist(M, method = "jaccard", binary = TRUE))
  )
  
  fmla_str <- if (!is.null(Z)) {
    paste0(mat_name, " | (", paste(colnames(Zdf), collapse = " + "), ") ~ trt + Y")
  } else { paste0(mat_name, " ~ trt + Y") }
  fmla_call <- parse(text = fmla_str)[[1]]
  
  tm_fit <- system.time({
    res <- permanovaFL(formula = fmla_call, data = get(meta_name, envir = .GlobalEnv),
                       dist.list = dist.list, test.mediation = TRUE,
                       n.perm.max = n.perm.max, seed = seed,
                       square.dist = square.dist, center.dist = center.dist)
  })
  
  global_p <- as.numeric(res$med.p.permanova.omni)
  runtime_sec <- as.numeric((proc.time() - t0_all)[["elapsed"]])
  
  list(runtime_sec = runtime_sec, runtime_fit_sec = as.numeric(tm_fit[["elapsed"]]),
       global_p = global_p, p_each = res$med.p.permanova, fit = res)
}


# ── MODIMA ───────────────────────────────────────────────────────────────────

MODIMA_sim <- function(count_m, treat_cov, y,
                       nrep = 2000, adjust_cov = TRUE,
                       seed = NULL, dist_set = c("bray", "jaccard")) {
  t0_all <- proc.time()
  if (!exists("modima", mode = "function"))
    stop("modima() not found. Please source MODIMA.R first.")
  
  count_m <- as.matrix(count_m); storage.mode(count_m) <- "numeric"; n <- nrow(count_m)
  stopifnot(length(y) == n)
  
  if (is.data.frame(treat_cov) || is.matrix(treat_cov)) {
    tc <- as.data.frame(treat_cov); stopifnot(nrow(tc) == n)
    X <- tc[[1]]; Z <- if (ncol(tc) >= 2) tc[, -1, drop = FALSE] else NULL
  } else { X <- treat_cov; Z <- NULL }
  X <- as.numeric(X); Y <- as.numeric(y)
  
  # Adjust for confounders via FWL residualization
  if (!is.null(Z) && isTRUE(adjust_cov)) {
    Zdf <- as.data.frame(Z)
    W_x <- stats::model.matrix(~ ., data = Zdf)
    W_y <- stats::model.matrix(~ X + ., data = data.frame(X = X, Zdf))
    X_use <- stats::lm.fit(W_x, X)$residuals
    Y_use <- stats::lm.fit(W_y, Y)$residuals
  } else { X_use <- X; Y_use <- Y }
  
  dist_X <- stats::dist(X_use); dist_Y <- stats::dist(Y_use)
  if (any(rowSums(count_m) <= 0)) stop("Some samples have zero library size.")
  
  dist_set <- unique(dist_set)
  D_list <- list()
  if ("bray" %in% dist_set) D_list$BC <- vegan::vegdist(count_m, method = "bray")
  if ("jaccard" %in% dist_set) D_list$JAC <- vegan::vegdist(count_m, method = "jaccard", binary = TRUE)
  
  if (!is.null(seed)) set.seed(seed)
  
  tm_test <- system.time({
    p_each <- vapply(names(D_list), function(nm) {
      modima(exposure = dist_X, mediator = D_list[[nm]], response = dist_Y, nrep = nrep)$p.value
    }, numeric(1))
  })
  
  names(p_each) <- names(D_list)
  global_p <- min(min(p_each) * length(p_each), 1)  # Bonferroni
  runtime_sec <- as.numeric((proc.time() - t0_all)[["elapsed"]])
  
  list(runtime_sec = runtime_sec, runtime_test_sec = as.numeric(tm_test[["elapsed"]]),
       global_p = as.numeric(global_p), p_each = as.numeric(p_each))
}


# ── MedTest ──────────────────────────────────────────────────────────────────

Medtest_sim <- function(count_m, treat_cov, y,
                        nperm = 2000, adjust_cov = TRUE,
                        seed = NULL, dist_set = c("bray", "jaccard")) {
  t0_all <- proc.time()
  if (!exists("MedOmniTest", mode = "function"))
    stop("MedOmniTest() not found. Please source MedTest.R first.")
  
  count_m <- as.matrix(count_m); storage.mode(count_m) <- "numeric"; n <- nrow(count_m)
  stopifnot(length(y) == n)
  
  if (is.data.frame(treat_cov) || is.matrix(treat_cov)) {
    tc <- as.data.frame(treat_cov); stopifnot(nrow(tc) == n)
    X <- tc[[1]]; Z <- if (ncol(tc) >= 2) tc[, -1, drop = FALSE] else NULL
  } else { X <- treat_cov; Z <- NULL }
  X <- as.numeric(X); Y <- as.numeric(y)
  
  if (!is.null(Z) && isTRUE(adjust_cov)) {
    Zdf <- as.data.frame(Z)
    W_x <- stats::model.matrix(~ ., data = Zdf)
    W_y <- stats::model.matrix(~ X + ., data = data.frame(X = X, Zdf))
    X_use <- stats::lm.fit(W_x, X)$residuals
    Y_use <- stats::lm.fit(W_y, Y)$residuals; Z_use <- NULL
  } else { X_use <- X; Y_use <- Y; Z_use <- NULL }
  
  if (any(rowSums(count_m) <= 0)) stop("Some samples have zero library size.")
  
  dist_set <- unique(dist_set)
  m.list <- list()
  if ("bray" %in% dist_set) m.list$BC <- as.matrix(vegan::vegdist(count_m, method = "bray"))
  if ("jaccard" %in% dist_set) m.list$JAC <- as.matrix(vegan::vegdist(count_m, method = "jaccard", binary = TRUE))
  
  if (!is.null(seed)) set.seed(seed)
  
  tm_test <- system.time({
    rslt <- MedOmniTest(x = X_use, y = Y_use, m.list = m.list, z = Z_use, nperm = nperm)
  })
  
  runtime_sec <- as.numeric((proc.time() - t0_all)[["elapsed"]])
  list(runtime_sec = runtime_sec, global_p = as.numeric(rslt$permP),
       p_each = setNames(as.numeric(rslt$margPs), names(m.list)), fit = rslt)
}

##############################################################################
# SECTION 3: CAMRA
##############################################################################
CAMRA <- function(count_m, treat_cov, y,
                  sudo = 0.5, cov_ad = NULL,
                  FDR_level = 0.05,
                  pre_filter = FALSE,
                  CClasso = FALSE,
                  cov_true = NULL,
                  method = c("hdmt", "sbmh"),
                  seed = 42) {
  
  set.seed(seed)
  method <- match.arg(method)
  t0 <- proc.time()[["elapsed"]]
  
  select_otu <- seq_len(ncol(count_m))
  
  if (pre_filter) {
    select_otu <- pre_filter_fun(
      count_matrix = count_m,
      treat_cov    = treat_cov,
      y            = y,
      const        = 2,
      seed         = seed,
      sudo         = sudo,
      cov_ad       = cov_ad
    )
    
    select_otu <- sort(unique(select_otu))
    select_otu <- select_otu[select_otu >= 1 & select_otu <= ncol(count_m)]
  }
  
  res1 <- recover_l_PALM(
    count_m,
    treat_cov,
    cov_ad = cov_ad
  )
  
  res2 <- recover_r(
    count_m,
    treat_cov,
    y,
    cov_ad   = cov_ad,
    CClasso  = CClasso,
    cov_true = cov_true,
    sudo     = sudo
  )
  
  p1 <- res1$p
  p2 <- res2$p
  
  p_matrix <- cbind(p1, p2)

  p_vec_all <- rep(1, length(p1))
  idx_detected <- integer(0)
  globalp.perm <- NA_real_
  
  if (length(select_otu) == 0) {
    
    idx_detected <- integer(0)
    globalp.perm <- NA_real_
    
  } else if (method == "hdmt") {
    
    tmp_q <- try(
      p_mediation_hdmt_fdr(
        p_matrix[select_otu, 1],
        p_matrix[select_otu, 2],
        exact_p = 0
      ),
      silent = TRUE
    )
    
    if (inherits(tmp_q, "try-error")) {
      p_vec_all[select_otu] <- NA_real_
      idx_detected <- integer(0)
      globalp.perm <- NA_real_
    } else {
      p_vec_all[select_otu] <- as.numeric(tmp_q)
      idx_sub <- which(tmp_q <= FDR_level)
      idx_detected <- select_otu[idx_sub]
      globalp.perm <- min(p_vec_all, na.rm = TRUE)
    }
    
  } else if (method == "sbmh") {
    
    tmp_q <- try(
      MultiMed::medTest.SBMH(
        p_matrix[select_otu, 1],
        p_matrix[select_otu, 2],
        MCP.type = "FDR",
        t1 = FDR_level / 2,
        t2 = FDR_level / 2
      ),
      silent = TRUE
    )
    
    if (inherits(tmp_q, "try-error")) {
      p_vec_all[select_otu] <- NA_real_
      idx_detected <- integer(0)
      globalp.perm <- NA_real_
    } else {
      p_vec_all[select_otu] <- as.numeric(tmp_q)
      idx_sub <- which(tmp_q <= FDR_level)
      idx_detected <- select_otu[idx_sub]
      globalp.perm <- min(p_vec_all, na.rm = TRUE)
    }
  }
  
  runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
  
  return(list(
    idx_detected  = idx_detected,
    fdr_value     = p_vec_all,
    runtime_sec   = runtime_sec,
    global_p      = globalp.perm,
    beta_l        = res1$beta_l,
    beta_r        = res2$beta_r,
    taxa_detected = colnames(count_m)[idx_detected],
    p_matrix      = p_matrix,
    p_left        = p1,
    p_right       = p2
  ))
}


##############################################################################
# SECTION 4: Data generation from AA template
##############################################################################

generate_data_from_AA <- function(
    num, p, num1_A, num1_B, num2, AA_real_sorted,
    beta_treat = log(5), beta_y = 1, d = 0.8, sigma_y = 5,
    seed = 5, row_id = NULL,
    lib_size_lambda = 2e7, lib_size = 1e6,
    dm_conc = 1e6, dirichlet_eps = 1e-10, beta_TX = 5) {
  
  stopifnot(num >= 1, p >= 1, num1_A >= 0, num1_A <= p, num1_B >= 0, num1_B <= p)
  stopifnot(num2 >= 0, num2 <= min(num1_A, num1_B), num1_A + num1_B - num2 <= p)
  stopifnot(d >= 0, d <= 1)
  if (nrow(AA_real_sorted) < num) stop("Template has fewer rows than num")
  if (ncol(AA_real_sorted) < p)   stop("Template has fewer columns than p")
  
  set.seed(seed)
  
  abs_base <- if (is.null(row_id)) AA_real_sorted[1:num, 1:p, drop = FALSE]
  else { stopifnot(length(row_id) == num); AA_real_sorted[row_id, 1:p, drop = FALSE] }
  abs_base <- as.matrix(abs_base)
  
  # Construct signal sets
  allowed_taxa <- seq_len(p)
  overlap      <- if (num2 > 0) sample(allowed_taxa, num2) else integer(0)
  remain       <- setdiff(allowed_taxa, overlap)
  treat_only   <- if (num1_A - num2 > 0) sample(remain, num1_A - num2) else integer(0)
  remain2      <- setdiff(remain, treat_only)
  outcome_only <- if (num1_B - num2 > 0) sample(remain2, num1_B - num2) else integer(0)
  S_treat   <- sort(c(overlap, treat_only))
  S_outcome <- sort(c(overlap, outcome_only))
  
  X <- sample(rep(0:1, length.out = num))
  
  # Spike-in exposure effects
  abs_true <- abs_base; up_id <- dn_id <- integer(0)
  if (length(S_treat) > 0) {
    k_up <- round(d * length(S_treat))
    up_id <- if (k_up > 0) sample(S_treat, k_up) else integer(0)
    dn_id <- setdiff(S_treat, up_id); mult <- exp(beta_treat)
    rows_t1 <- which(X == 1)
    if (length(rows_t1) > 0 && length(up_id) > 0)
      abs_true[rows_t1, up_id] <- abs_true[rows_t1, up_id, drop = FALSE] * mult
    rows_t0 <- which(X == 0)
    if (length(rows_t0) > 0 && length(dn_id) > 0)
      abs_true[rows_t0, dn_id] <- abs_true[rows_t0, dn_id, drop = FALSE] * mult
  }
  
  # Outcome-side coefficients
  beta_out <- rep(0, p)
  if (length(S_outcome) > 0) {
    m_out <- length(S_outcome); mags <- runif(m_out, 0, beta_y)
    k_pos <- round(d * m_out)
    pos_ids <- if (k_pos > 0) sample(S_outcome, k_pos) else integer(0)
    neg_ids <- setdiff(S_outcome, pos_ids)
    beta_out[pos_ids] <-  mags[match(pos_ids, S_outcome)]
    beta_out[neg_ids] <- -mags[match(neg_ids, S_outcome)]
  }
  
  log_abs <- log(abs_true + 1); log_abs_filled <- scale(log_abs, scale = FALSE)
  Y <- as.numeric(log_abs_filled %*% beta_out) + beta_TX * X + rnorm(num, sd = sigma_y)
  
  # Sequencing counts via multinomial sampling
  lib_size_vec <- if (!is.null(lib_size_lambda)) rpois(num, lambda = lib_size_lambda)
  else rep(as.integer(lib_size), num)
  
  source_AA <- pmax(abs_true, 0); rs <- rowSums(source_AA)
  count <- matrix(0L, nrow = num, ncol = p)
  for (i in seq_len(num)) {
    depth_i <- as.integer(lib_size_vec[i])
    if (!is.finite(depth_i) || depth_i <= 0) { count[i, ] <- 0L; next }
    base_prob <- if (is.finite(rs[i]) && rs[i] > 0) source_AA[i, ] / rs[i] else rep(1/p, p)
    if (is.finite(dm_conc) && dm_conc > 0) {
      prob_adj <- pmax(base_prob, dirichlet_eps); prob_adj <- prob_adj / sum(prob_adj)
      g <- rgamma(p, shape = dm_conc * prob_adj, rate = 1); prob <- g / sum(g)
    } else { prob <- base_prob }
    count[i, ] <- as.integer(rmultinom(1, size = depth_i, prob = prob))
  }
  colnames(count) <- colnames(source_AA); rownames(count) <- rownames(source_AA)
  
  list(abs_true = abs_true, Y = Y, treat = X, count = count,
       idx1 = sort(overlap), idx2 = setdiff(seq_len(p), sort(overlap)),
       lib_size = lib_size_vec, abs_base = abs_base, log_abs_filled = log_abs_filled,
       sets = list(treat = S_treat, outcome = S_outcome, overlap = overlap),
       effects = list(treat_up_taxa = up_id, treat_dn_taxa = dn_id,
                      beta_treat = beta_treat, mult = exp(beta_treat),
                      beta_outcome_vec = beta_out, beta_TX = beta_TX),
       params = list(num = num, p = p, num1_A = num1_A, num1_B = num1_B, num2 = num2,
                     beta_treat = beta_treat, beta_y = beta_y, d = d,
                     sigma_y = sigma_y, lib_size_lambda = lib_size_lambda, seed = seed))
}


##############################################################################
# SECTION 5: Global-level simulation runner
##############################################################################
CAMRA_HDMT_sim <- function(count_m, treat_cov, y) {
  CAMRA(count_m, treat_cov, y, CClasso = FALSE, method = "hdmt")
}

CAMRA_SBMH_sim <- function(count_m, treat_cov, y) {
  CAMRA(count_m, treat_cov, y, CClasso = FALSE, method = "sbmh")
}
                 
runone_simulation_Global <- function(
    n, p, num1_A, num1_B, num2,
    beta_treat, beta_outcome, d,
    template = "GALAXYMicrobLiver_study",
    template_dir = ".", save_dir = ".",
    seed = 1, safe = FALSE, save_rds = TRUE) {
  
  # Input validation
  stopifnot(n >= 1, p >= 1, num1_A >= 0, num1_A <= p, num1_B >= 0, num1_B <= p)
  stopifnot(num2 <= min(num1_A, num1_B), num1_A + num1_B - num2 <= p, d >= 0, d <= 1)
  
  # Load template
  template <- match.arg(template, c("GALAXYMicrobLiver_study"))
  f <- file.path(template_dir, paste0(template, ".RData"))
  if (!file.exists(f)) stop("Template file not found: ", f)
  e <- new.env(parent = emptyenv()); load(f, envir = e)
  template_mat <- as.matrix(get("AA_real", envir = e))
  if (n > nrow(template_mat)) stop("n exceeds template rows")
  
  # Generate simulated data
  sim_data <- generate_data_from_AA(
    num = n, p = p, num1_A = num1_A, num1_B = num1_B, num2 = num2,
    beta_treat = beta_treat, beta_y = beta_outcome, d = d,
    row_id = sample.int(nrow(template_mat), n),
    seed = seed, AA_real_sorted = template_mat
  )
  
  y         <- as.vector(sim_data$Y)
  treat_cov <- as.vector(sim_data$treat)
  count_m   <- as.matrix(sim_data$count)
  idx_true  <- sort(sim_data$sets$overlap)
  
  # ── Define method set (CMM excluded when p != 200 due to scalability) ──
  if (p == 200) {
  method_map <- list(
    CAMRA_HDMT   = "CAMRA_HDMT_sim",
    CAMRA_SBMH   = "CAMRA_SBMH_sim",
    LDM          = "ldm_sim",
    CMM          = "ccmm_sim",
    permanovaFL  = "permanovaFL_sim",
    MODIMA       = "MODIMA_sim",
    MedTest      = "Medtest_sim"
  )
} else {
  method_map <- list(
    CAMRA_HDMT   = "CAMRA_HDMT_sim",
    CAMRA_SBMH   = "CAMRA_SBMH_sim",
    LDM          = "ldm_sim",
    permanovaFL  = "permanovaFL_sim",
    MODIMA       = "MODIMA_sim",
    MedTest      = "Medtest_sim"
  )
}
  
  methods <- names(method_map)
  res_methods <- vector("list", length(methods)); names(res_methods) <- methods
  err_methods <- setNames(rep(NA_character_, length(methods)), methods)
  
  # ── Run each method ──
  for (m in methods) {
    fn_name <- method_map[[m]]
    if (!exists(fn_name, mode = "function")) {
      msg <- paste0("Method function not found: ", fn_name)
      if (isTRUE(safe)) {
        err_methods[[m]] <- msg
        res_methods[[m]] <- list(global_p = NA_real_, runtime_sec = NA_real_)
        next
      } else { stop(msg) }
    }
    
    fun <- get(fn_name, mode = "function")
    if (isTRUE(safe)) {
      res_methods[[m]] <- tryCatch(
        fun(count_m, treat_cov, y),
        error = function(e) {
          err_methods[[m]] <<- conditionMessage(e)
          list(global_p = NA_real_, runtime_sec = NA_real_)
        }
      )
    } else {
      res_methods[[m]] <- fun(count_m, treat_cov, y)
    }
    
    # Enforce required fields
    if (!is.list(res_methods[[m]]) ||
        !all(c("global_p", "runtime_sec") %in% names(res_methods[[m]]))) {
      msg <- paste0("Method '", fn_name, "' output must contain global_p and runtime_sec.")
      if (isTRUE(safe)) {
        err_methods[[m]] <- paste(err_methods[[m]], msg, sep = " | ")
        res_methods[[m]] <- list(global_p = NA_real_, runtime_sec = NA_real_)
      } else { stop(msg) }
    }
  }
  
  # ── Build 2 x K summary matrix ──
  pvals <- vapply(methods, function(m) as.numeric(res_methods[[m]]$global_p), numeric(1))
  rtsec <- vapply(methods, function(m) as.numeric(res_methods[[m]]$runtime_sec), numeric(1))
  summary_mat <- rbind(global_p = pvals, runtime_sec = rtsec)
  colnames(summary_mat) <- methods
  
  out_list <- list(summary_mat = summary_mat, idx_true = idx_true)
  
  if (isTRUE(save_rds)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    save_filename <- sprintf(
      "template_%s_n_%d_p%d_d%s_num1A_%d_num1B_%d_num2_%d_seed_%d.rds",
      gsub("_study$", "", template), n, p, as.character(d),
      num1_A, num1_B, num2, seed
    )
    saveRDS(out_list, file.path(save_dir, save_filename))
  }
  
  out_list
}

##############################################################################
# Example
##############################################################################

kk <- runone_simulation_Global(
  n = 200, 
  p = 200, 
  num1_A = 10, 
  num1_B = 10, 
  num2 = 0,
  beta_treat = log(5), 
  beta_outcome = 1, 
  d = 0.5,
  template = "GALAXYMicrobLiver_study",
  template_dir = ".", 
  save_dir = ".",
  seed = 1)
kk

