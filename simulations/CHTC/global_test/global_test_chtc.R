##############################################################################
# Global-level mediation test simulation
# ───────────────────────────────────────
# Benchmark of global mediation tests:
#   CAMRA, LDM-med, CMM, permanovaFL, MODIMA, MedTest
#
# Usage (CHTC):
#   Rscript global_level_simulation.R <template> <n> <p> <num1_A> <num1_B>
#                                      <num2> <d> <seed>
##############################################################################

# ── Command-line arguments ──────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Need 8 parameters: template, n, p, num1_A, num1_B, num2, d, seed")
}

template_in <- args[1]             # e.g., "GALAXYMicrobLiver_study"
n_in        <- as.numeric(args[2]) # sample size: 200/400/800
p_in        <- as.numeric(args[3]) # number of taxa: 200/400
num1_A_in   <- as.numeric(args[4]) # number of exposure-associated taxa: 10 or 0
num1_B_in   <- as.numeric(args[5]) # number of outcome-associated taxa: 10 or 0
num2_in     <- as.numeric(args[6]) # overlap (true mediators): 0 or 3
d_in        <- as.numeric(args[7]) # positive-effect proportion: 0.5 or 0.9
seed_in     <- as.numeric(args[8]) # replicate seed: 1-500

# ── Load packages ───────────────────────────────────────────────────────────
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
# SECTION 2: CAMRA core functions
##############################################################################

recover_l_PALM <- function(count_m, treat_cov, cov_ad = NULL,
                           prev.filter = 0, eps_p = 1e-10) {
  count_m <- as.matrix(count_m)
  storage.mode(count_m) <- "numeric"
  n <- nrow(count_m); p <- ncol(count_m)
  
  orig_taxa <- colnames(count_m)
  if (is.null(orig_taxa)) orig_taxa <- paste0("O", seq_len(p))
  colnames(count_m) <- orig_taxa
  
  stopifnot(length(treat_cov) == n)
  treat_cov <- matrix(as.numeric(treat_cov), ncol = 1)
  colnames(treat_cov) <- "treat"
  rn <- paste0("T", seq_len(n))
  rownames(count_m) <- rn; rownames(treat_cov) <- rn
  
  if (!is.null(cov_ad)) {
    cov_ad <- as.matrix(cov_ad); stopifnot(nrow(cov_ad) == n)
    storage.mode(cov_ad) <- "numeric"
    rownames(cov_ad) <- rn
    colnames(cov_ad) <- paste0("Cov", seq_len(ncol(cov_ad)))
  }
  
  result1 <- PALM::palm(
    rel.abd = count_m, covariate.interest = treat_cov,
    covariate.adjust = cov_ad, prev.filter = prev.filter
  )
  r1 <- result1$treat
  
  p_full <- rep(1, p); z_full <- rep(0, p); beta_full <- rep(0, p)
  names(p_full) <- names(z_full) <- names(beta_full) <- orig_taxa
  
  feat <- as.character(r1$feature)
  if (length(feat) > 0) {
    idx <- match(feat, orig_taxa); ok <- which(!is.na(idx))
    p_kept <- as.numeric(r1$pval); beta_kept <- as.numeric(r1$coef)
    p_full[idx[ok]] <- p_kept[ok]; beta_full[idx[ok]] <- beta_kept[ok]
    p_adj <- pmax(p_kept[ok], eps_p)
    z_full[idx[ok]] <- stats::qnorm(1 - p_adj / 2) * sign(beta_kept[ok])
  }
  
  list(p = p_full, z = z_full, beta_l = beta_full, feature_kept = feat)
}

recover_r <- function(count_matrix, treat_cov, y, sudo = 0.5, cov_ad = NULL,
                      CClasso = FALSE, cov_true = NULL) {
  logdata  <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
  por_data <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
  
  if (CClasso) {
    est_cov <- {
      res_cov <- fastCCLasso(count_matrix, isCnt = TRUE)
      diag(sqrt(res_cov$cov_diag)) %*% res_cov$rho %*% diag(sqrt(res_cov$cov_diag))
    }
  } else {
    est_cov <- SparCC.count(count_matrix)$cov.w
  }
  if (!is.null(cov_true)) est_cov <- cov_true
  
  p <- ncol(count_matrix); n <- nrow(count_matrix)
  ilr_basis <- compositions::ilrBase(por_data)
  lasso_data_ilr <- as.matrix(compositions::ilr(por_data))
  R2 <- ilr_basis
  Z_ilr <- lasso_data_ilr %*% solve(t(R2) %*% est_cov %*% R2) %*% t(R2) %*% est_cov
  
  # FWL residualization: treatment (+ confounders) unpenalized
  if (is.null(dim(treat_cov))) {
    treat_df <- data.frame(treat = as.numeric(treat_cov))
  } else {
    treat_df <- as.data.frame(treat_cov)
    if (nrow(treat_df) != n) stop("treat_cov nrow mismatch with count_matrix")
  }
  Zdf <- if (!is.null(cov_ad)) {
    cov_df <- as.data.frame(cov_ad)
    if (nrow(cov_df) != n) stop("cov_ad nrow mismatch with count_matrix")
    cbind(treat_df, cov_df)
  } else { treat_df }
  
  Z0 <- model.matrix(~ ., data = Zdf)
  qrZ <- qr(Z0)
  y_tilde <- as.numeric(y - qr.fitted(qrZ, y))
  X_tilde <- Z_ilr - qr.fitted(qrZ, Z_ilr)
  
  outRidge <- hdi::ridge.proj(x = X_tilde, y = y_tilde)
  all_p  <- outRidge$pval
  beta_r <- as.vector(outRidge$bhat)
  z      <- as.vector(qnorm(1 - all_p / 2) * sign(beta_r))
  
  list(p = all_p, z = z, beta_r = beta_r,
       y_tilde = y_tilde, X_tilde = X_tilde, Z0 = Z0)
}

pre_filter_fun <- function(count_matrix, treat_cov, y,
                           const = 2, seed = 42, sudo = 0.5,
                           cov_ad = NULL, adaptive_L = FALSE) {
  set.seed(seed)
  count_matrix <- as.matrix(count_matrix)
  storage.mode(count_matrix) <- "numeric"
  n <- nrow(count_matrix); p <- ncol(count_matrix)
  
  if (length(treat_cov) != n) stop("treat_cov length != nrow(count_matrix)")
  y <- as.numeric(y); if (length(y) != n) stop("y length != nrow(count_matrix)")
  if (!is.null(cov_ad)) {
    cov_ad <- as.matrix(cov_ad); storage.mode(cov_ad) <- "numeric"
    if (nrow(cov_ad) != n) stop("nrow(cov_ad) != nrow(count_matrix)")
  }
  
  logdata <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
  logdata[logdata < (-10)] <- (-10)
  
  est_cov <- SparCC.count(count_matrix)$cov.w
  por_data <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
  R2 <- compositions::ilrBase(por_data)
  Z_ilr <- (logdata %*% R2) %*% solve(t(R2) %*% est_cov %*% R2) %*% t(R2) %*% est_cov
  
  treat_vec <- as.numeric(treat_cov)
  Z0 <- cbind(treat = treat_vec)
  if (!is.null(cov_ad)) {
    Z0 <- cbind(Z0, cov_ad); colnames(Z0) <- make.names(colnames(Z0), unique = TRUE)
  }
  
  X <- cbind(Z_ilr, Z0)
  pZ <- ncol(Z_ilr); p0 <- ncol(Z0)
  pf <- c(rep(1, pZ), rep(0, p0))
  
  if (isTRUE(adaptive_L)) {
    cvfit <- glmnet::cv.glmnet(x = X, y = y, alpha = 1, penalty.factor = pf,
                               nfolds = 5, type.measure = "mse", standardize = TRUE)
    b <- as.matrix(coef(cvfit, s = "lambda.min"))
    beta_Z <- as.numeric(b)[-1][1:pZ]
    selection_set <- which(beta_Z != 0)
    if (length(selection_set) == 0) selection_set <- order(abs(beta_Z), decreasing = TRUE)[1]
    return(sort(unique(as.integer(selection_set))))
  }
  
  K <- max(1L, min(pZ, floor(const * n / log(max(n, 3)))))
  fit <- glmnet::glmnet(x = X, y = y, alpha = 1, penalty.factor = pf,
                        dfmax = min(K + p0, pZ + p0), nlambda = 500,
                        lambda.min.ratio = 1e-6, standardize = TRUE)
  B <- as.matrix(fit$beta)
  dfZ <- colSums(B[1:pZ, , drop = FALSE] != 0)
  idx <- which(dfZ >= K)[1]; if (is.na(idx)) idx <- ncol(B)
  beta_Z_full <- as.numeric(B[1:pZ, idx])
  ord <- order(abs(beta_Z_full), decreasing = TRUE)
  sort(unique(as.integer(ord[seq_len(min(K, length(ord)))])))
}

p_mediation_maxp <- function(p_alpha, p_beta,
                             pi_alpha0 = NULL, pi_beta0 = NULL,
                             pi_method = c("cp4p", "JC"),
                             weight_method = c("maxp", "product", "indenp")) {
  stopifnot(length(p_alpha) == length(p_beta))
  pi_method     <- match.arg(pi_method)
  weight_method <- match.arg(weight_method)
  
  mix_weights_product <- function(pi_alpha0, pi_beta0) {
    eps <- 1e-8
    pa <- min(max(pi_alpha0, eps), 1 - 1e-6)
    pb <- min(max(pi_beta0,  eps), 1 - 1e-6)
    pi0 <- max(1 - (1 - pa) * (1 - pb), 1e-6)
    c(w00 = (pa * pb) / pi0, w10 = ((1 - pa) * pb) / pi0,
      w01 = (pa * (1 - pb)) / pi0, pi0 = pi0)
  }
  
  mix_weights_maxp <- function(p_alpha, p_beta, pi_alpha0, pi_beta0, pi_method) {
    p_max <- pmax(p_alpha, p_beta)
    if (pi_method == "cp4p") {
      obj <- cp4p::estim.pi0(p_max)
      pi0_hat <- if (!is.null(obj$pi0)) as.numeric(obj$pi0) else mean(unlist(obj))
    } else {
      pi0_hat <- miMediation:::.pi0_JC(qnorm(1 - p_max))
    }
    clip01 <- function(x) min(max(x, 1e-6), 1 - 1e-6)
    pi_alpha0 <- clip01(pi_alpha0); pi_beta0 <- clip01(pi_beta0); pi0_hat <- clip01(pi0_hat)
    w <- pmax(c(w00 = (pi_alpha0 + pi_beta0 - pi0_hat) / pi0_hat,
                w10 = (pi0_hat - pi_alpha0) / pi0_hat,
                w01 = (pi0_hat - pi_beta0)  / pi0_hat), 0)
    w <- w / sum(w); names(w) <- c("w00", "w10", "w01"); w
  }
  
  # Step 1: Clean p-values for null proportion estimation
  p_alpha_pi0 <- p_alpha; p_beta_pi0 <- p_beta
  p_alpha_pi0[!is.finite(p_alpha_pi0)] <- 1; p_beta_pi0[!is.finite(p_beta_pi0)] <- 1
  p_alpha_pi0 <- pmin(pmax(p_alpha_pi0, 0), 1); p_beta_pi0 <- pmin(pmax(p_beta_pi0, 0), 1)
  
  # Step 2: Estimate marginal null proportions
  if (is.null(pi_alpha0) || is.null(pi_beta0)) {
    if (pi_method == "cp4p") {
      grab <- function(x) if (!is.null(x$pi0)) as.numeric(x$pi0) else mean(unlist(x), na.rm = TRUE)
      pi_alpha0 <- grab(cp4p::estim.pi0(p_alpha_pi0))
      pi_beta0  <- grab(cp4p::estim.pi0(p_beta_pi0))
    } else {
      pi_alpha0 <- miMediation:::.pi0_JC(qnorm(1 - p_alpha_pi0))
      pi_beta0  <- miMediation:::.pi0_JC(qnorm(1 - p_beta_pi0))
    }
  }
  eps <- 1e-8
  pi_alpha0 <- min(max(pi_alpha0, eps), 1 - eps)
  pi_beta0  <- min(max(pi_beta0,  eps), 1 - eps)
  
  # Step 3: Valid positions
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  out  <- rep(NA_real_, length(p_alpha))
  if (!any(keep)) return(out)
  p_a <- pmin(pmax(p_alpha[keep], eps), 1 - eps)
  p_b <- pmin(pmax(p_beta[keep],  eps), 1 - eps)
  
  # Step 4: Estimate mixture weights
  if (weight_method == "maxp") {
    w <- mix_weights_maxp(p_a, p_b, pi_alpha0, pi_beta0, pi_method)
  } else if (weight_method == "product") {
    w_raw <- mix_weights_product(pi_alpha0, pi_beta0)
    w_vec <- pmax(w_raw[c("w00", "w10", "w01")], 0); w <- w_vec / sum(w_vec)
    names(w) <- c("w00", "w10", "w01")
  } else {
    w_vec <- c(w00 = pi_alpha0 * pi_beta0, w10 = (1 - pi_alpha0) * pi_beta0,
               w01 = pi_alpha0 * (1 - pi_beta0))
    w_vec <- pmax(w_vec, 0); w <- w_vec / sum(w_vec)
  }
  w00 <- as.numeric(w["w00"]); w10 <- as.numeric(w["w10"]); w01 <- as.numeric(w["w01"])
  
  # Step 5: Grenander-based alternative CDF estimation
  estimate_F1_grenander <- function(p, pi0, eps = 1e-8) {
    p <- p[is.finite(p)]; p <- pmin(pmax(p, 0), 1)
    n <- length(p); stopifnot(n > 0); pi0 <- min(max(pi0, 1e-6), 1 - 1e-6)
    x <- sort(unique(c(0, sort(p), 1))); Fn <- ecdf(p); y <- Fn(x)
    dx <- diff(x); keep <- dx > eps; xR <- x[-1][keep]; dx <- dx[keep]
    yL <- y[-length(y)][keep]; yR <- y[-1][keep]
    s <- (yR - yL) / dx; s_hat <- -Iso::pava(-s, w = dx)
    f1_hat <- pmax((s_hat - pi0) / (1 - pi0), 0)
    area <- sum(f1_hat * dx)
    if (area <= 0) return(function(t) rep(0, length(t)))
    f1_hat <- f1_hat / area
    x_knots <- c(0, xR); F1_cum <- c(0, cumsum(f1_hat * dx))
    function(t) {
      t <- pmin(pmax(t, 0), 1)
      pmin(pmax(approx(x_knots, F1_cum, xout = t, method = "linear",
                       ties = "ordered", rule = 2)$y, 0), 1)
    }
  }
  
  F1a <- estimate_F1_grenander(p_a, pi_alpha0)
  F1b <- estimate_F1_grenander(p_b, pi_beta0)
  
  # Step 6: Compute mixture-null p-values
  t <- pmax(p_a, p_b)
  p_mix <- w00 * (t^2) + w10 * (t * F1a(t)) + w01 * (t * F1b(t))
  out[keep] <- pmin(pmax(p_mix, 0), 1)
  out
}

p_mediation_hdmt_fdr <- function(p_alpha, p_beta, exact_p = 1) {
  stopifnot(length(p_alpha) == length(p_beta))
  n <- length(p_alpha); out <- rep(NA_real_, n)
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  if (!any(keep)) return(out)
  pa <- pmin(pmax(p_alpha[keep], 0), 1); pb <- pmin(pmax(p_beta[keep], 0), 1)
  input <- cbind(pa, pb)
  nullprop <- HDMT::null_estimation(input)
  fdr <- HDMT::fdr_est(nullprop$alpha00, nullprop$alpha01, nullprop$alpha10,
                       nullprop$alpha1, nullprop$alpha2,
                       input_pvalues = input, exact = exact_p)
  out[keep] <- fdr; out
}

CAMRA <- function(count_m, treat_cov, y,
                  sudo = 0.5, cov_ad = NULL, FDR_level = 0.05,
                  pre_filter = FALSE, CClasso = FALSE,
                  cov_true = NULL, seed = 42) {
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  
  select_otu <- seq_len(ncol(count_m))
  if (pre_filter) {
    select_otu <- pre_filter_fun(count_matrix = count_m, treat_cov = treat_cov,
                                 y = y, const = 2, seed = seed, sudo = sudo, cov_ad = cov_ad)
  }
  
  res1 <- recover_l_PALM(count_m, treat_cov, cov_ad = cov_ad)
  res2 <- recover_r(count_m, treat_cov, y, cov_ad = cov_ad,
                    CClasso = CClasso, cov_true = cov_true, sudo = sudo)
  
  p1 <- res1$p; p2 <- res2$p
  p_matrix <- cbind(p1, p2)
  
  rawp.perm <- p_mediation_maxp(p1, p2, pi_method = "cp4p", weight_method = "product")
  p_vec <- p.adjust(rawp.perm, method = "BH")
  p_vec_all <- p_vec
  
  if (pre_filter) {
    p_vec_f <- p.adjust(rawp.perm[select_otu], method = "BH")
    p_vec_all[select_otu] <- p_vec_f
    p_vec_all[-select_otu] <- 1
  }
  
  # HDMT-based FDR with fallback
  tmp_locfdr <- try(
    p_mediation_hdmt_fdr(p_matrix[select_otu, 1], p_matrix[select_otu, 2], exact_p = 1),
    silent = TRUE
  )
  if (inherits(tmp_locfdr, "try-error")) {
    tmp_locfdr <- try(
      p_mediation_hdmt_fdr(p_matrix[select_otu, 1], p_matrix[select_otu, 2], exact_p = 0),
      silent = TRUE
    )
  }
  
  if (inherits(tmp_locfdr, "try-error")) {
    idx_detected <- which(p_vec_all < FDR_level)
  } else {
    p_vec_all[select_otu] <- tmp_locfdr
    idx_detected <- select_otu[which(tmp_locfdr <= FDR_level)]
  }
  
  globalp.perm <- min(p_vec_all, na.rm = TRUE)
  runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
  
  list(idx_detected = idx_detected, fdr_value = p_vec_all,
       runtime_sec = runtime_sec, global_p = globalp.perm,
       beta_l = res1$beta_l, beta_r = res2$beta_r,
       taxa_detected = colnames(count_m)[idx_detected], p_matrix = p_matrix)
}


##############################################################################
# SECTION 3: Competing method wrappers
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
      LDM         = "ldm_sim",
      CMM         = "ccmm_sim",
      permanovaFL = "permanovaFL_sim",
      MODIMA      = "MODIMA_sim",
      MedTest     = "Medtest_sim"
    )
  } else {
    method_map <- list(
      LDM         = "ldm_sim",
      permanovaFL = "permanovaFL_sim",
      MODIMA      = "MODIMA_sim",
      MedTest     = "Medtest_sim"
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
# SECTION 6: CHTC entry point
##############################################################################

tryCatch({
  result_list <- runone_simulation_Global(
    n = n_in, p = p_in,
    num1_A = num1_A_in, num1_B = num1_B_in, num2 = num2_in,
    beta_treat = log(5), beta_outcome = 1, d = d_in,
    template = template_in,
    template_dir = ".", save_dir = ".",
    seed = seed_in, save_rds = FALSE, safe = TRUE
  )
  
  output_filename <- sprintf(
    "template_%s_n_%d_p%d_d%s_num1A_%d_num1B_%d_num2_%d_seed_%d.rds",
    gsub("_study$", "", template_in),
    n_in, p_in, as.character(d_in),
    num1_A_in, num1_B_in, num2_in, seed_in
  )
  
  saveRDS(result_list, file = output_filename)
  cat(sprintf("Success! Results saved to: %s\n", output_filename))
}, error = function(e) {
  cat("Error occurred during simulation:\n")
  message(e)
  quit(status = 1)
})
