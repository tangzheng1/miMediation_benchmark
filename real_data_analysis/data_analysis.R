load("real_data.RData")

fastCCLasso <- function(xx, isCnt = FALSE, pseudo = 0.5, k_cv = 3,
                        lam_min_ratio = 1e-4, k_max = 20,
                        n_boot = 100, aa = NULL, bb = NULL) {
  n <- nrow(xx)
  p <- ncol(xx)
  
  if (isCnt) {
    xx <- xx + pseudo
    xx <- xx / rowSums(xx)
  }
  
  xx2  <- log(xx) - rowMeans(log(xx))
  vxx2 <- stats::var(xx2)
  
  if (is.null(aa)) aa <- rep(1, p)
  if (is.null(bb)) bb <- 1 / diag(vxx2)
  
  # --- Golden-section search for lambda (log10 scale) ---
  xx_w <- vxx2 * (aa * rep(bb, each = p) + bb * rep(aa, each = p)) / 2
  diag(xx_w) <- 0
  lam_max  <- max(abs(xx_w))
  lam_int2 <- log10(lam_max * c(lam_min_ratio, 1))
  a1 <- lam_int2[1]
  b1 <- lam_int2[2]
  
  lams  <- NULL
  fvals <- NULL
  
  # Two initial trial points
  a2  <- a1 + 0.382 * (b1 - a1)
  b2  <- a1 + 0.618 * (b1 - a1)
  fb2 <- cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2,
                       aa = aa, bb = bb)
  lams  <- c(lams, b2)
  fvals <- c(fvals, fb2)
  fa2 <- cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2,
                       aa = aa, bb = bb)
  lams  <- c(lams, a2)
  fvals <- c(fvals, fa2)
  
  # Convergence tolerances
  err_lam2 <- 1e-1 * max(1, lam_int2)
  err_fval <- 1e-4
  err <- b1 - a1
  k   <- 0
  
  while (err > err_lam2 && k < k_max) {
    fval_max <- max(fa2, fb2)
    if (fa2 > fb2) {
      a1  <- a2;  a2  <- b2;  fa2 <- fb2
      b2  <- a1 + 0.618 * (b1 - a1)
      fb2 <- cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2,
                           aa = aa, bb = bb)
      lams  <- c(lams, b2)
      fvals <- c(fvals, fb2)
    } else {
      b1  <- b2;  b2  <- a2;  fb2 <- fa2
      a2  <- a1 + 0.382 * (b1 - a1)
      fa2 <- cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2,
                           aa = aa, bb = bb)
      lams  <- c(lams, a2)
      fvals <- c(fvals, fa2)
    }
    fval_min <- min(fa2, fb2)
    k   <- k + 1
    err <- b1 - a1
    if (abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) break
  }
  
  info_cv <- list(lams = lams, fvals = fvals, k = k + 2,
                  lam_int = 10^c(a1, b1))
  
  lambda  <- 10^((a2 + b2) / 2)
  fit_res <- fastcclasso_sub(lambda = lambda, SS2 = vxx2, aa = aa, bb = bb)
  
  # Bootstrap inference
  sigma_mod <- boot_fastCCLasso(
    xx2 = xx2, sigma_hat = fit_res$sigma,
    lambda = lambda, aa = aa, bb = bb,
    n_boot = n_boot, max_iter = 200, stop_eps = 1e-6
  )
  
  list(
    rho         = sigma_mod$cor_w,
    cov_diag    = sigma_mod$var_w,
    lambda_best = lambda,
    info_cv     = info_cv,
    p_vals      = sigma_mod$p_vals
  )
}


# ── CV loss for a single lambda ──────────────────────────────────────────────

cvfastCCLasso <- function(lambda, k_cv, xx2, aa, bb) {
  n   <- nrow(xx2)
  p   <- ncol(xx2)
  n_b <- floor(n / k_cv)
  cv.loss <- 0
  
  for (k in 1:k_cv) {
    ite    <- (n_b * (k - 1) + 1):(n_b * k)
    vxx2te <- stats::var(xx2[ite, ])
    vxx2tr <- stats::var(xx2[-ite, ])
    out    <- fastcclasso_sub(lambda = lambda, SS2 = vxx2tr, aa = aa, bb = bb)
    mm     <- out$sigma - out$ww - rep(out$ww, each = p) - vxx2te
    cv.loss <- cv.loss + mean(mm^2 * aa * rep(bb, each = p))
  }
  cv.loss
}


# ── Core solver for a single lambda ─────────────────────────────────────────

fastcclasso_sub <- function(lambda, SS2, aa, bb, k_max = 200, x_tol = 1e-4) {
  p    <- ncol(SS2)
  cc   <- 1 / (aa * sum(bb) + bb * sum(aa))
  aa2  <- aa * cc
  bb2  <- bb * cc
  cab1 <- 1 + sum(aa * bb2)
  caa  <- sum(aa * aa2)
  cbb  <- sum(bb * bb2)
  aabb <- aa * rep(bb, each = p) + bb * rep(aa, each = p)
  lambda2 <- 2 * lambda / aabb
  ss2  <- rowSums(SS2 * aabb)
  sigma <- SS2
  ww   <- colMeans(sigma) - mean(sigma) / 2
  k    <- 0
  err  <- 1
  
  while (err > x_tol && k < k_max) {
    # Update ww
    xx  <- rowSums(sigma * aabb) - ss2
    ax1 <- sum(aa2 * xx)
    bx1 <- sum(bb2 * xx)
    ww2 <- xx * cc + (aa2 * (cbb * ax1 - cab1 * bx1) +
                        bb2 * (caa * bx1 - cab1 * ax1)) / (cab1^2 - caa * cbb)
    
    # Update sigma with soft-thresholding
    sigma2 <- SS2 + ww2 + rep(ww2, each = p)
    oo     <- diag(sigma2)
    sigma2 <- (sigma2 > lambda2) * (sigma2 - lambda2) +
      (sigma2 < -lambda2) * (sigma2 + lambda2)
    diag(sigma2) <- oo
    
    # Check convergence
    err   <- max(abs(sigma2 - sigma) / (abs(sigma) + 1))
    k     <- k + 1
    sigma <- sigma2
  }
  
  list(sigma = sigma, ww = ww2)
}


# ── Bootstrap inference for fastCCLasso ──────────────────────────────────────

boot_fastCCLasso <- function(xx2, sigma_hat, lambda, aa, bb,
                             n_boot = 100, max_iter = 200, stop_eps = 1e-6) {
  n <- nrow(xx2)
  p <- ncol(xx2)
  
  cors_boot <- matrix(0, nrow = p * (p - 1) / 2, ncol = n_boot + 1)
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1)
  cors_mat  <- matrix(0, p, p)
  ind_low   <- lower.tri(cors_mat)
  
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = TRUE),
                     ncol = n_boot)
  
  for (k in 1:n_boot) {
    ind_samp <- sam_boot[, k]
    S_samp   <- stats::var(xx2[ind_samp, ])
    cov_est  <- fastcclasso_sub(lambda, SS2 = S_samp, aa = aa, bb = bb,
                                k_max = 200, x_tol = stop_eps)
    vars_boot[, k] <- diag(cov_est$sigma)
    Is <- 1 / sqrt(vars_boot[, k])
    cor_est <- Is * cov_est$sigma * rep(Is, each = p)
    cors_boot[, k] <- cor_est[ind_low]
  }
  
  # Original estimate
  vars_boot[, n_boot + 1] <- diag(sigma_hat)
  Is <- 1 / sqrt(vars_boot[, n_boot + 1])
  cor_est <- Is * sigma_hat * rep(Is, each = p)
  cors_boot[, n_boot + 1] <- cor_est[ind_low]
  
  # Aggregate bootstrap results
  vars2    <- rowMeans(vars_boot)
  cors2mod <- rowMeans(cors_boot)
  
  cors2_mat <- diag(p)
  cors2_mat[ind_low] <- cors2mod
  cors2_mat <- t(cors2_mat)
  cors2_mat[ind_low] <- cors2mod
  
  # Two-sided p-values via t-distribution
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2)
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals)
  
  pval_mat <- diag(p)
  pval_mat[ind_low] <- p_vals
  pval_mat <- t(pval_mat)
  pval_mat[ind_low] <- p_vals
  
  list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat)
}


# ── SparCC (count-based) ────────────────────────────────────────────────────

SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  p <- ncol(x)
  n <- nrow(x)
  
  # Dirichlet prior: add 1 as pseudo-count
  
  x <- x + 1
  
  cov.w <- cor.w <- matrix(0, p, p)
  indLow <- lower.tri(cov.w, diag = TRUE)
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax)
  
  for (i in 1:imax) {
    # Draw fractions from Dirichlet posterior
    y <- t(apply(x, 1, function(xi) gtools::rdirichlet(n = 1, alpha = xi)))
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin)
    covs[, i] <- cov_cor$cov.w[indLow]
    cors[, i] <- cov_cor$cor.w[indLow]
  }
  
  # Median across posterior samples
  cov.w[indLow] <- apply(covs, 1, median)
  cor.w[indLow] <- apply(cors, 1, median)
  
  cov.w <- cov.w + t(cov.w)
  diag(cov.w) <- diag(cov.w) / 2
  cor.w <- cor.w + t(cor.w)
  diag(cor.w) <- 1
  
  list(cov.w = cov.w, cor.w = cor.w)
}


# ── SparCC (fraction-based) ─────────────────────────────────────────────────
#
# Core SparCC algorithm for known fractions.
#
# @param x     n x p fraction matrix (rows sum to ~1)
# @param kmax  max iteration steps
# @param alpha threshold for strong correlation exclusion
# @param Vmin  minimum variance floor
# @return list(cov.w, cor.w)

SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  x <- log(x)
  p <- ncol(x)
  
  # Variation matrix: T0 = var(log(xi / xj))
  TT <- stats::var(x)
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT
  
  # Basic SparCC variance and correlation
  rowT0 <- rowSums(T0)
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2)) / (p - 2)
  var.w[var.w < Vmin] <- Vmin
  
  Is    <- sqrt(1 / var.w)
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5
  cor.w[cor.w <= -1] <- -1
  cor.w[cor.w >=  1] <-  1
  
  # Iterative exclusion of strong correlations
  Lmat <- diag(rep(p - 2, p)) + 1
  rp   <- NULL
  cp   <- rep(TRUE, p)
  k    <- 0
  
  while (k < kmax && sum(cp) > 3) {
    T02 <- T0
    curr_cor.w <- cor.w
    diag(curr_cor.w) <- 0
    if (!is.null(rp)) curr_cor.w[rp] <- 0
    
    # Find strongest remaining correlation
    n_rp <- which.max(abs(curr_cor.w))
    
    if (abs(curr_cor.w[n_rp]) >= alpha) {
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)))
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1
      
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp)
      rp   <- c(rp, n_rp)
      T02[rp] <- 0
      cp <- (diag(Lmat) > 0)
      
      # Re-estimate variance and correlation
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]))
      var.w[var.w <= Vmin] <- Vmin
      
      Is    <- sqrt(1 / var.w)
      cor.w <- (var.w + rep(var.w, each = p) - T0) *
        Is * rep(Is, each = p) * 0.5
      cor.w[cor.w <= -1] <- -1
      cor.w[cor.w >=  1] <-  1
    } else {
      break
    }
    k <- k + 1
  }
  
  Is    <- sqrt(var.w)
  cov.w <- cor.w * Is * rep(Is, each = p)
  
  list(cov.w = cov.w, cor.w = cor.w)
}


# ── Dispatcher for correlation methods ───────────────────────────────────────

Callallmethods <- function(method, xMat, cv_k, lambda_min_ratio = 1e-4,
                           Edge_eps = 1e-4) {
  p <- ncol(xMat)
  S <- var(log(xMat) - rowMeans(log(xMat)))
  lambda_max <- max(max(S - diag(p)), -min(S - diag(p)))
  lambda_min <- lambda_min_ratio * lambda_max
  lambda_int <- c(lambda_min, lambda_max)
  
  if (method == "fastCCLasso") {
    begin_time <- proc.time()
    result     <- fastCCLasso(xx = xMat, lam_min_ratio = lambda_min_ratio,
                              k_cv = cv_k, k_max = 20)
    end_time   <- proc.time()
    result_cor <- result$rho
    
  } else if (method == "SparCC") {
    begin_time <- proc.time()
    result     <- compute_corr_mod(fracs = xMat, iter = 10, th = 0.1)
    end_time   <- proc.time()
    result_cor <- result$Cor.mat
    
  } else if (method == "CCLasso") {
    begin_time <- proc.time()
    result     <- cclasso(xMat, counts = FALSE, pseudo = 0.5, k_cv = cv_k,
                          lam_int = lambda_int, k_max = 20)
    end_time   <- proc.time()
    result_cor <- result$cor_w
    
  } else if (method == "COAT") {
    begin_time <- proc.time()
    result     <- coat(xMat, nFoler = cv_k, soft = 1)
    end_time   <- proc.time()
    result_cor <- result$corr
    
  } else {
    return(message("Method not recognized. Choose: fastCCLasso, SparCC, CCLasso, COAT"))
  }
  
  result_cor[abs(result_cor) < Edge_eps] <- 0
  
  list(
    runtime   = as.numeric((end_time - begin_time)[3]),
    est_lower = result_cor[lower.tri(result_cor)],
    cor_est   = result_cor
  )
}

recover_l_PALM <- function(count_m, treat_cov, cov_ad = NULL,
                           prev.filter = 0, eps_p = 1e-10) {
  count_m <- as.matrix(count_m)
  storage.mode(count_m) <- "numeric"
  n <- nrow(count_m)
  p <- ncol(count_m)
  
  # Assign internal taxa names to preserve ordering
  orig_taxa <- paste0("O", seq_len(p))
  colnames(count_m) <- orig_taxa
  
  stopifnot(length(treat_cov) == n)
  treat_cov <- matrix(as.numeric(treat_cov), ncol = 1)
  colnames(treat_cov) <- "treat"
  
  rn <- paste0("T", seq_len(n))
  rownames(count_m)  <- rn
  rownames(treat_cov) <- rn
  
  if (!is.null(cov_ad)) {
    cov_ad <- data.frame(cov_ad)
    stopifnot(nrow(cov_ad) == n)
    rownames(cov_ad) <- rn
    colnames(cov_ad) <- paste0("Cov", seq_len(ncol(cov_ad)))
  }
  
  # Run PALM
  result1 <- PALM::palm(
    rel.abd            = count_m,
    covariate.interest = treat_cov,
    covariate.adjust   = cov_ad,
    prev.filter        = prev.filter
  )
  
  r1 <- result1$treat
  
  # Initialize full-length output vectors (filtered taxa retain defaults)
  p_full    <- rep(1, p)
  z_full    <- rep(0, p)
  beta_full <- rep(0, p)
  names(p_full) <- names(z_full) <- names(beta_full) <- orig_taxa
  
  # Align kept features back to original order
  feat <- as.character(r1$feature)
  if (length(feat) > 0) {
    idx <- match(feat, orig_taxa)
    ok  <- which(!is.na(idx))
    
    p_kept    <- as.numeric(r1$pval)
    beta_kept <- as.numeric(r1$coef)
    
    p_full[idx[ok]]    <- p_kept[ok]
    beta_full[idx[ok]] <- beta_kept[ok]
    
    p_adj <- pmax(p_kept[ok], eps_p)
    z_full[idx[ok]] <- stats::qnorm(1 - p_adj / 2) * sign(beta_kept[ok])
  }
  
  list(p = p_full, z = z_full, beta_l = beta_full, feature_kept = feat)
}


# ── Step B: Taxon -> Outcome via ridge projection ───────────────────────────

recover_r <- function(count_matrix, treat_cov, y, sudo = 0.5, cov_ad = NULL,
                      CClasso = FALSE, cov_true = NULL) {
  
  # 1) Compositional transform + covariance estimation
  logdata  <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
  por_data <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
  
  CClasso_core <- function(count_m) {
    res_cov <- fastCCLasso(count_m, isCnt = TRUE)
    diag(sqrt(res_cov$cov_diag)) %*% res_cov$rho %*% diag(sqrt(res_cov$cov_diag))
  }
  
  if (CClasso) {
    est_cov <- CClasso_core(count_matrix)
  } else {
    res_l   <- SparCC.count(count_matrix)
    est_cov <- res_l$cov.w
  }
  
  if (!is.null(cov_true)) est_cov <- cov_true
  
  # 2) Build Z_ilr (ILR-based compositional predictor)
  p <- ncol(count_matrix)
  n <- nrow(count_matrix)
  
  ilr_basis      <- compositions::ilrBase(por_data)
  lasso_data_ilr <- as.matrix(compositions::ilr(por_data))
  R2 <- ilr_basis
  
  Z_ilr <- lasso_data_ilr %*%
    solve(t(R2) %*% est_cov %*% R2) %*% t(R2) %*% est_cov
  
  # 3) FWL residualization: treatment (+ confounders) are unpenalized
  if (is.null(dim(treat_cov))) {
    treat_df <- data.frame(treat = as.numeric(treat_cov))
  } else {
    treat_df <- as.data.frame(treat_cov)
    if (nrow(treat_df) != n)
      stop("nrow(treat_cov) != nrow(count_matrix)")
  }
  
  if (!is.null(cov_ad)) {
    cov_df <- as.data.frame(cov_ad)
    if (nrow(cov_df) != n)
      stop("nrow(cov_ad) != nrow(count_matrix)")
    Zdf <- cbind(treat_df, cov_df)
  } else {
    Zdf <- treat_df
  }
  
  Z0  <- model.matrix(~ ., data = Zdf)
  qrZ <- qr(Z0)
  
  y_hat   <- qr.fitted(qrZ, y)
  y_tilde <- as.numeric(y - y_hat)
  X_hat   <- qr.fitted(qrZ, Z_ilr)
  X_tilde <- Z_ilr - X_hat
  
  # 4) Ridge projection on residualized high-dimensional part
  outRidge <- hdi::ridge.proj(x = X_tilde, y = y_tilde)
  
  all_p  <- outRidge$pval
  beta_r <- as.vector(outRidge$bhat)
  z      <- as.vector(qnorm(1 - all_p / 2) * sign(beta_r))
  
  list(
    p        = all_p,
    z        = z,
    beta_r   = beta_r,
    y_tilde  = y_tilde,
    X_tilde  = X_tilde,
    Z0       = Z0,
    X_doubel = Z_ilr
  )
}


# ── Pre-screening via penalized regression ──────────────────────────────────

pre_filter_fun <- function(count_matrix, treat_cov, y,
                           const = 2, seed = 42, sudo = 0.5,
                           cov_ad = NULL, adaptive_L = FALSE) {
  set.seed(seed)
  
  count_matrix <- as.matrix(count_matrix)
  storage.mode(count_matrix) <- "numeric"
  n <- nrow(count_matrix)
  p <- ncol(count_matrix)
  
  if (length(treat_cov) != n) stop("treat_cov length != nrow(count_matrix)")
  y <- as.numeric(y)
  if (length(y) != n) stop("y length != nrow(count_matrix)")
  
  if (!is.null(cov_ad)) {
    cov_ad <- as.matrix(cov_ad)
    storage.mode(cov_ad) <- "numeric"
    if (nrow(cov_ad) != n) stop("nrow(cov_ad) != nrow(count_matrix)")
  }
  
  # Log-ratio transformation
  logdata <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
  logdata[logdata < -10] <- -10
  
  # Covariance estimation via SparCC
  res_l   <- SparCC.count(count_matrix)
  est_cov <- res_l$cov.w
  
  por_data  <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
  ilr_basis <- compositions::ilrBase(por_data)
  R2 <- ilr_basis
  
  Z_ilr <- (logdata %*% R2) %*%
    solve(t(R2) %*% est_cov %*% R2) %*% t(R2) %*% est_cov
  
  # Design matrix: mediators (penalized) + treatment/confounders (unpenalized)
  treat_vec <- as.numeric(treat_cov)
  Z0 <- cbind(treat = treat_vec)
  if (!is.null(cov_ad)) {
    Z0 <- cbind(Z0, cov_ad)
    colnames(Z0) <- make.names(colnames(Z0), unique = TRUE)
  }
  
  X  <- cbind(Z_ilr, Z0)
  pZ <- ncol(Z_ilr)   # number of mediators (penalized)
  p0 <- ncol(Z0)      # number of unpenalized variables
  pf <- c(rep(1, pZ), rep(0, p0))
  
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required for pre_filter_fun().")
  
  # --- Adaptive mode: use cv.glmnet lambda.min ---
  if (isTRUE(adaptive_L)) {
    cvfit <- glmnet::cv.glmnet(
      x = X, y = y, alpha = 1, penalty.factor = pf,
      nfolds = 5, type.measure = "mse", standardize = TRUE
    )
    b        <- as.matrix(coef(cvfit, s = "lambda.min"))
    beta_all <- as.numeric(b)[-1]
    beta_Z   <- beta_all[1:pZ]
    
    selection_set <- which(beta_Z != 0)
    if (length(selection_set) == 0) {
      ord <- order(abs(beta_Z), decreasing = TRUE)
      selection_set <- ord[1]
    }
    return(sort(unique(as.integer(selection_set))))
  }
  
  # --- Fixed-K mode: const * n / log(n) ---
  K_raw <- floor(const * n / log(max(n, 3)))
  K     <- max(1L, min(pZ, K_raw))
  
  fit <- glmnet::glmnet(
    x = X, y = y, alpha = 1, penalty.factor = pf,
    dfmax = min(K + p0, pZ + p0),
    nlambda = 500, lambda.min.ratio = 1e-6, standardize = TRUE
  )
  
  B   <- as.matrix(fit$beta)
  dfZ <- colSums(B[1:pZ, , drop = FALSE] != 0)
  
  idx <- which(dfZ >= K)[1]
  if (is.na(idx)) idx <- ncol(B)
  
  beta_Z_full <- as.numeric(B[1:pZ, idx])
  ord  <- order(abs(beta_Z_full), decreasing = TRUE)
  keep <- ord[seq_len(min(K, length(ord)))]
  
  sort(unique(as.integer(keep)))
}


# ── Mediation p-value: mixture null (maxP / product) ────────────────────────

p_mediation_maxp <- function(p_alpha, p_beta,
                             pi_alpha0 = NULL, pi_beta0 = NULL,
                             pi_method = c("JC", "cp4p"),
                             weight_method = c("maxp", "product", "indenp")) {
  stopifnot(length(p_alpha) == length(p_beta))
  pi_method     <- match.arg(pi_method)
  weight_method <- match.arg(weight_method)
  
  # --- Internal helpers ---
  
  mix_weights_product <- function(pi_alpha0, pi_beta0) {
    eps <- 1e-8
    pa  <- min(max(pi_alpha0, eps), 1 - 1e-6)
    pb  <- min(max(pi_beta0,  eps), 1 - 1e-6)
    pi0 <- max(1 - (1 - pa) * (1 - pb), 1e-6)
    w00 <- (pa * pb) / pi0
    w10 <- ((1 - pa) * pb) / pi0
    w01 <- (pa * (1 - pb)) / pi0
    c(w00 = w00, w10 = w10, w01 = w01, pi0 = pi0)
  }
  
  # Jin-Cai pi0 estimator via characteristic function
  from_miMed_estpi0 <- function(z) {
    xi <- (0:100) / 100
    tmax <- sqrt(log(length(z)))
    tt <- seq(0, tmax, 0.05)
    epsest <- NULL
    for (j in seq_along(tt)) {
      t <- tt[j]
      f <- exp((t * xi)^2 / 2)
      w <- 1 - abs(xi)
      co <- numeric(101)
      for (i in 1:101) co[i] <- mean(cos(t * xi[i] * z))
      epshat <- sum(w * f * co) / sum(w)
      epsest <- c(epsest, epshat)
    }
    tmp <- min(epsest)
    if (tmp > 1) tmp <- 1
    tmp
  }
  
  mix_weights_maxp <- function(p_alpha, p_beta,
                               pi_alpha0 = NULL, pi_beta0 = NULL,
                               pi_method = c("cp4p", "JC")) {
    pi_method <- match.arg(pi_method)
    
    if (is.null(pi_alpha0) || is.null(pi_beta0)) {
      if (pi_method == "cp4p") {
        pa <- cp4p::estim.pi0(p_alpha)
        pb <- cp4p::estim.pi0(p_beta)
        grab <- function(x) {
          if (!is.null(x$pi0)) as.numeric(x$pi0)
          else mean(unlist(x), na.rm = TRUE)
        }
        pi_alpha0 <- grab(pa)
        pi_beta0  <- grab(pb)
      } else {
        z1 <- qnorm(1 - p_alpha)
        z2 <- qnorm(1 - p_beta)
        pi_alpha0 <- from_miMed_estpi0(z1)
        pi_beta0  <- from_miMed_estpi0(z2)
      }
    }
    
    p_max <- pmax(p_alpha, p_beta)
    if (pi_method == "cp4p") {
      obj     <- cp4p::estim.pi0(p_max)
      pi0_hat <- if (!is.null(obj$pi0)) as.numeric(obj$pi0)
      else mean(unlist(obj))
    } else {
      zmax    <- qnorm(1 - p_max)
      pi0_hat <- from_miMed_estpi0(zmax)
    }
    
    clip01 <- function(x) min(max(x, 1e-6), 1 - 1e-6)
    pi_alpha0 <- clip01(pi_alpha0)
    pi_beta0  <- clip01(pi_beta0)
    pi0_hat   <- clip01(pi0_hat)
    
    w00 <- (pi_alpha0 + pi_beta0 - pi0_hat) / pi0_hat
    w10 <- (pi0_hat - pi_alpha0) / pi0_hat
    w01 <- (pi0_hat - pi_beta0)  / pi0_hat
    
    w <- pmax(c(w00, w10, w01), 0)
    w <- w / sum(w)
    names(w) <- c("w00", "w10", "w01")
    w
  }
  
  # --- 1. Clean p-values for pi0 estimation ---
  p_alpha_pi0 <- p_alpha
  p_beta_pi0  <- p_beta
  p_alpha_pi0[!is.finite(p_alpha_pi0)] <- 1
  p_beta_pi0[!is.finite(p_beta_pi0)]   <- 1
  p_alpha_pi0 <- pmin(pmax(p_alpha_pi0, 0), 1)
  p_beta_pi0  <- pmin(pmax(p_beta_pi0,  0), 1)
  
  # --- 2. Estimate marginal pi0 if not provided ---
  if (is.null(pi_alpha0) || is.null(pi_beta0)) {
    if (pi_method == "cp4p") {
      pa <- cp4p::estim.pi0(p_alpha_pi0)
      pb <- cp4p::estim.pi0(p_beta_pi0)
      grab <- function(x) {
        if (!is.null(x$pi0)) as.numeric(x$pi0)
        else mean(unlist(x), na.rm = TRUE)
      }
      pi_alpha0 <- grab(pa)
      pi_beta0  <- grab(pb)
    } else {
      z1 <- qnorm(1 - p_alpha_pi0)
      z2 <- qnorm(1 - p_beta_pi0)
      pi_alpha0 <- from_miMed_estpi0(z1)
      pi_beta0  <- from_miMed_estpi0(z2)
    }
  }
  
  eps <- 1e-8
  pi_alpha0 <- min(max(pi_alpha0, eps), 1 - eps)
  pi_beta0  <- min(max(pi_beta0,  eps), 1 - eps)
  
  # --- 3. Keep only finite pairs ---
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  out  <- rep(NA_real_, length(p_alpha))
  if (!any(keep)) return(out)
  
  p_a <- pmin(pmax(p_alpha[keep], eps), 1 - eps)
  p_b <- pmin(pmax(p_beta[keep],  eps), 1 - eps)
  
  # --- 4. Compute mixture weights ---
  if (weight_method == "maxp") {
    w <- mix_weights_maxp(p_a, p_b,
                          pi_alpha0 = pi_alpha0,
                          pi_beta0  = pi_beta0,
                          pi_method = pi_method)
  } else if (weight_method == "product") {
    w_raw <- mix_weights_product(pi_alpha0, pi_beta0)
    w_vec <- pmax(w_raw[c("w00", "w10", "w01")], 0)
    w     <- w_vec / sum(w_vec)
    names(w) <- c("w00", "w10", "w01")
  } else if (weight_method == "indenp") {
    w_vec <- c(
      w00 = pi_alpha0 * pi_beta0,
      w10 = (1 - pi_alpha0) * pi_beta0,
      w01 = pi_alpha0 * (1 - pi_beta0)
    )
    w_vec <- pmax(w_vec, 0)
    w     <- w_vec / sum(w_vec)
  }
  
  w00 <- as.numeric(w["w00"])
  w10 <- as.numeric(w["w10"])
  w01 <- as.numeric(w["w01"])
  
  # --- 5. Grenander-estimated alternative CDFs ---
  estimate_F1_grenander <- function(p, pi0, eps = 1e-8) {
    p <- p[is.finite(p)]
    p <- pmin(pmax(p, 0), 1)
    n <- length(p)
    stopifnot(n > 0)
    pi0 <- min(max(pi0, 1e-6), 1 - 1e-6)
    
    x  <- sort(unique(c(0, sort(p), 1)))
    Fn <- ecdf(p)
    y  <- Fn(x)
    dx <- diff(x)
    keep_dx <- dx > eps
    xL <- x[-length(x)][keep_dx]
    xR <- x[-1][keep_dx]
    yL <- y[-length(y)][keep_dx]
    yR <- y[-1][keep_dx]
    dx <- xR - xL
    
    s     <- (yR - yL) / dx
    s_hat <- -Iso::pava(-s, w = dx)
    
    f1_hat <- pmax((s_hat - pi0) / (1 - pi0), 0)
    area   <- sum(f1_hat * dx)
    if (area <= 0) return(function(t) rep(0, length(t)))
    f1_hat <- f1_hat / area
    
    x_knots <- c(0, xR)
    F1_cum  <- c(0, cumsum(f1_hat * dx))
    
    function(t) {
      t <- pmin(pmax(t, 0), 1)
      v <- approx(x_knots, F1_cum, xout = t, method = "linear",
                  ties = "ordered", rule = 2)$y
      pmin(pmax(v, 0), 1)
    }
  }
  
  F1a <- estimate_F1_grenander(p_a, pi_alpha0)
  F1b <- estimate_F1_grenander(p_b, pi_beta0)
  
  # --- 6. Mixture null p-value ---
  t     <- pmax(p_a, p_b)
  p_mix <- w00 * (t^2) + w10 * (t * F1a(t)) + w01 * (t * F1b(t))
  p_mix <- pmin(pmax(p_mix, 0), 1)
  
  out[keep] <- p_mix
  out
}


# ── HDMT-based local FDR estimation ─────────────────────────────────────────
#

p_mediation_hdmt_fdr <- function(p_alpha, p_beta, exact_p = 1) {
  stopifnot(length(p_alpha) == length(p_beta))
  n   <- length(p_alpha)
  out <- rep(NA_real_, n)
  
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  if (!any(keep)) return(out)
  
  pa    <- pmin(pmax(p_alpha[keep], 0), 1)
  pb    <- pmin(pmax(p_beta[keep],  0), 1)
  input <- cbind(pa, pb)
  
  nullprop <- HDMT::null_estimation(input)
  fdr <- HDMT::fdr_est(
    nullprop$alpha00, nullprop$alpha01, nullprop$alpha10,
    nullprop$alpha1, nullprop$alpha2,
    input_pvalues = input, exact = exact_p
  )
  
  out[keep] <- fdr
  out
}

# ── CAMRA main function ─────────────────────────────────────────────────────

CAMRA <- function(mediators,
                  treatment,
                  outcome,
                  confounders = NULL,
                  pseudo      = 0.5,
                  FDR_level   = 0.05,
                  hdmt.exact  = 0,
                  screen      = FALSE,
                  const       = 2,
                  CClasso     = FALSE,
                  seed        = 42) {
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  
  # Default: keep all taxa
  select_otu <- seq_len(ncol(mediators))
  
  # Optional pre-screening
  if (isTRUE(screen)) {
    select_otu <- pre_filter_fun(
      count_matrix = mediators,
      treat_cov    = treatment,
      y            = outcome,
      const        = const,
      seed         = seed,
      sudo         = pseudo,
      cov_ad       = confounders
    )
  }
  
  # Step A: exposure -> taxon
  res1 <- recover_l_PALM(
    count_m   = mediators,
    treat_cov = treatment,
    cov_ad    = confounders
  )
  
  # Step B: taxon -> outcome
  res2 <- recover_r(
    count_matrix = mediators,
    treat_cov    = treatment,
    y            = outcome,
    cov_ad       = confounders,
    CClasso      = CClasso,
    sudo         = pseudo
  )
  
  p1 <- res1$p
  p2 <- res2$p
  z1 <- res1$z
  z2 <- as.vector(res2$z)
  
  p_matrix <- cbind(p1, p2)
  rownames(p_matrix) <- colnames(mediators)
  
  # Combine p-values via mixture null (product weights)
  rawp.perm <- p_mediation_maxp(p1, p2, pi_method = "cp4p",
                                weight_method = "product")
  p_vec <- p.adjust(rawp.perm, method = "BH")
  
  rawp.perm.rm <- na.omit(rawp.perm)
  rawp.perm.rm[rawp.perm.rm < 1e-10] <- 1e-10
  
  p_vec_f   <- NULL
  p_vec_all <- p_vec
  
  # If screening, apply BH correction within selected set
  if (screen) {
    p_vec_f <- p.adjust(rawp.perm[select_otu], method = "BH")
    p_vec_all <- p_vec
    p_vec_all[select_otu]  <- p_vec_f
    p_vec_all[-select_otu] <- 1
  }
  
  # HDMT-based FDR estimation (try exact_p = hdmt.exact first, then fallback)
  stopifnot(hdmt.exact %in% c(0, 1))
  exact_first  <- as.integer(hdmt.exact)
  exact_second <- 1L - exact_first
  
  tmp_locfdr <- try(
    p_mediation_hdmt_fdr(
      p_matrix[select_otu, 1],
      p_matrix[select_otu, 2],
      exact_p = exact_first
    ),
    silent = TRUE
  )
  
  if (inherits(tmp_locfdr, "try-error")) {
    tmp_locfdr <- try(
      p_mediation_hdmt_fdr(
        p_matrix[select_otu, 1],
        p_matrix[select_otu, 2],
        exact_p = exact_second
      ),
      silent = TRUE
    )
  }
  
  # Fallback to BH q-values if HDMT fails
  if (inherits(tmp_locfdr, "try-error")) {
    selected_values <- p_vec_all
    idx_detected    <- which(selected_values < FDR_level)
    locfdr <- rep(NA_real_, length(select_otu))
  } else {
    p_vec_all[select_otu] <- tmp_locfdr
    idx_sub      <- which(tmp_locfdr <= FDR_level)
    idx_detected <- select_otu[idx_sub]
  }
  
  globalp     <- min(p_vec_all, na.rm = TRUE)
  runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
  
  list(
    idx_detected  = idx_detected,
    qval.med      = p_vec_all,
    runtime_sec   = runtime_sec,
    global_p      = globalp,
    pval.alpha    = p1,
    pval.beta     = p2,
    beta_l        = res1$beta_l,
    beta_r        = res2$beta_r,
    taxa_detected = colnames(mediators)[idx_detected],
    p_matrix      = p_matrix,
    filter_index  = select_otu
  )
}

# =============================================================================
# LDM-med 
# =============================================================================
ldm_med_sim <- function(count_m, treat_cov, y,
                        fdr_nominal = 0.05, seed = 42) {
  t0 <- proc.time()
  library(LDM)
  
  # Assign to .GlobalEnv with random names to avoid LDM formula-parsing bug
  rand_id   <- paste0(sample(letters, 12), collapse = "")
  mat_name  <- paste0("RA_mat_", rand_id)
  meta_name <- paste0("meta_", rand_id)
  
  assign(mat_name,  data.frame(as.matrix(count_m)), envir = .GlobalEnv)
  assign(meta_name, data.frame(trt = treat_cov, outcome = y), envir = .GlobalEnv)
  
  on.exit({
    if (exists(mat_name,  envir = .GlobalEnv))
      rm(list = mat_name,  envir = .GlobalEnv)
    if (exists(meta_name, envir = .GlobalEnv))
      rm(list = meta_name, envir = .GlobalEnv)
  }, add = TRUE)
  
  # Build formula as a parsed expression so LDM can find the objects
  fmla <- parse(text = paste0(mat_name, " ~ trt + outcome"))[[1]]
  
  set.seed(seed)
  res <- LDM::ldm(
    formula        = fmla,
    data           = get(meta_name, envir = .GlobalEnv),
    fdr.nominal    = fdr_nominal,
    test.mediation = TRUE
  )
  
  # Per-taxon omnibus p-values (2 rows: trt, outcome)
  P <- as.matrix(res$p.otu.omni)
  rn <- rownames(P)
  if (!is.null(rn) && all(c("trt", "outcome") %in% rn)) {
    p_EM <- as.numeric(P["trt", ])
    p_MY <- as.numeric(P["outcome", ])
  } else {
    p_EM <- as.numeric(P[1, ])
    p_MY <- as.numeric(P[2, ])
  }
  
  # Joint mediation p-value: max(p_exposure->mediator, p_mediator->outcome)
  p_med <- pmax(p_EM, p_MY)
  taxa_names <- colnames(P)
  if (is.null(taxa_names))
    taxa_names <- paste0("taxon_", seq_len(ncol(P)))
  names(p_med) <- taxa_names
  
  # Detected mediators from LDM's built-in FDR procedure
  det     <- res$med.detected.otu.omni
  det_sel <- rep(FALSE, ncol(P))
  
  if (is.logical(det)) {
    det_sel <- det
  } else if (is.numeric(det)) {
    valid <- as.integer(det)
    valid <- valid[valid >= 1 & valid <= ncol(P)]
    if (length(valid) > 0) det_sel[valid] <- TRUE
  } else if (is.character(det)) {
    det_sel <- taxa_names %in% det
  }
  
  discoveries  <- which(det_sel)
  otu_detected <- taxa_names[det_sel]
  global_p     <- res$med.p.global.omni
  runtime_sec  <- as.numeric((proc.time() - t0)["elapsed"])
  
  list(
    discoveries  = discoveries,
    otu_detected = otu_detected,
    p_med        = p_med,
    p_otu        = res$p.otu.omni,
    beta         = res$beta,
    global_p     = global_p,
    runtime_sec  = runtime_sec,
    ldm_obj      = res
  )
}

# =============================================================================
# CCM 
# =============================================================================

ccmm_sim <- function(count1, treat1, y1, sudo_count = 0.5,
                     method = c("boot", "normal")) {
  method <- match.arg(method)
  t0 <- proc.time()[["elapsed"]]
  treat1_vec <- as.vector(treat1)
  
  # Helper: derive p-values from confidence intervals
  get_p_from_ci <- function(IDEs, p_matrix_cmm, method = "fdr",
                            ci_level = 0.95) {
    ci_lower <- p_matrix_cmm[1, ]
    ci_upper <- p_matrix_cmm[2, ]
    z_alpha2 <- qnorm(1 - (1 - ci_level) / 2)
    se     <- (ci_upper - ci_lower) / (2 * z_alpha2)
    bad_se <- !is.finite(se) | se <= 0
    z_val  <- IDEs / se
    z_val[bad_se] <- NA_real_
    p_val <- 2 * pnorm(-abs(z_val))
    p_adj <- p.adjust(p_val, method = method)
    data.frame(IDE = IDEs, SE = se, z = z_val,
               p_value = p_val, p_adj = p_adj)
  }
  
  realdata_prop <- (count1 + sudo_count) / rowSums(count1 + sudo_count)
  M <- realdata_prop
  
  if (method == "boot") {
    res_ccmm   <- ccmm::ccmm(y = as.numeric(y1), M = M,
                             tr = treat1_vec, n.boot = 500)
    res_ccmm_p <- get_p_from_ci(res_ccmm$IDEs, res_ccmm$IDE.CIs)
    p_adj_cmm  <- res_ccmm_p$p_adj
    idx_cmm    <- which(p_adj_cmm < 0.05)
    global_p   <- ifelse(
      res_ccmm$TIDE.CI[1] > 0 | res_ccmm$TIDE.CI[2] < 0, 1e-6, 1
    )
  } else {
    res_ccmm  <- ccmm::ccmm(y1, realdata_prop, treat1_vec,
                            method.est.cov = "normal")
    se        <- sqrt(res_ccmm$Var.IDEs)
    z_val     <- res_ccmm$IDEs / se
    p_raw     <- 2 * pnorm(-abs(z_val))
    p_adj_cmm <- p.adjust(p_raw, method = "BH")
    idx_cmm   <- which(p_adj_cmm < 0.05)
    se_tide   <- sqrt(res_ccmm$Var.TIDE)
    global_p  <- 2 * pnorm(-abs(res_ccmm$TIDE / se_tide))
  }
  
  runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
  list(discoveries = idx_cmm, p_med = p_adj_cmm,
       runtime_sec = runtime_sec, global_p = global_p)
}


# =============================================================================
# CRAmed 
# =============================================================================

CRAmed_fullp <- function(M_mat, Y, Exposure,
                         FDR = 0.05, n.times = 100, prefilter = TRUE,
                         n.perm = 100, CI = FALSE,
                         modely = "gaussian", modelm = "ZINB",
                         method = "BH") {
  library(glmnet)
  library(MASS)
  
  # --- Helper: build p-value table ---
  build_p_table <- function(M_mat, selection_set, pvl_lst, final.p,
                            method, mediation_set = integer(0)) {
    p    <- ncol(M_mat)
    taxa <- colnames(M_mat)
    if (is.null(taxa)) taxa <- paste0("taxa_", seq_len(p))
    
    # Non-selected taxa default to p = 1
    p_d0crt <- rep(1, p)
    p_model <- rep(1, p)
    q_d0crt <- rep(1, p)
    q_model <- rep(1, p)
    p_joint <- rep(1, p)
    
    selection_set <- sort(unique(as.integer(selection_set)))
    selection_set <- selection_set[selection_set >= 1 & selection_set <= p]
    
    if (length(selection_set) > 0) {
      p_d0crt[selection_set] <- pvl_lst[selection_set]
      p_model[selection_set] <- as.numeric(final.p)
      q_d0crt[selection_set] <- p.adjust(p_d0crt[selection_set], method = method)
      q_model[selection_set] <- p.adjust(p_model[selection_set], method = method)
      p_joint[selection_set] <- pmax(q_d0crt[selection_set],
                                     q_model[selection_set])
    }
    
    p_table <- data.frame(
      taxa = taxa, in_select = FALSE, is_mediator = FALSE,
      p_d0crt = p_d0crt, p_model = p_model,
      q_d0crt = q_d0crt, q_model = q_model, p_joint = p_joint,
      stringsAsFactors = FALSE
    )
    if (length(selection_set) > 0) p_table$in_select[selection_set] <- TRUE
    if (length(mediation_set) > 0) {
      mediation_set <- mediation_set[mediation_set >= 1 & mediation_set <= p]
      p_table$is_mediator[mediation_set] <- TRUE
    }
    p_table
  }
  
  # --- Conditional generators (ZINB, NB, ZIP) ---
  
  Creat_condition_zinb <- function(M_mat, indx, Exposure, n.times) {
    library(pscl)
    n <- nrow(M_mat)
    Y_data <- as.data.frame(cbind(Exposure, M_mat[, indx]))
    colnames(Y_data) <- c("Exposure", "Mediator")
    
    zinb.fit <- try(
      zeroinfl(Mediator ~ Exposure | Exposure,
               data = Y_data, dist = "negbin", link = "logit"),
      silent = TRUE
    )
    
    if (inherits(zinb.fit, "try-error")) {
      zinb.fit <- try(glm.nb(Mediator ~ Exposure, data = Y_data), silent = TRUE)
      alpha <- summary(zinb.fit)$coefficients[1:2, 1]
      phi   <- zinb.fit$theta
      M_bar <- zinb.fit$fitted.values
      M_res <- zinb.fit$residuals
      x     <- cbind(1, Exposure)
      lamda <- exp(x %*% matrix(alpha))
      
      samp_M <- matrix(NA, n, n.times)
      for (j in 1:n.times) {
        set.seed(j)
        for (i in 1:n) samp_M[i, j] <- rnbinom(1, mu = lamda[i], size = phi)
      }
      res_M_samp <- lapply(seq_len(n.times), function(j) {
        Yd  <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
        fit <- try(glm.nb(Mediator ~ Exposure, data = Yd), silent = TRUE)
        data.frame(t(fit$residuals))
      })
    } else {
      alpha <- summary(zinb.fit)$coefficients$count[1:2, 1]
      gamma <- summary(zinb.fit)$coefficients$zero[1:2, 1]
      phi   <- zinb.fit$theta
      M_bar <- zinb.fit$fitted.values
      M_res <- zinb.fit$residuals
      
      x       <- cbind(1, Exposure)
      logit.p <- x %*% matrix(gamma)
      p0      <- 1 / (1 + exp(-logit.p))
      lamda   <- exp(x %*% matrix(alpha))
      
      samp_M <- matrix(NA, n, n.times)
      for (j in 1:n.times) {
        set.seed(j)
        Z1 <- rbinom(n, 1, p0)
        for (i in 1:n) {
          samp_M[i, j] <- if (Z1[i] == 1) 0
          else rnbinom(1, mu = lamda[i], size = phi)
        }
      }
      res_M_samp <- lapply(seq_len(n.times), function(j) {
        if (sum(samp_M[, j]) == 0) return(data.frame(t(rep(0, n))))
        Yd  <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
        fit <- try(
          zeroinfl(Mediator ~ Exposure | Exposure,
                   data = Yd, dist = "negbin", link = "logit"),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) {
          fit <- try(glm.nb(Mediator ~ Exposure, data = Yd), silent = TRUE)
          if (inherits(fit, "try-error")) {
            fit <- try(
              zeroinfl(Mediator ~ Exposure | Exposure,
                       data = Yd, dist = "poisson", link = "logit"),
              silent = TRUE
            )
          }
        }
        data.frame(t(fit$residuals))
      })
    }
    list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp)
  }
  
  Creat_condition_nb <- function(M_mat, indx, Exposure, n.times) {
    n <- nrow(M_mat)
    Y_data <- data.frame(Exposure = Exposure, Mediator = M_mat[, indx])
    nb.fit <- try(glm.nb(Mediator ~ Exposure, data = Y_data), silent = TRUE)
    
    alpha <- summary(nb.fit)$coefficients[1:2, 1]
    phi   <- nb.fit$theta
    M_bar <- nb.fit$fitted.values
    M_res <- nb.fit$residuals
    x     <- cbind(1, Exposure)
    lamda <- exp(x %*% matrix(alpha))
    
    samp_M <- matrix(NA, n, n.times)
    for (j in 1:n.times) {
      set.seed(j)
      for (i in 1:n) samp_M[i, j] <- rnbinom(1, mu = lamda[i], size = phi)
    }
    res_M_samp <- lapply(seq_len(n.times), function(j) {
      if (sum(samp_M[, j]) == 0) return(data.frame(t(rep(0, n))))
      Yd  <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
      fit <- try(glm.nb(Mediator ~ Exposure, data = Yd), silent = TRUE)
      data.frame(t(fit$residuals))
    })
    list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp)
  }
  
  Creat_condition_zip <- function(M_mat, indx, Exposure, n.times) {
    library(pscl)
    n <- nrow(M_mat)
    Y_data <- data.frame(Exposure = Exposure, Mediator = M_mat[, indx])
    zip.fit <- try(
      zeroinfl(Mediator ~ Exposure | Exposure,
               data = Y_data, dist = "poisson", link = "logit"),
      silent = TRUE
    )
    
    alpha   <- summary(zip.fit)$coefficients$count[1:2, 1]
    gamma   <- summary(zip.fit)$coefficients$zero[1:2, 1]
    M_bar   <- zip.fit$fitted.values
    M_res   <- zip.fit$residuals
    x       <- cbind(1, Exposure)
    logit.p <- x %*% matrix(gamma)
    p0      <- 1 / (1 + exp(-logit.p))
    lamda   <- exp(x %*% matrix(alpha))
    
    samp_M <- matrix(NA, n, n.times)
    for (j in 1:n.times) {
      set.seed(j)
      Z1 <- rbinom(n, 1, p0)
      for (i in 1:n) {
        samp_M[i, j] <- if (Z1[i] == 1) 0 else rpois(1, lambda = lamda[i])
      }
    }
    res_M_samp <- lapply(seq_len(n.times), function(j) {
      if (sum(samp_M[, j]) == 0) return(data.frame(t(rep(0, n))))
      Yd  <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
      fit <- try(
        zeroinfl(Mediator ~ Exposure | Exposure,
                 data = Yd, dist = "poisson", link = "logit"),
        silent = TRUE
      )
      data.frame(t(fit$residuals))
    })
    list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp)
  }
  
  # --- Main body ---
  p <- ncol(M_mat)
  n <- nrow(M_mat)
  
  # Prefilter / selection set
  if (prefilter) {
    set.seed(123)
    cv_lasso <- cv.glmnet(cbind(M_mat, Exposure), Y, alpha = 1,
                          family = modely, dfmax = as.integer(p / 2))
    lamb      <- cv_lasso$lambda.min
    opt_model <- glmnet(cbind(M_mat, Exposure), Y, alpha = 1,
                        lambda = lamb, family = modely,
                        dfmax = as.integer(p / 2))
    residualsy    <- as.numeric(Y - predict(opt_model, cbind(M_mat, Exposure)))
    beta_fit      <- opt_model$beta[-(1 + p)]
    beta_2        <- opt_model$beta[(1 + p)]
    selection_set <- which(beta_fit != 0)
  } else {
    y.data    <- data.frame(Y = Y, Exposure = Exposure, M_mat)
    lm.fit    <- lm(Y ~ ., data = y.data)
    residualsy <- lm.fit$residuals
    beta_fit  <- summary(lm.fit)$coefficients[-1, 1]
    beta_2    <- beta_fit["Exposure"]
    selection_set <- seq_len(p)
  }
  
  if (length(selection_set) == 0) selection_set <- integer(0)
  
  # d0CRT p-values for taxa in selection_set
  pvl_lst <- rep(1, p)
  nde.p   <- 1
  
  if (p == 1) {
    selection_set <- 1
    if (modelm == "ZINB") Cond_M <- Creat_condition_zinb(M_mat, 1, Exposure, n.times)
    if (modelm == "ZIP")  Cond_M <- Creat_condition_zip(M_mat, 1, Exposure, n.times)
    if (modelm == "NB")   Cond_M <- Creat_condition_nb(M_mat, 1, Exposure, n.times)
    
    M_res_ob <- Cond_M$res_m
    data0    <- data.frame(Y = Y, Exposure = Exposure)
    
    if (modely == "binomial") {
      eps_res <- Y - 1 / (1 + exp(-predict(glm(Y ~ Exposure,
                                               family = modely, data = data0))))
    } else {
      eps_res <- Y - predict(lm(Y ~ Exposure, data = data0))
    }
    
    imp_obe <- abs(mean(M_res_ob * eps_res)) / mean(M_res_ob^2)
    
    library(plyr)
    list.s   <- unlist(lapply(Cond_M$res_m_samp, length))
    list.sel <- which(list.s == 1)
    M_res_sample <- if (length(list.sel) != 0) {
      as.matrix(do.call(rbind.fill, Cond_M$res_m_samp[-list.sel]))
    } else {
      as.matrix(do.call(rbind.fill, Cond_M$res_m_samp))
    }
    
    var_lst_sample <- apply(M_res_sample, 1, function(v) mean(unlist(v)^2))
    t_lst <- abs(M_res_sample %*% eps_res / n) / var_lst_sample
    pvl_lst[1] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
    
  } else {
    # NDE p-value via exposure permutation
    Cond_E <- list(res_e_samp = list())
    for (j in 1:n.times) {
      set.seed(j)
      indx.e <- sample(1:n)
      Cond_E$res_e_samp[[j]] <- as.data.frame(t(Exposure[indx.e]))
    }
    
    set.seed(123)
    cv_lasso_null <- cv.glmnet(cbind(M_mat), Y, alpha = 1,
                               family = modely, dfmax = as.integer(p / 2))
    lamb_null      <- cv_lasso_null$lambda.min
    model_res_null <- glmnet(cbind(M_mat), Y, alpha = 1,
                             lambda = lamb_null, family = modely,
                             dfmax = as.integer(p / 2))
    
    if (modely == "binomial") {
      eps_res <- Y - 1 / (1 + exp(-predict(model_res_null, cbind(M_mat))))
    } else {
      eps_res <- as.numeric(Y - predict(model_res_null, cbind(M_mat)))
    }
    
    imp_exp <- abs(mean(Exposure * eps_res)) / mean(Exposure^2)
    
    library(plyr)
    list.s   <- unlist(lapply(Cond_E$res_e_samp, length))
    list.sel <- which(list.s == 1)
    E_res_sample <- if (length(list.sel) != 0) {
      as.matrix(do.call(rbind.fill, Cond_E$res_e_samp[-list.sel]))
    } else {
      as.matrix(do.call(rbind.fill, Cond_E$res_e_samp))
    }
    
    var_lst_sample <- apply(E_res_sample, 1, function(v) mean(unlist(v)^2))
    t_lst <- abs(E_res_sample %*% eps_res / n) / var_lst_sample
    nde.p <- mean(c(1, ifelse(t_lst >= imp_exp, 1, 0)))
    
    # Taxa-level d0CRT p-values (only for selection_set)
    if (length(selection_set) > 0) {
      for (j in seq_along(selection_set)) {
        indx <- selection_set[j]
        
        if (modelm == "ZINB") Cond_M <- Creat_condition_zinb(M_mat, indx, Exposure, n.times)
        if (modelm == "ZIP")  Cond_M <- Creat_condition_zip(M_mat, indx, Exposure, n.times)
        if (modelm == "NB")   Cond_M <- Creat_condition_nb(M_mat, indx, Exposure, n.times)
        
        M_res_ob <- Cond_M$res_m
        
        set.seed(123)
        cv_lasso_null2 <- cv.glmnet(
          cbind(M_mat[, -indx, drop = FALSE], Exposure), Y,
          alpha = 1, family = modely, dfmax = as.integer(p / 2)
        )
        lamb_null2      <- cv_lasso_null2$lambda.min
        model_res_null2 <- glmnet(
          cbind(M_mat[, -indx, drop = FALSE], Exposure), Y,
          alpha = 1, lambda = lamb_null2, family = modely,
          dfmax = as.integer(p / 2)
        )
        
        if (modely == "binomial") {
          eps_res2 <- Y - 1 / (1 + exp(-predict(
            model_res_null2,
            cbind(M_mat[, -indx, drop = FALSE], Exposure)
          )))
        } else {
          eps_res2 <- as.numeric(Y - predict(
            model_res_null2,
            cbind(M_mat[, -indx, drop = FALSE], Exposure)
          ))
        }
        
        imp_obe <- abs(mean(M_res_ob * eps_res2)) / mean(M_res_ob^2)
        
        list.s   <- unlist(lapply(Cond_M$res_m_samp, length))
        list.sel <- which(list.s == 1)
        M_res_sample <- if (length(list.sel) != 0) {
          as.matrix(do.call(rbind.fill, Cond_M$res_m_samp[-list.sel]))
        } else {
          as.matrix(do.call(rbind.fill, Cond_M$res_m_samp))
        }
        
        var_lst_sample <- apply(M_res_sample, 1, function(v) mean(unlist(v)^2))
        t_lst <- abs(M_res_sample %*% eps_res2 / n) / var_lst_sample
        pvl_lst[indx] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
      }
    }
  }
  
  # Offset (log library size)
  tij_mat <- log(rowSums(M_mat))
  
  # --- Model-specific mediator testing (ZINB / ZIP / NB) ---
  # Only process selection_set; initialize output vectors
  aic <- bic <- final.p <- niea <- niep <- nie <- numeric(0)
  alpha.p <- gamma.p <- numeric(0)
  residualsm <- list()
  alpha_mat <- gamma_mat <- NULL
  
  if (length(selection_set) > 0) {
    ns <- length(selection_set)
    aic <- bic <- final.p <- niea <- niep <- nie <- rep(NA_real_, ns)
    alpha.p <- gamma.p <- rep(NA_real_, ns)
    residualsm <- vector("list", ns)
    alpha_mat  <- matrix(NA_real_, ns, 2)
    gamma_mat  <- matrix(NA_real_, ns, 2)
  }
  
  mediation_set <- integer(0)
  nie_keep <- niea_keep <- niep_keep <- numeric(0)
  niepval_keep <- nieap_keep <- niepp_keep <- numeric(0)
  
  # Note: Only ZINB branch shown in detail below; ZIP and NB follow the same
  
  # pattern with minor model differences. See original code for full branches.
  
  if (modelm == "ZINB" && length(selection_set) > 0) {
    for (j in seq_along(selection_set)) {
      Y_data <- data.frame(
        Exposure = Exposure,
        Mediator = as.vector(M_mat[, selection_set[j]])
      )
      
      if (sum(Y_data$Mediator == 0) == 0) {
        nb.fit <- glm.nb(Mediator ~ Exposure, data = Y_data)
        residualsm[[j]] <- nb.fit$residuals
        final.p[j] <- summary(nb.fit)$coefficients[2, 4]
        alpha_mat[j, ] <- summary(nb.fit)$coefficients[1:2, 1]
      } else {
        zinb.fit <- if (sum(is.infinite(tij_mat)) != 0) {
          try(zeroinfl(Mediator ~ Exposure | Exposure,
                       data = Y_data, dist = "negbin", link = "logit"),
              silent = TRUE)
        } else {
          try(zeroinfl(Mediator ~ Exposure | Exposure,
                       offset = tij_mat, data = Y_data,
                       dist = "negbin", link = "logit"),
              silent = TRUE)
        }
        
        if (inherits(zinb.fit, "try-error")) {
          final.p[j] <- NA
          aic[j] <- bic[j] <- NA
          residualsm[[j]] <- NA
        } else {
          aic[j] <- AIC(zinb.fit)
          bic[j] <- BIC(zinb.fit)
          residualsm[[j]] <- zinb.fit$residuals
          alpha_mat[j, ] <- summary(zinb.fit)$coefficients$count[1:2, 1]
          gamma_mat[j, ] <- summary(zinb.fit)$coefficients$zero[1:2, 1]
          alpha.p[j] <- summary(zinb.fit)$coefficients$count[2, 4]
          gamma.p[j] <- summary(zinb.fit)$coefficients$zero[2, 4]
          
          cov.mat <- matrix(c(
            vcov(zinb.fit)[2, 2], vcov(zinb.fit)[2, 4],
            vcov(zinb.fit)[4, 2], vcov(zinb.fit)[4, 4]
          ), 2, 2)
          
          wald.t <- try(
            t(matrix(c(alpha_mat[j, 2], gamma_mat[j, 2]))) %*%
              ginv(cov.mat) %*%
              matrix(c(alpha_mat[j, 2], gamma_mat[j, 2])),
            silent = TRUE
          )
          final.p[j] <- if (inherits(wald.t, "try-error")) NA
          else (1 - pchisq(as.numeric(wald.t), df = 2))
        }
      }
      
      # NIE computation
      beta.val  <- beta_fit[selection_set[j]]
      alpha.val <- alpha_mat[j, ]
      gamma.val <- gamma_mat[j, ]
      
      ones_mat <- as.matrix(data.frame(rep(1, n)))
      one_one  <- as.matrix(data.frame(1, rep(1, n)))
      
      niea[j] <- mean(beta.val *
                        (1 / (1 + exp(ones_mat %*% matrix(gamma.val[1])))) *
                        (exp(one_one %*% matrix(alpha.val[1:2]) + tij_mat) -
                           exp(ones_mat %*% matrix(alpha.val[1]) + tij_mat)))
      
      niep[j] <- mean(beta.val *
                        (exp(one_one %*% matrix(alpha.val[1:2]) + tij_mat)) *
                        ((1 / (1 + exp(one_one %*% matrix(gamma.val[1:2])))) -
                           (1 / (1 + exp(ones_mat %*% matrix(gamma.val[1]))))))
      
      nie[j] <- niea[j] + niep[j]
    }
  }
  
  # Selection based on joint p within selection_set
  index.p <- if (length(selection_set) > 0) pvl_lst[selection_set] else numeric(0)
  bind.p  <- rbind(p.adjust(index.p, method), p.adjust(final.p, method))
  joint.p <- if (length(selection_set) > 0) apply(bind.p, 2, max) else numeric(0)
  
  index.mi      <- which(joint.p <= FDR)
  mediation_set <- if (length(index.mi) > 0) selection_set[index.mi] else integer(0)
  
  nie_keep     <- if (length(index.mi) > 0) nie[index.mi] else numeric(0)
  niea_keep    <- if (length(index.mi) > 0) niea[index.mi] else numeric(0)
  niep_keep    <- if (length(index.mi) > 0) niep[index.mi] else numeric(0)
  niepval_keep <- if (length(index.mi) > 0) joint.p[index.mi] else numeric(0)
  
  p_table <- build_p_table(M_mat, selection_set, pvl_lst, final.p,
                           method, mediation_set)
  
  list(
    Mediators  = mediation_set,
    NDE        = beta_2,
    NIE        = nie_keep,
    NIEA       = niea_keep,
    NIEP       = niep_keep,
    NDE.pval   = nde.p,
    NIE.pval   = niepval_keep,
    AIC = aic, BIC = bic,
    residualsy = residualsy,
    residualsm = residualsm,
    p_table    = p_table
  )
}

CRAmed_sim <- function(count1, treat1, y1) {
  t0 <- proc.time()
  
  results.infant <- CRAmed_fullp(
    M_mat    = count1,
    Y        = y1,
    Exposure = as.matrix(treat1),
    n.times  = 100,
    n.perm   = 100,
    CI       = FALSE
  )
  
  p_joint <- setNames(results.infant$p_table$p_joint,
                      results.infant$p_table$taxa)
  
  t1      <- proc.time() - t0
  runtime <- unname(t1["elapsed"])
  
  list(
    discoveries = results.infant$Mediators,
    p_med       = p_joint,
    runtime_sec = runtime
  )
}


# =============================================================================
# MarZIC
# =============================================================================

Mar_sim <- function(count1, treat1, y1, alpha = 0.05,
                    adjust_method = "fdr", num_cores = 1) {
  t0 <- proc.time()
  stopifnot(nrow(count1) == length(treat1), nrow(count1) == length(y1))
  
  MicrobData <- as.matrix(count1)
  mode(MicrobData) <- "numeric"
  orig_taxa <- colnames(MicrobData)
  if (is.null(orig_taxa))
    orig_taxa <- paste0("Taxon", seq_len(ncol(MicrobData)))
  colnames(MicrobData) <- orig_taxa
  MicrobData[!is.finite(MicrobData)] <- 0
  
  # Library size as covariate
  libsize <- pmax(1, rowSums(MicrobData))
  CovData <- data.frame(Y = as.numeric(y1), X = as.numeric(treat1),
                        libsize = libsize)
  
  # Remove invalid rows
  keep_row    <- is.finite(CovData$Y) & is.finite(CovData$X) &
    (rowSums(MicrobData) > 0)
  MicrobData2 <- MicrobData[keep_row, , drop = FALSE]
  CovData2    <- CovData[keep_row, , drop = FALSE]
  
  # Remove all-zero taxa and normalize
  keep_col    <- colSums(MicrobData2) > 0
  MicrobData2 <- MicrobData2[, keep_col, drop = FALSE]
  MicrobData2 <- t(apply(MicrobData2, 1, function(x) x / sum(x)))
  
  res <- MarZIC::MarZIC(
    MicrobData = MicrobData2, CovData = CovData2,
    lib_name = "libsize", y_name = "Y", x_name = "X",
    x4_inter = FALSE, x5_inter = FALSE, conf_name = NULL,
    taxa_of_interest = "all", num_cores = num_cores,
    adjust_method = adjust_method,
    zero_prop_NIE2 = 0.001, zero_count_NIE2 = 1,
    taxDropThresh = 0.999, taxDropCount = 1, SDThresh = 0.001
  )
  
  # Map adjusted p-values back to full taxa
  kept_taxa   <- colnames(MicrobData2)
  p_med       <- res$NIE_save$`p value adj`
  p_full      <- rep(1, length(orig_taxa))
  names(p_full) <- orig_taxa
  idx_match   <- match(kept_taxa, orig_taxa)
  p_med_clean <- as.numeric(p_med)
  p_med_clean[!is.finite(p_med_clean) | is.na(p_med_clean)] <- 1
  p_full[idx_match] <- p_med_clean
  
  sig_idx     <- which(p_full < alpha)
  runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
  
  list(discoveries = sig_idx, p_med = p_full, runtime_sec = runtime_sec)
}


# =============================================================================
# microHIMA
# =============================================================================

hima_discrete_q <- function(p_vec, taxa_names = NULL,
                            alpha_grid = c(0.01, 0.02, 0.05, 0.10, 0.20),
                            eps = 1e-4, simes = FALSE) {
  p_vec <- as.numeric(p_vec)
  m     <- length(p_vec)
  if (is.null(taxa_names)) {
    taxa_names <- names(p_vec)
    if (is.null(taxa_names)) taxa_names <- paste0("taxon_", seq_len(m))
  }
  
  q_disc <- rep(1, m)
  bad    <- !is.finite(p_vec) | p_vec < 0 | p_vec > 1
  q_disc[bad] <- NA_real_
  
  ok <- which(!bad)
  if (length(ok) == 0) {
    return(list(
      table = data.frame(Index = taxa_names, p = p_vec, q_value = q_disc),
      selected = setNames(vector("list", length(alpha_grid)),
                          as.character(alpha_grid))
    ))
  }
  
  p_ok        <- p_vec[ok]
  hom         <- hommel::hommel(p_ok, simes = simes)
  alpha_grid  <- sort(unique(alpha_grid))
  selected_list <- setNames(vector("list", length(alpha_grid)),
                            as.character(alpha_grid))
  
  for (a in alpha_grid) {
    sel <- which(p_ok < a)
    if (length(sel) == 0) {
      selected_list[[as.character(a)]] <- character(0)
      next
    }
    N1 <- hommel::discoveries(hom, sel, incremental = TRUE, alpha = a)
    L  <- length(sel)
    N2 <- numeric(L)
    if (L >= 2) N2[2:L] <- N1[1:(L - 1)]
    N0 <- N1 - N2
    ID_local  <- sel[which(N0 > 0)]
    ID_global <- ok[ID_local]
    
    selected_list[[as.character(a)]] <- taxa_names[ID_global]
    newly <- ID_global[is.finite(q_disc[ID_global]) & q_disc[ID_global] >= 1]
    if (length(newly) > 0) q_disc[newly] <- a - eps
  }
  
  list(
    table    = data.frame(Index = taxa_names, p = p_vec, q_value = q_disc),
    selected = selected_list
  )
}


HIMA_micro_sim1 <- function(count1, X, Y, COV = NULL,
                            verbose = TRUE, parallel = FALSE, ncore = 1) {
  pseudo   <- 0.5
  OTU_comp <- sweep(count1 + pseudo, 1, rowSums(count1 + pseudo), "/")
  
  HIMA_recover <- function(X, OTU, Y, COV = NULL, FDRcut = 0.05,
                           verbose = FALSE, parallel = FALSE, ncore = 1) {
    X     <- matrix(X, ncol = 1)
    M_raw <- as.matrix(OTU)
    M_ID_name <- colnames(M_raw)
    if (is.null(M_ID_name)) M_ID_name <- seq_len(ncol(M_raw))
    if (!is.null(COV)) { COV <- as.matrix(COV); X <- cbind(X, COV) }
    X <- scale(X)
    Y <- Y - mean(Y)
    M <- M_raw
    n <- nrow(M)
    d <- ncol(M)
    M1 <- t(t(M_raw[, 1]))
    
    if (verbose) message("Step 1: ILR Transformation + De-biased Lasso ...")
    HIMA:::checkParallel("hima_microbiome", parallel, ncore, verbose)
    
    library(foreach)
    results_loop <- foreach(k = seq_len(d), .combine = rbind) %dopar% {
      M <- M_raw; M[, 1] <- M[, k]; M[, k] <- M1
      MT <- matrix(0, n, d - 1)
      for (i in 1:n) {
        for (j in 1:(d - 1)) {
          C_1 <- sqrt((d - j) / (d - j + 1))
          C_2 <- prod(M[i, (j + 1):d]^(1 / (d - j)))
          MT[i, j] <- C_1 * log(M[i, j] / C_2)
        }
      }
      MT <- scale(MT)
      MX <- cbind(MT, X)
      fit.dlasso <- HIMA:::DLASSO_fun(MX, Y)
      beta_est <- fit.dlasso[1]; beta_se <- fit.dlasso[2]
      P_b <- 2 * (1 - pnorm(abs(beta_est / beta_se), 0, 1))
      lm.out    <- summary(stats::lm(MT[, 1] ~ X))
      alpha_est <- lm.out$coefficients[2, 1]
      alpha_se  <- lm.out$coefficients[2, 2]
      P_a <- 2 * (1 - pnorm(abs(alpha_est / alpha_se), 0, 1))
      c(beta_est, beta_se, alpha_est, alpha_se, max(P_a, P_b), P_a, P_b)
    }
    
    if (is.null(dim(results_loop)))
      results_loop <- matrix(results_loop, nrow = 1)
    P_raw_DLASSO <- results_loop[, 5]
    
    pmax_res <- hima_discrete_q(P_raw_DLASSO)
    list(q_value = pmax_res$table$q_value, select_list = pmax_res$selected)
  }
  
  t0  <- proc.time()
  fit <- HIMA_recover(X = X, OTU = OTU_comp, Y = Y, COV = COV,
                      verbose = verbose, parallel = parallel, ncore = ncore)
  runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
  
  list(index = NULL, q_value = fit$q_value, runtime_sec = runtime_sec)
}


# =============================================================================
# multimedia
# =============================================================================

multimedia_sim <- function(count1, treat1, y1, q_value = 0.05,
                           add_lib_size = TRUE, pseudo = 0.5,
                           alpha_out = 1) {
  t0 <- proc.time()
  stopifnot(length(treat1) == nrow(count1), length(y1) == nrow(count1))
  
  orig_names_full <- colnames(count1)
  if (is.null(orig_names_full))
    orig_names_full <- paste0("V", seq_len(ncol(count1)))
  
  X0 <- as.matrix(count1)
  storage.mode(X0) <- "numeric"
  
  # Drop all-zero samples
  rs <- rowSums(X0)
  keep_samp <- rs > 0
  if (!all(keep_samp)) {
    X0     <- X0[keep_samp, , drop = FALSE]
    treat1 <- treat1[keep_samp]
    y1     <- y1[keep_samp]
    rs     <- rs[keep_samp]
  }
  
  # Drop all-zero / zero-variance taxa
  keep_taxa <- (colSums(X0) > 0) & (apply(X0, 2, var) > 0)
  if (sum(keep_taxa) == 0) {
    runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
    p_full <- rep(1, length(orig_names_full))
    names(p_full) <- orig_names_full
    return(list(discoveries_default = integer(0), q_value_disc = p_full,
                runtime_sec = runtime_sec, discoveries_by_q = list()))
  }
  
  X_use          <- X0[, keep_taxa, drop = FALSE]
  orig_names_use <- orig_names_full[keep_taxa]
  
  # CLR transformation
  X_use <- X_use + pseudo
  prop  <- sweep(X_use, 1, rowSums(X_use), "/")
  logp  <- log(prop)
  clrM  <- logp - rowMeans(logp)
  
  safe_names <- make.names(orig_names_use, unique = TRUE)
  colnames(clrM) <- safe_names
  mediator <- safe_names
  
  # Build data frame
  treat01 <- as.numeric(treat1)
  stopifnot(all(treat01 %in% c(0, 1)), length(unique(treat01)) >= 2)
  treat12 <- treat01 + 1  # 0/1 -> 1/2
  
  df <- data.frame(treatment = treat12, outcome = as.numeric(y1),
                   check.names = FALSE)
  if (add_lib_size) df$lib_size <- as.numeric(log(rs + 1))
  df <- cbind(df, as.data.frame(clrM, check.names = FALSE))
  
  # Fit multimedia model
  exper <- multimedia::mediation_data(
    df, outcomes = "outcome", treatments = "treatment",
    mediator = mediator,
    pretreatments = if (add_lib_size) "lib_size" else NULL
  )
  
  med_est <- multimedia::lm_model()
  out_est <- if (!add_lib_size) multimedia::lm_model()
  else multimedia::glmnet_model(alpha = alpha_out)
  
  fit <- multimedia::multimedia(
    exper, outcome_estimator = out_est, mediation_estimator = med_est
  ) |> multimedia::estimate(exper)
  
  # Null contrast + FDR at multiple levels
  contrast <- multimedia::null_contrast(
    fit, exper, nullification = "T->M", f = multimedia::indirect_pathwise
  )
  
  q_grid <- c(0.01, 0.02, 0.05, 0.10, 0.20)
  eps    <- 1e-4
  
  fdr_list <- lapply(q_grid, function(q) {
    fdr <- multimedia::fdr_summary(contrast, effect = "indirect_pathwise",
                                   q_value = q)
    fdr[fdr$source == "real", , drop = FALSE]
  })
  
  # Discoveries at each q level
  discoveries_by_q <- setNames(vector("list", length(q_grid)),
                               paste0("q=", q_grid))
  for (i in seq_along(q_grid)) {
    sel_mediator <- fdr_list[[i]]$mediator[fdr_list[[i]]$keep]
    idx_safe <- match(sel_mediator, safe_names)
    idx_safe <- idx_safe[!is.na(idx_safe)]
    discoveries_by_q[[i]] <- which(keep_taxa)[idx_safe]
  }
  
  # Construct discrete q-values
  q_disc_safe <- rep(1, length(safe_names))
  names(q_disc_safe) <- safe_names
  for (i in seq_along(q_grid)) {
    q <- q_grid[i]
    sel_mediator <- fdr_list[[i]]$mediator[fdr_list[[i]]$keep]
    idx_safe <- match(sel_mediator, safe_names)
    idx_safe <- idx_safe[!is.na(idx_safe)]
    to_set   <- idx_safe[q_disc_safe[idx_safe] == 1]
    if (length(to_set)) q_disc_safe[to_set] <- max(q - eps, 0)
  }
  
  # Map back to full taxa
  q_disc_full <- rep(1, length(orig_names_full))
  names(q_disc_full) <- orig_names_full
  q_disc_full[keep_taxa] <- as.numeric(q_disc_safe[safe_names])
  q_disc_full[is.na(q_disc_full)] <- 1
  
  discoveries_default <- which(q_disc_full <= q_value)
  runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
  
  list(
    discoveries_default = discoveries_default,
    q_value_disc        = q_disc_full,
    runtime_sec         = runtime_sec,
    discoveries_by_q    = discoveries_by_q
  )
}

###############################################################################
# ========================= 1. Run all methods ================================
###############################################################################

# ── 1.1 CAMRA ────────────────────────────────────────────────────────────────

res1_prop <- CAMRA(
  outcome   = y1,
  treatment = treat1,
  mediators = count1,
  screen    = FALSE,
  CClasso   = FALSE,
  seed      = 123
)

# Extract results
idx_PRO  <- res1_prop$idx_detected
otu_PRO  <- colnames(count1)[idx_PRO]
PRO_beta <- rbind(res1_prop$beta_l, res1_prop$beta_r)
IDE_PRO  <- res1_prop$beta_l * res1_prop$beta_r
colnames(PRO_beta) <- colnames(count1)
rownames(PRO_beta) <- c("treat", "outcome")

cat("\nCAMRA detected", length(idx_PRO), "mediators\n")
if (length(idx_PRO) > 0) print(PRO_beta[, idx_PRO])

# One-sided p-values (for volcano plots)
res_l <- recover_l_PALM(count1, treat1)
selected_values1 <- res_l$p

res_r <- recover_r(count_matrix = count1, treat_cov = treat1, y = y1)
selected_values2 <- res_r$p


# ── 1.2 LDM-med ─────────────────────────────────────────────────────────────

res_ldm <- ldm_med_sim(count1, treat1, y1, fdr_nominal = 0.05, seed = 42)

p_LDM    <- res_ldm$p_otu
LDM_beta <- res_ldm$beta
otu_LDM  <- res_ldm$otu_detected
idx_LDM  <- res_ldm$discoveries

cat("LDM-med detected", length(idx_LDM), "mediators\n")


# ── 1.3 CRAmed ──────────────────────────────────────────────────────────────

res_CRA <- CRAmed_sim(count1, treat1, y1)
idx_CRA <- res_CRA$discoveries
otu_CRA <- colnames(count1)[idx_CRA]

cat("CRAmed detected", length(idx_CRA), "mediators\n")


# ── 1.4 microHIMA ───────────────────────────────────────────────────────────

res_hima1 <- HIMA_micro_sim1(count1, treat1, y1)
otu_hima1 <- colnames(count1)[which(res_hima1$q_value < 0.05)]

cat("microHIMA detected", length(otu_hima1), "mediators\n")


# ── 1.5 multimedia ──────────────────────────────────────────────────────────

res_multimedia <- multimedia_sim(count1, treat1, y1)
otu_multimedia <- colnames(count1)[which(res_multimedia$q_value_disc < 0.05)]

cat("multimedia detected", length(otu_multimedia), "mediators\n")


# ── 1.6 MarZIC ──────────────────────────────────────────────────────────────

MicrobData <- count1
CovData <- data.frame(
  Y       = y1,
  X       = treat1,
  libsize = pmax(1, rowSums(MicrobData))
)

res_mar <- MarZIC::MarZIC(
  MicrobData = MicrobData,
  CovData    = CovData,
  lib_name   = "libsize",
  y_name     = "Y",
  x_name     = "X",
  x4_inter   = FALSE,
  x5_inter   = FALSE,
  conf_name  = NULL,
  taxa_of_interest = "all",
  num_cores       = 5,
  adjust_method   = "fdr",
  zero_prop_NIE2  = 0.001,
  zero_count_NIE2 = 1,
  taxDropThresh   = 0.999,
  taxDropCount    = 1,
  SDThresh        = 0.001
)

pvals   <- res_mar$NIE_save$`p value adj`
otu_Mar <- colnames(count1)[!is.na(pvals) & !is.nan(pvals) & pvals < 0.05]
idx_Mar <- match(otu_Mar, colnames(count1))

cat("MarZIC detected", length(otu_Mar), "mediators\n")


# ── 1.7 Global mediation tests ──────────────────────────────────────────────

library(matrixStats)

permanovaFL_res <- permanovaFL_sim(count1, treat1, y1)
MODIMA_res      <- MODIMA_sim(count1, treat1, y1)
Medtest_res     <- Medtest_sim(count1, treat1, y1)
cmm_res         <- ccmm_sim(count1, treat1, y1)

cat("\n--- Global mediation p-values ---\n")
cat("CCM:          ", cmm_res$global_p, "\n")
cat("LDM-med:      ", res_ldm$global_p, "\n")
cat("PERMANOVA-FL: ", permanovaFL_res$global_p, "\n")
cat("MODIMA:       ", MODIMA_res$global_p, "\n")
cat("MedTest:      ", Medtest_res$global_p, "\n")


###############################################################################
# ========================= 2. Shared plot settings ===========================
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(forcats)
  library(scales)
  library(ComplexUpset)
  library(compositions)
})

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

base_fs   <- 22
theme_big <- theme_bw(base_size = base_fs) +
  theme(
    axis.title   = element_text(size = rel(1.05)),
    axis.text    = element_text(size = rel(0.95)),
    legend.title = element_text(size = rel(1.00)),
    legend.text  = element_text(size = rel(0.95)),
    strip.text   = element_text(size = rel(0.95)),
    plot.title   = element_text(size = rel(1.10))
  )

treat_lab   <- c("CHN", "USA")
cut_median  <- 25
legend_name <- "BMI group"

# Helper: format OTU labels for display
make_otu_label <- function(otu) {
  otu <- as.character(otu)
  has_bar <- grepl("\\|", otu)
  left    <- sub("\\|.*$", "", otu)
  right   <- sub("^.*\\|", "", otu)
  ifelse(has_bar & left != right, paste0(left, " | ", right), otu)
}


###############################################################################
# ========================= 3.1 Panel C: Treatment -> OTU =====================
###############################################################################

{
  otu_tab <- count1
  otu_sel <- intersect(otu_PRO, colnames(otu_tab))
  
  if (!is.null(names(treat1)) && all(rownames(otu_tab) %in% names(treat1))) {
    treat_use <- treat1[rownames(otu_tab)]
  } else {
    stopifnot(length(treat1) == nrow(otu_tab))
    treat_use <- treat1
  }
  
  meta2 <- data.frame(
    sample = rownames(otu_tab),
    treat  = treat_use,
    stringsAsFactors = FALSE
  )
  if (is.numeric(meta2$treat) || all(meta2$treat %in% c(0, 1))) {
    meta2$treat <- factor(meta2$treat, levels = c(0, 1), labels = treat_lab)
  } else {
    meta2$treat <- factor(meta2$treat)
  }
  
  lv <- levels(meta2$treat)
  
  cols_fill  <- c("#E69F00", "#56B4E9")
  cols_color <- c("#B05D00", "#1F78B4")
  pal_fill   <- setNames(cols_fill[seq_along(lv)], lv)
  pal_color  <- setNames(cols_color[seq_along(lv)], lv)
  
  # ---- Prevalence (non-zero proportion) by treatment group ----
  prev_df <- (otu_tab[, idx_PRO, drop = FALSE] > 0) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "otu", values_to = "detected") %>%
    left_join(meta2, by = "sample") %>%
    group_by(treat, otu) %>%
    summarise(prevalence = mean(detected, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    mutate(otu_label = make_otu_label(otu))
  
  prev_df2 <- prev_df %>%
    group_by(otu_label) %>%
    mutate(.ord = mean(prevalence, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(otu_label = fct_reorder(otu_label, .ord),
           panel = "")
  
  p_prev <- ggplot(prev_df2, aes(x = treat, y = prevalence, fill = treat)) +
    geom_col(width = 0.75) +
    scale_y_continuous(limits = c(0, 0.25),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = pal_fill, guide = "none") +
    labs(x = NULL, y = "Prevalence (Non-zero proportion)", fill = NULL) +
    theme_bw(base_size = 12) + theme_big +
    theme(
      strip.position   = "top",
      strip.placement  = "outside",
      strip.background = element_rect(fill = "grey95", color = "grey40"),
      legend.position  = "none",
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  
  ggsave("figure/treat_prev.png", plot = p_prev,
         width = 8, height = 7, dpi = 300, device = ragg::agg_png)
  
  # ---- Abundance (log10 RA) by treatment group ----
  pseudo    <- 0.5
  abund_mat <- log10(sweep(count1 + pseudo, 1, rowSums(count1 + pseudo), "/"))
  otu_clr   <- abund_mat
  
  clr_long <- otu_clr[, idx_PRO, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "otu", values_to = "clr_abund")
  
  cnt_long <- count1[, idx_PRO, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "otu", values_to = "otu_count")
  
  ab_df <- clr_long %>%
    left_join(cnt_long, by = c("sample", "otu")) %>%
    filter(!is.na(otu_count) & otu_count > 0) %>%
    left_join(meta2, by = "sample") %>%
    mutate(otu_label = make_otu_label(otu))
  
  ab_df <- ab_df %>%
    mutate(otu_label = factor(otu_label,
                              levels = c("Acidaminococcus intestini")))
  
  p_abund2 <- ggplot(ab_df, aes(x = treat, y = clr_abund, fill = treat)) +
    geom_violin(trim = TRUE, alpha = 0.45, color = "grey25") +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.25,
                 color = "grey25") +
    geom_jitter(width = 0.08, size = 0.45, alpha = 0.25, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6,
                 fill = "white", color = "black", stroke = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                 color = "black", linewidth = 0.35) +
    scale_fill_manual(values = pal_fill) +
    labs(x = "Treat", y = expression(log[10]("RA")), fill = "Treat") +
    theme_bw(base_size = 12) + theme_big +
    theme(
      legend.position  = "none",
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      strip.text       = element_blank(),
      strip.background = element_blank()
    )
  
  ggsave("figure/treat_compare.png", plot = p_abund2,
         width = 8, height = 7, dpi = 300, device = ragg::agg_png)
}


###############################################################################
# ========================= 3.2 Panel D: OTU -> Outcome ======================
###############################################################################

{
  otu_tab <- count1[, idx_PRO, drop = FALSE]
  otu_sel <- colnames(otu_tab)
  
  if (!is.null(names(y1)) && all(rownames(otu_tab) %in% names(y1))) {
    y_use <- y1[rownames(otu_tab)]
  } else {
    stopifnot(length(y1) == nrow(otu_tab))
    y_use <- y1
  }
  y_use <- as.numeric(y_use)
  
  bw_group <- ifelse(y_use <= cut_median,
                     paste0("\u2264 ", cut_median),
                     paste0("> ", cut_median))
  
  meta_bw <- data.frame(
    sample   = rownames(otu_tab),
    y        = y_use,
    bw_group = factor(bw_group,
                      levels = c(paste0("\u2264 ", cut_median),
                                 paste0("> ", cut_median))),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(bw_group))
  
  meta_bw <- meta_bw %>%
    mutate(bw_group_n = {
      tt <- table(bw_group)
      factor(bw_group,
             levels = levels(bw_group),
             labels = paste0(levels(bw_group), " (n=",
                             as.integer(tt[levels(bw_group)]), ")"))
    })
  
  otu_tab2 <- otu_tab[meta_bw$sample, otu_sel, drop = FALSE]
  lv_bw    <- levels(meta_bw$bw_group_n)
  cols_bw  <- c("#008837", "#7B3294")
  pal_bw   <- setNames(cols_bw[seq_along(lv_bw)], lv_bw)
  
  # ---- Prevalence by BMI group ----
  prev_df <- (otu_tab2 > 0) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "otu", values_to = "detected") %>%
    left_join(meta_bw, by = "sample") %>%
    group_by(bw_group_n, otu) %>%
    summarise(prevalence = mean(detected), n = n(), .groups = "drop") %>%
    mutate(otu_label = make_otu_label(otu))
  
  prev_df_bw2 <- prev_df %>%
    group_by(otu_label) %>%
    mutate(.ord = mean(prevalence, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(otu_label = fct_reorder(otu_label, .ord),
           panel = "")
  
  p_prev_bw <- ggplot(prev_df_bw2,
                      aes(x = bw_group_n, y = prevalence,
                          fill = bw_group_n)) +
    geom_col(width = 0.75) +
    scale_y_continuous(limits = c(0, 0.25),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = pal_bw, name = legend_name) +
    labs(x = NULL, y = "Prevalence (Non-zero proportion)") +
    theme_bw(base_size = 12) + theme_big +
    theme(
      strip.position   = "top",
      strip.placement  = "outside",
      strip.background = element_rect(fill = "grey95", color = "grey40"),
      legend.position  = "none",
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  
  ggsave("figure/bw_median_prev.png", plot = p_prev_bw,
         width = 8, height = 7, dpi = 300, device = ragg::agg_png)
  
  # ---- Abundance by BMI group (log10 RA) ----
  pseudo     <- 0.5
  count_all2 <- count1[meta_bw$sample, , drop = FALSE]
  abund_all  <- log10(sweep(count_all2 + pseudo, 1,
                            rowSums(count_all2 + pseudo), "/"))
  abund_mat  <- abund_all[, idx_PRO, drop = FALSE]
  value_nm   <- "clr_abund"
  ylab       <- expression(log[10]("RA"))
  
  count_long <- otu_tab2 %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "otu", values_to = "otu_count")
  
  ab_df <- abund_mat %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "otu", values_to = value_nm) %>%
    left_join(meta_bw, by = "sample") %>%
    left_join(count_long, by = c("sample", "otu")) %>%
    filter(!is.na(otu_count) & otu_count > 0) %>%
    mutate(otu_label = make_otu_label(otu))
  
  ab_df$bw_group_n <- factor(ab_df$bw_group_n, levels = lv_bw)
  
  p_abund_bw <- ggplot(ab_df,
                       aes(x = bw_group_n, y = .data[[value_nm]],
                           fill = bw_group_n)) +
    geom_violin(trim = TRUE, alpha = 0.45, color = "grey25") +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.25,
                 color = "grey25") +
    geom_jitter(width = 0.08, size = 0.45, alpha = 0.25, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6,
                 fill = "white", color = "black", stroke = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                 color = "black", linewidth = 0.35) +
    scale_fill_manual(values = pal_bw) +
    labs(y = ylab, fill = "Group") +
    theme_bw(base_size = 12) + theme_big +
    theme(
      legend.position  = "none",
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      strip.text       = element_blank(),
      strip.background = element_blank()
    )
  
  # Wilcoxon test for group comparison
  library(rstatix)
  library(ggpubr)
  
  ng <- nlevels(factor(ab_df$bw_group_n))
  if (ng == 2) {
    p_df <- ab_df %>%
      group_by(otu_label) %>%
      wilcox_test(as.formula(paste0(value_nm, " ~ bw_group_n"))) %>%
      ungroup() %>%
      mutate(
        group1 = levels(factor(ab_df$bw_group_n))[1],
        group2 = levels(factor(ab_df$bw_group_n))[2]
      )
  } else {
    p_df <- ab_df %>%
      group_by(otu_label) %>%
      pairwise_wilcox_test(
        as.formula(paste0(value_nm, " ~ bw_group_n")),
        p.adjust.method = "BH"
      ) %>%
      ungroup()
  }
  
  ggsave("figure/median_outcome_compare.png", plot = p_abund_bw,
         width = 8, height = 7, dpi = 300, device = ragg::agg_png)
}


###############################################################################
# ========================= 3.3 Panel A: UpSet plot ===========================
###############################################################################

{
  # Normalize taxa names for cross-method comparison
  normalize_taxa_sets <- function(taxa_sets) {
    normalize_taxa <- function(x) {
      x <- trimws(as.character(x))
      x <- gsub("^(k|p|c|o|f|g|s)__", "", x, ignore.case = TRUE)
      x <- gsub("^(kingdom|phylum|class|order|family|genus|species)[:\\.]",
                "", x, ignore.case = TRUE)
      x <- gsub("[\\.|:|\\||_|;]+", " ", x)
      x <- gsub("\\s+", " ", x)
      x <- tolower(trimws(x))
      x[!duplicated(x)]
    }
    stopifnot(is.list(taxa_sets))
    out <- lapply(taxa_sets, normalize_taxa)
    nm  <- names(out)
    if (!is.null(nm) && anyDuplicated(nm)) {
      out2 <- lapply(unique(nm), function(k)
        unique(unlist(out[nm == k], use.names = FALSE)))
      names(out2) <- unique(nm)
      out <- out2
    }
    out
  }
  
  plot_taxa_upset_two_color <- function(
    taxa_sets,
    highlight_set    = "CAMRA",
    other_color      = "gray",
    highlight_color  = "gray",
    min_intersection = 1,
    max_intersections = 40,
    sort_sets             = c("descending", "ascending", FALSE),
    sort_intersections    = c("descending", "ascending", FALSE),
    sort_intersections_by = c("cardinality", "degree", "ratio"),
    title = NULL,
    quiet_warnings = TRUE
  ) {
    sort_sets             <- match.arg(sort_sets)
    sort_intersections    <- match.arg(sort_intersections)
    sort_intersections_by <- match.arg(sort_intersections_by)
    
    stopifnot(is.list(taxa_sets), length(taxa_sets) >= 2)
    if (is.null(names(taxa_sets)) || any(names(taxa_sets) == ""))
      names(taxa_sets) <- paste0("Set", seq_along(taxa_sets))
    
    # Deduplicate set names
    nm <- names(taxa_sets)
    if (anyDuplicated(nm)) {
      taxa_sets <- lapply(unique(nm), function(k)
        unique(unlist(taxa_sets[nm == k], use.names = FALSE)))
      names(taxa_sets) <- unique(nm)
    }
    
    taxa_sets <- lapply(taxa_sets, function(x) {
      x <- as.character(x)
      unique(x[!is.na(x) & nzchar(x)])
    })
    
    set_names <- names(taxa_sets)
    all_taxa  <- sort(unique(unlist(taxa_sets, use.names = FALSE)))
    stopifnot(length(all_taxa) > 0)
    
    df <- data.frame(taxon = all_taxa, check.names = FALSE)
    for (s in set_names) df[[s]] <- df$taxon %in% taxa_sets[[s]]
    
    set_sizes_plot <- ComplexUpset::upset_set_size(
      mapping = ggplot2::aes(fill = I(other_color))
    ) + ggplot2::labs(y = "Number of selected mediators", x = NULL)
    
    queries <- list()
    if (highlight_set %in% set_names) {
      queries <- list(
        ComplexUpset::upset_query(
          set = highlight_set, fill = highlight_color,
          only_components = "overall_sizes"
        )
      )
    }
    
    make_plot <- function() {
      big_themes <- ComplexUpset::upset_modify_themes(list(
        default = ggplot2::theme(
          text       = element_text(size = 28),
          plot.title = element_text(size = 28),
          axis.title = element_text(size = 20),
          axis.text  = element_text(size = 18),
          strip.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text  = element_text(size = 17)
        ),
        intersections_matrix = ggplot2::theme(
          axis.text.y  = element_text(size = 18),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
        ),
        overall_sizes = ggplot2::theme(
          axis.text  = element_text(size = 18),
          axis.title = element_text(size = 20)
        ),
        `Intersection size` = ggplot2::theme(
          axis.text  = element_text(size = 18),
          axis.title = element_text(size = 20)
        )
      ))
      
      base_ann <- list(
        "Intersection size" =
          ComplexUpset::intersection_size(bar_number_threshold = 1) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
          labs(x = NULL) +
          theme(
            axis.title.x = element_blank(),
            text       = element_text(size = 24),
            axis.title = element_text(size = 24),
            axis.text  = element_text(size = 18)
          )
      )
      
      ComplexUpset::upset(
        df,
        intersect          = set_names,
        min_size           = min_intersection,
        n_intersections    = max_intersections,
        sort_sets          = sort_sets,
        sort_intersections = sort_intersections,
        sort_intersections_by = sort_intersections_by,
        wrap             = TRUE,
        set_sizes        = set_sizes_plot,
        queries          = queries,
        base_annotations = base_ann,
        themes           = big_themes,
        stripes = ComplexUpset::upset_stripes(
          geom   = geom_segment(size = 0),
          colors = c("white", "white")
        )
      )
    }
    
    p <- if (quiet_warnings) suppressWarnings(make_plot()) else make_plot()
    if (!is.null(title)) p <- p + ggtitle(title)
    p
  }
  
  taxa_sets <- list(
    CAMRA      = otu_PRO,
    `LDM-med`  = otu_LDM,
    microHIMA  = otu_hima1,
    CRAmed     = otu_CRA,
    multimedia = otu_multimedia,
    MarZIC     = otu_Mar
  )
  
  taxa_sets_std <- normalize_taxa_sets(taxa_sets)
  
  p_upset <- plot_taxa_upset_two_color(
    taxa_sets_std,
    highlight_set    = "CAMRA",
    min_intersection = 1
  )
  print(p_upset)
  
  ggsave("figure/upset.png", plot = p_upset,
         width = 15, height = 7, dpi = 300, device = ragg::agg_png)
}


###############################################################################
# ========================= 3.4 Panel B: Volcano plots ========================
###############################################################################

{
  alpha   <- 0.05
  use_fdr <- FALSE
  p_floor <- 1e-5
  
  plot_volcano_idxPRO_red <- function(p_kk, beta_kk, idx_PRO, xlab, outfile,
                                      idx_circle = NULL,
                                      circle_color = "yellow2",
                                      circle_size = 4.2,
                                      circle_stroke = 1.2) {
    stopifnot(length(p_kk) == length(beta_kk))
    
    df <- tibble::tibble(
      idx  = seq_along(beta_kk),
      beta = as.numeric(beta_kk),
      p    = as.numeric(p_kk)
    ) %>%
      dplyr::mutate(
        p    = ifelse(is.na(p) | !is.finite(p), NA_real_, p),
        beta = ifelse(is.na(beta) | !is.finite(beta), NA_real_, beta),
        p_use = if (use_fdr) p.adjust(p, method = "BH") else p,
        p_use = pmax(p_use, p_floor, na.rm = FALSE),
        neglog10p = -log10(p_use),
        group = dplyr::if_else(idx %in% idx_PRO, "Sig", "Not Sig")
      )
    
    g <- ggplot(df, aes(x = beta, y = neglog10p)) +
      geom_point(
        data = filter(df, group == "Not Sig"),
        aes(color = group), size = 1.6, alpha = 0.25, na.rm = TRUE
      ) +
      geom_point(
        data = filter(df, group == "Sig"),
        aes(color = group), size = 1.9, alpha = 0.90, na.rm = TRUE
      ) +
      scale_color_manual(
        values = c("Not Sig" = "grey75", "Sig" = "red"),
        breaks = c("Not Sig", "Sig"), name = NULL
      ) +
      labs(
        x = xlab,
        y = if (use_fdr) expression(-log[10]("FDR (BH)"))
        else expression(-log[10]("P-value"))
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position   = c(0.95, 0.60),
        legend.background = element_blank(),
        axis.line  = element_line(linewidth = 0.9),
        axis.ticks = element_line(linewidth = 0.8)
      ) +
      coord_cartesian(clip = "off") +
      theme_big
    
    # Optional: circle highlighted taxa
    if (!is.null(idx_circle)) {
      idx_circle <- unique(as.integer(idx_circle))
      g <- g +
        geom_point(
          data = filter(df, idx %in% idx_circle),
          aes(x = beta, y = neglog10p),
          shape = 21, stroke = circle_stroke, size = circle_size,
          fill = NA, color = circle_color,
          inherit.aes = FALSE, na.rm = TRUE
        )
    }
    
    print(g)
    ggsave(outfile, plot = g,
           width = 8, height = 7, dpi = 300, device = ragg::agg_png)
    invisible(g)
  }
  
  # CAMRA volcano: treatment side
  plot_volcano_idxPRO_red(
    p_kk    = selected_values1,
    beta_kk = PRO_beta[1, ],
    idx_PRO = idx_PRO,
    xlab    = "Effect estimate",
    outfile = "figure/treat_vocano.png"
  )
  
  # CAMRA volcano: outcome side
  plot_volcano_idxPRO_red(
    p_kk    = selected_values2,
    beta_kk = PRO_beta[2, ],
    idx_PRO = idx_PRO,
    xlab    = "Effect estimate",
    outfile = "figure/outcome_vocano.png"
  )
  
  # LDM volcano: treatment side
  plot_volcano_idxPRO_red(
    p_kk    = p_LDM[1, ],
    beta_kk = LDM_beta[1, ],
    idx_PRO = idx_LDM,
    xlab    = "Effect estimate",
    outfile = "figure/treat_vocano_LDM.png"
  )
  
  # LDM volcano: outcome side
  plot_volcano_idxPRO_red(
    p_kk    = p_LDM[2, ],
    beta_kk = LDM_beta[2, ],
    idx_PRO = idx_LDM,
    xlab    = "Effect estimate",
    outfile = "figure/outcome_vocano_LDM.png"
  )
}