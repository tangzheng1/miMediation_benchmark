##############################################################################
# Load packages, data and souce
##############################################################################
pkgs <- c(
  "glmnet", "pscl", "plyr", "hdi", "compositions",
  "Iso", "cp4p", "HDMT", "tidyverse", "LDM",
  "harmonicmeanp", "precrec", "multimedia",
  "ccmm", "MarZIC", "HIMA", "PALM", "gtools",
  "MultiMed", "permute", "vegan", "matrixStats", "energy"
)
for (pkg in pkgs) library(pkg, character.only = TRUE)

source("LDM_fun.R")
load("GALAXYMicrobLiver_study.RData")

##############################################################################
# SECTION 1: SparCC covariance estimation
##############################################################################

SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  p <- ncol(x)
  n <- nrow(x)
  x <- x + 1  # add pseudocount for Dirichlet posterior
  
  cov.w <- cor.w <- matrix(0, p, p)
  indLow <- lower.tri(cov.w, diag = TRUE)
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax)
  
  for (i in 1:imax) {
    # Draw fractions from Dirichlet posterior
    y <- t(apply(x, 1, function(row) gtools::rdirichlet(n = 1, alpha = row)))
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin)
    covs[, i] <- cov_cor$cov.w[indLow]
    cors[, i] <- cov_cor$cor.w[indLow]
  }
  
  # Median across posterior samples
  cov.w[indLow] <- apply(covs, 1, median)
  cor.w[indLow] <- apply(cors, 1, median)
  cov.w <- cov.w + t(cov.w); diag(cov.w) <- diag(cov.w) / 2
  cor.w <- cor.w + t(cor.w); diag(cor.w) <- 1
  
  list(cov.w = cov.w, cor.w = cor.w)
}

SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  x <- log(x)
  p <- ncol(x)
  
  # T0 = variation matrix: var(log(xi/xj))
  TT <- stats::var(x)
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT
  
  # Basic SparCC: estimate basis variances
  rowT0 <- rowSums(T0)
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2)) / (p - 2)
  var.w[var.w < Vmin] <- Vmin
  
  # Compute initial correlations
  Is <- sqrt(1 / var.w)
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5
  cor.w[cor.w <= -1] <- -1
  cor.w[cor.w >=  1] <-  1
  
  # Left matrix of estimation equation
  Lmat <- diag(rep(p - 2, p)) + 1
  rp <- NULL
  cp <- rep(TRUE, p)
  
  # Iteratively exclude highly correlated pairs
  k <- 0
  while (k < kmax && sum(cp) > 3) {
    T02 <- T0
    curr_cor.w <- cor.w
    diag(curr_cor.w) <- 0
    if (!is.null(rp)) curr_cor.w[rp] <- 0
    
    n_rp <- which.max(abs(curr_cor.w))
    
    if (abs(curr_cor.w[n_rp]) >= alpha) {
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)))
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp)
      rp <- c(rp, n_rp)
      T02[rp] <- 0
      cp <- (diag(Lmat) > 0)
      
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]))
      var.w[var.w <= Vmin] <- Vmin
      
      Is <- sqrt(1 / var.w)
      cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5
      cor.w[cor.w <= -1] <- -1
      cor.w[cor.w >=  1] <-  1
    } else {
      break
    }
    k <- k + 1
  }
  
  Is <- sqrt(var.w)
  cov.w <- cor.w * Is * rep(Is, each = p)
  
  list(cov.w = cov.w, cor.w = cor.w)
}


##############################################################################
# SECTION 2: CRAmed (with full p-value table)
##############################################################################

CRAmed_fullp <- function(
    M_mat, Y, Exposure,
    FDR = 0.05,
    n.times = 100,
    prefilter = TRUE,
    n.perm = 100,
    CI = FALSE,
    modely = "gaussian",
    modelm = "ZINB",
    method = "BH"
) {
  library(glmnet)
  library(MASS)
  
  # ── Helper: build p-value table (non-selected taxa default to 1) ──
  build_p_table <- function(M_mat, selection_set, pvl_lst, final.p, method,
                            mediation_set = integer(0)) {
    p <- ncol(M_mat)
    taxa <- colnames(M_mat)
    if (is.null(taxa)) taxa <- paste0("taxa_", seq_len(p))
    
    p_d0crt <- p_model <- q_d0crt <- q_model <- p_joint <- rep(1, p)
    
    selection_set <- sort(unique(as.integer(selection_set)))
    selection_set <- selection_set[selection_set >= 1 & selection_set <= p]
    
    if (length(selection_set) > 0) {
      p_d0crt[selection_set] <- pvl_lst[selection_set]
      p_model[selection_set] <- as.numeric(final.p)
      q_d0crt[selection_set] <- p.adjust(p_d0crt[selection_set], method = method)
      q_model[selection_set] <- p.adjust(p_model[selection_set], method = method)
      p_joint[selection_set] <- pmax(q_d0crt[selection_set], q_model[selection_set])
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
  
  # ── Conditional generators for mediator models ──────────────────────
  
  # ZINB conditional generator
  Creat_condition_zinb <- function(M_mat, indx, Exposure, n.times) {
    library(pscl)
    n <- nrow(M_mat)
    Y_data <- data.frame(Exposure = Exposure, Mediator = M_mat[, indx])
    
    zinb.fit <- try(
      zeroinfl(Mediator ~ Exposure | Exposure,
               data = Y_data, dist = "negbin", link = "logit"),
      silent = TRUE
    )
    
    if (inherits(zinb.fit, "try-error")) {
      # Fall back to NB model
      zinb.fit <- try(glm.nb(Mediator ~ Exposure, data = Y_data), silent = TRUE)
      alpha <- summary(zinb.fit)$coefficients[1:2, 1]
      phi <- zinb.fit$theta
      M_bar <- zinb.fit$fitted.values
      M_res <- zinb.fit$residuals
      x <- cbind(1, Exposure)
      lamda <- exp(x %*% matrix(alpha))
      
      samp_M <- matrix(NA, n, n.times)
      for (j in 1:n.times) {
        set.seed(j)
        for (i in 1:n) samp_M[i, j] <- rnbinom(1, mu = lamda[i], size = phi)
      }
      
      res_M_samp <- lapply(seq_len(n.times), function(j) {
        Yd <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
        fit <- try(glm.nb(Mediator ~ Exposure, data = Yd), silent = TRUE)
        data.frame(t(fit$residuals))
      })
    } else {
      alpha <- summary(zinb.fit)$coefficients$count[1:2, 1]
      gamma <- summary(zinb.fit)$coefficients$zero[1:2, 1]
      phi <- zinb.fit$theta
      M_bar <- zinb.fit$fitted.values
      M_res <- zinb.fit$residuals
      
      x <- cbind(1, Exposure)
      logit.p <- x %*% matrix(gamma)
      p0 <- 1 / (1 + exp(-logit.p))
      lamda <- exp(x %*% matrix(alpha))
      
      samp_M <- matrix(NA, n, n.times)
      for (j in 1:n.times) {
        set.seed(j)
        Z1 <- rbinom(n, 1, p0)
        for (i in 1:n) {
          samp_M[i, j] <- if (Z1[i] == 1) 0 else rnbinom(1, mu = lamda[i], size = phi)
        }
      }
      
      res_M_samp <- lapply(seq_len(n.times), function(j) {
        if (sum(samp_M[, j]) == 0) return(data.frame(t(rep(0, n))))
        Yd <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
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
  
  # NB conditional generator
  Creat_condition_nb <- function(M_mat, indx, Exposure, n.times) {
    n <- nrow(M_mat)
    Y_data <- data.frame(Exposure = Exposure, Mediator = M_mat[, indx])
    nb.fit <- try(glm.nb(Mediator ~ Exposure, data = Y_data), silent = TRUE)
    
    alpha <- summary(nb.fit)$coefficients[1:2, 1]
    phi <- nb.fit$theta
    M_bar <- nb.fit$fitted.values
    M_res <- nb.fit$residuals
    
    x <- cbind(1, Exposure)
    lamda <- exp(x %*% matrix(alpha))
    
    samp_M <- matrix(NA, n, n.times)
    for (j in 1:n.times) {
      set.seed(j)
      for (i in 1:n) samp_M[i, j] <- rnbinom(1, mu = lamda[i], size = phi)
    }
    
    res_M_samp <- lapply(seq_len(n.times), function(j) {
      if (sum(samp_M[, j]) == 0) return(data.frame(t(rep(0, n))))
      Yd <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
      fit <- try(glm.nb(Mediator ~ Exposure, data = Yd), silent = TRUE)
      data.frame(t(fit$residuals))
    })
    
    list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp)
  }
  
  # ZIP conditional generator
  Creat_condition_zip <- function(M_mat, indx, Exposure, n.times) {
    library(pscl)
    n <- nrow(M_mat)
    Y_data <- data.frame(Exposure = Exposure, Mediator = M_mat[, indx])
    zip.fit <- try(
      zeroinfl(Mediator ~ Exposure | Exposure,
               data = Y_data, dist = "poisson", link = "logit"),
      silent = TRUE
    )
    
    alpha <- summary(zip.fit)$coefficients$count[1:2, 1]
    gamma <- summary(zip.fit)$coefficients$zero[1:2, 1]
    M_bar <- zip.fit$fitted.values
    M_res <- zip.fit$residuals
    
    x <- cbind(1, Exposure)
    logit.p <- x %*% matrix(gamma)
    p0 <- 1 / (1 + exp(-logit.p))
    lamda <- exp(x %*% matrix(alpha))
    
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
      Yd <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
      fit <- try(
        zeroinfl(Mediator ~ Exposure | Exposure,
                 data = Yd, dist = "poisson", link = "logit"),
        silent = TRUE
      )
      data.frame(t(fit$residuals))
    })
    
    list(mean_m = M_bar, res_m = M_res, res_m_samp = res_M_samp)
  }
  
  # ── Main CRAmed body ───────────────────────────────────────────────
  p <- ncol(M_mat)
  n <- nrow(M_mat)
  
  # Pre-filter via Lasso or use all taxa
  if (prefilter) {
    set.seed(123)
    cv_lasso <- cv.glmnet(cbind(M_mat, Exposure), Y,
                          alpha = 1, family = modely, dfmax = as.integer(p / 2))
    opt_model <- glmnet(cbind(M_mat, Exposure), Y,
                        alpha = 1, lambda = cv_lasso$lambda.min,
                        family = modely, dfmax = as.integer(p / 2))
    residualsy <- as.numeric(Y - predict(opt_model, cbind(M_mat, Exposure)))
    beta_fit <- opt_model$beta[-(1 + p)]
    beta_2 <- opt_model$beta[(1 + p)]
    selection_set <- which(beta_fit != 0)
  } else {
    y.data <- data.frame(Y = Y, Exposure = Exposure, M_mat)
    lm.fit <- lm(Y ~ ., data = y.data)
    residualsy <- lm.fit$residuals
    beta_fit <- summary(lm.fit)$coefficients[-1, 1]
    beta_2 <- beta_fit["Exposure"]
    selection_set <- seq_len(p)
  }
  
  if (length(selection_set) == 0) selection_set <- integer(0)
  
  # Default p-values for all taxa (1 = non-significant)
  pvl_lst <- rep(1, p)
  nde.p <- 1
  tij_mat <- log(rowSums(M_mat))
  
  # ── Compute d0CRT p-values ─────────────────────────────────────────
  if (p == 1) {
    selection_set <- 1
    Cond_M <- switch(modelm,
                     ZINB = Creat_condition_zinb(M_mat, 1, Exposure, n.times),
                     ZIP  = Creat_condition_zip(M_mat, 1, Exposure, n.times),
                     NB   = Creat_condition_nb(M_mat, 1, Exposure, n.times)
    )
    M_res_ob <- Cond_M$res_m
    data0 <- data.frame(Y = Y, Exposure = Exposure)
    if (modely == "binomial") {
      eps_res <- Y - 1 / (1 + exp(-predict(lm(Y ~ Exposure, family = modely, data = data0))))
    } else {
      eps_res <- Y - predict(lm(Y ~ Exposure, family = modely, data = data0))
    }
    imp_obe <- abs(mean(M_res_ob * eps_res)) / mean(M_res_ob^2)
    
    list.s <- unlist(lapply(Cond_M$res_m_samp, length))
    list.sel <- which(list.s == 1)
    M_res_sample <- if (length(list.sel) != 0) {
      as.matrix(do.call(rbind.fill, Cond_M$res_m_samp[-list.sel]))
    } else {
      as.matrix(do.call(rbind.fill, Cond_M$res_m_samp))
    }
    var_lst_sample <- apply(M_res_sample, 1, function(v) mean((unlist(v))^2))
    t_lst <- abs(M_res_sample %*% eps_res / n) / var_lst_sample
    pvl_lst[1] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
    
  } else {
    # NDE p-value via exposure residual permutation
    Cond_E <- list(res_e_samp = list())
    for (j in 1:n.times) {
      set.seed(j)
      indx.e <- sample(1:n)
      Cond_E$res_e_samp[[j]] <- as.data.frame(t(Exposure[indx.e]))
    }
    
    set.seed(123)
    cv_lasso_null <- cv.glmnet(M_mat, Y, alpha = 1, family = modely,
                               dfmax = as.integer(p / 2))
    model_res_null <- glmnet(M_mat, Y, alpha = 1,
                             lambda = cv_lasso_null$lambda.min,
                             family = modely, dfmax = as.integer(p / 2))
    
    if (modely == "binomial") {
      eps_res <- Y - 1 / (1 + exp(-predict(model_res_null, M_mat)))
    } else {
      eps_res <- as.numeric(Y - predict(model_res_null, M_mat))
    }
    imp_exp <- abs(mean(Exposure * eps_res)) / mean(Exposure^2)
    
    list.s <- unlist(lapply(Cond_E$res_e_samp, length))
    list.sel <- which(list.s == 1)
    E_res_sample <- if (length(list.sel) != 0) {
      as.matrix(do.call(rbind.fill, Cond_E$res_e_samp[-list.sel]))
    } else {
      as.matrix(do.call(rbind.fill, Cond_E$res_e_samp))
    }
    var_lst_sample <- apply(E_res_sample, 1, function(v) mean((unlist(v))^2))
    t_lst <- abs(E_res_sample %*% eps_res / n) / var_lst_sample
    nde.p <- mean(c(1, ifelse(t_lst >= imp_exp, 1, 0)))
    
    # Taxon-level d0CRT p-values (only for selected taxa)
    if (length(selection_set) > 0) {
      for (j in seq_along(selection_set)) {
        indx <- selection_set[j]
        Cond_M <- switch(modelm,
                         ZINB = Creat_condition_zinb(M_mat, indx, Exposure, n.times),
                         ZIP  = Creat_condition_zip(M_mat, indx, Exposure, n.times),
                         NB   = Creat_condition_nb(M_mat, indx, Exposure, n.times)
        )
        M_res_ob <- Cond_M$res_m
        
        set.seed(123)
        cv_null2 <- cv.glmnet(cbind(M_mat[, -indx, drop = FALSE], Exposure), Y,
                              alpha = 1, family = modely, dfmax = as.integer(p / 2))
        model_null2 <- glmnet(cbind(M_mat[, -indx, drop = FALSE], Exposure), Y,
                              alpha = 1, lambda = cv_null2$lambda.min,
                              family = modely, dfmax = as.integer(p / 2))
        
        if (modely == "binomial") {
          eps_res2 <- Y - 1 / (1 + exp(-predict(model_null2, cbind(M_mat[, -indx, drop = FALSE], Exposure))))
        } else {
          eps_res2 <- as.numeric(Y - predict(model_null2, cbind(M_mat[, -indx, drop = FALSE], Exposure)))
        }
        imp_obe <- abs(mean(M_res_ob * eps_res2)) / mean(M_res_ob^2)
        
        list.s <- unlist(lapply(Cond_M$res_m_samp, length))
        list.sel <- which(list.s == 1)
        M_res_sample <- if (length(list.sel) != 0) {
          as.matrix(do.call(rbind.fill, Cond_M$res_m_samp[-list.sel]))
        } else {
          as.matrix(do.call(rbind.fill, Cond_M$res_m_samp))
        }
        var_lst_sample <- apply(M_res_sample, 1, function(v) mean((unlist(v))^2))
        t_lst <- abs(M_res_sample %*% eps_res2 / n) / var_lst_sample
        pvl_lst[indx] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
      }
    }
  }
  
  # ── Fit mediator models and compute NIE ────────────────────────────
  aic <- bic <- final.p <- niea <- niep <- nie <- alpha.p <- gamma.p <- numeric(0)
  residualsm <- list()
  alpha_mat <- gamma_mat <- NULL
  
  if (length(selection_set) > 0) {
    ns <- length(selection_set)
    aic <- bic <- final.p <- niea <- niep <- nie <- rep(NA_real_, ns)
    alpha.p <- gamma.p <- rep(NA_real_, ns)
    residualsm <- vector("list", ns)
    alpha_mat <- matrix(NA_real_, ns, 2)
    gamma_mat <- matrix(NA_real_, ns, 2)
  }
  
  mediation_set <- integer(0)
  
  # ── Model-specific fitting (ZINB / ZIP / NB) ──────────────────────
  fit_mediator_models <- function(modelm) {
    if (length(selection_set) == 0) return(NULL)
    
    for (j in seq_along(selection_set)) {
      Y_data <- data.frame(Exposure = Exposure,
                           Mediator = as.vector(M_mat[, selection_set[j]]))
      
      if (modelm == "NB" || sum(Y_data$Mediator == 0) == 0) {
        # Pure NB model (no zero inflation needed)
        nb.fit <- glm.nb(Mediator ~ Exposure, data = Y_data)
        residualsm[[j]] <<- nb.fit$residuals
        final.p[j] <<- summary(nb.fit)$coefficients[2, 4]
        alpha_mat[j, ] <<- summary(nb.fit)$coefficients[1:2, 1]
        alpha.p[j] <<- summary(nb.fit)$coefficients[2, 4]
        if (modelm == "NB") {
          aic[j] <<- AIC(nb.fit)
          bic[j] <<- BIC(nb.fit)
        }
      } else {
        # Zero-inflated model
        dist <- if (modelm == "ZINB") "negbin" else "poisson"
        offset_arg <- if (sum(is.infinite(tij_mat)) != 0) NULL else tij_mat
        
        zi.fit <- try(
          if (is.null(offset_arg)) {
            zeroinfl(Mediator ~ Exposure | Exposure,
                     data = Y_data, dist = dist, link = "logit")
          } else {
            zeroinfl(Mediator ~ Exposure | Exposure,
                     offset = offset_arg, data = Y_data, dist = dist, link = "logit")
          },
          silent = TRUE
        )
        
        if (inherits(zi.fit, "try-error")) {
          final.p[j] <<- NA
          aic[j] <<- bic[j] <<- NA
          residualsm[[j]] <<- NA
        } else {
          aic[j] <<- AIC(zi.fit)
          bic[j] <<- BIC(zi.fit)
          residualsm[[j]] <<- zi.fit$residuals
          
          alpha_mat[j, ] <<- summary(zi.fit)$coefficients$count[1:2, 1]
          gamma_mat[j, ] <<- summary(zi.fit)$coefficients$zero[1:2, 1]
          alpha.p[j] <<- summary(zi.fit)$coefficients$count[2, 4]
          gamma.p[j] <<- summary(zi.fit)$coefficients$zero[2, 4]
          
          # Wald test for joint significance
          cov.mat <- matrix(c(
            vcov(zi.fit)[2, 2], vcov(zi.fit)[2, 4],
            vcov(zi.fit)[4, 2], vcov(zi.fit)[4, 4]
          ), 2, 2)
          
          wald.t <- try(
            t(matrix(c(alpha_mat[j, 2], gamma_mat[j, 2]))) %*%
              ginv(cov.mat) %*% matrix(c(alpha_mat[j, 2], gamma_mat[j, 2])),
            silent = TRUE
          )
          final.p[j] <<- if (inherits(wald.t, "try-error")) NA else
            (1 - pchisq(as.numeric(wald.t), df = 2))
        }
      }
      
      # Compute NIE components
      beta.val <- beta_fit[selection_set[j]]
      alpha.val <- alpha_mat[j, ]
      gamma.val <- gamma_mat[j, ]
      ones <- rep(1, nrow(M_mat))
      
      if (modelm == "NB") {
        nie[j] <<- mean(beta.val * (
          exp(cbind(1, ones) %*% matrix(alpha.val[1:2]) + tij_mat) -
            exp(cbind(ones) %*% matrix(alpha.val[1]) + tij_mat)
        ))
      } else {
        # ZINB / ZIP: decompose into abundance and prevalence components
        niea[j] <<- mean(beta.val *
                           (1 / (1 + exp(cbind(ones) %*% matrix(gamma.val[1])))) *
                           (exp(cbind(1, ones) %*% matrix(alpha.val[1:2]) + tij_mat) -
                              exp(cbind(ones) %*% matrix(alpha.val[1]) + tij_mat)))
        
        niep[j] <<- mean(beta.val *
                           (exp(cbind(1, ones) %*% matrix(alpha.val[1:2]) + tij_mat)) *
                           ((1 / (1 + exp(cbind(1, ones) %*% matrix(gamma.val[1:2])))) -
                              (1 / (1 + exp(cbind(ones) %*% matrix(gamma.val[1]))))))
        
        nie[j] <<- niea[j] + niep[j]
      }
    }
  }
  
  fit_mediator_models(modelm)
  
  # ── Selection based on joint p-value ───────────────────────────────
  index.p <- if (length(selection_set) > 0) pvl_lst[selection_set] else numeric(0)
  bind.p <- rbind(p.adjust(index.p, method), p.adjust(final.p, method))
  joint.p <- if (length(selection_set) > 0) apply(bind.p, 2, max) else numeric(0)
  
  index.mi <- which(joint.p <= FDR)
  mediation_set <- if (length(index.mi) > 0) selection_set[index.mi] else integer(0)
  
  # Extract results for selected mediators
  nie_keep  <- if (length(index.mi) > 0) nie[index.mi] else numeric(0)
  niea_keep <- if (length(index.mi) > 0) niea[index.mi] else numeric(0)
  niep_keep <- if (length(index.mi) > 0) niep[index.mi] else numeric(0)
  niepval_keep <- if (length(index.mi) > 0) joint.p[index.mi] else numeric(0)
  
  # Build full p-value table
  p_table <- build_p_table(M_mat, selection_set, pvl_lst, final.p, method,
                           mediation_set)
  
  list(
    Mediators = mediation_set,
    NDE = beta_2, NIE = nie_keep, NIEA = niea_keep, NIEP = niep_keep,
    NDE.pval = nde.p, NIE.pval = niepval_keep,
    AIC = aic, BIC = bic,
    residualsy = residualsy, residualsm = residualsm,
    p_table = p_table
  )
}


##############################################################################
# SECTION 3: CAMRA core functions
##############################################################################

recover_l_PALM <- function(count_m, treat_cov, cov_ad = NULL,
                           prev.filter = 0, eps_p = 1e-10) {
  count_m <- as.matrix(count_m)
  storage.mode(count_m) <- "numeric"
  n <- nrow(count_m)
  p <- ncol(count_m)
  
  orig_taxa <- paste0("O", seq_len(p))
  colnames(count_m) <- orig_taxa
  rn <- paste0("T", seq_len(n))
  rownames(count_m) <- rn
  
  stopifnot(length(treat_cov) == n)
  treat_cov <- matrix(as.numeric(treat_cov), ncol = 1)
  colnames(treat_cov) <- "treat"
  rownames(treat_cov) <- rn
  
  if (!is.null(cov_ad)) {
    cov_ad <- as.matrix(cov_ad)
    stopifnot(nrow(cov_ad) == n)
    storage.mode(cov_ad) <- "numeric"
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
  
  # Initialize full-length outputs
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

recover_r <- function(count_matrix, treat_cov, y, sudo = 0.5, cov_ad = NULL,
                      CClasso = FALSE, cov_true = NULL) {
  
  # Compositional transform + covariance estimation
  logdata  <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
  por_data <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
  
  if (CClasso) {
    est_cov <- {
      res_cov <- fastCCLasso(count_matrix, isCnt = TRUE)
      diag(sqrt(res_cov$cov_diag)) %*% res_cov$rho %*% diag(sqrt(res_cov$cov_diag))
    }
  } else {
    res_l  <- SparCC.count(count_matrix)
    est_cov <- res_l$cov.w
  }
  
  if (!is.null(cov_true)) est_cov <- cov_true
  
  # Build PALAR-transformed predictors Z_ilr
  p <- ncol(count_matrix)
  n <- nrow(count_matrix)
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
  } else {
    treat_df
  }
  
  Z0 <- model.matrix(~ ., data = Zdf)
  
  # QR-based residualization
  qrZ <- qr(Z0)
  y_tilde <- as.numeric(y - qr.fitted(qrZ, y))
  X_tilde <- Z_ilr - qr.fitted(qrZ, Z_ilr)
  
  # De-biased ridge projection on residualized high-dimensional part
  outRidge <- hdi::ridge.proj(x = X_tilde, y = y_tilde)
  
  all_p  <- outRidge$pval
  beta_r <- as.vector(outRidge$bhat)
  z      <- as.vector(qnorm(1 - all_p / 2) * sign(beta_r))
  
  list(p = all_p, z = z, beta_r = beta_r,
       y_tilde = y_tilde, X_tilde = X_tilde,
       Z0 = Z0, X_doubel = Z_ilr)
}

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
  
  # Log-ratio transform
  logdata <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
  logdata[logdata < (-10)] <- (-10)
  
  # Covariance estimation via SparCC
  res_l  <- SparCC.count(count_matrix)
  est_cov <- res_l$cov.w
  
  por_data <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
  ilr_basis <- compositions::ilrBase(por_data)
  R2 <- ilr_basis
  Z_ilr <- (logdata %*% R2) %*%
    solve(t(R2) %*% est_cov %*% R2) %*% t(R2) %*% est_cov
  
  # Build design matrix: mediators (penalized) + treatment/confounders (unpenalized)
  treat_vec <- as.numeric(treat_cov)
  Z0 <- cbind(treat = treat_vec)
  if (!is.null(cov_ad)) {
    Z0 <- cbind(Z0, cov_ad)
    colnames(Z0) <- make.names(colnames(Z0), unique = TRUE)
  }
  
  X <- cbind(Z_ilr, Z0)
  pZ <- ncol(Z_ilr)
  p0 <- ncol(Z0)
  pf <- c(rep(1, pZ), rep(0, p0))  # penalty factors: 1=penalized, 0=unpenalized
  
  if (isTRUE(adaptive_L)) {
    cvfit <- glmnet::cv.glmnet(
      x = X, y = y, alpha = 1, penalty.factor = pf,
      nfolds = 5, type.measure = "mse", standardize = TRUE
    )
    b <- as.matrix(coef(cvfit, s = "lambda.min"))
    beta_Z <- as.numeric(b)[-1][1:pZ]
    selection_set <- which(beta_Z != 0)
    if (length(selection_set) == 0) {
      selection_set <- order(abs(beta_Z), decreasing = TRUE)[1]
    }
    return(sort(unique(as.integer(selection_set))))
  }
  
  # Fixed-K selection rule: K = const * n / log(n)
  K_raw <- floor(const * n / log(max(n, 3)))
  K <- max(1L, min(pZ, K_raw))
  
  fit <- glmnet::glmnet(
    x = X, y = y, alpha = 1, penalty.factor = pf,
    dfmax = min(K + p0, pZ + p0), nlambda = 500,
    lambda.min.ratio = 1e-6, standardize = TRUE
  )
  
  B <- as.matrix(fit$beta)
  dfZ <- colSums(B[1:pZ, , drop = FALSE] != 0)
  idx <- which(dfZ >= K)[1]
  if (is.na(idx)) idx <- ncol(B)
  
  beta_Z_full <- as.numeric(B[1:pZ, idx])
  ord <- order(abs(beta_Z_full), decreasing = TRUE)
  keep <- ord[seq_len(min(K, length(ord)))]
  
  sort(unique(as.integer(keep)))
}

p_mediation_maxp <- function(p_alpha, p_beta,
                             pi_alpha0 = NULL, pi_beta0 = NULL,
                             pi_method = c("cp4p", "JC"),
                             weight_method = c("maxp", "product", "indenp")) {
  stopifnot(length(p_alpha) == length(p_beta))
  pi_method     <- match.arg(pi_method)
  weight_method <- match.arg(weight_method)
  
  # ── Helper: product-based mixture weights ──
  mix_weights_product <- function(pi_alpha0, pi_beta0) {
    eps <- 1e-8
    pa <- min(max(pi_alpha0, eps), 1 - 1e-6)
    pb <- min(max(pi_beta0,  eps), 1 - 1e-6)
    pi0 <- max(1 - (1 - pa) * (1 - pb), 1e-6)
    w00 <- (pa * pb) / pi0
    w10 <- ((1 - pa) * pb) / pi0
    w01 <- (pa * (1 - pb)) / pi0
    c(w00 = w00, w10 = w10, w01 = w01, pi0 = pi0)
  }
  
  # ── Helper: max-p-based mixture weights ──
  mix_weights_maxp <- function(p_alpha, p_beta, pi_alpha0, pi_beta0, pi_method) {
    p_max <- pmax(p_alpha, p_beta)
    if (pi_method == "cp4p") {
      pi0_hat <- {
        obj <- cp4p::estim.pi0(p_max)
        if (!is.null(obj$pi0)) as.numeric(obj$pi0) else mean(unlist(obj))
      }
    } else {
      zmax <- qnorm(1 - p_max)
      pi0_hat <- miMediation:::.pi0_JC(zmax)
    }
    
    clip01 <- function(x) min(max(x, 1e-6), 1 - 1e-6)
    pi_alpha0 <- clip01(pi_alpha0)
    pi_beta0  <- clip01(pi_beta0)
    pi0_hat   <- clip01(pi0_hat)
    
    # MaxP decomposition (conditional on "mediation null")
    w00 <- (pi_alpha0 + pi_beta0 - pi0_hat) / pi0_hat
    w10 <- (pi0_hat - pi_alpha0) / pi0_hat
    w01 <- (pi0_hat - pi_beta0)  / pi0_hat
    
    w <- pmax(c(w00, w10, w01), 0)
    w <- w / sum(w)
    names(w) <- c("w00", "w10", "w01")
    w
  }
  
  # ── Step 1: Clean p-values for null proportion estimation ──
  p_alpha_pi0 <- p_alpha; p_beta_pi0 <- p_beta
  p_alpha_pi0[!is.finite(p_alpha_pi0)] <- 1
  p_beta_pi0[!is.finite(p_beta_pi0)]   <- 1
  p_alpha_pi0 <- pmin(pmax(p_alpha_pi0, 0), 1)
  p_beta_pi0  <- pmin(pmax(p_beta_pi0, 0), 1)
  
  # ── Step 2: Estimate marginal null proportions ──
  if (is.null(pi_alpha0) || is.null(pi_beta0)) {
    if (pi_method == "cp4p") {
      pa <- cp4p::estim.pi0(p_alpha_pi0)
      pb <- cp4p::estim.pi0(p_beta_pi0)
      grab <- function(x) {
        if (!is.null(x$pi0)) as.numeric(x$pi0) else mean(unlist(x), na.rm = TRUE)
      }
      pi_alpha0 <- grab(pa); pi_beta0 <- grab(pb)
    } else {
      z1 <- qnorm(1 - p_alpha_pi0)
      z2 <- qnorm(1 - p_beta_pi0)
      pi_alpha0 <- miMediation:::.pi0_JC(z1)
      pi_beta0  <- miMediation:::.pi0_JC(z2)
    }
  }
  
  eps <- 1e-8
  pi_alpha0 <- min(max(pi_alpha0, eps), 1 - eps)
  pi_beta0  <- min(max(pi_beta0,  eps), 1 - eps)
  
  # ── Step 3: Valid positions & cleaned p-values ──
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  out  <- rep(NA_real_, length(p_alpha))
  if (!any(keep)) return(out)
  
  p_a <- pmin(pmax(p_alpha[keep], eps), 1 - eps)
  p_b <- pmin(pmax(p_beta[keep],  eps), 1 - eps)
  
  # ── Step 4: Estimate mixture weights ──
  if (weight_method == "maxp") {
    w <- mix_weights_maxp(p_a, p_b, pi_alpha0, pi_beta0, pi_method)
  } else if (weight_method == "product") {
    w_raw <- mix_weights_product(pi_alpha0, pi_beta0)
    w_vec <- pmax(w_raw[c("w00", "w10", "w01")], 0)
    w <- w_vec / sum(w_vec)
    names(w) <- c("w00", "w10", "w01")
  } else {
    w_vec <- c(w00 = pi_alpha0 * pi_beta0,
               w10 = (1 - pi_alpha0) * pi_beta0,
               w01 = pi_alpha0 * (1 - pi_beta0))
    w_vec <- pmax(w_vec, 0)
    w <- w_vec / sum(w_vec)
  }
  
  w00 <- as.numeric(w["w00"])
  w10 <- as.numeric(w["w10"])
  w01 <- as.numeric(w["w01"])
  
  # ── Step 5: Grenander-based alternative CDF estimation ──
  estimate_F1_grenander <- function(p, pi0, eps = 1e-8) {
    p <- p[is.finite(p)]
    p <- pmin(pmax(p, 0), 1)
    n <- length(p); stopifnot(n > 0)
    pi0 <- min(max(pi0, 1e-6), 1 - 1e-6)
    
    x <- sort(unique(c(0, sort(p), 1)))
    Fn <- ecdf(p); y <- Fn(x)
    dx <- diff(x); keep <- dx > eps
    xR <- x[-1][keep]; dx <- dx[keep]
    yL <- y[-length(y)][keep]; yR <- y[-1][keep]
    
    s <- (yR - yL) / dx
    s_hat <- -Iso::pava(-s, w = dx)   # enforce monotone decreasing
    f1_hat <- pmax((s_hat - pi0) / (1 - pi0), 0)
    area <- sum(f1_hat * dx)
    if (area <= 0) return(function(t) rep(0, length(t)))
    f1_hat <- f1_hat / area
    
    x_knots <- c(0, xR)
    F1_cum  <- c(0, cumsum(f1_hat * dx))
    function(t) {
      t <- pmin(pmax(t, 0), 1)
      pmin(pmax(approx(x_knots, F1_cum, xout = t, method = "linear",
                       ties = "ordered", rule = 2)$y, 0), 1)
    }
  }
  
  F1a <- estimate_F1_grenander(p_a, pi_alpha0)
  F1b <- estimate_F1_grenander(p_b, pi_beta0)
  
  # ── Step 6: Compute mixture-null p-values ──
  t <- pmax(p_a, p_b)
  p_mix <- w00 * (t^2) + w10 * (t * F1a(t)) + w01 * (t * F1b(t))
  p_mix <- pmin(pmax(p_mix, 0), 1)
  
  out[keep] <- p_mix
  out
}

p_mediation_hdmt_fdr <- function(p_alpha, p_beta, exact_p = 1) {
  stopifnot(length(p_alpha) == length(p_beta))
  n   <- length(p_alpha)
  out <- rep(NA_real_, n)
  
  keep <- is.finite(p_alpha) & is.finite(p_beta)
  if (!any(keep)) return(out)
  
  pa <- pmin(pmax(p_alpha[keep], 0), 1)
  pb <- pmin(pmax(p_beta[keep],  0), 1)
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
# SECTION 4: Competing method wrappers
##############################################################################

# ── LDM-med ─────────────────────────────────────────────────────────────────
ldm_sim <- function(count_m, treat_cov, y) {
  # Generate random names to avoid global-env conflicts in parallel
  rand_id <- paste0(sample(letters, 12), collapse = "")
  mat_name  <- paste0("M_mat_", rand_id)
  meta_name <- paste0("meta_", rand_id)
  
  assign(mat_name, as.matrix(count_m), envir = .GlobalEnv)
  meta_df <- data.frame(trt = treat_cov, Y = y)
  assign(meta_name, meta_df, envir = .GlobalEnv)
  
  on.exit({
    if (exists(mat_name, envir = .GlobalEnv)) rm(list = mat_name, envir = .GlobalEnv)
    if (exists(meta_name, envir = .GlobalEnv)) rm(list = meta_name, envir = .GlobalEnv)
  }, add = TRUE)
  
  fmla_call <- parse(text = paste0(mat_name, " ~ trt + Y"))[[1]]
  
  tm <- system.time({
    library(permute); library(vegan); library(matrixStats)
    res <- ldm_new(formula = fmla_call,
                   data = get(meta_name, envir = .GlobalEnv),
                   seed = 67817, fdr.nominal = 0.05,
                   test.mediation = TRUE)
  })
  runtime_sec <- unname(tm["elapsed"])
  p_global <- res$med.p.global.omni
  
  # Extract per-taxon mediation p-values (maxP approach)
  P <- as.matrix(res$p.otu.omni)
  p_joint <- MultiMed::medTest.SBMH(P[1, ], P[2, ],
                                    MCP.type = "FDR", t1 = 0.05 / 2, t2 = 0.05 / 2)
  
  rn <- rownames(P)
  if (!is.null(rn) && all(c("trt", "Y") %in% rn)) {
    p_EM <- as.numeric(P["trt", ]); p_MY <- as.numeric(P["Y", ])
  } else {
    p_EM <- as.numeric(P[1, ]); p_MY <- as.numeric(P[2, ])
  }
  
  p_med <- pmax(p_EM, p_MY)
  taxa_names <- colnames(P)
  if (is.null(taxa_names)) taxa_names <- paste0("taxon_", seq_len(ncol(P)))
  names(p_med) <- taxa_names
  
  # Identify discoveries
  det <- res$med.detected.otu.omni
  discoveries <- integer(ncol(P))
  names(discoveries) <- taxa_names
  det_sel <- rep(FALSE, ncol(P))
  
  if (is.logical(det)) {
    det_sel <- det; discoveries <- as.integer(det)
  } else if (is.numeric(det)) {
    det <- as.integer(det)
    valid_det <- det[det >= 1 & det <= ncol(P)]
    if (length(valid_det) > 0) { det_sel[valid_det] <- TRUE; discoveries[valid_det] <- 1L }
  } else if (is.character(det)) {
    det_sel <- taxa_names %in% det; discoveries[det_sel] <- 1L
  }
  
  discoveries <- which(discoveries == 1)
  
  list(discoveries = discoveries, p_med = p_joint,
       runtime_sec = runtime_sec, global_p = p_global)
}

# ── CRAmed wrapper ──────────────────────────────────────────────────────────
CRAmed_sim <- function(count1, treat1, y1) {
  t0 <- proc.time()
  results <- CRAmed_fullp(
    M_mat = count1, Y = y1, Exposure = as.matrix(treat1),
    n.times = 100, n.perm = 100, CI = FALSE
  )
  p_joint <- setNames(results$p_table$p_joint, results$p_table$taxa)
  runtime <- unname((proc.time() - t0)["elapsed"])
  list(discoveries = results$Mediators, p_med = p_joint, runtime_sec = runtime)
}


# ── MarZIC ──────────────────────────────────────────────────────────────────
Mar_sim <- function(count1, treat1, y1, alpha = 0.05,
                    adjust_method = "fdr", num_cores = 1) {
  t0 <- proc.time()
  stopifnot(nrow(count1) == length(treat1), nrow(count1) == length(y1))
  
  MicrobData <- as.matrix(count1)
  mode(MicrobData) <- "numeric"
  orig_taxa <- colnames(MicrobData)
  if (is.null(orig_taxa)) orig_taxa <- paste0("Taxon", seq_len(ncol(MicrobData)))
  colnames(MicrobData) <- orig_taxa
  
  MicrobData[!is.finite(MicrobData)] <- 0
  
  # Library size as covariate
  libsize <- pmax(1, rowSums(MicrobData))
  CovData <- data.frame(Y = as.numeric(y1), X = as.numeric(treat1), libsize = libsize)
  
  # Remove rows with invalid Y/X or all-zero counts
  keep_row <- is.finite(CovData$Y) & is.finite(CovData$X) & (rowSums(MicrobData) > 0)
  MicrobData2 <- MicrobData[keep_row, , drop = FALSE]
  CovData2    <- CovData[keep_row, , drop = FALSE]
  
  # Remove all-zero taxa
  keep_col <- colSums(MicrobData2) > 0
  MicrobData2 <- MicrobData2[, keep_col, drop = FALSE]
  MicrobData2 <- t(apply(MicrobData2, 1, function(x) x / sum(x)))
  
  # Run MarZIC without taxon filtering
  res <- MarZIC::MarZIC(
    MicrobData = MicrobData2, CovData = CovData2,
    lib_name = "libsize", y_name = "Y", x_name = "X",
    x4_inter = FALSE, x5_inter = FALSE, conf_name = NULL,
    taxa_of_interest = "all", num_cores = num_cores,
    adjust_method = adjust_method,
    zero_prop_NIE2 = 0.001, zero_count_NIE2 = 1,
    taxDropThresh = 0.999, taxDropCount = 1, SDThresh = 0.001
  )
  
  # Extract adjusted p-values
  NIE_tbl <- as.data.frame(res$NIE_save %||% res$NIE)
  pcol <- if ("Adjusted p value" %in% colnames(NIE_tbl)) "Adjusted p value" else {
    cand <- grep("p\\s*value|pval|adj", colnames(NIE_tbl), ignore.case = TRUE, value = TRUE)
    if (length(cand) == 0) stop("Cannot locate p-value column in MarZIC output")
    cand[1]
  }
  
  # Map back to original taxa order
  p_raw <- as.numeric(NIE_tbl[[pcol]])
  taxa_in_tbl <- rownames(NIE_tbl)
  if (is.null(taxa_in_tbl)) taxa_in_tbl <- NIE_tbl[[1]]
  
  pvals <- rep(NA_real_, length(orig_taxa)); names(pvals) <- orig_taxa
  m <- match(orig_taxa, taxa_in_tbl)
  pvals[!is.na(m)] <- p_raw[m[!is.na(m)]]
  
  # Restore full p-values (filtered taxa get 1)
  kept_taxa <- colnames(MicrobData2)
  p_med <- res$NIE_save$`p value adj`
  p_full <- rep(1, length(orig_taxa)); names(p_full) <- orig_taxa
  idx_match <- match(kept_taxa, orig_taxa)
  p_med_clean <- as.numeric(p_med)
  p_med_clean[!is.finite(p_med_clean) | is.na(p_med_clean)] <- 1
  p_full[idx_match] <- p_med_clean
  
  sig_idx <- which(!is.na(pvals) & pvals < alpha)
  runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
  
  list(discoveries = sig_idx, p_med = p_full, runtime_sec = runtime_sec)
}


# ── microHIMA ───────────────────────────────────────────────────────────────
hima_discrete_q <- function(p_vec, taxa_names = NULL,
                            alpha_grid = c(0.01, 0.02, 0.05, 0.10, 0.20),
                            eps = 1e-4, simes = FALSE) {
  p_vec <- as.numeric(p_vec)
  m <- length(p_vec)
  if (is.null(taxa_names)) {
    taxa_names <- names(p_vec)
    if (is.null(taxa_names)) taxa_names <- paste0("taxon_", seq_len(m))
  }
  
  q_disc <- rep(1, m)
  bad <- !is.finite(p_vec) | p_vec < 0 | p_vec > 1
  q_disc[bad] <- NA_real_
  
  ok <- which(!bad)
  if (length(ok) == 0) {
    return(list(
      table = data.frame(Index = taxa_names, p = p_vec, q_value = q_disc),
      selected = setNames(vector("list", length(alpha_grid)), as.character(alpha_grid))
    ))
  }
  
  p_ok <- p_vec[ok]
  hom <- hommel::hommel(p_ok, simes = simes)
  alpha_grid <- sort(unique(alpha_grid))
  selected_list <- setNames(vector("list", length(alpha_grid)), as.character(alpha_grid))
  
  for (alpha in alpha_grid) {
    set <- which(p_ok < alpha)
    if (length(set) == 0) { selected_list[[as.character(alpha)]] <- character(0); next }
    
    N1 <- hommel::discoveries(hom, set, incremental = TRUE, alpha = alpha)
    L <- length(set)
    N2 <- numeric(L)
    if (L >= 2) N2[2:L] <- N1[1:(L - 1)]
    N0 <- N1 - N2
    ID_local  <- set[which(N0 > 0)]
    ID_global <- ok[ID_local]
    
    selected_list[[as.character(alpha)]] <- taxa_names[ID_global]
    newly <- ID_global[is.finite(q_disc[ID_global]) & q_disc[ID_global] >= 1]
    if (length(newly) > 0) q_disc[newly] <- alpha - eps
  }
  
  list(
    table = data.frame(Index = taxa_names, p = p_vec, q_value = q_disc),
    selected = selected_list
  )
}

HIMA_micro_sim1 <- function(count1, X, Y, COV = NULL,
                            verbose = TRUE, parallel = FALSE, ncore = 1) {
  pseudo <- 0.5
  OTU_comp <- sweep(count1 + pseudo, 1, rowSums(count1 + pseudo), "/")
  
  # Modified HIMA function with q-value recovery
  HIMA_recover <- function(X, OTU, Y, COV = NULL, FDRcut = 0.05,
                           verbose = FALSE, parallel = FALSE, ncore = 1) {
    X <- matrix(X, ncol = 1)
    M_raw <- as.matrix(OTU)
    M_ID_name <- colnames(M_raw)
    if (is.null(M_ID_name)) M_ID_name <- seq_len(ncol(M_raw))
    if (!is.null(COV)) { COV <- as.matrix(COV); X <- cbind(X, COV) }
    X <- scale(X); Y <- Y - mean(Y)
    M <- M_raw; n <- dim(M)[1]; d <- dim(M)[2]
    M1 <- t(t(M_raw[, 1]))
    
    if (verbose) message("Step 1: ILR Transformation and De-biased Lasso estimates ...")
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
      lm.out <- summary(stats::lm(MT[, 1] ~ X))
      alpha_est <- lm.out$coefficients[2, 1]; alpha_se <- lm.out$coefficients[2, 2]
      P_a <- 2 * (1 - pnorm(abs(alpha_est / alpha_se), 0, 1))
      c(beta_est, beta_se, alpha_est, alpha_se, max(P_a, P_b), P_a, P_b)
    }
    
    if (is.null(dim(results_loop))) results_loop <- matrix(results_loop, nrow = 1)
    P_raw_DLASSO <- results_loop[, 5]
    
    # Discrete q-values via hommel procedure
    pmax_res <- hima_discrete_q(P_raw_DLASSO)
    list(q_value = pmax_res$table$q_value, select_list = pmax_res$selected)
  }
  
  t0 <- proc.time()
  fit <- HIMA_recover(X = X, OTU = OTU_comp, Y = Y, COV = COV,
                      verbose = verbose, parallel = parallel, ncore = ncore)
  runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
  
  list(index = NULL, q_value = fit$q_value, runtime_sec = runtime_sec)
}


# ── multimedia ──────────────────────────────────────────────────────────────

multimedia_sim <- function(count1, treat1, y1, q_value = 0.05,
                           add_lib_size = TRUE, pseudo = 0.5, alpha_out = 1) {
  t0 <- proc.time()
  stopifnot(length(treat1) == nrow(count1), length(y1) == nrow(count1))
  
  orig_names_full <- colnames(count1)
  if (is.null(orig_names_full)) orig_names_full <- paste0("V", seq_len(ncol(count1)))
  
  X0 <- as.matrix(count1); storage.mode(X0) <- "numeric"
  
  # Drop all-zero samples
  rs <- rowSums(X0); keep_samp <- rs > 0
  if (!all(keep_samp)) {
    X0 <- X0[keep_samp, , drop = FALSE]
    treat1 <- treat1[keep_samp]; y1 <- y1[keep_samp]; rs <- rs[keep_samp]
  }
  
  # Drop all-zero / zero-variance taxa
  keep_taxa <- (colSums(X0) > 0) & (apply(X0, 2, var) > 0)
  if (sum(keep_taxa) == 0) {
    runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
    p_full <- rep(1, length(orig_names_full)); names(p_full) <- orig_names_full
    return(list(discoveries_default = integer(0), q_value_disc = p_full,
                runtime_sec = runtime_sec, discoveries_by_q = list()))
  }
  
  X_use <- X0[, keep_taxa, drop = FALSE]
  orig_names_use <- orig_names_full[keep_taxa]
  
  # CLR transformation
  X_use <- X_use + pseudo
  prop <- sweep(X_use, 1, rowSums(X_use), "/")
  logp <- log(prop)
  clrM <- logp - rowMeans(logp)
  
  safe_names <- make.names(orig_names_use, unique = TRUE)
  colnames(clrM) <- safe_names
  mediator <- safe_names
  
  # Build data frame
  treat01 <- as.numeric(treat1)
  stopifnot(all(treat01 %in% c(0, 1)), length(unique(treat01)) >= 2)
  treat12 <- treat01 + 1  # 0/1 -> 1/2
  
  df <- data.frame(treatment = treat12, outcome = as.numeric(y1), check.names = FALSE)
  if (add_lib_size) df$lib_size <- as.numeric(log(rs + 1))
  df <- cbind(df, as.data.frame(clrM, check.names = FALSE))
  
  # multimedia mediation_data
  exper <- multimedia::mediation_data(
    df, outcomes = "outcome", treatments = "treatment",
    mediator = mediator,
    pretreatments = if (add_lib_size) "lib_size" else NULL
  )
  
  # Fit models
  med_est <- multimedia::lm_model()
  out_est <- if (!add_lib_size) multimedia::lm_model() else
    multimedia::glmnet_model(alpha = alpha_out)
  
  fit <- multimedia::multimedia(
    exper, outcome_estimator = out_est, mediation_estimator = med_est
  ) |> multimedia::estimate(exper)
  
  # Null contrast + FDR at multiple levels
  contrast <- multimedia::null_contrast(
    fit, exper, nullification = "T->M", f = multimedia::indirect_pathwise
  )
  
  q_grid <- c(0.01, 0.02, 0.05, 0.10, 0.20)
  eps <- 1e-4
  
  fdr_list <- lapply(q_grid, function(q) {
    fdr <- multimedia::fdr_summary(contrast, effect = "indirect_pathwise", q_value = q)
    fdr[fdr$source == "real", , drop = FALSE]
  })
  
  # Discoveries at each q level
  discoveries_by_q <- setNames(vector("list", length(q_grid)), paste0("q=", q_grid))
  for (i in seq_along(q_grid)) {
    sel_mediator <- fdr_list[[i]]$mediator[fdr_list[[i]]$keep]
    idx_safe <- match(sel_mediator, safe_names)
    idx_safe <- idx_safe[!is.na(idx_safe)]
    discoveries_by_q[[i]] <- which(keep_taxa)[idx_safe]
  }
  
  # Construct discrete q-values
  q_disc_safe <- rep(1, length(safe_names)); names(q_disc_safe) <- safe_names
  for (i in seq_along(q_grid)) {
    q <- q_grid[i]
    sel_mediator <- fdr_list[[i]]$mediator[fdr_list[[i]]$keep]
    idx_safe <- match(sel_mediator, safe_names)
    idx_safe <- idx_safe[!is.na(idx_safe)]
    to_set <- idx_safe[q_disc_safe[idx_safe] == 1]
    if (length(to_set)) q_disc_safe[to_set] <- max(q - eps, 0)
  }
  
  # Map back to full taxa
  q_disc_full <- rep(1, length(orig_names_full)); names(q_disc_full) <- orig_names_full
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


##############################################################################
# SECTION 5: Data generation from AA template
##############################################################################

generate_data_from_AA <- function(
    num, p, num1_A, num1_B, num2,
    AA_real_sorted,
    beta_treat = log(5), beta_y = 1, d = 0.8,
    sigma_y = 5, seed = 5, row_id = NULL,
    lib_size_lambda = 2e7, lib_size = 1e6,
    dm_conc = 1e6, dirichlet_eps = 1e-10, beta_TX = 5
) {
  # Input validation
  stopifnot(is.finite(num), num >= 1, is.finite(p), p >= 1)
  stopifnot(num1_A >= 0, num1_A <= p, num1_B >= 0, num1_B <= p)
  stopifnot(num2 >= 0, num2 <= min(num1_A, num1_B))
  stopifnot(num1_A + num1_B - num2 <= p)
  stopifnot(d >= 0, d <= 1)
  if (nrow(AA_real_sorted) < num) stop("Template has fewer rows than num")
  if (ncol(AA_real_sorted) < p)   stop("Template has fewer columns than p")
  
  set.seed(seed)
  
  # Extract baseline AA from template
  if (is.null(row_id)) {
    abs_base <- AA_real_sorted[1:num, 1:p, drop = FALSE]
  } else {
    stopifnot(length(row_id) == num)
    abs_base <- AA_real_sorted[row_id, 1:p, drop = FALSE]
  }
  abs_base <- as.matrix(abs_base)
  
  # Construct signal sets: S_treat, S_outcome, overlap
  allowed_taxa <- seq_len(p)
  overlap      <- if (num2 > 0) sample(allowed_taxa, num2) else integer(0)
  remain       <- setdiff(allowed_taxa, overlap)
  treat_only   <- if (num1_A - num2 > 0) sample(remain, num1_A - num2) else integer(0)
  remain2      <- setdiff(remain, treat_only)
  outcome_only <- if (num1_B - num2 > 0) sample(remain2, num1_B - num2) else integer(0)
  
  S_treat   <- sort(c(overlap, treat_only))
  S_outcome <- sort(c(overlap, outcome_only))
  
  # Binary treatment (balanced)
  X <- sample(rep(0:1, length.out = num))
  
  # Spike-in exposure effects on AA
  abs_true <- abs_base
  up_id <- dn_id <- integer(0)
  
  if (length(S_treat) > 0) {
    k_up <- round(d * length(S_treat))
    up_id <- if (k_up > 0) sample(S_treat, k_up) else integer(0)
    dn_id <- setdiff(S_treat, up_id)
    mult <- exp(beta_treat)
    
    rows_t1 <- which(X == 1)
    if (length(rows_t1) > 0 && length(up_id) > 0)
      abs_true[rows_t1, up_id] <- abs_true[rows_t1, up_id, drop = FALSE] * mult
    rows_t0 <- which(X == 0)
    if (length(rows_t0) > 0 && length(dn_id) > 0)
      abs_true[rows_t0, dn_id] <- abs_true[rows_t0, dn_id, drop = FALSE] * mult
  }
  
  # Outcome-side coefficients: magnitude ~ U(0, beta_y), sign controlled by d
  beta_out <- rep(0, p)
  if (length(S_outcome) > 0) {
    m_out <- length(S_outcome)
    mags  <- runif(m_out, 0, beta_y)
    k_pos <- round(d * m_out)
    pos_ids <- if (k_pos > 0) sample(S_outcome, k_pos) else integer(0)
    neg_ids <- setdiff(S_outcome, pos_ids)
    beta_out[pos_ids] <-  mags[match(pos_ids, S_outcome)]
    beta_out[neg_ids] <- -mags[match(neg_ids, S_outcome)]
  }
  
  # Generate outcome: Y = beta_TX * X + sum(beta_k * centered_log_AA_k) + noise
  log_abs <- log(abs_true + 1)
  log_abs_filled <- scale(log_abs, scale = FALSE)
  Y <- as.numeric(log_abs_filled %*% beta_out) + beta_TX * X + rnorm(num, sd = sigma_y)
  
  # Generate sequencing counts via multinomial sampling
  lib_size_vec <- if (!is.null(lib_size_lambda)) {
    rpois(num, lambda = lib_size_lambda)
  } else {
    rep(as.integer(lib_size), num)
  }
  
  source_AA <- pmax(abs_true, 0)
  rs <- rowSums(source_AA)
  
  count <- matrix(0L, nrow = num, ncol = p)
  for (i in seq_len(num)) {
    depth_i <- as.integer(lib_size_vec[i])
    if (!is.finite(depth_i) || depth_i <= 0) { count[i, ] <- 0L; next }
    
    base_prob <- if (is.finite(rs[i]) && rs[i] > 0) source_AA[i, ] / rs[i] else rep(1 / p, p)
    
    # Dirichlet overdispersion
    if (is.finite(dm_conc) && dm_conc > 0) {
      prob_adj <- pmax(base_prob, dirichlet_eps)
      prob_adj <- prob_adj / sum(prob_adj)
      g <- rgamma(p, shape = dm_conc * prob_adj, rate = 1)
      prob <- g / sum(g)
    } else {
      prob <- base_prob
    }
    count[i, ] <- as.integer(rmultinom(1, size = depth_i, prob = prob))
  }
  colnames(count) <- colnames(source_AA)
  rownames(count) <- rownames(source_AA)
  
  list(
    abs_true = abs_true, Y = Y, treat = X, count = count,
    idx1 = sort(overlap), idx2 = setdiff(seq_len(p), sort(overlap)),
    lib_size = lib_size_vec,
    abs_base = abs_base, log_abs_filled = log_abs_filled,
    sets = list(treat = S_treat, outcome = S_outcome, overlap = overlap),
    effects = list(
      treat_up_taxa = up_id, treat_dn_taxa = dn_id,
      beta_treat = beta_treat, mult = exp(beta_treat),
      beta_outcome_vec = beta_out, beta_TX = beta_TX
    ),
    params = list(
      num = num, p = p, num1_A = num1_A, num1_B = num1_B, num2 = num2,
      beta_treat = beta_treat, beta_y = beta_y, d = d,
      sigma_y = sigma_y, lib_size_lambda = lib_size_lambda, seed = seed
    )
  )
}


##############################################################################
# SECTION 6: Simulation runner
##############################################################################

runone_simulation <- function(
    n, p, num1_A, num1_B, num2,
    beta_treat, beta_outcome, d,
    template = "GALAXYMicrobLiver_study",
    template_dir = ".",
    save_dir = ".",
    seed = 1,
    safe = FALSE,
    save_rds = TRUE
) {
  # Input validation
  stopifnot(is.finite(n), n >= 1, is.finite(p), p >= 1)
  stopifnot(num1_A >= 0, num1_A <= p, num1_B >= 0, num1_B <= p)
  stopifnot(num2 <= min(num1_A, num1_B))
  stopifnot(num1_A + num1_B - num2 <= p)
  stopifnot(d >= 0, d <= 1)
  
  # Load template
  template <- match.arg(template, c("GALAXYMicrobLiver_study"))
  f <- file.path(template_dir, paste0(template, ".RData"))
  if (!file.exists(f)) stop("Template file not found: ", f)
  
  e <- new.env(parent = emptyenv())
  load(f, envir = e)
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
  p_use     <- ncol(count_m)
  idx_true  <- sort(sim_data$sets$overlap)
  
  # ── Helper: align p-value vector to length p ──
  p_to_lenp <- function(pvec, method = "UNKNOWN") {
    out <- rep(NA_real_, p_use)
    if (is.null(pvec)) return(out)
    pvec <- as.numeric(pvec)
    if (length(pvec) == p_use) return(pvec)
    
    warning(sprintf("[%s] length(pvec)=%d != p=%d, aligning.", method, length(pvec), p_use))
    if (!is.null(names(pvec)) && !is.null(colnames(count_m))) {
      mm <- match(colnames(count_m), names(pvec))
      out[!is.na(mm)] <- as.numeric(pvec[mm[!is.na(mm)]])
      return(out)
    }
    L <- min(p_use, length(pvec))
    out[seq_len(L)] <- pvec[seq_len(L)]
    out
  }
  
  # ── Run methods ──
  methods <- c("CAMRA", "HIMA", "LDM", "MarZIC", "multimedia", "CRAmed")
  p_mat <- matrix(NA_real_, nrow = length(methods), ncol = p_use,
                  dimnames = list(methods, colnames(count_m)))
  runtime <- setNames(rep(NA_real_, length(methods)), methods)
  
  safe_eval <- function(expr) {
    if (!isTRUE(safe)) return(eval(expr))
    tryCatch(eval(expr), error = function(e) NULL)
  }
  
  run_one_method <- function(method_name, expr) {
    res <- safe_eval(expr)
    if (is.null(res)) {
      p_mat[method_name, ] <<- rep(NA_real_, p_use)
      runtime[method_name] <<- NA_real_
      return(invisible(NULL))
    }
    p_mat[method_name, ] <<- p_to_lenp(res[[2]], method = method_name)
    rt <- suppressWarnings(as.numeric(res[[3]]))
    runtime[method_name] <<- if (length(rt) == 1 && is.finite(rt)) rt else NA_real_
    invisible(NULL)
  }
  
  run_one_method("CAMRA-HDMT",  quote(CAMRA(count_m, treat_cov, y, CClasso = FALSE, method = "hdmt")))
  run_one_method("CAMRA-SBMH",  quote(CAMRA(count_m, treat_cov, y, CClasso = FALSE, method = "sbmh")))
  run_one_method("HIMA",       quote(HIMA_micro_sim1(count_m, treat_cov, y)))
  run_one_method("LDM",        quote(ldm_sim(count_m, treat_cov, y)))
  run_one_method("MarZIC",     quote(Mar_sim(count_m, treat_cov, y)))
  run_one_method("multimedia", quote(multimedia_sim(count_m, treat_cov, y)))
  run_one_method("CRAmed",     quote(CRAmed_sim(count_m, treat_cov, y)))
  
  out_list <- list(p_mat = p_mat, idx_true = idx_true, runtime = runtime)
  
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

kk <- runone_simulation(
  n = 200, 
  p = 200, 
  num1_A = 10, 
  num1_B = 10, 
  num2 = 7,
  beta_treat = log(5), 
  beta_outcome = 1, 
  d = 0.5,
  template = "GALAXYMicrobLiver_study",
  template_dir = ".",
  save_dir = ".",
  seed = 1
)
kk
