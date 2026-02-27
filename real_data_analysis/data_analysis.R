load("real_data.RData")
pkgs <- c(
  "glmnet", "pscl", "plyr", "hdi", "compositions",
  "Iso", "cp4p", "HDMT", "tidyverse", "LDM",
  "harmonicmeanp", "precrec", "multimedia",
  "ccmm", "MarZIC", "HIMA", "PALM", "gtools",
  "MultiMed", "permute", "vegan", "matrixStats", "energy"
)
for (pkg in pkgs) library(pkg, character.only = TRUE)
source("MODIMA.R")
source("MedTest.R")

############################============1.set all method funciton
####1.1 CAMRA
{
  ############computation for correlation
  {
    
    {
      
      fastCCLasso <- function(xx, isCnt = FALSE, pseudo = 0.5, k_cv = 3, 
                              lam_min_ratio = 1E-4, k_max = 20, n_boot=100, aa=NULL, bb=NULL) {
        n = nrow(xx);
        p = ncol(xx);
        if(isCnt) {
          xx = xx + pseudo;
          xx = xx / rowSums(xx);
        };
        xx2 = log(xx) - rowMeans(log(xx));
        vxx2 = stats::var(xx2);
        
        if(is.null(aa)){
          aa = rep(1, p);
        }
        if(is.null(bb)){
          bb = 1 / diag(vxx2);
        }		
        
        #-Golden section method for the selection of lambda (log10 scale)
        xx = vxx2 * (aa * rep(bb, each = p) + bb * rep(aa, each = p))/2;
        diag(xx) = 0;
        lam_max = max(abs(xx));
        lam_int2 = log10(lam_max * c(lam_min_ratio, 1));
        a1 = lam_int2[1]; 
        b1 = lam_int2[2];
        
        #-Store lambda and corresponding cross validation's loss
        lams = NULL; 
        fvals = NULL;
        #-Two trial points in first 
        a2 = a1 + 0.382 * (b1 - a1); 
        b2 = a1 + 0.618 * (b1 - a1);
        fb2 = cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2, 
                            aa = aa, bb = bb);
        lams = c(lams, b2); 
        fvals = c(fvals, fb2);
        fa2 = cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2, 
                            aa = aa, bb = bb);
        lams = c(lams, a2); 
        fvals = c(fvals, fa2);
        #Error tolerance for convergence
        err_lam2 = 1e-1 * max(1, lam_int2);
        err_fval = 1e-4;
        err = b1 - a1;
        k = 0;
        
        while(err > err_lam2 && k < k_max) {
          fval_max = max(fa2, fb2);
          if(fa2 > fb2) {
            a1 = a2;
            a2 = b2;
            fa2 = fb2;
            b2 = a1 + 0.618 * (b1 - a1);
            fb2 = cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2, 
                                aa = aa, bb = bb);
            lams = c(lams, b2); 
            fvals = c(fvals, fb2);
          } else {
            b1 = b2;
            b2 = a2;
            fb2 = fa2;
            a2 = a1 + 0.382 * (b1 - a1);
            fa2 = cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2, 
                                aa = aa, bb = bb);
            lams = c(lams, a2);
            fvals = c(fvals, fa2);
          };
          fval_min = min(fa2, fb2);
          k = k + 1;
          err = b1 - a1;
          if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
            break;
          };
        };
        info_cv = list(lams = lams, fvals = fvals, k = k + 2, 
                       lam_int = 10^c(a1, b1));
        #if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
        #	cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
        #		10^b1, "])\n", sep = "");
        #};
        lambda = 10^((a2 + b2)/2);
        fit_res = fastcclasso_sub(lambda = lambda, SS2 = vxx2, aa = aa, bb = bb);
        
        sigma_mod <- boot_fastCCLasso(xx2=xx2,sigma_hat= fit_res$sigma,
                                      lambda=lambda,aa=aa,bb=bb,
                                      n_boot = n_boot, max_iter=200,
                                      stop_eps=1e-6)
        
        return(list(rho=sigma_mod$cor_w,
                    cov_diag=sigma_mod$var_w,
                    lambda_best=lambda,
                    info_cv=info_cv,
                    p_vals=sigma_mod$p_vals))
      };
      #---------------------------------------
      #cross validation's loss of fastcclasso for single lambda
      cvfastCCLasso <- function(lambda, k_cv, xx2, aa, bb) {
        n = nrow(xx2);
        p = ncol(xx2);
        n_b = floor(n / k_cv);
        cv.loss = 0;
        for(k in 1:k_cv) {
          ite = (n_b * (k - 1) + 1):(n_b * k);
          vxx2te = stats::var(xx2[ite, ]);
          vxx2tr = stats::var(xx2[-ite, ]);
          out = fastcclasso_sub(lambda = lambda, SS2 = vxx2tr, aa = aa, bb = bb);
          mm = out$sigma - out$ww - rep(out$ww, each = p) - vxx2te;
          cv.loss = cv.loss + mean(mm^2 * aa * rep(bb, each = p));
        };
        return(cv.loss);
      };
      #---------------------------------------
      #fastcclasso for single lambda
      fastcclasso_sub <- function(lambda, SS2, aa, bb, k_max = 200, x_tol = 1E-4) {
        p = ncol(SS2);
        cc = 1 / (aa * sum(bb) + bb * sum(aa));
        aa2 = aa * cc;
        bb2 = bb * cc;
        cab1 = 1 + sum(aa * bb2);
        caa = sum(aa * aa2);
        cbb = sum(bb * bb2);
        aabb = aa * rep(bb, each = p) + bb * rep(aa, each = p);
        lambda2 = 2 * lambda / aabb;
        ss2 = rowSums(SS2 * aabb);
        sigma = SS2;
        ww = colMeans(sigma) - mean(sigma)/2;
        k = 0;
        err = 1;
        while(err > x_tol && k < k_max) {
          # Update ww
          xx = rowSums(sigma * aabb) - ss2;
          ax1 = sum(aa2 * xx);
          bx1 = sum(bb2 * xx);
          ww2 = xx * cc + (aa2 * (cbb * ax1 - cab1 * bx1) + bb2 * (caa * bx1 -
                                                                     cab1 * ax1)) / (cab1^2 - caa * cbb);
          # Update sigma
          sigma2 = SS2 + ww2 + rep(ww2, each = p);
          oo = diag(sigma2);
          sigma2 = (sigma2 > lambda2) * (sigma2 - lambda2) + 
            (sigma2 < - lambda2) * (sigma2 + lambda2);
          diag(sigma2) = oo;
          # Check convergence
          err = max(abs(sigma2 - sigma)/(abs(sigma) + 1)); 
          k = k + 1;
          sigma = sigma2;
        };
        #if(k >= k_max) {
        #	cat("WARNING of fastcclasso_sub:\n", "\tMaximum Iteration:", k_max, 
        #		"&& Relative error:", err, "!\n");
        #};
        return(list(sigma = sigma, ww = ww2));
      };
      #---------------------------------------
      boot_fastCCLasso <- function(xx2,sigma_hat,lambda,aa,bb, 
                                   n_boot = 100,
                                   max_iter=200, 
                                   stop_eps=1e-6) {
        n <- nrow(xx2);
        p <- ncol(xx2);
        
        # Store the result of bootstrap
        cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
        vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
        cors_mat <- matrix(0, p, p);
        ind_low <- lower.tri(cors_mat);
        
        # Bootstrap procedure
        sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
                           ncol = n_boot);
        for(k in 1:n_boot) {
          ind_samp <- sam_boot[, k];
          S_samp <- stats::var(xx2[ind_samp,])
          cov_est<- fastcclasso_sub(lambda, SS2=S_samp, 
                                    aa=aa, bb=bb,
                                    k_max = 200, x_tol = stop_eps);
          vars_boot[, k] <- diag(cov_est$sigma);
          Is <- 1 / sqrt(vars_boot[, k]);
          cor_est <- Is * cov_est$sigma * rep(Is, each = p);
          cors_boot[, k] <- cor_est[ind_low];
        }
        
        vars_boot[, n_boot + 1] <- diag(sigma_hat);
        Is <- 1 / sqrt(vars_boot[, n_boot + 1]);
        cor_est <- Is * sigma_hat * rep(Is, each = p);  
        cors_boot[, n_boot + 1] <- cor_est[ind_low]; 
        
        # Variance estimation via bootstrap
        vars2 <- rowMeans(vars_boot);
        cors2mod <- rowMeans(cors_boot);
        cors2_mat <- diag(p);
        cors2_mat[ind_low] <- cors2mod;
        cors2_mat <- t(cors2_mat);
        cors2_mat[ind_low] <- cors2mod;
        # P-values with null distribution of correlation estimations of absolute data
        p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
        p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
        pval_mat <- diag(p);
        pval_mat[ind_low] <- p_vals;
        pval_mat <- t(pval_mat);
        pval_mat[ind_low] <- p_vals;
        
        return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
      }
      #-------------------------------------------------------------------------------
      # call all methods
      Callallmethods <- function(method,xMat,cv_k, lambda_min_ratio=1e-4,
                                 Edge_eps=1e-4){
        # xMat: compositional data
        p <- dim(xMat)[2];
        S <- var(log(xMat) - rowMeans(log(xMat)));
        lambda_max <- max(max(S - diag(p)), -min(S - diag(p)));
        lambda_min <- lambda_min_ratio * lambda_max;
        lambda_int <- c(lambda_min,lambda_max);
        
        # method  
        if(method=="fastCCLasso"){
          begin_time <- proc.time();
          result <- fastCCLasso(xx=xMat,lam_min_ratio = lambda_min_ratio, 
                                k_cv = cv_k, k_max = 20);
          end_time <- proc.time();
          result_cor <- result$rho;
        }else if(method=="SparCC"){
          begin_time <- proc.time();
          result <- compute_corr_mod(fracs=xMat, iter=10, th=0.1);
          end_time <- proc.time();
          result_cor <- result$Cor.mat;
        }else if(method=="CCLasso"){    
          begin_time <- proc.time();
          result <- cclasso(xMat, counts = FALSE, pseudo = 0.5, k_cv = cv_k, 
                            lam_int = lambda_int, k_max = 20);
          end_time <- proc.time();	
          result_cor <- result$cor_w;
        }else if(method=="COAT"){
          begin_time <- proc.time();
          result <- coat(xMat, nFoler = cv_k, soft = 1);
          end_time <- proc.time();
          result_cor <- result$corr;
        }else{
          return(message("This method is not exist, please check it! "));
        } 
        result_cor[abs(result_cor)< Edge_eps] <- 0 ;
        
        return(list( runtime=as.numeric((end_time-begin_time)[3]),
                     est_lower=result_cor[lower.tri(result_cor)] ,
                     cor_est=result_cor));
      }
      
    }
    SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
      # dimension for w (latent variables)
      p <- ncol(x);
      n <- nrow(x);
      # posterior distribution (alpha)
      x <- x + 1;
      # store generate data
      y <- matrix(0, n, p);
      # store covariance/correlation matrix
      cov.w <- cor.w <- matrix(0, p, p);
      indLow <- lower.tri(cov.w, diag = T);
      # store covariance/correlation for several posterior samples
      covs <- cors <- matrix(0, p * (p + 1) / 2, imax);
      for(i in 1:imax) {
        # generate fractions from posterior distribution
        y <- t(apply(x, 1, function(x) 
          gtools::rdirichlet(n = 1, alpha = x)));
        # estimate covariance/correlation
        cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin);
        # store variance/correlation only low triangle 
        covs[, i] <- cov_cor$cov.w[indLow];
        cors[, i] <- cov_cor$cor.w[indLow];
      }
      # calculate median for several posterior samples
      cov.w[indLow] <- apply(covs, 1, median); 
      cor.w[indLow] <- apply(cors, 1, median);
      #
      cov.w <- cov.w + t(cov.w);
      diag(cov.w) <- diag(cov.w) / 2;
      cor.w <- cor.w + t(cor.w);
      diag(cor.w) <- 1;
      #
      return(list(cov.w = cov.w, cor.w = cor.w));
    }
    #-------------------------------------------------------------------------------
    # SparCC for fractions known
    #   function: SparCC.frac
    #   input:
    #          x ------ nxp fraction data matrix, row is sample, col is variable
    #       kmax ------ max iteration steps for SparCC. default is 10
    #      alpha ------ the threshold for strong correlation. default is 0.1
    #       Vmin ------ minimal variance if negative variance appears. default is 1e-4
    #   output: a list structure
    #      cov.w ------ covariance estimation
    #      cor.w ------ correlation estimation
    SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
      # Log transformation
      x <- log(x);
      p <- ncol(x);
      # T0 = var(log(xi/xj)) variation matrix
      TT <- stats::var(x);
      T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT;
      # Variance and correlation coefficients for Basic SparCC  
      rowT0 <- rowSums(T0);
      var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2);
      var.w[var.w < Vmin] <- Vmin;
      #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
      #  sqrt(outer(var.w, var.w, "*")) / 2;
      Is <- sqrt(1/var.w);
      cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5;
      # Truncated correlation in [-1, 1]
      cor.w[cor.w <= - 1] <- - 1; 
      cor.w[cor.w >= 1] <- 1;
      # Left matrix of estimation equation
      Lmat <- diag(rep(p - 2, p)) + 1; 
      # Remove pairs
      rp <- NULL;
      # Left components
      cp <- rep(TRUE, p);
      # Do loops until max iteration or only 3 components left
      k <- 0;  
      while(k < kmax && sum(cp) > 3) {
        # Left T0 = var(log(xi/xj)) after removing pairs
        T02 <- T0;
        # Store current correlation to find the strongest pair
        curr_cor.w <- cor.w;
        # Remove diagonal
        diag(curr_cor.w) <- 0;
        # Remove removed pairs
        if(!is.null(rp)) {
          curr_cor.w[rp] <- 0;
        }
        # Find the strongest pair in vector form
        n_rp <- which.max(abs(curr_cor.w));
        # Remove the pair if geater than alpha
        if(abs(curr_cor.w[n_rp]) >= alpha) {
          # Which pair in matrix form
          t_id <- c(arrayInd(n_rp, .dim = c(p, p)));
          Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1;
          # Update remove pairs
          n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp);
          rp <- c(rp, n_rp);
          # Update T02
          T02[rp] <- 0;
          # Which component left
          cp <- (diag(Lmat) > 0);
          # Update variance and truncated lower by Vmin
          var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]));
          var.w[var.w <= Vmin] <- Vmin;
          # Update correlation matrix and truncated by [-1, 1]
          #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
          #  sqrt(outer(var.w, var.w, "*")) / 2;    
          Is <- sqrt(1/var.w);
          cor.w <- (var.w + rep(var.w, each = p) - T0) * 
            Is * rep(Is, each = p) * 0.5;
          # Truncated correlation in [-1, 1]
          cor.w[cor.w <= - 1] <- - 1;
          cor.w[cor.w >= 1] <- 1;
        }
        else {
          break;
        }
        # 
        k <- k + 1;
      }
      # Covariance
      Is <- sqrt(var.w);
      cov.w <- cor.w * Is * rep(Is, each = p);
      #
      return(list(cov.w = cov.w, cor.w = cor.w));
    }
    #-------------------------------------------------------------------------------
    
  }
  
  
  
  ###############--------internal function for CAMRA 
  {
    
    
    #####CAMRA
    {
      recover_l_PALM <- function(count_m, treat_cov, cov_ad = NULL,
                                 prev.filter = 0, eps_p = 1e-10) {
        count_m <- as.matrix(count_m)
        storage.mode(count_m) <- "numeric"
        n <- nrow(count_m)
        p <- ncol(count_m)
        
        # ---- keep original taxa order ----
        orig_taxa <- colnames(count_m)
        orig_taxa <- paste0("O", seq_len(p))
        colnames(count_m) <- orig_taxa
        
        # ---- treat ----
        stopifnot(length(treat_cov) == n)
        treat_cov <- matrix(as.numeric(treat_cov), ncol = 1)
        colnames(treat_cov) <- "treat"
        
        # ---- rownames (optional, but keep consistent) ----
        rn <- paste0("T", seq_len(n))
        rownames(count_m) <- rn
        rownames(treat_cov) <- rn
        
        if (!is.null(cov_ad)) {
          cov_ad <- data.frame(cov_ad)
          stopifnot(nrow(cov_ad) == n)
          rownames(cov_ad) <- rn
          colnames(cov_ad) <- paste0("Cov", seq_len(ncol(cov_ad)))
        }
        
        # ---- run PALM once ----
        result1 <- PALM::palm(
          rel.abd            = count_m,
          covariate.interest = treat_cov,
          covariate.adjust   = cov_ad,
          prev.filter        = prev.filter
        )
        
        r1 <- result1$treat
        
        # ---- initialize full-length outputs (filtered taxa will remain defaults) ----
        p_full    <- rep(1, p)  
        z_full    <- rep(0, p)  
        beta_full <- rep(0, p)  
        
        names(p_full) <- names(z_full) <- names(beta_full) <- orig_taxa
        
        # ---- align kept features back to original order ----
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
        
        
        return(list(
          p      = p_full,
          z      = z_full,
          beta_l = beta_full,
          feature_kept = feat
        ))
      }
      
      recover_r <- function(count_m, treat_cov, y, sudo = 0.5, cov_ad = NULL,
                            CClasso = FALSE, cov_true = NULL) {
        
        ## ===== 1) compositional transform + covariance estimate =====
        logdata  <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
        por_data <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
        
        CClasso_core <- function(count_m) {
          res_cov <- fastCCLasso(count_m, isCnt = TRUE)
          diag(sqrt(res_cov$cov_diag)) %*% res_cov$rho %*% diag(sqrt(res_cov$cov_diag))
        }
        
        if (CClasso) {
          est_cov <- CClasso_core(count_matrix)
        } else {
          res_l  <- SparCC.count(count_matrix)
          est_cov <- res_l$cov.w
        }
        
        if (!is.null(cov_true)) {
          est_cov <- cov_true
        }
        
        ## ===== 2) build Z_ilr 
        p <- ncol(count_matrix)
        n <- nrow(count_matrix)
        
        ilr_basis <- compositions::ilrBase(por_data)
        lasso_data_ilr <- as.matrix(compositions::ilr(por_data))
        
        R2 <- ilr_basis
        
        # Z_ilr: n x p  (as in your code)
        Z_ilr <- lasso_data_ilr %*% solve(t(R2) %*% est_cov %*% R2) %*% t(R2) %*% est_cov
        
        
        ## ===== 3) FWL residualization: treat (+ cov_ad) unpenalized =====
        # --- build unpenalized design matrix Z0 with intercept ---
        # allow treat_cov to be vector / matrix / data.frame
        if (is.null(dim(treat_cov))) {
          treat_df <- data.frame(treat = as.numeric(treat_cov))
        } else {
          treat_df <- as.data.frame(treat_cov)
          if (nrow(treat_df) != n) stop("The number of rows in treat_cov is inconsistent with count_matrix.")
        }
        
        if (!is.null(cov_ad)) {
          cov_df <- as.data.frame(cov_ad)
          if (nrow(cov_df) != n) stop("The number of rows in cov_ad is inconsistent with count_matrix.")
          Zdf <- cbind(treat_df, cov_df)
        } else {
          Zdf <- treat_df
        }
        
        
        Z0 <- model.matrix(~ ., data = Zdf)
        
        # --- QR-based residualization (stable) ---
        qrZ <- qr(Z0)
        
        y_hat   <- qr.fitted(qrZ, y)
        y_tilde <- as.numeric(y - y_hat)
        
        X_hat   <- qr.fitted(qrZ, Z_ilr)
        X_tilde <- Z_ilr - X_hat
        
        ## ===== 4) ridge.proj only on penalized part (residualized high-dim) =====
        outRidge <- hdi::ridge.proj(x = X_tilde, y = y_tilde)
        
        all_p  <- outRidge$pval
        beta_r <- outRidge$bhat
        z      <- qnorm(1 - all_p / 2) * sign(beta_r)
        
        beta_r <- as.vector(beta_r)
        z <- as.vector(z)
        return(list(
          p      = all_p,     
          z      = z,
          beta_r = beta_r,    
          y_tilde = y_tilde,
          X_tilde = X_tilde,
          Z0      = Z0,
          X_doubel = Z_ilr
        ))
      }
      
      pre_filter_fun <- function(count_matrix, treat_cov, y,
                                 const = 2,
                                 seed = 42,
                                 sudo = 0.5,
                                 cov_ad = NULL,
                                 adaptive_L = FALSE) {
        
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
        
        
        logdata <- log((count_matrix + sudo) / rowSums(count_matrix + sudo))
        logdata[logdata < (-10)] <- (-10)
        
        # ---- covariance estimation (SparCC) ----
        res_l  <- SparCC.count(count_matrix)
        est_cov <- res_l$cov.w
        
        por_data <- (count_matrix + sudo) / rowSums(count_matrix + sudo)
        
        
        ilr_basis <- compositions::ilrBase(por_data)
        R2 <- ilr_basis
        
        
        Z_ilr <- (logdata %*% R2) %*%
          solve(t(R2) %*% est_cov %*% R2) %*%
          t(R2) %*% est_cov
        
        
        treat_vec <- as.numeric(treat_cov)
        Z0 <- cbind(treat = treat_vec)
        if (!is.null(cov_ad)) {
          Z0 <- cbind(Z0, cov_ad)
          colnames(Z0) <- make.names(colnames(Z0), unique = TRUE)
        }
        
        
        X <- cbind(Z_ilr, Z0)
        
        pZ <- ncol(Z_ilr)   # mediators (penalized)
        p0 <- ncol(Z0)      # unpenalized: treat + cov
        pf <- c(rep(1, pZ), rep(0, p0))
        
        if (!requireNamespace("glmnet", quietly = TRUE)) {
          stop("Package 'glmnet' is required for pre_filter_fun().")
        }
        
        
        
        ##############Use flexiable filter
        if (isTRUE(adaptive_L)) {
          cvfit <- glmnet::cv.glmnet(
            x = X, y = y,
            alpha = 1,
            penalty.factor = pf,
            nfolds = 5,
            type.measure = "mse",
            standardize = TRUE
          )
          
          
          b <- as.matrix(coef(cvfit, s = "lambda.min"))  
          
          beta_all <- as.numeric(b)[-1]  
          beta_Z   <- beta_all[1:pZ]     
          
          selection_set <- which(beta_Z != 0)
          
          if (length(selection_set) == 0) {
            ord <- order(abs(beta_Z), decreasing = TRUE)
            selection_set <- ord[1]
          }
          
          return(sort(unique(as.integer(selection_set))))
        }
        
        # ============================
        # Fixed-K: your original rule
        # ============================
        K_raw <- floor(const * n / log(max(n, 3)))
        K     <- max(1L, min(pZ, K_raw))
        
        fit <- glmnet::glmnet(
          x = X, y = y,
          alpha = 1,
          penalty.factor = pf,
          dfmax = min(K + p0, pZ + p0),  # K mediators + all unpenalized vars
          nlambda = 500, lambda.min.ratio = 1e-6,
          standardize = TRUE
        )
        
        B  <- as.matrix(fit$beta)  # (pZ+p0) x L (no intercept)
        dfZ <- colSums(B[1:pZ, , drop = FALSE] != 0)
        
        idx <- which(dfZ >= K)[1]
        if (is.na(idx)) idx <- ncol(B)
        
        beta_Z_full <- as.numeric(B[1:pZ, idx])
        ord  <- order(abs(beta_Z_full), decreasing = TRUE)
        keep <- ord[seq_len(min(K, length(ord)))]
        
        sort(unique(as.integer(keep)))
      }
      
      
      
      p_mediation_maxp <- function(p_alpha, p_beta,
                                   pi_alpha0 = NULL, pi_beta0 = NULL,
                                   pi_method = c("JC","cp4p"),
                                   weight_method = c("maxp","product","indenp")) {
        stopifnot(length(p_alpha) == length(p_beta))
        pi_method     <- match.arg(pi_method)
        weight_method <- match.arg(weight_method)
        
        
        ####
        {
          mix_weights_product <- function(pi_alpha0, pi_beta0) {
            eps <- 1e-8
            pa <- min(max(pi_alpha0, eps), 1 - 1e-6)
            pb <- min(max(pi_beta0,  eps), 1 - 1e-6)
            pi0 <- 1 - (1 - pa) * (1 - pb)        
            pi0 <- max(pi0, 1e-6)                 
            
            w00 <- (pa * pb) / pi0
            w10 <- ((1 - pa) * pb) / pi0
            w01 <- (pa * (1 - pb)) / pi0
            c(w00 = w00, w10 = w10, w01 = w01, pi0 = pi0)
          }
          
          #############
          from_miMed_estpi0 <- function(z)
          {
            
            xi = c(0:100)/100
            tmax = sqrt(log(length(z)))
            tt = seq(0, tmax, 0.05)
            epsest = NULL
            for (j in 1:length(tt)) {
              t = tt[j]
              f = t * xi
              f = exp(f^2/2)
              w = (1 - abs(xi))
              co = 0 * xi
              for (i in 1:101) {
                co[i] = mean(cos(t * xi[i] * z))
              }
              epshat = sum(w * f * co)/sum(w)
              epsest = c(epsest, epshat)
            }
            tmp = min(epsest)
            if (tmp > 1) 
              tmp = 1
            return(tmp)
            
          }
          
          mix_weights_maxp <- function(p_alpha, p_beta,
                                       pi_alpha0 = NULL, pi_beta0 = NULL,
                                       pi_method = c("cp4p","JC")  ##default cp4p
          ) {
            pi_method <- match.arg(pi_method)
            
            
            if (is.null(pi_alpha0) || is.null(pi_beta0)) {
              if (pi_method == "cp4p") {
                pa <- cp4p::estim.pi0(p_alpha); pb <- cp4p::estim.pi0(p_beta)
                grab <- function(x) if (!is.null(x$pi0)) as.numeric(x$pi0)
                else mean(unlist(x), na.rm = TRUE)
                pi_alpha0 <- grab(pa); pi_beta0 <- grab(pb)
              } else {
                z1 <- qnorm(1 - p_alpha); z2 <- qnorm(1 - p_beta)
                pi_alpha0 <- from_miMed_estpi0(z1)
                pi_beta0  <- from_miMed_estpi0(z2)
              }
            }
            
            
            p_max <- pmax(p_alpha, p_beta)
            if (pi_method == "cp4p") {
              pi0_hat <- {
                obj <- cp4p::estim.pi0(p_max)
                if (!is.null(obj$pi0)) as.numeric(obj$pi0) else mean(unlist(obj))
              }
            } else {
              zmax   <- qnorm(1 - p_max)
              pi0_hat <- from_miMed_estpi0(zmax)
            }
            
            
            clip01 <- function(x) min(max(x, 1e-6), 1 - 1e-6)
            pi_alpha0 <- clip01(pi_alpha0); pi_beta0 <- clip01(pi_beta0); pi0_hat <- clip01(pi0_hat)
            
            w00 <- (pi_alpha0 + pi_beta0 - pi0_hat) / pi0_hat
            w10 <- (pi0_hat - pi_alpha0) / pi0_hat
            w01 <- (pi0_hat - pi_beta0)  / pi0_hat
            
            w <- pmax(c(w00, w10, w01), 0); w <- w / sum(w)
            names(w) <- c("w00","w10","w01")
            w
          }
        }
        ## ---------- 1. Clean p for pi0 estimation -
        p_alpha_pi0 <- p_alpha
        p_beta_pi0  <- p_beta
        
        bad_a <- !is.finite(p_alpha_pi0)
        bad_b <- !is.finite(p_beta_pi0)
        p_alpha_pi0[bad_a] <- 1
        p_beta_pi0 [bad_b] <- 1
        
        p_alpha_pi0 <- pmin(pmax(p_alpha_pi0, 0), 1)
        p_beta_pi0  <- pmin(pmax(p_beta_pi0,  0), 1)
        
        ## ---------- 2. Estimate marginal pi0 if not given ----------
        if (is.null(pi_alpha0) || is.null(pi_beta0)) {
          if (pi_method == "cp4p") {
            pa <- cp4p::estim.pi0(p_alpha_pi0)
            pb <- cp4p::estim.pi0(p_beta_pi0)
            grab <- function(x) {
              if (!is.null(x$pi0)) as.numeric(x$pi0) else mean(unlist(x), na.rm = TRUE)
            }
            pi_alpha0 <- grab(pa)
            pi_beta0  <- grab(pb)
          } else {  # "JC"
            z1 <- qnorm(1 - p_alpha_pi0)
            z2 <- qnorm(1 - p_beta_pi0)
            pi_alpha0 <- from_miMed_estpi0(z1)
            pi_beta0  <- from_miMed_estpi0(z2)
          }
        }
        
        eps <- 1e-8
        pi_alpha0 <- min(max(pi_alpha0, eps), 1 - eps)
        pi_beta0  <- min(max(pi_beta0,  eps), 1 - eps)
        
        ## ---------- 3. Use only finite pairs to compute p_mix ----------
        keep <- is.finite(p_alpha) & is.finite(p_beta)
        out  <- rep(NA_real_, length(p_alpha))
        if (!any(keep)) return(out)
        
        p_a <- p_alpha[keep]
        p_b <- p_beta[keep]
        p_a <- pmin(pmax(p_a, eps), 1 - eps)
        p_b <- pmin(pmax(p_b, eps), 1 - eps)
        
        ## ---------- 4.  Weights ---------- w00, w10, w01 ----------
        if (weight_method == "maxp") {
          w <- mix_weights_maxp(p_a, p_b,
                                pi_alpha0 = pi_alpha0,
                                pi_beta0  = pi_beta0,
                                pi_method = pi_method)
        } else if (weight_method == "product") {
          
          w_raw <- mix_weights_product(pi_alpha0, pi_beta0) 
          w_vec <- pmax(w_raw[c("w00","w10","w01")], 0)
          w     <- w_vec / sum(w_vec)
          names(w) <- c("w00","w10","w01")
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
        
        ## ---------- 5. Grenander  F1α, F1β ----------
        estimate_F1_grenander <- function(p, pi0, eps = 1e-8) {
          p <- p[is.finite(p)]
          p <- pmin(pmax(p, 0), 1)
          n <- length(p); stopifnot(n > 0)
          pi0 <- min(max(pi0, 1e-6), 1 - 1e-6)
          
          x <- sort(unique(c(0, sort(p), 1)))
          Fn <- ecdf(p); y <- Fn(x)
          dx <- diff(x); keep <- dx > eps
          xL <- x[-length(x)][keep]; xR <- x[-1][keep]
          yL <- y[-length(y)][keep]; yR <- y[-1][keep]; dx <- xR - xL
          
          s <- (yR - yL) / dx                
          s_hat <- -Iso::pava(-s, w = dx)      
          f1_hat <- pmax((s_hat - pi0) / (1 - pi0), 0)
          area <- sum(f1_hat * dx)
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
        
        ## ---------- 6. mixture null
        t     <- pmax(p_a, p_b)
        p_mix <- w00 * (t^2) + w10 * (t * F1a(t)) + w01 * (t * F1b(t))
        p_mix <- pmin(pmax(p_mix, 0), 1)
        
        out[keep] <- p_mix
        out
      }
      
      p_mediation_hdmt_fdr <- function(p_alpha, p_beta ,exact_p =1) {
        stopifnot(length(p_alpha) == length(p_beta))
        n   <- length(p_alpha)
        out <- rep(NA_real_, n)
        
        keep <- is.finite(p_alpha) & is.finite(p_beta)
        if (!any(keep)) return(out)
        
        pa <- pmin(pmax(p_alpha[keep], 0), 1)
        pb <- pmin(pmax(p_beta [keep], 0), 1)
        input <- cbind(pa, pb)
        
        nullprop <- HDMT::null_estimation(input)  
        fdr <- HDMT::fdr_est(nullprop$alpha00, nullprop$alpha01, nullprop$alpha10,
                             nullprop$alpha1, nullprop$alpha2,
                             input_pvalues = input, exact = exact_p )  
        
        out[keep] <- fdr
        
        out
      }
      
      
      
      
    }
    
  }
  
  
  CAMRA <- function(mediators,
                    treatment,
                    outcome,
                    confounders= NULL ,
                    pseudo=0.5,
                    FDR_level =0.05,
                    hdmt.exact = 0,
                    screen = FALSE,        #### Prefilter function
                    const =2,              #### The number of filter bacteria, const * n/log(n)
                    CClasso = FALSE,       #### if use CClasso to compute correlation,default is Rsparcc
                    seed=42
  )
  {
    set.seed(seed)
    t0 <- proc.time()[["elapsed"]]
    
    ## Default keep all taxa
    select_otu <- c(1:ncol(mediators))
    ## Optional pre-screening
    
    if (isTRUE(screen)) {
      select_otu <- pre_filter_fun(
        count_matrix = mediators,
        treat_cov    = treatment,
        y            = outcome,
        const        = const ,
        seed         = seed,
        sudo         = pseudo,
        cov_ad       = confounders
      )
    }
    
    
    ## Step A: exposure->taxon
    res1 <- recover_l_PALM(count_m=mediators,
                           treat_cov = treatment,
                           cov_ad = confounders)
    ## Step B: taxon->outcome 
    res2 <- recover_r(count_m=mediators,
                      treat_cov = treatment,
                      y =outcome,
                      cov_ad = confounders,
                      CClasso = CClasso,
                      sudo=pseudo)
    
    
    
    p1 <- res1$p
    p2 <- res2$p
    z1 <- res1$z
    z2 <- as.vector(res2$z)
    
    p_matrix <- cbind(p1,p2)
    rownames(p_matrix) <- colnames(mediators)
    rawp.perm <- p_mediation_maxp(p1,p2,pi_method="cp4p",weight_method = "product")   
    p_vec <- p.adjust(rawp.perm, method = "BH")
    
    rawp.perm.rm = na.omit(rawp.perm)
    L = length(rawp.perm.rm)
    rawp.perm.rm[rawp.perm.rm< 1e-10] <- 1e-10
    
    
    p_vec_f <- NULL
    p_vec_all <- p_vec
    
    if(screen)
    {
      p_vec_f <- p.adjust(rawp.perm[select_otu], method = "BH")
      p_vec_all <- p_vec
      p_vec_all[select_otu] <- p_vec_f 
      p_vec_all[-select_otu] <- 1
    }
    
    
    
    # user-facing option:
    #   hdmt.exact = 0 (default) means try exact_p = 0 first
    #   hdmt.exact = 1 means try exact_p = 1 first
    stopifnot(hdmt.exact %in% c(0, 1))
    
    exact_first  <- as.integer(hdmt.exact)          # ensure 0/1
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
    
    
    
    ## If HDMT fails, fall back to BH q-values (possibly screened)
    
    if (inherits(tmp_locfdr, "try-error")) {
      selected_values <- p_vec_all
      idx_detected    <- which(selected_values <FDR_level)  
      
      locfdr <- rep(NA_real_, length(select_otu))
    } else {
      p_vec_all[select_otu] <- tmp_locfdr
      idx_sub <- which(tmp_locfdr<= FDR_level)
      idx_detected <- select_otu[idx_sub]
    }
    
    
    globalp <- min(p_vec_all, na.rm = TRUE)   
    
    runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
    
    
    return(list(idx_detected=idx_detected,
                qval.med =p_vec_all,
                runtime_sec = runtime_sec,
                global_p = globalp,
                pval.alpha = p1,
                pval.beta =p2,
                beta_l = res1$beta_l,
                beta_r = res2$beta_r,
                taxa_detected = colnames(mediators)[idx_detected],
                p_matrix =p_matrix,
                filter_index = select_otu))
  }
  
}


####1.2 other methods
{
  # ── LDM ──────────────────────────────────────────────────────────────────────
  ldm_sim <- function(count_m, treat_cov, y,
                      fdr_nominal = 0.05, seed = 67817) {
    t0 <- proc.time()
    library(LDM)
    
    rand_id   <- paste0(sample(letters, 12), collapse = "")
    mat_name  <- paste0("RA_mat_", rand_id)
    meta_name <- paste0("meta_", rand_id)
    
    assign(mat_name,  as.matrix(count_m), envir = .GlobalEnv)
    assign(meta_name, data.frame(trt = treat_cov, Y = y), envir = .GlobalEnv)
    
    on.exit({
      if (exists(mat_name,  envir = .GlobalEnv)) rm(list = mat_name,  envir = .GlobalEnv)
      if (exists(meta_name, envir = .GlobalEnv)) rm(list = meta_name, envir = .GlobalEnv)
    }, add = TRUE)
    
    fmla <- parse(text = paste0(mat_name, " ~ trt + Y"))[[1]]
    
    res <- LDM::ldm(
      formula        = fmla,
      data           = get(meta_name, envir = .GlobalEnv),
      seed           = seed,
      fdr.nominal    = fdr_nominal,
      test.mediation = TRUE
    )
    

    P <- as.matrix(res$p.otu.omni)
    taxa_names <- colnames(P)
    if (is.null(taxa_names)) taxa_names <- paste0("taxon_", seq_len(ncol(P)))
    
    # Joint p-value via MultiMed::medTest.SBMH (same as paper)
    p_joint <- MultiMed::medTest.SBMH(
      P[1, ], P[2, ],
      MCP.type = "FDR", t1 = 0.05 / 2, t2 = 0.05 / 2
    )
    
    # Detected mediators from LDM built-in procedure
    det     <- res$med.detected.otu.omni
    det_sel <- rep(FALSE, ncol(P))
    if (is.logical(det))        det_sel <- det
    else if (is.numeric(det))   { v <- as.integer(det); det_sel[v[v >= 1 & v <= ncol(P)]] <- TRUE }
    else if (is.character(det)) det_sel <- taxa_names %in% det
    
    runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
    
    list(
      discoveries  = which(det_sel),
      otu_detected = taxa_names[det_sel],
      p_med        = p_joint,
      p_otu        = res$p.otu.omni,
      beta         = res$beta,
      global_p     = res$med.p.global.omni,
      runtime_sec  = runtime_sec
    )
  }
  
  # ── CMM ──────────────────────────────────────────────────────────────────────
  ccmm_sim <- function(count1, treat1, y1, sudo_count = 0.5,
                       method = c("boot", "normal")) {
    set.seed(123)
    method <- match.arg(method)
    t0 <- proc.time()[["elapsed"]]
    treat1_vec <- as.vector(treat1)
    
    # Helper: extract p-values from confidence intervals
    get_p_from_ci <- function(IDEs, p_matrix_cmm, method = "fdr", ci_level = 0.95) {
      ci_lower <- p_matrix_cmm[1, ]; ci_upper <- p_matrix_cmm[2, ]
      z_alpha2 <- qnorm(1 - (1 - ci_level) / 2)
      se <- (ci_upper - ci_lower) / (2 * z_alpha2)
      bad_se <- !is.finite(se) | se <= 0
      z_val <- IDEs / se; z_val[bad_se] <- NA_real_
      p_val <- 2 * pnorm(-abs(z_val))
      p_adj <- p.adjust(p_val, method = method)
      data.frame(IDE = IDEs, SE = se, z = z_val, p_value = p_val, p_adj = p_adj)
    }
    
    realdata_prop <- (count1 + sudo_count) / rowSums(count1 + sudo_count)
    M <- realdata_prop
    
    if (method == "boot") {
      res_ccmm <- ccmm::ccmm(y = as.numeric(y1), M = M, tr = treat1_vec, n.boot = 500)
      res_ccmm_p <- get_p_from_ci(res_ccmm$IDEs, res_ccmm$IDE.CIs)
      p_adj_cmm <- res_ccmm_p$p_adj
      idx_cmm <- which(p_adj_cmm < 0.05)
      global_p <- ifelse(res_ccmm$TIDE.CI[1] > 0 | res_ccmm$TIDE.CI[2] < 0, 1e-6, 1)
    } else {
      res_ccmm <- ccmm::ccmm(y1, realdata_prop, treat1_vec, method.est.cov = "normal")
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
  
  
  # ── CRAmed wrapper ──────────────────────────────────────────────────────────
  
  {
    {
      
      ####Redefine CRAMed
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
        
        # -----------------------------
        # helper: build p-value table (non-selection taxa set to 1)
        # -----------------------------
        build_p_table <- function(M_mat, selection_set, pvl_lst, final.p, method, mediation_set = integer(0)) {
          p <- ncol(M_mat)
          taxa <- colnames(M_mat)
          if (is.null(taxa)) taxa <- paste0("taxa_", seq_len(p))
          
          # default: not in selection_set => 1
          p_d0crt <- rep(1, p)
          p_model <- rep(1, p)
          q_d0crt <- rep(1, p)
          q_model <- rep(1, p)
          p_joint <- rep(1, p)
          
          selection_set <- sort(unique(as.integer(selection_set)))
          selection_set <- selection_set[selection_set >= 1 & selection_set <= p]
          
          if (length(selection_set) > 0) {
            # pvl_lst is indexed by taxa id (original indices)
            p_d0crt[selection_set] <- pvl_lst[selection_set]
            
            # final.p is ordered by selection_set (j=1..len(selection_set))
            p_model[selection_set] <- as.numeric(final.p)
            
            # adjusted (within selection_set only), same as original logic
            q_d0crt[selection_set] <- p.adjust(p_d0crt[selection_set], method = method)
            q_model[selection_set] <- p.adjust(p_model[selection_set], method = method)
            
            # joint p used for selection
            p_joint[selection_set] <- pmax(q_d0crt[selection_set], q_model[selection_set])
          }
          
          p_table <- data.frame(
            taxa = taxa,
            in_select = FALSE,
            is_mediator = FALSE,
            p_d0crt = p_d0crt,
            p_model = p_model,
            q_d0crt = q_d0crt,
            q_model = q_model,
            p_joint = p_joint,
            stringsAsFactors = FALSE
          )
          if (length(selection_set) > 0) p_table$in_select[selection_set] <- TRUE
          if (length(mediation_set) > 0) {
            mediation_set <- mediation_set[mediation_set >= 1 & mediation_set <= p]
            p_table$is_mediator[mediation_set] <- TRUE
          }
          p_table
        }
        
        # -----------------------------
        # Conditional generator functions (original)
        # -----------------------------
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
            phi <- zinb.fit$theta
            M_bar <- zinb.fit$fitted.values
            M_res <- zinb.fit$residuals
            x <- cbind(1, Exposure)
            lamda <- exp(x %*% matrix(alpha))
            
            samp_M <- matrix(NA, n, n.times)
            for (j in 1:n.times) {
              set.seed(j)
              for (i in 1:n) {
                samp_M[i, j] <- rnbinom(1, mu = lamda[i], size = phi)
              }
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
                if (Z1[i] == 1) samp_M[i, j] <- 0
                if (Z1[i] == 0) samp_M[i, j] <- rnbinom(1, mu = lamda[i], size = phi)
              }
            }
            
            res_M_samp <- lapply(seq_len(n.times), function(j) {
              if (sum(samp_M[, j]) == 0) return(data.frame(t(rep(0, n)))) # fallback
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
            for (i in 1:n) {
              samp_M[i, j] <- rnbinom(1, mu = lamda[i], size = phi)
            }
          }
          
          res_M_samp <- lapply(seq_len(n.times), function(j) {
            if (sum(samp_M[, j]) == 0) return(data.frame(t(rep(0, n))))
            Yd <- data.frame(Exposure = Exposure, Mediator = samp_M[, j])
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
              if (Z1[i] == 1) samp_M[i, j] <- 0
              if (Z1[i] == 0) samp_M[i, j] <- rpois(1, lambda = lamda[i])
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
        
        # -----------------------------
        # main
        # -----------------------------
        p <- ncol(M_mat)
        n <- nrow(M_mat)
        
        # prefilter / selection_set
        if (prefilter) {
          set.seed(123)
          cv_lasso <- cv.glmnet(cbind(M_mat, Exposure), Y,
                                alpha = 1, family = modely, dfmax = as.integer(p/2))
          lamb <- cv_lasso$lambda.min
          opt_model <- glmnet(cbind(M_mat, Exposure), Y,
                              alpha = 1, lambda = lamb, family = modely, dfmax = as.integer(p/2))
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
          selection_set <- seq_len(p)  # no prefilter => treat as all selected
        }
        
        
        if (length(selection_set) == 0) selection_set <- integer(0)
        
        # Step: compute d0CRT p-values for taxa in selection_set
        pvl_lst <- rep(1, p)  # default 1 for all taxa (including non-selection)
        nde.p <- 1            # default
        
        if (p == 1) {
      
          selection_set <- 1
          
          # Mediator conditional
          if (modelm == "ZINB") Cond_M <- Creat_condition_zinb(M_mat, indx = 1, Exposure, n.times)
          if (modelm == "ZIP")  Cond_M <- Creat_condition_zip(M_mat, indx = 1, Exposure, n.times)
          if (modelm == "NB")   Cond_M <- Creat_condition_nb(M_mat, indx = 1, Exposure, n.times)
          
          M_res_ob <- Cond_M$res_m
          data0 <- data.frame(Y = Y, Exposure = Exposure)
          
          if (modely == "binomial") {
            eps_res <- Y - 1 / (1 + exp(-predict(lm(Y ~ Exposure, family = modely, data = data0))))
          } else {
            eps_res <- Y - predict(lm(Y ~ Exposure, family = modely, data = data0))
          }
          
          imp_obe <- abs(mean(M_res_ob * eps_res)) / mean(M_res_ob^2)
          
          library(plyr)
          list.s <- unlist(lapply(Cond_M$res_m_samp, length))
          list.sel <- which(list.s == 1)
          if (length(list.sel) != 0) {
            M_res_sample <- as.matrix(do.call(rbind.fill, Cond_M$res_m_samp[-list.sel]))
          } else {
            M_res_sample <- as.matrix(do.call(rbind.fill, Cond_M$res_m_samp))
          }
          
          var_lst_sample <- apply(M_res_sample, 1, function(v) mean((unlist(v))^2))
          t_lst <- abs(M_res_sample %*% eps_res / n) / var_lst_sample
          pvl_lst[1] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
        } else {
          # nde p-value part (Exposure residual permutation)
          Cond_E <- list(res_e_samp = list())
          for (j in 1:n.times) {
            set.seed(j)
            indx.e <- sample(1:n)
            Cond_E$res_e_samp[[j]] <- as.data.frame(t(Exposure[indx.e]))
          }
          
          set.seed(123)
          cv_lasso_null <- cv.glmnet(cbind(M_mat), Y,
                                     alpha = 1, family = modely, dfmax = as.integer(p/2))
          lamb_null <- cv_lasso_null$lambda.min
          model_res_null <- glmnet(cbind(M_mat), Y,
                                   alpha = 1, lambda = lamb_null, family = modely, dfmax = as.integer(p/2))
          
          if (modely == "binomial") {
            eps_res <- Y - 1 / (1 + exp(-predict(model_res_null, cbind(M_mat))))
          } else {
            eps_res <- as.numeric(Y - predict(model_res_null, cbind(M_mat)))
          }
          
          imp_exp <- abs(mean(Exposure * eps_res)) / mean(Exposure^2)
          
          library(plyr)
          list.s <- unlist(lapply(Cond_E$res_e_samp, length))
          list.sel <- which(list.s == 1)
          if (length(list.sel) != 0) {
            E_res_sample <- as.matrix(do.call(rbind.fill, Cond_E$res_e_samp[-list.sel]))
          } else {
            E_res_sample <- as.matrix(do.call(rbind.fill, Cond_E$res_e_samp))
          }
          
          var_lst_sample <- apply(E_res_sample, 1, function(v) mean((unlist(v))^2))
          t_lst <- abs(E_res_sample %*% eps_res / n) / var_lst_sample
          nde.p <- mean(c(1, ifelse(t_lst >= imp_exp, 1, 0)))
          
          # taxa-level d0CRT p-values only for selection_set; others remain 1
          if (length(selection_set) > 0) {
            for (j in seq_along(selection_set)) {
              indx <- selection_set[j]
              
              if (modelm == "ZINB") Cond_M <- Creat_condition_zinb(M_mat, indx, Exposure, n.times)
              if (modelm == "ZIP")  Cond_M <- Creat_condition_zip(M_mat, indx, Exposure, n.times)
              if (modelm == "NB")   Cond_M <- Creat_condition_nb(M_mat, indx, Exposure, n.times)
              
              M_res_ob <- Cond_M$res_m
              
              set.seed(123)
              cv_lasso_null2 <- cv.glmnet(cbind(M_mat[, -indx, drop = FALSE], Exposure), Y,
                                          alpha = 1, family = modely, dfmax = as.integer(p/2))
              lamb_null2 <- cv_lasso_null2$lambda.min
              model_res_null2 <- glmnet(cbind(M_mat[, -indx, drop = FALSE], Exposure), Y,
                                        alpha = 1, lambda = lamb_null2, family = modely, dfmax = as.integer(p/2))
              
              if (modely == "binomial") {
                eps_res2 <- Y - 1 / (1 + exp(-predict(model_res_null2, cbind(M_mat[, -indx, drop = FALSE], Exposure))))
              } else {
                eps_res2 <- as.numeric(Y - predict(model_res_null2, cbind(M_mat[, -indx, drop = FALSE], Exposure)))
              }
              
              imp_obe <- abs(mean(M_res_ob * eps_res2)) / mean(M_res_ob^2)
              
              library(plyr)
              list.s <- unlist(lapply(Cond_M$res_m_samp, length))
              list.sel <- which(list.s == 1)
              if (length(list.sel) != 0) {
                M_res_sample <- as.matrix(do.call(rbind.fill, Cond_M$res_m_samp[-list.sel]))
              } else {
                M_res_sample <- as.matrix(do.call(rbind.fill, Cond_M$res_m_samp))
              }
              
              var_lst_sample <- apply(M_res_sample, 1, function(v) mean((unlist(v))^2))
              t_lst <- abs(M_res_sample %*% eps_res2 / n) / var_lst_sample
              pvl_lst[indx] <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))
            }
          }
        }
        
        # offset
        tij_mat <- log(rowSums(M_mat))
        
        # -----------------------------
        # model-specific (ZINB / ZIP / NB)
        # -----------------------------
        aic <- bic <- numeric(0)
        final.p <- numeric(0)
        niea <- niep <- nie <- numeric(0)
        alpha.p <- gamma.p <- numeric(0)
        residualsm <- list()
        alpha_mat <- gamma_mat <- NULL
        
        # in original: only for selection_set
        if (length(selection_set) > 0) {
          aic <- rep(NA_real_, length(selection_set))
          bic <- rep(NA_real_, length(selection_set))
          final.p <- rep(NA_real_, length(selection_set))
          niea <- rep(NA_real_, length(selection_set))
          niep <- rep(NA_real_, length(selection_set))
          nie  <- rep(NA_real_, length(selection_set))
          alpha.p <- rep(NA_real_, length(selection_set))
          gamma.p <- rep(NA_real_, length(selection_set))
          residualsm <- vector("list", length(selection_set))
          alpha_mat <- matrix(NA_real_, length(selection_set), 2)
          gamma_mat <- matrix(NA_real_, length(selection_set), 2)
        }
        
        mediation_set <- integer(0)
        nie_keep <- niea_keep <- niep_keep <- niepval_keep <- nieap_keep <- niepp_keep <- numeric(0)
        
        if (modelm == "ZINB") {
          if (length(selection_set) > 0) {
            for (j in seq_along(selection_set)) {
              Y_data <- data.frame(Exposure = Exposure, Mediator = as.vector(M_mat[, selection_set[j]]))
              
              if (sum(Y_data$Mediator == 0) == 0) {
                nb.fit <- glm.nb(Mediator ~ Exposure, data = Y_data)
                residualsm[[j]] <- nb.fit$residuals
                final.p[j] <- summary(nb.fit)$coefficients[2, 4]
                alpha_mat[j, ] <- summary(nb.fit)$coefficients[1:2, 1]
              } else {
                if (sum(is.infinite(tij_mat)) != 0) {
                  zinb.fit <- try(
                    zeroinfl(Mediator ~ Exposure | Exposure,
                             data = Y_data, dist = "negbin", link = "logit"),
                    silent = TRUE
                  )
                } else {
                  zinb.fit <- try(
                    zeroinfl(Mediator ~ Exposure | Exposure,
                             offset = tij_mat, data = Y_data, dist = "negbin", link = "logit"),
                    silent = TRUE
                  )
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
                    t(matrix(c(alpha_mat[j, 2], gamma_mat[j, 2]))) %*% ginv(cov.mat) %*%
                      matrix(c(alpha_mat[j, 2], gamma_mat[j, 2])),
                    silent = TRUE
                  )
                  final.p[j] <- if (inherits(wald.t, "try-error")) NA else (1 - pchisq(as.numeric(wald.t), df = 2))
                }
              }
              
              beta.val <- beta_fit[selection_set[j]]
              alpha.val <- alpha_mat[j, ]
              gamma.val <- gamma_mat[j, ]
              
              # NIE components (original formula)
              niea[j] <- mean(beta.val *
                                (1 / (1 + exp(as.matrix(data.frame(rep(1, nrow(M_mat)))) %*% matrix(gamma.val[c(1)])))) *
                                (exp(as.matrix(data.frame(1, rep(1, nrow(M_mat)))) %*% matrix(alpha.val[1:2]) + tij_mat) -
                                   exp(as.matrix(data.frame(rep(1, nrow(M_mat)))) %*% matrix(alpha.val[c(1)]) + tij_mat)))
              
              niep[j] <- mean(beta.val *
                                (exp(as.matrix(data.frame(1, rep(1, nrow(M_mat)))) %*% matrix(alpha.val[1:2]) + tij_mat)) *
                                ((1 / (1 + exp(as.matrix(data.frame(1, rep(1, nrow(M_mat)))) %*% matrix(gamma.val[c(1, 2)])))) -
                                   (1 / (1 + exp(as.matrix(data.frame(rep(1, nrow(M_mat)))) %*% matrix(gamma.val[c(1)]))))))
              
              nie[j] <- niea[j] + niep[j]
            }
          }
          
          # selection based on joint p within selection_set
          index.p <- if (length(selection_set) > 0) pvl_lst[selection_set] else numeric(0)
          
          bind.p <- rbind(p.adjust(index.p, method), p.adjust(final.p, method))
          joint.p <- if (length(selection_set) > 0) apply(bind.p, 2, max) else numeric(0)
          
          index.mi <- which(joint.p <= FDR)
          mediation_set <- if (length(index.mi) > 0) selection_set[index.mi] else integer(0)
          
          # keep only selected mediators for NIE outputs (original)
          nie_keep  <- if (length(index.mi) > 0) nie[index.mi] else numeric(0)
          niea_keep <- if (length(index.mi) > 0) niea[index.mi] else numeric(0)
          niep_keep <- if (length(index.mi) > 0) niep[index.mi] else numeric(0)
          niepval_keep <- if (length(index.mi) > 0) joint.p[index.mi] else numeric(0)
          
          # NIEA.pval / NIEP.pval (original)
          bind.p2 <- rbind(p.adjust(index.p, method), p.adjust(alpha.p, method))
          joint.p2 <- if (length(selection_set) > 0) apply(bind.p2, 2, max) else numeric(0)
          nieap_keep <- if (length(index.mi) > 0) joint.p2[index.mi] else numeric(0)
          
          bind.p3 <- rbind(p.adjust(index.p, method), p.adjust(gamma.p, method))
          joint.p3 <- if (length(selection_set) > 0) apply(bind.p3, 2, max) else numeric(0)
          niepp_keep <- if (length(index.mi) > 0) joint.p3[index.mi] else numeric(0)
          
          # Build p_table (non-selection taxa set to 1)
          p_table <- build_p_table(M_mat, selection_set, pvl_lst, final.p, method, mediation_set)
          
          # CI part (kept as original behavior: only returns CI for selected mediators; CI code太长且与p值无关，这里不改动结构)
          out <- list(
            Mediators = mediation_set,
            NDE = beta_2,
            NIE = nie_keep, NIEA = niea_keep, NIEP = niep_keep,
            NDE.pval = nde.p,
            NIE.pval = niepval_keep, NIEA.pval = nieap_keep, NIEP.pval = niepp_keep,
            AIC = aic, BIC = bic,
            residualsy = residualsy,
            residualsm = residualsm,
            p_table = p_table
          )
          
          if (CI) {
            
            warning(" no CI")
          }
          
          return(out)
        }
        
        if (modelm == "ZIP") {
          if (length(selection_set) > 0) {
            for (j in seq_along(selection_set)) {
              Y_data <- data.frame(Exposure = Exposure, Mediator = as.vector(M_mat[, selection_set[j]]))
              
              if (sum(Y_data$Mediator == 0) == 0) {
                nb.fit <- glm.nb(Mediator ~ Exposure, data = Y_data)
                residualsm[[j]] <- nb.fit$residuals
                final.p[j] <- summary(nb.fit)$coefficients[2, 4]
                alpha_mat[j, ] <- summary(nb.fit)$coefficients[1:2, 1]
              } else {
                zip.fit <- try(
                  zeroinfl(Mediator ~ Exposure | Exposure,
                           offset = tij_mat, data = Y_data,
                           dist = "poisson", link = "logit"),
                  silent = TRUE
                )
                
                if (inherits(zip.fit, "try-error")) {
                  final.p[j] <- NA
                  aic[j] <- bic[j] <- NA
                  residualsm[[j]] <- NA
                } else {
                  aic[j] <- AIC(zip.fit)
                  bic[j] <- BIC(zip.fit)
                  residualsm[[j]] <- zip.fit$residuals
                  
                  alpha_mat[j, ] <- summary(zip.fit)$coefficients$count[1:2, 1]
                  gamma_mat[j, ] <- summary(zip.fit)$coefficients$zero[1:2, 1]
                  alpha.p[j] <- summary(zip.fit)$coefficients$count[2, 4]
                  gamma.p[j] <- summary(zip.fit)$coefficients$zero[2, 4]
                  
                  cov.mat <- matrix(c(
                    vcov(zip.fit)[2, 2], vcov(zip.fit)[2, 4],
                    vcov(zip.fit)[4, 2], vcov(zip.fit)[4, 4]
                  ), 2, 2)
                  
                  wald.t <- try(
                    t(matrix(c(alpha_mat[j, 2], gamma_mat[j, 2]))) %*% ginv(cov.mat) %*%
                      matrix(c(alpha_mat[j, 2], gamma_mat[j, 2])),
                    silent = TRUE
                  )
                  final.p[j] <- if (inherits(wald.t, "try-error")) NA else (1 - pchisq(as.numeric(wald.t), df = 2))
                }
              }
              
              beta.val <- beta_fit[selection_set[j]]
              alpha.val <- alpha_mat[j, ]
              gamma.val <- gamma_mat[j, ]
              
              niea[j] <- mean(beta.val *
                                (1 / (1 + exp(as.matrix(data.frame(rep(1, nrow(M_mat)))) %*% matrix(gamma.val[c(1)])))) *
                                (exp(as.matrix(data.frame(1, rep(1, nrow(M_mat)))) %*% matrix(alpha.val[1:2]) + tij_mat) -
                                   exp(as.matrix(data.frame(rep(1, nrow(M_mat)))) %*% matrix(alpha.val[c(1)]) + tij_mat)))
              
              niep[j] <- mean(beta.val *
                                (exp(as.matrix(data.frame(1, rep(1, nrow(M_mat)))) %*% matrix(alpha.val[1:2]) + tij_mat)) *
                                ((1 / (1 + exp(as.matrix(data.frame(1, rep(1, nrow(M_mat)))) %*% matrix(gamma.val[c(1, 2)])))) -
                                   (1 / (1 + exp(as.matrix(data.frame(rep(1, nrow(M_mat)))) %*% matrix(gamma.val[c(1)]))))))
              
              nie[j] <- niea[j] + niep[j]
            }
          }
          
          index.p <- if (length(selection_set) > 0) pvl_lst[selection_set] else numeric(0)
          
          bind.p <- rbind(p.adjust(index.p, method), p.adjust(final.p, method))
          joint.p <- if (length(selection_set) > 0) apply(bind.p, 2, max) else numeric(0)
          
          index.mi <- which(joint.p <= FDR)
          mediation_set <- if (length(index.mi) > 0) selection_set[index.mi] else integer(0)
          
          nie_keep  <- if (length(index.mi) > 0) nie[index.mi] else numeric(0)
          niea_keep <- if (length(index.mi) > 0) niea[index.mi] else numeric(0)
          niep_keep <- if (length(index.mi) > 0) niep[index.mi] else numeric(0)
          niepval_keep <- if (length(index.mi) > 0) joint.p[index.mi] else numeric(0)
          
          bind.p2 <- rbind(p.adjust(index.p, method), p.adjust(alpha.p, method))
          joint.p2 <- if (length(selection_set) > 0) apply(bind.p2, 2, max) else numeric(0)
          nieap_keep <- if (length(index.mi) > 0) joint.p2[index.mi] else numeric(0)
          
          bind.p3 <- rbind(p.adjust(index.p, method), p.adjust(gamma.p, method))
          joint.p3 <- if (length(selection_set) > 0) apply(bind.p3, 2, max) else numeric(0)
          niepp_keep <- if (length(index.mi) > 0) joint.p3[index.mi] else numeric(0)
          
          p_table <- build_p_table(M_mat, selection_set, pvl_lst, final.p, method, mediation_set)
          
          out <- list(
            Mediators = mediation_set,
            NDE = beta_2,
            NIE = nie_keep, NIEA = niea_keep, NIEP = niep_keep,
            NDE.pval = nde.p,
            NIE.pval = niepval_keep, NIEA.pval = nieap_keep, NIEP.pval = niepp_keep,
            AIC = aic, BIC = bic,
            residualsy = residualsy,
            residualsm = residualsm,
            p_table = p_table
          )
          
         
          
          return(out)
        }
        
        if (modelm == "NB") {
          if (length(selection_set) > 0) {
            for (j in seq_along(selection_set)) {
              Y_data <- data.frame(Exposure = Exposure, Mediator = as.vector(M_mat[, selection_set[j]]))
              nb.fit <- glm.nb(Mediator ~ Exposure, data = Y_data)
              
              residualsm[[j]] <- nb.fit$residuals
              aic[j] <- AIC(nb.fit)
              bic[j] <- BIC(nb.fit)
              
              final.p[j] <- summary(nb.fit)$coefficients[2, 4]
              alpha_mat[j, ] <- summary(nb.fit)$coefficients[1:2, 1]
              alpha.p[j] <- summary(nb.fit)$coefficients[2, 4]
              
              beta.val <- beta_fit[selection_set[j]]
              alpha.val <- alpha_mat[j, ]
              nie[j] <- mean(beta.val * (exp(as.matrix(data.frame(1, rep(1, nrow(M_mat)))) %*% matrix(alpha.val[1:2]) + tij_mat) -
                                           exp(as.matrix(data.frame(rep(1, nrow(M_mat)))) %*% matrix(alpha.val[c(1)]) + tij_mat)))
            }
          }
          
          index.p <- if (length(selection_set) > 0) pvl_lst[selection_set] else numeric(0)
          
          bind.p <- rbind(p.adjust(index.p, method), p.adjust(final.p, method))
          joint.p <- if (length(selection_set) > 0) apply(bind.p, 2, max) else numeric(0)
          
          index.mi <- which(joint.p <= FDR)
          mediation_set <- if (length(index.mi) > 0) selection_set[index.mi] else integer(0)
          
          nie_keep <- if (length(index.mi) > 0) nie[index.mi] else numeric(0)
          niepval_keep <- if (length(index.mi) > 0) joint.p[index.mi] else numeric(0)
          
          p_table <- build_p_table(M_mat, selection_set, pvl_lst, final.p, method, mediation_set)
          
          out <- list(
            Mediators = mediation_set,
            NDE = beta_2,
            NIE = nie_keep,
            NDE.pval = nde.p,
            NIE.pval = niepval_keep,
            AIC = aic, BIC = bic,
            residualsy = residualsy,
            residualsm = residualsm,
            p_table = p_table
          )
          
          if (CI) {
            warning(" no CI")
          }
          
          return(out)
        }
        
        stop("modelm must be one of: 'ZINB', 'ZIP', 'NB'.")
      }
      
      CRAmed_sim <- function(count1, treat1, y1) {
        t0 <- proc.time() 
        
        results.infant <- CRAmed_fullp(
          M_mat   = count1,
          Y       = y1,
          Exposure= as.matrix(treat1),
          n.times = 100,
          n.perm  = 100,
          CI      = FALSE
        )
        
        ss <- results.infant$p_table
        idx1 <- results.infant$Mediators
        
        #####
        p_joint <- setNames(results.infant$p_table$p_joint, results.infant$p_table$taxa)
        
        t1 <- proc.time() - t0
        runtime <- unname(t1["elapsed"])   # 
        
        taxa <- results.infant$p_table$taxa
        
        return(list(
          discoveries    = idx1,
          p_med = p_joint,
          runtime_sec = runtime
        ))
      }
      
    }
  }
  
  
  # ── MarZIC ──────────────────────────────────────────────────────────────────
  Mar_sim <- function(count1, treat1, y1, alpha = 0.05,
                      adjust_method = "fdr", num_cores = 1) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(123)
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
                             add_lib_size = TRUE,   
                             pseudo = 0.5,          
                             alpha_out = 1) {     
    
    set.seed(42)
    t0 <- proc.time()
    stopifnot(length(treat1) == nrow(count1),
              length(y1)     == nrow(count1))
    
    ## ---- original names ----
    orig_names_full <- colnames(count1)
    if (is.null(orig_names_full)) orig_names_full <- paste0("V", seq_len(ncol(count1)))
    
    ## ---- numeric matrix ----
    X0 <- as.matrix(count1)
    storage.mode(X0) <- "numeric"
    
    ## ---- drop all-zero samples ----
    rs <- rowSums(X0)
    keep_samp <- rs > 0
    if (!all(keep_samp)) {
      X0     <- X0[keep_samp, , drop = FALSE]
      treat1 <- treat1[keep_samp]
      y1     <- y1[keep_samp]
      rs     <- rs[keep_samp]
    }
    
    ## ---- drop all-zero / zero-variance taxa ----
    keep_taxa <- (colSums(X0) > 0) & (apply(X0, 2, var) > 0)
    
    if (sum(keep_taxa) == 0) {
      runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
      p_full <- rep(1, length(orig_names_full)); names(p_full) <- orig_names_full
      return(list(discoveries = integer(0), p_med = p_full, runtime_sec = runtime_sec))
    }
    
    X_use <- X0[, keep_taxa, drop = FALSE]
    orig_names_use <- orig_names_full[keep_taxa]
    
    ## ---- CLR transform (microbiome as continuous mediators) ----
    X_use <- X_use + pseudo
    prop  <- sweep(X_use, 1, rowSums(X_use), "/")
    logp  <- log(prop)
    clrM  <- logp - rowMeans(logp)
    
    ## ---- safe names for modeling (kept taxa only) ----
    safe_names <- make.names(orig_names_use, unique = TRUE)
    colnames(clrM) <- safe_names
    mediator <- safe_names
    

    if (add_lib_size) {
      lib_size <- log(rs + 1)
    } else {
      lib_size <- NULL
    }
    
    ## ---- treatment: force numeric 0/1, then shift to 1/2 ----
    treat01 <- treat1
    if (is.factor(treat01)) treat01 <- as.character(treat01)
    treat01 <- as.numeric(treat01)
    
    if (!all(treat01 %in% c(0, 1))) stop("treat1 must be coded as 0/1.")
    if (length(unique(treat01)) < 2) stop("treat1 has <2 levels after filtering samples.")
    
    treat12 <- treat01 + 1  # 0/1 -> 1/2
    df <- data.frame(
      treatment = treat12,
      outcome   = as.numeric(y1),
      check.names = FALSE
    )
    
    if (add_lib_size) df$lib_size <- as.numeric(lib_size)
    
    df <- cbind(df, as.data.frame(clrM, check.names = FALSE))
    
    ## ---- mediation_data ----
    exper <- multimedia::mediation_data(
      df,
      outcomes   = "outcome",
      treatments = "treatment",
      mediator  = mediator,
      pretreatments = if (add_lib_size) "lib_size" else NULL
    )
    
    ## ---- choose estimators ----
    med_est <- multimedia::lm_model()
  
    out_est <- if (!add_lib_size) {
      multimedia::lm_model()
    } else {
      multimedia::glmnet_model(alpha = alpha_out)
    }
    
    # out_est <- multimedia::lm_model()
    ## ---- fit + estimate ----
    fit <- multimedia::multimedia(
      exper,
      outcome_estimator   = out_est,
      mediation_estimator = med_est
    ) |>
      multimedia::estimate(exper)
    
    ## ---- null contrast + FDR summary ----
    contrast <- multimedia::null_contrast(
      fit, exper,
      nullification = "T->M",
      f = multimedia::indirect_pathwise
    )
    
    
    # 
    # 
    sig_res <- multimedia::fdr_summary(
      contrast,
      effect  = "indirect_pathwise",
      q_value = q_value
    )
    
    
    kk <- sig_res$fdr_hat
    
  
    q_grid <- c(0.01, 0.02, 0.05, 0.10, 0.20)
    eps    <- 1e-4
    q_grid <- sort(unique(q_grid))
    stopifnot(all(q_grid > 0), all(q_grid < 1))
    

    fdr_list <- lapply(q_grid, function(q) {
      fdr <- multimedia::fdr_summary(
        contrast,
        effect  = "indirect_pathwise",
        q_value = q
      )
      fdr[fdr$source == "real", , drop = FALSE]   
    })
    

    discoveries_by_q <- setNames(vector("list", length(q_grid)), paste0("q=", q_grid))
    for (i in seq_along(q_grid)) {
      fdr_i <- fdr_list[[i]]
      sel_mediator <- fdr_i$mediator[fdr_i$keep]  
      idx_safe <- match(sel_mediator, safe_names)
      idx_safe <- idx_safe[!is.na(idx_safe)]
 
      idx_full <- which(keep_taxa)[idx_safe]
      discoveries_by_q[[i]] <- idx_full
    }

    q_disc_safe <- rep(1, length(safe_names))
    names(q_disc_safe) <- safe_names
    
    for (i in seq_along(q_grid)) {
      q <- q_grid[i]
      fdr_i <- fdr_list[[i]]
      sel_mediator <- fdr_i$mediator[fdr_i$keep]
      idx_safe <- match(sel_mediator, safe_names)
      idx_safe <- idx_safe[!is.na(idx_safe)]
  
      to_set <- idx_safe[q_disc_safe[idx_safe] == 1]
      if (length(to_set)) q_disc_safe[to_set] <- max(q - eps, 0)
    }
    
 
    q_disc_full <- rep(1, length(orig_names_full))
    names(q_disc_full) <- orig_names_full
    q_disc_full[keep_taxa] <- as.numeric(q_disc_safe[safe_names])
    q_disc_full[is.na(q_disc_full)] <- 1
    

    {
      key_count <- make.names(colnames(count1))
      key_fdr   <- make.names(fdr_i$mediator)
      
     
      fdr_map <- tapply(fdr_i$fdr_hat, key_fdr, function(x) min(x, na.rm = TRUE))
 
      fdr_hat_vec <- unname(fdr_map[key_count])
      names(fdr_hat_vec) <- colnames(count1)
 
      sum(!is.na(fdr_hat_vec))
      head(fdr_hat_vec)
    }

    discoveries_default <- which(q_disc_full <= q_value)
    
    runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])
    
    return(list(
      discoveries_default = discoveries_default,
      q_value_disc        = fdr_hat_vec,
      runtime_sec         = runtime_sec,
      discoveries_by_q    = discoveries_by_q
    ))
    
  }
  
}




###############################============2. run all models

{
  res1_prop <- CAMRA( outcome= y1,
                      treatment= treat1,
                      mediators = count1,  #[,idx_1]
                      screen = FALSE,
                      CClasso = FALSE,
                      seed = 123)
  
  


  {

    {
      idx_PRO <- res1_prop$idx_detected
      otu_PRO <- colnames(count1)[idx_PRO] #colnames(count1)[idx_1][idx_PRO]

      PRO_beta <- rbind(res1_prop$beta_l,res1_prop$beta_r)
      IDE_PRO <- res1_prop$beta_l*res1_prop$beta_r
      colnames(PRO_beta) <- colnames(count1)
      rownames(PRO_beta) <- c("treat",'outcome')
    }
    PRO_beta[,idx_PRO]
    
  }
  
  ########results for one side
  res_l <- recover_l_PALM(count1,
                          treat1)
  selected_values1 <- res_l$p

  res_r <-  recover_r(
    y= y1,
    treat_cov = treat1,
    count_matrix = count1)
  
  selected_values2 <- res_r$p

}

#################LDM
{

  
  
  {
  res.ldm.med <- ldm_sim(count1,treat1,y1,seed = 42)
  p_LDM  <- res.ldm.med$p_otu
  LDM_beta <- res.ldm.med$beta
  otu_LDM <- res.ldm.med$otu_detected
  idx_LDM <- which(colnames(count1) %in% otu_LDM)
  }
}




res_CRA <-CRAmed_sim(count1,treat1,y1)
idx_CRA <- res_CRA[[1]]
otu_CRA <-  colnames(count1)[idx_CRA] 


res_hima1 <- HIMA_micro_sim1(count1,treat1,y1)
otu_hima1 <- colnames(count1)[which(res_hima1[[2]] < 0.05)]


res_multimedia<- multimedia_sim(count1,treat1,y1)
# res_multimedia <- multimedia_sim(count1,treat1,y1)
otu_multimedia <- colnames(count1)[which(res_multimedia[[2]] < 0.05)]

############MarZIC
{

  
  res_mar  <- Mar_sim(count1,treat1,y1,num_cores = 5)
  otu_mar <- colnames(count1)[which(res_mar[[2]] < 0.05)]
}


##################global test
{
  library(matrixStats)
  permanovaFL_res <- permanovaFL_sim(count1, treat1, y1)
  MODIMA_res <- MODIMA_sim(count1, treat1, y1)
  Medtest_res <- Medtest_sim(count1, treat1, y1)
  cmm_res <- ccmm_sim(count1, treat1, y1,method="normal")
  
  cmm_res$global_p
  res.ldm.med$med.p.global.omni
  permanovaFL_res$global_p
  MODIMA_res$global_p
  Medtest_res$global_p
}




###################============3. plot for real data
base_fs <- 22
library(ggplot2)
theme_big <- theme_bw(base_size = base_fs) +
  theme(
    axis.title   = element_text(size = rel(1.05)),
    axis.text    = element_text(size = rel(0.95)),
    legend.title = element_text(size = rel(1.00)),
    legend.text  = element_text(size = rel(0.95)),
    strip.text   = element_text(size = rel(0.95)),
    plot.title   = element_text(size = rel(1.10))
  )
################========= 3.1 panel c
treat_lab <- c("CHN","USA")
otu_tab <- count1

######### treat -> OTU
{
  make_otu_label <- function(otu) {
    otu <- as.character(otu)
    
    has_bar <- grepl("\\|", otu)
    
    left  <- sub("\\|.*$", "", otu)
    right <- sub("^.*\\|", "", otu)

    out <- ifelse(has_bar & left != right,
                  paste0(left, " | ", right),
                  otu)
    
    out
  }
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(compositions)  # clr()

  
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
  

  stopifnot(is.factor(meta2$treat))
  lv <- levels(meta2$treat)   
  

  cols_fill  <- c("#E69F00", "#56B4E9")
  cols_color <- c("#B05D00", "#1F78B4")
  
  pal_fill  <- setNames(cols_fill[seq_along(lv)],  lv)
  pal_color <- setNames(cols_color[seq_along(lv)], lv)

  stopifnot(all(lv %in% names(pal_fill)))
  

  prev_df <- (otu_tab[, idx_PRO, drop = FALSE] > 0) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::pivot_longer(cols = -sample, names_to = "otu", values_to = "detected") %>%
    dplyr::left_join(meta2, by = "sample") %>%
    dplyr::group_by(treat, otu) %>%
    dplyr::summarise(
      prevalence = mean(detected, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(otu_label = make_otu_label(otu))
  library(forcats)
  prev_df2 <- prev_df %>%
    group_by(otu_label) %>%
    mutate(.ord = mean(prevalence, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      otu_label = fct_reorder(otu_label, .ord),
      panel = ""   
    )
  prev_df2<- prev_df2 %>%
    dplyr::mutate(
      otu_label = factor(
        otu_label
      )
    )

  p_prev <- ggplot(prev_df2, aes(x = treat, y = prevalence, fill = treat)) +
    geom_col(width = 0.75) +
    # facet_grid(panel ~ otu_label) +
    scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = pal_fill, guide = "none") +   
    labs(x = NULL, y = "Prevalence (Non-zero proportion)", fill = NULL) +
    theme_bw(base_size = 12) +
    theme_big +
    theme(
      strip.position   = "top",
      strip.placement  = "outside",
      strip.background = element_rect(fill = "grey95", color = "grey40"),
      
      legend.position  = "none",       
      
      axis.text.x      = element_blank(), 
      axis.ticks.x     = element_blank(),  
      axis.line.x      = element_blank()   
    )
  
  p_prev
  
  ggsave(
    "figure/treat_prev.png",
    plot   = p_prev,
    width  = 8,
    height = 7,
    dpi    = 300,
    device = ragg::agg_png
  )
  

  pseudo <- 0.5

  
  abund_mat <- log10(sweep(count1+ pseudo, 1, rowSums(count1+ pseudo), "/")) 
  otu_clr <- abund_mat
  ## 1) CLR(long)
  clr_long <- otu_clr[, idx_PRO, drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "otu",
      values_to = "clr_abund"
    )
  

  cnt_long <- count1[, idx_PRO, drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "otu",
      values_to = "otu_count"
    )
  
  ab_df <- clr_long %>%
    dplyr::left_join(cnt_long, by = c("sample", "otu")) %>%
    dplyr::filter(!is.na(.data$otu_count) & .data$otu_count > 0) %>%
    dplyr::left_join(meta2, by = "sample") %>%
    dplyr::mutate(otu_label = make_otu_label(.data$otu))
  

  

  ab_df <- ab_df %>%
    dplyr::mutate(
      otu_label = factor(
        otu_label,
        levels = c( "Acidaminococcus intestini")  
      )
    )
  
  
  p_abund2 <- ggplot(ab_df, aes(x = treat, y = clr_abund, fill = treat)) +
    geom_violin(trim = TRUE, alpha = 0.45, color = "grey25") +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.25, color = "grey25") +

    geom_jitter(
      width = 0.08,
      size  = 0.45,
      alpha = 0.25,
      color = "black"
    ) +
    
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size  = 2.6,
      fill  = "white",
      color = "black",
      stroke = 0.8
    ) +
    
 
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.45,
      color = "black",
      linewidth = 0.35
    ) +
    
    # facet_wrap(~ otu_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = pal_fill) +
    # scale_color_manual(values = pal_color) +
    labs(x = "Treat", y =  expression(log[10]("RA")), fill = "Treat") +
    theme_bw(base_size = 12) +
    theme_big +
    theme(
      # legend.position = c(0.95, 0.95),
      legend.position = "none",
      legend.background = element_rect(fill = "white", colour = "grey80"),
      axis.text.x  = element_blank(), 
      axis.ticks.x = element_blank(),  
      axis.title.x = element_blank(),   
      strip.text       = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
  
  
  p_abund2
  
  
  
  ggsave(
    "figure/treat_compare.png",
    plot   = p_abund2,
    width  = 8,
    height = 7,
    dpi    = 300,
    device = ragg::agg_png
  )
  
}

cut_median <- 25  
legend_name <- "BMI group"

########## OTU -> outcome
{
  {
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(compositions)   
    
    
    otu_tab <- count1[, idx_PRO, drop = FALSE]
    otu_sel <- colnames(otu_tab)
    

    if (!is.null(names(y1)) && all(rownames(otu_tab) %in% names(y1))) {
      y_use <- y1[rownames(otu_tab)]
    } else {
      stopifnot(length(y1) == nrow(otu_tab))
      y_use <- y1
    }
    y_use <- as.numeric(y_use)
    
   
    bw_group <- ifelse(y_use <= cut_median, paste0("≤ ", cut_median), paste0("> ", cut_median))
    
    meta_bw <- data.frame(
      sample   = rownames(otu_tab),
      y        = y_use,
      bw_group = factor(bw_group, levels = c(paste0("≤ ", cut_median), paste0("> ", cut_median))),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(bw_group))
    
   
    meta_bw <- meta_bw %>%
      mutate(
        bw_group_n = {
          tt <- table(bw_group)
          factor(
            bw_group,
            levels = levels(bw_group),
            labels = paste0(levels(bw_group), " (n=", as.integer(tt[levels(bw_group)]), ")")
          )
        }
      )
    

    otu_tab2 <- otu_tab[meta_bw$sample, otu_sel, drop = FALSE]

    lv_bw <- levels(meta_bw$bw_group_n)   
    
    cols_bw <- c( "#008837" ,"#7B3294")     
    pal_bw  <- setNames(cols_bw[seq_along(lv_bw)], lv_bw)
    
    prev_df <- (otu_tab2 > 0) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample") %>%
      tidyr::pivot_longer(cols = -sample, names_to = "otu", values_to = "detected") %>%
      dplyr::left_join(meta_bw, by = "sample") %>%
      dplyr::group_by(bw_group_n, otu) %>%
      dplyr::summarise(
        prevalence = mean(detected),  
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(otu_label = make_otu_label(otu))
    
    ord <- prev_df %>%
      dplyr::group_by(.data$otu_label) %>%
      dplyr::summarise(prev_mean = mean(.data$prevalence), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$prev_mean)) %>%
      dplyr::pull(.data$otu_label)
    
    prev_df <- prev_df %>%
      mutate(otu_label = factor(otu_label, levels = ord))
    
    
    library(forcats)
    
    prev_df_bw2 <- prev_df %>%
      group_by(otu_label) %>%
      mutate(.ord = mean(prevalence, na.rm = TRUE)) %>%   
      ungroup() %>%
      mutate(
        otu_label = fct_reorder(otu_label, .ord),
        panel = ""               
      )
 
    p_prev_bw <- ggplot(prev_df_bw2, aes(x = bw_group_n, y = prevalence, fill = bw_group_n)) +
      geom_col(width = 0.75) +
      # facet_grid(panel ~ otu_label) +                    
      scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0.05))) +
      scale_fill_manual(values = pal_bw, name = legend_name) +
      labs(x = NULL, y = "Prevalence (Non-zero proportion)") +
      theme_bw(base_size = 12) + theme_big +
      theme(
        strip.position   = "top",
        strip.placement  = "outside",
        strip.background = element_rect(fill = "grey95", color = "grey40"),
        
        legend.position  = "none",      
        axis.text.x      = element_blank(),  
        axis.ticks.x     = element_blank(),  
        axis.line.x      = element_blank()  
      )
    
    p_prev_bw
    
    
    print(p_prev_bw)
 
    ggsave(
      "figure/bw_median_prev.png",
      plot   = p_prev_bw,
      width  = 8,
      height = 7,
      dpi    = 300,
      device = ragg::agg_png
    )
    

    transform <- "clr"   # "log1p" 或 "clr"
    
 
    count_long <- otu_tab2 %>%
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      pivot_longer(-sample, names_to = "otu", values_to = "otu_count")
    
    if (transform == "log1p") {
      abund_mat <- log10((otu_tab2 + 0.5)/rowSums((otu_tab2 + 0.5)))
      
      value_nm  <- "log1p_abund"
      ylab      <-  expression(log[10]("RA"))
      
    } else if (transform == "clr") {
      pseudo <- 0.5
      
      count_all2 <- count1[meta_bw$sample, , drop = FALSE]
      abund_all  <- log10(sweep(count_all2 + pseudo, 1, rowSums(count_all2 + pseudo), "/"))
      
      abund_mat <- abund_all[, idx_PRO, drop = FALSE]
      value_nm  <- "clr_abund"
      ylab      <- expression(log[10]("RA"))
      
    } else {
      stop("")
    }
    
    ab_df <- abund_mat %>%
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      pivot_longer(-sample, names_to = "otu", values_to = value_nm) %>%
      left_join(meta_bw, by = "sample") %>%
      left_join(count_long, by = c("sample","otu")) %>%  
      filter(!is.na(otu_count) & otu_count > 0) %>%      
      mutate(otu_label = make_otu_label(otu))
    
    ab_df$bw_group_n <- factor(ab_df$bw_group_n, levels = lv_bw)
    
    p_abund_bw <- ggplot(ab_df, aes(x = bw_group_n, y = .data[[value_nm]], fill = bw_group_n)) +
      geom_violin(trim = TRUE, alpha = 0.45, color = "grey25") +
      geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.25, color = "grey25") +
      geom_jitter(width = 0.08, size = 0.45, alpha = 0.25, color = "black") +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6,
                   fill = "white", color = "black", stroke = 0.8) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                   color = "black", linewidth = 0.35) +
      # facet_wrap(~ otu_label, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = pal_bw) +
      labs(y = ylab, fill = "Group") +
      theme_bw(base_size = 12) +
      theme_big +
      theme(
        legend.position = "none",
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text       = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank()
        
      )

    {
      library(dplyr)
      library(rstatix)
      library(ggpubr)
      

      mid_y_df <- ab_df %>%
        dplyr::group_by(otu_label) %>%
        dplyr::summarise(
          y_mid = (min(.data[[value_nm]], na.rm = TRUE) + max(.data[[value_nm]], na.rm = TRUE)) / 2,
          .groups = "drop"
        )
      

      ng <- nlevels(factor(ab_df$bw_group_n))
      
      if (ng == 2) {
        p_df <- ab_df %>%
          dplyr::group_by(otu_label) %>%
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
     
      p_df2 <- p_df %>%
        left_join(mid_y_df, by = "otu_label") %>%
        mutate(
          y.position = y_mid,
          p.label = rstatix::p_format(p, digits = 2)  
        )
      
      }
    print(p_abund_bw)
    
    

    ggsave(
      "figure/median_outcome_compare.png",
      plot   = p_abund_bw,
      width  = 8,
      height = 7,
      dpi    = 300,
      device = ragg::agg_png
    )
  }
  
}



###############======================3.2 panel a
{
  plot_taxa_upset_two_color <- function(
    taxa_sets,
    highlight_set   = "CAMRA",
    other_color     = "gray",
    highlight_color = "gray",
    min_intersection = 1,
    max_intersections = 40,
    sort_sets = c("descending", "ascending", FALSE),
    sort_intersections = c("descending", "ascending", FALSE),
    sort_intersections_by = c("cardinality", "degree", "ratio"),
    title = NULL,
    quiet_warnings = TRUE
  ) {
   
    
    sort_sets <- match.arg(sort_sets)
    sort_intersections <- match.arg(sort_intersections)
    sort_intersections_by <- match.arg(sort_intersections_by)
    

    if (!is.list(taxa_sets) || length(taxa_sets) < 2) {
      stop("  ")
    }
    if (is.null(names(taxa_sets)) || any(names(taxa_sets) == "")) {
      names(taxa_sets) <- paste0("Set", seq_along(taxa_sets))
    }
    

    nm <- names(taxa_sets)
    if (anyDuplicated(nm)) {
      taxa_sets <- lapply(unique(nm), function(k) {
        unique(unlist(taxa_sets[nm == k], use.names = FALSE))
      })
      names(taxa_sets) <- unique(nm)
    }
    
    taxa_sets <- lapply(taxa_sets, function(x) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(x)]
      unique(x)
    })
    
    set_names <- names(taxa_sets)
    all_taxa <- sort(unique(unlist(taxa_sets, use.names = FALSE)))
    if (length(all_taxa) == 0) stop("   ")
    
    df <- data.frame(taxon = all_taxa, check.names = FALSE)
    for (s in set_names) df[[s]] <- df$taxon %in% taxa_sets[[s]]
    
  
    set_sizes_plot <- ComplexUpset::upset_set_size(
      mapping = ggplot2::aes(fill = I(other_color))
    ) +
      ggplot2::labs(
        y = "Number of selected mediators",
        x = NULL
      )
   
    queries <- list()
    if (highlight_set %in% set_names) {
      queries <- list(
        ComplexUpset::upset_query(
          set = highlight_set,
          fill = highlight_color,
          only_components = "overall_sizes"
        )
      )
    }
    
    make_plot <- function() {
      
      
      big_themes <- ComplexUpset::upset_modify_themes(list(
        default = ggplot2::theme(
          text        = ggplot2::element_text(size = 28),
          plot.title  = ggplot2::element_text(size = 28),
          axis.title  = ggplot2::element_text(size = 20),
          axis.text   = ggplot2::element_text(size = 18),
          strip.text  = ggplot2::element_text(size = 18),
          legend.title= ggplot2::element_text(size = 18),
          legend.text = ggplot2::element_text(size = 17)
        ),

        intersections_matrix = ggplot2::theme(
          # axis.text.x = ggplot2::element_text(size = 16),
          axis.text.y = ggplot2::element_text(size = 18),
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank()   
        ),
        
     
        overall_sizes = ggplot2::theme(
          axis.text = ggplot2::element_text(size = 18),
          axis.title= ggplot2::element_text(size = 20)
        ),
        

        `Intersection size` = ggplot2::theme(
          axis.text = ggplot2::element_text(size = 18),
          axis.title= ggplot2::element_text(size = 20)
        )
      ))
      
     
      base_ann <- list(
        "Intersection size" =
          ComplexUpset::intersection_size(bar_number_threshold = 1) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
          ggplot2::labs(x = NULL) +
          ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            text        = ggplot2::element_text(size = 24),
            axis.title  = ggplot2::element_text(size = 24),
            axis.text   = ggplot2::element_text(size = 18)
          )
      )
      
      p <- ComplexUpset::upset(
        df,
        intersect = set_names,
        min_size = min_intersection,
        n_intersections = max_intersections,
        sort_sets = sort_sets,
        sort_intersections = sort_intersections,
        sort_intersections_by = sort_intersections_by,
        wrap = TRUE,
        set_sizes = set_sizes_plot,
        queries = queries,
        base_annotations = base_ann,
        themes = big_themes,
        stripes = ComplexUpset::upset_stripes(
          geom   = ggplot2::geom_segment(size = 0),   
          colors = c("white", "white") )               
        
      )
      
      if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
      
      p
    }
    
    
    if (quiet_warnings) suppressWarnings(make_plot()) else make_plot()
  }
  
  
  
  taxa_sets <- list(
    CAMRA     = otu_PRO,
    "LDM-med"        = otu_LDM,
    microHIMA      = otu_hima1,
    CRAmed        = otu_CRA,
    multimedia = otu_multimedia,
    MarZIC     = otu_Mar
  )
  normalize_taxa_sets <- function(taxa_sets) {
    normalize_taxa <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      
      
      x <- gsub("^(k|p|c|o|f|g|s)__", "", x, ignore.case = TRUE)  
      x <- gsub("^(kingdom|phylum|class|order|family|genus|species)[:\\.]", "", x,
                ignore.case = TRUE)
  
      x <- gsub("[\\.|:|\\||_|;]+", " ", x)
      
  
      x <- gsub("\\s+", " ", x)

      x <- tolower(trimws(x))
      
   
      x[!duplicated(x)]
    }
    
    stopifnot(is.list(taxa_sets))
    out <- lapply(taxa_sets, normalize_taxa)
    
   
    nm <- names(out)
    if (!is.null(nm) && anyDuplicated(nm)) {
      out2 <- lapply(unique(nm), function(k) unique(unlist(out[nm == k], use.names = FALSE)))
      names(out2) <- unique(nm)
      out <- out2
    }
    out
  }
  
  taxa_sets
  
  taxa_sets_std <- normalize_taxa_sets(taxa_sets)
  
  p <- plot_taxa_upset_two_color(
    taxa_sets_std,
    highlight_set = "CAMRA",
    min_intersection = 1,  
  )
  print(p)

  library(ggplot2)
  ggsave(
    "figure/upset.png",
    plot   = p,
    width  = 15,
    height = 7,
    dpi    = 300,
    device = ragg::agg_png
  )
  
}

###############======================3.3 panel b

##############CAMRA vocano
{
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(scales)
  })
  
  if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
  
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
    
    g <- ggplot2::ggplot(df, ggplot2::aes(x = beta, y = neglog10p)) +
   
      ggplot2::geom_point(
        data = dplyr::filter(df, group == "Not Sig"),
        mapping = ggplot2::aes(x = beta, y = neglog10p, color = group),
        size = 1.6, alpha = 0.45, na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = dplyr::filter(df, group == "Sig"),
        mapping = ggplot2::aes(x = beta, y = neglog10p, color = group),
        size = 1.9, alpha = 0.90, na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c("Not Sig" = "grey65", "Sig" = "red"),
        breaks = c("Not Sig", "Sig"),
        name = NULL
      ) +
      ggplot2::labs(
        x = xlab,
        y = if (use_fdr) expression(-log[10]("FDR (BH)")) else expression(-log[10]("P-value"))
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        legend.position = c(0.95, 0.60),
        legend.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(linewidth = 0.9),
        axis.ticks = ggplot2::element_line(linewidth = 0.8)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      theme_big
   
    if (!is.null(idx_circle)) {
      idx_circle <- unique(as.integer(idx_circle))
      g <- g +
        ggplot2::geom_point(
          data = dplyr::filter(df, idx %in% idx_circle),
          mapping = ggplot2::aes(x = beta, y = neglog10p),
          shape = 21,
          stroke = circle_stroke,
          size = circle_size,
          fill = NA,
          color = circle_color,
          inherit.aes = FALSE,
          na.rm = TRUE
        )
    }
    
    print(g)
    
    ggplot2::ggsave(
      outfile,
      plot   = g,
      width  = 8,
      height = 7,
      dpi    = 300,
      device = ragg::agg_png
    )
    
    invisible(g)
  }
  # ---- treat side ----
  plot_volcano_idxPRO_red(
    p_kk   = selected_values1,
    beta_kk= PRO_beta[1,],
    idx_PRO= idx_PRO,
    xlab   = "Effect estimate",
    outfile= "figure/treat_vocano.png"
  )
  
  # ---- outcome side ----
  plot_volcano_idxPRO_red(
    p_kk   = selected_values2,
    beta_kk= PRO_beta[2,],
    idx_PRO= idx_PRO,
    xlab   = "Effect estimate ",
    outfile= "figure/outcome_vocano.png"
  )
  
}

######################LDM vocano
{
  
  p_kk <- p_LDM[1,]
  beta_kk <- LDM_beta[1,]
  
  plot_volcano_idxPRO_red(
    p_kk   = p_LDM[1,],
    beta_kk= LDM_beta[1,],
    idx_PRO= idx_LDM,
    # idx_circle = idx_PRO,
    xlab   = "Effect estimate",
    outfile= "figure/treat_vocano_LDM.png"
  )
  
  plot_volcano_idxPRO_red(
    p_kk   = p_LDM[2,],
    beta_kk= LDM_beta[2,],
    idx_PRO= idx_LDM,
    # idx_circle = idx_PRO,
    xlab   = "Effect estimate",
    outfile= "figure/outcome_vocano_LDM.png"
  )
  
}



############################============= Table S4

{

  q_cut <- 0.05
  
  taxa <- colnames(count1)
  res_hima1[[2]][res_hima1[[2]] < 0.05] <- 0.05
  q_tab <- data.frame(
    Species    = taxa,
    CAMRA      = as.numeric(res1_prop$qval.med),
    CRAmed     = as.numeric(res_CRA[[2]]),
    `LDM-med`  = as.numeric(res.ldm.med$p_med),
    MarZIC     = as.numeric(res_mar$p_med),
    microHIMA  = as.numeric(res_hima1[[2]]),
    multimedia = as.numeric(res_multimedia[[2]]),
    check.names = FALSE
  )
  q_tab$Species <- as.character(q_tab$Species)
  
  key <- tolower(trimws(q_tab$Species))
  q_tab <- q_tab[order(key), , drop = FALSE]

  stopifnot(nrow(q_tab) == length(taxa))
  
  keep <- apply(q_tab[, -1, drop = FALSE], 1, function(x) any(!is.na(x) & x <= q_cut))
  q_tab_keep <- q_tab[keep, , drop = FALSE]
  
  q_min <- apply(q_tab_keep[, 2, drop = FALSE], 1, function(x) min(x, na.rm = TRUE))
  q_tab_keep <- q_tab_keep[order(q_min), , drop = FALSE]
  
  library(knitr)
  library(kableExtra)
  
  tex <- knitr::kable(
    q_tab_keep,
    format    = "latex",
    booktabs  = TRUE,
    row.names = FALSE,                
    digits    = 4,
    align     = c("l", rep("r", ncol(q_tab_keep) - 1)),
    escape    = FALSE,
    col.names = c(
      "Species",
      "\\textbf{CAMRA}",
      "\\textbf{CRAmed}",
      "\\textbf{LDM-med}",
      "\\textbf{MarZIC}",
      "\\textbf{microHIMA}",
      "\\textbf{multimedia}"
    ),
    caption = paste0(
      "\\textbf{Table S6:} Taxon-level mediation $q$-values from different methods in the real-data analysis. ",
      "The table lists taxa that were selected by at least one method. ",
      "We report $q$-values for taxa that were declared significant by at least one method at $q \\le ", q_cut, "$."
    ),
    label = "tab:S6_qvalues"
  ) |>
    kableExtra::kable_styling(position = "center", latex_options = "hold_position")
  
  cat(tex)
}
