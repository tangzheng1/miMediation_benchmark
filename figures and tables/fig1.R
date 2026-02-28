rm(list = ls())
library(data.table)
library(vegan)
library(energy) 
library(ggplot2)
library(dplyr)

# --- read files: downloaded from https://zenodo.org/records/14280080 ---

X0 <- fread("GALAXY_mOTUs_v25.tsv")   # taxa table
O0 <- fread("GALAXY_load.tsv")        # load table
sid_X <- X0[[1]]
sid_O <- O0[[1]]
common <- intersect(sid_X, sid_O)

X <- as.matrix(X0[sid_X %in% common, -1, with=FALSE])
O <- O0[sid_O %in% common, ][match(common, sid_O), ][[2]]  # load column

rs <- rowSums(X)
X_ra <- X / rs

eps <- 1e-8
logO <- log(O)

set.seed(1)
A_real <- O * X_ra
O_perm <- sample(O)                 # breaks dependence
logO_perm <- log(O_perm)
A_perm  <- O_perm * X_ra

########### PANEL A
set.seed(1)
top_taxa <- 50     # use top taxa by mean RA
n_pairs  <- 1000    # number of random pairs 
prev_min <- 0.10   # keep taxa present (>0) in >=10% samples
q_cut    <- 0.10   # BH threshold for “significant pairs”
eps_local <- eps  

# ---- select taxa ----
taxa_mean <- colMeans(X_ra, na.rm = TRUE)
taxa_prev <- colMeans(X_ra > 0, na.rm = TRUE)

idx <- order(taxa_mean, decreasing = TRUE)
idx <- idx[taxa_prev[idx] >= prev_min]
idx <- idx[1:min(top_taxa, length(idx))]

if (length(idx) < 10) stop("Too few taxa after filtering; lower prev_min or increase top_taxa.")

# ---- sample random pairs among selected taxa ----
pairs <- replicate(n_pairs, sample(idx, 2), simplify = FALSE)

# ---- function: slope + p-value for one pair ----
fit_pair <- function(jk, logO_vec) {
  j <- jk[1]; k <- jk[2]
  r <- log((X_ra[, j] + eps_local) / (X_ra[, k] + eps_local))
  fit <- lm(r ~ logO_vec)
  s <- summary(fit)$coefficients
  c(slope = s[2,1], p = s[2,4])
}

# ---- compute for real and permuted baseline ----
res_real <- t(vapply(pairs, fit_pair, numeric(2), logO_vec = logO))
res_perm <- t(vapply(pairs, fit_pair, numeric(2), logO_vec = logO_perm))

# ---- BH adjustment & summary stats ----
q_real <- p.adjust(res_real[, "p"], method = "BH")
q_perm <- p.adjust(res_perm[, "p"], method = "BH")

summ <- data.frame(
  setting = c("Real", "Permuted"),
  sd_slope = c(sd(res_real[, "slope"]), sd(res_perm[, "slope"])),
  iqr_abs  = c(IQR(abs(res_real[, "slope"])), IQR(abs(res_perm[, "slope"]))),
  n_sig    = c(sum(q_real < q_cut), sum(q_perm < q_cut))
)

# ---- plotting data ----
dfB <- rbind(
  data.frame(slope = res_real[, "slope"], setting = "Real"),
  data.frame(slope = res_perm[, "slope"], setting = "Permuted")
)

# ---- slope density plot ----
pal <- c("Real" = "#D55E00", "Permuted" = "#0072B2")

pA <- ggplot(dfB, aes(x = slope, fill = setting, color = setting)) +
  geom_density(alpha = 0.35, linewidth = 0.9, adjust = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  theme_classic(base_size = 12) +
  labs(
    #title = "(B) Log-ratio slopes vs log(load)",
    x = "Estimated load effect on log-ratio",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme(
    legend.position = c(0.80, 0.80),      # inside panel (x,y in [0,1])
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 11)
  )

pA <- pA +
  scale_fill_manual(
    values = pal,
    breaks = c("Real", "Permuted"),
    labels = c("Real load", "Permuted load")
  ) +
  scale_color_manual(
    values = pal,
    breaks = c("Real", "Permuted"),
    labels = c("Real load", "Permuted load")
  )

pA



########### PANEL B

# --- Global AA summary (mean log-AA across taxa) ---
m_real <- rowMeans(log(A_real + eps))
m_ind  <- rowMeans(log(A_perm + eps))

dfC <- bind_rows(
  data.frame(setting = "Real AA",            logload = logO,      m = m_real),
  data.frame(setting = "RA×permuted load",   logload = logO_perm, m = m_ind)
)

# Put panels in narrative order: Real first, baseline second
dfC$setting <- factor(dfC$setting, levels = c("Real AA", "RA×permuted load"))

# Common axes across panels 
xlim <- quantile(dfC$logload, c(0.01, 0.99), na.rm = TRUE)
ylim <- quantile(dfC$m,       c(0.01, 0.99), na.rm = TRUE)

# Correlation labels with 3 decimals
annoC <- dfC %>%
  group_by(setting) %>%
  summarise(r = cor(m, logload), .groups = "drop") %>%
  mutate(
    label = paste0("cor = ", sprintf("%.3f", r)),
    # place Real at top-right; baseline slightly left to reduce overlap
    x = ifelse(setting == "Real AA", xlim[2], xlim[1] + 0.65*(xlim[2]-xlim[1])),
    y = ylim[2]
  )

binsC <- 35

pB <- ggplot(dfC, aes(logload, m)) +
  geom_hex(bins = binsC) +
  scale_fill_viridis_c(option = "D", end = 0.95, name = "count") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "grey15") +
  facet_wrap(~setting, nrow = 1) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  geom_label(
    data = annoC,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.05,
    size = 3.2, label.size = 0.2
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "log(load)",
    y = "Mean log(AA) across taxa"
  )

pB







