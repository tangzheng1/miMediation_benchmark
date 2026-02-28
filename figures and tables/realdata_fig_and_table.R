suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(compositions)
  library(forcats)
  library(rstatix)
  library(ggpubr)
  library(scales)
  library(ComplexUpset)
  library(knitr)
  library(kableExtra)
})

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

# -- Output directory --
if (!dir.exists("figure")) dir.create("figure", recursive = TRUE)

# -- Global theme --
BASE_FONT_SIZE <- 22

theme_big <- theme_bw(base_size = BASE_FONT_SIZE) +
  theme(
    axis.title   = element_text(size = rel(1.05)),
    axis.text    = element_text(size = rel(0.95)),
    legend.title = element_text(size = rel(1.00)),
    legend.text  = element_text(size = rel(0.95)),
    strip.text   = element_text(size = rel(0.95)),
    plot.title   = element_text(size = rel(1.10))
  )

# Format OTU label: show "Genus | Species" when the two differ, else as-is.
make_otu_label <- function(otu) {
  otu   <- as.character(otu)
  left  <- sub("\\|.*$", "", otu)
  right <- sub("^.*\\|", "", otu)
  has_bar <- grepl("\\|", otu)
  ifelse(has_bar & left != right, paste0(left, " | ", right), otu)
}

# Build a metadata data.frame that maps sample names to a treatment factor.

build_treat_meta <- function(otu_tab, treat, labels = c("CHN", "USA")) {
  if (!is.null(names(treat)) && all(rownames(otu_tab) %in% names(treat))) {
    treat_use <- treat[rownames(otu_tab)]
  } else {
    stopifnot(length(treat) == nrow(otu_tab))
    treat_use <- treat
  }
  meta <- data.frame(
    sample = rownames(otu_tab),
    treat  = treat_use,
    stringsAsFactors = FALSE
  )
  if (is.numeric(meta$treat) || all(meta$treat %in% c(0, 1))) {
    meta$treat <- factor(meta$treat, levels = c(0, 1), labels = labels)
  } else {
    meta$treat <- factor(meta$treat)
  }
  meta
}

# Compute log10 relative abundance with pseudo-count.

compute_log10_ra <- function(count_mat, pseudo = 0.5) {
  log10(sweep(count_mat + pseudo, 1, rowSums(count_mat + pseudo), "/"))
}

#' Normalise taxa names for cross-method comparison (UpSet plot).
normalize_taxa <- function(x) {
  x <- as.character(trimws(x))
  x <- gsub("^(k|p|c|o|f|g|s)__", "", x, ignore.case = TRUE)
  x <- gsub("^(kingdom|phylum|class|order|family|genus|species)[:\\.]", "",
            x, ignore.case = TRUE)
  x <- gsub("[\\.|:|\\||_|;]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- tolower(trimws(x))
  unique(x)
}

# Apply normalize_taxa() to every element of a named list of taxon vectors.
normalize_taxa_sets <- function(taxa_sets) {
  stopifnot(is.list(taxa_sets))
  out <- lapply(taxa_sets, normalize_taxa)
  nm  <- names(out)
  if (!is.null(nm) && anyDuplicated(nm)) {
    out <- lapply(unique(nm), function(k)
      unique(unlist(out[nm == k], use.names = FALSE)))
    names(out) <- unique(nm)
  }
  out
}

# Colour palettes used throughout.
pal_treat <- function(levels) {
  cols_fill  <- c("#E69F00", "#56B4E9")
  cols_color <- c("#B05D00", "#1F78B4")
  list(
    fill  = setNames(cols_fill[seq_along(levels)],  levels),
    color = setNames(cols_color[seq_along(levels)], levels)
  )
}

pal_bw_group <- function(levels) {
  cols <- c("#008837", "#7B3294")
  setNames(cols[seq_along(levels)], levels)
}

# =============================================================================
# fig5c_prevalence — Prevalence bar plots (Treatment & Outcome sides)
# =============================================================================
# Requires: count1, treat1, y1, idx_PRO, otu_PRO
# Outputs:  figure/treat_prev.png, figure/bw_median_prev.png

# -- Configuration --
TREAT_LABELS <- c("CHN", "USA")
OUTCOME_CUT  <- 25
LEGEND_NAME  <- "BMI group"

# =========================================================================== #
#  1. Treatment → OTU prevalence                                              #
# =========================================================================== #
compute_prevalence <- function(otu_mat, group_df, group_col) {
  (otu_mat > 0) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "otu", values_to = "detected") %>%
    left_join(group_df, by = "sample") %>%
    group_by(!!sym(group_col), otu) %>%
    summarise(prevalence = mean(detected, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    mutate(otu_label = make_otu_label(otu))
}

plot_prevalence_bar <- function(prev_df, x_var, fill_var, pal, y_max = 0.25) {
  prev_df <- prev_df %>%
    group_by(otu_label) %>%
    mutate(.ord = mean(prevalence, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(otu_label = fct_reorder(otu_label, .ord))

  ggplot(prev_df, aes(x = .data[[x_var]], y = prevalence, fill = .data[[fill_var]])) +
    geom_col(width = 0.75) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "Prevalence (Non-zero proportion)") +
    theme_bw(base_size = 12) + theme_big +
    theme(
      legend.position = "none",
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank()
    )
}

# -- Treatment side --
meta_treat <- build_treat_meta(count1, treat1, labels = TREAT_LABELS)
pal        <- pal_treat(levels(meta_treat$treat))

prev_treat <- compute_prevalence(
  count1[, idx_PRO, drop = FALSE], meta_treat, "treat"
)

p_prev_treat <- plot_prevalence_bar(prev_treat, "treat", "treat", pal$fill)
ggsave("figure/treat_prev.png", p_prev_treat,
       width = 8, height = 7, dpi = 300, device = ragg::agg_png)

# -- Outcome (BMI) side --
y_use <- if (!is.null(names(y1)) && all(rownames(count1) %in% names(y1))) {
  y1[rownames(count1)]
} else {
  y1
}
y_use <- as.numeric(y_use)

meta_bw <- data.frame(
  sample   = rownames(count1),
  y        = y_use,
  bw_group = factor(
    ifelse(y_use <= OUTCOME_CUT,
           paste0("\u2264 ", OUTCOME_CUT),
           paste0("> ", OUTCOME_CUT)),
    levels = c(paste0("\u2264 ", OUTCOME_CUT), paste0("> ", OUTCOME_CUT))
  ),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(bw_group))

# Append sample sizes to labels
meta_bw <- meta_bw %>%
  mutate(bw_group_n = {
    tt <- table(bw_group)
    factor(bw_group,
           levels = levels(bw_group),
           labels = paste0(levels(bw_group), " (n=", as.integer(tt[levels(bw_group)]), ")"))
  })

pal_bw <- pal_bw_group(levels(meta_bw$bw_group_n))

prev_bw <- compute_prevalence(
  count1[meta_bw$sample, idx_PRO, drop = FALSE], meta_bw, "bw_group_n"
)

p_prev_bw <- plot_prevalence_bar(prev_bw, "bw_group_n", "bw_group_n", pal_bw)
ggsave("figure/bw_median_prev.png", p_prev_bw,
       width = 8, height = 7, dpi = 300, device = ragg::agg_png)

message("Done: prevalence plots saved.")
                 
# =============================================================================
# fig5c_abundance — Abundance violin + box plots (Treatment & Outcome)
# =============================================================================
# Requires: count1, treat1, y1, idx_PRO, meta_treat, meta_bw 
# Outputs:  figure/treat_compare.png, figure/median_outcome_compare.png

PSEUDO <- 0.5

# --------------------------------------------------------------------------- #
#  Shared plotting function                                                    #
# --------------------------------------------------------------------------- #
plot_abundance <- function(ab_df, x_var, fill_var, value_col, pal,
                           ylab = expression(log[10]("RA"))) {
  ggplot(ab_df, aes(x = .data[[x_var]], y = .data[[value_col]],
                     fill = .data[[fill_var]])) +
    geom_violin(trim = TRUE, alpha = 0.45, color = "grey25") +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.25, color = "grey25") +
    geom_jitter(width = 0.08, size = 0.45, alpha = 0.25, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6,
                 fill = "white", color = "black", stroke = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                 color = "black", linewidth = 0.35) +
    scale_fill_manual(values = pal) +
    labs(y = ylab) +
    theme_bw(base_size = 12) + theme_big +
    theme(
      legend.position  = "none",
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      strip.text       = element_blank(),
      strip.background = element_blank()
    )
}

# --------------------------------------------------------------------------- #
#  Build long-form abundance data for selected OTUs                            #
# --------------------------------------------------------------------------- #
build_abundance_long <- function(count_mat, idx_sel, meta_df,
                                 group_col, pseudo = 0.5) {
  log_ra <- compute_log10_ra(count_mat, pseudo = pseudo)
  abund_sel <- log_ra[, idx_sel, drop = FALSE]

  cnt_long <- count_mat[, idx_sel, drop = FALSE] %>%
    as.data.frame() %>% rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "otu", values_to = "otu_count")

  abund_sel %>%
    as.data.frame() %>% rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "otu", values_to = "clr_abund") %>%
    left_join(cnt_long, by = c("sample", "otu")) %>%
    filter(!is.na(otu_count) & otu_count > 0) %>%
    left_join(meta_df, by = "sample") %>%
    mutate(otu_label = make_otu_label(otu))
}

# =========================================================================== #
#  1. Treatment side                                                          #
# =========================================================================== #
ab_treat <- build_abundance_long(count1, idx_PRO, meta_treat, "treat",
                                  pseudo = PSEUDO)

pal_t <- pal_treat(levels(meta_treat$treat))

p_abund_treat <- plot_abundance(ab_treat, "treat", "treat", "clr_abund", pal_t$fill)
ggsave("figure/treat_compare.png", p_abund_treat,
       width = 8, height = 7, dpi = 300, device = ragg::agg_png)

# =========================================================================== #
#  2. Outcome (BMI) side                                                      #
# =========================================================================== #
ab_bw <- build_abundance_long(
  count1[meta_bw$sample, , drop = FALSE],
  idx_PRO, meta_bw, "bw_group_n", pseudo = PSEUDO
)
ab_bw$bw_group_n <- factor(ab_bw$bw_group_n, levels = levels(meta_bw$bw_group_n))

p_abund_bw <- plot_abundance(ab_bw, "bw_group_n", "bw_group_n", "clr_abund",
                              pal_bw)

# Optional: Wilcoxon test annotation
if (nlevels(factor(ab_bw$bw_group_n)) == 2) {
  p_df <- ab_bw %>%
    group_by(otu_label) %>%
    wilcox_test(clr_abund ~ bw_group_n) %>%
    ungroup()
  message("Wilcoxon p-values:\n", paste(capture.output(print(p_df)), collapse = "\n"))
}

ggsave("figure/median_outcome_compare.png", p_abund_bw,
       width = 8, height = 7, dpi = 300, device = ragg::agg_png)

message("Done: abundance plots saved.")

# =============================================================================
# fig5a_upset — UpSet plot comparing selected taxa across methods
# =============================================================================
# Requires: otu_PRO, otu_LDM, otu_hima1, otu_CRA, otu_multimedia, otu_Mar
# Outputs:  figure/upset.png

# --------------------------------------------------------------------------- #
#  UpSet plotting function                                                     #
# --------------------------------------------------------------------------- #
plot_taxa_upset <- function(taxa_sets,
                            highlight_set    = "CAMRA",
                            bar_color        = "gray",
                            highlight_color  = "gray",
                            min_intersection = 1,
                            max_intersections = 40,
                            sort_sets          = "descending",
                            sort_intersections  = "descending",
                            sort_intersections_by = "cardinality",
                            title = NULL) {

  stopifnot(is.list(taxa_sets), length(taxa_sets) >= 2)
  if (is.null(names(taxa_sets)) || any(names(taxa_sets) == ""))
    names(taxa_sets) <- paste0("Set", seq_along(taxa_sets))

  # De-duplicate set names
  nm <- names(taxa_sets)
  if (anyDuplicated(nm)) {
    taxa_sets <- lapply(unique(nm), function(k)
      unique(unlist(taxa_sets[nm == k], use.names = FALSE)))
    names(taxa_sets) <- unique(nm)
  }
  taxa_sets <- lapply(taxa_sets, function(x) unique(x[!is.na(x) & nzchar(x)]))

  set_names <- names(taxa_sets)
  all_taxa  <- sort(unique(unlist(taxa_sets, use.names = FALSE)))
  stopifnot(length(all_taxa) > 0)

  # Binary membership matrix
  df <- data.frame(taxon = all_taxa, check.names = FALSE)
  for (s in set_names) df[[s]] <- df$taxon %in% taxa_sets[[s]]

  # Set-size sidebar
  set_sizes_plot <- ComplexUpset::upset_set_size(
    mapping = aes(fill = I(bar_color))
  ) + labs(y = "Number of selected mediators", x = NULL)

  # Highlight query
  queries <- if (highlight_set %in% set_names) {
    list(ComplexUpset::upset_query(
      set = highlight_set, fill = highlight_color,
      only_components = "overall_sizes"
    ))
  } else list()

  # Theme overrides
  big_themes <- ComplexUpset::upset_modify_themes(list(
    default = theme(
      text         = element_text(size = 28),
      plot.title   = element_text(size = 28),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 18),
      strip.text   = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 17)
    ),
    intersections_matrix = theme(
      axis.text.y  = element_text(size = 18),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ),
    overall_sizes = theme(
      axis.text  = element_text(size = 18),
      axis.title = element_text(size = 20)
    ),
    `Intersection size` = theme(
      axis.text  = element_text(size = 18),
      axis.title = element_text(size = 20)
    )
  ))

  base_ann <- list(
    "Intersection size" =
      ComplexUpset::intersection_size(bar_number_threshold = 1) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL) +
      theme(axis.title.x = element_blank(),
            text       = element_text(size = 24),
            axis.title = element_text(size = 24),
            axis.text  = element_text(size = 18))
  )

  p <- suppressWarnings(ComplexUpset::upset(
    df, intersect = set_names,
    min_size = min_intersection, n_intersections = max_intersections,
    sort_sets = sort_sets, sort_intersections = sort_intersections,
    sort_intersections_by = sort_intersections_by,
    wrap = TRUE, set_sizes = set_sizes_plot,
    queries = queries, base_annotations = base_ann, themes = big_themes,
    stripes = ComplexUpset::upset_stripes(
      geom = geom_segment(size = 0), colors = c("white", "white")
    )
  ))

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

# =========================================================================== #
#  Build input and plot                                                       #
# =========================================================================== #
taxa_sets <- list(
  CAMRA      = otu_PRO,
  `LDM-med`  = otu_LDM,
  microHIMA  = otu_hima1,
  CRAmed     = otu_CRA,
  multimedia = otu_multimedia,
  MarZIC     = otu_Mar
)

taxa_sets_std <- normalize_taxa_sets(taxa_sets)

p_upset <- plot_taxa_upset(taxa_sets_std, highlight_set = "CAMRA",
                           min_intersection = 1)
print(p_upset)

ggsave("figure/upset.png", p_upset,
       width = 15, height = 7, dpi = 300, device = ragg::agg_png)

message("Done: UpSet plot saved.")

# =============================================================================
# fig5b_volcano — Volcano plots for CAMRA and LDM-med
# =============================================================================
# Requires: selected_values1/2, PRO_beta, idx_PRO,
#           p_LDM, LDM_beta, idx_LDM
# Outputs:  figure/treat_volcano.png,   figure/outcome_volcano.png
#           figure/treat_volcano_LDM.png, figure/outcome_volcano_LDM.png

# -- Configuration --
ALPHA   <- 0.05
USE_FDR <- FALSE
P_FLOOR <- 1e-5

# --------------------------------------------------------------------------- #
#  Volcano plotting function                                                   #
# --------------------------------------------------------------------------- #
plot_volcano <- function(p_values, betas, idx_sig, xlab, outfile,
                         idx_circle = NULL,
                         circle_color  = "yellow2",
                         circle_size   = 4.2,
                         circle_stroke = 1.2) {

  stopifnot(length(p_values) == length(betas))

  df <- tibble(
    idx  = seq_along(betas),
    beta = as.numeric(betas),
    p    = as.numeric(p_values)
  ) %>%
    mutate(
      p    = ifelse(is.na(p) | !is.finite(p), NA_real_, p),
      beta = ifelse(is.na(beta) | !is.finite(beta), NA_real_, beta),
      p_use     = if (USE_FDR) p.adjust(p, "BH") else p,
      p_use     = pmax(p_use, P_FLOOR, na.rm = FALSE),
      neglog10p = -log10(p_use),
      group     = if_else(idx %in% idx_sig, "Sig", "Not Sig")
    )

  g <- ggplot(df, aes(x = beta, y = neglog10p)) +
    geom_point(data = filter(df, group == "Not Sig"),
               aes(color = group), size = 1.6, alpha = 0.45, na.rm = TRUE) +
    geom_point(data = filter(df, group == "Sig"),
               aes(color = group), size = 1.9, alpha = 0.90, na.rm = TRUE) +
    scale_color_manual(
      values = c("Not Sig" = "grey65", "Sig" = "red"),
      breaks = c("Not Sig", "Sig"), name = NULL
    ) +
    labs(
      x = xlab,
      y = if (USE_FDR) expression(-log[10]("FDR (BH)"))
          else         expression(-log[10]("P-value"))
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

  # Optional circle overlay
 if (!is.null(idx_circle)) {
    g <- g + geom_point(
      data = filter(df, idx %in% as.integer(idx_circle)),
      aes(x = beta, y = neglog10p),
      shape = 21, stroke = circle_stroke, size = circle_size,
      fill = NA, color = circle_color, inherit.aes = FALSE, na.rm = TRUE
    )
  }

  print(g)
  ggsave(outfile, g, width = 8, height = 7, dpi = 300, device = ragg::agg_png)
  invisible(g)
}

# =========================================================================== #
#  CAMRA volcano plots                                                        #
# =========================================================================== #
plot_volcano(selected_values1, PRO_beta[1, ], idx_PRO,
             xlab = "Effect estimate", outfile = "figure/treat_volcano.png")

plot_volcano(selected_values2, PRO_beta[2, ], idx_PRO,
             xlab = "Effect estimate", outfile = "figure/outcome_volcano.png")

# =========================================================================== #
#  LDM-med volcano plots                                                      #
# =========================================================================== #
plot_volcano(p_LDM[1, ], LDM_beta[1, ], idx_LDM,
             xlab = "Effect estimate", outfile = "figure/treat_volcano_LDM.png")

plot_volcano(p_LDM[2, ], LDM_beta[2, ], idx_LDM,
             xlab = "Effect estimate", outfile = "figure/outcome_volcano_LDM.png")

message("Done: volcano plots saved.")

# =============================================================================
# TableS4 — Taxon-level mediation q-values across methods
# =============================================================================
# Requires: count1, res1_prop, res_CRA, res.ldm.med, res_mar, res_hima1,
#           res_multimedia
# Outputs:  LaTeX table printed to console

Q_CUT <- 0.05

# -- Clamp microHIMA q-values at 0.05 floor --
res_hima1[[2]][res_hima1[[2]] < 0.05] <- 0.05

# -- Assemble q-value table --
taxa  <- colnames(count1)

q_tab <- data.frame(
  Species    = as.character(taxa),
  CAMRA      = as.numeric(res1_prop$qval.med),
  CRAmed     = as.numeric(res_CRA[[2]]),
  `LDM-med`  = as.numeric(res.ldm.med$p_med),
  MarZIC     = as.numeric(res_mar$p_med),
  microHIMA  = as.numeric(res_hima1[[2]]),
  multimedia = as.numeric(res_multimedia[[2]]),
  check.names = FALSE
)

# Sort alphabetically by species name
q_tab <- q_tab[order(tolower(trimws(q_tab$Species))), , drop = FALSE]
stopifnot(nrow(q_tab) == length(taxa))

# Keep rows where at least one method is significant
keep    <- apply(q_tab[, -1, drop = FALSE], 1,
                 function(x) any(!is.na(x) & x <= Q_CUT))
q_tab_keep <- q_tab[keep, , drop = FALSE]

# Sort by CAMRA q-value
q_min <- apply(q_tab_keep[, 2, drop = FALSE], 1, min, na.rm = TRUE)
q_tab_keep <- q_tab_keep[order(q_min), , drop = FALSE]

# -- Render LaTeX --
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
    "\\textbf{CAMRA}", "\\textbf{CRAmed}", "\\textbf{LDM-med}",
    "\\textbf{MarZIC}", "\\textbf{microHIMA}", "\\textbf{multimedia}"
  ),
  caption = paste0(
    "\\textbf{Table S6:} Taxon-level mediation $q$-values from different ",
    "methods in the real-data analysis. The table lists taxa that were ",
    "selected by at least one method. We report $q$-values for taxa that ",
    "were declared significant by at least one method at $q \\le ", Q_CUT, "$."
  ),
  label = "tab:S6_qvalues"
) |>
  kableExtra::kable_styling(position = "center", latex_options = "hold_position")

cat(tex)
message("\nDone: Table S6 LaTeX printed to console.")
