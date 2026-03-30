#!/usr/bin/env Rscript
# =============================================================================
# antidefense_enrichment_analysis.R
#
# Tests whether CRISPR spacers preferentially target Anti_CRISPR systems
# over other anti-defense system classes, in both phages and plasmids.
#
# Tests performed:
#   - Fisher's exact: Anti_CRISPR vs all others combined
#   - Chi-square: across all classes with n >= 5
#   - Binomial: per-class enrichment vs genome-wide background rate
#   - Combined dot plot of hit rates by class (phage + plasmid)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ============================================================
# CONFIG
# ============================================================
PHAGE_OVERLAPS   <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/spacer_antidef_overlap/spacer_antidef_overlaps.tsv"
PHAGE_ANTIDEF    <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/antidef_ridgelines/derep_phage_antidef_systems_nt_coords.tsv"
PLASMID_OVERLAPS <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/spacer_antidef_overlap/spacer_antidef_overlaps.tsv"
PLASMID_ANTIDEF  <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/plasmid_antidef_ridgelines/plasmid_antidef_systems_nt_coords.tsv"
OUTDIR           <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/selfishgene_stats/antidefense_enrichment"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper: compute hit rates per class
# ============================================================
compute_enrichment <- function(overlaps, antidef) {
  hit_sys_ids <- unique(overlaps$sys_id)
  antidef %>%
    group_by(type) %>%
    summarise(
      n_total = n(),
      n_hit   = sum(sys_id %in% hit_sys_ids),
      .groups = "drop"
    ) %>%
    mutate(
      n_not_hit = n_total - n_hit,
      hit_rate  = n_hit / n_total
    ) %>%
    arrange(desc(hit_rate))
}

# ============================================================
# Helper: run Fisher's exact (Anti_CRISPR vs all others)
# ============================================================
run_fisher <- function(enrichment, label) {
  acr   <- enrichment %>% filter(type == "Anti_CRISPR")
  other <- enrichment %>% filter(type != "Anti_CRISPR") %>%
    summarise(n_hit = sum(n_hit), n_not_hit = sum(n_not_hit))

  ct <- matrix(
    c(acr$n_hit, acr$n_not_hit, other$n_hit, other$n_not_hit),
    nrow = 2,
    dimnames = list(c("Hit", "Not hit"), c("Anti_CRISPR", "Other"))
  )

  cat("\n=== Fisher's exact test —", label, "===\n")
  cat("Contingency table:\n")
  print(ct)

  result <- fisher.test(ct, alternative = "greater")
  cat("OR:      ", round(result$estimate, 3), "\n")
  cat("95% CI:  [", round(result$conf.int[1], 1), ", Inf]\n")
  cat("p-value: ", signif(result$p.value, 3), "\n")

  invisible(list(ct = ct, result = result))
}

# ============================================================
# Helper: chi-square across classes with n_total >= 5
# ============================================================
run_chisq <- function(enrichment, label) {
  cat("\n=== Chi-square test —", label, "(n_total >= 5) ===\n")
  mat <- enrichment %>%
    filter(n_total >= 5) %>%
    select(type, n_hit, n_not_hit) %>%
    column_to_rownames("type")
  print(chisq.test(mat))
}

# ============================================================
# Helper: per-class binomial test vs background rate
# ============================================================
run_binomial <- function(enrichment, label) {
  # Background rate = overall hit rate across ALL systems
  bg_rate <- sum(enrichment$n_hit) / sum(enrichment$n_total)
  cat("\n=== Binomial tests —", label,
      "(background rate =", round(bg_rate * 100, 2), "%) ===\n")

  enrichment %>%
    filter(n_total >= 3) %>%
    rowwise() %>%
    mutate(
      binom_p   = binom.test(n_hit, n_total, p = bg_rate,
                             alternative = "greater")$p.value,
      binom_p   = signif(binom_p, 3),
      enrichment_fold = round(hit_rate / bg_rate, 2)
    ) %>%
    ungroup() %>%
    select(type, n_total, n_hit, hit_rate, enrichment_fold, binom_p) %>%
    mutate(hit_rate = round(hit_rate * 100, 1)) %>%
    print(n = Inf)
}

# ============================================================
# Load data
# ============================================================
message("Loading phage data...")
phage_overlaps <- read_tsv(PHAGE_OVERLAPS, col_types = cols(.default = "c"),
                            show_col_types = FALSE)
phage_antidef  <- read_tsv(PHAGE_ANTIDEF,  col_types = cols(.default = "c"),
                            show_col_types = FALSE) %>%
  mutate(norm_mid = as.numeric(norm_mid))

message("Loading plasmid data...")
plasmid_overlaps <- read_tsv(PLASMID_OVERLAPS, col_types = cols(.default = "c"),
                              show_col_types = FALSE)
plasmid_antidef  <- read_tsv(PLASMID_ANTIDEF,  col_types = cols(.default = "c"),
                              show_col_types = FALSE) %>%
  mutate(norm_mid = as.numeric(norm_mid))

# ============================================================
# Compute enrichment
# ============================================================
phage_enrichment   <- compute_enrichment(phage_overlaps,   phage_antidef)
plasmid_enrichment <- compute_enrichment(plasmid_overlaps, plasmid_antidef)

cat("\n=== Phage hit rates by class ===\n")
print(phage_enrichment, n = Inf)

cat("\n=== Plasmid hit rates by class ===\n")
print(plasmid_enrichment, n = Inf)

# ============================================================
# Statistical tests
# ============================================================
phage_fisher   <- run_fisher(phage_enrichment,   "Phage")
plasmid_fisher <- run_fisher(plasmid_enrichment, "Plasmid")

run_chisq(phage_enrichment,   "Phage")
run_chisq(plasmid_enrichment, "Plasmid")

run_binomial(phage_enrichment,   "Phage")
run_binomial(plasmid_enrichment, "Plasmid")

# ============================================================
# Plot 1: Combined dot plot — hit rate by class, phage and plasmid
# ============================================================
plot_df <- bind_rows(
  phage_enrichment   %>% mutate(dataset = "Phage"),
  plasmid_enrichment %>% mutate(dataset = "Plasmid")
) %>%
  filter(n_total >= 3) %>%
  mutate(
    type = str_replace(type, "Anti_", "Anti-"),
    type = fct_reorder(type, hit_rate, .fun = max)
  )

p <- ggplot(plot_df,
            aes(x = hit_rate, y = type,
                color = dataset, size = n_total)) +
  geom_point(alpha = 0.85) +
  geom_vline(
    data = plot_df %>%
      group_by(dataset) %>%
      summarise(bg = sum(n_hit) / sum(n_total), .groups = "drop"),
    aes(xintercept = bg, color = dataset),
    linetype = "dashed", linewidth = 0.7, alpha = 0.6
  ) +
  scale_x_continuous(
    name   = "Fraction of systems hit by spacers",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_size_continuous(
    name   = "Total systems",
    range  = c(2, 10),
    breaks = c(5, 20, 50, 100, 200)
  ) +
  scale_color_manual(
    values = c("Phage" = "#1565C0", "Plasmid" = "#C62828"),
    name   = "Dataset"
  ) +
  labs(
    y        = NULL,
    title    = "CRISPR spacer targeting rate by anti-defense system class",
    subtitle = "Dashed lines = background hit rate per dataset  |  Point size = total systems"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    legend.position  = "right"
  )

ggsave(file.path(OUTDIR, "antidefense_hit_rate_by_class.pdf"),
       p, width = 9, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "antidefense_hit_rate_by_class.png"),
       p, width = 9, height = 6, dpi = 300, bg = "white")
message("Written: antidefense_hit_rate_by_class.pdf/.png")



# ============================================================
# Save results tables
# ============================================================
write_tsv(phage_enrichment   %>% mutate(dataset = "Phage"),
          file.path(OUTDIR, "phage_enrichment_by_class.tsv"))
write_tsv(plasmid_enrichment %>% mutate(dataset = "Plasmid"),
          file.path(OUTDIR, "plasmid_enrichment_by_class.tsv"))

# Fisher results summary
fisher_summary <- tibble(
  dataset     = c("Phage", "Plasmid"),
  OR          = c(phage_fisher$result$estimate,
                  plasmid_fisher$result$estimate),
  CI_lower    = c(phage_fisher$result$conf.int[1],
                  plasmid_fisher$result$conf.int[1]),
  p_value     = c(phage_fisher$result$p.value,
                  plasmid_fisher$result$p.value)
)
write_tsv(fisher_summary, file.path(OUTDIR, "fisher_acr_vs_other.tsv"))
message("Written: fisher_acr_vs_other.tsv")

cat("\n=== Summary ===\n")
cat("Phage   — Anti_CRISPR hit rate: ",
    round(phage_enrichment %>%
            filter(type == "Anti_CRISPR") %>% pull(hit_rate) * 100, 1), "%\n")
cat("Plasmid — Anti_CRISPR hit rate: ",
    round(plasmid_enrichment %>%
            filter(type == "Anti_CRISPR") %>% pull(hit_rate) * 100, 1), "%\n")
cat("Phage   — Fisher OR: ",
    round(phage_fisher$result$estimate, 1),
    "  p =", signif(phage_fisher$result$p.value, 3), "\n")
cat("Plasmid — Fisher OR: ",
    round(plasmid_fisher$result$estimate, 1),
    "  p =", signif(plasmid_fisher$result$p.value, 3), "\n")
message("Done.")


## ADDITIONAL PLOTS
# Plot 2:
# Wilson confidence intervals for proportions
wilson_ci <- function(n_hit, n_total, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- n_hit / n_total
  center <- (p + z^2 / (2 * n_total)) / (1 + z^2 / n_total)
  margin <- (z * sqrt(p * (1 - p) / n_total + z^2 / (4 * n_total^2))) /
    (1 + z^2 / n_total)
  tibble(estimate = center, lower = center - margin, upper = center + margin)
}

bar_df <- bind_rows(
  phage_enrichment   %>% mutate(dataset = "Phage"),
  plasmid_enrichment %>% mutate(dataset = "Plasmid")
) %>%
  filter(n_total >= 3) %>%
  rowwise() %>%
  mutate(wilson_ci(n_hit, n_total)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type, estimate, .fun = max))

p_bar <- ggplot(bar_df, aes(x = type, y = estimate,
                            fill = dataset, group = dataset)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.8), width = 0.25,
                linewidth = 0.6) +
  scale_y_continuous(name   = "Fraction of systems hit",
                     labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05)) +
  scale_fill_manual(values = c("Phage" = "#1565C0", "Plasmid" = "#C62828"),
                    name = NULL) +
  labs(x = NULL,
       title    = "CRISPR spacer targeting rate by anti-defense class",
       subtitle = "Wilson 95% confidence intervals  |  n >= 3 systems") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    legend.position  = "top"
  )

ggsave(file.path(OUTDIR, "antidefense_hit_rate_bar.pdf"),
       p_bar, width = 10, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "antidefense_hit_rate_bar.png"),
       p_bar, width = 10, height = 6, dpi = 300, bg = "white")
message("Written: antidefense_hit_rate_bar.pdf/.png")

# Plot 3
lollipop_df <- bind_rows(
phage_enrichment   %>% mutate(dataset = "Phage"),
plasmid_enrichment %>% mutate(dataset = "Plasmid")
) %>%
  filter(n_total >= 3) %>%
  mutate(type = fct_reorder(type, hit_rate, .fun = max))

p_lollipop <- ggplot(lollipop_df,
                     aes(x = hit_rate, y = type, color = dataset)) +
  geom_segment(aes(x = 0, xend = hit_rate,
                   y = type, yend = type),
               position = position_dodge(0.6),
               linewidth = 0.8, alpha = 0.6) +
  geom_point(position = position_dodge(0.6),
             size = 4, alpha = 0.9) +
  geom_text(aes(label = paste0(n_hit, "/", n_total)),
            position = position_dodge(0.6),
            hjust = -0.3, size = 2.8) +
  scale_x_continuous(name   = "Fraction of systems hit",
                     labels = percent_format(accuracy = 1),
                     limits = c(0, 1.15)) +
  scale_color_manual(values = c("Phage" = "#1565C0", "Plasmid" = "#C62828"),
                     name = NULL) +
  labs(y = NULL,
       title    = "CRISPR spacer targeting rate by anti-defense class",
       subtitle = "Labels show hits/total systems") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(face = "bold"),
    legend.position  = "top"
  )

ggsave(file.path(OUTDIR, "antidefense_hit_rate_lollipop.pdf"),
       p_lollipop, width = 9, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "antidefense_hit_rate_lollipop.png"),
       p_lollipop, width = 9, height = 6, dpi = 300, bg = "white")
message("Written: antidefense_hit_rate_lollipop.pdf/.png")

# Plot 4
heatmap_df <- bind_rows(
  phage_enrichment   %>% mutate(dataset = "Phage"),
  plasmid_enrichment %>% mutate(dataset = "Plasmid")
) %>%
  mutate(
    type    = fct_reorder(type, hit_rate, .fun = max),
    dataset = factor(dataset, levels = c("Phage", "Plasmid")),
    label   = paste0(n_hit, "/", n_total)
  )

p_heatmap <- ggplot(heatmap_df, aes(x = dataset, y = type, fill = hit_rate)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label,
                color = hit_rate > 0.5),
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low      = "#F5F5F5",
    mid      = "#FFAB40",
    high     = "#B71C1C",
    midpoint = 0.5,
    limits   = c(0, 1),
    labels   = percent_format(accuracy = 1),
    name     = "Hit rate"
  ) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"),
                     guide = "none") +
  scale_x_discrete(name = NULL, expand = c(0, 0)) +
  scale_y_discrete(name = NULL, expand = c(0, 0)) +
  labs(
    title    = "CRISPR spacer targeting rate by anti-defense class",
    subtitle = "Cell labels show hits / total systems"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(face = "bold"),
    axis.ticks   = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(OUTDIR, "antidefense_hit_rate_heatmap.pdf"),
       p_heatmap, width = 7,
       height = max(4, 0.4 * n_distinct(heatmap_df$type) + 2),
       device = "pdf")
ggsave(file.path(OUTDIR, "antidefense_hit_rate_heatmap.png"),
       p_heatmap, width = 7,
       height = max(4, 0.4 * n_distinct(heatmap_df$type) + 2),
       dpi = 300, bg = "white")
message("Written: antidefense_hit_rate_heatmap.pdf/.png")

# Plot 5
fold_df <- bind_rows(
  phage_enrichment %>% mutate(
    dataset = "Phage",
    bg_rate = sum(n_hit) / sum(n_total)
  ),
  plasmid_enrichment %>% mutate(
    dataset = "Plasmid",
    bg_rate = sum(n_hit) / sum(n_total)
  )
) %>%
  filter(n_total >= 3) %>%
  mutate(
    fold      = hit_rate / bg_rate,
    type      = fct_reorder(type, fold, .fun = max),
    log2_fold = log2(fold + 0.01)   # pseudocount for 0 hit rates
  )

p_fold <- ggplot(fold_df, aes(x = log2_fold, y = type, fill = log2_fold)) +
  geom_col(alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "grey40",
             linetype = "dashed") +
  facet_wrap(~ dataset, scales = "free_x") +
  scale_fill_gradient2(low  = "#1565C0", mid = "#F5F5F5",
                       high = "#B71C1C", midpoint = 0,
                       guide = "none") +
  scale_x_continuous(name = "log2 fold enrichment over background") +
  labs(y = NULL,
       title    = "Anti-defense class enrichment over background spacer hit rate",
       subtitle = "Background = overall hit rate across all systems per dataset") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

ggsave(file.path(OUTDIR, "antidefense_fold_enrichment.pdf"),
       p_fold, width = 11, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "antidefense_fold_enrichment.png"),
       p_fold, width = 11, height = 6, dpi = 300, bg = "white")
message("Written: antidefense_fold_enrichment.pdf/.png")

