#!/usr/bin/env Rscript
# =============================================================================
# permutation_test_acr_outofplace.R
#
# Tests whether spacers preferentially target Anti_CRISPR systems that are
# "out of place" — i.e. located outside the 95% positional distribution
# of Anti_CRISPR systems (below 2.5th or above 97.5th percentile).
#
# Design:
#   1. Define out-of-place threshold from background Anti_CRISPR positions
#   2. For each Anti_CRISPR system, classify as in-place or out-of-place
#   3. Observed: fraction of spacer-HIT systems that are out-of-place
#   4. Null: randomly draw N systems from background, compute fraction
#      out-of-place — repeat 10,000 times
#   5. P-value: fraction of null iterations >= observed fraction
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

set.seed(42)

# ============================================================
# CONFIG
# ============================================================
OVERLAPS_TSV   <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/spacer_antidef_overlap/spacer_antidef_overlaps.tsv"
ANTIDEF_COORDS <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/plasmid_antidef_ridgelines/plasmid_antidef_systems_nt_coords.tsv"
OUTDIR         <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/spacer_antidef_overlap/permutationtest"

N_PERM       <- 10000
TARGET_CLASS <- "Anti_CRISPR"
LOWER_PCTILE <- 0.025   # 2.5th percentile
UPPER_PCTILE <- 0.975   # 97.5th percentile

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load data
# ============================================================
message("Loading data...")

antidef <- read_tsv(ANTIDEF_COORDS, col_types = cols(.default = "c"),
                    show_col_types = FALSE) %>%
  mutate(norm_mid = as.numeric(norm_mid)) %>%
  filter(!is.na(norm_mid), type == TARGET_CLASS)

overlaps <- read_tsv(OVERLAPS_TSV, col_types = cols(.default = "c"),
                     show_col_types = FALSE) %>%
  filter(type == TARGET_CLASS) %>%
  mutate(
    sstart        = as.integer(sstart),
    send          = as.integer(send),
    genome_length = as.integer(genome_length),
    hit_min       = pmin(sstart, send),
    hit_max       = pmax(sstart, send),
    hit_mid_norm  = (hit_min + hit_max) / 2 / genome_length
  ) %>%
  filter(!is.na(hit_mid_norm), hit_mid_norm >= 0, hit_mid_norm <= 1)

# ============================================================
# Define out-of-place threshold from background distribution
# ============================================================
lower_bound <- quantile(antidef$norm_mid, LOWER_PCTILE)
upper_bound <- quantile(antidef$norm_mid, UPPER_PCTILE)

message("\n=== Out-of-place threshold (95% interval) ===")
message("  Lower bound (2.5th pctile): ", round(lower_bound * 100, 2), "%")
message("  Upper bound (97.5th pctile): ", round(upper_bound * 100, 2), "%")

# Classify each Anti_CRISPR system as in-place or out-of-place
antidef <- antidef %>%
  mutate(
    out_of_place = norm_mid < lower_bound | norm_mid > upper_bound
  )

n_total_systems   <- nrow(antidef)
n_outofplace_bg   <- sum(antidef$out_of_place)
pct_outofplace_bg <- n_outofplace_bg / n_total_systems

message("\n=== Background Anti_CRISPR systems ===")
message("  Total systems:        ", n_total_systems)
message("  Out-of-place systems: ", n_outofplace_bg,
        " (", round(pct_outofplace_bg * 100, 1), "%)")

# ============================================================
# Find which systems were hit by spacers (unique sys_id level)
# ============================================================
hit_sys_ids <- unique(overlaps$sys_id)
n_hit       <- length(hit_sys_ids)

antidef_hit <- antidef %>%
  filter(sys_id %in% hit_sys_ids)

n_hit_outofplace   <- sum(antidef_hit$out_of_place)
pct_hit_outofplace <- n_hit_outofplace / nrow(antidef_hit)

message("\n=== Spacer-hit Anti_CRISPR systems ===")
message("  Total systems hit:          ", n_hit)
message("  Out-of-place systems hit:   ", n_hit_outofplace,
        " (", round(pct_hit_outofplace * 100, 1), "%)")
message("  In-place systems hit:       ", n_hit - n_hit_outofplace,
        " (", round((1 - pct_hit_outofplace) * 100, 1), "%)")

# ============================================================
# Permutation test
# Null: randomly draw n_hit systems from background, what fraction
#       are out-of-place?
# One-tailed: is observed fraction > null fraction?
# ============================================================
message("\nRunning permutation test (", N_PERM, " iterations)...")

# Background out-of-place flag vector to resample from
bg_outofplace <- antidef$out_of_place

null_fractions <- replicate(N_PERM, {
  sampled <- sample(bg_outofplace, size = n_hit, replace = TRUE)
  mean(sampled)
})

# One-tailed p-value: fraction of null >= observed
p_value <- mean(null_fractions >= pct_hit_outofplace)

null_mean <- mean(null_fractions)
null_sd   <- sd(null_fractions)
z_score   <- (pct_hit_outofplace - null_mean) / null_sd

message("  Null mean fraction out-of-place: ",
        round(null_mean * 100, 2), "%")
message("  Null SD:                         ",
        round(null_sd * 100, 4), "%")
message("  Observed fraction out-of-place:  ",
        round(pct_hit_outofplace * 100, 2), "%")
message("  Z-score:                         ", round(z_score, 3))
message("  P-value (one-tailed, higher):    ", round(p_value, 4))

# Also run Fisher's exact test
contingency <- matrix(
  c(n_hit_outofplace,
    n_hit - n_hit_outofplace,
    n_outofplace_bg - n_hit_outofplace,
    n_total_systems - n_outofplace_bg - (n_hit - n_hit_outofplace)),
  nrow = 2,
  dimnames = list(
    c("Out-of-place", "In-place"),
    c("Spacer-hit", "Not-hit")
  )
)

#### Alternatives

# One-tailed p-value in the LOWER direction
p_value_lower <- mean(null_fractions <= pct_hit_outofplace)
cat("P-value (one-tailed, lower):", round(p_value_lower, 4), "\n")

fisher_lower <- fisher.test(ct, alternative = "less")
cat("Fisher OR:", round(fisher_lower$estimate, 3), "\n")
cat("Fisher p:", signif(fisher_lower$p.value, 3), "\n")

#### 

message("\n  Contingency table:")
print(contingency)

fisher_result <- fisher.test(contingency, alternative = "greater")
message("  Fisher's exact test (one-tailed):")
message("    OR:      ", round(fisher_result$estimate, 3))
message("    p-value: ", signif(fisher_result$p.value, 3))

# ============================================================
# Plots
# ============================================================

# Plot 1: Null distribution with observed fraction marked
null_df <- tibble(null_frac = null_fractions * 100)

p1 <- ggplot(null_df, aes(x = null_frac)) +
  geom_histogram(bins = 60, fill = "#90A4AE", color = "white",
                 alpha = 0.8, linewidth = 0.2) +
  geom_vline(xintercept = pct_hit_outofplace * 100,
             color = "#B71C1C", linewidth = 1.2) +
  geom_vline(xintercept = null_mean * 100,
             color = "#37474F", linewidth = 0.8, linetype = "dashed") +
  annotate("text",
           x = pct_hit_outofplace * 100 + 0.3,
           y = N_PERM * 0.06,
           label = paste0("Observed\n",
                          round(pct_hit_outofplace * 100, 1), "%"),
           color = "#B71C1C", hjust = 0, size = 3.5, fontface = "bold") +
  scale_x_continuous(
    name = "% of sampled Anti_CRISPR systems that are out-of-place"
  ) +
  scale_y_continuous(name = "Count") +
  labs(
    title    = "Permutation test: do spacers preferentially target out-of-place Anti_CRISPR?",
    subtitle = paste0(N_PERM, " resamplings  |  ",
                      "n systems hit = ", n_hit, "  |  ",
                      "p = ", ifelse(p_value == 0,
                                     paste0("< ", 1/N_PERM), round(p_value, 4)),
                      "  |  Z = ", round(z_score, 2),
                      "  |  Fisher OR = ", round(fisher_result$estimate, 2),
                      "  p = ", signif(fisher_result$p.value, 3))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40")
  )

# Plot 2: Background Anti_CRISPR distribution with out-of-place
# shading and spacer-hit systems marked
p2 <- ggplot(antidef, aes(x = norm_mid)) +
  # Shade out-of-place regions
  annotate("rect", xmin = 0, xmax = lower_bound,
           ymin = -Inf, ymax = Inf,
           fill = "#FFCDD2", alpha = 0.4) +
  annotate("rect", xmin = upper_bound, xmax = 1,
           ymin = -Inf, ymax = Inf,
           fill = "#FFCDD2", alpha = 0.4) +
  geom_density(fill = "#90A4AE", color = "#37474F",
               alpha = 0.4, linewidth = 0.9) +
  # Mark all systems as points
  geom_rug(data = antidef %>% filter(!out_of_place),
           color = "#37474F", alpha = 0.5, linewidth = 0.4) +
  geom_rug(data = antidef %>% filter(out_of_place),
           color = "#B71C1C", alpha = 0.8, linewidth = 0.6) +
  # Mark spacer-hit systems
  geom_point(data = antidef_hit,
             aes(x = norm_mid, y = 0,
                 color = out_of_place, shape = out_of_place),
             size = 3, alpha = 0.8, position = position_jitter(height = 0.05)) +
  geom_vline(xintercept = c(lower_bound, upper_bound),
             linetype = "dashed", color = "#B71C1C", linewidth = 0.8) +
  scale_x_continuous(
    name   = "Normalized genome position",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1), breaks = seq(0, 1, 0.1)
  ) +
  scale_y_continuous(name = "Density") +
  scale_color_manual(
    values = c("FALSE" = "#37474F", "TRUE" = "#B71C1C"),
    labels = c("FALSE" = "In-place (hit)", "TRUE" = "Out-of-place (hit)"),
    name = "Spacer-hit systems"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "In-place (hit)", "TRUE" = "Out-of-place (hit)"),
    name = "Spacer-hit systems"
  ) +
  annotate("text", x = lower_bound / 2, y = Inf,
           label = "Out of\nplace", vjust = 1.5, size = 3,
           color = "#B71C1C", fontface = "italic") +
  annotate("text", x = (upper_bound + 1) / 2, y = Inf,
           label = "Out of\nplace", vjust = 1.5, size = 3,
           color = "#B71C1C", fontface = "italic") +
  labs(
    title    = "Anti_CRISPR system positions with out-of-place threshold",
    subtitle = paste0("Red shading = outside 95% interval  |  ",
                      "Triangles = out-of-place systems hit by spacers  |  ",
                      n_hit_outofplace, "/", n_hit,
                      " (", round(pct_hit_outofplace * 100, 1),
                      "%) spacer-hit systems are out-of-place")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    legend.position  = "top"
  )

# ============================================================
# Save
# ============================================================
ggsave(file.path(OUTDIR, "acr_outofplace_permutation_null.pdf"),
       p1, width = 10, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "acr_outofplace_distribution.pdf"),
       p2, width = 10, height = 6, device = "pdf")

message("\nWritten: acr_outofplace_permutation_null.pdf")
message("Written: acr_outofplace_distribution.pdf")

# Save results
results <- tibble(
  lower_bound_pctile     = lower_bound,
  upper_bound_pctile     = upper_bound,
  n_total_systems        = n_total_systems,
  n_outofplace_bg        = n_outofplace_bg,
  pct_outofplace_bg      = pct_outofplace_bg,
  n_systems_hit          = n_hit,
  n_hit_outofplace       = n_hit_outofplace,
  pct_hit_outofplace     = pct_hit_outofplace,
  null_mean              = null_mean,
  null_sd                = null_sd,
  z_score                = z_score,
  p_permutation          = ifelse(p_value == 0, 1/N_PERM, p_value),
  fisher_OR              = fisher_result$estimate,
  fisher_p               = fisher_result$p.value,
  fisher_p_lower         = fisher_lower$p.value,
  fisher_OR_lower        = fisher_lower$estimate
)

write_tsv(results, file.path(OUTDIR, "acr_outofplace_results.tsv"))
message("Written: acr_outofplace_results.tsv")

cat("\n=== Final Results ===\n")
cat("  95% interval: [",
    round(lower_bound * 100, 1), "%, ",
    round(upper_bound * 100, 1), "%]\n", sep = "")
cat("  Background out-of-place:        ",
    round(pct_outofplace_bg * 100, 1), "%  (", n_outofplace_bg,
    "/", n_total_systems, ")\n")
cat("  Spacer-hit out-of-place:        ",
    round(pct_hit_outofplace * 100, 1), "%  (", n_hit_outofplace,
    "/", n_hit, ")\n")
cat("  Z-score:                        ", round(z_score, 3), "\n")
cat("  Permutation p (one-tailed):     ",
    ifelse(p_value == 0, paste0("< ", 1/N_PERM), round(p_value, 4)), "\n")
cat("  Fisher's exact OR:              ", round(fisher_result$estimate, 3), "\n")
cat("  Fisher's exact p (one-tailed):  ", signif(fisher_result$p.value, 3), "\n")
cat("  Fisher's exact OR (lower):      ", round(fisher_lower$estimate, 3), "\n")
cat("  Fisher's exact p (lower; one-tailed):   ", signif(fisher_lower$p.value, 3), "\n")
message("Done.")



