#!/usr/bin/env Rscript
# =============================================================================
# mystery_orientation_crossref.R
#
# Cross-references mystery-oriented genomes (dnaapler found no TerL)
# against out-of-place Anti_CRISPR systems to test whether positional
# outliers are orientation artifacts.
# =============================================================================

library(tidyverse)

# ============================================================
# CONFIG
# ============================================================
MYSTERY_GENOMES <- "/tmp/mystery_oriented_genomes.txt"
ANTIDEF_COORDS  <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/antidef_ridgelines/derep_phage_antidef_systems_nt_coords.tsv"
OUTDIR          <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/spacer_antidef_overlap"

LOWER_PCTILE <- 0.025
UPPER_PCTILE <- 0.975
TARGET_CLASS <- "Anti_CRISPR"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load data
# ============================================================
mystery_genomes <- read_lines(MYSTERY_GENOMES) %>% str_trim()
message("Mystery-oriented genomes: ", length(mystery_genomes))

antidef <- read_tsv(ANTIDEF_COORDS, col_types = cols(.default = "c"),
                    show_col_types = FALSE) %>%
  mutate(norm_mid = as.numeric(norm_mid)) %>%
  filter(!is.na(norm_mid), type == TARGET_CLASS)

message("Anti_CRISPR systems: ", nrow(antidef))

# ============================================================
# Classify systems
# ============================================================
lower_bound <- quantile(antidef$norm_mid, LOWER_PCTILE)
upper_bound <- quantile(antidef$norm_mid, UPPER_PCTILE)

message("95% interval: [", round(lower_bound*100,1), "%, ",
        round(upper_bound*100,1), "%]")

antidef <- antidef %>%
  mutate(
    out_of_place   = norm_mid < lower_bound | norm_mid > upper_bound,
    mystery_orient = replicon %in% mystery_genomes
  )

# ============================================================
# Results
# ============================================================
cat("\n=== Cross-tabulation ===\n")
print(table(out_of_place   = antidef$out_of_place,
            mystery_orient = antidef$mystery_orient))

cat("\n=== Fraction of out-of-place systems from mystery-oriented genomes ===\n")
antidef %>%
  filter(out_of_place) %>%
  count(mystery_orient) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print()

cat("\n=== Out-of-place systems from mystery-oriented genomes ===\n")
antidef %>%
  filter(out_of_place, mystery_orient) %>%
  select(sys_id, replicon, subtype, norm_mid) %>%
  print(n = Inf)

cat("\n=== Out-of-place systems NOT from mystery-oriented genomes ===\n")
antidef %>%
  filter(out_of_place, !mystery_orient) %>%
  select(sys_id, replicon, subtype, norm_mid) %>%
  print(n = Inf)

# Fisher's exact test — are out-of-place systems enriched in mystery genomes?
ct <- table(out_of_place   = antidef$out_of_place,
            mystery_orient = antidef$mystery_orient)
fisher_result <- fisher.test(ct, alternative = "greater")
cat("\n=== Fisher's exact test (out-of-place ~ mystery orientation) ===\n")
cat("  OR:      ", round(fisher_result$estimate, 3), "\n")
cat("  p-value: ", signif(fisher_result$p.value, 3), "\n")

write_tsv(
  antidef %>% select(sys_id, replicon, subtype, norm_mid,
                     out_of_place, mystery_orient),
  file.path(OUTDIR, "acr_mystery_orientation_crossref.tsv")
)
message("\nWritten: acr_mystery_orientation_crossref.tsv")
message("Done.")