#!/usr/bin/env Rscript
# =============================================================================
# spacer_antidef_overlap_derep.R
#
# Overlaps spacer BLAST hits against dereplicated Pa phages with
# DefenseFinder anti-defense system coordinates. Plots:
#   1. Distribution of ALL anti-defense systems (background)
#   2. Distribution of spacer-targeted anti-defense systems (signal)
#   3. Overlay comparison
#   4. Enrichment by subtype
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(scales)
  library(patchwork)
})

# ============================================================
# CONFIG
# ============================================================
BLAST_TSV    <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/blast_results/spacers_vs_orient_plasmids.blast6.tsv"
ANTIDEF_COORDS <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/plasmid_antidef_ridgelines/plasmid_antidef_systems_nt_coords.tsv"
OUTDIR       <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/spacer_antidef_overlap"

PIDENT_MIN   <- 90    # % identity cutoff
LEN_MIN      <- 20    # minimum alignment length (nt)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load BLAST results
# ============================================================
message("Loading BLAST results...")
blast_cols <- c("qseqid", "sseqid", "sstart", "send", "sstrand",
                "pident", "length", "mismatch", "gapopen", "qlen")

blast <- read_tsv(BLAST_TSV, col_names = blast_cols,
                  col_types = cols(.default = "c"),
                  show_col_types = FALSE) %>%
  mutate(
    pident  = as.numeric(pident),
    length  = as.integer(length),
    sstart  = as.integer(sstart),
    send    = as.integer(send),
    hit_min = pmin(sstart, send),
    hit_max = pmax(sstart, send)
  ) %>%
  filter(pident >= PIDENT_MIN, length >= LEN_MIN)

message("  Hits after filter: ", nrow(blast))
message("  Unique spacers: ", n_distinct(blast$qseqid))
message("  Unique phage accessions: ", n_distinct(blast$sseqid))

# ============================================================
# Load anti-defense system coordinates
# ============================================================
message("\nLoading anti-defense system coordinates...")
antidef <- read_tsv(ANTIDEF_COORDS, col_types = cols(.default = "c"),
                   show_col_types = FALSE) %>%
  mutate(
    sys_nt_begin  = as.integer(sys_nt_begin),
    sys_nt_end    = as.integer(sys_nt_end),
    genome_length = as.integer(genome_length),
    norm_mid      = as.numeric(norm_mid)
  )

message("  Anti-defense systems: ", nrow(antidef))
message("  Unique phages with anti-defense: ", n_distinct(antidef$replicon))

# ============================================================
# Overlap: flag BLAST hits that land within anti-defense loci
# ============================================================
message("\nFlagging spacer hits overlapping anti-defense loci...")

# For each anti-defense system, find BLAST hits where:
#   sseqid == replicon AND hit overlaps [sys_nt_begin, sys_nt_end]
overlapping_hits <- blast %>%
  inner_join(
    antidef %>% select(sys_id, replicon, type, subtype,
                       sys_nt_begin, sys_nt_end, genome_length, norm_mid),
    by = c("sseqid" = "replicon"),
    relationship = "many-to-many"
  ) %>%
  filter(
    hit_min <= sys_nt_end,
    hit_max >= sys_nt_begin
  )

message("  Spacer hits overlapping anti-defense loci: ", nrow(overlapping_hits))
message("  Unique spacers hitting anti-defense: ",
        n_distinct(overlapping_hits$qseqid))
message("  Unique anti-defense systems hit: ",
        n_distinct(overlapping_hits$sys_id))
message("  Unique phages: ", n_distinct(overlapping_hits$sseqid))

# Anti-defense systems that are hit by at least one spacer
antidef_hit <- antidef %>%
  filter(sys_id %in% overlapping_hits$sys_id)

message("  Anti-defense systems with >= 1 spacer hit: ", nrow(antidef_hit),
        " / ", nrow(antidef), " total (",
        round(100 * nrow(antidef_hit) / nrow(antidef), 1), "%)")

# ============================================================
# Save overlap table
# ============================================================
write_tsv(overlapping_hits, file.path(OUTDIR, "spacer_antidef_overlaps.tsv"))
write_tsv(antidef_hit,      file.path(OUTDIR, "antidef_systems_hit_by_spacers.tsv"))
message("\nWritten: spacer_antidef_overlaps.tsv")
message("Written: antidef_systems_hit_by_spacers.tsv")

# ============================================================
# Plots
# ============================================================
theme_clean <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      axis.title       = element_text(size = 11),
      legend.position  = "right"
    )
}

x_scale <- scale_x_continuous(
  name   = "Normalized genome position",
  labels = percent_format(accuracy = 1),
  limits = c(0, 1), breaks = seq(0, 1, 0.1)
)

# Plot 1: Background — all anti-defense system positions
p1 <- ggplot(antidef, aes(x = norm_mid)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "#90A4AE", color = "white",
                 alpha = 0.8, linewidth = 0.2) +
  geom_density(color = "#37474F", linewidth = 1.0) +
  x_scale +
  scale_y_continuous(name = "Density") +
  labs(
    title    = "Background: all anti-defense system positions",
    subtitle = paste0("n = ", nrow(antidef), " systems  |  ",
                      n_distinct(antidef$replicon), " phage genomes")
  ) +
  theme_clean()

# Plot 2: Spacer-targeted anti-defense system positions
p2 <- ggplot(antidef_hit, aes(x = norm_mid)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "#EF5350", color = "white",
                 alpha = 0.8, linewidth = 0.2) +
  geom_density(color = "#B71C1C", linewidth = 1.0) +
  x_scale +
  scale_y_continuous(name = "Density") +
  labs(
    title    = "Spacer-targeted anti-defense system positions",
    subtitle = paste0("n = ", nrow(antidef_hit), " systems hit  |  ",
                      n_distinct(antidef_hit$replicon), " phage genomes")
  ) +
  theme_clean()

# Plot 3: Overlay — background vs spacer-targeted
overlay_df <- bind_rows(
  antidef     %>% mutate(group = "All anti-defense systems"),
  antidef_hit %>% mutate(group = "Spacer-targeted systems")
)

p3 <- ggplot(overlay_df, aes(x = norm_mid, fill = group, color = group)) +
  geom_density(alpha = 0.25, linewidth = 0.9) +
  x_scale +
  scale_y_continuous(name = "Density") +
  scale_fill_manual(values  = c("All anti-defense systems" = "#90A4AE",
                                "Spacer-targeted systems"  = "#EF5350"),
                    name = NULL) +
  scale_color_manual(values = c("All anti-defense systems" = "#37474F",
                                "Spacer-targeted systems"  = "#B71C1C"),
                     name = NULL) +
  labs(
    title    = "Anti-defense system positions: background vs spacer-targeted",
    subtitle = paste0("Background n=", nrow(antidef),
                      "  |  Spacer-targeted n=", nrow(antidef_hit))
  ) +
  theme_clean() +
  theme(legend.position = "top")

# Plot 4: Enrichment by subtype — fraction of each subtype hit by spacers
enrichment_df <- antidef %>%
  group_by(subtype) %>%
  summarise(
    n_total = n(),
    n_hit   = sum(sys_id %in% antidef_hit$sys_id),
    .groups = "drop"
  ) %>%
  mutate(
    pct_hit  = n_hit / n_total,
    subtype  = fct_reorder(subtype, pct_hit)
  ) %>%
  filter(n_total >= 2)

p4 <- ggplot(enrichment_df, aes(x = pct_hit, y = subtype,
                                 fill = pct_hit)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = paste0(n_hit, "/", n_total)),
            hjust = -0.15, size = 3.2) +
  scale_x_continuous(
    name   = "Fraction of systems hit by spacers",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.15)
  ) +
  scale_fill_gradient(low = "#FFCDD2", high = "#B71C1C", guide = "none") +
  labs(y = NULL,
       title    = "Fraction of anti-defense systems targeted by spacers",
       subtitle = "By subtype  |  n >= 2 systems") +
  theme_clean() +
  theme(legend.position = "none")

# ============================================================
# Save plots
# ============================================================
ggsave(file.path(OUTDIR, "antidef_background_positions.pdf"),
       p1, width = 9, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "antidef_spacer_targeted_positions.pdf"),
       p2, width = 9, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "antidef_overlay_background_vs_targeted.pdf"),
       p3, width = 9, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "antidef_enrichment_by_subtype.pdf"),
       p4, width = 9, height = max(5, 0.3 * nrow(enrichment_df) + 2),
       device = "pdf")

message("Written: antidef_background_positions.pdf")
message("Written: antidef_spacer_targeted_positions.pdf")
message("Written: antidef_overlay_background_vs_targeted.pdf")
message("Written: antidef_enrichment_by_subtype.pdf")

# ============================================================
# Summary stats
# ============================================================
cat("\n=== Summary ===\n")
cat("  BLAST hits (filtered):          ", nrow(blast), "\n")
cat("  Anti-defense systems (total):   ", nrow(antidef), "\n")
cat("  Anti-defense systems hit:       ", nrow(antidef_hit), "\n")
cat("  Fraction hit:                   ",
    round(100 * nrow(antidef_hit) / nrow(antidef), 1), "%\n")
cat("  Unique spacers hitting antidef: ",
    n_distinct(overlapping_hits$qseqid), "\n")
cat("  Output dir:                     ", normalizePath(OUTDIR), "\n")
message("Done.")
