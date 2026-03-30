#!/usr/bin/env Rscript
# =============================================================================
# plot_terL_positions.R
#
# Concatenates hmmsearch --tblout results for TerL HMM profiles,
# extracts genomic coordinates from the Prodigal annotation embedded
# in each hit line, normalizes by genome length, and plots position
# distribution across normalized phage genome length.
#
# Input: directory of *_terl_hits.tsv files from run_terl_hmmsearch.sh
# =============================================================================

library(tidyverse)
library(scales)
library(Biostrings)

# ============================================================
# CONFIG
# ============================================================
RESULTS_DIR <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/terl_hmmsearch_results_oriented"
OUTDIR       <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/terl_plots"
MIN_GENOME   <- 2000    # nt
EVALUE_CUTOFF <- 1e-5

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parse hmmsearch tblout format
# Columns (space-delimited, comments start with #):
#  1: target name (e.g. AB008550.1_3)
#  2: target accession
#  3: query name (HMM profile name)
#  4: query accession (e.g. PF04466.22)
#  5: E-value (full sequence)
#  6: score (full sequence)
#  ... rest of fields ...
#  rest of line after fields: Prodigal annotation
#    # nt_start # nt_end # strand # ID=...
# ============================================================
parse_tblout <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  # Keep only non-comment, non-empty lines
  hits  <- lines[!startsWith(lines, "#") & nchar(trimws(lines)) > 0]
  if (length(hits) == 0) return(NULL)

  map_dfr(hits, function(line) {
    # Split on whitespace for the first 18 fields, rest is description
    fields <- str_split(line, "\\s+")[[1]]
    if (length(fields) < 6) return(NULL)

    target_name <- fields[1]
    query_name  <- fields[3]
    query_acc   <- fields[4]
    evalue      <- as.numeric(fields[5])
    score       <- as.numeric(fields[6])

    # Extract Prodigal coords from description: # nt_start # nt_end # strand
    nt_start <- as.integer(str_match(line, "# (\\d+) # \\d+ # -?1")[, 2])
    nt_end   <- as.integer(str_match(line, "# \\d+ # (\\d+) # -?1")[, 2])

    # Extract accession from target_name: everything before last _DIGITS
    accession <- str_extract(target_name, "^.+?(?=_\\d+$)")

    tibble(
      target_name = target_name,
      accession   = accession,
      hmm_profile = query_name,
      hmm_acc     = query_acc,
      evalue      = evalue,
      score       = score,
      nt_start    = nt_start,
      nt_end      = nt_end
    )
  })
}

# ============================================================
# Load all tblout files
# ============================================================
message("Loading hmmsearch results from: ", RESULTS_DIR)

tsv_files <- list.files(RESULTS_DIR, pattern = "_terl_hits\\.tsv$",
                        full.names = TRUE)
message("  Files found: ", length(tsv_files))

all_hits <- map_dfr(tsv_files, function(f) {
  tryCatch(parse_tblout(f), error = function(e) NULL)
}) %>%
  filter(!is.na(evalue), evalue <= EVALUE_CUTOFF)

message("  Total hits after E-value filter (<= ", EVALUE_CUTOFF, "): ",
        nrow(all_hits))
message("  Unique accessions with TerL hit: ",
        n_distinct(all_hits$accession))
message("  HMM profiles matched: ",
        paste(sort(unique(all_hits$hmm_profile)), collapse = ", "))


DF_RESULTS <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/terLanalysis/defensefinder_results"

genome_lengths <- map_dfr(
  list.files(DF_RESULTS, pattern = "\\.prt$", recursive = TRUE,
             full.names = TRUE),
  function(f) {
    lines <- readLines(f, warn = FALSE)
    headers <- lines[startsWith(lines, ">")]
    if (length(headers) == 0) return(NULL)
    ends <- as.integer(str_match(headers, "# \\d+ # (\\d+) # -?1")[, 2])
    acc  <- basename(dirname(f))
    tibble(accession = acc, genome_length = max(ends, na.rm = TRUE))
  }
) %>% filter(is.finite(genome_length), genome_length >= MIN_GENOME)

# ============================================================
# Keep best hit per genome (lowest e-value) and normalize
# ============================================================
terl_positions <- all_hits %>%
  left_join(genome_lengths, by = "accession") %>%
  filter(!is.na(genome_length)) %>%
  group_by(accession) %>%
  slice_min(evalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    norm_mid = (nt_start + nt_end) / 2 / genome_length
  ) %>%
  filter(norm_mid >= 0, norm_mid <= 1)

message("\nGenomes with TerL position: ", nrow(terl_positions))
message("norm_mid summary:")
print(summary(terl_positions$norm_mid))

# ============================================================
# Save concatenated results
# ============================================================
write_tsv(terl_positions, file.path(OUTDIR, "terl_positions.tsv"))
message("\nWritten: terl_positions.tsv")

# ============================================================
# Plots
# ============================================================
theme_terl <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      axis.title       = element_text(size = 11)
    )
}

# Plot 1: Density of TerL positions across all genomes
p1 <- ggplot(terl_positions, aes(x = norm_mid)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, fill = "#2196F3", color = "white",
                 alpha = 0.8, linewidth = 0.2) +
  geom_density(color = "#0D47A1", linewidth = 1.0) +
  geom_rug(alpha = 0.3, linewidth = 0.3, color = "#2196F3") +
  scale_x_continuous(
    name   = "Normalized genome position",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1), breaks = seq(0, 1, 0.1)
  ) +
  scale_y_continuous(name = "Density") +
  labs(
    title    = "TerL (large terminase) position across normalized phage genome length",
    subtitle = paste0("n = ", nrow(terl_positions), " genomes  |  ",
                      "best hit per genome  |  E-value <= ", EVALUE_CUTOFF)
  ) +
  theme_terl()

# Plot 2: By HMM profile — are different TerL families at different positions?
p2 <- ggplot(terl_positions, aes(x = norm_mid, fill = hmm_profile,
                                  color = hmm_profile)) +
  geom_density(alpha = 0.25, linewidth = 0.9) +
  scale_x_continuous(
    name   = "Normalized genome position",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1), breaks = seq(0, 1, 0.1)
  ) +
  scale_y_continuous(name = "Density") +
  scale_fill_brewer(palette = "Set1", name = "HMM profile") +
  scale_color_brewer(palette = "Set1", name = "HMM profile") +
  labs(
    title    = "TerL position by HMM profile",
    subtitle = paste0("n = ", nrow(terl_positions), " genomes")
  ) +
  theme_terl() +
  theme(legend.position = "right")

# Plot 3: Score vs position — are high-confidence hits clustered?
p3 <- ggplot(terl_positions, aes(x = norm_mid, y = score,
                                  color = hmm_profile)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = median(terl_positions$norm_mid),
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  annotate("text",
           x = median(terl_positions$norm_mid) + 0.02,
           y = max(terl_positions$score) * 0.95,
           hjust = 0, size = 3.5, color = "grey40",
           label = paste0("median = ",
                          round(median(terl_positions$norm_mid) * 100, 1), "%")) +
  scale_x_continuous(
    name   = "Normalized genome position",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1), breaks = seq(0, 1, 0.1)
  ) +
  scale_y_continuous(name = "HMM score") +
  scale_color_brewer(palette = "Set1", name = "HMM profile") +
  labs(
    title    = "TerL HMM score vs. genomic position",
    subtitle = "Dashed line = median position"
  ) +
  theme_terl() +
  theme(legend.position = "right")

# ============================================================
# Save plots
# ============================================================
ggsave(file.path(OUTDIR, "terl_position_density.pdf"),
       p1, width = 9, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "terl_position_by_profile.pdf"),
       p2, width = 9, height = 6, device = "pdf")
ggsave(file.path(OUTDIR, "terl_score_vs_position.pdf"),
       p3, width = 9, height = 6, device = "pdf")

message("Written: terl_position_density.pdf")
message("Written: terl_position_by_profile.pdf")
message("Written: terl_score_vs_position.pdf")

cat("\n=== Summary ===\n")
cat("  Genomes with TerL hit:  ", nrow(terl_positions), "\n")
cat("  Genomes without TerL:   ",
    nrow(genome_lengths) - nrow(terl_positions), "\n")
cat("  Median TerL position:   ",
    round(median(terl_positions$norm_mid) * 100, 1), "%\n")
cat("  IQR:                    ",
    round(IQR(terl_positions$norm_mid) * 100, 1), "%\n")
cat("  Output dir:             ", normalizePath(OUTDIR), "\n")
message("Done.")
