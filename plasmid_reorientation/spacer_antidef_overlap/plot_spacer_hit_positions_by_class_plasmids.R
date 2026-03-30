#!/usr/bin/env Rscript
# =============================================================================
# plot_spacer_hit_positions_by_class.R
#
# Ridgeline plot showing where spacers hit within anti-defense systems,
# by anti-defense class. Mirrors the background binline ridgeline plot
# but uses spacer protospacer midpoint positions instead of system midpoints.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(scales)
})

# ============================================================
# CONFIG
# ============================================================
OVERLAPS_TSV   <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/spacer_antidef_overlap/spacer_antidef_overlaps.tsv"
ANTIDEF_COORDS <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/plasmid_antidef_ridgelines/plasmid_antidef_systems_nt_coords.tsv"
OUTDIR         <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/spacer_antidef_overlap"

BINS         <- 250
RIDGE_SCALE  <- 1.0
MIN_HITS     <- 3    # minimum spacer hits to show a class

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load data
# ============================================================
message("Loading overlaps...")
overlaps <- read_tsv(OVERLAPS_TSV, col_types = cols(.default = "c"),
                     show_col_types = FALSE) %>%
  mutate(
    sstart        = as.integer(sstart),
    send          = as.integer(send),
    genome_length = as.integer(genome_length),
    sys_nt_begin  = as.integer(sys_nt_begin),
    sys_nt_end    = as.integer(sys_nt_end),
    hit_min       = pmin(sstart, send),
    hit_max       = pmax(sstart, send),
    # Normalized midpoint of the spacer hit
    hit_mid_norm  = (hit_min + hit_max) / 2 / genome_length
  ) %>%
  filter(!is.na(hit_mid_norm), hit_mid_norm >= 0, hit_mid_norm <= 1)

message("  Overlap hits loaded: ", nrow(overlaps))

# Load background for overlay
antidef <- read_tsv(ANTIDEF_COORDS, col_types = cols(.default = "c"),
                    show_col_types = FALSE) %>%
  mutate(norm_mid = as.numeric(norm_mid)) %>%
  filter(!is.na(norm_mid))

# ============================================================
# Prepare class-level factors — match order from background plot
# ============================================================
class_order <- antidef %>%
  count(type, sort = TRUE) %>%
  pull(type)

# Filter to classes with enough spacer hits
hit_counts <- overlaps %>%
  count(type, name = "n_hits") %>%
  filter(n_hits >= MIN_HITS)

message("  Classes with >= ", MIN_HITS, " spacer hits: ",
        nrow(hit_counts))
message("  ", paste(hit_counts$type, collapse = ", "))

plot_df <- overlaps %>%
  filter(type %in% hit_counts$type) %>%
  mutate(
    type_plot = factor(type, levels = rev(class_order[class_order %in% type]))
  )

# ============================================================
# Plot 1: Spacer hit positions by anti-defense CLASS (binline)
# ============================================================
nlev <- nlevels(plot_df$type_plot)
h    <- max(7, 0.45 * nlev + 3)

p1 <- ggplot(plot_df, aes(x = hit_mid_norm, y = type_plot,
                           fill = type_plot)) +
  geom_density_ridges(
    stat = "binline",
    aes(height = after_stat(count)),
    bins = BINS,
    scale = RIDGE_SCALE,
    rel_min_height = 0,
    alpha = 0.85,
    color = "white",
    linewidth = 0.2
  ) +
  scale_x_continuous(
    name   = "Normalized genome position (binned)",
    limits = c(0, 1), breaks = seq(0, 1, 0.1),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    y        = NULL,
    title    = "Spacer hit positions within anti-defense systems — by class",
    subtitle = paste0("n = ", nrow(plot_df), " spacer hits  |  ",
                      n_distinct(plot_df$sseqid), " phage genomes  |  ",
                      "min ", MIN_HITS, " hits per class")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    axis.title       = element_text(size = 11)
  )

ggsave(file.path(OUTDIR, "RIDGE_spacer_hits_by_class_binline.png"),
       p1, width = 13, height = h, dpi = 300, bg = "white", device = "png")
message("Written: RIDGE_spacer_hits_by_class_binline.png")

# ============================================================
# Plot 2: Same but kernel density
# ============================================================
p2 <- ggplot(plot_df, aes(x = hit_mid_norm, y = type_plot,
                           fill = type_plot)) +
  geom_density_ridges(
    stat = "density",
    aes(height = after_stat(count)),
    scale = RIDGE_SCALE,
    rel_min_height = 0,
    alpha = 0.85,
    color = "white",
    linewidth = 0.2,
    trim = TRUE
  ) +
  scale_x_continuous(
    name   = "Normalized genome position",
    limits = c(0, 1), breaks = seq(0, 1, 0.1),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    y        = NULL,
    title    = "Spacer hit positions within anti-defense systems — by class",
    subtitle = paste0("n = ", nrow(plot_df), " spacer hits  |  ",
                      n_distinct(plot_df$sseqid), " phage genomes")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    axis.title       = element_text(size = 11)
  )

ggsave(file.path(OUTDIR, "RIDGE_spacer_hits_by_class_kernel.png"),
       p2, width = 13, height = h, dpi = 300, bg = "white", device = "png")
message("Written: RIDGE_spacer_hits_by_class_kernel.png")

# ============================================================
# Plot 3: Subtype level (if enough hits)
# ============================================================
subtype_counts <- overlaps %>%
  count(subtype, type, name = "n_hits") %>%
  filter(n_hits >= MIN_HITS)

if (nrow(subtype_counts) > 0) {
  subtype_order <- antidef %>%
    count(subtype, sort = TRUE) %>%
    pull(subtype)

  plot_sub <- overlaps %>%
    filter(subtype %in% subtype_counts$subtype) %>%
    mutate(
      subtype_plot = factor(subtype,
                            levels = rev(subtype_order[subtype_order %in% subtype]))
    )

  nlev_s <- nlevels(plot_sub$subtype_plot)
  h_s    <- max(7, 0.40 * nlev_s + 3)

  p3 <- ggplot(plot_sub, aes(x = hit_mid_norm, y = subtype_plot,
                              fill = subtype_plot)) +
    geom_density_ridges(
      stat = "binline",
      aes(height = after_stat(count)),
      bins = BINS,
      scale = RIDGE_SCALE,
      rel_min_height = 0,
      alpha = 0.85,
      color = "white",
      linewidth = 0.2
    ) +
    scale_x_continuous(
      name   = "Normalized genome position (binned)",
      limits = c(0, 1), breaks = seq(0, 1, 0.1),
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      y        = NULL,
      title    = "Spacer hit positions within anti-defense systems — by subtype",
      subtitle = paste0("n = ", nrow(plot_sub), " spacer hits  |  ",
                        "min ", MIN_HITS, " hits per subtype")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      axis.title       = element_text(size = 11)
    )

  ggsave(file.path(OUTDIR, "RIDGE_spacer_hits_by_subtype_binline.png"),
         p3, width = 13, height = h_s, dpi = 300, bg = "white", device = "png")
  message("Written: RIDGE_spacer_hits_by_subtype_binline.png")
}

# ============================================================
# Summary table
# ============================================================
summary_tbl <- overlaps %>%
  group_by(type, subtype) %>%
  summarise(
    n_spacer_hits   = n(),
    n_unique_spacers = n_distinct(qseqid),
    n_unique_phages  = n_distinct(sseqid),
    n_systems_hit    = n_distinct(sys_id),
    median_norm_pos  = round(median(hit_mid_norm, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_spacer_hits))

write_tsv(summary_tbl, file.path(OUTDIR, "spacer_hits_summary_by_subtype.tsv"))
message("Written: spacer_hits_summary_by_subtype.tsv")

print(summary_tbl, n = 30)
message("Done.")
