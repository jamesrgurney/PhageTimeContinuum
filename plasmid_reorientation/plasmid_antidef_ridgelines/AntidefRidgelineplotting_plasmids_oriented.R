#!/usr/bin/env Rscript
# =============================================================================
# AntidefRidgelineplotting_plasmids.R
#
# Plots anti-defense system positions along normalized plasmid genome length.
# Uses:
#   - ALLplasmid_defense_finder_systems.tsv  (concatenated systems)
#   - Per-plasmid .prt files for true nucleotide coordinates
#
# Directory structure:
#   DEFENSEFINDER_DIR/
#     ALLplasmid_defense_finder_systems.tsv
#     GCA_XXXX.part_ACCESSION/
#       GCA_XXXX.part_ACCESSION.prt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(scales)
})

# ============================================================
# CONFIG — edit these
# ============================================================
DEFENSEFINDER_DIR <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/defensefinder_results/"
SYSTEMS_TSV       <- file.path(DEFENSEFINDER_DIR, "ALLplasmid_defense_finder_systems.tsv")
OUTDIR            <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/plasmid_antidef_ridgelines"
MIN_GENOME_SIZE   <- 2000
TOP_N_SUBTYPES    <- 30
RIDGE_SCALE       <- 1.2
BINS              <- 250
KERNEL_BW         <- NA_real_

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parse one .prt file -> hit_id, nt_start, nt_end
# ============================================================
parse_prt <- function(prt_path) {
  lines   <- readLines(prt_path, warn = FALSE)
  headers <- lines[startsWith(lines, ">")]
  if (length(headers) == 0) return(tibble())
  tibble(raw = headers) %>%
    mutate(
      hit_id   = str_extract(raw, "(?<=>)\\S+"),
      nt_start = as.integer(str_match(raw, "# (\\d+) # \\d+ # -?1")[, 2]),
      nt_end   = as.integer(str_match(raw, "# \\d+ # (\\d+) # -?1")[, 2])
    ) %>%
    select(hit_id, nt_start, nt_end) %>%
    filter(!is.na(nt_start), !is.na(nt_end))
}

# ============================================================
# Load systems TSV
# ============================================================
message("Loading systems TSV: ", SYSTEMS_TSV)

systems <- read_tsv(SYSTEMS_TSV, col_types = cols(.default = "c"),
                    show_col_types = FALSE) %>%
  rename_with(str_trim)

message("  Systems loaded: ", nrow(systems))
message("  Unique accessions: ", n_distinct(systems$accession))
message("  Subtypes: ", paste(sort(unique(systems$subtype)), collapse = ", "))

# Extract clean accession from full GCA string: everything after .part_
systems <- systems %>%
  mutate(
    replicon = str_extract(accession, "(?<=\\.part_).+$"),
    replicon = if_else(is.na(replicon), accession, replicon)
  )

# ============================================================
# Build .prt coordinate lookup across all accession directories
# ============================================================
message("\nScanning .prt files in: ", DEFENSEFINDER_DIR)

acc_dirs <- list.dirs(DEFENSEFINDER_DIR, full.names = TRUE, recursive = FALSE)
acc_dirs <- acc_dirs[!str_detect(basename(acc_dirs), "\\.log$")]
message("  Subdirectories found: ", length(acc_dirs))

prt_all <- map_dfr(acc_dirs, function(d) {
  dir_name <- basename(d)
  prt_path <- file.path(d, paste0(dir_name, "_reoriented.prt"))
  if (!file.exists(prt_path)) return(NULL)
  tryCatch(parse_prt(prt_path), error = function(e) NULL)
})

message("  Total gene entries parsed from .prt files: ", nrow(prt_all))

# Genome length per accession = max nt_end across all its genes
# hit_id: GCA_XXXX.part_ACCESSION_GENENUM -> strip trailing _DIGITS for prefix
prt_all <- prt_all %>%
  mutate(acc_prefix = str_extract(hit_id, "^.+?(?=_\\d+$)"))

genome_lengths <- prt_all %>%
  group_by(acc_prefix) %>%
  summarise(genome_length = max(nt_end, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    replicon = str_extract(acc_prefix, "(?<=\\.part_).+$"),
    replicon = if_else(is.na(replicon), acc_prefix, replicon)
  )

message("  Genome lengths computed for: ", nrow(genome_lengths), " accessions")

# ============================================================
# Join coordinates to systems via protein_in_syst
# ============================================================
message("\nJoining coordinates to systems...")

systems_with_coords <- systems %>%
  separate_rows(protein_in_syst, sep = ",") %>%
  mutate(protein_in_syst = str_trim(protein_in_syst)) %>%
  left_join(
    prt_all %>% select(hit_id, nt_start, nt_end),
    by = c("protein_in_syst" = "hit_id")
  ) %>%
  group_by(sys_id, accession, replicon, type, subtype) %>%
  summarise(
    sys_nt_begin = min(nt_start, na.rm = TRUE),
    sys_nt_end   = max(nt_end,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sys_nt_begin = if_else(is.infinite(sys_nt_begin), NA_integer_,
                           as.integer(sys_nt_begin)),
    sys_nt_end   = if_else(is.infinite(sys_nt_end),   NA_integer_,
                           as.integer(sys_nt_end))
  ) %>%
  left_join(genome_lengths %>% select(replicon, genome_length),
            by = "replicon") %>%
  filter(
    !is.na(sys_nt_begin), !is.na(sys_nt_end),
    !is.na(genome_length), genome_length >= MIN_GENOME_SIZE
  ) %>%
  mutate(norm_mid = (sys_nt_begin + sys_nt_end) / 2 / genome_length) %>%
  filter(norm_mid >= 0, norm_mid <= 1)

message("Systems with coordinates: ", nrow(systems_with_coords))
message("Unique plasmids: ", n_distinct(systems_with_coords$replicon))

if (nrow(systems_with_coords) == 0)
  stop("No systems with coordinates — check .prt paths and protein_in_syst join")

spread <- IQR(systems_with_coords$norm_mid, na.rm = TRUE)
message("norm_mid IQR: ", round(spread, 4))
if (spread < 0.01)
  warning("Very low IQR — check .prt coordinate parsing")

# ============================================================
# Prepare factors
# ============================================================
top_subtypes <- systems_with_coords %>%
  count(subtype, sort = TRUE) %>%
  slice_head(n = TOP_N_SUBTYPES) %>%
  pull(subtype)

plot_df <- systems_with_coords %>%
  mutate(
    subtype_plot = if_else(subtype %in% top_subtypes, subtype, "other"),
    type_plot    = as.character(type)
  ) %>%
  mutate(
    subtype_plot = factor(subtype_plot,
                          levels = c(top_subtypes[top_subtypes %in% subtype],
                                     "other")),
    type_plot    = factor(type_plot,
                          levels = names(sort(table(type_plot),
                                             decreasing = TRUE)))
  )

# ============================================================
# Plot helpers
# ============================================================
theme_ridge <- function() {
  theme_minimal(base_size = 13) +
    theme(
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      axis.title       = element_text(size = 11)
    )
}

save_png <- function(p, fname, w = 13, nlev = 10) {
  h <- max(7, 0.35 * nlev + 3)
  ggsave(file.path(OUTDIR, fname), plot = p,
         width = w, height = h, dpi = 300, bg = "white", device = "png")
  message("  Written: ", fname)
}

make_kernel_ridge <- function(df, yvar, title, fname) {
  if (!nrow(df)) return(invisible(NULL))
  nlev <- n_distinct(df[[yvar]])
  sub  <- paste0("n = ", nrow(df), " systems  |  ",
                 n_distinct(df$replicon), " plasmid genomes")
  base <- ggplot(df, aes(x = norm_mid, y = .data[[yvar]],
                         fill = .data[[yvar]]))
  p <- if (is.na(KERNEL_BW)) {
    base + geom_density_ridges(stat = "density",
                               aes(height = after_stat(count)),
                               scale = RIDGE_SCALE, rel_min_height = 0,
                               alpha = 0.85, color = "white",
                               linewidth = 0.2, trim = TRUE)
  } else {
    base + geom_density_ridges(stat = "density", bw = KERNEL_BW,
                               aes(height = after_stat(count)),
                               scale = RIDGE_SCALE, rel_min_height = 0,
                               alpha = 0.85, color = "white",
                               linewidth = 0.2, trim = TRUE)
  }
  p + scale_x_continuous(name = "Normalized genome position",
                          limits = c(0, 1), breaks = seq(0, 1, 0.1),
                          labels = percent_format(accuracy = 1)) +
    labs(y = NULL, title = title, subtitle = sub) +
    theme_ridge() -> p
  save_png(p, fname, nlev = nlev)
}

make_binline_ridge <- function(df, yvar, title, fname) {
  if (!nrow(df)) return(invisible(NULL))
  nlev <- n_distinct(df[[yvar]])
  sub  <- paste0("n = ", nrow(df), " systems  |  ",
                 n_distinct(df$replicon), " plasmid genomes")
  p <- ggplot(df, aes(x = norm_mid, y = .data[[yvar]],
                      fill = .data[[yvar]])) +
    geom_density_ridges(stat = "binline",
                        aes(height = after_stat(count)),
                        bins = BINS, scale = RIDGE_SCALE,
                        rel_min_height = 0, alpha = 0.85,
                        color = "white", linewidth = 0.2) +
    scale_x_continuous(name = "Normalized genome position (binned)",
                        limits = c(0, 1), breaks = seq(0, 1, 0.1),
                        labels = percent_format(accuracy = 1)) +
    labs(y = NULL, title = title, subtitle = sub) +
    theme_ridge()
  save_png(p, fname, nlev = nlev)
}

# ============================================================
# Generate plots
# ============================================================
message("\nGenerating plots...")

sub_df <- plot_df %>%
  filter(!is.na(subtype_plot)) %>%
  add_count(subtype_plot, name = "n_sub") %>%
  filter(n_sub >= 3) %>%
  droplevels()

make_kernel_ridge(sub_df, "subtype_plot",
  "Anti-defense subtype — position across normalized plasmid genome length",
  "RIDGE_plasmid_subtype_kernel.png")

make_binline_ridge(sub_df, "subtype_plot",
  "Anti-defense subtype — binned position across normalized plasmid genome length",
  "RIDGE_plasmid_subtype_binline.png")

make_kernel_ridge(plot_df %>% filter(!is.na(type_plot)), "type_plot",
  "Anti-defense class — position across normalized plasmid genome length",
  "RIDGE_plasmid_class_kernel.png")

make_binline_ridge(plot_df %>% filter(!is.na(type_plot)), "type_plot",
  "Anti-defense class — binned position across normalized plasmid genome length",
  "RIDGE_plasmid_class_binline.png")

write_tsv(systems_with_coords,
          file.path(OUTDIR, "plasmid_antidef_systems_nt_coords.tsv"))
message("  Written: plasmid_antidef_systems_nt_coords.tsv")

cat("\n=== Summary ===\n")
cat("  Systems in TSV:        ", nrow(systems), "\n")
cat("  Systems with coords:   ", nrow(systems_with_coords), "\n")
cat("  Unique plasmids:       ", n_distinct(systems_with_coords$replicon), "\n")
cat("  Subtypes:              ", n_distinct(systems_with_coords$subtype), "\n")
cat("  norm_mid IQR:          ", round(spread, 4), "\n")
cat("  Output dir:            ", normalizePath(OUTDIR), "\n")
message("Done.")
