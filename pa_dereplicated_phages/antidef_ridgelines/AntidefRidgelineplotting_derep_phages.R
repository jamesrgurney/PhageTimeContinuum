#!/usr/bin/env Rscript
# =============================================================================
# AntidefRidgelineplotting_derep_phages.R
#
# Plots anti-defense system positions along normalized phage genome length
# using reoriented (dnaapler) genomes and per-genome .prt files.
#
# Directory structure:
#   DEFENSEFINDER_DIR/
#     ALLderep_defense_finder_systems.tsv
#     ACCESSION/
#       ACCESSION_reoriented.prt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(scales)
})

# ============================================================
# CONFIG — edit these
# ============================================================
DEFENSEFINDER_DIR <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/defensefinder_results"
SYSTEMS_TSV       <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/ALLderep_defense_finder_systems.tsv"
OUTDIR            <- "/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/antidef_ridgelines"
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
message("  Subtypes: ", paste(sort(unique(systems$subtype)), collapse = ", "))

# Extract clean accession from sys_id: everything before _SUBTYPENAME_NUM
# e.g. AB008550_reoriented_acrif3_1 -> AB008550_reoriented -> AB008550
# accession column should be present; use it directly
# No accession column — extract from protein_in_syst (first gene ID before last _DIGITS)
systems <- systems %>%
  mutate(
    first_gene = str_extract(protein_in_syst, "^[^,]+"),
    replicon   = str_extract(first_gene, "^.+?(?=_\\d+$)")
  )

message("  Unique accessions: ", n_distinct(systems$replicon))

# ============================================================
# Build .prt coordinate lookup from per-genome _reoriented.prt files
# ============================================================
message("\nScanning _reoriented.prt files in: ", DEFENSEFINDER_DIR)

prt_files <- list.files(DEFENSEFINDER_DIR,
                        pattern = "_reoriented\\.prt$",
                        recursive = TRUE,
                        full.names = TRUE)
prt_files <- prt_files[!str_detect(basename(prt_files), "^\\.")]

message("  .prt files found: ", length(prt_files))

prt_all <- map_dfr(prt_files, function(f) {
  tryCatch(parse_prt(f), error = function(e) NULL)
})

message("  Total gene entries parsed: ", nrow(prt_all))

# Genome length = max nt_end per accession
# hit_id format: AB008550_reoriented_1 -> prefix = AB008550_reoriented
prt_all <- prt_all %>%
  mutate(acc_prefix = str_extract(hit_id, "^.+?(?=_\\d+$)"))

genome_lengths <- prt_all %>%
  group_by(acc_prefix) %>%
  summarise(genome_length = max(nt_end, na.rm = TRUE), .groups = "drop") %>%
  mutate(replicon = str_remove(acc_prefix, "_reoriented$"))

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
  group_by(sys_id, replicon, type, subtype) %>%
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
message("Unique phages: ", n_distinct(systems_with_coords$replicon))

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
                 n_distinct(df$replicon), " phage genomes (reoriented)")
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
                 n_distinct(df$replicon), " phage genomes (reoriented)")
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
  "Anti-defense subtype — position across normalized phage genome length",
  "RIDGE_derep_phage_subtype_kernel.png")

make_binline_ridge(sub_df, "subtype_plot",
  "Anti-defense subtype — binned position across normalized phage genome length",
  "RIDGE_derep_phage_subtype_binline.png")

make_kernel_ridge(plot_df %>% filter(!is.na(type_plot)), "type_plot",
  "Anti-defense class — position across normalized phage genome length",
  "RIDGE_derep_phage_class_kernel.png")

make_binline_ridge(plot_df %>% filter(!is.na(type_plot)), "type_plot",
  "Anti-defense class — binned position across normalized phage genome length",
  "RIDGE_derep_phage_class_binline.png")

write_tsv(systems_with_coords,
          file.path(OUTDIR, "derep_phage_antidef_systems_nt_coords.tsv"))
message("  Written: derep_phage_antidef_systems_nt_coords.tsv")

cat("\n=== Summary ===\n")
cat("  Systems in TSV:        ", nrow(systems), "\n")
cat("  Systems with coords:   ", nrow(systems_with_coords), "\n")
cat("  Unique phages:         ", n_distinct(systems_with_coords$replicon), "\n")
cat("  Subtypes:              ", n_distinct(systems_with_coords$subtype), "\n")
cat("  norm_mid IQR:          ", round(spread, 4), "\n")
cat("  Output dir:            ", normalizePath(OUTDIR), "\n")
message("Done.")
