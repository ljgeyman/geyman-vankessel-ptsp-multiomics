#!/usr/bin/env Rscript

# ==============================================================================
# SmcR TF ChIP-seq Quantitative Analysis: Dynamic FE vs. Paired CPM NCIS Pipeline
# ------------------------------------------------------------------------------
# Supports USE_FE_BWS = TRUE  (04b_bigwigs: pre-computed *_FE.bw)
#      and USE_FE_BWS = FALSE (04_bigwigs:  pairs *_IP.bw with *_Control.bw to
#                              compute true un-floored NCIS-calibrated log2(IP/Input))
#
# Step 1: Regulon & Direct Target Peak Identification (4 Tiers & 2 Modes)
# Step 2: Deep Background (NCIS/SES) Calibration Across Peak-Depleted Genome
# Step 3: Analysis 1A & 1B - Summit-Centered Calibrated Log2(IP/Input) Metagenes
# Step 4: Analysis 2A & 2B - Peak Occupancy [-150, +150 bp] with Normality-
#         Adaptive Omnibus (RM-ANOVA / Friedman), Post-Hoc, & CLD Lettering
# Step 5: Analysis 3A & 3B - Difference Metagenes (Delta Log2(IP/Input) + 95% CI)
# Step 6: DiffBind Differential Binding Analysis (Appended to Excel Workbook)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(DESeq2)
  library(openxlsx)
  library(ashr)
  library(zoo)
  library(Biostrings)
  library(lme4)
  library(lmerTest)
  library(broom)
  library(broom.mixed)
  library(DiffBind)
})

# ==============================================================================
# 1. PATHS & PARAMETERS
# ==============================================================================
# ---> GENOME & RNA-SEQ PATHS
ORIGINAL_REFSEQ_GTF <- "refs/refseq_vvul.gtf"
RNA_COUNTS          <- "analysis_results/counts/all_samples_counts.txt"
METADATA_FILE       <- "metadata.txt"

# ---> TF CHIP PEAK & MOTIF PATHS
PEAKS_FILE          <- "05_final_output/combined_macs2_peaks.txt"
SUMMITS_EXCEL       <- "05_final_output/macs2_consolidated_peaks.xlsx"
FIMO_FILE           <- "PTSP_meme_chip_analysis/meme_chip_results/fimo_out_1/fimo.tsv"

# ---> DIFFBIND PATHS & PARAMETERS
# Provide a DiffBind CSV sample sheet OR leave NULL to auto-construct from BAM & Peak directories
DIFFBIND_SAMPLESHEET <- NULL
CHIP_BAM_DIR         <- "02_alignments"
CHIP_PEAK_DIR        <- "03_macs2_peaks"

# ---> LOCI EXCLUSIONS (e.g. self-deleted genes or highly variable pseudo-genes)
EXCLUDED_LOCI       <- c("AOT11_RS09875", "AOT11_RS01725")

# ---> TF CHIP BIGWIG MODE SWITCH
# Set USE_FE_BWS = FALSE to pair un-floored *_IP.bw and *_Control.bw from 04_bigwigs
# (recommended to rescue global binding shifts from MACS2 FE flooring),
# or TRUE to use pre-computed *_FE.bw files from 04b_bigwigs.
USE_FE_BWS          <- FALSE
TF_BW_DIR           <- if (USE_FE_BWS) "04b_bigwigs" else "04_bigwigs"
PSEUDOCOUNT         <- 1e-3

# ---> DYNAMIC REGEX PATTERNS FOR SMCR TF CHIP BIGWIGS
if (USE_FE_BWS) {
  TF_DMSO_IP_PATTERN  <- "^(wt|dluxO).*DMSO.*_FE[.]bw$"
  TF_PTSP_IP_PATTERN  <- "^(wt|dluxO).*PTSP.*_FE[.]bw$"
  TF_DSMCR_IP_PATTERN <- "^(d?smcR).*DMSO.*_FE[.]bw$"
} else {
  TF_DMSO_IP_PATTERN  <- "^(wt).*DMSO.*IP[.]bw$" 
  TF_PTSP_IP_PATTERN  <- "^(wt).*PTSP.*IP[.]bw$" 
  TF_DSMCR_IP_PATTERN <- "^(d?smcR).*DMSO.*_IP[.]bw$"
}

# ---> DEEP BACKGROUND (NCIS / SES) CALIBRATION PARAMETERS
BG_TILE_WIDTH       <- 50   # Tile width (bp) across peak-depleted genome
Peak_Exclusion_Pad  <- 50   # Buffer (bp) around MACS2 peaks excluded from background

# ---> REGULON & WINDOW PARAMETERS
PADJ_CUTOFF         <- 0.05
LFC_CUTOFF          <- 1.0
PEAK_MAXGAP         <- 2000

# Summit-centered Metagene Window [-400, +400 bp] & Smoothing
UP_DIST             <- 400
DOWN_DIST           <- 400
SMOOTH_K            <- 51

# Peak Quantification Window around Summit/Motif Center: [-150, +150 bp]
PEAK_QUANT_HALFWIN  <- 400

# Bootstrap iterations for Analysis 3 Difference Metagenes
N_BOOT              <- 1000

OUT_DIR             <- if (USE_FE_BWS) "08_tf_chip_fe_calibrated" else "08_tf_chip_paired_cpm_ncis"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

assert_bws_exist <- function(paths, label) {
  if (length(paths) == 0) {
    stop(sprintf("[FATAL] No BigWig files found for '%s' in '%s'. Check directory path and filename regex.", label, TF_BW_DIR), call. = FALSE)
  }
  message(sprintf("[INFO] Found %d BigWig file(s) for %s:\n       %s", length(paths), label, paste(basename(paths), collapse = ", ")))
}

# ==============================================================================
# 2. ANNOTATIONS, DESEQ2 & REGULON CLASSIFICATION (2 MODES & 4 TIERS)
# ==============================================================================
message("[INFO] Parsing GTF annotations and running DESeq2 for regulon classification...")

ref_gtf   <- import(ORIGINAL_REFSEQ_GTF)
ref_genes <- ref_gtf[ref_gtf$type == "gene"]

raw_ids <- if ("locus_tag" %in% colnames(mcols(ref_genes))) ref_genes$locus_tag else ref_genes$gene_id
if (is(raw_ids, "List") || is(raw_ids, "CharacterList")) raw_ids <- unlist(as.list(raw_ids))
ref_genes$locus_tag <- str_remove(as.character(raw_ids), "^gene-")

raw_counts_full <- suppressWarnings(read_tsv(RNA_COUNTS, comment = "#", show_col_types = FALSE))
colnames(raw_counts_full) <- colnames(raw_counts_full) %>%
  str_remove("^.*bam[/\\\\.]") %>%
  str_remove("_sorted.*$")

count_matrix <- raw_counts_full %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

if (file.exists(METADATA_FILE)) {
  metadata <- read.table(METADATA_FILE, header = FALSE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
  metadata <- metadata[colnames(count_matrix), , drop = FALSE]
  metadata$genotype  <- tolower(sapply(strsplit(rownames(metadata), "_"), `[`, 1))
  metadata$condition <- tolower(sapply(strsplit(rownames(metadata), "_"), `[`, 2))
} else {
  sample_names <- colnames(count_matrix)
  metadata <- data.frame(
    genotype  = tolower(str_extract(sample_names, "^[a-zA-Z0-9]+")),
    condition = tolower(str_extract(sample_names, "(dmso|ptsp)")),
    row.names = sample_names
  )
}
metadata$group <- factor(paste(metadata$genotype, metadata$condition, sep = "_"))

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

get_shrunken_res <- function(dds_obj, num_group, den_group) {
  res <- lfcShrink(dds_obj, contrast = c("group", num_group, den_group), type = "ashr", quiet = TRUE)
  as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    filter(!gene_id %in% EXCLUDED_LOCI)
}

res_dSmcR <- get_shrunken_res(dds, "dsmcr_dmso", "dluxo_dmso")
res_inhib <- get_shrunken_res(dds, "dluxo_ptsp", "dluxo_dmso")

rna_res_all <- merge(res_dSmcR, res_inhib, by = "gene_id", suffixes = c("_dSmcR", "_Inhib")) %>%
  mutate(
    Mode = case_when(
      log2FoldChange_dSmcR >  LFC_CUTOFF & padj_dSmcR < PADJ_CUTOFF ~ "Repression",
      log2FoldChange_dSmcR < -LFC_CUTOFF & padj_dSmcR < PADJ_CUTOFF ~ "Activation",
      TRUE ~ "Independent"
    ),
    Mode_Label = factor(
      Mode,
      levels = c("Repression", "Activation", "Independent"),
      labels = c("Directly Repressed by SmcR", "Directly Activated by SmcR", "Independent")
    ),
    Tier = case_when(
      Mode == "Repression" & log2FoldChange_Inhib >  LFC_CUTOFF & padj_Inhib < PADJ_CUTOFF ~ "Repression_Sensitive",
      Mode == "Repression" ~ "Repression_Muted",
      Mode == "Activation" & log2FoldChange_Inhib < -LFC_CUTOFF & padj_Inhib < PADJ_CUTOFF ~ "Activation_Sensitive",
      Mode == "Activation" ~ "Activation_Muted",
      TRUE ~ "Not_Direct_Target"
    ),
    Tier_Label = factor(
      Tier,
      levels = c("Repression_Sensitive", "Repression_Muted", "Activation_Sensitive", "Activation_Muted", "Not_Direct_Target"),
      labels = c("Repressed: Sensitive", "Repressed: Muted", "Activated: Sensitive", "Activated: Muted", "Not Direct Target")
    ),
    Response_Magnitude = abs(log2FoldChange_Inhib),
    Baseline_Magnitude = abs(log2FoldChange_dSmcR)
  )

# ==============================================================================
# 3. MAP PEAKS TO REGULATED LOCI & ANCHOR AT EXACT BINDING SUMMIT / MOTIF CENTER
# ==============================================================================
message("[INFO] Mapping MACS2 peaks and FIMO motifs to define summit-anchored TF binding sites...")

peaks <- read.delim(PEAKS_FILE, stringsAsFactors = FALSE)
if (!"peak_id" %in% colnames(peaks)) colnames(peaks)[4] <- "peak_id"

peaks_gr <- makeGRangesFromDataFrame(
  peaks,
  keep.extra.columns = TRUE,
  seqnames.field     = "chrom",
  start.field        = "chromStart",
  end.field          = "chromEnd"
)

fimo <- read_tsv(FIMO_FILE, comment = "#", show_col_types = FALSE) %>%
  filter(!is.na(matched_sequence))

fimo_gr <- makeGRangesFromDataFrame(
  fimo,
  keep.extra.columns = TRUE,
  seqnames.field     = "sequence_name",
  start.field        = "start",
  end.field          = "stop"
)

gene_overlaps   <- findOverlaps(peaks_gr, ref_genes, maxgap = PEAK_MAXGAP)
annotated_peaks <- peaks[queryHits(gene_overlaps), ]
annotated_peaks$gene_id <- ref_genes$locus_tag[subjectHits(gene_overlaps)]

motif_overlaps <- findOverlaps(peaks_gr, fimo_gr)
overlap_df <- data.frame(
  peak_id          = peaks$peak_id[queryHits(motif_overlaps)],
  motif_start      = start(fimo_gr)[subjectHits(motif_overlaps)],
  motif_end        = end(fimo_gr)[subjectHits(motif_overlaps)],
  motif_strand     = as.character(strand(fimo_gr))[subjectHits(motif_overlaps)],
  motif_pvalue     = fimo_gr$`p-value`[subjectHits(motif_overlaps)],
  matched_sequence = fimo_gr$matched_sequence[subjectHits(motif_overlaps)]
) %>%
  mutate(motif_center = floor((motif_start + motif_end) / 2))

best_motifs <- overlap_df %>%
  group_by(peak_id) %>%
  slice_min(order_by = motif_pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

peaks_motifs <- inner_join(annotated_peaks, best_motifs, by = "peak_id")

master_df <- rna_res_all %>%
  filter(Mode != "Independent") %>%
  inner_join(peaks_motifs, by = "gene_id") %>%
  group_by(gene_id) %>%
  filter(!gene_id %in% EXCLUDED_LOCI) %>%
  slice_min(order_by = motif_pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Tier_Label  = droplevels(Tier_Label),
    Mode_Label  = droplevels(Mode_Label),
    peak_center = ifelse(!is.na(motif_center), motif_center, floor((chromStart + chromEnd) / 2))
  )

master_peaks_unique <- master_df %>%
  group_by(peak_id) %>%
  slice_max(order_by = Baseline_Magnitude, n = 1, with_ties = FALSE) %>%
  ungroup()

tf_summits_gr <- GRanges(
  seqnames = master_peaks_unique$chrom,
  ranges   = IRanges(start = master_peaks_unique$peak_center, end = master_peaks_unique$peak_center),
  strand   = ifelse(master_peaks_unique$motif_strand %in% c("+", "-"), master_peaks_unique$motif_strand, "+")
)
mcols(tf_summits_gr) <- master_peaks_unique %>%
  dplyr::select(peak_id, gene_id, Mode, Mode_Label, Tier, Tier_Label,
                log2FoldChange_dSmcR, log2FoldChange_Inhib, motif_pvalue, matched_sequence)

message(sprintf("[INFO] Mapped %d unique directly regulated SmcR peaks (linked to %d regulated genes).",
                length(tf_summits_gr), nrow(master_df)))

# ==============================================================================
# 4. DYNAMIC BIGWIG PAIRING & DEEP BACKGROUND (NCIS / SES) CALIBRATION
# ==============================================================================
message(sprintf("[INFO] Mode: USE_FE_BWS = %s | Scanning '%s'...", USE_FE_BWS, TF_BW_DIR))

tf_dmso_ip_bws  <- list.files(TF_BW_DIR, pattern = TF_DMSO_IP_PATTERN,  full.names = TRUE, ignore.case = TRUE)
tf_ptsp_ip_bws  <- list.files(TF_BW_DIR, pattern = TF_PTSP_IP_PATTERN,  full.names = TRUE, ignore.case = TRUE)
tf_dsmcr_ip_bws <- list.files(TF_BW_DIR, pattern = TF_DSMCR_IP_PATTERN, full.names = TRUE, ignore.case = TRUE)

# Exclude any RNAP tracks if present in the same directory
tf_dmso_ip_bws  <- tf_dmso_ip_bws[!grepl("(clpP|r179c)", basename(tf_dmso_ip_bws), ignore.case = TRUE)]
tf_ptsp_ip_bws  <- tf_ptsp_ip_bws[!grepl("(clpP|r179c)", basename(tf_ptsp_ip_bws), ignore.case = TRUE)]
tf_dsmcr_ip_bws <- tf_dsmcr_ip_bws[!grepl("(clpP|r179c)", basename(tf_dsmcr_ip_bws), ignore.case = TRUE)]

assert_bws_exist(tf_dmso_ip_bws,  "SmcR TF ChIP wt/dluxO DMSO")
assert_bws_exist(tf_ptsp_ip_bws,  "SmcR TF ChIP wt/dluxO PTSP")
assert_bws_exist(tf_dsmcr_ip_bws, "SmcR TF ChIP dsmcR DMSO")

all_ip_bws <- c(tf_dmso_ip_bws, tf_ptsp_ip_bws, tf_dsmcr_ip_bws)

# Helper to locate paired *_Control.bw when USE_FE_BWS == FALSE
resolve_control_bw <- function(ip_bw_path) {
  if (USE_FE_BWS) return(NA_character_)
  ctrl_path <- sub("_IP[.]bw$", "_Control.bw", ip_bw_path, ignore.case = TRUE)
  if (!file.exists(ctrl_path)) {
    stop(sprintf("[FATAL] Paired Input/Control BigWig not found for '%s'. Expected: '%s'",
                 basename(ip_bw_path), basename(ctrl_path)), call. = FALSE)
  }
  return(ctrl_path)
}

# Build peak-depleted background tiles across the genome:
# Exclude MACS2 called peaks (+ Peak_Exclusion_Pad) and only FIMO hits that fall inside called peaks
first_bw_info <- seqinfo(import(all_ip_bws[1]))
genome_gr     <- GRanges(seqnames = seqnames(first_bw_info), ranges = IRanges(start = 1, width = seqlengths(first_bw_info)))

bound_fimo_gr <- subsetByOverlaps(fimo_gr, peaks_gr)
all_exclusion_zones <- IRanges::reduce(c(
  resize(peaks_gr,      width = width(peaks_gr)      + 2 * Peak_Exclusion_Pad, fix = "center"),
  resize(bound_fimo_gr, width = width(bound_fimo_gr) + 2 * Peak_Exclusion_Pad, fix = "center")
))

bg_regions_gr <- GenomicRanges::setdiff(genome_gr, all_exclusion_zones)
bg_tiles_gr   <- unlist(tile(bg_regions_gr, width = BG_TILE_WIDTH))
bg_tiles_gr   <- bg_tiles_gr[width(bg_tiles_gr) == BG_TILE_WIDTH]

message(sprintf("[INFO] Constructed %d peak-depleted background tiles (%d bp each; %.2f Mb total background)...",
                length(bg_tiles_gr), BG_TILE_WIDTH, sum(width(bg_tiles_gr)) / 1e6))

# Extract width-weighted linear mean across regions for a single BigWig
get_tile_linear_means <- function(bw_path, gr) {
  s <- import(bw_path, which = gr)
  hits <- findOverlaps(gr, s)
  overlaps <- pintersect(gr[queryHits(hits)], s[subjectHits(hits)])
  
  df_ov <- data.frame(
    idx   = queryHits(hits),
    score = pmax(s$score[subjectHits(hits)], PSEUDOCOUNT),
    w     = width(overlaps)
  ) %>%
    group_by(idx) %>%
    summarise(mean_val = weighted.mean(score, w = w, na.rm = TRUE), .groups = "drop") %>%
    complete(idx = seq_along(gr), fill = list(mean_val = PSEUDOCOUNT))
  
  return(df_ov$mean_val)
}

# Compute raw uncalibrated log2(IP/Input) across GRanges for a given replicate
# - If USE_FE_BWS = TRUE:  log2(FE_score)
# - If USE_FE_BWS = FALSE: log2((IP_CPM + PSEUDOCOUNT) / (Control_CPM + PSEUDOCOUNT))
get_raw_log2_ip_over_input_regions <- function(ip_bw_path, ctrl_bw_path, gr) {
  ip_vals <- get_tile_linear_means(ip_bw_path, gr)
  if (USE_FE_BWS) {
    return(log2(pmax(ip_vals, PSEUDOCOUNT)))
  } else {
    ctrl_vals <- get_tile_linear_means(ctrl_bw_path, gr)
    return(log2((ip_vals + PSEUDOCOUNT) / (ctrl_vals + PSEUDOCOUNT)))
  }
}

# Compute NCIS / SES Deep Background Calibration Parameters per replicate
bg_calib_list <- lapply(all_ip_bws, function(ip_bw) {
  ctrl_bw <- resolve_control_bw(ip_bw)
  cond_group <- case_when(
    ip_bw %in% tf_dmso_ip_bws  ~ "wt_DMSO",
    ip_bw %in% tf_ptsp_ip_bws  ~ "wt_PTSP",
    ip_bw %in% tf_dsmcr_ip_bws ~ "dsmcR_DMSO",
    TRUE                       ~ "Unknown"
  )
  
  ip_bg_tiles <- get_tile_linear_means(ip_bw, bg_tiles_gr)
  
  if (USE_FE_BWS) {
    ctrl_bg_tiles <- rep(1.0, length(ip_bg_tiles))
    ratio_tiles   <- pmax(ip_bg_tiles, PSEUDOCOUNT)
  } else {
    ctrl_bg_tiles <- get_tile_linear_means(ctrl_bw, bg_tiles_gr)
    ratio_tiles   <- (ip_bg_tiles + PSEUDOCOUNT) / (ctrl_bg_tiles + PSEUDOCOUNT)
  }
  
  # Filter to bottom 90% of background tiles to exclude any unannotated weak binding sites
  q90      <- quantile(ratio_tiles, probs = 0.90, na.rm = TRUE)
  clean_ok <- which(ratio_tiles <= q90 & ip_bg_tiles > PSEUDOCOUNT)
  
  ip_bg_med    <- median(ip_bg_tiles[clean_ok], na.rm = TRUE)
  ctrl_bg_med  <- if (USE_FE_BWS) NA_real_ else median(ctrl_bg_tiles[clean_ok], na.rm = TRUE)
  ratio_bg_med <- median(ratio_tiles[clean_ok], na.rm = TRUE)
  log2_bg_med  <- median(log2(ratio_tiles[clean_ok]), na.rm = TRUE)
  log2_bg_mad  <- mad(log2(ratio_tiles[clean_ok]), constant = 1.4826, na.rm = TRUE)
  
  data.frame(
    ip_bw_path           = ip_bw,
    ctrl_bw_path         = ctrl_bw,
    sample_id            = basename(ip_bw),
    paired_control       = if (USE_FE_BWS) "Pre-computed FE" else basename(ctrl_bw),
    Condition            = cond_group,
    IP_BG_Median_CPM     = ip_bg_med,
    Input_BG_Median_CPM  = ctrl_bg_med,
    BG_IP_Input_Linear   = ratio_bg_med,
    BG_Log2_IP_Input_Med = log2_bg_med,
    BG_Log2_MAD          = log2_bg_mad,
    NCIS_Linear_Scale    = 1.0 / ratio_bg_med,
    Log2_Calib_Offset    = -log2_bg_med,
    stringsAsFactors     = FALSE
  )
}) %>% bind_rows()

bg_calib_params <- bg_calib_list

message("\n[INFO] Deep Background (NCIS/SES) Log2(IP/Input) Calibration Parameters:")
print(bg_calib_params %>%
        dplyr::select(sample_id, paired_control, Condition, IP_BG_Median_CPM, Input_BG_Median_CPM,
                      BG_IP_Input_Linear, BG_Log2_IP_Input_Med, NCIS_Linear_Scale, Log2_Calib_Offset))
write.xlsx(bg_calib_params, file.path(OUT_DIR, "00_TF_Deep_Background_Calibration_Params.xlsx"), overwrite = TRUE)

get_sample_calib_row <- function(ip_bw_path) {
  row <- bg_calib_params[bg_calib_params$ip_bw_path == ip_bw_path, ]
  if (nrow(row) == 0) stop(paste("BigWig not found in calibration table:", ip_bw_path))
  return(row[1, ])
}

# Diagnostic Plot: Background Log2(IP/Input) Tile Distributions Pre- vs Post-Calibration
bg_diag_df <- lapply(seq_len(nrow(bg_calib_params)), function(i) {
  ip_bw   <- bg_calib_params$ip_bw_path[i]
  ctrl_bw <- bg_calib_params$ctrl_bw_path[i]
  raw_lfc <- get_raw_log2_ip_over_input_regions(ip_bw, ctrl_bw, bg_tiles_gr)
  off     <- bg_calib_params$Log2_Calib_Offset[i]
  
  bind_rows(
    data.frame(Sample = bg_calib_params$sample_id[i], Condition = bg_calib_params$Condition[i],
               Stage = "1. Uncalibrated Background Log2(IP / Input)", Value = raw_lfc),
    data.frame(Sample = bg_calib_params$sample_id[i], Condition = bg_calib_params$Condition[i],
               Stage = "2. NCIS Calibrated Background Log2(IP / Input)", Value = raw_lfc + off)
  )
}) %>% bind_rows()

p_bg_diag <- ggplot(bg_diag_df, aes(x = Value, color = Sample, linetype = Condition)) +
  geom_density(linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~ Stage, scales = "free_y") +
  coord_cartesian(xlim = quantile(bg_diag_df$Value, c(0.01, 0.99), na.rm = TRUE)) +
  theme_bw(base_size = 12) +
  labs(
    title    = "Deep Background (NCIS/SES) Calibration of Log2(IP / Input) Across Non-Peak Genome Tiles",
    subtitle = sprintf("Mode: %s | %d non-peak 500-bp tiles (>%d bp from any MACS2 peak)",
                       ifelse(USE_FE_BWS, "Pre-computed FE BigWigs", "Paired IP & Control CPM BigWigs"),
                       length(bg_tiles_gr), Peak_Exclusion_Pad),
    x        = expression("Background Tile "*Log[2]*"(IP / Input)"),
    y        = "Density"
  ) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(OUT_DIR, "00_TF_Deep_Background_Density_Diagnostic.tiff"), p_bg_diag, width = 10, height = 5, dpi = 300)

# ==============================================================================
# 5. ANALYSIS 1A & 1B: SUMMIT-CENTERED CALIBRATED LOG2(IP/INPUT) METAGENES
#    (1A: 4 FUNCTIONAL TIERS | 1B: 2 REGULATION MODES)
# ==============================================================================
message("[INFO] Extracting summit-centered NCIS-calibrated Log2(IP/Input) metagene matrices [-400, +400 bp]...")

PAD_BP      <- floor(SMOOTH_K / 2)
EXT_UP      <- UP_DIST + PAD_BP
EXT_DOWN    <- DOWN_DIST + PAD_BP
x_positions <- seq(-EXT_UP, EXT_DOWN)

# Extract base-pair-resolved calibrated log2(IP/Input) matrix around summits
get_tf_calibrated_log2fe_matrix <- function(ip_bw_paths, gr, up_dist, dn_dist) {
  total_bins <- up_dist + dn_dist + 1
  n_loci     <- length(gr)
  if (length(ip_bw_paths) == 0 || n_loci == 0) {
    return(matrix(0, nrow = n_loci, ncol = total_bins))
  }
  
  summit_windows <- promoters(gr, upstream = up_dist, downstream = dn_dist + 1)
  accum_mat      <- matrix(0, nrow = n_loci, ncol = total_bins)
  
  for (ip_bw in ip_bw_paths) {
    cal_row  <- get_sample_calib_row(ip_bw)
    log2_off <- cal_row$Log2_Calib_Offset
    
    ip_data  <- import(ip_bw)
    ip_cov   <- coverage(ip_data, weight = "score")
    
    if (!USE_FE_BWS) {
      ctrl_data <- import(cal_row$ctrl_bw_path)
      ctrl_cov  <- coverage(ctrl_data, weight = "score")
    }
    
    mat_list <- lapply(seq_along(summit_windows), function(i) {
      chrom      <- as.character(seqnames(summit_windows)[i])
      start_pos  <- start(summit_windows)[i]
      end_pos    <- end(summit_windows)[i]
      strand_dir <- as.character(strand(summit_windows)[i])
      
      if (!chrom %in% names(ip_cov)) return(rep(0, total_bins))
      
      chr_ip  <- ip_cov[[chrom]]
      chr_len <- length(chr_ip)
      
      s <- max(1L, start_pos)
      e <- min(chr_len, end_pos)
      
      if (s <= e) {
        ip_vec <- as.numeric(chr_ip[s:e])
        if (USE_FE_BWS) {
          vec <- log2(pmax(ip_vec, PSEUDOCOUNT)) + log2_off
        } else {
          chr_ctrl <- ctrl_cov[[chrom]]
          ctrl_vec <- as.numeric(chr_ctrl[s:e])
          vec      <- log2((pmax(ip_vec, 0) + PSEUDOCOUNT) / (pmax(ctrl_vec, 0) + PSEUDOCOUNT)) + log2_off
        }
        if (start_pos < 1L) vec <- c(rep(0, 1L - start_pos), vec)
        if (end_pos > chr_len) vec <- c(vec, rep(0, end_pos - chr_len))
      } else {
        vec <- rep(0, total_bins)
      }
      
      if (strand_dir == "-") vec <- rev(vec)
      return(vec)
    })
    
    bw_mat <- do.call(rbind, mat_list)
    if (is.null(dim(bw_mat))) bw_mat <- matrix(bw_mat, nrow = 1)
    accum_mat <- accum_mat + bw_mat
  }
  
  return(accum_mat / length(ip_bw_paths))
}

tf_mat_dmso  <- get_tf_calibrated_log2fe_matrix(tf_dmso_ip_bws,  tf_summits_gr, EXT_UP, EXT_DOWN)
tf_mat_ptsp  <- get_tf_calibrated_log2fe_matrix(tf_ptsp_ip_bws,  tf_summits_gr, EXT_UP, EXT_DOWN)
tf_mat_dsmcr <- get_tf_calibrated_log2fe_matrix(tf_dsmcr_ip_bws, tf_summits_gr, EXT_UP, EXT_DOWN)

build_and_plot_tf_metagene <- function(group_vec, plot_title, out_filename, ncol_facet = 2, fig_w = 9.5, fig_h = 5.0) {
  grp_levels  <- levels(droplevels(factor(group_vec)))
  df_metagene <- data.frame()
  
  for (grp in grp_levels) {
    idx <- which(group_vec == grp)
    if (length(idx) == 0) next
    
    dmso_prof  <- apply(tf_mat_dmso[idx, , drop = FALSE],  2, median, na.rm = TRUE)
    ptsp_prof  <- apply(tf_mat_ptsp[idx, , drop = FALSE],  2, median, na.rm = TRUE)
    dsmcr_prof <- apply(tf_mat_dsmcr[idx, , drop = FALSE], 2, median, na.rm = TRUE)
    
    df_metagene <- bind_rows(
      df_metagene,
      data.frame(Position = x_positions, Signal = dmso_prof,  Plot_Group = "wt_DMSO",    Group_Label = grp, N_Peaks = length(idx)),
      data.frame(Position = x_positions, Signal = ptsp_prof,  Plot_Group = "wt_PTSP",    Group_Label = grp, N_Peaks = length(idx)),
      data.frame(Position = x_positions, Signal = dsmcr_prof, Plot_Group = "dsmcR_DMSO", Group_Label = grp, N_Peaks = length(idx))
    )
  }
  
  df_metagene <- df_metagene %>%
    group_by(Plot_Group, Group_Label) %>%
    arrange(Position) %>%
    mutate(Signal = rollmean(Signal, k = SMOOTH_K, fill = NA, align = "center")) %>%
    ungroup() %>%
    filter(!is.na(Signal) & Position >= -UP_DIST & Position <= DOWN_DIST) %>%
    mutate(
      Plot_Group  = factor(Plot_Group, levels = c("wt_DMSO", "wt_PTSP", "dsmcR_DMSO")),
      Group_Label = factor(Group_Label, levels = grp_levels),
      Facet_Label = factor(
        paste0(Group_Label, " (n=", N_Peaks, " peaks)"),
        levels = unique(paste0(Group_Label, " (n=", N_Peaks, " peaks)"))
      )
    )
  
  y_max        <- max(df_metagene$Signal, na.rm = TRUE)
  annotation_y <- y_max * 1.05
  gene_model_y <- min(df_metagene$Signal, na.rm = TRUE) * 0.9
  
  p_meta <- ggplot(df_metagene, aes(x = Position, y = Signal, group = Plot_Group, color = Plot_Group, linetype = Plot_Group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30", linewidth = 0.7, alpha = 0.6) +
    geom_ribbon(aes(ymin = -Inf, ymax = Signal, fill = Plot_Group), alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line(data = subset(df_metagene, Plot_Group == "dsmcR_DMSO"), linewidth = 0.8, alpha = 0.85) +
    geom_line(data = subset(df_metagene, Plot_Group != "dsmcR_DMSO"), linewidth = 0.95) +
    coord_cartesian(ylim = c(-6, 0.5)) +
    facet_wrap(~ Facet_Label, ncol = ncol_facet) +
    scale_color_manual(
      values = c("wt_DMSO" = "#4E79A7", "wt_PTSP" = "#B00006", "dsmcR_DMSO" = "grey35"),
      labels = c("wt_DMSO" = "wt + DMSO", "wt_PTSP" = "wt + PTSP", "dsmcR_DMSO" = expression(Delta*italic("smcR")*" + DMSO"))
    ) +
    scale_fill_manual(
      values = c("wt_DMSO" = "#4E79A7", "wt_PTSP" = "#B00006", "dsmcR_DMSO" = "grey35")
    ) +
    scale_linetype_manual(
      values = c("wt_DMSO" = "solid", "wt_PTSP" = "solid", "dsmcR_DMSO" = "dashed"),
      labels = c("wt_DMSO" = "wt + DMSO", "wt_PTSP" = "wt + PTSP", "dsmcR_DMSO" = expression(Delta*italic("smcR")*" + DMSO"))
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = 12),
      legend.position  = "bottom",
      axis.line        = element_line(colour = "black", linewidth = 0.5),
      axis.text        = element_text(face = "bold")
    ) +
    scale_x_continuous(
      breaks = c(-200, 0, 200),
      labels = c("-200 bp", "Binding Summit", "+200 bp"),
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(min(0, gene_model_y * 0.85), annotation_y * 1.05)) +
    labs(
      title    = plot_title,
      x        = "Position Relative to SmcR Binding Summit / Motif Center",
      y        = expression("Median Calibrated "*Log[2]*"(IP / Input)"),
      color    = "Condition",
      linetype = "Condition"
    )
  
  ggsave(file.path(OUT_DIR, out_filename), p_meta, width = fig_w, height = fig_h, dpi = 300)
  return(df_metagene)
}

message("[INFO] Generating Analysis 1B (2 Modes) and Analysis 1A (4 Tiers) TF Metagenes...")
df_tf_meta_2modes <- build_and_plot_tf_metagene(
  tf_summits_gr$Mode_Label,
  "SmcR TF Occupancy: Directly Repressed vs. Activated Loci (2 Modes)",
  "01B_TF_Calibrated_Log2FE_Metagene_2Modes.tiff",
  ncol_facet = 2, fig_w = 9.5, fig_h = 5.0
)

df_tf_meta_4tiers <- build_and_plot_tf_metagene(
  tf_summits_gr$Tier_Label,
  "SmcR TF Occupancy Across Functional Regulon Tiers (4 Groups)",
  "01A_TF_Calibrated_Log2FE_Metagene_4Tiers.tiff",
  ncol_facet = 2, fig_w = 10.0, fig_h = 7.5
)

# ==============================================================================
# 6. ANALYSIS 2A & 2B: PEAK WINDOW [-150, +150 bp] OCCUPANCY
#    (PAIRED POST-HOC WITH ASTERISK ANNOTATIONS)
# ==============================================================================
message("[INFO] Running Analysis 2A & 2B: Summit Window [-150, +150 bp] Paired Statistical Suite...")

peak_win_gr <- resize(tf_summits_gr, width = 2 * PEAK_QUANT_HALFWIN + 1, fix = "center")

get_tf_calibrated_region_log2fe <- function(ip_bw_paths, gr) {
  rowMeans(do.call(cbind, lapply(ip_bw_paths, function(ip_bw) {
    cal_row <- get_sample_calib_row(ip_bw)
    raw_lfc <- get_raw_log2_ip_over_input_regions(ip_bw, cal_row$ctrl_bw_path, gr)
    raw_lfc + cal_row$Log2_Calib_Offset
  })), na.rm = TRUE)
}

df_tf_peaks <- as.data.frame(mcols(tf_summits_gr)) %>%
  dplyr::mutate(
    Peak_wt_DMSO     = get_tf_calibrated_region_log2fe(tf_dmso_ip_bws,  peak_win_gr),
    Peak_wt_PTSP     = get_tf_calibrated_region_log2fe(tf_ptsp_ip_bws,  peak_win_gr),
    Peak_dsmcR_DMSO  = get_tf_calibrated_region_log2fe(tf_dsmcr_ip_bws, peak_win_gr),
    Peak_Shift_PTSP  = Peak_wt_PTSP    - Peak_wt_DMSO,
    Peak_Shift_dsmcR = Peak_dsmcR_DMSO - Peak_wt_DMSO
  )

evaluate_tf_peak_comparisons_pooled <- function(df_in, group_col, alpha = 0.05) {
  cond_levels <- c("wt + DMSO", "wt + PTSP", "dsmcR + DMSO")
  
  df_long <- df_in %>%
    dplyr::mutate(Group_Var = droplevels(factor(.data[[group_col]]))) %>%
    tidyr::pivot_longer(
      cols      = c(Peak_wt_DMSO, Peak_wt_PTSP, Peak_dsmcR_DMSO),
      names_to  = "Condition",
      values_to = "Calibrated_Log2_FE"
    ) %>%
    dplyr::mutate(
      Condition = factor(
        Condition,
        levels = c("Peak_wt_DMSO", "Peak_wt_PTSP", "Peak_dsmcR_DMSO"),
        labels = cond_levels
      )
    )
  
  grp_levels <- levels(df_long$Group_Var)
  
  pairwise_results <- list()
  summary_list     <- list()
  
  for (g in grp_levels) {
    sub_df  <- df_in %>% dplyr::filter(.data[[group_col]] == g)
    n_loci  <- nrow(sub_df)
    
    # wt PTSP vs wt DMSO (Paired)
    diff_wt <- sub_df$Peak_wt_PTSP - sub_df$Peak_wt_DMSO
    diff_norm_p <- if (n_loci >= 3 && n_loci <= 5000) shapiro.test(diff_wt)$p.value else 0.0
    is_diff_normal <- diff_norm_p >= alpha
    
    if (is_diff_normal) {
      t_res <- t.test(sub_df$Peak_wt_PTSP, sub_df$Peak_wt_DMSO, paired = TRUE)
      p_wt_ptsp_vs_dmso <- t_res$p.value
      test_wt_ptsp      <- "Paired Student's t-test"
    } else {
      w_res <- wilcox.test(sub_df$Peak_wt_PTSP, sub_df$Peak_wt_DMSO, paired = TRUE, exact = FALSE)
      p_wt_ptsp_vs_dmso <- w_res$p.value
      test_wt_ptsp      <- "Paired Wilcoxon Signed-Rank Test"
    }
    
    # dsmcR DMSO vs wt DMSO (Paired - exactly matched loci)
    diff_dsmc <- sub_df$Peak_dsmcR_DMSO - sub_df$Peak_wt_DMSO
    dsmc_norm_p <- if (n_loci >= 3 && n_loci <= 5000) shapiro.test(diff_dsmc)$p.value else 0.0
    is_dsmc_normal <- dsmc_norm_p >= alpha
    
    if (is_dsmc_normal) {
      t_dsmc <- t.test(sub_df$Peak_dsmcR_DMSO, sub_df$Peak_wt_DMSO, paired = TRUE)
      p_dsmc_vs_wt <- t_dsmc$p.value
      test_dsmc    <- "Paired Two-Sample t-test"
    } else {
      w_dsmc <- wilcox.test(sub_df$Peak_dsmcR_DMSO, sub_df$Peak_wt_DMSO, paired = TRUE, exact = FALSE)
      p_dsmc_vs_wt <- w_dsmc$p.value
      test_dsmc    <- "Paired Wilcoxon Signed-Rank Test"
    }
    
    pairwise_results[[g]] <- data.frame(
      Group           = g,
      Contrast        = c("wt + PTSP vs. wt + DMSO", "dsmcR + DMSO vs. wt + DMSO"),
      Raw_P           = c(p_wt_ptsp_vs_dmso, p_dsmc_vs_wt),
      Test_Used       = c(test_wt_ptsp, test_dsmc),
      stringsAsFactors = FALSE
    )
    
    summary_list[[g]] <- data.frame(
      Group                  = g,
      N_Peaks                = n_loci,
      Median_wt_DMSO         = median(sub_df$Peak_wt_DMSO, na.rm = TRUE),
      Median_wt_PTSP         = median(sub_df$Peak_wt_PTSP, na.rm = TRUE),
      Median_dsmcR_DMSO      = median(sub_df$Peak_dsmcR_DMSO, na.rm = TRUE),
      Median_Shift_PTSP      = median(sub_df$Peak_Shift_PTSP, na.rm = TRUE),
      Median_Shift_dsmcR     = median(sub_df$Peak_Shift_dsmcR, na.rm = TRUE),
      P_Val_wt_PTSP_vs_DMSO  = p_wt_ptsp_vs_dmso,
      P_Val_dsmcR_vs_wt      = p_dsmc_vs_wt,
      stringsAsFactors       = FALSE
    )
  }
  
  # Combine and apply Holm FWER correction across all tests
  pairwise_df <- dplyr::bind_rows(pairwise_results) %>%
    dplyr::mutate(
      Adjusted_P = p.adjust(Raw_P, method = "holm"),
      Significance = dplyr::case_when(
        Adjusted_P < 0.001 ~ "***",
        Adjusted_P < 0.01  ~ "**",
        Adjusted_P < 0.05  ~ "*",
        TRUE               ~ "ns"
      )
    )
  
  stats_table <- dplyr::bind_rows(summary_list) %>%
    dplyr::left_join(
      pairwise_df %>% 
        dplyr::select(Group, Contrast, Adjusted_P, Significance) %>% 
        tidyr::pivot_wider(names_from = Contrast, values_from = c(Adjusted_P, Significance)),
      by = "Group"
    )
  
  # Summary stats per cell for plotting bracket heights
  cell_summary_df <- df_long %>%
    dplyr::group_by(Group_Var, Condition) %>%
    dplyr::summarise(
      Max_Y  = max(Calibrated_Log2_FE, na.rm = TRUE),
      Q75_Y  = quantile(Calibrated_Log2_FE, 0.75, na.rm = TRUE),
      Median = median(Calibrated_Log2_FE, na.rm = TRUE),
      Mean   = mean(Calibrated_Log2_FE, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    stats_table       = stats_table,
    pairwise_table    = pairwise_df,
    cell_summary_df   = cell_summary_df,
    df_long           = df_long
  )
}

res_tf_2modes <- evaluate_tf_peak_comparisons_pooled(df_tf_peaks, "Mode_Label")
res_tf_4tiers <- evaluate_tf_peak_comparisons_pooled(df_tf_peaks, "Tier_Label")

tf_stats_2modes    <- res_tf_2modes$stats_table
tf_stats_4tiers    <- res_tf_4tiers$stats_table
tf_pairwise_2modes <- res_tf_2modes$pairwise_table
tf_pairwise_4tiers <- res_tf_4tiers$pairwise_table

message("\n[RESULT] Analysis 2B Post-Hoc Significance Summary (2 Modes):")
print(tf_stats_2modes %>% dplyr::select(Group, N_Peaks, `Significance_wt + PTSP vs. wt + DMSO`, `Significance_dsmcR + DMSO vs. wt + DMSO`))

message("\n[RESULT] Analysis 2A Post-Hoc Significance Summary (4 Tiers):")
print(tf_stats_4tiers %>% dplyr::select(Group, N_Peaks, `Significance_wt + PTSP vs. wt + DMSO`, `Significance_dsmcR + DMSO vs. wt + DMSO`))

plot_tf_peaks_stats <- function(eval_obj, plot_title, x_label, out_filename, fig_w = 7.2, fig_h = 5.8, angle_x = 0) {
  df_long  <- eval_obj$df_long
  cell_sum <- eval_obj$cell_summary_df
  pw_stats <- eval_obj$pairwise_table
  
  y_span  <- diff(range(df_long$Calibrated_Log2_FE, na.rm = TRUE))
  
  pw_map <- pw_stats %>%
    dplyr::mutate(Condition = dplyr::case_when(
      Contrast == "wt + PTSP vs. wt + DMSO"    ~ "wt + PTSP",
      Contrast == "dsmcR + DMSO vs. wt + DMSO" ~ "dsmcR + DMSO",
      TRUE                                     ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(Condition), Significance != "ns") %>%
    dplyr::rename(Group_Var = Group)
  
  sig_df <- cell_sum %>%
    dplyr::inner_join(pw_map, by = c("Group_Var", "Condition")) %>%
    dplyr::mutate(Label_Y = Max_Y + 0.05 * y_span)
  
  p_box <- ggplot(df_long, aes(x = Group_Var, y = Calibrated_Log2_FE, fill = Condition)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30", linewidth = 0.7, alpha = 0.6) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, position = position_dodge(width = 0.75)) +
    geom_point(
      aes(group = Condition),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      size = 1.7, alpha = 0.6, color = "grey15", show.legend = FALSE
    ) +
    geom_text(
      data     = sig_df,
      aes(x = Group_Var, y = Label_Y, label = Significance, group = Condition),
      position = position_dodge(width = 0.75),
      fontface = "bold",
      size     = 7,
      color    = "black",
      vjust    = 0
    ) +
    scale_fill_manual(
      values = c("wt + DMSO" = "#4E79A7", "wt + PTSP" = "#B00006", "dsmcR + DMSO" = "grey55"),
      labels = c("wt + DMSO" = "wt + DMSO", "wt + PTSP" = "wt + PTSP", "dsmcR + DMSO" = expression(Delta*italic("smcR")*" + DMSO"))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw(base_size = 13) +
    labs(
      title    = plot_title,
      subtitle = "Significance relative to wt + DMSO baseline (*p < 0.05, **p < 0.01, ***p < 0.001)",
      x        = x_label,
      y        = expression("Calibrated Peak "*Log[2]*"(IP / Input) [-150 to +150 bp]"),
      fill     = "Condition"
    ) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9.5, color = "grey25"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 0.5),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text.x      = element_text(face = "bold", angle = angle_x, hjust = ifelse(angle_x > 0, 1, 0.5)),
      axis.text.y      = element_text(face = "bold")
    )
  
  ggsave(file.path(OUT_DIR, out_filename), p_box, width = fig_w, height = fig_h, dpi = 300)
}

plot_tf_peaks_stats(
  res_tf_2modes,
  "SmcR Peak Occupancy [-150 to +150 bp]: Repressed vs. Activated (Paired Model)",
  "Regulation Mode",
  "02B_TF_Peak_Calibrated_Log2FE_Grouped_2Modes.tiff",
  fig_w = 7.2, fig_h = 5.8, angle_x = 0
)

plot_tf_peaks_stats(
  res_tf_4tiers,
  "SmcR Peak Occupancy [-150 to +150 bp] Across Functional Tiers (Paired Model)",
  "Functional Regulon Tier",
  "02A_TF_Peak_Calibrated_Log2FE_Grouped_4Tiers.tiff",
  fig_w = 8.5, fig_h = 5.8, angle_x = 15
)

# ==============================================================================
# 7. ANALYSIS 3A & 3B: DIFFERENCE METAGENE TRACKS (DELTA LOG2 IP/INPUT + 95% CI)
#    (ROUND 1: 4 TIERS | ROUND 2: 2 MODES)
# ==============================================================================
message(sprintf("[INFO] Running Analysis 3: Difference Metagene Tracks (%d Bootstrap Resamples)...", N_BOOT))

set.seed(42)

compute_diff_bootstrap <- function(mat_treat, mat_ctrl, x_pos, smooth_k, up_lim, dn_lim, n_boot = 500) {
  diff_mat <- mat_treat - mat_ctrl
  n_loci   <- nrow(diff_mat)
  
  smooth_row <- function(v) rollmean(v, k = smooth_k, fill = NA, align = "center")
  diff_mat_smooth <- t(apply(diff_mat, 1, smooth_row))
  if (n_loci == 1) diff_mat_smooth <- matrix(diff_mat_smooth, nrow = 1)
  
  keep_cols <- which(x_pos >= -up_lim & x_pos <= dn_lim)
  diff_crop <- diff_mat_smooth[, keep_cols, drop = FALSE]
  pos_crop  <- x_pos[keep_cols]
  
  obs_median <- apply(diff_crop, 2, median, na.rm = TRUE)
  
  if (n_loci >= 3) {
    boot_medians <- replicate(n_boot, {
      idx <- sample.int(n_loci, size = n_loci, replace = TRUE)
      apply(diff_crop[idx, , drop = FALSE], 2, median, na.rm = TRUE)
    })
    ci_low  <- apply(boot_medians, 1, quantile, probs = 0.025, na.rm = TRUE)
    ci_high <- apply(boot_medians, 1, quantile, probs = 0.975, na.rm = TRUE)
  } else {
    ci_low  <- obs_median
    ci_high <- obs_median
  }
  
  data.frame(
    Position = pos_crop,
    Delta    = obs_median,
    CI_Lower = ci_low,
    CI_Upper = ci_high
  )
}

build_and_plot_diff_tf_metagene <- function(group_vec, plot_title, out_filename, fig_w = 10, fig_h = 7.5) {
  grp_levels   <- levels(droplevels(factor(group_vec)))
  df_diff_meta <- data.frame()
  
  for (grp in grp_levels) {
    idx <- which(group_vec == grp)
    if (length(idx) == 0) next
    
    res_ptsp <- compute_diff_bootstrap(
      tf_mat_ptsp[idx, , drop = FALSE], tf_mat_dmso[idx, , drop = FALSE],
      x_positions, SMOOTH_K, UP_DIST, DOWN_DIST, N_BOOT
    ) %>%
      mutate(Contrast = "wt + PTSP vs. wt + DMSO", Group_Label = grp, N_Peaks = length(idx))
    
    res_dsmcr <- compute_diff_bootstrap(
      tf_mat_dsmcr[idx, , drop = FALSE], tf_mat_dmso[idx, , drop = FALSE],
      x_positions, SMOOTH_K, UP_DIST, DOWN_DIST, N_BOOT
    ) %>%
      mutate(Contrast = "dsmcR + DMSO vs. wt + DMSO", Group_Label = grp, N_Peaks = length(idx))
    
    df_diff_meta <- bind_rows(df_diff_meta, res_ptsp, res_dsmcr)
  }
  
  df_diff_meta <- df_diff_meta %>%
    mutate(
      Group_Label = factor(Group_Label, levels = grp_levels),
      Facet_Label = factor(
        paste0(Group_Label, " (n=", N_Peaks, " peaks)"),
        levels = unique(paste0(Group_Label, " (n=", N_Peaks, " peaks)"))
      ),
      Contrast    = factor(Contrast, levels = c("wt + PTSP vs. wt + DMSO", "dsmcR + DMSO vs. wt + DMSO"))
    )
  
  p_diff <- ggplot(df_diff_meta, aes(x = Position, y = Delta, color = Contrast, fill = Contrast, linetype = Contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7, alpha = 0.3) +
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.22, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.95) +
    coord_cartesian(ylim = c(-6, 0.5)) +
    facet_wrap(~ Facet_Label, ncol = 2) +
    scale_color_manual(
      values = c("wt + PTSP vs. wt + DMSO" = "#B00006", "dsmcR + DMSO vs. wt + DMSO" = "grey30"),
      labels = c(
        "wt + PTSP vs. wt + DMSO"    = "PTSP Effect (wt+PTSP - wt+DMSO)",
        "dsmcR + DMSO vs. wt + DMSO" = expression("Deletion Effect ("*Delta*italic("smcR")*" - wt+DMSO)")
      )
    ) +
    scale_fill_manual(
      values = c("wt + PTSP vs. wt + DMSO" = "#B00006", "dsmcR + DMSO vs. wt + DMSO" = "grey30")
    ) +
    scale_linetype_manual(
      values = c("wt + PTSP vs. wt + DMSO" = "solid", "dsmcR + DMSO vs. wt + DMSO" = "dashed"),
      labels = c(
        "wt + PTSP vs. wt + DMSO"    = "PTSP Effect (wt+PTSP - wt+DMSO)",
        "dsmcR + DMSO vs. wt + DMSO" = expression("Deletion Effect ("*Delta*italic("smcR")*" - wt+DMSO)")
      )
    ) +
    scale_x_continuous(
      breaks = c(-200, 0, 200),
      labels = c("-200 bp", "Binding Summit", "+200 bp"),
      expand = c(0, 0)
    ) +
    theme_bw(base_size = 14) +
    labs(
      title    = plot_title,
      x        = "Position Relative to SmcR Binding Summit / Motif Center",
      y        = expression(Delta*" Median Calibrated "*Log[2]*"(IP / Input)"),
      color    = "Contrast",
      linetype = "Contrast"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = 12),
      legend.position  = "bottom",
      axis.line        = element_line(colour = "black", linewidth = 0.5),
      axis.text        = element_text(face = "bold")
    )
  
  ggsave(file.path(OUT_DIR, out_filename), p_diff, width = fig_w, height = fig_h, dpi = 300)
}

build_and_plot_diff_tf_metagene(
  tf_summits_gr$Tier_Label,
  "Difference Metagene Tracks Across Functional Tiers (4 Groups, 95% Bootstrap CI)",
  "03A_TF_Difference_Metagene_Tracks_4Tiers.tiff",
  fig_w = 10, fig_h = 7.5
)

build_and_plot_diff_tf_metagene(
  tf_summits_gr$Mode_Label,
  "Difference Metagene Tracks: Repressed vs. Activated (2 Groups, 95% Bootstrap CI)",
  "03B_TF_Difference_Metagene_Tracks_2Modes.tiff",
  fig_w = 9.5, fig_h = 5.0
)

# ==============================================================================
# 7B. DIFFBIND DIFFERENTIAL BINDING ANALYSIS (FAST / TARGETED)
# ==============================================================================
message("[INFO] Running focused DiffBind differential binding analysis...")

diffbind_results_df <- data.frame()

tryCatch({
  if (!is.null(DIFFBIND_SAMPLESHEET) && file.exists(DIFFBIND_SAMPLESHEET)) {
    dba_samples <- read.csv(DIFFBIND_SAMPLESHEET, stringsAsFactors = FALSE)
  } else {
    all_bams  <- if (dir.exists(CHIP_BAM_DIR))  list.files(CHIP_BAM_DIR,  pattern = "\\.bam$", full.names = TRUE) else character(0)
    all_peaks <- if (dir.exists(CHIP_PEAK_DIR)) list.files(CHIP_PEAK_DIR, pattern = "\\.(narrowPeak|bed|xls)$", full.names = TRUE) else character(0)
    
    # Focus strictly on wt_DMSO and wt_PTSP samples for fast two-group differential testing
    target_calib <- bg_calib_params %>% filter(Condition %in% c("wt_DMSO", "wt_PTSP"))
    
    dba_samples <- lapply(seq_len(nrow(target_calib)), function(i) {
      raw_id    <- target_calib$sample_id[i]
      cond_grp  <- target_calib$Condition[i]
      core_pref <- sub("(_IP|_FE)?\\.bw$", "", raw_id, ignore.case = TRUE)
      rep_num   <- as.integer(coalesce(str_extract(core_pref, "(?<=rep)\\d+"), as.character(i)))
      
      ip_bam_hit <- all_bams[grepl(paste0("^", core_pref, "(_IP)?(\\.sorted|_sorted)?\\.bam$"), basename(all_bams), ignore.case = TRUE)]
      if (length(ip_bam_hit) == 0) ip_bam_hit <- all_bams[grepl(core_pref, basename(all_bams), ignore.case = TRUE) & !grepl("(Control|Input)", basename(all_bams), ignore.case = TRUE)]
      
      ctrl_pref    <- sub("(_Control|_Input)?\\.bw$", "", target_calib$paired_control[i], ignore.case = TRUE)
      ctrl_bam_hit <- all_bams[grepl(ctrl_pref, basename(all_bams), ignore.case = TRUE) & grepl("(Control|Input)", basename(all_bams), ignore.case = TRUE)]
      if (length(ctrl_bam_hit) == 0) ctrl_bam_hit <- all_bams[grepl(core_pref, basename(all_bams), ignore.case = TRUE) & grepl("(Control|Input)", basename(all_bams), ignore.case = TRUE)]
      
      peak_hit    <- all_peaks[grepl(core_pref, basename(all_peaks), ignore.case = TRUE)]
      peak_path   <- if (length(peak_hit) > 0) peak_hit[1] else PEAKS_FILE
      peak_caller <- if (grepl("narrowPeak$", peak_path, ignore.case = TRUE)) "narrow" else "bed"
      
      data.frame(
        SampleID   = core_pref,
        Tissue     = sapply(strsplit(cond_grp, "_"), `[`, 1),
        Factor     = "SmcR",
        Condition  = cond_grp,
        Treatment  = sapply(strsplit(cond_grp, "_"), `[`, 2),
        Replicate  = rep_num,
        bamReads   = if (length(ip_bam_hit) > 0) ip_bam_hit[1] else NA_character_,
        ControlID  = paste0(core_pref, "_Input"),
        bamControl = if (length(ctrl_bam_hit) > 0) ctrl_bam_hit[1] else NA_character_,
        Peaks      = peak_path,
        PeakCaller = peak_caller,
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows() %>% filter(!is.na(bamReads) & file.exists(bamReads))
  }
  
  if (nrow(dba_samples) < 2) {
    warning("[WARN] Fewer than 2 valid BAM files found for DiffBind. Skipping.")
  } else {
    if (all(is.na(dba_samples$bamControl) | !file.exists(dba_samples$bamControl))) {
      dba_samples$bamControl <- NULL
      dba_samples$ControlID  <- NULL
    }
    
    dba_obj <- dba(sampleSheet = dba_samples, config = data.frame(RunParallel = TRUE, cores = 4))
    
    # Fast C++ counting restricted strictly to your pre-defined target windows
    dba_obj <- dba.count(
      dba_obj,
      peaks                 = peak_win_gr,
      summits               = FALSE,
      score                 = DBA_SCORE_READS,
      bUseSummarizeOverlaps = FALSE,
      bParallel             = TRUE
    )
    
    # Fast library-size normalization (skips slow genome-wide background binning)
    dba_obj <- dba.normalize(dba_obj, method = DBA_DESEQ2, normalize = DBA_NORM_LIB, background = FALSE)
    
    # Single focused contrast: wt_PTSP vs wt_DMSO
    dba_obj <- dba.contrast(
      dba_obj,
      group1 = dba_obj$masks$wt_PTSP,
      group2 = dba_obj$masks$wt_DMSO,
      name1  = "wt_PTSP",
      name2  = "wt_DMSO"
    )
    
    dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2, bBlacklist = FALSE, bGreylist = FALSE)
    
    rep_gr <- dba.report(dba_obj, contrast = 1, th = 1, bCalled = FALSE, bCounts = TRUE)
    
    if (!is.null(rep_gr) && length(rep_gr) > 0) {
      hits <- findOverlaps(rep_gr, tf_summits_gr)
      rep_df <- as.data.frame(rep_gr) %>%
        mutate(
          Contrast   = "wt_PTSP_vs_wt_DMSO",
          peak_id    = NA_character_,
          gene_id    = NA_character_,
          Mode_Label = NA_character_,
          Tier_Label = NA_character_
        )
      if (length(hits) > 0) {
        rep_df$peak_id[queryHits(hits)]    <- as.character(tf_summits_gr$peak_id[subjectHits(hits)])
        rep_df$gene_id[queryHits(hits)]    <- as.character(tf_summits_gr$gene_id[subjectHits(hits)])
        rep_df$Mode_Label[queryHits(hits)] <- as.character(tf_summits_gr$Mode_Label[subjectHits(hits)])
        rep_df$Tier_Label[queryHits(hits)] <- as.character(tf_summits_gr$Tier_Label[subjectHits(hits)])
      }
      diffbind_results_df <- rep_df %>%
        relocate(Contrast, peak_id, gene_id, Mode_Label, Tier_Label)
      
      message(sprintf("[INFO] Focused DiffBind complete: %d peak intervals evaluated.", nrow(diffbind_results_df)))
    }
  }
}, error = function(e) {
  warning(sprintf("[WARN] DiffBind skipped due to error: %s", e$message))
})

# ==============================================================================
# 8. EXPORT CONSOLIDATED EXCEL WORKBOOK
# ==============================================================================
wb <- createWorkbook()

addWorksheet(wb, "Unique_Peak_Occupancy")
writeData(wb, "Unique_Peak_Occupancy", df_tf_peaks)

addWorksheet(wb, "Peak_Stats_CLD_2Modes")
writeData(wb, "Peak_Stats_CLD_2Modes", tf_stats_2modes)

addWorksheet(wb, "Peak_Stats_CLD_4Tiers")
writeData(wb, "Peak_Stats_CLD_4Tiers", tf_stats_4tiers)

addWorksheet(wb, "Peak_AllPairs_2Modes")
writeData(wb, "Peak_AllPairs_2Modes", tf_pairwise_2modes)

addWorksheet(wb, "Peak_AllPairs_4Tiers")
writeData(wb, "Peak_AllPairs_4Tiers", tf_pairwise_4tiers)

addWorksheet(wb, "Deep_BG_Calibration_Params")
writeData(wb, "Deep_BG_Calibration_Params", bg_calib_params)

addWorksheet(wb, "DiffBind_Results")
if (nrow(diffbind_results_df) > 0) {
  writeData(wb, "DiffBind_Results", diffbind_results_df)
} else {
  writeData(wb, "DiffBind_Results", data.frame(Note = "No DiffBind results generated. Verify CHIP_BAM_DIR or DIFFBIND_SAMPLESHEET paths."))
}

saveWorkbook(wb, file.path(OUT_DIR, "SmcR_TF_Calibrated_Occupancy_Summary.xlsx"), overwrite = TRUE)

message(paste0("\n[SUCCESS] SmcR TF ChIP-seq NCIS/SES calibrated analysis + DiffBind complete. Outputs saved to: ", OUT_DIR))
