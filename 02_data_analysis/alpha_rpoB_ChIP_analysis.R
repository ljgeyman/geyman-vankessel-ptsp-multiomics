#!/usr/bin/env Rscript

# ==============================================================================
# RNAP Dynamics Integration: Calibrated Log2(IP/Input)
# Dual-Round Execution: 4 Functional Tiers AND 2 Modes (Activation vs Repression)
# ------------------------------------------------------------------------------
# Core Analysis 1: Calibrated Log2(IP/Input) Metagene Profiling (4 Tiers & 2 Modes)
# Core Analysis 2: Promoter Occupancy [-100, +150 bp] Grouped Distributions (4 Tiers & 2 Modes)
# Core Analysis 3: Integration with RNA-seq (Calibrated RNAP Shift vs RNA LFC; 4 Tiers & 2 Modes)
# Follow-up 4:     Promoter [-100, +50 bp] vs Coding Body Decoupling (4 Tiers & 2 Modes)
# Follow-up 5:     Difference Metagene Tracks (Delta Calibrated Log2 IP/Input + 95% CI; 4 Tiers & 2 Modes)
# Follow-up 6:     Sequence Symmetry Integration (Motif Hamming Asymmetry vs RNAP Shift)
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
})

# ==============================================================================
# 1. PATHS & PARAMETERS
# ==============================================================================
# ---> GENOME & RNA-SEQ PATHS
ORIGINAL_REFSEQ_GTF <- "refs/refseq_vvul.gtf"
RNA_COUNTS          <- "analysis_results/counts/all_samples_counts.txt"
METADATA_FILE       <- "metadata.txt"

# ---> TF CHIP & MOTIF PATHS (Symmetry Pipeline Integration)
PEAKS_FILE          <- "05_final_output/combined_macs2_peaks.txt"
FIMO_FILE           <- "PTSP_meme_chip_analysis/meme_chip_results/fimo_out_1/fimo.tsv"

# ---> RNAP LINEAR FOLD-ENRICHMENT (IP/Input) BIGWIG PATHS
RNAP_FE_BW_DIR      <- "04b_bigwigs"
FE_PSEUDOCOUNT      <- 1e-3  # Guard against 0s prior to log2 transformation

# ---> LOCI EXCLUSIONS (e.g. self-deleted genes or highly variable pseudo-genes)
EXCLUDED_LOCI       <- c("AOT11_RS09875", "AOT11_RS01725")

# ---> REGULON & MOTIF PARAMETERS
PADJ_CUTOFF         <- 0.05
LFC_CUTOFF          <- 1.0
MOTIF_LENGTH        <- 21
PEAK_MAXGAP         <- 2000

# ---> METAGENE & PROMOTER WINDOW PARAMETERS
UP_DIST             <- 400
DOWN_DIST           <- 400
SMOOTH_K            <- 51   # Odd integer for centered rolling mean window (bp)

# Analysis 2 Promoter Window: [-100, +150] relative to TSS/ATG
PRO_UP_A2           <- 100
PRO_DOWN_A2         <- 150

# Follow-up 4 Decoupling Windows: Core Promoter [-100, +50] vs Downstream Body [+51, End]
DEC_PRO_UP          <- 100
DEC_PRO_DOWN        <- 50
MIN_DEC_GENE_LEN    <- 200  # Minimum gene length (bp) for promoter vs body decoupling

# Bootstrap iterations for Follow-up 5 Difference Metagene 95% CIs
N_BOOT              <- 1000

OUT_DIR             <- "07_rnap_dynamics_calibrated"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Helper to fail fast if BigWig regex patterns don't match any files
assert_bws_exist <- function(paths, label) {
  if (length(paths) == 0) {
    stop(sprintf("[FATAL] No BigWig files found for '%s'. Check directory path and filename regex.", label), call. = FALSE)
  }
  message(sprintf("[INFO] Found %d BigWig file(s) for %s.", length(paths), label))
}

# Helper for exact Hamming Distance (Fixed-length mismatch count)
hamming_distance <- function(s1, s2) {
  if (is.na(s1) || is.na(s2) || nchar(s1) != nchar(s2)) return(NA_integer_)
  sum(strsplit(s1, "")[[1]] != strsplit(s2, "")[[1]])
}

# ==============================================================================
# 2. ANNOTATIONS, DESEQ2 & TPM CALCULATION
# ==============================================================================
message("[INFO] Parsing GTF annotations and anchoring at RefSeq CDS start (ATG)...")

ref_gtf <- import(ORIGINAL_REFSEQ_GTF)
ref_genes <- ref_gtf[ref_gtf$type == "gene"]

raw_ids <- if ("locus_tag" %in% colnames(mcols(ref_genes))) ref_genes$locus_tag else ref_genes$gene_id
if (is(raw_ids, "List") || is(raw_ids, "CharacterList")) raw_ids <- unlist(as.list(raw_ids))
ref_genes$locus_tag <- str_remove(as.character(raw_ids), "^gene-")

# Anchor directly to the 5' start of the annotated gene (ATG / start codon)
final_tss_coords <- ifelse(strand(ref_genes) == "+", start(ref_genes), end(ref_genes))

tss_gr <- GRanges(
  seqnames   = seqnames(ref_genes),
  ranges     = IRanges(start = final_tss_coords, end = final_tss_coords),
  strand     = strand(ref_genes),
  locus_tag  = ref_genes$locus_tag,
  gene_start = start(ref_genes),
  gene_end   = end(ref_genes)
)

# Full gene body GRanges for TPM calibration curve
all_genes_body_gr <- ref_genes
names(all_genes_body_gr) <- ref_genes$locus_tag

message("[INFO] Importing RNA-seq counts, computing TPMs, and running DESeq2...")
raw_counts_full <- suppressWarnings(read_tsv(RNA_COUNTS, comment = "#", show_col_types = FALSE))
colnames(raw_counts_full) <- colnames(raw_counts_full) %>%
  str_remove("^.*bam[/\\\\.]") %>%
  str_remove("_sorted.*$")

gene_lengths_kb <- raw_counts_full$Length / 1000.0
names(gene_lengths_kb) <- raw_counts_full$Geneid

count_matrix <- raw_counts_full %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# Compute sample-level and condition-level TPMs for calibration
rpk_matrix <- sweep(count_matrix, 1, gene_lengths_kb[rownames(count_matrix)], FUN = "/")
tpm_matrix <- sweep(rpk_matrix, 2, colSums(rpk_matrix, na.rm = TRUE) / 1e6, FUN = "/")

# Align sample metadata (supports metadata.txt or automatic regex parsing)
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

# Compute condition-specific mean log2(TPM + 1)
tpm_df <- data.frame(
  gene_id              = rownames(tpm_matrix),
  TPM_WT_DMSO          = rowMeans(tpm_matrix[, metadata$group == "dluxo_dmso", drop = FALSE], na.rm = TRUE),
  TPM_WT_PTSP          = rowMeans(tpm_matrix[, metadata$group == "dluxo_ptsp", drop = FALSE], na.rm = TRUE),
  TPM_dsmcR_DMSO       = rowMeans(tpm_matrix[, metadata$group == "dsmcr_dmso", drop = FALSE], na.rm = TRUE)
) %>%
  mutate(
    log2TPM_WT_DMSO    = log2(TPM_WT_DMSO + 1),
    log2TPM_WT_PTSP    = log2(TPM_WT_PTSP + 1),
    log2TPM_dsmcR_DMSO = log2(TPM_dsmcR_DMSO + 1),
    log2TPM_Global     = log2(rowMeans(across(c(TPM_WT_DMSO, TPM_WT_PTSP, TPM_dsmcR_DMSO)), na.rm = TRUE) + 1)
  )

# Run DESeq2 with symmetry pipeline filtering (rowSums >= 10)
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
  left_join(tpm_df, by = "gene_id")

# ==============================================================================
# 3. SYMMETRY PIPELINE TARGET LOCUS MAPPING (BOTH 4-TIER & 2-MODE CLASSIFICATIONS)
# ==============================================================================
message("[INFO] Defining 4-Tier Regulon and 2-Mode Classification using Symmetry Pipeline logic...")

rna_res_all <- rna_res_all %>%
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
    Response_Magnitude     = abs(log2FoldChange_Inhib),
    Concordant_Sensitivity = log2FoldChange_Inhib * sign(log2FoldChange_dSmcR),
    Baseline_Magnitude     = abs(log2FoldChange_dSmcR)
  )

# Parse FIMO motifs and calculate asymmetry scores
fimo <- read_tsv(FIMO_FILE, comment = "#", show_col_types = FALSE)

calculate_global_asymmetry <- function(sequence) {
  if (is.na(sequence) || nchar(sequence) < MOTIF_LENGTH) return(NA_integer_)
  seq_upper  <- toupper(sequence)
  seq_len    <- nchar(seq_upper)
  half_len   <- floor(seq_len / 2)
  left_half  <- substr(seq_upper, 1, half_len)
  right_half <- substr(seq_upper, (seq_len - half_len) + 1, seq_len)
  right_rc   <- as.character(reverseComplement(DNAString(right_half)))
  return(hamming_distance(left_half, right_rc))
}

calculate_palindromic_anchor_asymmetry <- function(sequence) {
  if (is.na(sequence) || nchar(sequence) < 18) return(NA_integer_)
  seq_upper    <- toupper(sequence)
  left_anchor  <- substr(seq_upper, 4, 7)
  right_anchor <- substr(seq_upper, 15, 18)
  right_rc     <- as.character(reverseComplement(DNAString(right_anchor)))
  return(hamming_distance(left_anchor, right_rc))
}

calculate_meme_anchor_asymmetry <- function(sequence) {
  if (is.na(sequence) || nchar(sequence) < 17) return(NA_integer_)
  seq_upper    <- toupper(sequence)
  left_anchor  <- substr(seq_upper, 4, 7)
  right_anchor <- substr(seq_upper, 14, 17)
  right_rc     <- as.character(reverseComplement(DNAString(right_anchor)))
  return(hamming_distance(left_anchor, right_rc))
}

fimo <- fimo %>%
  filter(!is.na(matched_sequence)) %>%
  rowwise() %>%
  mutate(
    Global_Asymmetry             = calculate_global_asymmetry(matched_sequence),
    Palindromic_Anchor_Asymmetry = calculate_palindromic_anchor_asymmetry(matched_sequence),
    MEME_Anchor_Asymmetry        = calculate_meme_anchor_asymmetry(matched_sequence)
  ) %>%
  ungroup()

fimo_gr <- makeGRangesFromDataFrame(
  fimo,
  keep.extra.columns = TRUE,
  seqnames.field     = "sequence_name",
  start.field        = "start",
  end.field          = "stop"
)

# Load MACS2 peaks and map to genes within 2000 bp (maxgap = 2000)
peaks <- read.delim(PEAKS_FILE, stringsAsFactors = FALSE)
if (!"peak_id" %in% colnames(peaks)) colnames(peaks)[4] <- "peak_id"

peaks_gr <- makeGRangesFromDataFrame(
  peaks,
  keep.extra.columns = TRUE,
  seqnames.field     = "chrom",
  start.field        = "chromStart",
  end.field          = "chromEnd"
)

gene_overlaps <- findOverlaps(peaks_gr, ref_genes, maxgap = PEAK_MAXGAP)
annotated_peaks <- peaks[queryHits(gene_overlaps), ]
annotated_peaks$gene_id <- ref_genes$locus_tag[subjectHits(gene_overlaps)]

motif_overlaps <- findOverlaps(peaks_gr, fimo_gr)
overlap_df <- data.frame(
  peak_id                      = peaks$peak_id[queryHits(motif_overlaps)],
  Global_Asymmetry             = fimo_gr$Global_Asymmetry[subjectHits(motif_overlaps)],
  Palindromic_Anchor_Asymmetry = fimo_gr$Palindromic_Anchor_Asymmetry[subjectHits(motif_overlaps)],
  MEME_Anchor_Asymmetry        = fimo_gr$MEME_Anchor_Asymmetry[subjectHits(motif_overlaps)],
  motif_pvalue                 = fimo_gr$`p-value`[subjectHits(motif_overlaps)],
  matched_sequence             = fimo_gr$matched_sequence[subjectHits(motif_overlaps)]
)

best_motifs <- overlap_df %>%
  group_by(peak_id) %>%
  slice_min(order_by = motif_pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

peaks_motifs <- inner_join(annotated_peaks, best_motifs, by = "peak_id")

# Master Regulated Dataset
master_df <- rna_res_all %>%
  filter(Mode != "Independent") %>%
  inner_join(peaks_motifs, by = "gene_id") %>%
  group_by(gene_id) %>%
  filter(!gene_id %in% EXCLUDED_LOCI) %>%
  slice_min(order_by = motif_pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Motif_Score_log10 = -log10(motif_pvalue),
    Tier_Label        = droplevels(Tier_Label),
    Mode_Label        = droplevels(Mode_Label)
  )

# Strict Promoter-Proximal Collapsed Dataset 
master_df_unique <- master_df %>%
  group_by(peak_id) %>%
  slice_max(order_by = Baseline_Magnitude, n = 1, with_ties = FALSE) %>%
  ungroup()

# Attach coordinates to regulated GRanges (anchored at TSS/ATG)
regulated_gr <- tss_gr[tss_gr$locus_tag %in% master_df_unique$gene_id]
matched_idx  <- match(regulated_gr$locus_tag, master_df_unique$gene_id)
mcols(regulated_gr) <- cbind(
  mcols(regulated_gr),
  master_df_unique[matched_idx, c("peak_id", "Mode", "Mode_Label", "Tier", "Tier_Label",
                                  "log2FoldChange_dSmcR", "padj_dSmcR",
                                  "log2FoldChange_Inhib", "padj_Inhib",
                                  "Response_Magnitude", "Baseline_Magnitude",
                                  "Global_Asymmetry", "Palindromic_Anchor_Asymmetry",
                                  "MEME_Anchor_Asymmetry", "Motif_Score_log10", "matched_sequence")]
)

message(sprintf("[INFO] Symmetry pipeline mapped %d regulated genes across %d unique peaks.",
                nrow(master_df), nrow(master_df_unique)))
message("\n--- 4-Tier Breakdown ---")
print(table(regulated_gr$Tier_Label))
message("\n--- 2-Mode Breakdown ---")
print(table(regulated_gr$Mode_Label))

# ==============================================================================
# 4. LOCATE FE BIGWIGS & BUILD CROSS-SAMPLE TPM CALIBRATION CURVE
# ==============================================================================
message("[INFO] Locating Linear Fold-Enrichment (IP/Input) BigWigs in ", RNAP_FE_BW_DIR, "...")

rnap_dmso_fe_bws  <- list.files(RNAP_FE_BW_DIR, pattern = "(wt).*DMSO.*[.]bw$", full.names = TRUE)
rnap_ptsp_fe_bws  <- list.files(RNAP_FE_BW_DIR, pattern = "(wt).*PTSP.*[.]bw$", full.names = TRUE)
rnap_dsmcr_fe_bws <- list.files(RNAP_FE_BW_DIR, pattern = "smcR.*DMSO.*[.]bw$", full.names = TRUE)

assert_bws_exist(rnap_dmso_fe_bws,  "Linear FE WT/dluxO DMSO")
assert_bws_exist(rnap_ptsp_fe_bws,  "Linear FE WT/dluxO PTSP")
assert_bws_exist(rnap_dsmcr_fe_bws, "Linear FE dsmcR DMSO")

all_fe_bws <- c(rnap_dmso_fe_bws, rnap_ptsp_fe_bws, rnap_dsmcr_fe_bws)

# Define Non-Regulated / Housekeeping Calibration Gene Set
bound_genes_all <- unique(annotated_peaks$gene_id)

non_reg_df <- rna_res_all %>%
  filter(
    !gene_id %in% bound_genes_all,
    !gene_id %in% EXCLUDED_LOCI,
    abs(log2FoldChange_dSmcR) < 0.5,
    abs(log2FoldChange_Inhib) < 0.5,
    !is.na(log2TPM_Global)
  )

if (nrow(non_reg_df) < 200) {
  non_reg_df <- rna_res_all %>%
    filter(Mode == "Independent", !gene_id %in% bound_genes_all, !is.na(log2TPM_Global), !gene_id %in% EXCLUDED_LOCI)
}

non_reg_gr <- all_genes_body_gr[names(all_genes_body_gr) %in% non_reg_df$gene_id]
non_reg_df <- non_reg_df[match(names(non_reg_gr), non_reg_df$gene_id), ]

message(sprintf("[INFO] Building continuous linear calibration curve across %d non-regulated genes...",
                length(non_reg_gr)))

get_single_bw_raw_log2fe <- function(bw_path, gr) {
  s <- import(bw_path, which = gr)
  hits <- findOverlaps(gr, s)
  overlaps <- pintersect(gr[queryHits(hits)], s[subjectHits(hits)])
  
  df_ov <- data.frame(
    idx     = queryHits(hits),
    log2_fe = log2(pmax(s$score[subjectHits(hits)], FE_PSEUDOCOUNT)),
    w       = width(overlaps)
  ) %>%
    group_by(idx) %>%
    summarise(mean_val = weighted.mean(log2_fe, w = w, na.rm = TRUE), .groups = "drop") %>%
    complete(idx = seq_along(gr), fill = list(mean_val = NA_real_))
  
  return(df_ov$mean_val)
}

calib_params <- lapply(all_fe_bws, function(bw) {
  cond_group <- case_when(
    bw %in% rnap_dmso_fe_bws  ~ "WT_DMSO",
    bw %in% rnap_ptsp_fe_bws  ~ "WT_PTSP",
    bw %in% rnap_dsmcr_fe_bws ~ "dsmcR_DMSO",
    TRUE                      ~ "Unknown"
  )
  
  x_tpm <- switch(
    cond_group,
    "WT_DMSO"    = non_reg_df$log2TPM_WT_DMSO,
    "WT_PTSP"    = non_reg_df$log2TPM_WT_PTSP,
    "dsmcR_DMSO" = non_reg_df$log2TPM_dsmcR_DMSO,
    non_reg_df$log2TPM_Global
  )
  
  y_raw_log2fe <- get_single_bw_raw_log2fe(bw, non_reg_gr)
  fit <- lm(y_raw_log2fe ~ x_tpm)
  
  data.frame(
    bw_path    = bw,
    sample_id  = basename(bw),
    Condition  = cond_group,
    alpha_i    = unname(coef(fit)[1]),
    beta_i     = unname(coef(fit)[2]),
    r_squared  = summary(fit)$r.squared,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

alpha_bar <- mean(calib_params$alpha_i, na.rm = TRUE)
beta_bar  <- mean(calib_params$beta_i,  na.rm = TRUE)

calib_params <- calib_params %>%
  mutate(
    alpha_ref    = alpha_bar,
    beta_ref     = beta_bar,
    scale_factor = beta_bar / beta_i,
    shift_offset = alpha_bar - (beta_bar / beta_i) * alpha_i
  )

message("\n[INFO] Cross-Sample TPM Calibration Parameters:")
print(calib_params %>% dplyr::select(sample_id, Condition, alpha_i, beta_i, r_squared, scale_factor, shift_offset))
write.xlsx(calib_params, file.path(OUT_DIR, "TPM_Calibration_Curve_Parameters.xlsx"), overwrite = TRUE)

calib_plot_df <- lapply(seq_len(nrow(calib_params)), function(i) {
  bw   <- calib_params$bw_path[i]
  cond <- calib_params$Condition[i]
  x_tpm <- switch(
    cond,
    "WT_DMSO"    = non_reg_df$log2TPM_WT_DMSO,
    "WT_PTSP"    = non_reg_df$log2TPM_WT_PTSP,
    "dsmcR_DMSO" = non_reg_df$log2TPM_dsmcR_DMSO,
    non_reg_df$log2TPM_Global
  )
  y_raw <- get_single_bw_raw_log2fe(bw, non_reg_gr)
  y_cal <- calib_params$shift_offset[i] + calib_params$scale_factor[i] * y_raw
  
  bind_rows(
    data.frame(Sample = calib_params$sample_id[i], Condition = cond, Stage = "1. Uncalibrated Log2(IP/Input)", log2TPM = x_tpm, Signal = y_raw),
    data.frame(Sample = calib_params$sample_id[i], Condition = cond, Stage = "2. Calibrated Log2(IP/Input)",   log2TPM = x_tpm, Signal = y_cal)
  )
}) %>% bind_rows()

p_calib <- ggplot(calib_plot_df, aes(x = log2TPM, y = Signal, color = Sample, linetype = Condition)) +
  geom_point(alpha = 0.08, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~ Stage, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Cross-Sample Reference Calibration Curve Across Non-Regulated Genes",
    x     = "RNA-seq Expression: Log2(TPM + 1)",
    y     = "Gene Body RNAP ChIP: Log2(IP/Input)"
  ) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(OUT_DIR, "00_TPM_Calibration_Curves_Diagnostic.tiff"), p_calib, width = 10, height = 5, dpi = 300)

get_calib_coeffs <- function(bw_path) {
  row <- calib_params[calib_params$bw_path == bw_path, ]
  if (nrow(row) == 0) stop(paste("BigWig not found in calibration table:", bw_path))
  list(scale = row$scale_factor[1], offset = row$shift_offset[1])
}

# ==============================================================================
# 5. CORE ANALYSIS 1: CALIBRATED LOG2(IP/INPUT) METAGENE PROFILING
#    (ROUND 1: 4 TIERS | ROUND 2: 2 MODES)
# ==============================================================================
message("[INFO] Extracting per-locus calibrated metagene matrices once for all regulated genes...")

PAD_BP      <- floor(SMOOTH_K / 2)
EXT_UP      <- UP_DIST + PAD_BP
EXT_DOWN    <- DOWN_DIST + PAD_BP
x_positions <- seq(-EXT_UP, EXT_DOWN)

get_calibrated_locus_matrix <- function(bw_paths, gr, up_dist, dn_dist) {
  total_bins <- up_dist + dn_dist + 1
  n_loci     <- length(gr)
  if (length(bw_paths) == 0 || n_loci == 0) {
    return(matrix(0, nrow = n_loci, ncol = total_bins))
  }
  
  promoter_windows <- promoters(gr, upstream = up_dist, downstream = dn_dist + 1)
  accum_mat <- matrix(0, nrow = n_loci, ncol = total_bins)
  
  for (bw in bw_paths) {
    cal      <- get_calib_coeffs(bw)
    bw_data  <- import(bw)
    cov_data <- coverage(bw_data, weight = "score")
    
    cal_cov_data <- lapply(cov_data, function(chr_rle) {
      raw_log2 <- log2(pmax(chr_rle, FE_PSEUDOCOUNT))
      cal$offset + cal$scale * raw_log2
    })
    
    mat_list <- lapply(seq_along(promoter_windows), function(i) {
      chrom      <- as.character(seqnames(promoter_windows)[i])
      start_pos  <- start(promoter_windows)[i]
      end_pos    <- end(promoter_windows)[i]
      strand_dir <- as.character(strand(promoter_windows)[i])
      
      if (!chrom %in% names(cal_cov_data)) return(rep(0, total_bins))
      
      chr_cov <- cal_cov_data[[chrom]]
      chr_len <- length(chr_cov)
      
      s <- max(1L, start_pos)
      e <- min(chr_len, end_pos)
      
      if (s <= e) {
        vec <- as.numeric(chr_cov[s:e])
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
  
  return(accum_mat / length(bw_paths))
}

# Compute the full gene x base-pair matrices once across all regulated_gr
full_mat_dmso  <- get_calibrated_locus_matrix(rnap_dmso_fe_bws,  regulated_gr, EXT_UP, EXT_DOWN)
full_mat_ptsp  <- get_calibrated_locus_matrix(rnap_ptsp_fe_bws,  regulated_gr, EXT_UP, EXT_DOWN)
full_mat_dsmcr <- get_calibrated_locus_matrix(rnap_dsmcr_fe_bws, regulated_gr, EXT_UP, EXT_DOWN)

# Modular function to build and plot calibrated metagenes for any grouping column
build_and_plot_metagene <- function(group_vec, plot_title, out_filename, ncol_facet = 2, fig_w = 10, fig_h = 7.5) {
  grp_levels  <- levels(droplevels(factor(group_vec)))
  df_metagene <- data.frame()
  
  for (grp in grp_levels) {
    idx <- which(group_vec == grp)
    if (length(idx) == 0) next
    
    dmso_prof  <- apply(full_mat_dmso[idx, , drop = FALSE],  2, median, na.rm = TRUE)
    ptsp_prof  <- apply(full_mat_ptsp[idx, , drop = FALSE],  2, median, na.rm = TRUE)
    dsmcr_prof <- apply(full_mat_dsmcr[idx, , drop = FALSE], 2, median, na.rm = TRUE)
    
    df_metagene <- bind_rows(
      df_metagene,
      data.frame(Position = x_positions, Signal = dmso_prof,  Plot_Group = "WT_DMSO",    Group_Label = grp, N_Loci = length(idx)),
      data.frame(Position = x_positions, Signal = ptsp_prof,  Plot_Group = "WT_PTSP",    Group_Label = grp, N_Loci = length(idx)),
      data.frame(Position = x_positions, Signal = dsmcr_prof, Plot_Group = "dsmcR_DMSO", Group_Label = grp, N_Loci = length(idx))
    )
  }
  
  df_metagene <- df_metagene %>%
    group_by(Plot_Group, Group_Label) %>%
    arrange(Position) %>%
    mutate(Signal = rollmean(Signal, k = SMOOTH_K, fill = NA, align = "center")) %>%
    ungroup() %>%
    filter(!is.na(Signal) & Position >= -UP_DIST & Position <= DOWN_DIST) %>%
    mutate(
      Plot_Group  = factor(Plot_Group, levels = c("WT_DMSO", "WT_PTSP", "dsmcR_DMSO")),
      Group_Label = factor(Group_Label, levels = grp_levels),
      Facet_Label = factor(
        paste0(Group_Label, " (n=", N_Loci, ")"),
        levels = unique(paste0(Group_Label, " (n=", N_Loci, ")"))
      )
    )
  
  y_max        <- max(df_metagene$Signal, na.rm = TRUE)
  annotation_y <- y_max * 1.05
  gene_model_y <- min(df_metagene$Signal, na.rm = TRUE) * 0.9
  
  p_meta <- ggplot(df_metagene, aes(x = Position, y = Signal, group = Plot_Group, color = Plot_Group, linetype = Plot_Group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8, alpha = 0.3) +
    geom_ribbon(aes(ymin = -Inf, ymax = Signal, fill = Plot_Group), alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line(data = subset(df_metagene, Plot_Group == "dsmcR_DMSO"), linewidth = 0.8, alpha = 0.85) +
    geom_line(data = subset(df_metagene, Plot_Group != "dsmcR_DMSO"), linewidth = 0.9) +
    facet_wrap(~ Facet_Label, ncol = ncol_facet, scales = "free_y") +
    scale_color_manual(
      values = c("WT_DMSO" = "#4E79A7", "WT_PTSP" = "#B00006", "dsmcR_DMSO" = "grey35"),
      labels = c("WT_DMSO" = "WT + DMSO", "WT_PTSP" = "WT + PTSP", "dsmcR_DMSO" = expression(Delta*italic("smcR")*" + DMSO"))
    ) +
    scale_fill_manual(
      values = c("WT_DMSO" = "#4E79A7", "WT_PTSP" = "#B00006", "dsmcR_DMSO" = "grey35")
    ) +
    scale_linetype_manual(
      values = c("WT_DMSO" = "solid", "WT_PTSP" = "solid", "dsmcR_DMSO" = "dashed"),
      labels = c("WT_DMSO" = "WT + DMSO", "WT_PTSP" = "WT + PTSP", "dsmcR_DMSO" = expression(Delta*italic("smcR")*" + DMSO"))
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
      labels = c("-200 bp", "ORF Start", "+200 bp"),
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(min(0, gene_model_y * 0.85), annotation_y * 1.05)) +
    labs(
      title    = plot_title,
      x        = "Position Relative to ORF Start",
      y        = expression("Median Calibrated "*Log[2]*"(IP / Input)"),
      color    = "Condition",
      linetype = "Condition"
    )
  
  ggsave(file.path(OUT_DIR, out_filename), p_meta, width = fig_w, height = fig_h, dpi = 300)
  return(df_metagene)
}

message("[INFO] Running Analysis 1 (Round 1: 4 Tiers & Round 2: 2 Modes)...")
df_meta_4tiers <- build_and_plot_metagene(
  regulated_gr$Tier_Label,
  "Calibrated RNAP Occupancy Across Functional Regulon Tiers (4 Groups)",
  "01A_Calibrated_Log2FE_Metagene_4Tiers.tiff",
  ncol_facet = 2, fig_w = 10, fig_h = 7.5
)

df_meta_2modes <- build_and_plot_metagene(
  regulated_gr$Mode_Label,
  "Calibrated RNAP Occupancy: Directly Repressed vs. Activated (2 Groups)",
  "01B_Calibrated_Log2FE_Metagene_2Modes.tiff",
  ncol_facet = 2, fig_w = 9.5, fig_h = 5.0
)

# ==============================================================================
# 6. CORE ANALYSIS 2: PROMOTER [-100, +150 bp] CALIBRATED LOG2(IP/INPUT)
#    (POOLED OMNIBUS & PAIRED/UNPAIRED POST-HOC WITH ASTERISK ANNOTATIONS)
# ==============================================================================
message("[INFO] Running Analysis 2: Pooled Paired/Unpaired Tests with Asterisk Annotations...")

tss_pos <- start(regulated_gr)
strands <- as.character(strand(regulated_gr))

pro_start_a2 <- pmax(1L, ifelse(strands == "+", tss_pos - PRO_UP_A2,   tss_pos - PRO_DOWN_A2))
pro_end_a2   <- ifelse(strands == "+",          tss_pos + PRO_DOWN_A2, tss_pos + PRO_UP_A2)

promoter_a2_gr <- regulated_gr
ranges(promoter_a2_gr) <- IRanges(start = pro_start_a2, end = pro_end_a2)

get_calibrated_region_log2fe <- function(bw_paths, gr) {
  rowMeans(do.call(cbind, lapply(bw_paths, function(bw) {
    cal       <- get_calib_coeffs(bw)
    raw       <- get_single_bw_raw_log2fe(bw, gr)
    cal_val   <- cal$offset + cal$scale * raw
    cal_floor <- cal$offset + cal$scale * log2(FE_PSEUDOCOUNT)
    ifelse(is.na(cal_val), cal_floor, cal_val)
  })), na.rm = TRUE)
}

df_promoter <- as.data.frame(mcols(regulated_gr)) %>%
  dplyr::mutate(
    gene_id          = locus_tag,
    Prom_WT_DMSO     = get_calibrated_region_log2fe(rnap_dmso_fe_bws,  promoter_a2_gr),
    Prom_WT_PTSP     = get_calibrated_region_log2fe(rnap_ptsp_fe_bws,  promoter_a2_gr),
    Prom_dsmcR_DMSO  = get_calibrated_region_log2fe(rnap_dsmcr_fe_bws, promoter_a2_gr),
    Prom_Shift_PTSP  = Prom_WT_PTSP - Prom_WT_DMSO,
    Prom_Shift_dsmcR = Prom_dsmcR_DMSO - Prom_WT_DMSO
  )

evaluate_promoter_comparisons_pooled <- function(df_in, group_col, alpha = 0.05) {
  cond_levels <- c("WT + DMSO", "WT + PTSP", "dsmcR + DMSO")
  
  df_long <- df_in %>%
    dplyr::mutate(Group_Var = droplevels(factor(.data[[group_col]]))) %>%
    tidyr::pivot_longer(
      cols      = c(Prom_WT_DMSO, Prom_WT_PTSP, Prom_dsmcR_DMSO),
      names_to  = "Condition",
      values_to = "Calibrated_Log2_FE"
    ) %>%
    dplyr::mutate(
      Condition = factor(
        Condition,
        levels = c("Prom_WT_DMSO", "Prom_WT_PTSP", "Prom_dsmcR_DMSO"),
        labels = cond_levels
      )
    )
  
  grp_levels <- levels(df_long$Group_Var)
  
  # 1. Perform within-group paired tests (WT + PTSP vs. WT + DMSO) and unpaired for dsmcR
  pairwise_results <- list()
  summary_list     <- list()
  
  for (g in grp_levels) {
    sub_df  <- df_in %>% dplyr::filter(.data[[group_col]] == g)
    n_loci  <- nrow(sub_df)
    
    # WT PTSP vs WT DMSO (Paired)
    diff_wt <- sub_df$Prom_WT_PTSP - sub_df$Prom_WT_DMSO
    diff_norm_p <- if (n_loci >= 3 && n_loci <= 5000) shapiro.test(diff_wt)$p.value else 0.0
    is_diff_normal <- diff_norm_p >= alpha
    
    if (is_diff_normal) {
      t_res <- t.test(sub_df$Prom_WT_PTSP, sub_df$Prom_WT_DMSO, paired = TRUE)
      p_wt_ptsp_vs_dmso <- t_res$p.value
      test_wt_ptsp      <- "Paired Student's t-test"
    } else {
      w_res <- wilcox.test(sub_df$Prom_WT_PTSP, sub_df$Prom_WT_DMSO, paired = TRUE, exact = FALSE)
      p_wt_ptsp_vs_dmso <- w_res$p.value
      test_wt_ptsp      <- "Paired Wilcoxon Signed-Rank Test"
    }
    
    # dsmcR DMSO vs WT DMSO (Paired)
    diff_dsmc <- sub_df$Prom_dsmcR_DMSO - sub_df$Prom_WT_DMSO
    dsmc_norm_p <- if (n_loci >= 3 && n_loci <= 5000) shapiro.test(diff_dsmc)$p.value else 0.0
    is_dsmc_normal <- dsmc_norm_p >= alpha
    
    if (is_dsmc_normal) {
      t_dsmc <- t.test(sub_df$Prom_dsmcR_DMSO, sub_df$Prom_WT_DMSO, paired = TRUE)
      p_dsmc_vs_wt <- t_dsmc$p.value
      test_dsmc    <- "Paired Two-Sample t-test"
    } else {
      w_dsmc <- wilcox.test(sub_df$Prom_dsmcR_DMSO, sub_df$Prom_WT_DMSO, paired = TRUE, exact = FALSE)
      p_dsmc_vs_wt <- w_dsmc$p.value
      test_dsmc    <- "Paired Wilcoxon Signed-Rank Test"
    }
    
    pairwise_results[[g]] <- data.frame(
      Group           = g,
      Contrast        = c("WT + PTSP vs. WT + DMSO", "dsmcR + DMSO vs. WT + DMSO"),
      Raw_P           = c(p_wt_ptsp_vs_dmso, p_dsmc_vs_wt),
      Test_Used       = c(test_wt_ptsp, test_dsmc),
      stringsAsFactors = FALSE
    )
    
    summary_list[[g]] <- data.frame(
      Group                  = g,
      N_Loci                 = n_loci,
      Median_WT_DMSO         = median(sub_df$Prom_WT_DMSO, na.rm = TRUE),
      Median_WT_PTSP         = median(sub_df$Prom_WT_PTSP, na.rm = TRUE),
      Median_dsmcR_DMSO      = median(sub_df$Prom_dsmcR_DMSO, na.rm = TRUE),
      P_Val_WT_PTSP_vs_DMSO  = p_wt_ptsp_vs_dmso,
      P_Val_dsmcR_vs_WT      = p_dsmc_vs_wt,
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

res_a2_4tiers <- evaluate_promoter_comparisons_pooled(df_promoter, "Tier_Label")
res_a2_2modes <- evaluate_promoter_comparisons_pooled(df_promoter, "Mode_Label")

prom_stats_4tiers    <- res_a2_4tiers$stats_table
prom_stats_2modes    <- res_a2_2modes$stats_table

message("\n[RESULT] Analysis 2A Post-Hoc Significance Summary (4 Tiers):")
print(prom_stats_4tiers %>% dplyr::select(Group, N_Loci, `Significance_WT + PTSP vs. WT + DMSO`, `Significance_dsmcR + DMSO vs. WT + DMSO`))

message("\n[RESULT] Analysis 2B Post-Hoc Significance Summary (2 Modes):")
print(prom_stats_2modes %>% dplyr::select(Group, N_Loci, `Significance_WT + PTSP vs. WT + DMSO`, `Significance_dsmcR + DMSO vs. WT + DMSO`))

# Plotting Function refactored for asterisks instead of CLD
plot_promoter_grouped_stats <- function(eval_obj, plot_title, x_label, out_filename, fig_w = 8.5, fig_h = 5.8, angle_x = 15) {
  df_long  <- eval_obj$df_long
  cell_sum <- eval_obj$cell_summary_df
  pw_stats <- eval_obj$pairwise_table
  
  y_span  <- diff(range(df_long$Calibrated_Log2_FE, na.rm = TRUE))
  
  # Map pairwise significance brackets exactly to the box being compared against WT+DMSO
  pw_map <- pw_stats %>%
    dplyr::mutate(Condition = dplyr::case_when(
      Contrast == "WT + PTSP vs. WT + DMSO"    ~ "WT + PTSP",
      Contrast == "dsmcR + DMSO vs. WT + DMSO" ~ "dsmcR + DMSO",
      TRUE                                     ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(Condition), Significance != "ns") %>%
    dplyr::rename(Group_Var = Group)
  
  # Inner join guarantees text is only plotted over significant boxes
  sig_df <- cell_sum %>%
    dplyr::inner_join(pw_map, by = c("Group_Var", "Condition")) %>%
    dplyr::mutate(Label_Y = Max_Y + 0.05 * y_span)
  
  p_box <- ggplot(df_long, aes(x = Group_Var, y = Calibrated_Log2_FE, fill = Condition)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7, alpha = 0.5) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, position = position_dodge(width = 0.75)) +
    geom_point(
      aes(group = Condition),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      size = 1.6, alpha = 0.55, color = "grey15", show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("WT + DMSO" = "#4E79A7", "WT + PTSP" = "#B00006", "dsmcR + DMSO" = "grey55"),
      labels = c("WT + DMSO" = "WT + DMSO", "WT + PTSP" = "WT + PTSP", "dsmcR + DMSO" = expression(Delta*italic("smcR")*" + DMSO"))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw(base_size = 13) +
    labs(
      title    = plot_title,
      subtitle = "Significance relative to WT + DMSO baseline (*p < 0.05, **p < 0.01, ***p < 0.001)",
      x        = x_label,
      y        = expression("Calibrated Promoter "*Log[2]*"(IP / Input)"),
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

plot_promoter_grouped_stats(
  res_a2_4tiers,
  "Promoter RNAP Occupancy [-100 to +150 bp] Across Functional Tiers (Paired Model)",
  "Functional Regulon Tier",
  "02A_Promoter_Calibrated_Log2FE_Grouped_4Tiers.tiff",
  fig_w = 8.5, fig_h = 5.8, angle_x = 15
)

plot_promoter_grouped_stats(
  res_a2_2modes,
  "Promoter RNAP Occupancy [-100 to +150 bp]: Repressed vs. Activated (Paired Model)",
  "Regulation Mode",
  "02B_Promoter_Calibrated_Log2FE_Grouped_2Modes.tiff",
  fig_w = 7.2, fig_h = 5.8, angle_x = 0
)

# ==============================================================================
# 7. CORE ANALYSIS 3: INTEGRATION WITH RNA-SEQ (CALIBRATED SHIFT VS RNA LFC)
#    (ROUND 1: 4 TIERS | ROUND 2: 2 MODES)
# ==============================================================================
message("[INFO] Running Analysis 3: Integrating Calibrated Promoter Log2(IP/Input) Shift with RNA-seq...")

df_integ_ptsp <- df_promoter %>%
  mutate(
    RNAP_Calibrated_Shift = Prom_Shift_PTSP,
    RNA_Log2FC            = log2FoldChange_Inhib,
    Comparison            = "PTSP Effect (WT + Drug vs DMSO)"
  )

df_integ_del <- df_promoter %>%
  mutate(
    RNAP_Calibrated_Shift = Prom_Shift_dsmcR,
    RNA_Log2FC            = log2FoldChange_dSmcR,
    Comparison            = "Deletion Effect (dsmcR vs WT)"
  )

df_integration <- bind_rows(df_integ_ptsp, df_integ_del) %>%
  filter(is.finite(RNAP_Calibrated_Shift) & is.finite(RNA_Log2FC))

tier_palette <- c(
  "Repressed: Sensitive" = "#D55E00",
  "Repressed: Muted"     = "#E69F00",
  "Activated: Sensitive" = "#009E73",
  "Activated: Muted"     = "#56B4E9"
)

mode_palette <- c(
  "Directly Repressed by SmcR" = "firebrick",
  "Directly Activated by SmcR" = "forestgreen"
)

# Plot 3A: 4 Tiers
p_integ_4tiers <- ggplot(df_integ_ptsp, aes(x = RNAP_Calibrated_Shift, y = RNA_Log2FC, color = Tier_Label, shape = Tier_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(aes(fill = Tier_Label), method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = tier_palette, name = "Regulon Tier") +
  scale_fill_manual(values = tier_palette,  name = "Regulon Tier") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Regulon Tier") +
  theme_bw(base_size = 13) +
  labs(
    title = "Promoter RNAP Shift vs. Transcript Abundance (PTSP - 4 Tiers)",
    x     = expression(Delta*" Calibrated Promoter "*Log[2]*"(IP / Input) [PTSP - DMSO]"),
    y     = "RNA-seq Expression Shift (Log2FC PTSP / DMSO)"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text        = element_text(face = "bold"),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black")
  ) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

ggsave(file.path(OUT_DIR, "03A_RNAP_vs_RNAseq_Integration_4Tiers.tiff"), p_integ_4tiers, width = 6.5, height = 6, dpi = 300)

# Plot 3B: 2 Modes
p_integ_2modes <- ggplot(df_integ_ptsp, aes(x = RNAP_Calibrated_Shift, y = RNA_Log2FC, color = Mode_Label, shape = Mode_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(aes(fill = Mode_Label), method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.9) +
  scale_color_manual(values = mode_palette, name = "Regulation Mode") +
  scale_fill_manual(values = mode_palette,  name = "Regulation Mode") +
  theme_bw(base_size = 13) +
  labs(
    title = "Promoter RNAP Shift vs. Transcript Abundance (PTSP - 2 Modes)",
    x     = expression(Delta*" Calibrated Promoter "*Log[2]*"(IP / Input) [PTSP - DMSO]"),
    y     = "RNA-seq Expression Shift (Log2FC PTSP / DMSO)"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text        = element_text(face = "bold"),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(OUT_DIR, "03B_RNAP_vs_RNAseq_Integration_2Modes.tiff"), p_integ_2modes, width = 6.0, height = 5.5, dpi = 300)

compute_cor_stats <- function(df_in, group_col) {
  df_in %>%
    group_by(Comparison, Group = .data[[group_col]]) %>%
    summarise(
      N_Loci    = n(),
      Pearson_R = cor(RNAP_Calibrated_Shift, RNA_Log2FC, method = "pearson", use = "complete.obs"),
      P_value   = ifelse(n() >= 3, cor.test(RNAP_Calibrated_Shift, RNA_Log2FC, method = "pearson")$p.value, NA_real_),
      R_squared = ifelse(n() >= 3, summary(lm(RNA_Log2FC ~ RNAP_Calibrated_Shift))$r.squared, NA_real_),
      .groups   = "drop"
    )
}

cor_stats_4tiers <- compute_cor_stats(df_integration, "Tier_Label")
cor_stats_2modes <- compute_cor_stats(df_integration, "Mode_Label")

message("\n[RESULT] Analysis 3A RNAP vs RNA-seq Correlations (4 Tiers):")
print(cor_stats_4tiers)
message("\n[RESULT] Analysis 3B RNAP vs RNA-seq Correlations (2 Modes):")
print(cor_stats_2modes)

# ==============================================================================
# 8. FOLLOW-UP 4: PROMOTER [-100, +50 bp] VS CODING BODY DECOUPLING (UNPAIRED)
#    (ROUND 1: 4 TIERS | ROUND 2: 2 Modes)
# ==============================================================================
message("[INFO] Running Follow-up 4: Promoter [-100, +50 bp] vs Coding Body [+51, End] Unpaired Decoupling...")

gene_lens <- abs(regulated_gr$gene_end - regulated_gr$gene_start) + 1
dec_gr    <- regulated_gr[gene_lens >= MIN_DEC_GENE_LEN]

tss_dec     <- start(dec_gr)
strands_dec <- as.character(strand(dec_gr))

core_pro_start <- pmax(1L, ifelse(strands_dec == "+", tss_dec - DEC_PRO_UP,   tss_dec - DEC_PRO_DOWN))
core_pro_end   <- ifelse(strands_dec == "+",          tss_dec + DEC_PRO_DOWN, tss_dec + DEC_PRO_UP)

cds_body_start <- ifelse(strands_dec == "+", tss_dec + DEC_PRO_DOWN + 1L, dec_gr$gene_start)
cds_body_end   <- ifelse(strands_dec == "+", dec_gr$gene_end,             tss_dec - DEC_PRO_DOWN - 1L)

core_pro_gr <- dec_gr
ranges(core_pro_gr) <- IRanges(start = core_pro_start, end = core_pro_end)

cds_body_gr <- dec_gr
ranges(cds_body_gr) <- IRanges(start = cds_body_start, end = cds_body_end)

df_decoupling <- as.data.frame(mcols(dec_gr)) %>%
  mutate(
    gene_id               = locus_tag,
    CorePro_DMSO          = get_calibrated_region_log2fe(rnap_dmso_fe_bws,  core_pro_gr),
    CorePro_PTSP          = get_calibrated_region_log2fe(rnap_ptsp_fe_bws,  core_pro_gr),
    CorePro_dsmcR         = get_calibrated_region_log2fe(rnap_dsmcr_fe_bws, core_pro_gr),
    CDSBody_DMSO          = get_calibrated_region_log2fe(rnap_dmso_fe_bws,  cds_body_gr),
    CDSBody_PTSP          = get_calibrated_region_log2fe(rnap_ptsp_fe_bws,  cds_body_gr),
    CDSBody_dsmcR         = get_calibrated_region_log2fe(rnap_dsmcr_fe_bws, cds_body_gr),
    Delta_Pro_PTSP        = CorePro_PTSP  - CorePro_DMSO,
    Delta_Body_PTSP       = CDSBody_PTSP  - CDSBody_DMSO,
    Delta_Pro_dsmcR       = CorePro_dsmcR - CorePro_DMSO,
    Delta_Body_dsmcR      = CDSBody_dsmcR - CDSBody_DMSO,
    Clearance_Shift_PTSP  = Delta_Body_PTSP  - Delta_Pro_PTSP,
    Clearance_Shift_dsmcR = Delta_Body_dsmcR - Delta_Pro_dsmcR
  )

plot_decoupling <- function(df_in, group_col, plot_title, out_filename, fig_w = 9, fig_h = 7) {
  df_dec_long <- df_in %>%
    mutate(Group_Var = .data[[group_col]]) %>%
    pivot_longer(
      cols      = c(Delta_Pro_PTSP, Delta_Body_PTSP),
      names_to  = "Compartment",
      values_to = "Delta_Calibrated_Log2FE"
    ) %>%
    mutate(
      Compartment = factor(
        Compartment,
        levels = c("Delta_Pro_PTSP", "Delta_Body_PTSP"),
        labels = c("Core Promoter\n[-100, +50 bp]", "Coding Body\n[+51 bp, End]")
      )
    )
  
  p_dec <- ggplot(df_dec_long, aes(x = Compartment, y = Delta_Calibrated_Log2FE, fill = Compartment)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(width = 0.45, alpha = 0.85, outlier.shape = NA) +
    geom_jitter(width = 0.14, size = 1.8, alpha = 0.65, color = "grey15") +
    facet_wrap(~ Group_Var, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c("Core Promoter\n[-100, +50 bp]" = "#4E79A7", "Coding Body\n[+51 bp, End]" = "#F28E2B")) +
    theme_bw(base_size = 13) +
    labs(
      title = plot_title,
      x     = "Genomic Region",
      y     = expression(Delta*" Calibrated "*Log[2]*"(IP / Input) [PTSP - DMSO]")
    ) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      strip.text      = element_text(face = "bold"),
      legend.position = "none",
      axis.text       = element_text(face = "bold")
    )
  
  ggsave(file.path(OUT_DIR, out_filename), p_dec, width = fig_w, height = fig_h, dpi = 300)
}

plot_decoupling(
  df_decoupling, "Tier_Label",
  "Promoter vs. Coding Body RNAP Decoupling Across Functional Tiers (4 Groups)",
  "04A_Promoter_vs_Body_Decoupling_4Tiers.tiff",
  fig_w = 9, fig_h = 7
)

plot_decoupling(
  df_decoupling, "Mode_Label",
  "Promoter vs. Coding Body RNAP Decoupling: Repressed vs. Activated (2 Groups)",
  "04B_Promoter_vs_Body_Decoupling_2Modes.tiff",
  fig_w = 8.5, fig_h = 4.5
)

compute_decoupling_stats <- function(df_in, group_col, alpha = 0.05) {
  df_in %>%
    group_by(Group = .data[[group_col]]) %>%
    summarise(
      N_Loci                  = n(),
      Median_Delta_Pro_PTSP   = median(Delta_Pro_PTSP, na.rm = TRUE),
      Median_Delta_Body_PTSP  = median(Delta_Body_PTSP, na.rm = TRUE),
      Median_Clearance_Diff   = median(Clearance_Shift_PTSP, na.rm = TRUE),
      Cor_Pro_vs_Body         = cor(Delta_Pro_PTSP, Delta_Body_PTSP, use = "complete.obs"),
      Shapiro_Min_P           = min(
        ifelse(n() >= 3, shapiro.test(Delta_Pro_PTSP)$p.value, 1.0),
        ifelse(n() >= 3, shapiro.test(Delta_Body_PTSP)$p.value, 1.0),
        na.rm = TRUE
      ),
      Selected_Test           = ifelse(Shapiro_Min_P >= alpha, "Unpaired Two-Sample t-test", "Unpaired Mann-Whitney U"),
      P_Unpaired_Pro_vs_Body  = ifelse(
        Shapiro_Min_P >= alpha,
        t.test(Delta_Body_PTSP, Delta_Pro_PTSP, paired = FALSE)$p.value,
        wilcox.test(Delta_Body_PTSP, Delta_Pro_PTSP, paired = FALSE, exact = FALSE)$p.value
      ),
      .groups                 = "drop"
    )
}

decoupling_stats_4tiers <- compute_decoupling_stats(df_decoupling, "Tier_Label")
decoupling_stats_2modes <- compute_decoupling_stats(df_decoupling, "Mode_Label")

message("\n[RESULT] Follow-up 4A Promoter vs. Body Unpaired Decoupling (4 Tiers):")
print(decoupling_stats_4tiers)
message("\n[RESULT] Follow-up 4B Promoter vs. Body Unpaired Decoupling (2 Modes):")
print(decoupling_stats_2modes)

# ==============================================================================
# 9. FOLLOW-UP 5: DIFFERENCE METAGENE TRACKS (DELTA LOG2 IP/INPUT + 95% BOOTSTRAP CI)
#    (ROUND 1: 4 TIERS | ROUND 2: 2 MODES)
# ==============================================================================
message(sprintf("[INFO] Running Follow-up 5: Difference Metagene Tracks (%d Bootstrap Resamples)...", N_BOOT))

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

build_and_plot_diff_metagene <- function(group_vec, plot_title, out_filename, fig_w = 10, fig_h = 7.5) {
  grp_levels   <- levels(droplevels(factor(group_vec)))
  df_diff_meta <- data.frame()
  
  for (grp in grp_levels) {
    idx <- which(group_vec == grp)
    if (length(idx) == 0) next
    
    res_ptsp <- compute_diff_bootstrap(
      full_mat_ptsp[idx, , drop = FALSE], full_mat_dmso[idx, , drop = FALSE],
      x_positions, SMOOTH_K, UP_DIST, DOWN_DIST, N_BOOT
    ) %>%
      mutate(Contrast = "WT + PTSP vs. WT + DMSO", Group_Label = grp, N_Loci = length(idx))
    
    res_dsmcr <- compute_diff_bootstrap(
      full_mat_dsmcr[idx, , drop = FALSE], full_mat_dmso[idx, , drop = FALSE],
      x_positions, SMOOTH_K, UP_DIST, DOWN_DIST, N_BOOT
    ) %>%
      mutate(Contrast = "dsmcR + DMSO vs. WT + DMSO", Group_Label = grp, N_Loci = length(idx))
    
    df_diff_meta <- bind_rows(df_diff_meta, res_ptsp, res_dsmcr)
  }
  
  df_diff_meta <- df_diff_meta %>%
    mutate(
      Group_Label = factor(Group_Label, levels = grp_levels),
      Facet_Label = factor(
        paste0(Group_Label, " (n=", N_Loci, ")"),
        levels = unique(paste0(Group_Label, " (n=", N_Loci, ")"))
      ),
      Contrast    = factor(Contrast, levels = c("WT + PTSP vs. WT + DMSO", "dsmcR + DMSO vs. WT + DMSO"))
    )
  
  p_diff <- ggplot(df_diff_meta, aes(x = Position, y = Delta, color = Contrast, fill = Contrast, linetype = Contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7, alpha = 0.3) +
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.22, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.95) +
    coord_cartesian(ylim = c(-1, 1.25)) +
    facet_wrap(~ Facet_Label, ncol = 2) +
    scale_color_manual(
      values = c("WT + PTSP vs. WT + DMSO" = "#B00006", "dsmcR + DMSO vs. WT + DMSO" = "grey30"),
      labels = c(
        "WT + PTSP vs. WT + DMSO"    = "PTSP Effect (WT+PTSP - WT+DMSO)",
        "dsmcR + DMSO vs. WT + DMSO" = expression("Deletion Effect ("*Delta*italic("smcR")*" - WT+DMSO)")
      )
    ) +
    scale_fill_manual(
      values = c("WT + PTSP vs. WT + DMSO" = "#B00006", "dsmcR + DMSO vs. WT + DMSO" = "grey30")
    ) +
    scale_linetype_manual(
      values = c("WT + PTSP vs. WT + DMSO" = "solid", "dsmcR + DMSO vs. WT + DMSO" = "dashed"),
      labels = c(
        "WT + PTSP vs. WT + DMSO"    = "PTSP Effect (WT+PTSP - WT+DMSO)",
        "dsmcR + DMSO vs. WT + DMSO" = expression("Deletion Effect ("*Delta*italic("smcR")*" - WT+DMSO)")
      )
    ) +
    scale_x_continuous(
      breaks = c(-200, 0, 200),
      labels = c("-200 bp", "ORF Start", "+200 bp"),
      expand = c(0, 0)
    ) +
    theme_bw(base_size = 14) +
    labs(
      title    = plot_title,
      x        = "Position Relative to ORF Start",
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

build_and_plot_diff_metagene(
  regulated_gr$Tier_Label,
  "Difference Metagene Tracks Across Functional Tiers (4 Groups, 95% Bootstrap CI)",
  "05A_Difference_Metagene_Tracks_4Tiers.tiff",
  fig_w = 10, fig_h = 7.5
)

build_and_plot_diff_metagene(
  regulated_gr$Mode_Label,
  "Difference Metagene Tracks: Repressed vs. Activated (2 Groups, 95% Bootstrap CI)",
  "05B_Difference_Metagene_Tracks_2Modes.tiff",
  fig_w = 9.5, fig_h = 5.0
)

# ==============================================================================
# 10. FOLLOW-UP 6: INTEGRATING SEQUENCE ASYMMETRY WITH CALIBRATED RNAP OCCUPANCY
# ==============================================================================
message("[INFO] Running Follow-up 6: Testing Motif Asymmetry vs. Promoter RNAP Occupancy Shift...")

df_sym_rnap <- df_promoter %>%
  mutate(
    RNAP_Prom_Response_Mag = abs(Prom_Shift_PTSP),
    RNAP_Prom_Baseline_Mag = abs(Prom_Shift_dsmcR),
    RNAP_Concordant_Shift  = Prom_Shift_PTSP * sign(Prom_Shift_dsmcR)
  )

df_sym_rnap_unique <- df_sym_rnap %>%
  group_by(peak_id) %>%
  slice_max(order_by = Baseline_Magnitude, n = 1, with_ties = FALSE) %>%
  ungroup()

lmer_rnap_asym <- tryCatch(
  lmer(
    RNAP_Prom_Response_Mag ~ Global_Asymmetry * Mode + RNAP_Prom_Baseline_Mag + Motif_Score_log10 + (1 | peak_id),
    data = df_sym_rnap
  ),
  error = function(e) NULL
)

lm_rnap_asym_unique <- lm(
  RNAP_Prom_Response_Mag ~ Global_Asymmetry * Mode + RNAP_Prom_Baseline_Mag + Motif_Score_log10,
  data = df_sym_rnap_unique
)

lm_rnap_anchor_unique <- lm(
  RNAP_Prom_Response_Mag ~ Palindromic_Anchor_Asymmetry * Mode + RNAP_Prom_Baseline_Mag + Motif_Score_log10,
  data = df_sym_rnap_unique
)

message("\n[RESULT] Follow-up 6 Collapsed OLS Model: Motif Asymmetry Predicting Promoter RNAP Shift (N=92):")
print(summary(lm_rnap_asym_unique))

# Plot 6A: Colored by 4 Tiers
p_sym_rnap_4tiers <- ggplot(df_sym_rnap, aes(x = Global_Asymmetry, y = Prom_Shift_PTSP, color = Tier_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(width = 0.15, size = 2.6, alpha = 0.85) +
  geom_smooth(aes(group = Mode, fill = Mode), method = "lm", se = TRUE, alpha = 0.18, color = "black", linewidth = 0.8) +
  facet_wrap(~ Mode_Label, scales = "free_y") +
  scale_color_manual(values = tier_palette, name = "Regulon Tier") +
  scale_fill_manual(values = c("Repression" = "#D55E00", "Activation" = "#009E73"), guide = "none") +
  theme_bw(base_size = 13) +
  labs(
    title = "Motif Hamming Asymmetry vs. Calibrated Promoter RNAP Shift (4 Tiers)",
    x     = "Global Motif Hamming Asymmetry (Mismatches Across Half-Sites)",
    y     = expression(Delta*" Calibrated Promoter "*Log[2]*"(IP / Input) [PTSP - DMSO]")
  ) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text       = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "06A_Sequence_Asymmetry_vs_Promoter_RNAP_Shift_4Tiers.tiff"), p_sym_rnap_4tiers, width = 8.5, height = 5.5, dpi = 300)

# Plot 6B: Colored purely by 2 Modes (Activation vs Repression)
p_sym_rnap_2modes <- ggplot(df_sym_rnap, aes(x = Global_Asymmetry, y = Prom_Shift_PTSP, color = Mode_Label, fill = Mode_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(width = 0.15, size = 2.6, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.18, linewidth = 0.9) +
  facet_wrap(~ Mode_Label, scales = "free_y") +
  scale_color_manual(values = mode_palette, name = "Regulation Mode") +
  scale_fill_manual(values = mode_palette,  name = "Regulation Mode") +
  theme_bw(base_size = 13) +
  labs(
    title = "Motif Hamming Asymmetry vs. Calibrated Promoter RNAP Shift (2 Modes)",
    x     = "Global Motif Hamming Asymmetry (Mismatches Across Half-Sites)",
    y     = expression(Delta*" Calibrated Promoter "*Log[2]*"(IP / Input) [PTSP - DMSO]")
  ) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text       = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "06B_Sequence_Asymmetry_vs_Promoter_RNAP_Shift_2Modes.tiff"), p_sym_rnap_2modes, width = 8.5, height = 5.5, dpi = 300)

extract_rnap_model <- function(mod, label) {
  if (is.null(mod)) return(NULL)
  td <- if (inherits(mod, "merMod")) {
    broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE)
  } else {
    broom::tidy(mod, conf.int = TRUE)
  }
  td %>%
    mutate(Model = label) %>%
    dplyr::select(Model, Term = term, Estimate = estimate, Std_Error = std.error, Statistic = statistic, P_value = p.value, CI_Low = conf.low, CI_High = conf.high)
}

rnap_asym_models_df <- bind_rows(
  extract_rnap_model(lmer_rnap_asym,        "1. LMM Operon-Controlled (All Regulated Loci)"),
  extract_rnap_model(lm_rnap_asym_unique,   "2. Collapsed OLS - Global Asymmetry (Unique Peaks)"),
  extract_rnap_model(lm_rnap_anchor_unique, "3. Collapsed OLS - Palindromic Anchor (Unique Peaks)")
)

# ==============================================================================
# 11. EXPORT CONSOLIDATED EXCEL WORKBOOK (BOTH 4-TIER & 2-MODE SHEETS)
# ==============================================================================
wb <- createWorkbook()

addWorksheet(wb, "Locus_Promoter_Data")
writeData(wb, "Locus_Promoter_Data", df_sym_rnap)

addWorksheet(wb, "Prom_Stats_CLD_4Tiers")
writeData(wb, "Prom_Stats_CLD_4Tiers", prom_stats_4tiers)

addWorksheet(wb, "Prom_Stats_CLD_2Modes")
writeData(wb, "Prom_Stats_CLD_2Modes", prom_stats_2modes)

addWorksheet(wb, "RNA_Correlations_4Tiers")
writeData(wb, "RNA_Correlations_4Tiers", cor_stats_4tiers)

addWorksheet(wb, "RNA_Correlations_2Modes")
writeData(wb, "RNA_Correlations_2Modes", cor_stats_2modes)

addWorksheet(wb, "Promoter_Body_Decoupling")
writeData(wb, "Promoter_Body_Decoupling", df_decoupling)

addWorksheet(wb, "Decoupling_Stats_4Tiers")
writeData(wb, "Decoupling_Stats_4Tiers", decoupling_stats_4tiers)

addWorksheet(wb, "Decoupling_Stats_2Modes")
writeData(wb, "Decoupling_Stats_2Modes", decoupling_stats_2modes)

addWorksheet(wb, "Asymmetry_RNAP_Models")
writeData(wb, "Asymmetry_RNAP_Models", rnap_asym_models_df)

addWorksheet(wb, "TPM_Calibration_Params")
writeData(wb, "TPM_Calibration_Params", calib_params)

saveWorkbook(wb, file.path(OUT_DIR, "RNAP_Calibrated_Dynamics_Complete_Summary.xlsx"), overwrite = TRUE)

message(paste0("\n[SUCCESS] Both 4-Tier and 2-Mode analyses completed. Outputs saved to: ", OUT_DIR))
