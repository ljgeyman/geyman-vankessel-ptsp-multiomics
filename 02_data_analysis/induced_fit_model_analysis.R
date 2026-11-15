#!/usr/bin/env Rscript

# ==============================================================================
# Operon-Aware Modeling & Figure Generation
# ==============================================================================
if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
if (!requireNamespace("broom.mixed", quietly = TRUE)) install.packages("broom.mixed")

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(Biostrings)
  library(ggplot2)
  library(ashr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggseqlogo)
  library(cowplot)
  library(tidyr)
  library(lme4)      # Linear Mixed-Effects Models
  library(lmerTest)  # Satterthwaite p-values for LMMs
  library(broom)
  library(broom.mixed)
})

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
CONFIG <- list(
  counts_file     = "analysis_results/counts/all_samples_counts.txt",
  metadata_file   = "metadata.txt",
  fimo_file       = "PTSP_meme_chip_analysis/meme_chip_results/fimo_out_1/fimo.tsv",
  peaks_file      = "05_final_output/combined_macs2_peaks.txt",
  gtf_file        = "refs/refseq_vvul.gtf",
  out_dir         = "09_sequence_symmetry",
  
  genotype_col    = "genotype",
  condition_col   = "condition",
  
  wt_group        = "dluxO_dmso",
  mut_group       = "dsmcR_dmso",
  inhib_group     = "dluxO_ptsp", 
  
  peaks_chr_col   = "chrom",
  peaks_start_col = "chromStart",
  peaks_end_col   = "chromEnd",
  
  padj_cutoff     = 0.05,
  lfc_cutoff      = 1.0,
  motif_length    = 21,
  
  # Excluded loci (e.g., self-deleted targets)
  exclude_loci    = c("AOT11_RS09875")
)

dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)

clean_ids <- function(ids) { return(ids) }

# Helper for exact Hamming Distance (Fixed-length mismatch count)
hamming_distance <- function(s1, s2) {
  if (is.na(s1) || is.na(s2) || nchar(s1) != nchar(s2)) return(NA_integer_)
  sum(strsplit(s1, "")[[1]] != strsplit(s2, "")[[1]])
}

# ==============================================================================
# 2. DATA IMPORT & ALIGNMENT
# ==============================================================================
raw_counts <- read.table(CONFIG$counts_file, header = TRUE, comment.char = "#")
count_data <- raw_counts[, 7:ncol(raw_counts)]
rownames(count_data) <- raw_counts$Geneid

colnames(count_data) <- colnames(count_data) %>%
  str_remove("^.*\\.bam\\.") %>%  
  str_remove("_sorted.*$")          

metadata <- read.table(CONFIG$metadata_file, header = FALSE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

missing_in_meta <- setdiff(colnames(count_data), rownames(metadata))
if (length(missing_in_meta) > 0) stop("Counts columns missing from metadata: ", paste(missing_in_meta, collapse = ", "))

metadata <- metadata[colnames(count_data), , drop = FALSE]
metadata[[CONFIG$genotype_col]]  <- sapply(strsplit(rownames(metadata), "_"), `[`, 1)
metadata[[CONFIG$condition_col]] <- sapply(strsplit(rownames(metadata), "_"), `[`, 2)
metadata$group <- factor(paste(metadata[[CONFIG$genotype_col]], metadata[[CONFIG$condition_col]], sep = "_"))

# ==============================================================================
# 3. DESEQ2 EXECUTION (~ group design)
# ==============================================================================
cat("Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_data), colData = metadata, design = ~ group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

get_shrunken_res <- function(dds_obj, num_group, den_group) {
  if (!all(c(num_group, den_group) %in% levels(dds_obj$group))) {
    stop(paste("CRITICAL ERROR: Group(s) not found in DESeq2 object:", num_group, "or", den_group))
  }
  res <- lfcShrink(dds_obj, contrast = c("group", num_group, den_group), type = "ashr", quiet = TRUE)
  as.data.frame(res) %>% rownames_to_column("Gene") %>% mutate(Gene = clean_ids(Gene))
}

res_dSmcR <- get_shrunken_res(dds, CONFIG$mut_group, CONFIG$wt_group)
res_inhib <- get_shrunken_res(dds, CONFIG$inhib_group, CONFIG$wt_group)

rna_res <- merge(res_dSmcR, res_inhib, by = "Gene", suffixes = c("_dSmcR", "_Inhib"))
colnames(rna_res)[1] <- "gene_id"

# ==============================================================================
# 4. DEFINE REGULON TIERS & CONTINUOUS METRICS (CORRECTED)
# ==============================================================================
rna_res <- rna_res %>%
  mutate(
    Mode = case_when(
      log2FoldChange_dSmcR > CONFIG$lfc_cutoff & padj_dSmcR < CONFIG$padj_cutoff ~ "Repression",
      log2FoldChange_dSmcR < -CONFIG$lfc_cutoff & padj_dSmcR < CONFIG$padj_cutoff ~ "Activation",
      TRUE ~ "Independent"
    ),
    # Simplify Muted logic: If it's Repression/Activation and NOT Sensitive, it's Muted.
    Tier = case_when(
      Mode == "Repression" & log2FoldChange_Inhib > CONFIG$lfc_cutoff & padj_Inhib < CONFIG$padj_cutoff ~ "Repression_Sensitive",
      Mode == "Repression" ~ "Repression_Muted",
      
      Mode == "Activation" & log2FoldChange_Inhib < -CONFIG$lfc_cutoff & padj_Inhib < CONFIG$padj_cutoff ~ "Activation_Sensitive",
      Mode == "Activation" ~ "Activation_Muted",
      
      TRUE ~ "Not_Direct_Target"
    ),
    Response_Magnitude = abs(log2FoldChange_Inhib),
    Concordant_Sensitivity = log2FoldChange_Inhib * sign(log2FoldChange_dSmcR),
    Baseline_Magnitude = abs(log2FoldChange_dSmcR)
  ) %>%
  filter(Mode != "Independent")

# ==============================================================================
# 5. DUAL METRIC SCORING (GLOBAL VS. EXPLICIT ANCHORS)
# ==============================================================================
cat("Calculating Motif Asymmetry Metrics...\n")
fimo <- read_tsv(CONFIG$fimo_file, comment = "#", show_col_types = FALSE)

calculate_global_asymmetry <- function(sequence) {
  if (is.na(sequence) || nchar(sequence) < CONFIG$motif_length) return(NA_integer_)
  seq_upper <- toupper(sequence)
  seq_len <- nchar(seq_upper)
  half_len <- floor(seq_len / 2)
  left_half <- substr(seq_upper, 1, half_len)
  right_half <- substr(seq_upper, (seq_len - half_len) + 1, seq_len)
  right_rc <- as.character(reverseComplement(DNAString(right_half)))
  return(hamming_distance(left_half, right_rc))
}

calculate_palindromic_anchor_asymmetry <- function(sequence) {
  if (is.na(sequence) || nchar(sequence) < 18) return(NA_integer_)
  seq_upper <- toupper(sequence)
  left_anchor  <- substr(seq_upper, 4, 7)
  right_anchor <- substr(seq_upper, 15, 18)
  right_rc <- as.character(reverseComplement(DNAString(right_anchor)))
  return(hamming_distance(left_anchor, right_rc))
}

calculate_meme_anchor_asymmetry <- function(sequence) {
  if (is.na(sequence) || nchar(sequence) < 17) return(NA_integer_)
  seq_upper <- toupper(sequence)
  left_anchor  <- substr(seq_upper, 4, 7)
  right_anchor <- substr(seq_upper, 14, 17)
  right_rc <- as.character(reverseComplement(DNAString(right_anchor)))
  return(hamming_distance(left_anchor, right_rc))
}

fimo <- fimo %>%
  filter(!is.na(matched_sequence)) %>%
  rowwise() %>%
  mutate(
    Global_Asymmetry = calculate_global_asymmetry(matched_sequence),
    Palindromic_Anchor_Asymmetry = calculate_palindromic_anchor_asymmetry(matched_sequence),
    MEME_Anchor_Asymmetry = calculate_meme_anchor_asymmetry(matched_sequence)
  ) %>%
  ungroup()

fimo_gr <- makeGRangesFromDataFrame(
  fimo,
  keep.extra.columns = TRUE,
  seqnames.field = "sequence_name",
  start.field = "start",
  end.field = "stop" 
)

# ==============================================================================
# 6. ANNOTATE PEAKS & CONSTRUCT DATASETS (MASTER N=146 & UNIQUE N=92)
# ==============================================================================
cat("Loading GTF and annotating peaks...\n")
peaks <- read.delim(CONFIG$peaks_file, stringsAsFactors = FALSE) 
if (!"peak_id" %in% colnames(peaks)) colnames(peaks)[4] <- "peak_id"

peaks_gr <- makeGRangesFromDataFrame(
  peaks, 
  keep.extra.columns = TRUE,
  seqnames.field = CONFIG$peaks_chr_col,
  start.field = CONFIG$peaks_start_col,
  end.field = CONFIG$peaks_end_col
)

gtf <- import(CONFIG$gtf_file)
genes <- gtf[gtf$type == "gene"]

raw_ids <- if ("locus_tag" %in% colnames(mcols(genes))) genes$locus_tag else genes$gene_id
if (is(raw_ids, "List") || is(raw_ids, "CharacterList")) raw_ids <- unlist(as.list(raw_ids))

gene_overlaps <- findOverlaps(peaks_gr, genes, maxgap = 2000)

annotated_peaks <- peaks[queryHits(gene_overlaps), ]
annotated_peaks$gene_id <- as.character(raw_ids[subjectHits(gene_overlaps)])
annotated_peaks$gene_id <- str_remove(annotated_peaks$gene_id, "^gene-")

motif_overlaps <- findOverlaps(peaks_gr, fimo_gr)

overlap_df <- data.frame(
  peak_id = peaks$peak_id[queryHits(motif_overlaps)],
  Global_Asymmetry = fimo_gr$Global_Asymmetry[subjectHits(motif_overlaps)],
  Palindromic_Anchor_Asymmetry = fimo_gr$Palindromic_Anchor_Asymmetry[subjectHits(motif_overlaps)],
  MEME_Anchor_Asymmetry = fimo_gr$MEME_Anchor_Asymmetry[subjectHits(motif_overlaps)],   motif_pvalue = fimo_gr$`p-value`[subjectHits(motif_overlaps)],
  matched_sequence = fimo_gr$matched_sequence[subjectHits(motif_overlaps)]
)

best_motifs <- overlap_df %>%
  group_by(peak_id) %>%
  slice_min(order_by = motif_pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

peaks_motifs <- inner_join(annotated_peaks, best_motifs, by = "peak_id")

# 1. Full Dataset (N = 146 genes across 92 unique peaks)
master_df <- inner_join(rna_res, peaks_motifs, by = "gene_id") %>%
  group_by(gene_id) %>%
  filter(!gene_id %in% CONFIG$exclude_loci) %>%
  slice_min(order_by = motif_pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Motif_Score_log10 = -log10(motif_pvalue))

# 2. Strict Promoter-Proximal Collapsed Dataset (N = 92 unique loci)
master_df_unique <- master_df %>%
  group_by(peak_id) %>%
  slice_max(order_by = Baseline_Magnitude, n = 1, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("\nDataset Summary: %d total genes mapped to %d unique binding loci.\n", 
            nrow(master_df), nrow(master_df_unique)))

# ==============================================================================
# 7. STATISTICAL MODEL SUITE
# ==============================================================================
cat("\n=== Fitting Regression Models ===\n")

# Model 1A: Naive OLS (N = 146)
lm_global <- lm(
  Response_Magnitude ~ Global_Asymmetry * Mode + Baseline_Magnitude + Motif_Score_log10 + log2(baseMean_dSmcR + 1), 
  data = master_df
)

# Model 1B: Linear Mixed-Effects Model (LMM) controlling for operon clustering (N = 146)
lmer_global <- lmer(
  Response_Magnitude ~ Global_Asymmetry * Mode + Baseline_Magnitude + Motif_Score_log10 + log2(baseMean_dSmcR + 1) + (1 | peak_id), 
  data = master_df
)

# Model 2A: Promoter-Proximal Collapsed OLS (N = 92) - Global Asymmetry
lm_global_unique <- lm(
  Response_Magnitude ~ Global_Asymmetry * Mode + Baseline_Magnitude + Motif_Score_log10 + log2(baseMean_dSmcR + 1), 
  data = master_df_unique
)

# Model 2B: Promoter-Proximal Collapsed OLS (N = 92) - Palindromic Anchor
lm_anchor_unique <- lm(
  Response_Magnitude ~ Palindromic_Anchor_Asymmetry * Mode + Baseline_Magnitude + Motif_Score_log10 + log2(baseMean_dSmcR + 1), 
  data = master_df_unique
)

cat("\n--- LMM Global Model (Polycistronic Operon Controlled) ---\n")
print(summary(lmer_global))

cat("\n--- Collapsed OLS Model (N=92 Unique Loci) ---\n")
print(summary(lm_global_unique))

# ==============================================================================
# 8. SLIDING WINDOW SCAN
# ==============================================================================
cat("\nExecuting Sliding Window Asymmetry Scan...\n")
window_size <- 6
scan_results <- data.frame()

for (start_pos in 1:(10 - window_size + 1)) {
  end_pos <- start_pos + window_size - 1
  right_start <- 22 - end_pos
  right_end   <- 22 - start_pos
  
  calc_dynamic_asymmetry <- function(seq) {
    if (is.na(seq) || nchar(seq) < 21) return(NA_integer_)
    seq_u <- toupper(seq)
    left_win  <- substr(seq_u, start_pos, end_pos)
    right_win <- substr(seq_u, right_start, right_end)
    right_rc  <- as.character(reverseComplement(DNAString(right_win)))
    return(hamming_distance(left_win, right_rc))
  }
  
  master_df$Temp_Asymmetry <- sapply(master_df$matched_sequence, calc_dynamic_asymmetry)
  
  lm_temp <- lm(
    Response_Magnitude ~ Temp_Asymmetry * Mode + Baseline_Magnitude + Motif_Score_log10 + log2(baseMean_dSmcR + 1), 
    data = master_df
  )
  
  summ <- summary(lm_temp)
  coefs <- summ$coefficients
  
  beta  <- if ("Temp_Asymmetry" %in% rownames(coefs)) coefs["Temp_Asymmetry", "Estimate"] else NA
  se    <- if ("Temp_Asymmetry" %in% rownames(coefs)) coefs["Temp_Asymmetry", "Std. Error"] else NA
  p_val <- if ("Temp_Asymmetry" %in% rownames(coefs)) coefs["Temp_Asymmetry", "Pr(>|t|)"] else NA
  
  scan_results <- rbind(scan_results, data.frame(
    Start_Pos         = start_pos,
    Left_Window       = paste0(start_pos, "-", end_pos),
    Right_Window      = paste0(right_start, "-", right_end),
    Asymmetry_Beta    = beta,
    Beta_SE           = se,
    Beta_CI_Lower     = beta - (1.96 * se),
    Beta_CI_Upper     = beta + (1.96 * se),
    Asymmetry_P_Value = p_val,
    Model_Adj_R2      = summ$adj.r.squared
  ))
}

# ==============================================================================
# 9. FIGURE PANEL GENERATION
# ==============================================================================
cat("\nBuilding Figure Panels...\n")

# PANEL A: Target Locus Collapse vs Baseline Knockout
panel_a <- ggplot(master_df, aes(x = log2FoldChange_dSmcR, y = log2FoldChange_Inhib, color = Global_Asymmetry)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_c(name = "Hamming\nAsymmetry") +
  theme_minimal(base_size = 11) +
  labs(
    title = "A. Response Scaling vs. Baseline Knockout Magnitude",
    x = "Deletion Effect",
    y = "Inhibitor Effect"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# PANEL B: Sliding Window Scan Track Plot
panel_b <- ggplot(scan_results, aes(x = factor(Left_Window, levels = Left_Window), y = Asymmetry_Beta, group = 1)) +
  geom_ribbon(aes(ymin = Beta_CI_Lower, ymax = Beta_CI_Upper), fill = "#4682B4", alpha = 0.25) +
  geom_line(color = "#4682B4", linewidth = 1) +
  ylim(-0.5, 0.5) +
  geom_point(aes(color = Asymmetry_P_Value < 0.05), size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_color_manual(values = c("FALSE" = "#4682B4", "TRUE" = "#D95F02"), guide = "none") +
  theme_minimal(base_size = 16) +
  labs(
    title = "B. Sliding Window Effect Sizes Across Binding Site",
    x = "4-bp Window Position (Left Half-Site)",
    y = "Regression Beta (95% CI)"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text = element_text(face = "bold"))

# PANEL C: Model Beta Forest Plot (Comparing Predictor Strengths across Formulations)
model_comparison_df <- data.frame(
  Model = c("Naive OLS\n(N=146)", "Linear Mixed\n(N=146, Operon-Controlled)", "Promoter-Proximal\n(N=92 Unique Loci)"),
  Beta_Baseline = c(coef(lm_global)["Baseline_Magnitude"], 
                    fixef(lmer_global)["Baseline_Magnitude"], 
                    coef(lm_global_unique)["Baseline_Magnitude"]),
  SE_Baseline   = c(summary(lm_global)$coefficients["Baseline_Magnitude", "Std. Error"],
                    summary(lmer_global)$coefficients["Baseline_Magnitude", "Std. Error"],
                    summary(lm_global_unique)$coefficients["Baseline_Magnitude", "Std. Error"]),
  Beta_Asym     = c(coef(lm_global)["Global_Asymmetry"], 
                    fixef(lmer_global)["Global_Asymmetry"], 
                    coef(lm_global_unique)["Global_Asymmetry"]),
  SE_Asym       = c(summary(lm_global)$coefficients["Global_Asymmetry", "Std. Error"],
                    summary(lmer_global)$coefficients["Global_Asymmetry", "Std. Error"],
                    summary(lm_global_unique)$coefficients["Global_Asymmetry", "Std. Error"])
)

forest_df <- pivot_longer(
  model_comparison_df, 
  cols = starts_with("Beta_"), 
  names_to = "Predictor", 
  values_to = "Beta"
) %>%
  mutate(
    SE = if_else(Predictor == "Beta_Baseline", SE_Baseline, SE_Asym),
    Predictor = if_else(Predictor == "Beta_Baseline", "Baseline Regulatory Mag.", "Global Sequence Asymmetry"),
    CI_Lower = Beta - 1.96 * SE,
    CI_Upper = Beta + 1.96 * SE
  )

panel_c <- ggplot(forest_df, aes(x = Beta, y = Model, color = Predictor, shape = Predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 6) +
  xlim(NA, 1) +
  scale_color_manual(values = c("Baseline Regulatory Mag." = "#2E8B57", "Global Sequence Asymmetry" = "#D95F02")) +
  theme_minimal(base_size = 16) +
  labs(
    title = "C. Predictor Consistency Across Model Formulations",
    x = "Effect Size Estimate (Beta ± 95% CI)",
    y = NULL,
    color = "Variable"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12), 
        legend.position = "bottom",
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text = element_text(face = "bold")
  )

# PANEL D: Sequence Logos across Functional Tiers
get_tier_seqs <- function(df, target_tier) {
  seqs <- df %>% filter(Tier == target_tier, !is.na(matched_sequence)) %>% pull(matched_sequence)
  if (length(seqs) < 2) return(NULL)
  return(seqs)
}

# Panel D sequence extraction using collapsed dataset (N = 91)
seqs_rep_sens <- get_tier_seqs(master_df_unique, "Repression_Sensitive")
seqs_rep_mut  <- get_tier_seqs(master_df_unique, "Repression_Muted")
seqs_act_sens <- get_tier_seqs(master_df_unique, "Activation_Sensitive")
seqs_act_mut  <- get_tier_seqs(master_df_unique, "Activation_Muted")

plot_logo_clean <- function(seqs, title) {
  if (is.null(seqs)) return(ggplot() + theme_void())
  ggseqlogo(seqs, method = 'bits') + 
    ggtitle(paste0(title, " (n=", length(seqs), ")")) +
    theme_minimal(base_size = 14) +
    ylim(0.0, 2.0) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(face = "bold")
    )
}

p_logo1 <- plot_logo_clean(seqs_rep_sens, "Repression: Sensitive")
p_logo2 <- plot_logo_clean(seqs_rep_mut,  "Repression: Muted")
p_logo3 <- plot_logo_clean(seqs_act_sens, "Activation: Sensitive")
p_logo4 <- plot_logo_clean(seqs_act_mut,  "Activation: Muted")

panel_d_grid <- plot_grid(p_logo1, p_logo3, p_logo2, p_logo4, ncol = 2, align = "hv")
panel_d <- plot_grid(
  ggdraw() + draw_label("D. Motif Sequence Conservation Across Functional Tiers", fontface = "bold", x = 0, hjust = 0, size = 12),
  panel_d_grid,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# ==============================================================================
# 10. COMPOSITE PUBLICATION FIGURE ASSEMBLY
# ==============================================================================
cat("\nAssembling Composite Publication Figure...\n")

top_row <- plot_grid(panel_a, panel_b, ncol = 2, labels = NULL, rel_widths = c(1, 1))
bottom_row <- plot_grid(panel_c, panel_d, ncol = 2, labels = NULL, rel_widths = c(1, 1))

final_figure <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1.1))

# Save composite figures
ggsave(file.path(CONFIG$out_dir, "Figure_SmcR_PTSP_Asymmetry_Analysis.png"), final_figure, width = 12, height = 10, dpi = 300)
ggsave(file.path(CONFIG$out_dir, "Figure_SmcR_PTSP_Asymmetry_Analysis.pdf"), final_figure, width = 12, height = 10)

ggsave(file.path(CONFIG$out_dir, "Model_Coefficient_Effect_Size.tiff"), panel_c, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(CONFIG$out_dir, "Sliding_window_model.png"), panel_b, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(CONFIG$out_dir, "Site_PWM.png"), panel_d, width = 12, height = 6, dpi = 300)


# ==============================================================================
# 11. SUMMARY TABLE EXTRACTION Across All Statistical Models
# ==============================================================================
cat("\n=== Compiling Model Results Table ===\n")

extract_model_stats <- function(model, model_label) {
  # Dynamically handle lmer (merMod class) vs standard lm objects
  td <- if (inherits(model, "merMod")) {
    broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  } else {
    broom::tidy(model, conf.int = TRUE)
  }
  
  td %>%
    mutate(
      Model = model_label,
      p.value = ifelse(p.value < 0.0001, "< 0.0001", sprintf("%.4f", p.value))
    ) %>%
    select(
      Model, 
      Term = term, 
      `Beta (Estimate)` = estimate, 
      `Std. Error` = std.error, 
      `Statistic` = statistic, 
      `p-value` = p.value, 
      `95% CI Lower` = conf.low, 
      `95% CI Upper` = conf.high
    )
}

# Combine results from all 4 regression models
model_summary_table <- bind_rows(
  extract_model_stats(lm_global, "1A. Naive OLS (N=145)"),
  extract_model_stats(lmer_global, "1B. Linear Mixed Model (N=145)"),
  extract_model_stats(lm_global_unique, "2A. Collapsed OLS - Global Asym (N=91)"),
  extract_model_stats(lm_anchor_unique, "2B. Collapsed OLS - Anchor Asym (N=91)")
)

# Print clean table to console
print(as.data.frame(model_summary_table))

# Save table to disk for manuscript submission
readr::write_csv(model_summary_table, file.path(CONFIG$out_dir, "Table_S1_Regression_Model_Summary.csv"))
cat(sprintf("\nModel summary table successfully saved to '%s/Table_S1_Regression_Model_Summary.csv'.\n", CONFIG$out_dir))
