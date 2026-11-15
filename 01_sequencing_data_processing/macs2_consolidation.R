#!/usr/bin/env Rscript

# ==============================================================================
# Consolidate MACS2 Peak Files into a Multi-Sheet Excel Workbook with Full Stats
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(GenomicRanges)
})

# --- Paths Configuration ---
DIR_MACS <- "03_macs2_peaks"
DIR_OUT  <- "05_final_output"
COMBINED_TXT <- file.path(DIR_OUT, "combined_macs2_peaks.txt")
OUTPUT_XLSX  <- file.path(DIR_OUT, "macs2_consolidated_peaks.xlsx")

wb <- createWorkbook()

# ------------------------------------------------------------------------------
# 1. Load Master Table of All Individual Peak Calls
# ------------------------------------------------------------------------------
message("Loading all individual peak calls...")
if (!file.exists(COMBINED_TXT)) {
  stop(paste("Could not find:", COMBINED_TXT))
}

df_all <- read_tsv(COMBINED_TXT, col_types = cols()) %>%
  mutate(
    summit_genomic_coord = chromStart + peak_summit,
    peak_width = chromEnd - chromStart
  ) %>%
  relocate(peak_width, .after = chromEnd) %>%
  relocate(summit_genomic_coord, .after = peak_summit)

# Create GRanges object for spatial matching
gr_all <- GRanges(
  seqnames = df_all$chrom,
  ranges = IRanges(start = df_all$chromStart, end = df_all$chromEnd),
  mcols = df_all
)

# Sheet 1: Master Combined Table
addWorksheet(wb, sheetName = "All_Peaks_Combined")
writeData(wb, sheet = "All_Peaks_Combined", x = df_all)

# ------------------------------------------------------------------------------
# 2. Master Consensus Peak Boundaries (Annotated with Aggregate Stats)
# ------------------------------------------------------------------------------
consensus_peaks_bed <- file.path(DIR_MACS, "consensus_peaks_merged.bed")
if (file.exists(consensus_peaks_bed)) {
  message("Annotating Master Consensus Regions with peak metrics...")
  
  df_cons_raw <- read_tsv(consensus_peaks_bed, col_names = c("chrom", "start", "end"), col_types = cols())
  gr_cons <- GRanges(seqnames = df_cons_raw$chrom, ranges = IRanges(start = df_cons_raw$start, end = df_cons_raw$end))
  
  # Find overlapping peak calls from individual samples
  hits <- findOverlaps(gr_cons, gr_all)
  
  # Aggregate stats for each consensus region
  agg_stats <- data.frame(
    queryHits = queryHits(hits),
    fold_change = gr_all$mcols.fold_change[subjectHits(hits)],
    log10_pval = gr_all$mcols.log10_pval[subjectHits(hits)],
    log10_qval = gr_all$mcols.log10_qval[subjectHits(hits)],
    Strain = gr_all$mcols.Strain[subjectHits(hits)]
  ) %>%
    group_by(queryHits) %>%
    summarise(
      samples_overlapping = paste(unique(Strain), collapse = "; "),
      n_overlapping_calls = n(),
      mean_fold_change = mean(fold_change, na.rm = TRUE),
      max_fold_change = max(fold_change, na.rm = TRUE),
      mean_log10_qval = mean(log10_qval, na.rm = TRUE),
      max_log10_qval = max(log10_qval, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_consensus <- df_cons_raw %>%
    mutate(
      queryHits = row_number(),
      width = end - start
    ) %>%
    left_join(agg_stats, by = "queryHits") %>%
    select(-queryHits) %>%
    arrange(desc(mean_fold_change))
  
  addWorksheet(wb, sheetName = "Master_Consensus_Regions")
  writeData(wb, sheet = "Master_Consensus_Regions", x = df_consensus)
}

# ------------------------------------------------------------------------------
# 3. Master Reproducible Summits (Annotated with Peak Metrics)
# ------------------------------------------------------------------------------
consensus_summits_bed <- file.path(DIR_MACS, "consensus_summits_merged.bed")
if (file.exists(consensus_summits_bed)) {
  message("Annotating Master Reproducible Summits...")
  
  df_sum_raw <- read_tsv(consensus_summits_bed, col_names = c("chrom", "start", "end"), col_types = cols())
  gr_summits <- GRanges(seqnames = df_sum_raw$chrom, ranges = IRanges(start = df_sum_raw$start, end = df_sum_raw$end))
  
  hits_sum <- findOverlaps(gr_summits, gr_all)
  
  agg_sum_stats <- data.frame(
    queryHits = queryHits(hits_sum),
    fold_change = gr_all$mcols.fold_change[subjectHits(hits_sum)],
    log10_pval = gr_all$mcols.log10_pval[subjectHits(hits_sum)],
    log10_qval = gr_all$mcols.log10_qval[subjectHits(hits_sum)],
    Strain = gr_all$mcols.Strain[subjectHits(hits_sum)]
  ) %>%
    group_by(queryHits) %>%
    summarise(
      samples_overlapping = paste(unique(Strain), collapse = "; "),
      mean_fold_change = mean(fold_change, na.rm = TRUE),
      max_fold_change = max(fold_change, na.rm = TRUE),
      mean_log10_qval = mean(log10_qval, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_summits <- df_sum_raw %>%
    mutate(
      queryHits = row_number(),
      summit_center = floor((start + end) / 2),
      summit_cluster_width = end - start
    ) %>%
    left_join(agg_sum_stats, by = "queryHits") %>%
    select(-queryHits) %>%
    arrange(desc(mean_fold_change))
  
  addWorksheet(wb, sheetName = "Master_Reproducible_Summits")
  writeData(wb, sheet = "Master_Reproducible_Summits", x = df_summits)
}

# ------------------------------------------------------------------------------
# 4. Per-Condition Reproducible Summits (Annotated by Replicate)
# ------------------------------------------------------------------------------
reproducible_files <- list.files(DIR_MACS, pattern = "_reproducible_summits\\.bed$", full.names = TRUE)

for (bed_file in reproducible_files) {
  cond_name <- basename(bed_file) %>%
    str_remove("_reproducible_summits\\.bed")
  
  message(paste("Annotating condition:", cond_name))
  
  df_cond_raw <- read_tsv(bed_file, col_names = c("chrom", "start", "end"), col_types = cols())
  
  if (nrow(df_cond_raw) == 0) next
  
  gr_cond <- GRanges(seqnames = df_cond_raw$chrom, ranges = IRanges(start = df_cond_raw$start, end = df_cond_raw$end))
  
  # Match against individual replicate calls for this condition
  gr_cond_samples <- gr_all[grepl(paste0("^", cond_name, "_rep"), gr_all$mcols.Strain)]
  
  hits_cond <- findOverlaps(gr_cond, gr_cond_samples)
  
  agg_cond_stats <- data.frame(
    queryHits = queryHits(hits_cond),
    fold_change = gr_cond_samples$mcols.fold_change[subjectHits(hits_cond)],
    log10_pval = gr_cond_samples$mcols.log10_pval[subjectHits(hits_cond)],
    log10_qval = gr_cond_samples$mcols.log10_qval[subjectHits(hits_cond)],
    Strain = gr_cond_samples$mcols.Strain[subjectHits(hits_cond)]
  ) %>%
    group_by(queryHits) %>%
    summarise(
      replicates_present = paste(unique(Strain), collapse = "; "),
      mean_fold_change = mean(fold_change, na.rm = TRUE),
      mean_log10_qval = mean(log10_qval, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_cond <- df_cond_raw %>%
    mutate(
      queryHits = row_number(),
      summit_center = floor((start + end) / 2),
      width = end - start
    ) %>%
    left_join(agg_cond_stats, by = "queryHits") %>%
    select(-queryHits) %>%
    arrange(desc(mean_fold_change))
  
  sheet_name <- substr(paste0("Rep_", cond_name), 1, 31)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = df_cond)
}

# ------------------------------------------------------------------------------
# Save Excel File
# ------------------------------------------------------------------------------
saveWorkbook(wb, file = OUTPUT_XLSX, overwrite = TRUE)
message(paste0("[SUCCESS] Fully annotated Excel workbook saved at: ", OUTPUT_XLSX))
