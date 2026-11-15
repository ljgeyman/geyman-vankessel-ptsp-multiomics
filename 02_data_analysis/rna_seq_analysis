#!/usr/bin/env Rscript

# ==============================================================================
# Unified RNA-seq Master Pipeline: 
# DESeq2, Deming Potency + Global DBSCAN, Sieve, ORA, & Plotting
# ==============================================================================

cat("Checking and installing dependencies...\n")
cran_pkgs <- c("tidyverse", "openxlsx", "pheatmap", "ashr", "lmtest", 
               "sandwich", "deming", "dbscan", "ggpubr", "ggrepel", "ggbreak", "ggforce")
bioc_pkgs <- c("DESeq2", "clusterProfiler", "GO.db", "KEGGREST")

new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(new_cran)) install.packages(new_cran, repos = "http://cran.us.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(openxlsx)
  library(pheatmap)
  library(ashr)
  library(lmtest)
  library(sandwich)
  library(rlang)
  library(deming)
  library(dbscan)
  library(ggpubr)
  library(ggrepel)
  library(ggbreak)
  library(ggforce)
  library(clusterProfiler)
  library(GO.db)
  library(KEGGREST)
})

options(stringsAsFactors = FALSE)
set.seed(42)

# ==============================================================================
# 1. Master Configuration & Parameters
# ==============================================================================
CONFIG <- list(
  dataset_name      = "vvul_Transcriptomics",
  out_base_dir      = "analysis_results",
  counts_file       = "analysis_results/counts/all_samples_counts.txt",
  metadata_file     = "metadata.txt",
  
  # Reference paths
  gtf_file          = "refs/refseq_vvul.gtf",
  protein_fasta     = "refs/refseq_vvul.faa",
  eggnog_anno_file  = "refs/vvul_emapper_annotations.xlsx",
  vfdb_blast_file   = "refs/vvul_vfdb_results.tsv",
  
  master_excel_name = "Vvul_Analysis_Master.xlsx",
  
  # Experimental Design
  Target_Gene_Name  = "smcR",
  WT_Genotype_Name  = "dluxO",
  cond_treated      = "ptsp",
  cond_untreated    = "dmso",
  genotype_col      = "genotype",
  condition_col     = "treatment",
  
  # Exclusions for plotting / modeling
  target_locus      = "AOT11_RS09875",
  exclude_loci      = c("AOT11_RS09875", "AOT11_RS01725", "VP_RS12200", "VIBHAR_03459"),
  
  # Thresholds
  lfc_cutoff        = 1.0,
  padj_cutoff       = 0.05,
  vf_evalue_max     = 1e-5,
  vf_pident_min     = 50,
  vf_qcov_min       = 50,
  
  # Global DBSCAN Parameters
  dbscan_eps        = 0.5,
  dbscan_minPts     = 5,
  
  # ORA Settings
  run_eggnog        = FALSE, # Set to TRUE if emapper.py is installed in the environment                             
  eggnog_threads    = 8,
  ora_pval_cutoff   = 0.05,
  ora_qval_cutoff   = 0.20,
  ora_min_genes     = 5,                                 
  
  # Strain Aliases
  focal_strain      = "dluxO",  
  resistant_strain  = "r179c",  
  
  font_family       = if (.Platform$OS.type == "windows") "Arial" else "sans",
  plot_dpi          = 600,
  
  files = list(
    pca_pdf             = "PCA_Group_Means.pdf",
    vf_heatmap_pdf      = "Virulence_Factor_Heatmap.pdf",
    potency_all_pdf     = "Potency_Correlations_All.pdf",
    potency_prefix      = "Potency_DBSCAN_",
    ecdf_tiff           = "Incomplete_Penetrance_ECDF.tiff",
    violin_tiff         = "Incomplete_Penetrance_Horizontal_Violin.tiff",
    dir_volcano_tiff    = "dluxO_vs_dsmcR_Volcano_Directional.tiff",
    off_volcano_tiff    = "Vvul_Volcano_off_target_r179c.tiff"
  ),
  
  colors = list(
    potency_8way = c(
      "Discordant"                  = "#FF3D00",
      "Not Significant"             = "gray75",
      "Upregulated - Both"          = "#1E88E5",
      "Downregulated - Both"        = "#D96B6F",
      "Upregulated - PTSP Only"     = "#154375",
      "Downregulated - PTSP Only"   = "#7A1A1E",
      "Upregulated - Mutant Only"   = "#A0C4E2",
      "Downregulated - Mutant Only" = "#B00006"
    ),
    ecdf_groups = c(
      "Background"      = "#B0BEC5",
      "All Transcripts" = "#B0BEC5",
      "Missed Targets"  = "#1E88E5",
      "Shared Regulon"  = "#B00006",
      "Resistant"       = "#154375",
      "Deletion"        = "#7A1A1E"
    ),
    volcano_3way = c(
      "Downregulated"   = "#B00006",
      "Upregulated"     = "#154375",
      "Not Significant" = "gray75"
    )
  )
)

# ==============================================================================
# 2. Helper Functions
# ==============================================================================
clean_ids  <- function(x) gsub("^(gene-|CDS-|transcript-|rna-)", "", x)
safe_sheet <- function(name) stringr::str_trunc(name, 31, "right", "")

get_shrunken_res <- function(dds_obj, num_group, den_group, anno_df = NULL) {
  if (!all(c(num_group, den_group) %in% levels(dds_obj$group))) return(NULL)
  res <- lfcShrink(dds_obj, contrast = c("group", num_group, den_group), type = "ashr", quiet = TRUE)
  df <- as.data.frame(res) %>% rownames_to_column("Gene") %>% mutate(Gene = clean_ids(Gene))
  if (!is.null(anno_df) && nrow(anno_df) > 0) {
    df <- df %>% left_join(anno_df, by = "Gene") %>% relocate(Gene_Name, Product, .after = Gene)
  }
  return(df)
}

classify_8way <- function(df, lfc_del_col, padj_del_col, lfc_ptsp_col, padj_ptsp_col, lfc_cut, padj_cut) {
  relax_lfc <- 0.5
  relax_padj <- 0.15
  
  df %>%
    mutate(
      Reg_del = case_when(
        .data[[lfc_del_col]] >=  lfc_cut & .data[[padj_del_col]] < padj_cut ~ "Upregulated",
        .data[[lfc_del_col]] <= -lfc_cut & .data[[padj_del_col]] < padj_cut ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      Reg_ptsp = case_when(
        .data[[lfc_ptsp_col]] >=  lfc_cut & .data[[padj_ptsp_col]] < padj_cut ~ "Upregulated",
        .data[[lfc_ptsp_col]] <= -lfc_cut & .data[[padj_ptsp_col]] < padj_cut ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      Trend_del = case_when(
        .data[[lfc_del_col]] >=  relax_lfc & .data[[padj_del_col]] < relax_padj ~ "Upregulated",
        .data[[lfc_del_col]] <= -relax_lfc & .data[[padj_del_col]] < relax_padj ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      Trend_ptsp = case_when(
        .data[[lfc_ptsp_col]] >=  relax_lfc & .data[[padj_ptsp_col]] < relax_padj ~ "Upregulated",
        .data[[lfc_ptsp_col]] <= -relax_lfc & .data[[padj_ptsp_col]] < relax_padj ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      meta_regulation = case_when(
        Reg_del == "Upregulated"   & Reg_ptsp == "Upregulated"   ~ "Upregulated - Both",
        Reg_del == "Downregulated" & Reg_ptsp == "Downregulated" ~ "Downregulated - Both",
        Reg_del == "Upregulated"   & Trend_ptsp == "Upregulated"   ~ "Upregulated - Both",
        Reg_del == "Downregulated" & Trend_ptsp == "Downregulated" ~ "Downregulated - Both",
        Trend_del == "Upregulated"   & Reg_ptsp == "Upregulated"   ~ "Upregulated - Both",
        Trend_del == "Downregulated" & Reg_ptsp == "Downregulated" ~ "Downregulated - Both",
        Reg_del == "Not Significant" & Reg_ptsp == "Upregulated"   ~ "Upregulated - PTSP Only",
        Reg_del == "Not Significant" & Reg_ptsp == "Downregulated" ~ "Downregulated - PTSP Only",
        Reg_del == "Upregulated"   & Reg_ptsp == "Not Significant" ~ "Upregulated - Mutant Only",
        Reg_del == "Downregulated" & Reg_ptsp == "Not Significant" ~ "Downregulated - Mutant Only",
        (Reg_del == "Upregulated" & Reg_ptsp == "Downregulated") | (Reg_del == "Downregulated" & Reg_ptsp == "Upregulated") ~ "Discordant",
        TRUE ~ "Not Significant"
      ),
      meta_regulation = factor(meta_regulation, levels = names(CONFIG$colors$potency_8way))
    )
}

# ==============================================================================
# 3. Directory Setup & GTF Annotation
# ==============================================================================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
res_dir   <- file.path(CONFIG$out_base_dir, paste0(CONFIG$dataset_name, "_", timestamp))
plot_dir  <- file.path(res_dir, "Plots")
walk(c(res_dir, plot_dir), dir.create, recursive = TRUE, showWarnings = FALSE)
cat("Output directory created at:", res_dir, "\n")

cat("Parsing GTF for annotations...\n")
protein2locus <- tibble(protein_clean = character(), locus_clean = character())
anno_dict     <- tibble(Gene = character(), Gene_Name = character(), Product = character())

if (file.exists(CONFIG$gtf_file)) {
  gtf_data <- as.data.frame(rtracklayer::import(CONFIG$gtf_file))
  
  anno_dict <- gtf_data %>%
    dplyr::filter(type %in% c("gene", "CDS", "transcript", "mRNA")) %>%
    mutate(
      Raw_ID     = if ("locus_tag" %in% names(.)) locus_tag else gene_id,
      Clean_Gene = clean_ids(as.character(Raw_ID)),
      Gene_Name  = if ("gene" %in% names(.)) as.character(gene) else NA_character_,
      Product    = if ("product" %in% names(.)) as.character(product) else NA_character_
    ) %>%
    dplyr::filter(!is.na(Clean_Gene)) %>%
    group_by(Clean_Gene) %>%
    summarise(Gene_Name = paste(unique(na.omit(Gene_Name)), collapse = " | "), Product = paste(unique(na.omit(Product)), collapse = " | "), .groups = "drop") %>%
    dplyr::rename(Gene = Clean_Gene)
  
  if (all(c("protein_id", "locus_tag") %in% names(gtf_data))) {
    protein2locus <- gtf_data %>%
      dplyr::filter(type == "CDS", !is.na(protein_id), !is.na(locus_tag)) %>%
      transmute(
        protein_clean = clean_ids(as.character(protein_id)),
        locus_clean   = clean_ids(as.character(locus_tag))
      ) %>%
      distinct(protein_clean, .keep_all = TRUE)
  }
} else {
  warning("GTF file not found at ", CONFIG$gtf_file, ". Proceeding without GTF annotations.")
}

# ==============================================================================
# 4. Data Import, DESeq2 Execution, VST, PCA, TPM, & Z-Scores
# ==============================================================================
cat("Loading counts and metadata...\n")
raw_counts <- read.table(CONFIG$counts_file, header = TRUE, comment.char = "#")
count_data <- raw_counts[, 7:ncol(raw_counts)]
rownames(count_data) <- raw_counts$Geneid
colnames(count_data) <- colnames(count_data) %>% str_remove("^.*\\.bam\\.") %>% str_remove("_sorted.*$")

metadata <- read.table(CONFIG$metadata_file, header = FALSE, row.names = 1, sep = "\t")
missing_in_meta <- setdiff(colnames(count_data), rownames(metadata))
if (length(missing_in_meta) > 0) stop("Counts columns missing from metadata: ", paste(missing_in_meta, collapse = ", "))

metadata <- metadata[colnames(count_data), , drop = FALSE]
metadata[[CONFIG$genotype_col]]  <- sapply(strsplit(rownames(metadata), "_"), `[`, 1)
metadata[[CONFIG$condition_col]] <- sapply(strsplit(rownames(metadata), "_"), `[`, 2)
metadata$group <- factor(paste(metadata[[CONFIG$genotype_col]], metadata[[CONFIG$condition_col]], sep = "_"))

cat("Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_data), colData = metadata, design = ~ group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

cat("Calculating VST, PCA, TPM, and Z-Scores...\n")
vsd     <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

vst_mean_mat <- do.call(cbind, lapply(levels(dds$group), function(grp) {
  rowMeans(vst_mat[, rownames(metadata)[metadata$group == grp], drop = FALSE])
}))
colnames(vst_mean_mat) <- levels(dds$group)
rownames(vst_mean_mat) <- clean_ids(rownames(vst_mat))

master_vst_df <- as.data.frame(vst_mean_mat) %>% rownames_to_column("Gene")
vst_mean_mat  <- master_vst_df %>% column_to_rownames("Gene") %>% as.matrix()

gene_vars <- apply(vst_mean_mat, 1, var)
pca_res   <- prcomp(t(vst_mean_mat[order(gene_vars, decreasing = TRUE)[1:min(500, length(gene_vars))], ]), center = TRUE, scale. = FALSE)
pca_var   <- round((pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100, 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("Group") %>%
  mutate(Strain = sapply(strsplit(Group, "_"), `[`, 1), Treatment = sapply(strsplit(Group, "_"), `[`, 2)) %>%
  relocate(Strain, Treatment, .after = Group)

pdf(file.path(plot_dir, CONFIG$files$pca_pdf), width = 7, height = 6)
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Strain, shape = Treatment)) +
    geom_point(size = 5, alpha = 0.8) +
    labs(title = "PCA of Strain x Treatment (Group Means)", subtitle = "Top 500 most variable genes",
         x = paste0("PC1: ", pca_var[1], "% variance"), y = paste0("PC2: ", pca_var[2], "% variance")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
)
dev.off()

rpk    <- count_data / (raw_counts$Length / 1000)
tpm_df <- as.data.frame(t(t(rpk) / (colSums(rpk) / 1e6))) %>% rownames_to_column("Gene") %>% mutate(Gene = clean_ids(Gene))

zscore_mat <- t(apply(vst_mean_mat, 1, function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)))
master_zscore_df <- as.data.frame(zscore_mat) %>% rownames_to_column("Gene")

# ==============================================================================
# 5. Extract Dynamic Contrasts
# ==============================================================================
cat("Extracting contrasts for all strains...\n")
genos       <- unique(as.character(metadata[[CONFIG$genotype_col]]))
del_name    <- paste0("d", CONFIG$Target_Gene_Name)
all_results <- list()
drug_lfcs   <- list()

mutant_sig_res <- get_shrunken_res(dds, paste0(del_name, "_", CONFIG$cond_untreated), paste0(CONFIG$WT_Genotype_Name, "_", CONFIG$cond_untreated), anno_dict)
if (!is.null(mutant_sig_res)) all_results[[paste0(del_name, "_vs_WT_", CONFIG$cond_untreated)]] <- mutant_sig_res

for (g in genos) {
  res_df <- get_shrunken_res(dds, paste0(g, "_", CONFIG$cond_treated), paste0(g, "_", CONFIG$cond_untreated), anno_dict)
  if (!is.null(res_df)) {
    all_results[[paste0(g, "_Drug_Response")]] <- res_df
    drug_lfcs[[g]] <- res_df %>% dplyr::select(Gene, !!g := log2FoldChange)
  }
}

for (g in setdiff(genos, del_name)) {
  for (cond in c(CONFIG$cond_untreated, CONFIG$cond_treated)) {
    res_del <- get_shrunken_res(dds, paste0(g, "_", cond), paste0(del_name, "_", cond), anno_dict)
    if (!is.null(res_del)) all_results[[paste0(g, "_vs_", del_name, "_", cond)]] <- res_del
  }
  if (g != CONFIG$WT_Genotype_Name) {
    res_wt <- get_shrunken_res(dds, paste0(g, "_", CONFIG$cond_untreated), paste0(CONFIG$WT_Genotype_Name, "_", CONFIG$cond_untreated), anno_dict)
    if (!is.null(res_wt)) all_results[[paste0(g, "_vs_WT_", CONFIG$cond_untreated)]] <- res_wt
  }
}

direct_res <- get_shrunken_res(dds, paste0(CONFIG$WT_Genotype_Name, "_", CONFIG$cond_treated), paste0(del_name, "_", CONFIG$cond_untreated), anno_dict)
if (!is.null(direct_res)) all_results[[paste0(CONFIG$WT_Genotype_Name, "_", CONFIG$cond_treated, "_vs_", del_name, "_", CONFIG$cond_untreated)]] <- direct_res

# ==============================================================================
# 6. Potency Analysis, Deming Regression, & Global DBSCAN
# ==============================================================================
cat("Running Potency Analysis (Deming Regression + Global DBSCAN)...\n")
sig_genes_mutant <- mutant_sig_res %>% dplyr::filter(padj < CONFIG$padj_cutoff, abs(log2FoldChange) > CONFIG$lfc_cutoff) %>% pull(Gene)

norm_counts <- counts(dds, normalized = TRUE)
wt_samples  <- rownames(metadata[metadata[[CONFIG$genotype_col]] == CONFIG$WT_Genotype_Name & metadata[[CONFIG$condition_col]] == CONFIG$cond_untreated, ])
del_samples <- rownames(metadata[metadata[[CONFIG$genotype_col]] == del_name & metadata[[CONFIG$condition_col]] == CONFIG$cond_untreated, ])

baseline_rho <- NA_real_
if (length(wt_samples) >= 1 && length(del_samples) >= 2) {
  wt_mean        <- rowMeans(log2(norm_counts[, wt_samples, drop = FALSE] + 1))
  ceiling_matrix <- log2(norm_counts[, del_samples, drop = FALSE] + 1) - wt_mean
  ceiling_df     <- as.data.frame(ceiling_matrix) %>% rownames_to_column("Gene") %>%
    dplyr::filter(Gene %in% sig_genes_mutant, !(Gene %in% CONFIG$exclude_loci)) %>% drop_na()
  cor_mat      <- cor(ceiling_df %>% column_to_rownames("Gene"), method = "spearman")
  baseline_rho <- mean(cor_mat[upper.tri(cor_mat)])
}

potency_summary     <- tibble()
potency_data_sheets <- list()
potency_plots_list  <- list()

for (g in setdiff(genos, del_name)) {
  drug_key <- paste0(g, "_Drug_Response")
  if (!(g %in% names(drug_lfcs) && drug_key %in% names(all_results))) next
  current_del_res <- get_shrunken_res(dds, paste0(del_name, "_", CONFIG$cond_untreated), paste0(g, "_", CONFIG$cond_untreated), anno_dict)
  if (is.null(current_del_res)) next
  
  contrast_name <- paste0(del_name, "_vs_", g, "_", CONFIG$cond_untreated)
  drug_res_df <- all_results[[drug_key]]
  
  sig_pool <- union(
    current_del_res %>% dplyr::filter(padj < CONFIG$padj_cutoff, abs(log2FoldChange) > CONFIG$lfc_cutoff) %>% pull(Gene),
    drug_res_df     %>% dplyr::filter(padj < CONFIG$padj_cutoff, abs(log2FoldChange) > CONFIG$lfc_cutoff) %>% pull(Gene)
  )
  
  pot_df <- current_del_res %>%
    dplyr::select(Gene, any_of(c("Gene_Name", "Product")), baseMean, Del_LFC = log2FoldChange, Del_lfcSE = lfcSE, Del_padj = padj) %>%
    inner_join(drug_res_df %>% dplyr::select(Gene, Drug_LFC = log2FoldChange, Drug_lfcSE = lfcSE, Drug_padj = padj), by = "Gene") %>%
    dplyr::filter(Gene %in% sig_pool, !(Gene %in% CONFIG$exclude_loci)) %>%
    drop_na(Del_LFC, Drug_LFC, Del_lfcSE, Drug_lfcSE) %>%
    classify_8way("Del_LFC", "Del_padj", "Drug_LFC", "Drug_padj", CONFIG$lfc_cutoff, CONFIG$padj_cutoff)
  
  if (nrow(pot_df) < 5) next
  
  rho        <- cor(pot_df$Drug_LFC, pot_df$Del_LFC, method = "spearman")
  pot_deming <- deming(Drug_LFC ~ Del_LFC, data = pot_df, xstd = Del_lfcSE, ystd = Drug_lfcSE)
  
  d_int      <- unname(pot_deming$coefficients[1])
  d_slope    <- unname(pot_deming$coefficients[2])
  d_slope_se <- unname(pot_deming$se[2])
  pval_0     <- 2 * (1 - pnorm(abs(d_slope / d_slope_se)))
  pval_ident <- 2 * (1 - pnorm(abs((d_slope - 1) / d_slope_se)))
  
  # Global Unsegregated DBSCAN
  pot_df$DBSCAN_Cluster <- "Noise"
  
  if (nrow(pot_df) >= CONFIG$dbscan_minPts) {
    coords <- as.matrix(pot_df[, c("Del_LFC", "Drug_LFC")])
    db_res <- dbscan::dbscan(coords, eps = CONFIG$dbscan_eps, minPts = CONFIG$dbscan_minPts)
    has_clusters <- db_res$cluster > 0
    if (any(has_clusters)) {
      pot_df$DBSCAN_Cluster[has_clusters] <- paste0("Cluster_", db_res$cluster[has_clusters])
    }
    n_clusters <- max(db_res$cluster)
  } else {
    n_clusters <- 0
  }
  
  pot_df <- pot_df %>% arrange(DBSCAN_Cluster, meta_regulation)
  potency_data_sheets[[paste0("Potency_", g)]] <- pot_df
  
  n_noise  <- sum(pot_df$DBSCAN_Cluster == "Noise")
  axis_lim <- ceiling(max(abs(c(pot_df$Del_LFC, pot_df$Drug_LFC)), na.rm = TRUE))
  
  hull_df <- pot_df %>%
    dplyr::filter(DBSCAN_Cluster != "Noise") %>%
    group_by(DBSCAN_Cluster) %>%
    dplyr::filter(n() >= CONFIG$dbscan_minPts) %>%
    dplyr::slice(chull(Del_LFC, Drug_LFC)) %>%
    ungroup()
  
  p_pot <- ggplot(pot_df, aes(x = Del_LFC, y = Drug_LFC)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
    annotate("segment", x = -5, xend = 5, y = -5 + d_int, yend = 5 + d_int, color = "gray40", linetype = "dotted", linewidth = 1) +
    annotate("segment", x = -5, xend = 5, y = (d_slope * -5) + d_int, yend = (d_slope * 5) + d_int, color = "#B00006", linewidth = 1) +
    {if (nrow(hull_df) > 0) geom_mark_hull(
      data = dplyr::filter(pot_df, DBSCAN_Cluster %in% hull_df$DBSCAN_Cluster),
      aes(group = DBSCAN_Cluster, label = DBSCAN_Cluster),
      concavity = 2.5, expand = unit(2.5, "mm"), radius = unit(1.5, "mm"),
      color = "black", fill = "gray80", alpha = 0.2, linewidth = 0.8
    )} +
    geom_point(aes(color = meta_regulation, shape = DBSCAN_Cluster == "Noise"), size = 2.2, alpha = 0.85) +
    scale_color_manual(values = CONFIG$colors$potency_8way, name = "Regulation Status", drop = FALSE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "Clustered", `TRUE` = "DBSCAN Noise"), name = "DBSCAN") +
    coord_cartesian(xlim = c(-axis_lim, axis_lim), ylim = c(-axis_lim, axis_lim)) +
    labs(
      title    = paste0("Potency & Global DBSCAN Analysis: ", g),
      subtitle = sprintf("Spearman rho = %.3f | Deming Slope = %.2f | Global DBSCAN: %d clusters (noise = %d)",
                         rho, d_slope, n_clusters, n_noise),
      x        = paste0("LFC: ", del_name, " vs ", g, " (", CONFIG$cond_untreated, ")"),
      y        = paste0("LFC: ", g, " (", CONFIG$cond_treated, " vs ", CONFIG$cond_untreated, ")")
    ) +
    theme_classic(base_size = 13, base_family = CONFIG$font_family) +
    theme(
      axis.line        = element_line(linewidth = 1.2, color = "black"),
      axis.ticks       = element_line(linewidth = 1.2, color = "black"),
      axis.text        = element_text(face = "bold", color = "black"),
      axis.title       = element_text(face = "bold", size = 14),
      plot.title       = element_text(face = "bold", size = 15),
      plot.subtitle    = element_text(size = 10, color = "gray20"),
      legend.text      = element_text(size = 9, face = "bold"),
      legend.title     = element_text(size = 10, face = "bold"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  ggsave(file.path(plot_dir, paste0(CONFIG$files$potency_prefix, g, ".tiff")),
         plot = p_pot, bg = "transparent", width = 8.5, height = 6.5, dpi = CONFIG$plot_dpi)
  
  potency_plots_list[[g]] <- p_pot
  
  potency_summary <- bind_rows(potency_summary, tibble(
    Strain               = g,
    Matched_Contrast     = contrast_name,
    Spearman_Rho         = rho,
    Baseline_Rho         = baseline_rho,
    Deming_Slope         = d_slope,
    Deming_SE            = d_slope_se,
    Deming_P_Value_0     = pval_0,
    Deming_P_Value_Ident = pval_ident,
    Deming_Intercept     = d_int,
    DBSCAN_Clusters      = n_clusters,
    DBSCAN_Noise_Genes   = n_noise
  ))
}

if (length(potency_plots_list) > 0) {
  potency_pdf_path <- file.path(plot_dir, CONFIG$files$potency_all_pdf)
  pdf(potency_pdf_path, width = 8.5, height = 6.5, family = "Helvetica", useDingbats = FALSE)
  walk(potency_plots_list, print)
  dev.off()
}

# ==============================================================================
# 7. Transcriptomic "Sieve" Testing & Classification
# ==============================================================================
cat("Performing Classification and Transcriptomic Sieve testing...\n")
total_genes       <- nrow(dds)
sieve_summary     <- tibble()
sieve_data_sheets <- list()

for (g in setdiff(genos, del_name)) {
  if (!(g %in% names(drug_lfcs))) next
  
  mutant_vs_target_res <- get_shrunken_res(dds, paste0(g, "_", CONFIG$cond_untreated), paste0(del_name, "_", CONFIG$cond_untreated), anno_dict)
  if (is.null(mutant_vs_target_res)) next
  
  mutant_degs <- mutant_vs_target_res %>% dplyr::filter(padj < CONFIG$padj_cutoff, abs(log2FoldChange) > CONFIG$lfc_cutoff) %>% pull(Gene)
  drug_res    <- all_results[[paste0(g, "_Drug_Response")]]
  qsi_degs    <- drug_res %>% dplyr::filter(padj < CONFIG$padj_cutoff, abs(log2FoldChange) > CONFIG$lfc_cutoff) %>% pull(Gene)
  
  overlap_genes <- intersect(mutant_degs, qsi_degs)
  K_mut <- length(mutant_degs); n_qsi <- length(qsi_degs); k_ov <- length(overlap_genes)
  
  pval <- phyper(k_ov - 1, K_mut, total_genes - K_mut, n_qsi, lower.tail = FALSE)
  off_target_ratio <- ifelse(n_qsi > 0, (n_qsi - k_ov) / n_qsi, NA_real_)
  
  sieve_summary <- bind_rows(sieve_summary, tibble(
    Strain = g, Mutant_DEGs = K_mut, QSI_DEGs = n_qsi, Overlap = k_ov,
    Hypergeometric_P = pval, Off_Target_Ratio = off_target_ratio
  ))
  
  sieve_data_sheets[[paste0("Sieve_", g)]] <- drug_res %>%
    dplyr::filter(Gene %in% union(qsi_degs, mutant_degs)) %>%
    mutate(Classification = case_when(
      Gene %in% overlap_genes ~ "True Positive (Shared Regulon)",
      Gene %in% qsi_degs      ~ "False Positive (Off-Target QSI)",
      Gene %in% mutant_degs   ~ "False Negative (Missed Target)",
      TRUE                    ~ "Other"
    )) %>%
    dplyr::select(Gene, any_of(c("Gene_Name", "Product")), baseMean, log2FoldChange, lfcSE, pvalue, padj, Classification) %>%
    arrange(Classification, padj)
}

# ==============================================================================
# 8. Virulence Factor (VFDB) Heatmap Processing
# ==============================================================================
cat("Processing Virulence Factors...\n")
vf_heatmap_data <- tibble()
vf_zscore_data  <- tibble()

if (file.exists(CONFIG$vfdb_blast_file)) {
  vf_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "stitle", "vf_category")
  vf_raw <- tryCatch(read.delim(CONFIG$vfdb_blast_file, header = TRUE, col.names = vf_cols), error = function(e) NULL)
  
  if (!is.null(vf_raw) && "qlen" %in% colnames(vf_raw)) {
    vf_filtered <- vf_raw %>%
      mutate(qcov = (length / qlen) * 100, protein_clean = clean_ids(qseqid)) %>%
      dplyr::filter(evalue <= CONFIG$vf_evalue_max, pident >= CONFIG$vf_pident_min, qcov >= CONFIG$vf_qcov_min) %>%
      group_by(protein_clean) %>%
      dplyr::slice_max(bitscore, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(VF_Name = str_extract(stitle, "(?<=\\[)[^\\]]+(?=\\])") %>% coalesce(sseqid))
    
    vf_translated <- if (nrow(protein2locus) > 0) {
      vf_filtered %>% inner_join(protein2locus, by = "protein_clean") %>% dplyr::rename(Gene = locus_clean)
    } else {
      vf_filtered %>% dplyr::rename(Gene = protein_clean)
    }
    
    focal_drug_key <- paste0(CONFIG$focal_strain, "_Drug_Response")
    sig_focal_drug <- if (focal_drug_key %in% names(all_results)) {
      all_results[[focal_drug_key]] %>% dplyr::filter(padj < CONFIG$padj_cutoff) %>% pull(Gene)
    } else character(0)
    
    vf_translated   <- vf_translated %>% dplyr::filter(Gene %in% union(sig_genes_mutant, sig_focal_drug))
    vf_heatmap_data <- master_vst_df %>% inner_join(vf_translated %>% dplyr::select(Gene, VF_Name), by = "Gene") %>% drop_na() %>% distinct(Gene, .keep_all = TRUE)
    vf_zscore_data  <- master_zscore_df %>% inner_join(vf_translated %>% dplyr::select(Gene, VF_Name), by = "Gene") %>% drop_na() %>% distinct(Gene, .keep_all = TRUE)
    
    if (nrow(vf_heatmap_data) > 1) {
      mat <- vf_heatmap_data %>% dplyr::select(-Gene, -VF_Name) %>% as.matrix()
      rownames(mat) <- paste0(vf_heatmap_data$Gene, " (", vf_heatmap_data$VF_Name, ")")
      pdf(file.path(plot_dir, CONFIG$files$vf_heatmap_pdf), width = 8, height = max(6, nrow(mat) * 0.2))
      pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE, scale = "row", display_numbers = FALSE, main = "Virulence Factor Expression (Row-Scaled VST)")
      dev.off()
    }
  }
}

# ==============================================================================
# 9. Incomplete Penetrance Evaluation: ECDF & Horizontal Violin Plots
# ==============================================================================
sieve_focal_key <- paste0("Sieve_", CONFIG$focal_strain)
drug_focal_key  <- paste0(CONFIG$focal_strain, "_Drug_Response")
resist_drug_key <- paste0(CONFIG$resistant_strain, "_Drug_Response")
del_focal_key   <- paste0(del_name, "_vs_WT_", CONFIG$cond_untreated) # e.g. dsmcR_vs_WT_dmso

if (all(c(sieve_focal_key) %in% names(sieve_data_sheets)) && all(c(drug_focal_key, resist_drug_key, del_focal_key) %in% names(all_results))) {
  cat("Running Incomplete Penetrance ECDF & Rank-Sum Analysis...\n")
  sieve_data   <- sieve_data_sheets[[sieve_focal_key]]
  full_res     <- all_results[[drug_focal_key]]
  resistant_mt <- all_results[[resist_drug_key]]
  del_res      <- all_results[[del_focal_key]]
  
  background_data <- full_res %>%
    dplyr::filter(!(Gene %in% sieve_data$Gene)) %>%
    mutate(Classification = "Background") %>%
    dplyr::select(Gene, baseMean, log2FoldChange, lfcSE, pvalue, padj, Classification)
  
  r_mt_data <- resistant_mt %>%
    dplyr::filter(Gene %in% sieve_data$Gene) %>%
    mutate(Classification = "Resistant") %>%
    dplyr::select(Gene, baseMean, log2FoldChange, lfcSE, pvalue, padj, Classification)
    
  del_data <- del_res %>% 
    dplyr::filter(Gene %in% sieve_data$Gene) %>% 
    dplyr::filter(Gene != CONFIG$target_locus) %>%  # Exclude self (AOT11_RS09875)
    mutate(Classification = 'Deletion') %>% 
    dplyr::select(Gene, baseMean, log2FoldChange, lfcSE, pvalue, padj, Classification)
  
  eval_df <- bind_rows(sieve_data, background_data, r_mt_data, del_data) %>%
    dplyr::filter(Classification != "False Positive (Off-Target QSI)") %>%
    mutate(
      Group = case_when(
        Classification == "True Positive (Shared Regulon)" ~ "Shared Regulon",
        Classification == "False Negative (Missed Target)" ~ "Missed Targets",
        TRUE ~ Classification
      ),
      Group = factor(Group, levels = c("Background", "Missed Targets", "Shared Regulon", "Resistant", "Deletion")),
      Abs_Z = abs(log2FoldChange / lfcSE)
    ) %>%
    drop_na(Abs_Z)
  
  kw_res <- kruskal.test(Abs_Z ~ Group, data = eval_df)
  cat(sprintf("Kruskal-Wallis Global p-value: %g\n", kw_res$p.value))
  
  # --- ECDF Plot ---
  p_ecdf <- ggplot(eval_df, aes(x = Abs_Z, color = Group)) +
    stat_ecdf(geom = "step", linewidth = 1.2) +
    scale_color_manual(values = CONFIG$colors$ecdf_groups) +
    coord_cartesian(xlim = c(0, 7)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
    labs(x = "Absolute Pseudo-Wald Statistic (|Z|)", y = "Cumulative Proportion") +
    theme_classic(base_size = 14, base_family = CONFIG$font_family) +
    theme(
      legend.justification  = c(0.75, 0.25),
      legend.position       = c(0.95, 0.10),
      legend.title          = element_blank(),
      legend.text           = element_text(face = "bold", color = "black"),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.background     = element_rect(fill = "transparent", color = NA),
      axis.line             = element_line(linewidth = 1.5, color = "black"),
      axis.ticks            = element_line(linewidth = 1.5, color = "black"),
      axis.ticks.length     = unit(0.2, "cm"),
      axis.text             = element_text(face = "bold", color = "black"),
      axis.title.x          = element_text(size = 18, margin = margin(t = 15)),
      axis.title.y          = element_text(size = 18, margin = margin(r = 15)),
      plot.background       = element_rect(fill = "transparent", color = NA),
      panel.background      = element_rect(fill = "transparent", color = NA)
    )
  ggsave(file.path(plot_dir, CONFIG$files$ecdf_tiff), plot = p_ecdf, bg = "transparent", dpi = 400, width = 7, height = 5)
  
  # --- Horizontal Violin Plot ---
  plot_df_violin <- eval_df %>%
    mutate(
      Plot_Group = ifelse(as.character(Group) == "Background", "All Transcripts", as.character(Group)),
      Plot_Group = factor(Plot_Group, levels = c("Shared Regulon", "Missed Targets", "All Transcripts", "Resistant", "Deletion"))
    )
  
  true_max <- max(plot_df_violin$Abs_Z, na.rm = TRUE)
  step     <- max(0.5, true_max * 0.08)
  my_comparisons <- list(c("All Transcripts", "Missed Targets"), c("Missed Targets", "Shared Regulon"), c("All Transcripts", "Shared Regulon"))
  
  p_horizontal <- ggplot(plot_df_violin, aes(x = Plot_Group, y = Abs_Z, fill = Plot_Group)) +
    geom_violin(scale = "width", trim = FALSE, linewidth = 1.2, color = "black", draw_quantiles = c(0.5)) +
    scale_fill_manual(values = scales::alpha(CONFIG$colors$ecdf_groups, 0.8)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", p.adjust.method = "holm",
                       label = "p.signif", bracket.size = 1, size = 8,
                       label.y = c(true_max - (step * 3.5), true_max + (step * 3.0), true_max + (step * 4.5))) +
    coord_flip(ylim = c(-0.1, true_max + (step * 6.5))) +
    labs(x = NULL, y = "Absolute Pseudo-Wald Statistic (|Z|)") +
    theme_classic(base_size = 16, base_family = CONFIG$font_family) +
    theme(
      legend.position   = "none",
      axis.line         = element_line(linewidth = 1.5, color = "black"),
      axis.ticks        = element_line(linewidth = 1.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(face = "bold", color = "black"),
      axis.text.y       = element_text(size = 16),
      axis.title.x      = element_text(size = 18, margin = margin(t = 15)),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )
  ggsave(file.path(plot_dir, CONFIG$files$violin_tiff), plot = p_horizontal, bg = "transparent", width = 7, height = 5, dpi = CONFIG$plot_dpi)
}

# ==============================================================================
# 10. Volcano Plots (Directional Overlay & Off-Target)
# ==============================================================================
target_vol_key <- paste0(CONFIG$focal_strain, "_vs_", del_name, "_", CONFIG$cond_untreated)
if (all(c(target_vol_key, drug_focal_key) %in% names(all_results))) {
  cat("Generating Directional Volcano Plot...\n")
  res_df  <- all_results[[target_vol_key]] %>% drop_na(log2FoldChange, padj) %>% dplyr::filter(!(Gene %in% CONFIG$exclude_loci))
  ptsp_df <- all_results[[drug_focal_key]] %>% drop_na(log2FoldChange, padj) %>% dplyr::filter(!(Gene %in% CONFIG$exclude_loci)) %>%
    mutate(log2FoldChange = -log2FoldChange)
  
  vol_df <- res_df %>%
    inner_join(ptsp_df, by = "Gene", suffix = c("_del", "_ptsp")) %>%
    mutate(
      padj_plot_del      = ifelse(padj_del == 0, .Machine$double.xmin, padj_del),
      neg_log10_padj_del = -log10(padj_plot_del)
    ) %>%
    classify_8way("log2FoldChange_del", "padj_del", "log2FoldChange_ptsp", "padj_ptsp", CONFIG$lfc_cutoff, CONFIG$padj_cutoff) %>%
    arrange(meta_regulation)
  
  p_volcano <- ggplot(vol_df, aes(x = log2FoldChange_del, y = neg_log10_padj_del)) +
    geom_vline(xintercept = c(-CONFIG$lfc_cutoff, CONFIG$lfc_cutoff), linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_hline(yintercept = -log10(CONFIG$padj_cutoff), linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_point(aes(color = meta_regulation), alpha = 0.75, size = 1) +
    scale_color_manual(values = CONFIG$colors$potency_8way) +
    labs(x = expression("LFC: " * Delta * "luxO vs " * Delta * "smcR"), y = expression(-log[10] ~ "(adj-" ~ italic(p) * ")")) +
    scale_x_continuous(breaks = seq(-12, 12, by = 2), limits = c(-12, 12)) +
    theme_classic(base_size = 14, base_family = CONFIG$font_family) +
    theme(
      axis.line         = element_line(linewidth = 1.5, color = "black"),
      axis.ticks        = element_line(linewidth = 1.5, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(face = "bold", color = "black"),
      axis.title.x      = element_text(size = 18, margin = margin(t = 15)),
      axis.title.y      = element_text(size = 18, margin = margin(r = 15)),
      legend.title      = element_blank(),
      legend.text       = element_text(family = CONFIG$font_family, face = "bold", size = 11),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.position   = "inside",
      legend.justification = c(0.80, 0.20),
      legend.position.inside = c(0.85, 0.35),
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA)
    )
    
  # Apply y-axis break
  p_volcano_broken <- p_volcano + 
    scale_y_break(c(40.0, 60.0), scales = 0.1) +
    theme(
      axis.text.x.top  = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.line.x.top  = element_blank()
    )
  
  ggsave(file.path(plot_dir, CONFIG$files$dir_volcano_tiff), plot = p_volcano_broken, bg = "transparent", width = 7, height = 6, dpi = 300)
}

if (resist_drug_key %in% names(all_results)) {
  cat("Generating Off-Target Volcano Plot for:", resist_drug_key, "...\n")
  off_df <- all_results[[resist_drug_key]] %>%
    drop_na(log2FoldChange, padj) %>%
    dplyr::filter(!(Gene %in% CONFIG$exclude_loci)) %>%
    mutate(
      padj_plot      = ifelse(padj == 0, .Machine$double.xmin, padj),
      neg_log10_padj = -log10(padj_plot),
      Regulation     = case_when(
        log2FoldChange >=  CONFIG$lfc_cutoff & padj < CONFIG$padj_cutoff ~ "Upregulated",
        log2FoldChange <= -CONFIG$lfc_cutoff & padj < CONFIG$padj_cutoff ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      Regulation = factor(Regulation, levels = c("Downregulated", "Not Significant", "Upregulated")),
      Label      = if ("Gene_Name" %in% names(.) && !all(is.na(Gene_Name))) ifelse(is.na(Gene_Name) | Gene_Name == "", Gene, Gene_Name) else Gene
    )
  
  top_off_genes <- off_df %>%
    dplyr::filter(Regulation != "Not Significant") %>%
    group_by(Regulation) %>%
    dplyr::slice_min(order_by = padj, n = 5, with_ties = FALSE) %>%
    ungroup()
  
  p_off_volcano <- ggplot(off_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_vline(xintercept = c(-CONFIG$lfc_cutoff, CONFIG$lfc_cutoff), linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_hline(yintercept = -log10(CONFIG$padj_cutoff), linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_point(aes(color = Regulation), alpha = 0.75, size = 1) +
    scale_color_manual(values = CONFIG$colors$volcano_3way) +
    geom_text_repel(data = top_off_genes, aes(label = Label), family = CONFIG$font_family, fontface = "bold",
                    size = 4, box.padding = 0.4, point.padding = 0.3, max.overlaps = 20, segment.color = "black") +
    labs(x = expression("LFC: R179C PTSP"), y = expression(bold(-log[10] ~ "(adj-" ~ italic(p) * ")"))) +
    theme_classic(base_size = 14, base_family = CONFIG$font_family) +
    theme(
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(face = "bold", color = "black"),
      axis.line        = element_line(linewidth = 1.2, color = "black"),
      axis.ticks       = element_line(linewidth = 1.2, color = "black"),
      legend.position  = "none",
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  ggsave(file.path(plot_dir, CONFIG$files$off_volcano_tiff), plot = p_off_volcano, bg = "transparent", width = 7, height = 6, dpi = 300)
}

# ==============================================================================
# 11. Translated Over-Representation Analysis (ORA)
# ==============================================================================
cat("Running EggNOG Annotation Translation & Over-Representation Analysis (ORA)...\n")

if (CONFIG$run_eggnog) {
  if (file.exists(CONFIG$protein_fasta)) {
    cat("Executing emapper.py...\n")
    cmd <- sprintf("emapper.py -i %s -o vvul_eggnog --cpu %d", CONFIG$protein_fasta, CONFIG$eggnog_threads)
    system(cmd)
  }
}

ora_results_list <- list()

cog_dict <- c(
  "J" = "Translation, ribosomal structure and biogenesis",
  "A" = "RNA processing and modification",
  "K" = "Transcription",
  "L" = "Replication, recombination and repair",
  "B" = "Chromatin structure and dynamics",
  "D" = "Cell cycle control, cell division, chromosome partitioning",
  "Y" = "Nuclear structure",
  "V" = "Defense mechanisms",
  "T" = "Signal transduction mechanisms",
  "M" = "Cell wall/membrane/envelope biogenesis",
  "N" = "Cell motility",
  "Z" = "Cytoskeleton",
  "W" = "Extracellular structures",
  "U" = "Intracellular trafficking, secretion, and vesicular transport",
  "O" = "Posttranslational modification, protein turnover, chaperones",
  "X" = "Mobilome: transposons, prophages",
  "C" = "Energy production and conversion",
  "G" = "Carbohydrate transport and metabolism",
  "E" = "Amino acid transport and metabolism",
  "F" = "Nucleotide transport and metabolism",
  "H" = "Coenzyme transport and metabolism",
  "I" = "Lipid transport and metabolism",
  "P" = "Inorganic ion transport and metabolism",
  "Q" = "Secondary metabolites biosynthesis, transport and catabolism",
  "R" = "General function prediction only",
  "S" = "Function unknown"
)

if (file.exists(CONFIG$eggnog_anno_file) && nrow(protein2locus) > 0) {
  if (grepl("\\.xlsx?$", CONFIG$eggnog_anno_file, ignore.case = TRUE)) {
    preview_df <- openxlsx::read.xlsx(CONFIG$eggnog_anno_file, sheet = 1, rows = 1:10, colNames = FALSE)
    header_row <- which(apply(preview_df, 1, function(r) any(grepl("^#?query$", as.character(r), ignore.case = TRUE))))[1]
    if (is.na(header_row)) header_row <- 1
    eggnog_raw <- openxlsx::read.xlsx(CONFIG$eggnog_anno_file, sheet = 1, startRow = header_row)
  } else {
    eggnog_lines <- readLines(CONFIG$eggnog_anno_file)
    eggnog_raw   <- read.delim(text = eggnog_lines[!grepl("^##", eggnog_lines)], header = TRUE, stringsAsFactors = FALSE)
  }
  
  names(eggnog_raw) <- gsub("^#", "", trimws(names(eggnog_raw)))
  names(eggnog_raw)[1] <- "query"
  
  eggnog_mapped <- eggnog_raw %>%
    mutate(protein_clean = clean_ids(as.character(query))) %>%
    inner_join(protein2locus, by = "protein_clean") %>%
    dplyr::rename(Gene = locus_clean)
  
  cat(sprintf("  -> Mapped %d of %d EggNOG protein records to genomic locus tags.\n", nrow(eggnog_mapped), nrow(eggnog_raw)))
  
  cat("  -> Downloading Official NCBI COG Mapping Database...\n")
  cog_def_url <- "https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab"
  cog_mapping <- tryCatch({
    df <- read.delim(cog_def_url, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    df %>% dplyr::select(COG_ID = V1, func_category = V2)
  }, error = function(e) NULL)
  
  t2g_list <- list()
  t2n_list <- list()
  
  if ("COG_category" %in% names(eggnog_mapped)) {
    cog_df <- eggnog_mapped %>%
      dplyr::select(Gene, COG_category) %>%
      dplyr::filter(!is.na(COG_category) & COG_category != "" & COG_category != "-") %>%
      tidyr::separate_rows(COG_category, sep = "[,;\\s]+")
    
    if (!is.null(cog_mapping)) {
      cog_df <- cog_df %>%
        left_join(cog_mapping, by = c("COG_category" = "COG_ID")) %>%
        mutate(TERM = ifelse(!is.na(func_category), func_category, COG_category))
    } else {
      cog_df <- cog_df %>% mutate(TERM = COG_category)
    }
    
    t2g_list[["COG"]] <- cog_df %>%
      mutate(TERM_split = stringr::str_split(TERM, "")) %>%
      tidyr::unnest(TERM_split) %>%
      mutate(TERM_final = case_when(
        TERM_split %in% names(cog_dict) ~ paste0("[", TERM_split, "] ", cog_dict[TERM_split]),
        TRUE ~ TERM_split
      )) %>%
      dplyr::select(TERM = TERM_final, GENE = Gene) %>%
      distinct()
    
    t2n_list[["COG"]] <- data.frame(TERM = unique(t2g_list[["COG"]]$TERM), NAME = unique(t2g_list[["COG"]]$TERM))
  }
  
  if ("GOs" %in% names(eggnog_mapped)) {
    t2g_list[["GO"]] <- eggnog_mapped %>%
      dplyr::select(Gene, GOs) %>%
      dplyr::filter(!is.na(GOs) & GOs != "" & GOs != "-") %>%
      tidyr::separate_rows(GOs, sep = ",") %>%
      dplyr::select(TERM = GOs, GENE = Gene) %>%
      distinct()
    
    go_ids <- unique(t2g_list[["GO"]]$TERM)
    go_t2n <- tryCatch({
      res <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = go_ids, columns = "TERM", keytype = "GOID"))
      data.frame(TERM = res$GOID, NAME = coalesce(res$TERM, res$GOID))
    }, error = function(e) data.frame(TERM = go_ids, NAME = go_ids))
    t2n_list[["GO"]] <- go_t2n
  }
  
  if ("KEGG_Pathway" %in% names(eggnog_mapped)) {
    t2g_list[["KEGG_Pathway"]] <- eggnog_mapped %>%
      dplyr::select(Gene, KEGG_Pathway) %>%
      dplyr::filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "" & KEGG_Pathway != "-") %>%
      tidyr::separate_rows(KEGG_Pathway, sep = ",") %>%
      dplyr::filter(grepl("^ko", KEGG_Pathway)) %>%
      dplyr::select(TERM = KEGG_Pathway, GENE = Gene) %>%
      distinct()
    
    pathway_ids <- unique(t2g_list[["KEGG_Pathway"]]$TERM)
    path_t2n <- data.frame(TERM = pathway_ids, NAME = pathway_ids)
    try({
      kegg_paths <- KEGGREST::keggList("pathway")
      names(kegg_paths) <- gsub("path:map", "ko", names(kegg_paths))
      matched_names <- kegg_paths[pathway_ids]
      path_t2n$NAME <- ifelse(!is.na(matched_names), matched_names, pathway_ids)
    }, silent = TRUE)
    t2n_list[["KEGG_Pathway"]] <- path_t2n
  }
  
  universe_genes <- setdiff(clean_ids(rownames(dds)), CONFIG$exclude_loci)
  
  for (p_sheet in names(potency_data_sheets)) {
    strain_name <- gsub("^Potency_", "", p_sheet)
    pot_df      <- potency_data_sheets[[p_sheet]]
    
    valid_flags <- unique(as.character(pot_df$meta_regulation))
    valid_flags <- valid_flags[!valid_flags %in% c("Not Significant")]
    
    ora_plots <- list()
    
    for (flag in valid_flags) {
      flag_genes <- pot_df %>% dplyr::filter(meta_regulation == flag) %>% pull(Gene)
      
      if (length(flag_genes) < CONFIG$ora_min_genes) next
      
      for (db_name in names(t2g_list)) {
        t2g <- t2g_list[[db_name]]
        t2n <- t2n_list[[db_name]]
        if (is.null(t2g) || nrow(t2g) == 0) next
        
        enrich_res <- tryCatch({
          clusterProfiler::enricher(
            gene          = flag_genes,
            universe      = universe_genes,
            TERM2GENE     = t2g,
            TERM2NAME     = t2n,
            pvalueCutoff  = CONFIG$ora_pval_cutoff,
            pAdjustMethod = "BH",
            qvalueCutoff  = CONFIG$ora_qval_cutoff,
            minGSSize     = 2,
            maxGSSize     = 500
          )
        }, error = function(e) NULL)
        
        if (!is.null(enrich_res) && nrow(as.data.frame(enrich_res)) > 0) {
          res_df <- as.data.frame(enrich_res) %>%
            rowwise() %>%
            mutate(
              GR_num = as.numeric(strsplit(as.character(GeneRatio), "/")[[1]][1]),
              GR_den = as.numeric(strsplit(as.character(GeneRatio), "/")[[1]][2]),
              BR_num = as.numeric(strsplit(as.character(BgRatio), "/")[[1]][1]),
              BR_den = as.numeric(strsplit(as.character(BgRatio), "/")[[1]][2]),
              FoldEnrichment = (GR_num / GR_den) / (BR_num / BR_den)
            ) %>%
            ungroup() %>%
            dplyr::select(-GR_num, -GR_den, -BR_num, -BR_den) %>%
            mutate(Database = db_name, Meta_Regulation = flag, Strain = strain_name)
          
          key_name <- paste0(strain_name, "_", gsub("[^A-Za-z0-9]+", "_", flag), "_", db_name)
          ora_results_list[[key_name]] <- res_df
          
          plot_df <- res_df %>%
            arrange(p.adjust) %>%
            head(15) %>%
            mutate(Description = factor(Description, levels = rev(unique(Description))))
          
          p_ora <- ggplot(plot_df, aes(x = FoldEnrichment, y = Description, size = Count, color = p.adjust)) +
            geom_point() +
            scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
            scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
            labs(
              title    = paste0(db_name, " Enrichment: ", flag),
              subtitle = sprintf("Strain: %s | Input Genes: %d", strain_name, length(flag_genes)),
              x        = "Fold Enrichment",
              y        = ""
            ) +
            theme_classic(base_size = 12, base_family = CONFIG$font_family) +
            theme(
              axis.text.y = element_text(size = 9, face = "bold", color = "black"),
              axis.text.x = element_text(face = "bold", color = "black"),
              axis.title  = element_text(face = "bold"),
              plot.title  = element_text(face = "bold", size = 13),
              legend.title = element_text(face = "bold", size = 10),
              legend.text  = element_text(size = 9)
            )
          
          ora_plots[[paste0(flag, "_", db_name)]] <- p_ora
        }
      }
    }
    
    if (length(ora_plots) > 0) {
      ora_pdf_path <- file.path(plot_dir, paste0("ORA_Enrichment_", strain_name, ".pdf"))
      pdf(ora_pdf_path, width = 10, height = 7.5, family = "Helvetica", useDingbats = FALSE)
      purrr::walk(ora_plots, print)
      dev.off()
    }
  }
} else {
  cat("EggNOG annotation file or GTF protein mappings missing. Skipping ORA.\n")
}

# ==============================================================================
# 12. Master Excel Compilation
# ==============================================================================
cat("Compiling Master Excel Workbook...\n")
wb <- createWorkbook()

add_annotated_sheet <- function(wb_obj, sheet_title, df, anno_df = NULL) {
  s_name <- safe_sheet(sheet_title)
  addWorksheet(wb_obj, s_name)
  if (!is.null(anno_df) && nrow(anno_df) > 0 && "Gene" %in% names(df) && !("Gene_Name" %in% names(df))) {
    df <- df %>% left_join(anno_df, by = "Gene") %>% relocate(Gene_Name, Product, .after = Gene)
  }
  writeData(wb_obj, s_name, df)
}

add_annotated_sheet(wb, "PCA_Coordinates", pca_df)
add_annotated_sheet(wb, "Potency_Summary", potency_summary)
walk(names(potency_data_sheets), ~ add_annotated_sheet(wb, .x, potency_data_sheets[[.x]], anno_dict))
add_annotated_sheet(wb, "Sieve_Summary", sieve_summary)
walk(names(sieve_data_sheets), ~ add_annotated_sheet(wb, .x, sieve_data_sheets[[.x]], anno_dict))

if (nrow(vf_heatmap_data) > 0) {
  add_annotated_sheet(wb, "VF_Heatmap_VST", vf_heatmap_data, anno_dict)
  add_annotated_sheet(wb, "VF_Heatmap_Zscores", vf_zscore_data, anno_dict)
}

add_annotated_sheet(wb, "TPM_All_Samples", tpm_df, anno_dict)
add_annotated_sheet(wb, "Zscores_All_Samples", master_zscore_df, anno_dict)

if (length(ora_results_list) > 0) {
  addWorksheet(wb, "ORA_Enrichment")
  combined_ora <- bind_rows(ora_results_list) %>% relocate(Strain, Meta_Regulation, Database)
  writeData(wb, "ORA_Enrichment", combined_ora)
}

for (res_name in names(all_results)) {
  df_out <- all_results[[res_name]] %>%
    mutate(Status = case_when(
      padj < CONFIG$padj_cutoff & log2FoldChange >  CONFIG$lfc_cutoff ~ "UP",
      padj < CONFIG$padj_cutoff & log2FoldChange < -CONFIG$lfc_cutoff ~ "DOWN",
      TRUE ~ "NS"
    ))
  add_annotated_sheet(wb, res_name, df_out, anno_dict)
}

saveWorkbook(wb, file.path(res_dir, CONFIG$master_excel_name), overwrite = TRUE)
cat("Pipeline complete! All outputs and plots saved to:", res_dir, "\n")
