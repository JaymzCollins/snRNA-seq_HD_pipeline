# ============================================================
# PSEUDOBULK DIFFERENTIAL EXPRESSION, THRESHOLD SENSITIVITY,
# AND EXTERNAL VALIDATION
#
# Runs pseudobulk differential expression for both comparisons
# (genotype and treatment), assesses the impact of applying
# uniform vs. cell-type-specific significance thresholds, exports
# a supplementary Excel workbook of full DE results, generates
# volcano plots, and validates the disease signature against two
# independent external sources (Murillo et al., 2025; Lee et al.,
# 2020).
#
# Requires:
#   - merged_seurat_post_cluster_QC.rds, produced by
#     04_clustering_annotation.R
#     (PROJECT_ROOT/Analysis/merged_seurat_post_cluster_QC.rds)
#   - murillo_genotype_pseudobulk.csv, murillo_treatment_pseudobulk.csv
#     (Murillo et al. 2025 supplementary data, exported to CSV)
#   - lee2020_R62_zQ175_combined.csv (Lee et al. 2020 supplementary
#     Table S2, exported to CSV)
#
# Produces:
#   - DE_results_all_comparisons_postQC.rds / .csv
#   - threshold_sensitivity_strict_applied_to_others.csv
#   - threshold_sensitivity_standard_applied_to_MSN.csv
#   - Supplementary_DE_results_all_celltypes.xlsx
#   - Volcano plots (treatment comparison, all cell types)
#   - MSN_genotype_murillo_concordance.csv
#   - MSN_treatment_murillo_concordance.csv
#   - Lee2020_concordance_all.csv
#   - concordance_barplot_main_celltypes.png
# ============================================================

# ============================================================
# EDIT THIS: set to your own project root before running
# ============================================================
PROJECT_ROOT <- "/scratch/SCWF00151/shared/c.d21085818/Dissertation"

library(DESeq2)
library(Seurat)

# NOTE: Graph-based clustering (FindClusters - 04_clustering_annotation.R) 
# is not guaranteed to be perfectly reproducible across independent sessions, 
# even with a fixed random seed - see README for details. The checkpoint file
# this script produces (merged_seurat_post_cluster_QC.rds) is provided
# directly alongside this code archive to ensure exact reproducibility
# of all downstream results.
# Load the QC-corrected checkpoint (produced by 04_clustering_annotation.R)
merged_seurat <- readRDS(file.path(PROJECT_ROOT, "Analysis/merged_seurat_post_cluster_QC.rds"))

group_lookup <- c(
  "WT_Cas9_1_v2" = "WT_Cas9", "WT_Cas9_2_v2" = "WT_Cas9", "WT_Cas9_3_v2" = "WT_Cas9",
  "WT_Cas9_sg_1_v2" = "WT_Cas9_sg", "WT_Cas9_sg_2_v2" = "WT_Cas9_sg", "WT_Cas9_sg_3_v2" = "WT_Cas9_sg",
  "R61_Cas9_1_v2" = "R61_Cas9", "R61_Cas9_2_v2" = "R61_Cas9", "R61_Cas9_3_v2" = "R61_Cas9",
  "R61_Cas9_sg_1_v2" = "R61_Cas9_sg", "R61_Cas9_sg_2_v2" = "R61_Cas9_sg", "R61_Cas9_sg_3_v2" = "R61_Cas9_sg"
)

# ------------------------------------------------------------
# Comparisons to run: name -> c(numerator_group, denominator_group)
# ------------------------------------------------------------
comparisons <- list(
  treatment_R61 = c("R61_Cas9_sg", "R61_Cas9"),   # treatment effect within R6/1
  genotype_untreated = c("R61_Cas9", "WT_Cas9")   # disease signature, untreated arm
)

build_pseudobulk <- function(seurat_obj, cell_type_name) {
  obj <- subset(seurat_obj, subset = cell_type_provisional == cell_type_name)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  samples <- unique(obj$orig.ident)
  pb <- matrix(0, nrow = nrow(counts), ncol = length(samples))
  rownames(pb) <- rownames(counts)
  colnames(pb) <- samples
  for (s in samples) {
    cells <- colnames(obj)[obj$orig.ident == s]
    pb[, s] <- Matrix::rowSums(counts[, cells, drop = FALSE])
  }
  return(pb)
}

apply_global_gene_filter <- function(pb) {
  min_samples <- max(3, ceiling(ncol(pb) * 0.5))
  keep_genes <- rowSums(pb >= 10) >= min_samples
  return(pb[keep_genes, ])
}

run_DESeq2 <- function(pb, comparison) {
  groups <- group_lookup[colnames(pb)]
  keep <- names(groups)[groups %in% comparison]
  counts_subset <- pb[, keep, drop = FALSE]
  groups_subset <- factor(groups[keep], levels = comparison)
  if (any(table(groups_subset) < 3)) return(NULL)
  coldata <- data.frame(row.names = keep, group = groups_subset)
  group_counts <- table(groups_subset)
  keep_genes_stage2 <- rowSums(counts_subset >= 10) >= max(2, min(group_counts))
  counts_subset <- counts_subset[keep_genes_stage2, ]
  dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData = coldata, design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", comparison[1], comparison[2]))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  return(res)
}

# ------------------------------------------------------------
# Non-standard gene symbol filter (Gm-prefixed, numeric-prefixed,
# *Rik predicted/unnamed genes)
# ------------------------------------------------------------
is_standard_symbol <- function(genes) {
  !grepl("^Gm[0-9]", genes) & !grepl("^[0-9]", genes) & !grepl("Rik$", genes)
}

# ------------------------------------------------------------
# Cell-type-specific significance thresholds
# ------------------------------------------------------------
MSN_PADJ_THRESHOLD <- 0.01
MSN_LFC_THRESHOLD <- 0.585

OTHER_PADJ_THRESHOLD <- 0.05
OTHER_LFC_THRESHOLD <- log2(1.2)  # ≈ 0.263

all_cell_types <- setdiff(unique(merged_seurat$cell_type_provisional), "Low_quality")
cat("Cell types to process (Low_quality excluded):", paste(all_cell_types, collapse = ", "), "\n\n")

de_results_list <- list()
de_summary <- data.frame(comparison = character(), cell_type = character(), n_tested = integer(),
                         n_significant_total = integer(),
                         n_significant_standard_symbols = integer(),
                         n_excluded_nonstandard = integer(),
                         n_up = integer(), n_down = integer(),
                         padj_threshold_used = numeric(), lfc_threshold_used = numeric(),
                         stringsAsFactors = FALSE)

for (comp_name in names(comparisons)) {
  comparison <- comparisons[[comp_name]]
  cat("########## COMPARISON:", comp_name, "(", comparison[1], "vs", comparison[2], ") ##########\n\n")
  de_results_list[[comp_name]] <- list()
  
  for (ct in all_cell_types) {
    cat("=== Processing:", ct, "===\n")
    
    if (ct == "MSN") {
      padj_thresh <- MSN_PADJ_THRESHOLD
      lfc_thresh <- MSN_LFC_THRESHOLD
    } else {
      padj_thresh <- OTHER_PADJ_THRESHOLD
      lfc_thresh <- OTHER_LFC_THRESHOLD
    }
    
    pb <- build_pseudobulk(merged_seurat, ct)
    pb <- apply_global_gene_filter(pb)
    res <- run_DESeq2(pb, comparison)
    if (is.null(res)) { cat("  Skipped - insufficient replicates\n\n"); next }
    
    res$significant <- !is.na(res$padj) & res$padj < padj_thresh & abs(res$log2FoldChange) >= lfc_thresh
    res$direction <- ifelse(res$significant & res$log2FoldChange > 0, "Up_in_treated",
                            ifelse(res$significant & res$log2FoldChange < 0, "Down_in_treated", "NS"))
    res$standard_symbol <- is_standard_symbol(res$gene)
    
    de_results_list[[comp_name]][[ct]] <- res
    
    n_sig_total <- sum(res$significant, na.rm = TRUE)
    n_sig_standard <- sum(res$significant & res$standard_symbol, na.rm = TRUE)
    n_up <- sum(res$direction == "Up_in_treated" & res$standard_symbol, na.rm = TRUE)
    n_down <- sum(res$direction == "Down_in_treated" & res$standard_symbol, na.rm = TRUE)
    
    cat("  Thresholds used: padj <", padj_thresh, "| |log2FC| >=", round(lfc_thresh, 4), "\n")
    cat("  Genes tested:", nrow(res), "| Significant (total):", n_sig_total,
        "| Significant (standard symbols only):", n_sig_standard,
        "(", n_up, "up /", n_down, "down )\n\n")
    
    de_summary <- rbind(de_summary, data.frame(comparison = comp_name, cell_type = ct, n_tested = nrow(res),
                                               n_significant_total = n_sig_total,
                                               n_significant_standard_symbols = n_sig_standard,
                                               n_excluded_nonstandard = n_sig_total - n_sig_standard,
                                               n_up = n_up, n_down = n_down,
                                               padj_threshold_used = padj_thresh, lfc_threshold_used = round(lfc_thresh, 4)))
  }
}

de_summary <- de_summary[order(de_summary$comparison, -de_summary$n_significant_standard_symbols), ]
print(de_summary)

saveRDS(de_results_list, file.path(PROJECT_ROOT, "Analysis/DE_results_all_comparisons_postQC.rds"))
write.csv(de_summary, file.path(PROJECT_ROOT, "Analysis/DE_summary_all_comparisons_postQC.csv"), row.names = FALSE)


# ============================================================
# THRESHOLD SENSITIVITY ANALYSIS
# Tests whether cell-type-specific significance thresholds were
# necessary rather than arbitrary, by examining the impact of a
# single, uniform threshold applied in both directions. Re-uses
# de_results_list built above - no new DESeq2 runs required.
# ============================================================

# ------------------------------------------------------------
# Apply MSN's stricter threshold to every other cell
# type's results, in place of their own standard threshold.
# Restricted to cell types with >=10 significant genes under
# their actual threshold, matching the minimum applied for
# pathway enrichment eligibility.
# ------------------------------------------------------------
threshold_sensitivity_strict <- data.frame(
  comparison = character(), cell_type = character(),
  n_sig_actual_threshold = integer(), n_sig_msn_threshold = integer(),
  pct_reduction = numeric(), stringsAsFactors = FALSE
)

for (comp_name in names(de_results_list)) {
  for (ct in names(de_results_list[[comp_name]])) {
    if (ct == "MSN") next  # MSN already uses this threshold
    
    res <- de_results_list[[comp_name]][[ct]]
    
    n_actual <- sum(res$significant & res$standard_symbol, na.rm = TRUE)
    if (n_actual < 10) next  # eligibility cutoff
    
    n_msn_thresh <- sum(!is.na(res$padj) & res$padj < MSN_PADJ_THRESHOLD &
                          abs(res$log2FoldChange) >= MSN_LFC_THRESHOLD &
                          res$standard_symbol, na.rm = TRUE)
    
    pct_reduction <- round(100 * (n_actual - n_msn_thresh) / n_actual, 1)
    
    threshold_sensitivity_strict <- rbind(threshold_sensitivity_strict, data.frame(
      comparison = comp_name, cell_type = ct,
      n_sig_actual_threshold = n_actual, n_sig_msn_threshold = n_msn_thresh,
      pct_reduction = pct_reduction
    ))
  }
}

threshold_sensitivity_strict <- threshold_sensitivity_strict[order(-threshold_sensitivity_strict$pct_reduction), ]
cat("=== Table B1: MSN's stricter threshold applied to other cell types ===\n")
print(threshold_sensitivity_strict)

write.csv(threshold_sensitivity_strict,
          file.path(PROJECT_ROOT, "Analysis/threshold_sensitivity_strict_applied_to_others.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# Apply the standard, more lenient threshold to MSN,
# in place of its own stricter threshold.
# ------------------------------------------------------------
threshold_sensitivity_lenient <- data.frame(
  comparison = character(), n_sig_actual_threshold = integer(),
  n_sig_standard_threshold = integer(), pct_increase = numeric(),
  stringsAsFactors = FALSE
)

for (comp_name in names(de_results_list)) {
  res <- de_results_list[[comp_name]][["MSN"]]
  
  n_actual <- sum(res$significant & res$standard_symbol, na.rm = TRUE)
  
  n_standard_thresh <- sum(!is.na(res$padj) & res$padj < OTHER_PADJ_THRESHOLD &
                             abs(res$log2FoldChange) >= OTHER_LFC_THRESHOLD &
                             res$standard_symbol, na.rm = TRUE)
  
  pct_increase <- round(100 * (n_standard_thresh - n_actual) / n_actual, 1)
  
  threshold_sensitivity_lenient <- rbind(threshold_sensitivity_lenient, data.frame(
    comparison = comp_name, n_sig_actual_threshold = n_actual,
    n_sig_standard_threshold = n_standard_thresh, pct_increase = pct_increase
  ))
}

cat("\n=== Table B2: Standard threshold applied to MSN ===\n")
print(threshold_sensitivity_lenient)

write.csv(threshold_sensitivity_lenient,
          file.path(PROJECT_ROOT, "Analysis/threshold_sensitivity_standard_applied_to_MSN.csv"),
          row.names = FALSE)


# ============================================================
# EXPORT FULL DE RESULTS TO SUPPLEMENTARY EXCEL WORKBOOK
#
# Produces one tab per cell type per comparison, containing
# every tested gene (regardless of significance or symbol type),
# plus a README tab summarising sheet contents and gene counts.
# ============================================================

library(openxlsx)

COMPARISON_LABELS <- c(
  genotype_untreated = "Genotype",
  treatment_R61 = "Treatment"
)

COMPARISON_DESCRIPTIONS <- c(
  genotype_untreated = "R6/1 vs WT (untreated)",
  treatment_R61 = "R6/1 Cas9+sgRNA vs Cas9-only"
)

# ------------------------------------------------------------
# Build every data sheet, tracking a summary row for each as
# it's added, so the README index stays exactly in sync with
# the actual sheet contents
# ------------------------------------------------------------

wb <- createWorkbook()
sheet_index <- data.frame(
  Sheet = character(), Comparison = character(), Cell_Type = character(),
  Total_Genes_Tested = integer(), Significant_Genes = integer(),
  stringsAsFactors = FALSE
)

for (comp_name in names(de_results_list)) {
  for (ct in names(de_results_list[[comp_name]])) {
    
    res <- de_results_list[[comp_name]][[ct]]
    res <- res[order(res$padj), c("gene", "log2FoldChange", "padj", "significant")]
    if (nrow(res) == 0) next
    
    sheet_name <- substr(paste0(COMPARISON_LABELS[[comp_name]], "_", ct), 1, 31)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, res)
    addStyle(wb, sheet_name, createStyle(textDecoration = "bold"), rows = 1, cols = 1:4)
    setColWidths(wb, sheet_name, cols = 1:4, widths = c(18, 15, 15, 12))
    
    sheet_index <- rbind(sheet_index, data.frame(
      Sheet = sheet_name,
      Comparison = COMPARISON_DESCRIPTIONS[[comp_name]],
      Cell_Type = ct,
      Total_Genes_Tested = nrow(res),
      Significant_Genes = sum(res$significant, na.rm = TRUE)
    ))
    
    cat(sprintf("Added sheet: %-25s (%d genes, %d significant)\n",
                sheet_name, nrow(res), sum(res$significant, na.rm = TRUE)))
  }
}

# ------------------------------------------------------------
# README tab: description text, followed by the sheet index
# ------------------------------------------------------------

readme_text <- c(
  "Supplementary Table: Full Differential Expression Results",
  "",
  "This workbook contains complete pseudobulk differential expression results",
  "for every annotated cell type, across two comparisons:",
  "  - Genotype: R6/1 vs WT (untreated Cas9-only arm) - disease-associated signature",
  "  - Treatment: R6/1 Cas9+sgRNA vs R6/1 Cas9-only - treatment-associated response",
  "",
  "Each sheet contains every gene tested in that cell type/comparison (regardless",
  "of significance or gene symbol type), with the following columns:",
  "  gene            - gene symbol",
  "  log2FoldChange  - log2 fold change",
  "  padj            - Benjamini-Hochberg adjusted p-value",
  "  significant     - TRUE/FALSE, based on cell-type-specific significance",
  "                    thresholds applied in this analysis (MSN: padj<0.01,",
  "                    |log2FC|>=0.585; all other cell types: padj<0.05,",
  "                    |log2FC|>=0.263)",
  "",
  "Sheet index and gene counts are listed below."
)

addWorksheet(wb, "README", tabColour = "yellow")
writeData(wb, "README", data.frame(README = readme_text), colNames = FALSE)
setColWidths(wb, "README", cols = 1, widths = 90)

index_start_row <- length(readme_text) + 3
writeData(wb, "README", sheet_index, startRow = index_start_row)
addStyle(wb, "README", createStyle(textDecoration = "bold"), rows = index_start_row, cols = 1:5)

worksheetOrder(wb) <- c(which(names(wb) == "README"), which(names(wb) != "README"))

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveWorkbook(wb, file.path(PROJECT_ROOT, "Analysis/Supplementary_DE_results_all_celltypes.xlsx"),
             overwrite = TRUE)
cat("\nWorkbook saved with README tab (", nrow(sheet_index), "data sheets ).\n")



# ============================================================
# VOLCANO PLOTS - TREATMENT COMPARISON, ALL CELL TYPES
# Standard-symbol genes only; GFP excluded from the plot for
# visual clarity (reported separately in text as the most
# significantly differentially expressed gene where applicable)
# ============================================================
library(ggplot2)
library(ggrepel)

plot_volcano <- function(de_results_list, comparison, cell_type, comparison_label,
                         padj_thresh, lfc_thresh, n_labels = 8) {
  
  res <- de_results_list[[comparison]][[cell_type]]
  if (is.null(res)) {
    cat(cell_type, ": no results available, skipping\n")
    return(NULL)
  }
  
  res <- res[res$standard_symbol, ]
  
  # Report GFP's statistics before excluding it from the plot
  gfp_row <- res[res$gene == "GFP", ]
  if (nrow(gfp_row) > 0) {
    cat(cell_type, "- GFP log2FC:", round(gfp_row$log2FoldChange, 3),
        "| padj:", format.pval(gfp_row$padj), "\n")
  }
  
  # Calculate up/down counts BEFORE excluding GFP, so counts remain accurate
  res$category <- ifelse(res$significant & res$log2FoldChange > 0, "Up",
                         ifelse(res$significant & res$log2FoldChange < 0, "Down", "Not significant"))
  n_up <- sum(res$category == "Up")
  n_down <- sum(res$category == "Down")
  
  # Exclude GFP from the plot only
  res_plot <- res[res$gene != "GFP", ]
  res_plot$category <- factor(res_plot$category, levels = c("Down", "Not significant", "Up"))
  
  top_genes <- res_plot[res_plot$significant, ]
  top_genes <- top_genes[order(top_genes$padj), ]
  top_genes <- top_genes[1:min(n_labels, nrow(top_genes)), ]
  
  p <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = category)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_genes, aes(label = gene), size = 3.2, color = "black",
                    max.overlaps = 20, box.padding = 0.4) +
    scale_color_manual(values = c("Down" = "#2166AC", "Not significant" = "grey70", "Up" = "#B2182B")) +
    labs(title = paste0(cell_type, ": ", comparison_label),
         subtitle = paste0("Down: ", n_down, "    Up: ", n_up),
         x = expression(Log[2]~"Fold Change"), y = expression(-log[10]~"adjusted p-value"),
         color = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 11),
          legend.position = "bottom")
  
  return(p)
}

MSN_PADJ_THRESHOLD <- 0.01
MSN_LFC_THRESHOLD <- 0.585
OTHER_PADJ_THRESHOLD <- 0.05
OTHER_LFC_THRESHOLD <- log2(1.2)

all_cell_types <- names(de_results_list[["treatment_R61"]])

for (ct in all_cell_types) {
  if (ct == "MSN") {
    padj_thresh <- MSN_PADJ_THRESHOLD
    lfc_thresh <- MSN_LFC_THRESHOLD
  } else {
    padj_thresh <- OTHER_PADJ_THRESHOLD
    lfc_thresh <- OTHER_LFC_THRESHOLD
  }
  
  p <- plot_volcano(de_results_list, "treatment_R61", ct,
                    comparison_label = "R6/1 Cas9+sgRNA vs Cas9-only",
                    padj_thresh = padj_thresh, lfc_thresh = lfc_thresh)
  
  if (is.null(p)) next
  
  print(p)
  
  ggsave(file.path(PROJECT_ROOT, paste0("Analysis/volcano_", ct, "_treatment_R61.png")),
         plot = p, width = 7, height = 6, dpi = 300)
  
  cat("Saved:", ct, "\n")
}


# ============================================================
# EXTERNAL VALIDATION: MURILLO ET AL. (2025) AND LEE ET AL. (2020)
# Direction concordance + Fisher's exact test.
# ============================================================

# ------------------------------------------------------------
# STEP 1: Import Murillo et al. genotype and treatment pseudobulk
# data (exported from their supplementary Excel file to CSV)
# ------------------------------------------------------------
murillo_geno <- read.csv(file.path(PROJECT_ROOT, "Analysis/murillo_genotype_pseudobulk.csv"))
murillo_treat <- read.csv(file.path(PROJECT_ROOT, "Analysis/murillo_treatment_pseudobulk.csv"))

# ------------------------------------------------------------
# STEP 2: Import Lee et al. (2020) R6/2 and zQ175DN DEG data
# (exported from their supplementary Table S2 to CSV)
# ------------------------------------------------------------
lee2020_combined <- read.csv(file.path(PROJECT_ROOT, "Analysis/lee2020_R62_zQ175_combined.csv"))

# ------------------------------------------------------------
# STEP 3: MURILLO CONCORDANCE - GENOTYPE COMPARISON, MSN
# ------------------------------------------------------------
geno_msn_full <- de_results_list[["genotype_untreated"]][["MSN"]]
your_sig_genes <- geno_msn_full[geno_msn_full$significant, ]

merged <- merge(your_sig_genes[, c("gene", "log2FoldChange", "padj")],
                murillo_geno[, c("gene", "avg_log2FC", "p_val_adj", "Significance")],
                by = "gene")

merged$same_direction <- sign(merged$log2FoldChange) == sign(merged$avg_log2FC)
n_concordant_geno <- sum(merged$same_direction)
pct_concordant_geno <- round(100 * n_concordant_geno / nrow(merged), 1)

contingency_geno <- matrix(c(n_concordant_geno, nrow(merged) - n_concordant_geno,
                             nrow(merged) / 2, nrow(merged) / 2), nrow = 2, byrow = TRUE)
fisher_geno <- fisher.test(contingency_geno)

cat("=== Murillo - MSN Genotype ===\n")
cat("n =", nrow(merged), "| Concordant:", n_concordant_geno, "(", pct_concordant_geno, "%) | Fisher p =",
    format.pval(fisher_geno$p.value), "\n\n")

# ------------------------------------------------------------
# STEP 4: MURILLO CONCORDANCE - TREATMENT COMPARISON, MSN
# ------------------------------------------------------------
treat_msn_full <- de_results_list[["treatment_R61"]][["MSN"]]
your_sig_genes_treat <- treat_msn_full[treat_msn_full$significant, ]

merged_treat <- merge(your_sig_genes_treat[, c("gene", "log2FoldChange", "padj")],
                      murillo_treat[, c("gene", "avg_log2FC", "p_val_adj", "Significance")],
                      by = "gene")

merged_treat$same_direction <- sign(merged_treat$log2FoldChange) == sign(merged_treat$avg_log2FC)
n_concordant_treat <- sum(merged_treat$same_direction)
pct_concordant_treat <- round(100 * n_concordant_treat / nrow(merged_treat), 1)

contingency_treat <- matrix(c(n_concordant_treat, nrow(merged_treat) - n_concordant_treat,
                              nrow(merged_treat) / 2, nrow(merged_treat) / 2), nrow = 2, byrow = TRUE)
fisher_treat <- fisher.test(contingency_treat)

cat("=== Murillo - MSN Treatment ===\n")
cat("n =", nrow(merged_treat), "| Concordant:", n_concordant_treat, "(", pct_concordant_treat, "%) | Fisher p =",
    format.pval(fisher_treat$p.value), "\n\n")

# ------------------------------------------------------------
# STEP 5: LEE ET AL. CONCORDANCE - R6/2 AND zQ175DN MODELS,
# GENOTYPE COMPARISON, MULTIPLE CELL TYPES
# ------------------------------------------------------------
celltype_mapping_lee <- list(
  MSN = c("dSPN", "iSPN"),
  Astrocytes = "Astroglia",
  Microglia = "Microglia",
  Oligodendrocytes = "Oligodendrocyte",
  OPC = "OPC",
  PV_Th_Interneuron = "Pvalb_Th_IN",
  Sst_Npy_Interneuron = "Sst_Npy_IN",
  Ciliated_Ependymal = "Ciliated_Ependymal"
)

run_lee_concordance <- function(model_prefix) {
  results <- data.frame(model = character(), cell_type = character(),
                        n_compared = integer(), n_concordant = integer(),
                        pct_concordant = numeric(), fisher_p = numeric(),
                        stringsAsFactors = FALSE)
  
  for (my_ct in names(celltype_mapping_lee)) {
    their_sheets <- paste0(model_prefix, "_", celltype_mapping_lee[[my_ct]])
    geno_res <- de_results_list[["genotype_untreated"]][[my_ct]]
    if (is.null(geno_res)) next
    your_sig <- geno_res[geno_res$significant, c("gene", "log2FoldChange")]
    
    their_data <- lee2020_combined[lee2020_combined$sheet %in% their_sheets, ]
    if (nrow(their_data) == 0) next
    their_data <- aggregate(their_data$summary.logFC, by = list(gene = their_data$gene), FUN = mean)
    colnames(their_data) <- c("gene", "logFC_lee")
    
    merged_lee <- merge(your_sig, their_data, by = "gene")
    if (nrow(merged_lee) < 2) next
    
    merged_lee$same_direction <- sign(merged_lee$log2FoldChange) == sign(merged_lee$logFC_lee)
    n_concordant <- sum(merged_lee$same_direction)
    n_discordant <- sum(!merged_lee$same_direction)
    pct_concordant <- round(100 * n_concordant / nrow(merged_lee), 1)
    
    contingency <- matrix(c(n_concordant, n_discordant, nrow(merged_lee)/2, nrow(merged_lee)/2),
                          nrow = 2, byrow = TRUE)
    fisher_p <- fisher.test(contingency)$p.value
    
    cat(model_prefix, "-", my_ct, ": n =", nrow(merged_lee), "| concordant =", pct_concordant,
        "% | Fisher p =", format.pval(fisher_p), "\n")
    
    results <- rbind(results, data.frame(model = model_prefix, cell_type = my_ct,
                                         n_compared = nrow(merged_lee), n_concordant = n_concordant,
                                         pct_concordant = pct_concordant, fisher_p = fisher_p))
  }
  return(results)
}

cat("=== Lee et al. (2020) - R6/2 model ===\n")
r62_results <- run_lee_concordance("R62")

cat("\n=== Lee et al. (2020) - zQ175DN model ===\n")
zq175_results <- run_lee_concordance("zQ175")

lee2020_all_results <- rbind(r62_results, zq175_results)

# ------------------------------------------------------------
# Save all validation results
# ------------------------------------------------------------
write.csv(merged, file.path(PROJECT_ROOT, "Analysis/MSN_genotype_murillo_concordance.csv"), row.names = FALSE)
write.csv(merged_treat, file.path(PROJECT_ROOT, "Analysis/MSN_treatment_murillo_concordance.csv"), row.names = FALSE)
write.csv(lee2020_all_results, file.path(PROJECT_ROOT, "Analysis/Lee2020_concordance_all.csv"), row.names = FALSE)

cat("\n=== ALL VALIDATION RESULTS ===\n")
print(lee2020_all_results)


# ============================================================
# BAR CHART: DIRECTION CONCORDANCE ACROSS VALIDATION SOURCES
# MSN, Astrocytes, Oligodendrocytes, Ciliated Ependymal
# ============================================================
library(ggplot2)

concordance_data <- data.frame(
  cell_type = factor(rep(c("MSN", "Astrocytes", "Oligodendrocytes", "Ciliated Ependymal"), each = 3),
                     levels = c("MSN", "Astrocytes", "Oligodendrocytes", "Ciliated Ependymal")),
  source = factor(rep(c("Murillo (same dataset)", "Lee R6/2", "Lee zQ175DN"), times = 4),
                  levels = c("Murillo (same dataset)", "Lee R6/2", "Lee zQ175DN")),
  concordance = c(97.5, 95.9, 87.8,
                  NA, 100, 63.5,
                  NA, 98.7, 88.5,
                  NA, 76.8, 54.4)
)

p <- ggplot(concordance_data, aes(x = cell_type, y = concordance, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, na.rm = TRUE) +
  scale_fill_manual(values = c("Murillo (same dataset)" = "#4472C4",
                               "Lee R6/2" = "#ED7D31",
                               "Lee zQ175DN" = "#548235"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "Direction concordance of disease signature vs independent validation sources",
       x = NULL, y = "Direction concordance (%)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size = 11))

print(p)

ggsave(file.path(PROJECT_ROOT, "Analysis/concordance_barplot_main_celltypes.png"),
       plot = p, width = 11, height = 6.5, dpi = 300)