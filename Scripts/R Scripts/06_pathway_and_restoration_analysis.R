# ============================================================
# OFF-TARGET SPECIFICITY, INTERFERON HEATMAP, TREATMENT PATHWAY
# ENRICHMENT, AND RESTORATION ANALYSIS
#
# Tests whether the treatment vector itself (independent of
# disease genotype) produces transcriptional changes, visualises
# the resulting interferon-stimulated gene signature across all
# four experimental groups, runs pathway enrichment on the
# treatment comparison, and performs the restoration analysis
# (treated R6/1 vs WT) to determine whether treatment restores
# disease-associated gene expression toward wild-type levels.
#
# Requires:
#   - merged_seurat_post_cluster_QC.rds, produced by
#     04_clustering_annotation.R
#     (PROJECT_ROOT/Analysis/merged_seurat_post_cluster_QC.rds)
#   - DE_results_all_comparisons_postQC.rds, produced by
#     05_pseudobulk_DE_and_validation.R
#
# Produces:
#   - offtarget_WTsg_vs_WT_all_celltypes.rds / offtarget_summary_all_celltypes.csv
#   - heatmap_Oligodendrocytes_interferon_4groups_raw.png
#   - pathway_results_treatment_all_celltypes.rds + pathway bar charts
#   - restoration_results_all_celltypes.rds
#   - pathway_results_restoration_all_celltypes.rds + pathway bar charts
#   - restoration_summary_all_celltypes.csv
#   - shrinkage_summary_all_celltypes.csv
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
# of all downstream results; re-running Script 04 is not required and
# may not reproduce identical cluster assignments.
# Load the QC-corrected checkpoint (produced by 04_clustering_annotation.R)
merged_seurat <- readRDS(file.path(PROJECT_ROOT, "Analysis/merged_seurat_post_cluster_QC.rds"))

# Load the DE results (produced by 05_pseudobulk_DE_and_validation.R)
de_results_list <- readRDS(file.path(PROJECT_ROOT, "Analysis/DE_results_all_comparisons_postQC.rds"))

# ------------------------------------------------------------
# Shared helper functions, needed by the off-target check and
# the interferon heatmap below (also redefined further down,
# immediately ahead of the restoration analysis, for that
# section's own standalone clarity)
# ------------------------------------------------------------
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

# ============================================================
# OFF-TARGET SPECIFICITY CHECK - GROUP 2 vs GROUP 1
# WT + Cas9D10A + sgCTG-Mut+5 vs WT + Cas9D10A only
# Tests whether the treatment machinery itself (independent of
# disease genotype) produces transcriptional changes - run for
# every cell type, using each cell type's own established
# significance threshold
# ============================================================
group_lookup <- c(
  "WT_Cas9_1_v2" = "WT_Cas9", "WT_Cas9_2_v2" = "WT_Cas9", "WT_Cas9_3_v2" = "WT_Cas9",
  "WT_Cas9_sg_1_v2" = "WT_Cas9_sg", "WT_Cas9_sg_2_v2" = "WT_Cas9_sg", "WT_Cas9_sg_3_v2" = "WT_Cas9_sg"
)
is_standard_symbol <- function(genes) {
  !grepl("^Gm[0-9]", genes) & !grepl("^[0-9]", genes) & !grepl("Rik$", genes)
}
all_cell_types <- setdiff(unique(merged_seurat$cell_type_provisional), "Low_quality")
offtarget_results <- list()
offtarget_summary <- data.frame(cell_type = character(), n_tested = integer(),
                                n_significant = integer(), stringsAsFactors = FALSE)
for (ct in all_cell_types) {
  cat("\n=== Processing:", ct, "===\n")
  
  if (ct == "MSN") {
    padj_thresh <- 0.01
    lfc_thresh <- 0.585
  } else {
    padj_thresh <- 0.05
    lfc_thresh <- log2(1.2)
  }
  
  pb <- build_pseudobulk(merged_seurat, ct)
  pb <- apply_global_gene_filter(pb)
  
  wt_samples <- names(group_lookup)
  wt_samples <- intersect(wt_samples, colnames(pb))
  
  if (length(wt_samples) < 4) {
    cat("  Skipped - insufficient samples\n")
    next
  }
  
  pb_subset <- pb[, wt_samples, drop = FALSE]
  groups_subset <- factor(group_lookup[wt_samples], levels = c("WT_Cas9", "WT_Cas9_sg"))
  
  if (any(table(groups_subset) < 3)) {
    cat("  Skipped - insufficient replicates\n")
    next
  }
  
  coldata <- data.frame(row.names = wt_samples, group = groups_subset)
  group_counts <- table(groups_subset)
  keep_genes_stage2 <- rowSums(pb_subset >= 10) >= max(2, min(group_counts))
  pb_subset <- pb_subset[keep_genes_stage2, ]
  
  dds <- tryCatch({
    d <- DESeqDataSetFromMatrix(countData = pb_subset, colData = coldata, design = ~ group)
    DESeq(d)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(dds)) next
  
  res <- as.data.frame(results(dds, contrast = c("group", "WT_Cas9_sg", "WT_Cas9")))
  res$gene <- rownames(res)
  res$significant <- !is.na(res$padj) & res$padj < padj_thresh & abs(res$log2FoldChange) >= lfc_thresh
  res$standard_symbol <- is_standard_symbol(res$gene)
  
  n_sig <- sum(res$significant & res$standard_symbol, na.rm = TRUE)
  
  cat("  Genes tested:", nrow(res), "| Significant:", n_sig, "\n")
  
  offtarget_results[[ct]] <- res
  offtarget_summary <- rbind(offtarget_summary, data.frame(
    cell_type = ct, n_tested = nrow(res), n_significant = n_sig
  ))
}
print(offtarget_summary)
saveRDS(offtarget_results, file.path(PROJECT_ROOT, "Analysis/offtarget_WTsg_vs_WT_all_celltypes.rds"))
write.csv(offtarget_summary, file.path(PROJECT_ROOT, "Analysis/offtarget_summary_all_celltypes.csv"), row.names = FALSE)
for (ct in c("MSN", "PV_Th_Interneuron", "Endothelial", "Ciliated_Ependymal", "Microglia", "OPC", "Secretory_Ependymal", "Mural", "Sst_Npy_Interneuron")) {
  offtarget_ct <- offtarget_results[[ct]]
  if (is.null(offtarget_ct)) next
  offtarget_sig_genes <- offtarget_ct$gene[offtarget_ct$significant & offtarget_ct$standard_symbol]
  treatment_sig_genes <- de_results_list[["treatment_R61"]][[ct]]$gene[
    de_results_list[["treatment_R61"]][[ct]]$significant & de_results_list[["treatment_R61"]][[ct]]$standard_symbol]
  overlap <- intersect(offtarget_sig_genes, treatment_sig_genes)
  cat(ct, "- off-target:", length(offtarget_sig_genes), "| treatment:", length(treatment_sig_genes), "| overlap:", length(overlap), "\n")
  if (length(overlap) > 0) print(overlap)
}



# ============================================================
# HEATMAP: OLIGODENDROCYTE INTERFERON-STIMULATED GENES,
# ALL FOUR EXPERIMENTAL GROUPS
# Log-transformed, size-factor normalised expression (not
# VST/z-scored), across WT_Cas9, WT_Cas9_sg, R61_Cas9, and
# R61_Cas9_sg, for the eight interferon-stimulated genes
# identified in the Oligodendrocyte treatment comparison.
# Run after the off-target specificity check script, which
# already provides build_pseudobulk(), apply_global_gene_filter(),
# and merged_seurat.
# ============================================================
library(ggplot2)
library(reshape2)

group_lookup_full <- c(
  "WT_Cas9_1_v2" = "WT_Cas9", "WT_Cas9_2_v2" = "WT_Cas9", "WT_Cas9_3_v2" = "WT_Cas9",
  "WT_Cas9_sg_1_v2" = "WT_Cas9_sg", "WT_Cas9_sg_2_v2" = "WT_Cas9_sg", "WT_Cas9_sg_3_v2" = "WT_Cas9_sg",
  "R61_Cas9_1_v2" = "R61_Cas9", "R61_Cas9_2_v2" = "R61_Cas9", "R61_Cas9_3_v2" = "R61_Cas9",
  "R61_Cas9_sg_1_v2" = "R61_Cas9_sg", "R61_Cas9_sg_2_v2" = "R61_Cas9_sg", "R61_Cas9_sg_3_v2" = "R61_Cas9_sg"
)

interferon_genes <- c("Iigp1c", "H2-Q4", "Irgm2", "Trim12a", "Trim34a", "Nlrc5", "Ddx60", "Ifih1")

pb_oligo <- build_pseudobulk(merged_seurat, "Oligodendrocytes")
pb_oligo <- apply_global_gene_filter(pb_oligo)

coldata_full <- data.frame(row.names = colnames(pb_oligo),
                           group = factor(group_lookup_full[colnames(pb_oligo)],
                                          levels = c("WT_Cas9", "WT_Cas9_sg", "R61_Cas9", "R61_Cas9_sg")))

dds_vis <- DESeqDataSetFromMatrix(countData = pb_oligo, colData = coldata_full, design = ~ group)
dds_vis <- estimateSizeFactors(dds_vis)
norm_counts <- counts(dds_vis, normalized = TRUE)

genes_available <- intersect(interferon_genes, rownames(norm_counts))
missing_genes <- setdiff(interferon_genes, genes_available)
if (length(missing_genes) > 0) {
  cat("Note: genes not available post-filtering:", paste(missing_genes, collapse = ", "), "\n")
}

log_norm <- log2(norm_counts[genes_available, , drop = FALSE] + 1)

heat_df <- melt(log_norm, varnames = c("gene", "sample"), value.name = "log2_norm_expr")
heat_df$group <- group_lookup_full[as.character(heat_df$sample)]
heat_df$group <- factor(heat_df$group, levels = c("WT_Cas9", "WT_Cas9_sg", "R61_Cas9", "R61_Cas9_sg"))
heat_df$gene <- factor(heat_df$gene, levels = rev(interferon_genes[interferon_genes %in% genes_available]))

p <- ggplot(heat_df, aes(x = sample, y = gene, fill = log2_norm_expr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Log2 normalised\nexpression") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(title = "Oligodendrocytes: interferon-stimulated gene expression",
       subtitle = "All four experimental groups",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        strip.text.x = element_text(face = "bold", size = 9),
        panel.grid = element_blank())

print(p)

ggsave(file.path(PROJECT_ROOT, "Analysis/heatmap_Oligodendrocytes_interferon_4groups_raw.png"),
       plot = p, width = 10, height = 5, dpi = 300)

cat("Saved heatmap for", length(genes_available), "genes across all four groups.\n")



# ============================================================
# TOP 10 PATHWAY TERMS - TREATMENT COMPARISON, ALL CELL TYPES
# ============================================================
library(enrichR)
library(ggplot2)

websiteLive <- getOption("enrichR.live")
if (is.null(websiteLive) || !websiteLive) setEnrichrSite("Enrichr")

dbs_to_query <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                  "GO_Molecular_Function_2023", "KEGG_2019_Mouse")

MIN_SIG_GENES <- 10
all_cell_types <- names(de_results_list[["treatment_R61"]])

pathway_results_all <- list()

for (ct in all_cell_types) {
  res_ct <- de_results_list[["treatment_R61"]][[ct]]
  sig_genes <- res_ct$gene[res_ct$significant & res_ct$standard_symbol]
  
  if (length(sig_genes) < MIN_SIG_GENES) {
    cat(ct, ": skipped (", length(sig_genes), "significant genes, below", MIN_SIG_GENES, ")\n")
    next
  }
  
  cat("\n=== Running enrichment:", ct, "( n =", length(sig_genes), ") ===\n")
  
  enriched <- tryCatch({
    enrichr(sig_genes, dbs_to_query)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(enriched)) next
  
  all_terms <- do.call(rbind, lapply(names(enriched), function(db) {
    df <- enriched[[db]]
    if (nrow(df) == 0) return(NULL)
    df$database <- db
    df
  }))
  all_terms <- all_terms[order(all_terms$Adjusted.P.value), ]
  
  pathway_results_all[[ct]] <- all_terms
  
  # Build and save the bar chart
  top_terms <- head(all_terms, 10)
  top_terms$Term_clean <- gsub("\\s*\\(GO:[0-9]+\\)", "", top_terms$Term)
  top_terms$Term_clean <- make.unique(top_terms$Term_clean)
  top_terms$Term_clean <- sapply(top_terms$Term_clean, function(x) paste(strwrap(x, width = 28), collapse = "\n"))
  top_terms$Term_clean <- factor(top_terms$Term_clean,
                                 levels = top_terms$Term_clean[order(top_terms$Adjusted.P.value, decreasing = TRUE)])
  
  p <- ggplot(top_terms, aes(y = Term_clean, x = -log10(Adjusted.P.value))) +
    geom_col(fill = "#A6CEE3", width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(title = paste0(ct, ": R6/1 Cas9+sgRNA vs Cas9-only"),
         subtitle = "Top 10 pathway terms (all databases)",
         y = NULL, x = expression(-log[10]~"(adjusted P-value)")) +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(PROJECT_ROOT, paste0("Analysis/pathway_bar_", ct, "_treatment_R61_top10.png")),
         plot = p, width = 10, height = 7, dpi = 300)
  
  cat("  Saved:", ct, "\n")
}

saveRDS(pathway_results_all, file.path(PROJECT_ROOT, "Analysis/pathway_results_treatment_all_celltypes.rds"))


# ============================================================
# PRINT TOP 15 PATHWAY TERMS - ALL PROCESSED CELL TYPES
# ============================================================

for (ct in names(pathway_results_all)) {
  cat("\n========================================\n")
  cat(ct, "- Top 15 pathway terms\n")
  cat("========================================\n")
  print(head(pathway_results_all[[ct]][, c("database", "Term", "Overlap", "Adjusted.P.value")], 15))
}



# Confirm the actual current Astrocyte pathway enrichment result
astro_pathways <- pathway_results_all[["Astrocytes"]]
print(head(astro_pathways[, c("database", "Term", "Overlap", "Adjusted.P.value")], 10))




# ============================================================
# Restoration Analysis: R6/1 TREATED vs WT
# Direct test of whether treatment restores gene expression to
# WT levels, across all annotated cell types. Genes significant
# in the genotype comparison (R6/1 vs WT, untreated) that become
# NON-significant here are considered restored; those remaining
# significant are not.
# ============================================================

library(DESeq2)

group_lookup <- c(
  "WT_Cas9_1_v2" = "WT_Cas9", "WT_Cas9_2_v2" = "WT_Cas9", "WT_Cas9_3_v2" = "WT_Cas9",
  "R61_Cas9_sg_1_v2" = "R61_Cas9_sg", "R61_Cas9_sg_2_v2" = "R61_Cas9_sg", "R61_Cas9_sg_3_v2" = "R61_Cas9_sg"
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

is_standard_symbol <- function(genes) {
  !grepl("^Gm[0-9]", genes) & !grepl("^[0-9]", genes) & !grepl("Rik$", genes)
}

MSN_PADJ_THRESHOLD <- 0.01
MSN_LFC_THRESHOLD <- 0.585
OTHER_PADJ_THRESHOLD <- 0.05
OTHER_LFC_THRESHOLD <- log2(1.2)

all_cell_types <- setdiff(unique(merged_seurat$cell_type_provisional), "Low_quality")

restoration_results_all <- list()

for (ct in all_cell_types) {
  cat("=== Processing (restoration):", ct, "===\n")
  
  padj_thresh <- if (ct == "MSN") MSN_PADJ_THRESHOLD else OTHER_PADJ_THRESHOLD
  lfc_thresh <- if (ct == "MSN") MSN_LFC_THRESHOLD else OTHER_LFC_THRESHOLD
  
  pb <- build_pseudobulk(merged_seurat, ct)
  pb <- apply_global_gene_filter(pb)
  
  keep_samples <- intersect(names(group_lookup), colnames(pb))
  pb_subset <- pb[, keep_samples, drop = FALSE]
  groups_subset <- factor(group_lookup[keep_samples], levels = c("R61_Cas9_sg", "WT_Cas9"))
  
  if (any(table(groups_subset) < 3)) {
    cat("  Skipped - insufficient replicates\n\n")
    next
  }
  
  coldata <- data.frame(row.names = keep_samples, group = groups_subset)
  group_counts <- table(groups_subset)
  keep_genes_stage2 <- rowSums(pb_subset >= 10) >= max(2, min(group_counts))
  pb_subset <- pb_subset[keep_genes_stage2, ]
  
  dds <- DESeqDataSetFromMatrix(countData = pb_subset, colData = coldata, design = ~ group)
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, contrast = c("group", "R61_Cas9_sg", "WT_Cas9")))
  res$gene <- rownames(res)
  res$significant <- !is.na(res$padj) & res$padj < padj_thresh & abs(res$log2FoldChange) >= lfc_thresh
  res$standard_symbol <- is_standard_symbol(res$gene)
  
  restoration_results_all[[ct]] <- res
  
  cat("  Genes tested:", nrow(res), "| Significant:", sum(res$significant, na.rm = TRUE), "\n\n")
}

saveRDS(restoration_results_all, file.path(PROJECT_ROOT, "Analysis/restoration_results_all_celltypes.rds"))

cat("\nDone. Restoration comparison completed for:", paste(names(restoration_results_all), collapse = ", "), "\n")



# ============================================================
# PATHWAY ENRICHMENT - RESTORED GENES, ALL CELL TYPES
# For each cell type, identifies restored genes (genotype-
# significant genes that lose significance in the restoration
# comparison) and runs pathway enrichment on that gene set.
# ============================================================
library(enrichR)
library(ggplot2)

websiteLive <- getOption("enrichR.live")
if (is.null(websiteLive) || !websiteLive) setEnrichrSite("Enrichr")

dbs_to_query <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023",
                  "GO_Molecular_Function_2023", "KEGG_2019_Mouse")

MIN_SIG_GENES <- 10
pathway_results_restoration_all <- list()

for (ct in names(restoration_results_all)) {
  
  geno_res <- de_results_list[["genotype_untreated"]][[ct]]
  if (is.null(geno_res)) {
    cat(ct, ": no genotype comparison results available, skipping\n")
    next
  }
  
  geno_sig <- geno_res$gene[geno_res$significant & geno_res$standard_symbol]
  
  restoration_res <- restoration_results_all[[ct]]
  restored_genes <- restoration_res$gene[restoration_res$gene %in% geno_sig &
                                           !restoration_res$significant]
  
  if (length(restored_genes) < MIN_SIG_GENES) {
    cat(ct, ": skipped (", length(restored_genes), "restored genes, below", MIN_SIG_GENES, ")\n")
    next
  }
  
  cat("\n=== Running enrichment on restored genes:", ct, "( n =", length(restored_genes), ") ===\n")
  
  enriched <- tryCatch({
    enrichr(restored_genes, dbs_to_query)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(enriched)) next
  
  all_terms <- do.call(rbind, lapply(names(enriched), function(db) {
    df <- enriched[[db]]
    if (nrow(df) == 0) return(NULL)
    df$database <- db
    df
  }))
  all_terms <- all_terms[order(all_terms$Adjusted.P.value), ]
  
  pathway_results_restoration_all[[ct]] <- all_terms
  
  # Build and save the bar chart
  top_terms <- head(all_terms, 10)
  top_terms$Term_clean <- gsub("\\s*\\(GO:[0-9]+\\)", "", top_terms$Term)
  top_terms$Term_clean <- make.unique(top_terms$Term_clean)
  top_terms$Term_clean <- sapply(top_terms$Term_clean, function(x) paste(strwrap(x, width = 28), collapse = "\n"))
  top_terms$Term_clean <- factor(top_terms$Term_clean,
                                 levels = top_terms$Term_clean[order(top_terms$Adjusted.P.value, decreasing = TRUE)])
  
  p <- ggplot(top_terms, aes(y = Term_clean, x = -log10(Adjusted.P.value))) +
    geom_col(fill = "#A6CEE3", width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(title = paste0(ct, ": Restored genes (treated vs WT, no longer significant)"),
         subtitle = "Top 10 pathway terms (all databases)",
         y = NULL, x = expression(-log[10]~"(adjusted P-value)")) +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(PROJECT_ROOT, paste0("Analysis/pathway_bar_", ct, "_restored_genes_top10.png")),
         plot = p, width = 10, height = 7, dpi = 300)
  
  cat("  Saved:", ct, "\n")
}

saveRDS(pathway_results_restoration_all, file.path(PROJECT_ROOT, "Analysis/pathway_results_restoration_all_celltypes.rds"))

cat("\nDone. Restored-gene pathway enrichment completed for:",
    paste(names(pathway_results_restoration_all), collapse = ", "), "\n")


# ============================================================
# PRINT RESTORATION SUMMARY - ALL CELL TYPES
# Compare against genotype-significant genes (standard symbols
# only, matching convention used throughout this dissertation)
# to calculate n restored, n not restored, % restored per cell type
# ============================================================
for (ct in names(restoration_results_all)) {
  geno_res <- de_results_list[["genotype_untreated"]][[ct]]
  if (is.null(geno_res)) {
    cat(ct, ": no genotype comparison results available, skipping\n\n")
    next
  }
  
  geno_sig_genes <- geno_res$gene[geno_res$significant & geno_res$standard_symbol]
  
  restoration_res <- restoration_results_all[[ct]]
  restoration_tested <- restoration_res[restoration_res$gene %in% geno_sig_genes, ]
  
  n_tested <- nrow(restoration_tested)
  n_restored <- sum(!restoration_tested$significant, na.rm = TRUE)
  n_not_restored <- sum(restoration_tested$significant, na.rm = TRUE)
  pct_restored <- if (n_tested > 0) round(100 * n_restored / n_tested, 1) else NA
  
  cat("===", ct, "===\n")
  cat("Genotype-significant genes:", length(geno_sig_genes), "\n")
  cat("Tested in restoration comparison:", n_tested, "\n")
  cat("Restored:", n_restored, "(", pct_restored, "%)\n")
  cat("Not restored:", n_not_restored, "\n\n")
}


# ============================================================
# RESTORATION SUMMARY TABLE - ALL CELL TYPES
# Standard-symbol genes only, matching established convention
# ============================================================

restoration_summary_table <- data.frame(
  cell_type = character(), genotype_sig_genes = integer(),
  tested_in_restoration = integer(), n_restored = integer(),
  pct_restored = numeric(), stringsAsFactors = FALSE
)

for (ct in names(restoration_results_all)) {
  geno_res <- de_results_list[["genotype_untreated"]][[ct]]
  if (is.null(geno_res)) next
  
  geno_sig_genes <- geno_res$gene[geno_res$significant & geno_res$standard_symbol]
  
  restoration_res <- restoration_results_all[[ct]]
  restoration_tested <- restoration_res[restoration_res$gene %in% geno_sig_genes, ]
  
  n_tested <- nrow(restoration_tested)
  n_restored <- sum(!restoration_tested$significant, na.rm = TRUE)
  pct_restored <- if (n_tested > 0) round(100 * n_restored / n_tested, 1) else NA
  
  restoration_summary_table <- rbind(restoration_summary_table, data.frame(
    cell_type = ct,
    genotype_sig_genes = length(geno_sig_genes),
    tested_in_restoration = n_tested,
    n_restored = n_restored,
    pct_restored = pct_restored
  ))
}

restoration_summary_table <- restoration_summary_table[order(-restoration_summary_table$genotype_sig_genes), ]

print(restoration_summary_table)

write.csv(restoration_summary_table, file.path(PROJECT_ROOT, "Analysis/restoration_summary_all_celltypes.csv"), row.names = FALSE)



# ============================================================
# MEDIAN EFFECT-SIZE SHRINKAGE, ALL CELL TYPES
# For each cell type's restored genes: compare genotype log2FC
# against restoration log2FC to calculate % shrinkage, then
# summarise as a median per cell type
# ============================================================

shrinkage_summary_table <- data.frame(
  cell_type = character(), n_restored = integer(),
  median_pct_shrinkage = numeric(), n_low_shrinkage = integer(),
  stringsAsFactors = FALSE
)

for (ct in names(restoration_results_all)) {
  geno_res <- de_results_list[["genotype_untreated"]][[ct]]
  if (is.null(geno_res)) next
  
  geno_sig_genes <- geno_res$gene[geno_res$significant & geno_res$standard_symbol]
  
  restoration_res <- restoration_results_all[[ct]]
  restored <- restoration_res[restoration_res$gene %in% geno_sig_genes &
                                !restoration_res$significant, ]
  
  if (nrow(restored) == 0) next
  
  # Merge in genotype log2FC for each restored gene
  geno_lookup <- geno_res[geno_res$gene %in% restored$gene, c("gene", "log2FoldChange")]
  colnames(geno_lookup) <- c("gene", "log2FC_genotype")
  
  merged <- merge(restored[, c("gene", "log2FoldChange")], geno_lookup, by = "gene")
  colnames(merged)[2] <- "log2FC_restoration"
  
  merged$pct_shrinkage <- round(100 * (1 - abs(merged$log2FC_restoration) / abs(merged$log2FC_genotype)), 1)
  
  median_shrinkage <- round(median(merged$pct_shrinkage, na.rm = TRUE), 1)
  n_low <- sum(merged$pct_shrinkage < 20, na.rm = TRUE)
  
  shrinkage_summary_table <- rbind(shrinkage_summary_table, data.frame(
    cell_type = ct, n_restored = nrow(merged),
    median_pct_shrinkage = median_shrinkage,
    n_low_shrinkage = n_low
  ))
}

shrinkage_summary_table <- shrinkage_summary_table[order(-shrinkage_summary_table$n_restored), ]

print(shrinkage_summary_table)

write.csv(shrinkage_summary_table, file.path(PROJECT_ROOT, "Analysis/shrinkage_summary_all_celltypes.csv"), row.names = FALSE)