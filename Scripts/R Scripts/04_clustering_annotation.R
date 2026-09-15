# ============================================================
# CLUSTERING AND CELL TYPE ANNOTATION
#
# Clusters the integrated dataset, assigns provisional cell-type
# labels via marker-based module scoring, and performs full
# cluster-level quality control (marker specificity, doublet
# co-expression confirmation, cluster purity) to identify and
# exclude technically compromised clusters.
#
# Requires:
#   - merged_seurat_post_umap.rds, produced by
#     03_dimensionality_reduction.R
#     (PROJECT_ROOT/Analysis/merged_seurat_post_umap.rds)
#
# Produces:
#   - merged_seurat_post_clustering.rds
#   - merged_seurat_post_annotation.rds
#   - cluster_level_annotation_table.csv
#   - cluster_module_scores_all.csv
#   - cluster_size_distribution.csv
#   - cluster_marker_specificity_summary.csv
#   - cluster_top_named_marker.csv
#   - cluster_QC_summary.csv
#   - cluster_purity_summary.csv
#   - Final checkpoint: merged_seurat with Low_quality clusters
#     (9, 19, and confirmed doublet clusters) excluded, ready for
#     pseudobulk differential expression analysis
# ============================================================

# ============================================================
# EDIT THIS: set to your own project root before running
# ============================================================
PROJECT_ROOT <- "/scratch/SCWF00151/shared/c.d21085818/Dissertation"

library(Seurat)
library(ggplot2)

# Load the post-UMAP object (produced by 03_dimensionality_reduction.R)
merged_seurat <- readRDS(file.path(PROJECT_ROOT, "Analysis/merged_seurat_post_umap.rds"))

# ============================================================
# CLUSTERING
# ============================================================
# UPDATED: random.seed explicitly set on both calls, consistent with
# every other stochastic step in this pipeline (SCTransform, RunPCA,
# RunHarmony, RunUMAP). Without this, FindClusters falls back to its
# own internal default seed regardless of any earlier set.seed() call.
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:20,
                               random.seed = 123)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5,
                              random.seed = 123)
DimPlot(merged_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Cluster Assignments (UMAP)")
saveRDS(merged_seurat, file.path(PROJECT_ROOT, "Analysis/merged_seurat_post_clustering.rds"))
# ============================================================
# CELL TYPE ANNOTATION
# ============================================================
# Define marker panel (original, Section 2.10)
celltype_markers <- list(
  MSN = c("Ppp1r1b","Bcl11b","Pde1b","Drd2","Penk","Grik3","Ttc12","Gpr6","Adora2a","Drd1","Tac1","Pdyn","Ebf1","Slc35d3","Foxp2", "Cacng5", "Olfm3","Dcx"),
  Astrocytes = c("Slc1a2", "Gfap", "Slc6a11", "Grm3"),
  Microglia = c("Csf1r", "Trem2", "Gpr34", "Gal3st4"),
  OPC = c("Cspg4", "Gpr17", "Neu4", "Zfp488", "Olig2", "A930009A15Rik", "Olig1", "Rlbp1"),
  Oligodendrocytes = c("Plp1", "Mobp", "Prr5l", "Cdh20", "Mag", "Mog"),
  Cholinergic_Interneuron = "Chat",
  PV_Th_Interneuron = c("Hs3st2", "Cntnap4"),
  Sst_Npy_Interneuron = "Nos1",
  Endothelial = c("Atp13a5", "Gm5127", "Fbln5", "Zic3", "Ltbp4", "Slc2a1"),
  Mural = c("Car3", "Mylk", "Cald1", "Flna", "Lpl", "Pdgfrb"),
  Ciliated_Ependymal = c("Hydin", "Scgb1a1", "Sec14l3", "Cdc20b", "Cyp2f2"),
  Secretory_Ependymal = "Ttr",
  Cycling_cell = "Top2a"
)
# Validation marker panel - independent set used to cross-check and
# strengthen the original panel, added to increase robustness against
# ambient RNA contamination (Section 3.x)
proposed_new_markers <- list(
  MSN = c("Rgs9", "Meis2", "Rarb", "Camk4", "Arpp21", "Isl1"),
  Astrocytes = c("Aqp4", "Slc1a3", "Aldh1l1", "Gja1", "Sox9"),
  Microglia = c("P2ry12", "Cx3cr1", "C1qa", "C1qb", "Tmem119", "Hexb"),
  OPC = c("Pdgfra", "Sox10", "Vcan"),
  Oligodendrocytes = c("Sox10", "Cnp", "Opalin", "Aspa"),
  Endothelial = c("Cldn5", "Pecam1", "Flt1"),
  Ciliated_Ependymal = c("Foxj1", "Ccdc153", "Tmem212"),
  Cycling_cell = c("Mki67", "Pcna", "Ccnb2")
)
# Merge the two panels: categories in both are combined and deduplicated;
# categories only in the original panel (Cholinergic_Interneuron,
# PV_Th_Interneuron, Sst_Npy_Interneuron, Mural, Secretory_Ependymal)
# are carried through unchanged, since the validation panel doesn't cover them
all_categories <- union(names(celltype_markers), names(proposed_new_markers))
celltype_markers <- setNames(
  lapply(all_categories, function(ct) {
    unique(c(celltype_markers[[ct]], proposed_new_markers[[ct]]))
  }),
  all_categories
)
# Confirm the merge - inspect gene counts per category, and check for
# any silent duplication from case/typo mismatches between the two lists
cat("Marker counts per category after merge:\n")
for (ct in all_categories) print(paste0(ct, ": ", length(celltype_markers[[ct]]), " genes"))
print(celltype_markers)
# Score each cluster - seed set explicitly to 123 for reproducibility
DefaultAssay(merged_seurat) <- "SCT"
merged_seurat <- AddModuleScore(merged_seurat, features = celltype_markers, name = "celltype_score", seed = 123)
score_cols <- paste0("celltype_score", 1:length(celltype_markers))
cluster_scores <- aggregate(merged_seurat@meta.data[, score_cols],
                            by = list(cluster = merged_seurat$seurat_clusters), FUN = mean)
colnames(cluster_scores) <- c("cluster", names(celltype_markers))
cluster_scores$best_match <- apply(cluster_scores[, names(celltype_markers)], 1, function(x) names(celltype_markers)[which.max(x)])
print(cluster_scores[, c("cluster", "best_match")])
# Apply provisional cell_type labels from best_match, for a first look
provisional_labels <- setNames(cluster_scores$best_match, cluster_scores$cluster)
merged_seurat$cell_type_provisional <- unname(provisional_labels[as.character(merged_seurat$seurat_clusters)])
# Flip axis 1 to match Murillo et al.'s published UMAP orientation
merged_seurat$umap_flip_1 <- -Embeddings(merged_seurat, "umap")[,1]
merged_seurat[["umap"]]@cell.embeddings[,1] <- merged_seurat$umap_flip_1
library(ggplot2)
DimPlot(merged_seurat, reduction = "umap", group.by = "cell_type_provisional", 
        label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle("Cluster Annotation By Cell Type")
# Print cell numbers per cell type to console.
cat("Cell counts per cluster (post-annotation):\n")
print(table(merged_seurat$cell_type_provisional))
saveRDS(merged_seurat, file.path(PROJECT_ROOT, "Analysis/merged_seurat_post_annotation.rds"))
# ============================================================
# CLUSTER-LEVEL TABLE: ALL 25 CLUSTERS, CELL COUNT, AND
# PROVISIONAL CELL TYPE LABEL
# ============================================================
# Using merged_seurat (the post-annotation, pre-cluster-QC
# checkpoint) so all 25 original clusters and their initial
# provisional labels are visible, before any exclusions
cluster_table <- as.data.frame(table(merged_seurat$seurat_clusters,
                                     merged_seurat$cell_type_provisional))
colnames(cluster_table) <- c("cluster", "cell_type", "n_cells")
# Keep only non-zero rows (each cluster has exactly one assigned type)
cluster_table <- cluster_table[cluster_table$n_cells > 0, ]
cluster_table <- cluster_table[order(as.numeric(as.character(cluster_table$cluster))), ]
cat("=== All", nrow(cluster_table), "clusters, with cell count and provisional label ===\n")
print(cluster_table)
# Confirm this sums to your total pre-cluster-QC cell count
cat("\nTotal cells across all clusters:", sum(cluster_table$n_cells), "\n")
write.csv(cluster_table, file.path(PROJECT_ROOT, "Analysis/cluster_level_annotation_table.csv"), row.names = FALSE)
# ============================================================
# MODULE SCORES FOR ALL 25 CLUSTERS, ACROSS ALL CELL TYPE
# CATEGORIES, PLUS BEST-MATCH LABEL
# ============================================================
cluster_scores_full <- cluster_scores[order(as.numeric(as.character(cluster_scores$cluster))), ]
cat("=== Module scores, all 25 clusters, all categories ===\n")
print(cluster_scores_full)
write.csv(cluster_scores_full, file.path(PROJECT_ROOT, "Analysis/cluster_module_scores_all.csv"), row.names = FALSE)
# ============================================================
# CLUSTER SIZE DISTRIBUTION
#
# Basic sanity check on the clustering output before any further
# analysis: confirms the number of clusters produced, and the
# distribution of cell counts across them. Very small clusters
# (few cells) are flagged here as a first-pass signal for
# tracking into later QC steps (doublet detection, annotation
# confidence), though size alone is not sufficient evidence to
# act on at this stage.
# ============================================================
cat("Number of clusters:", length(levels(merged_seurat$seurat_clusters)), "\n\n")
cluster_sizes <- table(merged_seurat$seurat_clusters)
cluster_sizes_df <- data.frame(
  cluster = names(cluster_sizes),
  n_cells = as.integer(cluster_sizes)
)
cluster_sizes_df <- cluster_sizes_df[order(-cluster_sizes_df$n_cells), ]
cat("=== Cluster sizes, largest to smallest ===\n")
print(cluster_sizes_df)
cat("\nLargest cluster:", cluster_sizes_df$cluster[1], "(", cluster_sizes_df$n_cells[1], "cells )\n")
cat("Smallest cluster:", cluster_sizes_df$cluster[nrow(cluster_sizes_df)],
    "(", cluster_sizes_df$n_cells[nrow(cluster_sizes_df)], "cells )\n")
# Flag clusters under a minimum size threshold - not acted on here,
# just noted for awareness going into later QC steps
MIN_CLUSTER_SIZE <- 100
small_clusters <- cluster_sizes_df$cluster[cluster_sizes_df$n_cells < MIN_CLUSTER_SIZE]
cat("\nClusters below", MIN_CLUSTER_SIZE, "cells (flagged for awareness only):",
    paste(small_clusters, collapse = ", "), "\n")
# Visual: cluster sizes as a bar chart, ordered largest to smallest
library(ggplot2)
ggplot(cluster_sizes_df, aes(x = reorder(cluster, -n_cells), y = n_cells)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = MIN_CLUSTER_SIZE, linetype = "dashed", color = "red") +
  labs(title = "Cluster size distribution", x = "Cluster", y = "Number of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
write.csv(cluster_sizes_df, file.path(PROJECT_ROOT, "Analysis/cluster_size_distribution.csv"), row.names = FALSE)
# ============================================================
# MARKER GENE SPECIFICITY CHECK (all clusters)
#
# For each cluster, computes its top unbiased differentially
# expressed genes (FindMarkers, not drawn from any pre-selected
# marker panel) and evaluates how SPECIFIC its top marker is:
# a good marker should be expressed in most of the cluster's own
# cells (high pct.1) but rarely elsewhere (low pct.2), with a
# strong fold-change.
# ============================================================
DefaultAssay(merged_seurat) <- "SCT"
Idents(merged_seurat) <- "seurat_clusters"
all_cluster_markers <- list()
marker_specificity_summary <- data.frame(
  cluster = character(), n_cells = integer(),
  top_marker = character(), top_marker_log2FC = numeric(),
  top_marker_pct1 = numeric(), top_marker_pct2 = numeric(),
  specificity_gap = numeric(),
  stringsAsFactors = FALSE
)
for (cl in levels(merged_seurat$seurat_clusters)) {
  cat("Processing cluster:", cl, "\n")
  
  markers <- FindMarkers(
    merged_seurat,
    ident.1 = cl,
    group.by = "seurat_clusters",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5
  )
  markers <- markers[order(-markers$avg_log2FC), ]
  all_cluster_markers[[cl]] <- markers
  
  if (nrow(markers) == 0) {
    cat("  No markers passed thresholds for cluster", cl, "\n")
    next
  }
  
  top <- markers[1, ]
  marker_specificity_summary <- rbind(marker_specificity_summary, data.frame(
    cluster = cl,
    n_cells = sum(merged_seurat$seurat_clusters == cl),
    top_marker = rownames(top),
    top_marker_log2FC = round(top$avg_log2FC, 2),
    top_marker_pct1 = top$pct.1,
    top_marker_pct2 = top$pct.2,
    specificity_gap = round(top$pct.1 - top$pct.2, 3)
  ))
}
# Sort by specificity_gap ascending - clusters with the WEAKEST
# specificity (top marker not cleanly distinguishing the cluster
# from the rest of the dataset) appear first
marker_specificity_summary <- marker_specificity_summary[order(marker_specificity_summary$specificity_gap), ]
cat("\n=== Marker specificity summary, all clusters (weakest specificity first) ===\n")
print(marker_specificity_summary)
write.csv(marker_specificity_summary, 
          file.path(PROJECT_ROOT, "Analysis/cluster_marker_specificity_summary.csv"), 
          row.names = FALSE)
# ============================================================
# FOR EACH CLUSTER, FIND THE TOP-RANKED KNOWN (NAMED) GENE,
# SKIPPING Gm-PREFIXED / NUMERIC-PREFIXED / BC-PREFIXED /
# Rik-SUFFIXED PREDICTED GENES
# ============================================================
is_named_gene <- function(genes) {
  !grepl("^Gm[0-9]", genes) & !grepl("^[0-9]", genes) &
    !grepl("Rik$", genes) & !grepl("^BC[0-9]", genes) &
    !grepl("^ENSMUSG", genes)
}
top_named_marker_summary <- data.frame(cluster = character(), cell_type = character(),
                                       top_named_marker = character(), rank = integer(),
                                       avg_log2FC = numeric(), pct.1 = numeric(),
                                       pct.2 = numeric(), specificity_gap = numeric(),
                                       stringsAsFactors = FALSE)
for (cl in names(all_cluster_markers)) {
  markers <- all_cluster_markers[[cl]]
  if (nrow(markers) == 0) next
  
  markers$gene <- rownames(markers)
  named_markers <- markers[is_named_gene(markers$gene), ]
  
  if (nrow(named_markers) == 0) {
    cat("Cluster", cl, ": no named genes found in marker list\n")
    next
  }
  
  top_named <- named_markers[1, ]
  original_rank <- which(rownames(markers) == rownames(top_named))
  
  ct <- cluster_table$cell_type[cluster_table$cluster == cl]
  
  top_named_marker_summary <- rbind(top_named_marker_summary, data.frame(
    cluster = cl, cell_type = ifelse(length(ct) > 0, ct, NA),
    top_named_marker = rownames(top_named), rank = original_rank,
    avg_log2FC = round(top_named$avg_log2FC, 3),
    pct.1 = top_named$pct.1, pct.2 = top_named$pct.2,
    specificity_gap = round(top_named$pct.1 - top_named$pct.2, 3)
  ))
}
top_named_marker_summary <- top_named_marker_summary[order(as.numeric(as.character(top_named_marker_summary$cluster))), ]
cat("=== Top-ranked NAMED gene per cluster (rank shown = position in original ranking) ===\n")
print(top_named_marker_summary)
write.csv(top_named_marker_summary, file.path(PROJECT_ROOT, "Analysis/cluster_top_named_marker.csv"), row.names = FALSE)
# ============================================================
# Cluster 0 and 9: top 15 markers, full specificity picture
# (not just the single top gene) - checking whether weak
# specificity holds throughout, or was just an artifact of
# ranking by log2FC alone
# ============================================================
for (cl in c("0", "9")) {
  cat("\n=== Cluster", cl, ": top 15 markers, full specificity detail ===\n")
  markers <- all_cluster_markers[[cl]]
  top15 <- head(markers, 15)
  top15$specificity_gap <- round(top15$pct.1 - top15$pct.2, 3)
  print(top15[, c("avg_log2FC", "pct.1", "pct.2", "specificity_gap", "p_val_adj")])
}
# Cluster 9 percent.mt distribution vs rest of dataset
cat("Cluster 9 percent.mt summary:\n")
print(summary(merged_seurat$percent.mt[merged_seurat$seurat_clusters == "9"]))
cat("\nRest of dataset percent.mt summary:\n")
print(summary(merged_seurat$percent.mt[merged_seurat$seurat_clusters != "9"]))
VlnPlot(merged_seurat, features = "percent.mt", group.by = "seurat_clusters", pt.size = 0) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
  ggtitle("percent.mt per cluster (red = your QC threshold, 5%)")
# ============================================================
# CLUSTER 9: quality exclusion
#
# Cluster 9 (3,150 cells, previously labeled MSN) shows a
# mitochondrial signature markedly elevated relative to the rest
# of the dataset (median percent.mt 1.57% vs 0.25% dataset-wide,
# ~6.3x higher), and its top unbiased markers are dominated by
# mitochondrial and ribosomal/housekeeping genes rather than any
# specific cell-type signature - consistent with cellular stress
# or lower nuclear RNA quality rather than a distinct population.
# All cells individually pass the percent.mt < 5% QC threshold,
# so this pattern was only detectable at the
# cluster level, post-clustering.
# ============================================================
cat("Cells being reassigned from cluster 9:", sum(merged_seurat$seurat_clusters == "9"), "\n")
merged_seurat$cell_type_provisional[merged_seurat$seurat_clusters == "9"] <- "Low_quality"
cat("\nUpdated cell_type_provisional counts:\n")
print(table(merged_seurat$cell_type_provisional))
# ============================================================
# 04_cluster_QC_and_doublets.R
#
# Full, reproducible cluster QC pipeline: every diagnostic step
# that identifies suspicious clusters is included and printed in
# full for ALL clusters - nothing is pre-filtered or hidden.
# Confirmation (co-expression) and correction happen only after
# the full diagnostic picture is visible.
# ============================================================
# ============================================================
# SECTION 1: Cluster size distribution (all clusters)
# ============================================================
cat("=== Cluster size distribution ===\n")
print(table(merged_seurat$seurat_clusters))
# ============================================================
# SECTION 2: Doublet score per cluster (all clusters)
# ============================================================
doublet_by_cluster <- aggregate(
  doublet_score ~ seurat_clusters,
  data = merged_seurat@meta.data,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE))
)
doublet_by_cluster <- do.call(data.frame, doublet_by_cluster)
colnames(doublet_by_cluster) <- c("cluster", "mean_doublet_score", "median_doublet_score")
doublet_by_cluster <- doublet_by_cluster[order(-doublet_by_cluster$mean_doublet_score), ]
cat("\n=== Doublet score per cluster, all clusters, sorted highest first ===\n")
print(doublet_by_cluster)
VlnPlot(merged_seurat, features = "doublet_score", group.by = "seurat_clusters", pt.size = 0) +
  geom_hline(yintercept = median(merged_seurat$doublet_score, na.rm = TRUE), linetype = "dashed", color = "red") +
  ggtitle("Doublet score per cluster (red = dataset median)")
# ============================================================
# SECTION 3: nCount_RNA / nFeature_RNA per cluster (all clusters)
# ============================================================
qc_by_cluster <- aggregate(
  cbind(nCount_RNA, nFeature_RNA) ~ seurat_clusters,
  data = merged_seurat@meta.data,
  FUN = function(x) c(mean = mean(x), median = median(x))
)
qc_by_cluster <- do.call(data.frame, qc_by_cluster)
colnames(qc_by_cluster) <- c("cluster", "mean_nCount", "median_nCount", "mean_nFeature", "median_nFeature")
qc_by_cluster <- qc_by_cluster[order(-qc_by_cluster$mean_nCount), ]
cat("\n=== nCount/nFeature per cluster, all clusters, sorted highest nCount first ===\n")
print(qc_by_cluster)
VlnPlot(merged_seurat, features = "nCount_RNA", group.by = "seurat_clusters", pt.size = 0) +
  geom_hline(yintercept = median(merged_seurat$nCount_RNA), linetype = "dashed", color = "red") +
  ggtitle("nCount_RNA per cluster (red = dataset median)")
VlnPlot(merged_seurat, features = "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0) +
  geom_hline(yintercept = median(merged_seurat$nFeature_RNA), linetype = "dashed", color = "red") +
  ggtitle("nFeature_RNA per cluster (red = dataset median)")
# ============================================================
# SECTION 4: Cell type annotation panels, rebuilt for margin check
# ============================================================
celltype_markers_original <- list(
  MSN = c("Ppp1r1b","Bcl11b","Pde1b","Drd2","Penk","Grik3","Ttc12","Gpr6","Adora2a","Drd1","Tac1","Pdyn","Ebf1","Slc35d3","Foxp2","Cacng5","Olfm3","Dcx"),
  Astrocytes = c("Slc1a2","Gfap","Slc6a11","Grm3"),
  Microglia = c("Csf1r","Trem2","Gpr34","Gal3st4"),
  OPC = c("Cspg4","Gpr17","Neu4","Zfp488","Olig2","A930009A15Rik","Olig1","Rlbp1"),
  Oligodendrocytes = c("Plp1","Mobp","Prr5l","Cdh20","Mag","Mog"),
  Cholinergic_Interneuron = "Chat",
  PV_Th_Interneuron = c("Hs3st2","Cntnap4"),
  Sst_Npy_Interneuron = "Nos1",
  Endothelial = c("Atp13a5","Gm5127","Fbln5","Zic3","Ltbp4","Slc2a1"),
  Mural = c("Car3","Mylk","Cald1","Flna","Lpl","Pdgfrb"),
  Ciliated_Ependymal = c("Hydin","Scgb1a1","Sec14l3","Cdc20b","Cyp2f2"),
  Secretory_Ependymal = "Ttr",
  Cycling_cell = "Top2a"
)
proposed_new_markers <- list(
  MSN = c("Rgs9","Meis2","Rarb","Camk4","Arpp21","Isl1"),
  Astrocytes = c("Aqp4","Slc1a3","Aldh1l1","Gja1","Sox9"),
  Microglia = c("P2ry12","Cx3cr1","C1qa","C1qb","Tmem119","Hexb"),
  OPC = c("Pdgfra","Sox10","Vcan"),
  Oligodendrocytes = c("Sox10","Cnp","Opalin","Aspa"),
  Endothelial = c("Cldn5","Pecam1","Flt1"),
  Ciliated_Ependymal = c("Foxj1","Ccdc153","Tmem212"),
  Cycling_cell = c("Mki67","Pcna","Ccnb2")
)
all_categories <- union(names(celltype_markers_original), names(proposed_new_markers))
celltype_markers <- setNames(
  lapply(all_categories, function(ct) unique(c(celltype_markers_original[[ct]], proposed_new_markers[[ct]]))),
  all_categories
)
score_cols <- colnames(merged_seurat@meta.data)[grep("^celltype_score", colnames(merged_seurat@meta.data))]
stopifnot(length(score_cols) == length(celltype_markers))
cluster_scores <- aggregate(merged_seurat@meta.data[, score_cols],
                            by = list(cluster = merged_seurat$seurat_clusters), FUN = mean)
colnames(cluster_scores) <- c("cluster", names(celltype_markers))
cluster_scores$best_match <- apply(cluster_scores[, names(celltype_markers)], 1, function(x) names(celltype_markers)[which.max(x)])
cluster_scores$top_score <- apply(cluster_scores[, names(celltype_markers)], 1, max)
cluster_scores$second_score <- apply(cluster_scores[, names(celltype_markers)], 1, function(x) sort(x, decreasing = TRUE)[2])
cluster_scores$second_best_match <- apply(cluster_scores[, names(celltype_markers)], 1, function(x) names(celltype_markers)[order(-x)[2]])
cluster_scores$margin <- cluster_scores$top_score - cluster_scores$second_score
cluster_scores_sorted <- cluster_scores[order(cluster_scores$margin), ]
cat("\n=== Annotation confidence (margin), all clusters, sorted most ambiguous first ===\n")
print(cluster_scores_sorted[, c("cluster", "best_match", "second_best_match", "top_score", "second_score", "margin")])
# ============================================================
# SECTION 5: Consolidated QC summary table (combines Sections
# 1-4 into one master reference table, all clusters)
# ============================================================
qc_summary <- merge(qc_by_cluster[, c("cluster", "median_nCount", "median_nFeature")],
                    doublet_by_cluster[, c("cluster", "mean_doublet_score", "median_doublet_score")],
                    by = "cluster")
qc_summary <- merge(qc_summary, cluster_scores[, c("cluster", "best_match", "second_best_match", "margin")], by = "cluster")
qc_summary$n_cells <- as.integer(table(merged_seurat$seurat_clusters)[qc_summary$cluster])
qc_summary <- qc_summary[order(-qc_summary$mean_doublet_score), ]
cat("\n=== MASTER QC SUMMARY: all clusters, all metrics ===\n")
print(qc_summary)
write.csv(qc_summary, file.path(PROJECT_ROOT, "Analysis/cluster_QC_summary.csv"), row.names = FALSE)
# ============================================================
# SECTION 6: Marker co-expression confirmation
#
# Clusters tested here (3, 23, 24, 26, 27) were identified as
# suspicious:
#   - 23, 24: elevated doublet score AND elevated UMI/nFeature
#   - 26: close annotation margin (Oligo vs Microglia, same
#     pairing as 23), flagged by Section 4
#   - 3: single largest UMI outlier in the entire dataset
#     (Section 3), despite unremarkable doublet score
#   - 27: close annotation margin (Ciliated_Ependymal vs
#     Astrocytes), flagged by Section 4
# ============================================================
flagged_clusters <- c("3", "23", "24", "26", "27")
# ------------------------------------------------------------
# SECTION 6a: Unbiased top markers for the flagged clusters only
# ------------------------------------------------------------
DefaultAssay(merged_seurat) <- "SCT"
Idents(merged_seurat) <- "seurat_clusters"
all_cluster_markers <- list()
for (cl in flagged_clusters) {
  cat("Computing unbiased markers for cluster:", cl, "\n")
  markers <- FindMarkers(
    merged_seurat,
    ident.1 = cl,
    group.by = "seurat_clusters",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5
  )
  markers <- markers[order(-markers$avg_log2FC), ]
  all_cluster_markers[[cl]] <- markers
}
# ------------------------------------------------------------
# SECTION 6b: Co-expression check
#
# Explicitly sets DefaultAssay to SCT and filters each marker
# panel down to genes actually present in the object before
# testing, printing anything dropped rather than failing or
# silently mismatching (resolves the earlier "undefined columns
# selected" error, caused by a marker gene e.g. Scgb1a1 not being
# found in the expected assay/layer).
# ------------------------------------------------------------
check_coexpression_auto <- function(seurat_obj, cluster_id, cluster_scores_tbl, marker_list) {
  DefaultAssay(seurat_obj) <- "SCT"
  
  row <- cluster_scores_tbl[cluster_scores_tbl$cluster == cluster_id, ]
  type_A <- row$best_match
  type_B <- row$second_best_match
  markers_A_raw <- marker_list[[type_A]]
  markers_B_raw <- marker_list[[type_B]]
  
  shared_genes <- intersect(markers_A_raw, markers_B_raw)
  markers_A <- setdiff(markers_A_raw, shared_genes)
  markers_B <- setdiff(markers_B_raw, shared_genes)
  
  available_genes <- rownames(seurat_obj)
  missing_A <- setdiff(markers_A, available_genes)
  missing_B <- setdiff(markers_B, available_genes)
  if (length(missing_A) > 0) cat("Cluster", cluster_id, "-", type_A, "genes not in object, dropped:", paste(missing_A, collapse=", "), "\n")
  if (length(missing_B) > 0) cat("Cluster", cluster_id, "-", type_B, "genes not in object, dropped:", paste(missing_B, collapse=", "), "\n")
  markers_A <- intersect(markers_A, available_genes)
  markers_B <- intersect(markers_B, available_genes)
  
  if (length(shared_genes) > 0) {
    cat("Cluster", cluster_id, "- excluded shared genes (", type_A, "vs", type_B, "):",
        paste(shared_genes, collapse = ", "), "\n")
  }
  
  if (length(markers_A) == 0 | length(markers_B) == 0) {
    cat("Cluster", cluster_id, "- SKIPPED: one panel empty after filtering\n")
    return(data.frame(cluster = cluster_id, n_cells = NA, type_A = type_A, type_B = type_B,
                      pct_both = NA, pct_A_only = NA, pct_B_only = NA, pct_neither = NA))
  }
  
  cells <- subset(seurat_obj, subset = seurat_clusters == cluster_id)
  expr <- FetchData(cells, vars = c(markers_A, markers_B))
  stopifnot(all(c(markers_A, markers_B) %in% colnames(expr)))
  
  pos_A <- rowSums(expr[, markers_A, drop = FALSE] > 0) > 0
  pos_B <- rowSums(expr[, markers_B, drop = FALSE] > 0) > 0
  
  data.frame(
    cluster = cluster_id, n_cells = nrow(expr), type_A = type_A, type_B = type_B,
    pct_both = round(100 * mean(pos_A & pos_B), 1),
    pct_A_only = round(100 * mean(pos_A & !pos_B), 1),
    pct_B_only = round(100 * mean(!pos_A & pos_B), 1),
    pct_neither = round(100 * mean(!pos_A & !pos_B), 1)
  )
}
coexpression_results <- do.call(rbind, lapply(flagged_clusters, function(cl) {
  check_coexpression_auto(merged_seurat, cl, cluster_scores, celltype_markers)
}))
cat("\n=== Co-expression results, flagged clusters ===\n")
print(coexpression_results)
for (cl in flagged_clusters) {
  cat("\n--- Cluster", cl, ": top 10 unbiased markers ---\n")
  print(head(all_cluster_markers[[cl]], 10))
}
# ============================================================
# SECTION 7: Apply correction
# Threshold: >20% both-positive treated as doublet-confirmed
# (consistent with clusters 23 [88.1%] and 26 [60.6%]).
# ============================================================
DOUBLET_THRESHOLD_PCT <- 20
confirmed_doublet_clusters <- coexpression_results$cluster[!is.na(coexpression_results$pct_both) &
                                                             coexpression_results$pct_both > DOUBLET_THRESHOLD_PCT]
cat("\nClusters confirmed as doublet-dominated (>", DOUBLET_THRESHOLD_PCT, "% both-positive):",
    paste(confirmed_doublet_clusters, collapse = ", "), "\n")
merged_seurat$cell_type_provisional[merged_seurat$seurat_clusters %in% confirmed_doublet_clusters] <- "Low_quality"
cat("\nUpdated cell_type_provisional counts:\n")
print(table(merged_seurat$cell_type_provisional))
DimPlot(merged_seurat, reduction = "umap", group.by = "cell_type_provisional",
        label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle("Cell type annotation (confirmed doublet clusters excluded)")
# ============================================================
# SECTION 8: CLUSTER PURITY / CONTAMINATION SCORE
#
# Per-cell, within-cluster check - the most direct answer to "is
# this cluster contaminated with genes from other cell types."
# Unlike the annotation margin, which compares
# CLUSTER-LEVEL AVERAGE scores, this checks each INDIVIDUAL cell's
# own highest-scoring cell type, then asks: what fraction of cells
# in this cluster individually agree with the cluster's assigned
# label, vs. how many individually look more like a different
# cell type? A cluster can pass the margin check (clear average
# winner) while still containing a substantial contaminated
# minority - this check is what would catch that.
# ============================================================
score_cols <- colnames(merged_seurat@meta.data)[grep("^celltype_score", colnames(merged_seurat@meta.data))]
stopifnot(length(score_cols) == length(celltype_markers))
# For every single cell, determine which cell type it individually
# scores highest for (per-cell, not per-cluster average)
cell_scores_matrix <- merged_seurat@meta.data[, score_cols]
colnames(cell_scores_matrix) <- names(celltype_markers)
merged_seurat$individual_cell_best_match <- names(celltype_markers)[apply(cell_scores_matrix, 1, which.max)]
# ------------------------------------------------------------
# For each cluster: what % of its cells individually agree with
# the cluster's assigned cell_type_provisional label, and what's
# the largest single "contaminating" cell type among the rest?
# ------------------------------------------------------------
purity_summary <- data.frame(
  cluster = character(), n_cells = integer(),
  assigned_label = character(),
  pct_cells_matching_assigned_label = numeric(),
  top_contaminating_type = character(),
  pct_top_contaminant = numeric(),
  stringsAsFactors = FALSE
)
for (cl in levels(merged_seurat$seurat_clusters)) {
  cells_in_cluster <- merged_seurat$seurat_clusters == cl
  n <- sum(cells_in_cluster)
  if (n == 0) next
  
  assigned <- unique(merged_seurat$cell_type_provisional[cells_in_cluster])
  if (length(assigned) != 1) next  # skip if a cluster somehow has mixed labels
  
  individual_matches <- merged_seurat$individual_cell_best_match[cells_in_cluster]
  match_table <- sort(table(individual_matches), decreasing = TRUE)
  
  pct_matching <- round(100 * sum(individual_matches == assigned) / n, 1)
  
  non_assigned <- match_table[names(match_table) != assigned]
  top_contaminant <- if (length(non_assigned) > 0) names(non_assigned)[1] else NA
  pct_contaminant <- if (length(non_assigned) > 0) round(100 * non_assigned[1] / n, 1) else 0
  
  purity_summary <- rbind(purity_summary, data.frame(
    cluster = cl, n_cells = n, assigned_label = assigned,
    pct_cells_matching_assigned_label = pct_matching,
    top_contaminating_type = top_contaminant,
    pct_top_contaminant = pct_contaminant
  ))
}
purity_summary <- purity_summary[order(purity_summary$pct_cells_matching_assigned_label), ]
cat("=== Cluster purity / contamination summary, sorted LOWEST purity first ===\n")
print(purity_summary)
write.csv(purity_summary, file.path(PROJECT_ROOT, "Analysis/cluster_purity_summary.csv"), row.names = FALSE)
# ------------------------------------------------------------
# Visual: violin/distribution of the ASSIGNED type's own score,
# per cluster - shows the spread directly, not just the summary %
# ------------------------------------------------------------
for (ct in unique(purity_summary$assigned_label)) {
  clusters_of_this_type <- purity_summary$cluster[purity_summary$assigned_label == ct]
  score_col_index <- which(names(celltype_markers) == ct)
  if (length(score_col_index) == 0) next
  
  score_col_name <- score_cols[score_col_index]
  
  p <- VlnPlot(merged_seurat,
               features = score_col_name,
               idents = clusters_of_this_type,
               group.by = "seurat_clusters",
               pt.size = 0) +
    ggtitle(paste0(ct, " module score, own clusters (", paste(clusters_of_this_type, collapse=", "), ")"))
  print(p)
}
# Distribution of Cycling_cell score vs Oligodendrocyte score,
# individually, across cluster 19's cells - tests whether a small
# subset of high-Cycling_cell cells is pulling the cluster mean up
cluster19_scores <- merged_seurat@meta.data[merged_seurat$seurat_clusters == "19", score_cols]
colnames(cluster19_scores) <- names(celltype_markers)
cat("Cycling_cell score distribution, cluster 19:\n")
print(summary(cluster19_scores$Cycling_cell))
cat("\nOligodendrocytes score distribution, cluster 19:\n")
print(summary(cluster19_scores$Oligodendrocytes))
# How many cells have a NEAR-ZERO Cycling_cell score (i.e. the
# cluster mean is being driven by a minority of high scorers)?
cat("\n% of cluster 19 cells with Cycling_cell score below 0.05:",
    round(100 * mean(cluster19_scores$Cycling_cell < 0.05), 1), "%\n")
# ============================================================
# CLUSTER 19: quality exclusion
#
# Cluster 19 (916 cells, previously labeled Cycling_cell) shows:
#   - The Cycling_cell label was driven by a small subset of
#     outlier cells with strong scores pulling the cluster MEAN
#     upward, despite 87.7% of cells individually showing a
#     near-zero or negative Cycling_cell score (median -0.018)
#   - Only 9.1% individual-cell purity against its assigned label
#     (Section 8 purity check) - the lowest of any cluster in the
#     dataset

# ============================================================
cat("Cells being reassigned from cluster 19:", sum(merged_seurat$seurat_clusters == "19"), "\n")
merged_seurat$cell_type_provisional[merged_seurat$seurat_clusters == "19"] <- "Low_quality"
cat("\nUpdated cell_type_provisional counts:\n")
print(table(merged_seurat$cell_type_provisional))