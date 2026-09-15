# ============================================================
# STEP 10: Principal component analysis (Section 2.9)
# ============================================================
# npcs = 50: deliberately computes more components than will
# ultimately be used, to allow inspection of the variance-explained
# curve and confirm the number of components retained downstream
# (20, matching Murillo et al., 2025) is appropriate for this
# dataset, rather than assumed
set.seed(123)
merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = TRUE)
# Visualise variance explained per component to identify the "elbow".
# Confirmed a clear inflection between components 6-8, with 20
# components sitting comfortably beyond this point (Section 3.6)
ElbowPlot(merged_seurat, ndims = 50) +
  ggtitle("Variance Explained by Principal Component")

# ============================================================
# STEP 11: Harmony batch integration (Section 2.9)
# ============================================================
library(harmony)
# Harmony's iterative optimisation involves stochastic elements;
# seed set immediately before the call for reproducibility
set.seed(123)
# group.by.vars = "orig.ident": corrects for batch effects between
#   the 12 individual samples
# dims.use = 1:20: integrates using the 20 principal components
#   retained following elbow plot inspection
# theta = 1: diversity clustering penalty, matching the integration
#   parameter used by Murillo et al. (2025)
merged_seurat <- RunHarmony(
  merged_seurat,
  group.by.vars = "orig.ident",
  dims.use = 1:20,
  theta = 1
)
# Confirm the Harmony-corrected reduction was added successfully
Reductions(merged_seurat)
# Confirm dimensions: expect 102,575 cells x 20 components
dim(Embeddings(merged_seurat, reduction = "harmony"))
# ============================================================
# STEP 12: UMAP (Section 2.9)
# ============================================================
# UMAP's layout algorithm is explicitly randomised; seed set
# immediately before the call for reproducibility
library(ggplot2)
set.seed(123)
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:20)
# Visual check: cells from different samples should intermix within
# shared regions of the embedding, confirming successful integration
# rather than sample-specific clustering (Section 3.6)
DimPlot(merged_seurat, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample Distribution (UMAP)")
saveRDS(merged_seurat, "/scratch/SCWF00151/shared/c.d21085818/Dissertation/Analysis/merged_seurat_post_umap.rds")