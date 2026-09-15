# ============================================================
# MERGE AND NORMALISE - ALL 12 SAMPLES
#
# Merges the 12 QC-filtered Seurat objects into a single combined
# dataset, adds sample-level metadata (sex, genotype, treatment),
# and applies SCTransform normalisation.
#
# Requires:
#   - seurat_list_filtered.rds, produced by 01_load_and_QC_filter.R
#     (PROJECT_ROOT/Analysis/seurat_list_filtered.rds)
#
# Produces:
#   - merged_seurat_with_sex.rds: merged object with sex, genotype,
#     and treatment metadata added, prior to normalisation
#   - merged_seurat_post_sctransform.rds: final checkpoint for this
#     script, SCTransform-normalised and ready for dimensionality
#     reduction (03_dimensionality_reduction.R)
# ============================================================

# ============================================================
# EDIT THIS: set to your own project root before running
# ============================================================
PROJECT_ROOT <- "/scratch/SCWF00151/shared/c.d21085818/Dissertation"

library(Seurat)

# Load the QC-filtered Seurat objects (produced by 01_load_and_QC_filter.R)
seurat_list_filtered <- readRDS(file.path(PROJECT_ROOT, "Analysis/seurat_list_filtered.rds"))
samples <- names(seurat_list_filtered)

# ============================================================
# STEP 7: Merge all 12 filtered objects into one combined dataset
# ============================================================
merged_seurat <- merge(
  x = seurat_list_filtered[[1]],
  y = seurat_list_filtered[2:length(seurat_list_filtered)],
  add.cell.ids = samples,
  project = "HD_striatum_v2"
)
cat("Total cells in merged object:", ncol(merged_seurat), "\n")  # expect 102,575
cat("Total genes:", nrow(merged_seurat), "\n")                    # expect 55,406
# Consolidate the 12 per-sample count layers into a single unified layer
merged_seurat <- JoinLayers(merged_seurat)
Layers(merged_seurat[["RNA"]])  # expect just "counts"
# ============================================================
# STEP 8: Add independently-verified sample-level metadata
# (sex, genotype, treatment)
# ============================================================
# Sex determined via Ddx3y (Y-linked) and Xist (X-linked) read counts,
# extracted directly from each sample's aligned BAM file.
# Concordance with recorded sample metadata was 12/12.
sex_lookup <- c(
  "WT_Cas9_1_v2" = "Female", "WT_Cas9_2_v2" = "Female", "WT_Cas9_3_v2" = "Female",
  "WT_Cas9_sg_1_v2" = "Female", "WT_Cas9_sg_2_v2" = "Male", "WT_Cas9_sg_3_v2" = "Female",
  "R61_Cas9_1_v2" = "Female", "R61_Cas9_2_v2" = "Male", "R61_Cas9_3_v2" = "Male",
  "R61_Cas9_sg_1_v2" = "Female", "R61_Cas9_sg_2_v2" = "Female", "R61_Cas9_sg_3_v2" = "Male"
)
# unname() is required: without it, the resulting vector retains sample
# names rather than cell barcode names, causing Seurat's assignment to
# fail to match any cells (a "no cell overlap" error)
merged_seurat$sex <- unname(sex_lookup[merged_seurat$orig.ident])
table(merged_seurat$orig.ident, merged_seurat$sex)  # verify 12/12 correct
# Genotype and treatment, derived from orig.ident naming
# convention, using the same lookup-table approach as sex (rather than
# regex on merged_seurat$orig.ident directly) for consistency and to
# avoid silent mismatches if any sample name doesn't fit an assumed pattern
genotype_lookup <- c(
  "WT_Cas9_1_v2" = "WT", "WT_Cas9_2_v2" = "WT", "WT_Cas9_3_v2" = "WT",
  "WT_Cas9_sg_1_v2" = "WT", "WT_Cas9_sg_2_v2" = "WT", "WT_Cas9_sg_3_v2" = "WT",
  "R61_Cas9_1_v2" = "R6.1", "R61_Cas9_2_v2" = "R6.1", "R61_Cas9_3_v2" = "R6.1",
  "R61_Cas9_sg_1_v2" = "R6.1", "R61_Cas9_sg_2_v2" = "R6.1", "R61_Cas9_sg_3_v2" = "R6.1"
)
treatment_lookup <- c(
  "WT_Cas9_1_v2" = "Cas9_only", "WT_Cas9_2_v2" = "Cas9_only", "WT_Cas9_3_v2" = "Cas9_only",
  "WT_Cas9_sg_1_v2" = "Cas9_sgRNA", "WT_Cas9_sg_2_v2" = "Cas9_sgRNA", "WT_Cas9_sg_3_v2" = "Cas9_sgRNA",
  "R61_Cas9_1_v2" = "Cas9_only", "R61_Cas9_2_v2" = "Cas9_only", "R61_Cas9_3_v2" = "Cas9_only",
  "R61_Cas9_sg_1_v2" = "Cas9_sgRNA", "R61_Cas9_sg_2_v2" = "Cas9_sgRNA", "R61_Cas9_sg_3_v2" = "Cas9_sgRNA"
)
merged_seurat$genotype <- unname(genotype_lookup[merged_seurat$orig.ident])
merged_seurat$treatment <- unname(treatment_lookup[merged_seurat$orig.ident])
# Verify 12/12 correct assignment, matching the sex verification above
table(merged_seurat$orig.ident, merged_seurat$genotype)
table(merged_seurat$orig.ident, merged_seurat$treatment)
# ============================================================
# Checkpoint: object is fully loaded, filtered, merged, and annotated
# with sex, genotype, and treatment, ready for normalisation
# ============================================================
saveRDS(merged_seurat, file.path(PROJECT_ROOT, "Analysis/merged_seurat_with_sex.rds"))
# ============================================================
# STEP 9: SCTransform normalisation (Section 2.8)
# ============================================================
# Increases the future package's globals size limit. SCTransform
# distributes several GB of data to parallel workers for a dataset
# this size. The default 500MB limit causes this step to fail with a
# "maxSizeOfObjects" error if not raised beforehand.
options(future.globals.maxSize = 20000 * 1024^2)  # 20GB
# vst.flavor = "v2": the improved SCTransform statistical model
# vars.to.regress: removes the influence of mitochondrial percentage,
#   total UMI count, and sex before identifying variable genes and
#   computing normalised (Pearson residual) values, matching the
#   normalisation covariates used by Murillo et al. (2025)
# conserve.memory = TRUE: processes genes in smaller batches to
#   reduce peak memory usage, at some cost to speed - necessary at
#   this dataset size
set.seed(123)
merged_seurat <- SCTransform(
  merged_seurat,
  vst.flavor = "v2",
  vars.to.regress = c("percent.mt", "nCount_RNA", "sex"),
  verbose = TRUE,
  conserve.memory = TRUE
)
# Confirm SCTransform completed successfully: expect 45,580 genes
# retained (down from 55,406 pre-filtering), 3,000 of which are
# flagged as highly variable, and default assay switched to "SCT".
# Observed runtime: ~16 minutes with 4 CPU cores, ~31GB peak memory usage.
cat("Genes retained:", nrow(merged_seurat), "\n")
cat("Default assay:", DefaultAssay(merged_seurat), "\n")
saveRDS(merged_seurat, file.path(PROJECT_ROOT, "Analysis/merged_seurat_post_sctransform.rds"))