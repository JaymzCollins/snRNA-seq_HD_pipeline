# ============================================================
# QUALITY CONTROL FILTERING - ALL 12 SAMPLES
#
# Loads raw Cell Ranger outputs, calculates per-cell QC metrics
# (mitochondrial percentage, gene count), imports Scrublet doublet
# predictions, and applies combined per-cell filtering.
#
# Requires:
#   - Cell Ranger outputs in PROJECT_ROOT/Cellranger_outputs_v2/
#     (one subfolder per sample, each containing outs/filtered_feature_bc_matrix)
#   - Scrublet outputs in PROJECT_ROOT/Scrublet_outputs_v2/
#     (one CSV per sample, containing barcode, predicted_doublet,
#     and doublet_score columns)
#
# Produces:
#   - seurat_list_filtered: a list of 12 QC-filtered Seurat objects,
#     one per sample
#   - QC_filtering_summary_per_sample.csv: per-sample cell counts and
#     per-criterion exclusion counts before/after filtering
# ============================================================



# ============================================================
# EDIT THIS: set to your own project root before running
# ============================================================
PROJECT_ROOT <- "/scratch/SCWF00151/shared/c.d21085818/Dissertation"

# Ensure the Analysis output directory exists (created if missing)
dir.create(file.path(PROJECT_ROOT, "Analysis"), showWarnings = FALSE, recursive = TRUE)



# ============================================================
# Required Packages
# ============================================================
library(Seurat)
library(harmony)

# ============================================================
# STEP 1-2: Load raw Cell Ranger outputs into individual Seurat objects
# ============================================================
base_path <- file.path(PROJECT_ROOT, "Cellranger_outputs_v2")
samples <- c("WT_Cas9_1_v2", "WT_Cas9_2_v2", "WT_Cas9_3_v2",
             "WT_Cas9_sg_1_v2", "WT_Cas9_sg_2_v2", "WT_Cas9_sg_3_v2",
             "R61_Cas9_1_v2", "R61_Cas9_2_v2", "R61_Cas9_3_v2",
             "R61_Cas9_sg_1_v2", "R61_Cas9_sg_2_v2", "R61_Cas9_sg_3_v2")
seurat_list <- list()
for (sample in samples) {
  cat("Loading:", sample, "\n")
  data_dir <- file.path(base_path, sample, "outs", "filtered_feature_bc_matrix")
  counts <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)
  seurat_list[[sample]] <- seurat_obj
}
# Sanity check: expect 10,017 / 2,142 / 2,383 / 2,384 / 7,964 / 6,000 /
# 9,298 / 10,918 / 27,983 / 10,326 / 9,918 / 8,486 cells respectively,
# 55,406 genes for every sample
for (sample in samples) {
  cat(sample, ":", ncol(seurat_list[[sample]]), "cells,", nrow(seurat_list[[sample]]), "genes\n")
}
# ============================================================
# STEP 3: Calculate mitochondrial read percentage per cell
# ============================================================
# pattern = "^mt-": matches mouse mitochondrial gene naming convention
# (e.g. mt-Nd1, mt-Co1)
for (sample in samples) {
  seurat_list[[sample]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_list[[sample]], pattern = "^mt-"
  )
}
# Verification: confirms the calculation completed correctly
for (sample in samples) {
  cat(sample, "- Median mito%:", median(seurat_list[[sample]]$percent.mt),
      "| Max mito%:", max(seurat_list[[sample]]$percent.mt), "\n")
}
# ============================================================
# STEP 4: Global adaptive gene-count threshold (Section 3.4)
# ============================================================
# A fixed upper threshold (Murillo et al., 2025: >5000 genes) was found
# to produce inconsistent, sample-composition-dependent exclusion rates
# across this multi-cell-type dataset. A global adaptive threshold -
# calculated once from the pooled distribution of all 12 samples,
# rather than per-sample - was adopted instead, avoiding both the
# over-exclusion of high-complexity samples and the distortion caused
# by samples with atypical cell-type composition.
all_nfeat <- unlist(lapply(samples, function(s) seurat_list[[s]]$nFeature_RNA))
global_median <- median(all_nfeat)
global_mad <- mad(all_nfeat)
global_upper <- global_median + 3 * global_mad
global_lower <- max(0, global_median - 3 * global_mad)
# Confirmed values: median = 1748, MAD = 1455.9, upper = 6116, lower = 0
cat("Global median:", global_median, "\n")
cat("Global MAD:", global_mad, "\n")
cat("Global adaptive upper threshold:", round(global_upper, 0), "\n")
cat("Global adaptive lower threshold:", round(global_lower, 0), "\n")
# ============================================================
# STEP 5: Import Scrublet doublet predictions
# ============================================================
# Imports both the binary call (predicted_doublet) and
# the continuous doublet_score, matched by barcode for each. The
# continuous score is needed for downstream cluster-level QC checks
# (e.g. comparing doublet score distributions across clusters), which
# the binary call alone cannot support - a called doublet only tells
# you the value crossed Scrublet's (under-calling) automatic
# threshold, not how borderline nearby uncalled cells were.
for (sample in samples) {
  scrublet_path <- file.path(PROJECT_ROOT, "Scrublet_outputs_v2",
                             paste0(sample, "_scrublet.csv"))
  scrublet_results <- read.csv(scrublet_path)
  
  # Python's "True"/"False" strings are read as character by read.csv();
  # explicit conversion to logical is required before use
  scrublet_results$predicted_doublet <- as.logical(scrublet_results$predicted_doublet)
  
  match_idx <- match(colnames(seurat_list[[sample]]), scrublet_results$barcode)
  
  seurat_list[[sample]]$predicted_doublet <- scrublet_results$predicted_doublet[match_idx]
  seurat_list[[sample]]$doublet_score <- scrublet_results$doublet_score[match_idx]
  
  cat(sample, "- doublets matched:", sum(seurat_list[[sample]]$predicted_doublet, na.rm = TRUE),
      "out of", ncol(seurat_list[[sample]]), "cells",
      "| doublet_score NAs:", sum(is.na(seurat_list[[sample]]$doublet_score)), "\n")
}

# ============================================================
# STEP 6: Apply combined QC filter to all 12 samples
# ============================================================
# Criteria applied simultaneously:
#   - nFeature_RNA > 200: removes likely empty droplets/degraded nuclei
#   - nFeature_RNA < 6116: global adaptive upper threshold (Step 4)
#   - percent.mt < 5: standard mitochondrial QC threshold, matching
#     Murillo et al. (2025); contributed negligibly to filtering (0-2
#     cells per sample), consistent with the low mitochondrial content
#     typical of single-nucleus preparations
#   - predicted_doublet == FALSE: excludes Scrublet-flagged doublets
#
# Per-criterion exclusion counts are now tracked and saved to
# qc_summary, alongside the original before/after totals. This is
# purely additional bookkeeping - the filtering logic itself
# (the subset() call) is byte-for-byte unchanged, so results are
# identical.
qc_summary <- data.frame(
  sample = character(), cells_before = integer(),
  excl_nFeature_low = integer(), excl_nFeature_high = integer(),
  excl_percent_mt = integer(), excl_doublet = integer(),
  cells_after = integer(), pct_retained = numeric(),
  stringsAsFactors = FALSE
)
seurat_list_filtered <- list()
for (sample in samples) {
  obj <- seurat_list[[sample]]
  n_before <- ncol(obj)
  
  # Per-criterion exclusion counts (not mutually exclusive -
  # a cell can fail more than one criterion)
  fail_low <- sum(obj$nFeature_RNA <= 200)
  fail_high <- sum(obj$nFeature_RNA >= 6116)
  fail_mt <- sum(obj$percent.mt >= 5)
  fail_doublet <- sum(obj$predicted_doublet == TRUE, na.rm = TRUE)
  
  obj_filtered <- subset(obj, subset = nFeature_RNA > 200 &
                           nFeature_RNA < 6116 &
                           percent.mt < 5 &
                           predicted_doublet == FALSE)
  
  n_after <- ncol(obj_filtered)
  cat(sample, ": before =", n_before, "| after =", n_after,
      "| removed =", n_before - n_after,
      paste0("(", round(100*(n_before-n_after)/n_before, 1), "%)\n"))
  
  # Append this sample's row to qc_summary
  qc_summary <- rbind(qc_summary, data.frame(
    sample = sample, cells_before = n_before,
    excl_nFeature_low = fail_low, excl_nFeature_high = fail_high,
    excl_percent_mt = fail_mt, excl_doublet = fail_doublet,
    cells_after = n_after,
    pct_retained = round(100 * n_after / n_before, 1)
  ))
  
  seurat_list_filtered[[sample]] <- obj_filtered
}

# Expect total: 107,819 -> 102,575 cells (4.9% removed) across the dataset
# Save the per-sample QC summary table
print(qc_summary)
write.csv(qc_summary, file.path(PROJECT_ROOT, "Analysis/QC_filtering_summary_per_sample.csv"), row.names = FALSE)

# ============================================================
# Checkpoint: save the list of QC-filtered Seurat objects, so
# Script 02 (merge_normalise.R) can load this directly without
# needing to re-run this script in the same session
# ============================================================
saveRDS(seurat_list_filtered, file.path(PROJECT_ROOT, "Analysis/seurat_list_filtered.rds"))