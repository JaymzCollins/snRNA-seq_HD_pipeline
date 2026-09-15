#!/bin/bash
#SBATCH --job-name=mkref_Mus_custom_v2
#SBATCH --account=scwf00151_r_anney_149
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00

# Requires project_config.txt in the same directory this script is run from.
PROJECT_ROOT=$(grep -v '^#' Project_Config.txt | tail -n 1)

#SBATCH --output=PLACEHOLDER_OUT
#SBATCH --error=PLACEHOLDER_ERR

# Load Cell Ranger module
module load CellRanger/10.0.0

BASE=${PROJECT_ROOT}

# Move to Reference directory so the new reference folder is built here
cd ${PROJECT_ROOT}/Reference

# --------------------------------------------------------------------------
# Step 1: Filter the combined GTF to standard countable biotypes, following
# 10x Genomics' recommended reference-build practice (as used in their own
# GRCh38/GRCm39 reference construction scripts). Custom transgenes
# (Cas9D10A, GFP, huHTT) are all annotated as gene_biotype "protein_coding",
# so they are retained by this filter.
# --------------------------------------------------------------------------
cellranger mkgtf \
  ${BASE}/Mus_custom_v2.gtf \
  ${BASE}/Mus_custom_v2.filtered.gtf \
  --attribute=gene_biotype:protein_coding \
  --attribute=gene_biotype:lncRNA \
  --attribute=gene_biotype:antisense \
  --attribute=gene_biotype:IG_LV_gene \
  --attribute=gene_biotype:IG_V_gene \
  --attribute=gene_biotype:IG_V_pseudogene \
  --attribute=gene_biotype:IG_D_gene \
  --attribute=gene_biotype:IG_J_gene \
  --attribute=gene_biotype:IG_J_pseudogene \
  --attribute=gene_biotype:IG_C_gene \
  --attribute=gene_biotype:IG_C_pseudogene \
  --attribute=gene_biotype:TR_V_gene \
  --attribute=gene_biotype:TR_V_pseudogene \
  --attribute=gene_biotype:TR_D_gene \
  --attribute=gene_biotype:TR_J_gene \
  --attribute=gene_biotype:TR_J_pseudogene \
  --attribute=gene_biotype:TR_C_gene

# --------------------------------------------------------------------------
# Step 2: Build the corrected reference genome using the filtered GTF
# --------------------------------------------------------------------------
cellranger mkref \
  --genome=Mus_custom_v2 \
  --fasta=${BASE}/Mus_custom_v2.fa \
  --genes=${BASE}/Mus_custom_v2.filtered.gtf \
  --nthreads=8 \
  --memgb=64
