#!/bin/bash
#SBATCH --job-name=cellranger_WT_Cas9_sg_2_v2
#SBATCH --account=scwf00151_r_anney_149
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load CellRanger/10.0.0

# Requires Project_Config.txt in the same directory this script is run from.
PROJECT_ROOT=$(grep -v '^#' Project_Config.txt | tail -n 1)

mkdir -p ${PROJECT_ROOT}/Cellranger_outputs_v2
cd ${PROJECT_ROOT}/Cellranger_outputs_v2

cellranger count \
  --id=WT_Cas9_sg_2_v2 \
  --transcriptome=${PROJECT_ROOT}/Reference/Mus_custom_v2 \
  --fastqs=${PROJECT_ROOT}/cellranger_fastqs \
  --sample=WT_Cas9_sg_2 \
  --create-bam=true \
  --localcores=8 \
  --localmem=64
