#!/bin/bash
#SBATCH --job-name=cellranger_WT_Cas9_v2
#SBATCH --account=scwf00151_r_anney_149
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-11:00:00
#SBATCH --array=1-3
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

module load CellRanger/10.0.0

# Requires Project_Config.txt in the same directory this script is run from.
PROJECT_ROOT=$(grep -v '^#' Project_Config.txt | tail -n 1)

mkdir -p ${PROJECT_ROOT}/Cellranger_outputs_v2
cd ${PROJECT_ROOT}/Cellranger_outputs_v2

SAMPLE_ORIG="WT_Cas9_${SLURM_ARRAY_TASK_ID}"
SAMPLE_OUT="WT_Cas9_${SLURM_ARRAY_TASK_ID}_v2"

cellranger count \
  --id=${SAMPLE_OUT} \
  --transcriptome=${PROJECT_ROOT}/Reference/Mus_custom_v2 \
  --fastqs=${PROJECT_ROOT}/cellranger_fastqs \
  --sample=${SAMPLE_ORIG} \
  --create-bam=true \
  --localcores=8 \
  --localmem=64
