#!/bin/bash
# ==============================================================================
# Build corrected reference genome (Mus_custom_v2)
# Combines the standard Ensembl GRCm39 release 116 genome + annotation with
# the three custom transgene sequences (Cas9D10A, GFP, huHTT) extracted from
# the original (incomplete) Mus_custom reference.
#
# Requires project_config.txt in the same directory this script is run from,
# containing a single line: the root path of the project.
# ==============================================================================

set -euo pipefail

PROJECT_ROOT=$(grep -v '^#' Project_Config.txt | tail -n 1)
cd "${PROJECT_ROOT}"

echo "=== Step 1: Decompressing base FASTA and GTF ==="
gunzip -k Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip -k Mus_musculus.GRCm39.116.gtf.gz

echo "=== Step 2: Combining FASTA files ==="
cat Mus_musculus.GRCm39.dna.primary_assembly.fa custom_transgenes.fa > Mus_custom_v2.fa

echo "=== Step 3: Combining GTF files ==="
cat Mus_musculus.GRCm39.116.gtf custom_transgenes.gtf > Mus_custom_v2.gtf

echo "=== Step 4: Sanity checks ==="
echo "Total contigs in combined FASTA:"
grep -c ">" Mus_custom_v2.fa

echo ""
echo "Last 5 contigs (should show Cas9D10A, GFP, huHTT):"
grep ">" Mus_custom_v2.fa | tail -5

echo ""
echo "Checking for MT and Y chromosomes:"
grep ">" Mus_custom_v2.fa | grep -E "^>MT|^>Y" || echo "WARNING: MT or Y not found!"

echo ""
echo "=== Done. Combined files ready: Mus_custom_v2.fa, Mus_custom_v2.gtf ==="
