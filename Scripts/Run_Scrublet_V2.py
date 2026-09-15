#!/usr/bin/env python3
"""
Run_Scrublet_V2.py

Doublet detection for all 12 samples, re-run against Cell Ranger outputs
generated using the corrected reference genome (Mus_custom_v2).

- Uses Scrublet's automatic, formula-based expected doublet rate by default:
      expected_doublet_rate = min((n_cells / 1000) * 0.008, 0.15)
- Supports manual per-sample threshold overrides via MANUAL_THRESHOLDS,
  to be filled in only after visually inspecting a sample's doublet
  score histogram (see companion plots written to Scrublet_outputs_v2/).
- Skips samples that have already been processed (resumable).
- Can rebuild the summary CSV from existing per-sample output files
  without re-running Scrublet, via REBUILD_SUMMARY_ONLY.

Requires Project_Config.txt in the same directory this script is run from.
"""

import os
import sys
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import scrublet as scr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

with open("Project_Config.txt") as f:
    lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]
    PROJECT_ROOT = lines[0]

BASE = os.path.join(PROJECT_ROOT, "Cellranger_outputs_v2")
OUT_DIR = os.path.join(PROJECT_ROOT, "Scrublet_outputs_v2")

SAMPLES = [
    "WT_Cas9_1_v2", "WT_Cas9_2_v2", "WT_Cas9_3_v2",
    "WT_Cas9_sg_1_v2", "WT_Cas9_sg_2_v2", "WT_Cas9_sg_3_v2",
    "R61_Cas9_1_v2", "R61_Cas9_2_v2", "R61_Cas9_3_v2",
    "R61_Cas9_sg_1_v2", "R61_Cas9_sg_2_v2", "R61_Cas9_sg_3_v2",
]

MANUAL_THRESHOLDS = {
    # "SAMPLE_NAME_v2": 0.6,
}

REBUILD_SUMMARY_ONLY = False

os.makedirs(OUT_DIR, exist_ok=True)
SUMMARY_PATH = os.path.join(OUT_DIR, "scrublet_summary.csv")


def load_counts_matrix(sample):
    mtx_dir = os.path.join(BASE, sample, "outs", "filtered_feature_bc_matrix")
    matrix_path = os.path.join(mtx_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(mtx_dir, "barcodes.tsv.gz")
    counts_matrix = sio.mmread(matrix_path).T.tocsc()
    barcodes = pd.read_csv(barcodes_path, header=None)[0].values
    return counts_matrix, barcodes


def expected_doublet_rate(n_cells):
    return min((n_cells / 1000) * 0.008, 0.15)


def process_sample(sample):
    out_csv = os.path.join(OUT_DIR, f"{sample}_scrublet.csv")
    if os.path.exists(out_csv):
        print(f"[skip] {sample}: already processed ({out_csv} exists)")
        existing = pd.read_csv(out_csv)
        n_cells = len(existing)
        n_doublets = int(existing["predicted_doublet"].sum())
        return sample, n_cells, n_doublets

    print(f"[run]  {sample}: loading counts matrix...")
    counts_matrix, barcodes = load_counts_matrix(sample)
    n_cells = counts_matrix.shape[0]
    exp_rate = expected_doublet_rate(n_cells)
    print(f"[run]  {sample}: {n_cells} cells, expected doublet rate = {exp_rate:.4f}")

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=exp_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    if sample in MANUAL_THRESHOLDS:
        manual_thresh = MANUAL_THRESHOLDS[sample]
        print(f"[run]  {sample}: applying manual threshold override = {manual_thresh}")
        predicted_doublets = doublet_scores > manual_thresh

    fig = scrub.plot_histogram()[0]
    fig.savefig(os.path.join(OUT_DIR, f"{sample}_doublet_score_histogram.png"), dpi=150)
    plt.close(fig)

    result_df = pd.DataFrame({
        "barcode": barcodes,
        "doublet_score": doublet_scores,
        "predicted_doublet": predicted_doublets,
    })
    result_df.to_csv(out_csv, index=False)
    n_doublets = int(predicted_doublets.sum())
    print(f"[done] {sample}: {n_doublets} / {n_cells} predicted doublets")
    return sample, n_cells, n_doublets


def rebuild_summary_from_disk():
    rows = []
    for sample in SAMPLES:
        out_csv = os.path.join(OUT_DIR, f"{sample}_scrublet.csv")
        if not os.path.exists(out_csv):
            print(f"[warn] {sample}: no output file found, skipping in summary")
            continue
        df = pd.read_csv(out_csv)
        n_cells = len(df)
        n_doublets = int(df["predicted_doublet"].sum())
        note = "manual threshold" if sample in MANUAL_THRESHOLDS else ""
        rows.append({"sample": sample, "n_cells": n_cells, "n_doublets": n_doublets, "note": note})
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(SUMMARY_PATH, index=False)
    print(f"Summary rebuilt from disk -> {SUMMARY_PATH}")
    print(summary_df.to_string(index=False))


def main():
    if REBUILD_SUMMARY_ONLY:
        rebuild_summary_from_disk()
        return

    rows = []
    for sample in SAMPLES:
        sample_dir = os.path.join(BASE, sample, "outs", "filtered_feature_bc_matrix")
        if not os.path.isdir(sample_dir):
            print(f"[warn] {sample}: expected output directory not found ({sample_dir}), skipping")
            continue
        sample_name, n_cells, n_doublets = process_sample(sample)
        note = "manual threshold" if sample in MANUAL_THRESHOLDS else ""
        rows.append({"sample": sample_name, "n_cells": n_cells, "n_doublets": n_doublets, "note": note})

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(SUMMARY_PATH, index=False)
    print(f"\nSummary written -> {SUMMARY_PATH}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
