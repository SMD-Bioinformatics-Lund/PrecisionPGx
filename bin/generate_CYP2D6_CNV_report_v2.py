#!/usr/bin/env python

import argparse
import os
import pandas as pd
import numpy as np
import logging
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import re

# Read input CNV file
def read_cnv_file(cnv_file):
    """
    Reads a CNV file into a Pandas DataFrame.

    Parameters:
    cnv_file (str): Path to the CNV file.

    Returns:
    pd.DataFrame: DataFrame containing the CNV data.
    
    """

    df = pd.read_csv(cnv_file, dtype={0: str})

    return df

def gene_starts_and_ends(bed_file):
    """
    Reads a BED file and extracts the start and end positions of unique genes.

    Parameters:
    bed_file (str): Path to the BED file.

    Returns:
    pd.DataFrame: DataFrame with columns ["chr", "start", "end", "gene"].
    """
    # Read the BED file (assuming tab-separated values with columns: chr, start, end, gene)
    df_bed = pd.read_csv(bed_file, sep="\t", names=["chr", "start", "end", "gene"])

    gene_limits = []
    this_gene = None
    this_start = None

    for index, row in df_bed.iterrows():
        gene_name_parts = str(row["gene"]).split("_")

        # Handling gene naming convention
        if "/" in str(row["gene"]) and "ENS" not in str(row["gene"]):
            this_line_gene_name = f"{gene_name_parts[0]}_{gene_name_parts[1]}"
        else:
            this_line_gene_name = gene_name_parts[0]

        # If this is the first gene encountered
        if this_gene is None:
            this_gene = this_line_gene_name
            this_start = row["start"] + 1  # Convert 0-based start to 1-based

        # If a new gene is encountered, store previous gene details
        if this_line_gene_name != this_gene:
            gene_limits.append([df_bed.loc[index - 1, "chr"], this_start, df_bed.loc[index - 1, "end"], this_gene])
            #print(df_bed.loc[index - 1, "chr"],[row["chr"], this_start, df_bed.loc[index - 1, "end"],this_gene])

            # Reset values for the new gene
            this_gene = this_line_gene_name
            this_start = row["start"] + 1
            #print ("exception2")
            #print (this_gene)
            #print (this_start)

        # If it's the last row, store the final gene details
        if index == len(df_bed) - 1:
            gene_limits.append([row["chr"], this_start, row["end"], this_gene])

    return pd.DataFrame(gene_limits, columns=["chr", "start", "end", "gene"])

def plot_complete_panel(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, output_dir):
    """
    Generates plots for CNV analysis including Z-score, Copy Number, and Read Depth.
    
    Parameters:
    - data_frame: Pandas DataFrame containing CNV data (columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber)
    - gene_lims: DataFrame containing chromosome start and end points for genes
    - sample_name: Name of the sample for output file naming
    - z_threshold: Z-score threshold
    - cn_threshold: Copy number threshold
    - output_dir: Output directory to save the plots
    """
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Get unique chromosomes and first row indices where they appear
    panel_chrs = data_frame["Chr"].unique()
    panel_chrs_row_no = {chr_: data_frame.index[data_frame["Chr"] == chr_][0] for chr_ in panel_chrs}

    # Extract series for plotting
    y_series_z = data_frame["zscore"]
    y_series_cn = data_frame["copynumber"]
    y_series_depth = data_frame["Depth"]
    y_series_expDepth = data_frame["expDepth"]
    y_series_stdev = np.abs((data_frame["Depth"] / data_frame["prop"]) * data_frame["std"])

    # Identify points above thresholds (for red markers)
    threshold_series = (np.abs(y_series_z) > z_threshold) & (np.abs(y_series_cn - 2) > cn_threshold)

    # Create figure with 3 subplots
    fig, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)
    
    # Read Depth Plot
    axes[0].errorbar(range(len(y_series_depth)), y_series_expDepth, yerr=y_series_stdev, fmt='o', color="gray", label="Expected Depth")
    axes[0].scatter(range(len(y_series_depth)), y_series_depth, c=np.where(threshold_series, "red", "black"), s=5, label="Observed Depth")
    axes[0].axhline(y=20, linestyle="dashed", color="gray")  # Threshold line
    axes[0].set_ylabel("Read Depth")
    axes[0].legend()

    # Copy Number Plot
    axes[1].scatter(range(len(y_series_cn)), y_series_cn, c=np.where(threshold_series, "red", "black"), s=5, label="Copy Number")
    axes[1].axhline(y=2, linestyle="dashed", color="gray")
    axes[1].axhline(y=2 + cn_threshold, linestyle="dotted", color="gray")
    axes[1].axhline(y=2 - cn_threshold, linestyle="dotted", color="gray")
    axes[1].set_ylabel("Copy Number")

    # Z-score Plot
    axes[2].scatter(range(len(y_series_z)), y_series_z, c=np.where(threshold_series, "red", "black"), s=5, label="Z-score")
    axes[2].axhline(y=0, linestyle="dashed", color="gray")
    axes[2].axhline(y=z_threshold, linestyle="dotted", color="gray")
    axes[2].axhline(y=-z_threshold, linestyle="dotted", color="gray")
    axes[2].set_ylabel("Z-score")
    axes[2].set_xlabel("Positions covered by panel")

    # Add vertical lines for chromosomes
    for chr_, row_no in panel_chrs_row_no.items():
        for ax in axes:
            ax.axvline(x=row_no, linestyle="dotted", color="gray")
            ax.text(row_no, ax.get_ylim()[1] * 0.9, chr_, rotation=90, verticalalignment="top")

    # Title and save plot
    plt.suptitle(f"{sample_name} - Complete Panel CNV Analysis", fontsize=14)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{sample_name}_cnv_complete_panel.png")
    plt.savefig(output_path)
    plt.close()

    print(f"Plot saved at {output_path}")

def cnv_report(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, prop_threshold, bed_file, flag_bed_file, output_dir):
    """
    Generates a CNV report based on z-score and copy number thresholds.

    Parameters:
    - data_frame: Pandas DataFrame with columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: Pandas DataFrame containing chromosome start and end positions for genes
    - sample_name: Sample name prefix for output files
    - z_threshold: Threshold for the z-score
    - cn_threshold: Threshold for the copy number deviation from 2
    - prop_threshold: Proportion threshold for calling CNVs
    - bed_file: Path to BED file containing gene annotations
    - flag_bed_file: Path to flag BED file with problematic genomic regions
    - output_dir: Directory to save the CNV report

    Returns:
    - DataFrame of CNV regions
    """

    print(f"Creating CNV report for {sample_name}...")

    report_file = os.path.join(output_dir, f"{sample_name}_cnv_report_all.txt")
    filtered_report_file = os.path.join(output_dir, f"{sample_name}_cnv_report_filtered.txt")

    # Read BED and flag BED files
    bed_df = pd.read_csv(bed_file, sep="\t", names=["chr", "start_pos", "end_pos", "gene"])
    flag_bed_df = pd.read_csv(flag_bed_file, sep="\t", names=["chr", "start_pos", "end_pos", "annot"])

    # Identify positions exceeding both CNV and Z-score thresholds
    threshold_series = (data_frame["zscore"].abs() > z_threshold) & ((data_frame["copynumber"] - 2).abs() > cn_threshold)

    cnv_records = []

    if threshold_series.any():
        # Process each gene region from BED file
        for _, gene_row in bed_df.iterrows():
            chr_region = gene_row["chr"]
            start = gene_row["start_pos"]
            end = gene_row["end_pos"]
            gene_name = gene_row["gene"]

            # Filter CNV data for the current gene region
            gene_cnv_df = data_frame[(data_frame["Chr"] == chr_region) & 
                                     (data_frame["Position"] >= start) & 
                                     (data_frame["Position"] <= end)]
            
            if not gene_cnv_df.empty:
                # Calculate mean CN and Z-score
                mean_cn = gene_cnv_df["copynumber"].mean()
                mean_z = gene_cnv_df["zscore"].mean()
                mean_depth = gene_cnv_df["Depth"].mean()
                exp_mean_depth = gene_cnv_df["expDepth"].mean()
                
                # Determine CNV type
                if mean_cn > 2:
                    cnv_type = "DUP"
                elif mean_cn < 2:
                    cnv_type = "DEL"
                else:
                    cnv_type = "UNDETERMINED"

                # Count positions above threshold
                reg_length = len(gene_cnv_df)
                reg_pos_above = threshold_series.loc[gene_cnv_df.index].sum()
                reg_percent_above = reg_pos_above / reg_length

                # Flagging CNV
                reg_flag = "PASS" if reg_percent_above >= prop_threshold else "FAIL"
                reg_comment = "FAIL - Low proportion above threshold. " if reg_flag == "FAIL" else ""

                # Check for low mean expected depth
                if exp_mean_depth < 50:
                    reg_flag = "FAIL"
                    reg_comment += "FAIL - Cohort mean depth too low. "

                # Check for overlap with problematic regions
                flag_overlaps = flag_bed_df[(flag_bed_df["chr"] == chr_region) &
                                            ((flag_bed_df["start_pos"] <= end) & 
                                             (flag_bed_df["end_pos"] >= start))]

                if not flag_overlaps.empty:
                    reg_comment += "WARNING - Overlaps with problematic region: "
                    for _, flag_row in flag_overlaps.iterrows():
                        reg_comment += f"{flag_row['chr']}:{flag_row['start_pos']}-{flag_row['end_pos']} {flag_row['annot']} "

                # Save CNV data
                cnv_records.append([chr_region, start, end, gene_name, cnv_type, mean_depth, exp_mean_depth, 
                                    mean_cn, mean_z, reg_length, reg_pos_above, reg_percent_above, reg_flag, reg_comment.strip()])

    # Convert to DataFrame
    cnv_df = pd.DataFrame(cnv_records, columns=[
        "chr", "start_pos", "end_pos", "gene", "cnv", "mean_depth", "exp_mean_depth", "mean_cn", 
        "mean_z-score", "reg.length", "pos.above", "%above", "flag", "comment"
    ])

    # Save full CNV report
    cnv_df.to_csv(report_file, sep="\t", index=False)
    print(f"Report saved: {report_file}")

    # Save only PASS CNVs
    cnv_df_filtered = cnv_df[cnv_df["flag"] == "PASS"]
    cnv_df_filtered.to_csv(filtered_report_file, sep="\t", index=False)
    print(f"Filtered report saved: {filtered_report_file}")

    return cnv_df

def plot_complete_panel_all_pos(data_frame, gene_lims, sample_name, chr_length_file, z_threshold, cn_threshold, output_dir):
    """
    Plots CNV data across all chromosome positions.

    Parameters:
    - data_frame: Pandas DataFrame containing columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: Pandas DataFrame with gene start and end positions
    - sample_name: Name of the sample for plot title and output file naming
    - chr_length_file: Path to chromosome length file
    - z_threshold: Threshold for Z-score
    - cn_threshold: Threshold for Copy Number deviation from 2
    - output_dir: Directory to save the plots

    Saves:
    - A combined genome-wide CNV plot with Depth, CNV, and Z-score across all chromosomes.
    """

    # Load chromosome length file
    chr_lengths = pd.read_csv(chr_length_file, sep="\t", names=["chr", "length", "cumulative_length"])

    # Compute adjusted genome-wide positions
    adj_pos = data_frame.apply(lambda row: row["Position"] + 
                               chr_lengths.loc[chr_lengths["chr"] == row["Chr"], "cumulative_length"].values[0], axis=1)

    # Identify positions exceeding both CNV and Z-score thresholds
    threshold_series = (data_frame["zscore"].abs() > z_threshold) & ((data_frame["copynumber"] - 2).abs() > cn_threshold)

    # Prepare data for plotting
    x_series = adj_pos
    y_series_z = data_frame["zscore"]
    y_series_cn = data_frame["copynumber"]
    y_series_depth = data_frame["Depth"]
    y_series_expDepth = data_frame["expDepth"]
    y_series_stdev = np.abs((data_frame["Depth"] / data_frame["prop"]) * data_frame["std"])

    # Setup plot layout
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 8), sharex=True)

    # Plot Read Depth
    axes[0].errorbar(x_series, y_series_expDepth, yerr=y_series_stdev, fmt='o', markersize=2, color='gray', label="Expected Depth")
    axes[0].scatter(x_series, y_series_depth, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[0].axhline(y=20, linestyle='dashed', color='gray')
    axes[0].set_ylabel("Read Depth")
    axes[0].legend()

    # Plot Copy Number
    axes[1].scatter(x_series, y_series_cn, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[1].axhline(y=2, linestyle='dashed', color='gray')
    axes[1].axhline(y=2 + cn_threshold, linestyle='dotted', color='gray')
    axes[1].axhline(y=2 - cn_threshold, linestyle='dotted', color='gray')
    axes[1].set_ylabel("Copy Number")
    axes[1].set_ylim(0, max(y_series_cn.max(), 4))

    # Plot Z-score
    axes[2].scatter(x_series, y_series_z, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[2].axhline(y=0, linestyle='dashed', color='gray')
    axes[2].axhline(y=z_threshold, linestyle='dotted', color='gray')
    axes[2].axhline(y=-z_threshold, linestyle='dotted', color='gray')
    axes[2].set_ylabel("Z-score")
    axes[2].set_xlabel("Genomic Position")

    # Add chromosome separation lines
    for _, row in chr_lengths.iloc[1:].iterrows():
        for ax in axes:
            ax.axvline(x=row["cumulative_length"], linestyle='dotted', color='gray')

    # Add chromosome labels
    for i, row in chr_lengths.iloc[1:].iterrows():
        mid_chr_pos = row["cumulative_length"] - (row["length"] / 2)
        axes[-1].text(mid_chr_pos, z_threshold * 1.2, row["chr"], ha='center', fontsize=8)

    # Final plot settings
    fig.suptitle(f"{sample_name} - Complete Panel")
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    # Save plot
    output_file = os.path.join(output_dir, f"{sample_name}_cnv_whole_genome.png")
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved as {output_file}")
    plt.close()

def plot_chrs(data_frame, gene_lims, sample_name, chr_list, z_threshold, cn_threshold, output_dir):
    """
    Plots CNV data for each chromosome separately.

    Parameters:
    - data_frame: Pandas DataFrame containing columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: Pandas DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Name of the sample for plot title and output file naming
    - chr_list: List of chromosomes to plot
    - z_threshold: Threshold for Z-score
    - cn_threshold: Threshold for Copy Number deviation from 2
    - output_dir: Directory to save the plots

    Saves:
    - PNG plots for each chromosome.
    """

    print("Plotting chromosomes...")

    # Iterate through each chromosome in the list
    for chrom in chr_list:
        print(f"Processing {chrom}...")

        # Filter data for the current chromosome
        this_chr_df = data_frame[data_frame["Chr"] == chrom].copy()
        this_chr_gene_lims = gene_lims[gene_lims["chr"] == chrom].copy()

        if this_chr_df.empty:
            print(f"No data for chromosome {chrom}. Skipping...")
            continue

        # Create vertical line positions for genes
        gene_starts = []
        gene_ends = []
        gene_labels = []

        for _, row in this_chr_gene_lims.iterrows():
            if row["start"] in this_chr_df["Position"].values and row["end"] in this_chr_df["Position"].values:
                gene_starts.append(row["start"])
                gene_ends.append(row["end"])
                gene_labels.append(row["gene"])

        # Identify CNV and Z-score threshold exceedances
        threshold_series = (this_chr_df["zscore"].abs() > z_threshold) & ((this_chr_df["copynumber"] - 2).abs() > cn_threshold)

        # Extract values for plotting
        x_series = this_chr_df["Position"]
        y_series_z = this_chr_df["zscore"]
        y_series_cn = this_chr_df["copynumber"]
        y_series_depth = this_chr_df["Depth"]
        y_series_expDepth = this_chr_df["expDepth"]
        y_series_stdev = np.abs((this_chr_df["Depth"] / this_chr_df["prop"]) * this_chr_df["std"])

        # Setup plot layout
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 8), sharex=True)

        # Read Depth Plot
        axes[0].errorbar(x_series, y_series_expDepth, yerr=y_series_stdev, fmt='o', markersize=2, color='gray', label="Expected Depth")
        axes[0].scatter(x_series, y_series_depth, c=threshold_series.map({True: "red", False: "black"}), s=3)
        axes[0].axhline(y=20, linestyle='dashed', color='gray')
        axes[0].set_ylabel("Read Depth")
        axes[0].legend()

        # Copy Number Plot
        axes[1].scatter(x_series, y_series_cn, c=threshold_series.map({True: "red", False: "black"}), s=3)
        axes[1].axhline(y=2, linestyle='dashed', color='gray')
        axes[1].axhline(y=2 + cn_threshold, linestyle='dotted', color='gray')
        axes[1].axhline(y=2 - cn_threshold, linestyle='dotted', color='gray')
        axes[1].set_ylabel("Copy Number")
        axes[1].set_ylim(0, max(y_series_cn.max(), 4))

        # Z-score Plot
        axes[2].scatter(x_series, y_series_z, c=threshold_series.map({True: "red", False: "black"}), s=3)
        axes[2].axhline(y=0, linestyle='dashed', color='gray')
        axes[2].axhline(y=z_threshold, linestyle='dotted', color='gray')
        axes[2].axhline(y=-z_threshold, linestyle='dotted', color='gray')
        axes[2].set_ylabel("Z-score")
        axes[2].set_xlabel("Genomic Position")

        # Add vertical gene lines
        for ax in axes:
            for start, end in zip(gene_starts, gene_ends):
                ax.axvline(x=start, linestyle='solid', color='black')
                ax.axvline(x=end, linestyle='solid', color='black')

        # Add gene annotations
        for start, gene in zip(gene_starts, gene_labels):
            axes[2].text(start, z_threshold * 1.2, gene, ha='center', fontsize=8, rotation=45)

        # Final plot settings
        fig.suptitle(f"{sample_name} - Chromosome {chrom}")
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])

        # Save plot
        output_file = os.path.join(output_dir, f"{sample_name}_cnv_{chrom}.png")
        plt.savefig(output_file, dpi=300)
        print(f"Plot saved as {output_file}")
        plt.close()

def plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chrom, start_pos, end_pos, bed_file, output_dir):
    """
    Plots CNV data for a specific genomic region.

    Parameters:
    - data_frame: Pandas DataFrame containing columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Sample name for labeling and file naming
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - chrom: Chromosome of the region
    - start_pos: Start position of the region
    - end_pos: End position of the region
    - bed_file: Path to BED file
    - output_dir: Directory to save plots

    Saves:
    - PNG plot of the specified genomic region
    """
    print(f"Creating plot for region {chrom}:{start_pos}-{end_pos}")

    chr = chrom.replace("chr","")   

    # Load BED file containing annotation regions
    bed_df = pd.read_csv(bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "annot"])

    ## print(data_frame.head())

    # Filter CNV data for the given region
    this_reg_df = data_frame[
        (data_frame["Chr"] == chr) & 
        (data_frame["Position"] >= start_pos) & 
        (data_frame["Position"] <= end_pos)
    ].copy()


    if this_reg_df.empty:
        print(f"No data available for region {chrom}:{start_pos}-{end_pos}. Skipping...")
        return

    # Filter gene limits for the given region
    this_reg_gene_lims = gene_lims[
        (gene_lims["chr"] == chrom) & 
        (gene_lims["start"] < end_pos) & 
        (gene_lims["end"] > start_pos)
    ].copy()

    # Identify CNV and Z-score threshold exceedances
    threshold_series = (this_reg_df["zscore"].abs() > z_threshold) & ((this_reg_df["copynumber"] - 2).abs() > cn_threshold)

    # Extract values for plotting
    x_series = this_reg_df["Position"]
    y_series_z = this_reg_df["zscore"]
    y_series_cn = this_reg_df["copynumber"]
    y_series_depth = this_reg_df["Depth"]
    y_series_expDepth = this_reg_df["expDepth"]
    y_series_stdev = np.abs((this_reg_df["Depth"] / this_reg_df["prop"]) * this_reg_df["std"])

    # Setup plot layout
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 8), sharex=True)

    # Read Depth Plot
    axes[0].errorbar(x_series, y_series_expDepth, yerr=y_series_stdev, fmt='o', markersize=2, color='gray', label="Expected Depth")
    axes[0].scatter(x_series, y_series_depth, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[0].axhline(y=20, linestyle='dashed', color='gray')
    axes[0].set_ylabel("Read Depth")
    axes[0].legend()

    # Copy Number Plot
    axes[1].scatter(x_series, y_series_cn, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[1].axhline(y=2, linestyle='dashed', color='gray')
    axes[1].axhline(y=2 + cn_threshold, linestyle='dotted', color='gray')
    axes[1].axhline(y=2 - cn_threshold, linestyle='dotted', color='gray')
    axes[1].set_ylabel("Copy Number")
    axes[1].set_ylim(0, max(y_series_cn.max(), 4))

    # Z-score Plot
    axes[2].scatter(x_series, y_series_z, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[2].axhline(y=0, linestyle='dashed', color='gray')
    axes[2].axhline(y=z_threshold, linestyle='dotted', color='gray')
    axes[2].axhline(y=-z_threshold, linestyle='dotted', color='gray')
    axes[2].set_ylabel("Z-score")
    axes[2].set_xlabel("Genomic Position")

    # Final plot settings
    fig.suptitle(f"{sample_name} - {chrom}:{start_pos}-{end_pos}")
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    # Save plot
    output_file = os.path.join(output_dir, f"{sample_name}_{chrom}_{start_pos}_{end_pos}.png")
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved as {output_file}")
    plt.close()

def plot_cnv_regions(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, cnv_regions, bed_file, output_dir):
    """
    Plots CNV regions by generating zoomed-in plots for significant CNVs.

    Parameters:
    - data_frame: Pandas DataFrame containing CNV data
    - gene_lims: DataFrame containing gene start/end positions
    - sample_name: Sample name for labeling
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - cnv_regions: DataFrame with detected CNVs (chr, start_pos, end_pos, gene, etc.)
    - bed_file: Path to BED file
    - output_dir: Directory to save plots
    """
    print("Generating CNV region plots...")

    # Load BED file containing annotation regions
    bed_df = pd.read_csv(bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "gene"])

    # Get unique chromosomes with CNVs
    cnv_chrs = cnv_regions["chr"].unique()

    # Generate chromosome-wide plots
    for chrom in cnv_chrs:
        plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chrom, 0, float("inf"), bed_file, output_dir)

    # Determine zoomed-in CNV plot regions
    cnv_plot_regions = []
    
    for _, row in cnv_regions.iterrows():
        chrom = row["chr"]
        start_pos = row["start_pos"]
        end_pos = row["end_pos"]

        # Determine start and end positions for zoomed-in plot
        cnv_length = end_pos - start_pos
        cnv_length = max(cnv_length, 40)

        # Locate bed file regions covering the CNV
        bed_start_idx = bed_df[
            (bed_df["chr"] == chrom) & 
            (bed_df["start_pos"] <= start_pos) & 
            (bed_df["end_pos"] >= start_pos)
        ].index.min()

        bed_end_idx = bed_df[
            (bed_df["chr"] == chrom) & 
            (bed_df["start_pos"] <= end_pos) & 
            (bed_df["end_pos"] >= end_pos)
        ].index.max()

        # Adjust start and end positions
        plot_start_pos = bed_df.iloc[max(0, bed_start_idx - 3)]["start_pos"] if bed_start_idx is not None else start_pos - 2 * cnv_length
        plot_end_pos = bed_df.iloc[min(len(bed_df) - 1, bed_end_idx + 3)]["end_pos"] if bed_end_idx is not None else end_pos + 2 * cnv_length

        # Ensure plot regions do not overlap excessively
        if not cnv_plot_regions or (cnv_plot_regions[-1][1] < plot_start_pos):
            cnv_plot_regions.append((chrom, plot_start_pos, plot_end_pos))

    # Generate zoomed-in CNV plots
    for chrom, start_pos, end_pos in cnv_plot_regions:
        plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chrom, start_pos, end_pos, bed_file, output_dir)

def cyp2d6_report(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, prop_threshold, bed_file, flag_bed_file, output_dir):
    """
    Generates a CYP2D6 CNV report.

    Parameters:
    - data_frame: Pandas DataFrame containing CNV data (Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber)
    - gene_lims: DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Sample name for labeling
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - prop_threshold: Proportion threshold for CNV detection
    - bed_file: Path to BED file
    - flag_bed_file: Path to BED file containing flagged regions
    - output_dir: Directory to save reports and plots
    """
    print(f"Creating CYP2D6 CNV report for {sample_name}...")

    # Load BED and flagged regions data
    bed_df = pd.read_csv(bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "gene"])
    flag_bed_df = pd.read_csv(flag_bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "annot"])

    print(bed_df.head())
    print(flag_bed_df.head())
    print (data_frame.head())
    print (data_frame.tail())
    print (data_frame.shape)

    # Define the CYP2D6-D8 region boundaries
    cyp2d6_series =  ((data_frame["Chr"] == "22") & (data_frame["Position"] > 42123142) & (data_frame["Position"] < 42155251))
    #cyp2d6_series = (data_frame["Chr"] == int(22) & (data_frame["Position"] > int(42123142)) & (data_frame["Position"] < int(42155251)))

    # Identify CNV threshold exceedances
    threshold_series = (data_frame["zscore"].abs() > z_threshold) & ((data_frame["copynumber"] - 2).abs() > cn_threshold)

    report_file = os.path.join(output_dir, f"{sample_name}_CYP2D6_report_all.txt")
    results = []
    
    with open(report_file, "w") as f:
        f.write("chr\tstart_pos\tend_pos\tgene\tcnv\tmean_depth\texp_mean_depth\tmean_cn\tmean_z-score\treg.length\tpos.above\t%above\tflag\tcomment\n")

        result_list = [] 
        # Process CNVs within the CYP2D6 region
        for idx, row in data_frame.loc[cyp2d6_series].iterrows():
            chr_, pos = str("chr").replace("chr","")+str(row["Chr"]), row["Position"]
            
            #print (chr_, pos)

            # Find corresponding gene region
            gene_region = bed_df[(bed_df["chr"] == chr_) & (bed_df["start_pos"] < pos) & (bed_df["end_pos"] >= pos)]
            #gene_region = bed_df[(bed_df["chr"] == 22)]
            #print ("GENE_REGION")
            #print (gene_region)


            if gene_region.empty:
                continue

            start, end, gene_name = gene_region.iloc[0][["start_pos", "end_pos", "gene"]]
            #print("what is this")
            #print (chr_,start, end, gene_name)


            region_mask = (data_frame["Chr"] == chr_) & (data_frame["Position"] > start) & (data_frame["Position"] <= end)
            #region_mask = (data_frame["Chr"] == chr_)


            #print ("region_mask")
            #print(region_mask) 
            region_df = data_frame.loc[region_mask]
            
            #print ("REGION MASK")
            #print(region_df)

            mean_cn = region_df["copynumber"].mean()
            mean_z = region_df["zscore"].mean()
            mean_depth = region_df["Depth"].mean()
            mean_exp_depth = region_df["expDepth"].mean()

            
            #print(mean_depth,mean_z, mean_cn, mean_exp_depth)

            cnv_type = "DUP" if mean_cn > 2 else "DEL" if mean_cn < 2 else "UNDETERMINED"
            region_length = len(region_df)
            pos_above = region_df.loc[threshold_series & region_mask].shape[0]
            percent_above = pos_above / region_length if region_length else 0

            flag = "PASS" if percent_above >= prop_threshold else "FAIL"
            comment = "FAIL - Low proportion of positions above threshold. " if flag == "FAIL" else ""

            # Flag low expected depth
            if mean_exp_depth < 50:
                flag = "FAIL"
                comment += "FAIL - Cohort mean depth too low. "

            # Flag overlaps with flagged regions
            overlap_flags = flag_bed_df[
                (flag_bed_df["chr"] == chr_) &
                (flag_bed_df["start_pos"] < end) &
                (flag_bed_df["end_pos"] > start)
                ]
            
            if not overlap_flags.empty:
                comment += "WARNING - Overlap with problematic region: " + "; ".join(
                    [f"{r['chr']}:{r['start_pos']}-{r['end_pos']} {r['annot']}" for _, r in overlap_flags.iterrows()]
                )

            # Write to report file
            result_list.append([chr_, start, end, gene_name, cnv_type, mean_depth, mean_exp_depth, mean_cn, mean_z, region_length, pos_above, percent_above, flag, comment])

        unique_result_list =  [list(t) for t in set(tuple(sublist) for sublist in result_list)]
        unique_result_list = sorted(unique_result_list, key=lambda x: x[1])
        print (len(unique_result_list))

        for i in unique_result_list:
            f.write("\t".join([str(x) for x in i]) + "\n")

        #f.write(f"{chr_}\t{start}\t{end}\t{gene_name}\t{cnv_type}\t{mean_depth:.2f}\t{mean_exp_depth:.2f}\t{mean_cn:.2f}\t"
        #        f"{mean_z:.2f}\t{region_length}\t{pos_above}\t{percent_above:.2f}\t{flag}\t{comment}\n")

    print(f"CYP2D6 CNV report saved: {report_file}")

    # Load the report into a DataFrame
    cnv_regions = pd.read_csv(report_file, sep='\t')

    # Save a filtered report with only PASS flags
    filtered_report_file = os.path.join(output_dir, f"{sample_name}_CYP2D6_report_filtered.txt")
    cnv_regions_pass = cnv_regions[cnv_regions["flag"] == "PASS"]
    cnv_regions_pass.to_csv(filtered_report_file, sep='\t', index=False)
    print(f"Filtered report saved: {filtered_report_file}")

    # Create CYP2D6 vs CYP2D7 table
    element_df = cnv_regions.copy()
    #element_df["gene"] = element_df["gene"].apply(lambda x: x.split("_")[0])
    element_df["gene_new"] = element_df["gene"].apply(lambda x: x.split("/")[0])
    element_df["element"] = element_df["gene"].apply(lambda x: x.split("_")[-1])


    print (element_df["gene_new"].head() )
    print (element_df["element"].head() )


    elements = ["ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "ex8", "ex9", "REP"]
    cn_df = pd.DataFrame({"element": elements, "CYP2D6": 0.0, "CYP2D7": 0.0, "total": 0.0, "diff": 0.0, "ratio": 0.0, "CYP2D6_DEL_DUP": ""})


    for i, element in enumerate(elements):

        if element != "REP":
            print (i, element)
            d6_values = element_df[(element_df["element"] == element) & (element_df["gene_new"] == "CYP2D6")]["mean_cn"]
            d7_values = element_df[(element_df["element"] == element) & (element_df["gene_new"] == "CYP2D7")]["mean_cn"]
        else:
            print (i, element)
            d6_values = element_df[(element_df["element"] == "REP6") & (element_df["gene_new"] == "CYP2D6ds")]["mean_cn"]
            d7_values = element_df[(element_df["element"] == "REP7") & (element_df["gene_new"] == "CYP2D6-D7_REP7")]["mean_cn"]


        cn_df.at[i, "CYP2D6"]           = round(d6_values.mean(),1) if not d6_values.empty else 0.0
        cn_df.at[i, "CYP2D7"]           = round(d7_values.mean(),1) if not d7_values.empty else 0.0
        cn_df.at[i, "total"]            = round(cn_df.at[i, "CYP2D6"] + cn_df.at[i, "CYP2D7"],1)
        cn_df.at[i, "diff"]             = round(abs(cn_df.at[i, "CYP2D6"] - cn_df.at[i, "CYP2D7"]),1)
        cn_df.at[i, "ratio"]            = round(cn_df.at[i, "CYP2D6"] / cn_df.at[i, "CYP2D7"],1) if cn_df.at[i, "CYP2D7"] != 0 else np.nan
        cn_df.at[i, "CYP2D6_DEL_DUP"]   = "DUP" if cn_df.at[i, "ratio"] > 1 else ("DEL" if cn_df.at[i, "ratio"] < 1 else 0.0)

    # Save CN table
    print (cn_df)
    plot_CYP2D6_vs_CYP2D7 (cn_df, sample_name, output_dir)

    cn_table_file = os.path.join(output_dir, f"{sample_name}_CYP2D6_vs_CYP2D7_table.txt")
    cn_df.to_csv(cn_table_file, sep='\t', index=False)
    print(f"CYP2D6 vs CYP2D7 table saved: {cn_table_file}")

  

    return cnv_regions


def plot_CYP2D6_vs_CYP2D7(cn_df, sample_name, output_dir):
    # Find max and min copy numbers
    max_cn = max(
        np.sort(cn_df["CYP2D6"])[-1],
        np.sort(cn_df["CYP2D7"])[-1]
    )
    min_cn = min(
        np.sort(cn_df["CYP2D6"])[0],
        np.sort(cn_df["CYP2D7"])[0]
    )

    # Set Y-axis limits
    y_high_lim_cn = max_cn if max_cn > 4 else 4
    y_low_lim_cn = min_cn if min_cn < 0 else 0

    # Start plotting
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot CYP2D6
    ax.scatter(
        range(len(cn_df)), 
        cn_df["CYP2D6"], 
        label="CYP2D6", 
        color="black", 
        s=150, 
        alpha=1
    )

    # Plot CYP2D7
    ax.scatter(
        range(len(cn_df)), 
        cn_df["CYP2D7"], 
        label="CYP2D7", 
        color="grey", 
        s=120, 
        alpha=1
    )

    # X-ticks and labels
    ax.set_xticks(range(len(cn_df)))
    ax.set_xticklabels(cn_df["element"], rotation=90, fontsize=12)

    # Grid and axis settings
    ax.grid(True, which='major', axis='x', linestyle='--', alpha=0.7)
    ax.minorticks_on()
    ax.set_ylabel("Copy Number", fontsize=12)
    ax.set_ylim(y_low_lim_cn, y_high_lim_cn)
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', labelsize=12, direction='out')

    # Title and legend
    ax.set_title(f"{sample_name}: CYP2D6 vs CYP2D7 copy numbers", fontsize=14)
    ax.legend(fontsize=12, frameon=False)

    # Horizontal lines at 1.5 and 2.5
    ax.axhline(y=1.5, linestyle='dashed', color='grey')
    ax.axhline(y=2.5, linestyle='dashed', color='grey')

    # Save plot
    output_path = f"{output_dir}/{sample_name}_CYP2D6_vs_CYP2D7.png"
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    print (f"Saving CYP2D6 vs CYP2D7 plot as {output_path}")


def plot_all_genes(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, bed_file, output_dir):
    """
    Generates plots for all genes in the dataset.

    Parameters:
    - data_frame: Pandas DataFrame containing CNV data (Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber)
    - gene_lims: DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Sample name for labeling
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - bed_file: Path to BED file
    - output_dir: Directory to save plots
    """
    print("Generating plots for all genes...")
    

    for index, row in gene_lims.iterrows():
        chr_, start, end, gene = row["chr"], row["start"], row["end"], row["gene"]
        print(f"Plotting gene: {gene} ({chr_}:{start}-{end})")
        
        # Call the plot_region function for each gene
        plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chr_, start, end, bed_file, output_dir)

def setup_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def plot_region_with_layers(data_frame, gene_lims, sample_name, z_threshold,    cn_threshold, chr_name, start_pos, end_pos, bed_file, output_dir):

    # Log start
    """
    Creates a plot for a given region in a sample.

    Parameters:
    data_frame (pd.DataFrame): DataFrame containing the CNV data.
    gene_lims (pd.DataFrame): DataFrame containing the start and end positions of unique genes.
    sample_name (str): Name of the sample.
    z_threshold (float): Z-score threshold above or below which a marker is considered significant.
    cn_threshold (float): Copy number threshold above or below which a marker is considered significant.
    chr_name (str): Chromosome name.
    start_pos (int): Start position of the region.
    end_pos (int): End position of the region.
    bed_file (str): Path to the BED file.
    output_dir (str): Output directory for the plot.
    log_file (str): Log file to write to.

    Returns:
    None
    """

    # Read BED file
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start_pos', 'end_pos', 'annot'])

    print(bed_df.head())
    print (bed_df.shape)
    print (gene_lims.tail())
    print (gene_lims)

    # Filter region data
    this_reg_df = data_frame[(data_frame['Chr'] == chr_name) &
                             (data_frame['Position'] >= start_pos) &
                             (data_frame['Position'] <= end_pos)].reset_index(drop=True)

    this_reg_gene_lims = gene_lims[(gene_lims.iloc[:, 0] == chr_name) &
                                   (gene_lims.iloc[:, 1] < end_pos) &
                                   (gene_lims.iloc[:, 2] > start_pos)].reset_index(drop=True)

    print (this_reg_gene_lims.tail() )
    print (this_reg_gene_lims.shape)

    # Map gene limits to row indices
    gene_limits = []
    for _, row in this_reg_gene_lims.iterrows():
        start_idx = this_reg_df.index[this_reg_df['Position'] == row.iloc[1]].tolist()
        end_idx = this_reg_df.index[this_reg_df['Position'] == row.iloc[2]].tolist()

        start_idx = start_idx[0] if start_idx else 0
        end_idx = end_idx[0] if end_idx else len(this_reg_df) - 1

        gene_limits.append([start_idx, end_idx, row.iloc[3]])

    print (gene_limits)
    # Identify BED regions (REP/exon) for shading
    shaded_regions = []
    for _, row in bed_df.iterrows():
        if row['chr'] == chr_name:
            start_idx = this_reg_df.index[this_reg_df['Position'] == row['start_pos']].tolist()
            end_idx = this_reg_df.index[this_reg_df['Position'] == row['end_pos']].tolist()
            if start_idx and end_idx:
                shaded_regions.append((start_idx[0], end_idx[0], row['annot']))

    print (shaded_regions)
    # Threshold markers
    threshold_series = ((np.abs(this_reg_df['zscore']) > z_threshold) &
                        (np.abs(this_reg_df['copynumber'] - 2) > cn_threshold)).astype(int)

    color_map = mcolors.LinearSegmentedColormap.from_list('black_red', ['black', 'red'])

    # Plotting
    # fig, axs = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
    fig, axs = plt.subplots(3, 1, figsize=(16, 12))

    x_series = np.arange(len(this_reg_df))

    elements = ["dsREP","ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "ex8", "ex9","REP6","REP7"]

    # Add shaded regions (REP/exons)
    for start_idx, end_idx, annot in shaded_regions:
        
        annotation = annot.split('_')[-1]
       
        for ax in axs:
            if annotation in elements:
                ax.axvspan(start_idx, end_idx, color='lightgrey', alpha=0.75, linewidth=2, edgecolor='black', linestyle='dotted')
        for ax in axs[:3]:
            print (start_idx, end_idx, annotation)
            if annotation == 'REP6' or annotation == 'REP7' or annotation == 'dsREP':
                ax.annotate(annotation, xy=(start_idx, 1), 
                xycoords=('data', 'axes fraction'), 
                fontsize=10, color='black', 
                verticalalignment='top', 
                xytext=(0, -10), textcoords='offset points')


    # Z-score plot
    axs[0].scatter(x_series, this_reg_df['zscore'],
                   c=threshold_series, cmap=color_map, s=(threshold_series + 2) * 10)
    axs[0].axhline(0, linestyle='dashed', color='grey')
    axs[0].axhline(z_threshold, linestyle='dotted', color='grey')
    axs[0].axhline(-z_threshold, linestyle='dotted', color='grey')
    axs[0].set_ylabel("Z-score", fontsize=10)
    axs[0].set_ylim(-4, 4)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)



    # CN plot
    axs[1].scatter(x_series, this_reg_df['copynumber'],
                   c=threshold_series, cmap=color_map, s=(threshold_series + 2) * 10)
    axs[1].axhline(2, linestyle='dashed', color='grey')
    axs[1].axhline(2 + cn_threshold, linestyle='dotted', color='grey')
    axs[1].axhline(2 - cn_threshold, linestyle='dotted', color='grey')
    axs[1].set_ylim(0, max(4, np.ceil(this_reg_df['copynumber'].max())))
    axs[1].set_ylabel("Copy Number", fontsize=10)
    axs[1].set_ylim(0, 4)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    #axs[1].grid(True)

    # Depth plot
    axs[2].errorbar(x_series, this_reg_df['expDepth'],yerr=np.abs((this_reg_df['Depth'] / this_reg_df['prop']) * this_reg_df['std']),fmt='o', color='grey', label='Expected Depth', markersize=1,alpha=0.8)
  

    axs[2].scatter(x_series, this_reg_df['Depth'],
                   c=threshold_series, cmap=color_map, s=(threshold_series + 2) * 10)
    axs[2].axhline(20, linestyle='dashed', color='grey')
    axs[2].set_ylabel("Read Depth", fontsize=10)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)
    #axs[2].grid(True)

    # Add gene boundaries
    for start_idx, end_idx, gene in gene_limits:
        #print (start_idx, end_idx, gene)
        gene_names = gene.split('/')[0]
        pattern = r'^(CYP2D6|CYP2D7|CYP2D8P)$'
        if re.search(pattern, gene_names):
            label_names = gene_names
        else:
            label_names = ""
        for ax in axs:
            ax.axvline(start_idx, linestyle='solid', color='black')
            ax.axvline(end_idx, linestyle='solid', color='black')
        for ax in axs[:3]:
            ax.annotate(label_names, (start_idx, ax.get_ylim()[1]), fontsize=12, color='black')

    # Annotations
    axs[0].annotate(f"red markers = cn 2±{cn_threshold}, z-score ±{z_threshold}",
                    xy=(0, 1.15),  # 15% above the top
                    xycoords='axes fraction',
                    fontsize=12, ha='left', va='bottom')
    #                (0, axs[0].get_ylim()[1]*2), fontsize=15)                    
    axs[1].annotate("shaded areas = exons & REP-regions",
                    xy=(0, 1.15),  # 15% above the top
                    xycoords='axes fraction',
                    fontsize=12, ha='left', va='bottom')
    axs[2].annotate("grey markers = expected read depth",
                    xy=(0, 1.15),  # 15% above the top
                    xycoords='axes fraction',
                    fontsize=12, ha='left', va='bottom')

    # X-ticks as genomic positions
    tick_indices = np.linspace(0, len(this_reg_df) - 1, 10, dtype=int)
    tick_labels = [f"{chr_name}:{this_reg_df.iloc[idx]['Position']}" for idx in tick_indices]

    for ax in axs:
        ax.set_xticks(tick_indices)
        ax.set_xticklabels(tick_labels, rotation=0, fontsize=8)

    #axs[2].set_xticks(tick_indices)
    #axs[2].set_xticklabels(tick_labels, rotation=0, fontsize=8)
    axs[2].set_xlabel("Position")

    # Final plot
    fig.suptitle(f"{sample_name} {chr_name}:{start_pos}-{end_pos}", fontsize=16)
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.tight_layout(rect=[0, 0.01, 1, 0.95], h_pad=3)

    # Save plot
    output_path = os.path.join(output_dir, f"{sample_name}_{chr_name}_{start_pos}_{end_pos}.png")
    fig.savefig(output_path, dpi=300)


    plt.close(fig)

def main():
    """
    Main entry point for the script.

    This script takes in a sample name, input CNV file, input BED file, input flag BED file, z-score threshold, copy number threshold, proportion threshold, and chromosome lengths file as input and generates a report in the specified output folder.

    The report includes a CNV report, a plot of the CYP2D6 region, and plots for all genes in the dataset.

    The report directory is created if it does not exist, otherwise it is used as is.

    A log file is created in the report directory with the same name as the sample name.

    The script prints out the start time, the command line, and the end time to the log file.

    The script prints out the analysis steps to the log file.

    The script prints out the end time to the log file.
    """
    parser = argparse.ArgumentParser(description="Generate CYP2D6 CNV report from input files")
    parser.add_argument("-s","--input_sample_name", type=str, help="Prefix for output report and plots")
    parser.add_argument("-c","--input_cnv", type=str, help="Path to the input CNV file")
    parser.add_argument("-b","--input_bed", type=str, help="Path to the input BED file")
    parser.add_argument("-f","--input_flag_bed", type=str, help="Path to the input flag BED file")
    parser.add_argument("-z","--input_z_threshold", type=float, help="Threshold for the z-score")
    parser.add_argument("-cn","--input_cn_threshold", type=float, help="Threshold for the copy number deviation from 2")
    parser.add_argument("-p","--input_prop_threshold", type=float, help="Proportion threshold")
    parser.add_argument("-l","--chr_lengths", type=str, help="Path to chromosome lengths file")
    parser.add_argument("-o","--output_folder", type=str, help="Folder to store the report output")
    
    args = parser.parse_args()

    print("****************************************************")
    print(f"Running Python CYP2D6 CNV Report Script with arguments: {args}")

    # Create report directory
    report_dir = os.path.join(args.output_folder, f"{args.input_sample_name}_CYP2D6_report")
    if not os.path.exists(report_dir):
        print(f"\nCreating report directory: {report_dir}")
        os.makedirs(report_dir)
    else:
        print(f"\nUsing existing report§ directory: {report_dir}")

    log_file = os.path.join(report_dir, f"{args.input_sample_name}_CYP2D6_report_log.txt")
    print(f"\nCreating log file: {log_file}")

    # Setup logging
    setup_logging(log_file)

    logging.info("****************************************************")
    logging.info(f"Running Python generate_CYP2D6_CNV_report.py {args.input_sample_name} {args.input_cnv} {args.input_bed} {args.input_flag_bed} {args.input_z_threshold} {args.input_cn_threshold} {args.input_prop_threshold} {args.chr_lengths} {args.output_folder}")
    logging.info(f"Start time: {datetime.now()}")

    # Read input CNV file
    logging.info(f"\nReading input CNV file: {args.input_cnv}")
    df_cnv = read_cnv_file(args.input_cnv)

    # Read BED file
    logging.info(f"\nFinding gene starts and ends from input BED file: {args.input_bed}")
    ##file_name = args.input_bed
    ##df_bed = gene_starts_and_ends(file_name)
    

    # Extract gene limits
    df_gene_limits = gene_starts_and_ends(args.input_bed)

    # Create CYP2D6 report
    logging.info("\nCreating CYP2D6 report...")
    cyp2d6_regions = cyp2d6_report(df_cnv, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, args.input_prop_threshold, args.input_bed, args.input_flag_bed, report_dir)


    # Create CYP2D6 region plot
    logging.info("\nCreating CYP2D6 region plot...")
    #plot_region(df_cnv, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, "22", 42123143, 42155250, args.input_bed, report_dir)

    plot_region_with_layers(df_cnv, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, "22", 42123143, 42155250, args.input_bed, report_dir)

    # Plot all genes
    logging.info("\nPlotting all genes...")
    #plot_all_genes(df_cnv, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, args.input_bed, report_dir)

    # Final log statements
    logging.info("Analysis complete.")
    logging.info(f"End time: {datetime.now()}")

if __name__ == "__main__":
    main()
