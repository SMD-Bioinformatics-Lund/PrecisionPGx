#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import os
from datetime import datetime


def normgene(chr, pos):
    """Determine whether a genomic position should be used for normalization"""
    return chr != "chrX" and pos > 0


def read_depth_file(depth_file):
    """Read the depth file into a DataFrame"""
    return pd.read_csv(depth_file, sep="\t", names=["Chr", "Position", "Depth"])


def calculate_cnv(depth_df):
    """Calculate CNV metrics for a single sample"""
    # Normalize read depths
    norm_region = depth_df[
        (depth_df["Chr"] != ["chrX", "X"]) & (depth_df["Position"] > 0)
    ]
    norm_sum = norm_region["Depth"].sum()
    depth_df["prop"] = depth_df["Depth"] / norm_sum

    # Compute mean and standard deviation per position
    # TODO: Doing within the sample and not across multiple samples at the moment
    mean_prop = depth_df["prop"].mean()
    std_prop = depth_df["prop"].std()

    # Compute expected depth, z-score, and CNV
    depth_df["mean"] = mean_prop
    depth_df["std"] = std_prop
    depth_df["expDepth"] = depth_df["Depth"] / depth_df["prop"] * mean_prop
    depth_df["zscore"] = (depth_df["prop"] - mean_prop) / std_prop
    depth_df["copynumber"] = 2 * (depth_df["prop"] / mean_prop)

    return depth_df


def save_results(sample_id, cnv_df, output_file):
    """Save CNV results to a specified output file without creating directories"""
    cnv_df.to_csv(output_file, index=False)

    log_file = "calculate_CNV_log.txt"
    with open(log_file, "a") as log:
        log.write(
            f"{datetime.now()} - Processed sample {sample_id} and saved results to {output_file}\n"
        )

    print(f"CNV results saved for {sample_id} at {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate CNV from depth file.")
    parser.add_argument("-s", "--sample_id", type=str, help="Sample ID", required=True)
    parser.add_argument(
        "-d", "--depth_file", type=str, help="Path to depth file", required=True
    )
    parser.add_argument("-b", "--target_bed", type=str, help="Path to target BED file")
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        help="Path to output file (CSV format)",
        required=True,
    )

    args = parser.parse_args()

    # Read and process depth file
    depth_df = read_depth_file(args.depth_file)
    cnv_df = calculate_cnv(depth_df)

    # Save results
    save_results(args.sample_id, cnv_df, args.output_file)

    print("CNV analysis complete.")
