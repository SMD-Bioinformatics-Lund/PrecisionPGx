#!/usr/bin/env python3
import argparse
import pandas as pd
import datetime
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description="PharmCat Coverage Analysis")
    parser.add_argument("--sample", type=str, help="Path to sample file")
    parser.add_argument("--depths_file", type=str, help="Path to precomputed depths file for pharmcat positions")
    parser.add_argument("--coverage_threshold", type=int, help="Minimum coverage threshold")
    parser.add_argument("--output_folder", type=str, help="Path to PharmCat reports folder")
    return parser.parse_args()


def check_coverage(sample, output_folder, coverage_threshold, depths_file):
    df = pd.read_csv(depths_file, sep="\t", names=["chr", "pos", "depth"])

    low_coverage_positions = df[df["depth"] < coverage_threshold]
    low_cov_path = os.path.join(output_folder, f"{sample}.pharmcat.pos.low.coverage.txt")
    os.makedirs(os.path.dirname(low_cov_path), exist_ok=True)
    
    with open(low_cov_path, "w") as f:
        f.write(f"## Sample: {sample}\n")
        if low_coverage_positions.empty:
            f.write("## All positions above threshold.\n")
        else:
            f.write(f"## Positions below threshold ({coverage_threshold} reads):\n")
            low_coverage_positions.to_csv(f, sep="\t", index=False, header=False)


def main():
    args = parse_arguments()
    if not args.output_folder:
        args.output_folder = os.getcwd()
    
    check_coverage(args.sample, args.output_folder, args.coverage_threshold, args.depths_file)


if __name__ == "__main__":
    main()
