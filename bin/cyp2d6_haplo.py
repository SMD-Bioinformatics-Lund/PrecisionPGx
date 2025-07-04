#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import csv
import re

# Define argument parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="CYP2D6 diplotyping from variant VCF file.")
    parser.add_argument("--sample_id", type=str, help="Sample ID")
    parser.add_argument("--input_vcf", type=str, help="Path to input VCF file")
    parser.add_argument("--output_dir", type=str, help="Output directory")
    parser.add_argument("--idx", type=int, default=4, help="Index for analysis, default is 4")
    return parser.parse_args()

# Dictionary for identification of non-identical gene duplication
DICT4 = {
    "0202": "*4x1/*4x1",
    "1202": "*4x1/[*68+*4]x1",
    "1203": "*4x1/[*68+*4]x1",
    "1102": "*4x1/*10x1",
    "1111": "*4x1/*(other)x1",
    "2204": "[*68+*4]x1/[*68+*4]x1",
    "2103": "[*68+*4]x1/*10x1",
    "2112": "[*68+*4]x1/*(other)x1",
    "2113": "[*68+*4]x1/*(other)x1",
    "1112": "[*68+*4]x1/*(other)x1",
    "1212": "*4x2/*(other)x1"
}

# Function to read and reformat VCF file
def reformat_vcf(vcf_file):
    reformatted_data = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chrom, pos, id_, ref, alt, qual, filter_, info, format_, key_attributes = cols[:10]
            if chrom == "chr22" and 42126310 <= int(pos) <= 42132392:
                ad = key_attributes.split(":")[1]
                gt = key_attributes.split(":")[0]
                if gt == "0|1":
                    gt = "0/1"
                elif gt == "1|1":
                    gt = "1/1"
                reformatted_data.append(f"{pos} {pos} {ref} {alt} {gt} {ad}")
    return reformatted_data

# Function to initiate DataFrame M1
def initiate_m1(reformatted_vcf, dfgroup):
    M1 = pd.read_csv("Mtemplate2.csv")
    for i, entry in enumerate(reformatted_vcf):
        pos, _, ref, alt, gt, ad = entry.split(" ")
        M1.loc[i, 'Start'] = int(pos)
        M1.loc[i, 'Stop'] = int(pos)
        M1.loc[i, 'Ref'] = ref
        M1.loc[i, 'Alt'] = alt
        M1.loc[i, 'GT'] = gt
        M1.loc[i, 'AD'] = ad
        regex_pattern = f"{pos} {pos} {ref} {alt}"
        matching_rows = dfgroup[dfgroup.iloc[:, 1].str.contains(regex_pattern, regex=True, na=False)]
        for _, row in matching_rows.iterrows():
            M1.loc[i, row.iloc[0]] = 1
    M1 = M1[M1['Start'] != 11111111]  # Drop redundant rows
    return M1

# Function to analyze M1
def analyze_m1(M1, dfgroup):
    results = []
    for i in range(dfgroup.shape[0]):
        total_variants = M1.iloc[:, i + 8].sum()
        expected_variants = dfgroup.iloc[i, 1].count(",")
        haplotype = dfgroup.iloc[i, 0]
        if total_variants == expected_variants and "s" not in haplotype:
            results.append(f"*{haplotype} MATCH ({total_variants}/{expected_variants})")
        elif total_variants / expected_variants < 1:
            results.append(f"*{haplotype} Excluded ({total_variants}/{expected_variants})")
        elif total_variants == expected_variants and "s" in haplotype:
            results.append(f"Candidate suballele *{haplotype}")
    return results

if __name__ == "__main__":
    args = parse_arguments()
    print("=================")
    print(args.sample_id)
    print("idx:", args.idx)
    
    dfgroup = pd.read_csv("CYP2D6.haplotypes.grouped.modified.txt", sep="\t", header=None)
    reformatted_vcf = reformat_vcf(args.input_vcf)
    M1 = initiate_m1(reformatted_vcf, dfgroup)
    results = analyze_m1(M1, dfgroup)
    for result in results:
        print(result)
