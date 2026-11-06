#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import os
import sys

def calculate_sample_cn(depths_df, means_df):
    # Normalize read depths and add a proportion data field to the data frame
    norm_region = depths_df[ (depths_df["Chr"] != "chrX") & (depths_df["Position"] > 0) ]
    norm_sum = norm_region["Depth"].sum()
    depths_df["NormSum"] = norm_sum
    depths_df["Prop"] = depths_df["Depth"] / norm_sum

    #---------------------------------------------
    #Test with separate normalization for chr22 and the rest of the panel
    #norm_region_excl22 = depths_df[ (depths_df["Chr"] != "chrX") & (depths_df["Position"] > 0) & (depths_df["Chr"] != "chr22") ]
    #norm_sum_excl22 = norm_region_excl22["Depth"].sum()
    #depths_df.loc[depths_df["Chr"] != "chr22", "NormSum"] = norm_sum_excl22
    #depths_df.loc[depths_df["Chr"] != "chr22", "Prop"] = depths_df["Depth"] / norm_sum_excl22 
    #-------------------------------------------------
    
    #Filter the input sample dataframe. Calculations can only be made in positions included in the input mean depths file.
    #Add a temporary column combining chr and position in both dataframes
    depths_df["ChrPosition"] = depths_df["Chr"] + ":" + depths_df["Position"].astype(str)
    means_df["ChrPosition"] = means_df["Chr"] + ":" + means_df["Position"].astype(str)
    #Create a new dataframe for the copy number calculations keeping only the positions included in the input mean depths file
    cn_df = depths_df.loc[depths_df["ChrPosition"].isin(means_df["ChrPosition"])].reset_index(drop=True) 
    #In the case of using a mean depths file that includes other regions than the depths calculated for the sample, also filter the input mean depths
    overlap_means_df = means_df.loc[means_df["ChrPosition"].isin(cn_df["ChrPosition"])].reset_index(drop=True)

    #sample_cn_df and filtered_mean_prop_df should now hold the same positions
    if cn_df.empty:
        sys.stdout.write(f'There was no overlap between the sample depths and mean prop depths. Output will be empty.\n')
    elif not cn_df["ChrPosition"].equals(overlap_means_df["ChrPosition"]):
        sys.stdout.write(f'Sample depths and mean depths were not matched as expected. Output will be empty.\n')
        #Set the copy number frame to an empty data frame
        cn_df = pd.DataFrame(list(depths_df))
    else:
        #Calculate expected depth, copy number and z-score in each position for all samples
        #Add the MeanPropDepth and StDevPropDepth to the DataFrame (for transparency in calculations)
        cn_df["AppliedMeanPropDepth"] = overlap_means_df["MeanPropDepth"]
        cn_df["AppliedStDevPropDepth"] = overlap_means_df["StDevPropDepth"]
        #Calculate expected depth: NormSum * AppliedMeanPropDepth
        cn_df["ExpectedDepth"] = cn_df["NormSum"] * cn_df["AppliedMeanPropDepth"]
        #Calculate copy number (assuming the normal copy number is two): 2 * (PropDepth / AppliedMeanPropDepth)
        cn_df["CopyNumber"] = 2 * (cn_df["Prop"] / cn_df["AppliedMeanPropDepth"])
        #Calculate z-score as a measure of deviation: (PropDepth - AppliedMeanPropDepth) / AppliedStDevPropDepth
        cn_df["ZScore"] = (cn_df["Prop"] - cn_df["AppliedMeanPropDepth"]) / cn_df["AppliedStDevPropDepth"]

    #Remove the temporary column from the sample copy number dataframe
    results_df = cn_df.drop(["ChrPosition"], axis=1)

    return results_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate copy numbers from depth file using the CNVZ algorithm and precalculated expected depths.")
    parser.add_argument(
            "-s", "--sampleid", type=str, help="Sample id. Used to name output file: <sampleid.csv>", required=True
    )
    parser.add_argument(
            "-d", "--depths", type=str, help="Path to depth file: tab separated chr/pos/depth", required=True
    )
    parser.add_argument(
            "-m", "--mean_depths", type=str, help="Path to file with mean proportional depth and stdev for each position, to be used in the calculation. Comma-separated Chr,Position,MeanPropDepth,StDevPropDepth ", required=True
    )
    parser.add_argument(
        "-o", "--output_path", type=str, help="Path to output folder, will be created if it does not exist", required=True
    )

    args = parser.parse_args()

    # Create output folder if it does not exist
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    sys.stdout.write(f"CNV calculations using the CNVZ algorithm\nInput sample id: {args.sampleid}\nInput depths: {args.depths}\nInput mean depths: {args.depths}\nOutput path: {args.output_path}\n")

    # Read input files to data frames
    input_depths_df = pd.read_csv(args.depths, sep="\t", names=["Chr", "Position", "Depth"])
    input_means_df = pd.read_csv(args.mean_depths, sep=",", float_precision='round_trip')

    sys.stdout.write(f"Calculate copy number in each position based on input mean values.\n")
    
    # Process and calculate copy numbers for input sample
    sample_results_df = calculate_sample_cn(input_depths_df, input_means_df)

    # Save results to file
    sample_results_file = args.output_path + "/" + args.sampleid + ".cnvz.csv"
    sample_results_df.to_csv(sample_results_file, index=False)

    sys.stdout.write(f"Results saved as file: {sample_results_file}\nDone.\n")
