#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
import os

def calculate_means(depthfiles):
    #Create a list to hold all normalized depths
    norm_depths = []

    #Read depth data into dataframe and calculate the sample normalized depth (proportion = Prop) 
    for i in depthfiles:
        # Read the depth file as a dataframe
        in_depth_df = pd.read_csv(i, sep="\t", names=["Chr", "Position", "Depth"])
        #Add a temporary column combining chr and position
        in_depth_df["ChrPosition"] = in_depth_df["Chr"] + ":" + in_depth_df["Position"].astype(str)
        if i == depthfiles[0]:
            #The first input file is used as reference, so that mean values will only be calculated for sites included in this file
            ref_pos_df = in_depth_df
            depth_df = in_depth_df
        else:
            #Filter the dataframe keeping only the positions included in the first input depths file
            depth_df = in_depth_df.loc[in_depth_df["ChrPosition"].isin(ref_pos_df["ChrPosition"])].reset_index(drop=True)

        # Create a subset of the data frame with the regions used for normalization (excluding sex chromosomes and mitochondria)
        norm_region = depth_df[ (depth_df["Chr"] != "chrX") & (depth_df["Position"] > 0) ]
        # Add the normalized depths to the list
        norm_depths.append(depth_df["Depth"] / norm_region["Depth"].sum())
        
        #---------------------------------------------
        #Test with separate normalization for chr22 and the rest of the panel
        #norm_region_excl22 = depth_df[ (depth_df["Chr"] != "chrX") & (depth_df["Position"] > 0) & (depth_df["Chr"] != "chr22") ]
        # Temporary dataframe for the normalized depths of this sample
        #this_norm_depth_df = depth_df[['Chr', 'Position']].copy()
        #Add normalized depths, normalized with the panel sum excluding chr22
        #this_norm_depth_df["NormDepth"] = depth_df["Depth"] / norm_region_excl22["Depth"].sum()
        #For chr22 add the normalized depths using the complete panelsum
        #this_norm_depth_df.loc[this_norm_depth_df["Chr"] == "chr22", "NormDepth"] = depth_df["Depth"] / norm_region["Depth"].sum()
        # Add the normalized depths to the list
        #norm_depths.append(this_norm_depth_df["NormDepth"])
        #-------------------------------------------------

    #Calculate the mean normalized depth across all samples for each position
    #Convert the list of normalized depths to a dataframe to facilitate calculation of mean values
    #Using np.transpose so that positions are rows and samples are columns 
    norm_depths_df = pd.DataFrame(np.transpose(norm_depths))

    #Create a new dataframe to hold chr, pos and mean normalized depth across all samples
    mean_prop_df = pd.DataFrame(data=ref_pos_df, columns=["Chr", "Position"])
    mean_prop_df["MeanPropDepth"] = norm_depths_df.mean(axis=1)
    mean_prop_df["StDevPropDepth"] = norm_depths_df.std(axis=1)

    return mean_prop_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate mean depths from depth files.")
    parser.add_argument(
            "-d", "--depths", type=str, help="List of depth files", required=True
    )
    parser.add_argument(
            "-p", "--output_prefix", type=str, help="output_file_prefix", required=True
    )
    parser.add_argument(
        "-o", "--output_path", type=str, help="Path to output folder, will be created if it does not exist", required=True
    )

    sys.stdout.write(f"Reading input.\n")
    args = parser.parse_args()
    
    sys.stdout.write(f"Output path: {args.output_path}\n")
    sys.stdout.write(f"Depth files: {args.depths}\n")

    #Create output folder if it does not exist
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # Read input file list
    depthfiles = list(args.depths.split(","))
    #depthfiles = list(args.depths.replace("[","").replace("]","").split(","))
    
    # Process and calculate copy numbers for all listed depth files
    sys.stdout.write(f"Calculate cohort mean values.\n")
    means_df = calculate_means(depthfiles)

    sys.stdout.write(f"Saving output files.\n")
    #Save the mean values applied in the calculations
    means_df.to_csv(args.output_path + "/" + args.output_prefix + ".mean_proportional_depths.csv", index=False)

    sys.stdout.write(f"Done.\n")
