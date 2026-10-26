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
import gzip

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

# Read input VCF file
def get_vafs(vcf_file, sample_name, output_dir):
    """
    Reads a gzipped VCF file and creates a Pandas DataFrame with vafs.

    Parameters:
    vcf_file (str): Path to the VCF file.

    Returns:
    pd.DataFrame: DataFrame containing the vaf data.
    """

    vcfin = gzip.open(vcf_file, 'rt')
    
    vaf_list = [] 
    for line in vcfin:
        #Skip header lines
        if not("#" in line):
            #Create columns array
            thiscolumns = line.strip().split("\t")
            #Create values array
            thisvalues = thiscolumns[9].split(":")

            #Don't consider calls with GT = 0/0
            if not(thisvalues[0] == "0/0"):
                #Create FORMAT array
                thisformat = thiscolumns[8].split(":")
                #Get the index for the DP and AD fields
                dpidx = 0
                adidx = 0
                i = 0
                #Iterate the FORMAT fields
                for f in thisformat:
                    if f == "DP":
                        dpidx = i
                    if f == "AD":
                        adidx = i
                    i += 1
                thisaltad = thisvalues[adidx].split(",")
                #Add values to list
                vaf_list.append((thiscolumns[0], thiscolumns[1], thiscolumns[3], thiscolumns[4], thiscolumns[5], thiscolumns[6], thiscolumns[8], thiscolumns[9], int(thisaltad[-1])/int(thisvalues[dpidx])))

    #Save vaf list
    vaf_file = os.path.join(output_dir, f"{sample_name}.vafs.txt")

    # Convert to DataFrame
    vaf_df = pd.DataFrame(vaf_list, columns=["chr", "pos", "ref", "alt", "qual", "filter", "format", "sample", "vaf"])
    
    #Debug prints
    #print(vaf_df.dtypes)

    #Change data types
    vaf_df["pos"] = pd.to_numeric(vaf_df["pos"], errors='coerce')
    vaf_df["vaf"] = pd.to_numeric(vaf_df["vaf"], errors='coerce')
    
    #Debug prints
    #print(vaf_df.dtypes)

    #Save vaf list
    vaf_df.to_csv(vaf_file, sep="\t", index=False)
    
    return vaf_df

# Get the ranges of unique genes from the bed file
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
            
            #Convert to 1-based
            #this_start = row["start"] + 1  # Convert 0-based start to 1-based
            
            #Keep 0-based
            this_start = row["start"]

        # If a new gene is encountered, store previous gene details
        if this_line_gene_name != this_gene:
            gene_limits.append([df_bed.loc[index - 1, "chr"], this_start, df_bed.loc[index - 1, "end"], this_gene])
            
            # Reset values for the new gene
            this_gene = this_line_gene_name
            
            #Convert to 1-based
            #this_start = row["start"] + 1
            
            #Keep 0-based
            this_start = row["start"]
            
        # If it's the last row, store the final gene details
        if index == len(df_bed) - 1:
            gene_limits.append([row["chr"], this_start, row["end"], this_gene])

    return pd.DataFrame(gene_limits, columns=["chr", "start", "end", "gene"])

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

def plot_region_with_layers(data_frame, vafs, gene_lims, sample_name, z_threshold, cn_threshold, chr_name, start_pos, end_pos, bed_file, plotlabel, output_dir):

    # Log start
    """
    Creates a plot for a given region in a sample.

    Parameters:
    data_frame (pd.DataFrame): DataFrame containing the CNV data.
    vafs (pd:DataFrame): DataFrame containing variant info including vafs
    gene_lims (pd.DataFrame): DataFrame containing the start and end positions of unique genes.
    sample_name (str): Name of the sample.
    z_threshold (float): Z-score threshold above or below which a marker is considered significant.
    cn_threshold (float): Copy number threshold above or below which a marker is considered significant.
    chr_name (str): Chromosome name.
    start_pos (int): Start position of the region.
    end_pos (int): End position of the region.
    bed_file (str): Path to the BED file.
    plotlabel (str): Used in plot headings and output file names
    output_dir (str): Output directory for the plot.
    log_file (str): Log file to write to.

    Returns:
    None
    """

    # Read BED file
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start_pos', 'end_pos', 'annot'])

    #Debug prints
    #print ("bed:")
    #print(bed_df.head())
    #print (bed_df.shape)
    #print ("gene_lims")
    #print (gene_lims.tail())
    #print (gene_lims)
    #print(vafs.dtypes)
    #print(data_frame.dtypes)

    # Filter region data
    this_reg_df = data_frame[(data_frame['Chr'] == chr_name) &
                             (data_frame['Position'] >= start_pos) &
                             (data_frame['Position'] <= end_pos)].reset_index(drop=True)

    #Set copy number and Z-score to NaN for positions where the Expected depth is not greater than 50, as these positions are not reliable
    #For example if expected depth is 0 and observed depth is 1, the calculated copy number will be infinite
    this_reg_df["CopyNumber"] = this_reg_df[this_reg_df.ExpectedDepth > 50]["CopyNumber"]
    this_reg_df["ZScore"] = this_reg_df[this_reg_df.ExpectedDepth > 50]["ZScore"]
    
    this_reg_vaf_df = vafs[(vafs['chr'] == chr_name) &
                             (vafs['pos'] >= start_pos) &
                             (vafs['pos'] <= end_pos)].reset_index(drop=True)

    #Debug prints
    #print ("this_reg_df:")
    #print (this_reg_df)

    this_reg_gene_lims = gene_lims[(gene_lims.iloc[:, 0] == chr_name) &
                                   (gene_lims.iloc[:, 1] < end_pos) &
                                   (gene_lims.iloc[:, 2] > start_pos)].reset_index(drop=True)

    #Debug prints
    #print ("this_reg_gene_lims:")
    #print (this_reg_gene_lims )
    #print (this_reg_gene_lims.shape)

    # Map gene limits to row indices
    gene_limits = []
    for _, row in this_reg_gene_lims.iterrows():
        start_idx = this_reg_df.index[this_reg_df['Position'] == row.iloc[1]].tolist()
        end_idx = this_reg_df.index[this_reg_df['Position'] == row.iloc[2]].tolist()

        start_idx = start_idx[0] if start_idx else 0
        end_idx = end_idx[0] if end_idx else len(this_reg_df) - 1

        gene_limits.append([start_idx, end_idx, row.iloc[3]])

    #Debug prints
    #print ("gene_limits")
    #print (gene_limits)
    
    # Identify BED regions (REP/exon) for shading
    bed_regions = []
    for _, row in bed_df.iterrows():
        if row['chr'] == chr_name:
            start_idx = this_reg_df.index[this_reg_df['Position'] == row['start_pos']].tolist()
            end_idx = this_reg_df.index[this_reg_df['Position'] == row['end_pos']].tolist()
            
            #Debug prints
            #print("bed_df row:")
            #print(row)
            #print("start_idx:")
            #print(start_idx)
            #print("end_idx:")
            #print(end_idx)

            if start_idx and end_idx:
                bed_regions.append((start_idx[0], end_idx[0], row['annot']))

    #Debug prints
    #print ("bed_regions")
    #print (bed_regions)
    
    #Figure with read depth/copy number/Z-score for every position in the region
    #----------------------------------------------------------------------------

    # Threshold markers
    threshold_series = ((np.abs(this_reg_df['ZScore']) > z_threshold) &
                        (np.abs(this_reg_df['CopyNumber'] - 2) > cn_threshold)).astype(int)

    color_map = mcolors.LinearSegmentedColormap.from_list('black_red', ['black', 'red'])

    # Plotting
    fig, axs = plt.subplots(4, 1, figsize=(16, 12))

    x_series = np.arange(len(this_reg_df))
    vaf_x_series = [this_reg_df.loc[this_reg_df['Position']==x].index[0] for x in this_reg_vaf_df['pos']]

    elements = ["REP", "ex", "ts"]
    ignoreelements = [ "dsex" ]

    # Add shaded regions (REP/exons)
    for start_idx, end_idx, annot in bed_regions:
        
        annotation = annot.split('_')[-1]
       
        for ax in axs:
            #if annotation in elements:
            if any(w in annotation for w in elements) and not(any(w in annotation for w in ignoreelements)):
                ax.axvspan(start_idx, end_idx, color='lightgrey', alpha=0.75, linewidth=2, edgecolor='black', linestyle='dotted')
        for ax in axs[:3]:
            #Debug prints
            #print (start_idx, end_idx, annotation)

            if any(w in annotation for w in elements) and not(any(w in annotation for w in ignoreelements)):
                ax.annotate(annotation, xy=(start_idx, 1), 
                xycoords=('data', 'axes fraction'), 
                fontsize=10, color='black', 
                verticalalignment='top', 
                xytext=(0, -10), textcoords='offset points', rotation=90)
        
        #Place annotations at the bottom for the VAF graph
        for ax in axs[3:]:
            #Debug prints
            #print (start_idx, end_idx, annotation)

            if any(w in annotation for w in elements) and not(any(w in annotation for w in ignoreelements)):
                ax.annotate(annotation, xy=(start_idx, 0.1),
                xycoords=('data', 'axes fraction'),
                fontsize=10, color='black',
                verticalalignment='bottom',
                xytext=(0, -10), textcoords='offset points', rotation=90)

    # Z-score plot
    axs[0].scatter(x_series, this_reg_df['ZScore'],
                   c=threshold_series, cmap=color_map, s=(threshold_series + 2) * 10)
    axs[0].axhline(0, linestyle='dashed', color='grey')
    axs[0].axhline(z_threshold, linestyle='dotted', color='grey')
    axs[0].axhline(-z_threshold, linestyle='dotted', color='grey')
    axs[0].set_ylabel("z-score", fontsize=12)
    axs[0].set_ylim(min(-4, np.floor(this_reg_df['ZScore'].min())), max(4, np.ceil(this_reg_df['ZScore'].max())))
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    # CN plot
    axs[1].scatter(x_series, this_reg_df['CopyNumber'],
                   c=threshold_series, cmap=color_map, s=(threshold_series + 2) * 10)
    axs[1].axhline(2, linestyle='dashed', color='grey')
    axs[1].axhline(2 + cn_threshold, linestyle='dotted', color='grey')
    axs[1].axhline(2 - cn_threshold, linestyle='dotted', color='grey')
    axs[1].set_ylim(0, max(4, np.ceil(this_reg_df['CopyNumber'].max())))
    axs[1].set_ylabel("copy number", fontsize=12)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    # Depth plot
    axs[2].errorbar(x_series, this_reg_df['ExpectedDepth'],yerr=np.abs((this_reg_df['Depth'] / this_reg_df['Prop']) * this_reg_df['AppliedStDevPropDepth']),fmt='o', color='grey', label='expected read depth', markersize=1,alpha=0.8)
    axs[2].scatter(x_series, this_reg_df['Depth'], label='sample read depth', c=threshold_series, cmap=color_map, s=(threshold_series + 2) * 10, zorder=2.5)
    axs[2].legend()
    axs[2].axhline(20, linestyle='dashed', color='grey')
    axs[2].set_ylabel("read depth", fontsize=12)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)

    # VAF plot
    axs[3].scatter(vaf_x_series, this_reg_vaf_df['vaf'], c='black', s=20)
    axs[3].axhline(0.5, linestyle='dashed', color='grey')
    axs[3].axhline(1.0, linestyle='dashed', color='grey')
    axs[3].axhline(0.66, linestyle='dashed', color='lightgrey')
    axs[3].axhline(0.33, linestyle='dashed', color='lightgrey')
    axs[3].set_ylabel("variant allele frequency", fontsize=12)
    axs[3].set_ylim(0, 1.05)
    axs[3].spines['top'].set_visible(False)
    axs[3].spines['right'].set_visible(False)

    # Add gene boundaries
    for start_idx, end_idx, gene in gene_limits:
        #Don't draw boundaries for downstream, upstream or rs-SNPs
        ignoregenes=["ds", "us", "rs", "-D"]
        if any(w in gene for w in ignoregenes):
            label_names = ""
        else:
            label_names = gene
        for ax in axs:
            ax.axvline(start_idx, linestyle='solid', color='black')
            ax.axvline(end_idx, linestyle='solid', color='black')
        for ax in axs[:4]:
            ax.annotate(label_names, (start_idx, ax.get_ylim()[1]), fontsize=12, color='black')

    # Annotations
    axs[0].annotate(f"red markers indicate deviation from cn and z-score thresholds (dotted lines cn 2±{cn_threshold}, z-score ±{z_threshold})",
                     xy=(0, 1.10),  # 10% above the top
                     xycoords='axes fraction',
                     fontsize=10, ha='left', va='bottom')

    # X-ticks as genomic positions
    tick_indices = np.linspace(0, len(this_reg_df) - 1, 10, dtype=int)
    tick_labels = [f"{chr_name}:{this_reg_df.iloc[idx]['Position']}" for idx in tick_indices]

    for ax in axs:
        ax.set_xticks(tick_indices)
        ax.set_xticklabels(tick_labels, rotation=0, fontsize=8)

    # X-axis label only for bottom panel
    axs[3].set_xlabel("position", fontsize=12)

    # Final plot
    fig.suptitle(f"{sample_name} {plotlabel} {chr_name}:{start_pos}-{end_pos}", fontsize=16)
    fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], h_pad=2)

    # Save plot
    output_path = os.path.join(output_dir, f"{sample_name}.{plotlabel}.{chr_name}_{start_pos}_{end_pos}.png")
    fig.savefig(output_path, dpi=300)

    plt.close(fig)

    #Figure with mean read depth/copy number for all annotated elements in the bed file
    #----------------------------------------------------------------------------
    #Calculate mean cn for all bed elements in the region
    elements = ["REP", "ex", "ts"]
    ignoreelements = [ "dsex" , "dsREP"]

    element_means = []
    means_report = []
    for start_idx, end_idx, annot in bed_regions:
        if any(w in annot for w in elements) and not(any(w in annot for w in ignoreelements)):
            element_means.append((start_idx, end_idx, annot, np.nanmean(this_reg_df['Depth'].loc[start_idx:end_idx]), np.nanmean(this_reg_df['ExpectedDepth'].loc[start_idx:end_idx]), np.nanmean(this_reg_df['CopyNumber'].loc[start_idx:end_idx])))
            means_report.append((chr_name, this_reg_df['Position'][start_idx], this_reg_df['Position'][end_idx], annot, np.round(np.nanmean(this_reg_df['Depth'].loc[start_idx:end_idx]), 2), np.round(np.nanmean(this_reg_df['ExpectedDepth'].loc[start_idx:end_idx]), 2), np.round(np.nanmean(this_reg_df['CopyNumber'].loc[start_idx:end_idx]), 2)))

    #Do not try to create mean value reports or figures if there was no annotations with the above elements
    if len(element_means) == 0:
        return

    #Debug prints
    #print("element_means")
    #print(element_means)

    #Save text report
    report_file = os.path.join(output_dir, f"{sample_name}.{plotlabel}.exon_means.txt")

    # Convert to DataFrame
    means_report_df = pd.DataFrame(means_report, columns=["chr", "start_pos", "end_pos", "annotation", "mean_depth", "exp_mean_depth", "mean_copy_number"])

    # Save full CNV report
    means_report_df.to_csv(report_file, sep="\t", index=False)


    #For CYP2D6, highlight elements that deviate more than 25% from the calculated copy number of CYP2D8
    if ("CYP2D6" in [x[2] for x in gene_limits]) and ("CYP2D8P" in [x[2] for x in gene_limits]):
        cn_colors = []
        #Mean copy number of CYP2D8
        cyp2d8cn = [x[5] for x in element_means][[x[2] for x in element_means].index("CYP2D8P_ENSG00000226450_ts")]
        for x in element_means:
            if ((x[5] > (1.25 * cyp2d8cn)) or (x[5] < (0.75 * cyp2d8cn))):
                cn_colors.append("red")
            else:
                cn_colors.append("black")
    else:
        cn_colors = ["black"]*len(element_means)
    
    # Mean value plotting
    meanfig, meanaxs = plt.subplots(2, 1, figsize=(16, 12))
    meanx_series = np.arange(len(element_means)) 
    
    #Debug prints
    #print("length meanxseries:")
    #print(len(meanx_series))
    #print("length element_means:")
    #print(len(element_means)) 
    #print("length sliced element_means:")
    #print(len([x[2] for x in element_means]))
    #print("shape meanxseries:")
    #print(np.shape(meanx_series))
    #print("shape sliced element_means:")
    #print(np.shape([x[2] for x in element_means]))
    #print("sliced elementmeans")
    #print([x[2] for x in element_means])
    #print("sliced elementmeans no square brackets")
    #print((x[2] for x in element_means))

    # CN plot
    meanaxs[0].scatter(meanx_series, [x[5] for x in element_means], s=400, c=cn_colors,edgecolors='none', zorder=2.5)
    meanaxs[0].axhline(2, linestyle='dashed', color='grey')
    if ("CYP2D6" in [x[2] for x in gene_limits]) and ("CYP2D8P" in [x[2] for x in gene_limits]):    
        meanaxs[0].axhline(cyp2d8cn, linestyle='dotted', color='grey')
    meanaxs[0].set_ylim(0, max(4, np.ceil(max([x[5] for x in element_means]))))
    meanaxs[0].set_ylabel("mean copy number", fontsize=14)
    meanaxs[0].spines['top'].set_visible(False)
    meanaxs[0].spines['right'].set_visible(False)
    meanaxs[0].set_xticks([])

    # Depth plot
    meanaxs[1].errorbar(meanx_series, [x[4] for x in element_means],fmt='o', color='grey', label='expected depth', markersize=20,alpha=0.5,zorder=2.5)
    meanaxs[1].scatter(meanx_series, [x[3] for x in element_means], s=400, label="sample depth", c='black',edgecolors='none',zorder=2.8)
    meanaxs[1].legend()
    meanaxs[1].axhline(20, linestyle='dashed', color='grey')
    meanaxs[1].set_ylabel("mean read depth", fontsize=14)
    meanaxs[1].spines['top'].set_visible(False)
    meanaxs[1].spines['right'].set_visible(False)
    meanaxs[1].set_xticks([])
    meanaxs[1].set_xlabel("annotated exons", fontsize=14)

    # Add shaded regions (REP/exons)
    thisn = 0
    for start_idx, end_idx, annot, depth, expdepth, meancn in element_means:
        annotation = annot.split('_')[-1]
        for meanax in meanaxs:
            if any(w in annotation for w in elements) and not(any(w in annotation for w in ignoreelements)):
                meanax.axvline(thisn, linestyle='solid', color='lightgrey', linewidth=24, alpha=0.75)
        for meanax in meanaxs[:2]:
            if any(w in annotation for w in elements) and not(any(w in annotation for w in ignoreelements)):
                    meanax.annotate(annotation, xy=(thisn, 1),
                            xycoords=('data', 'axes fraction'),
                            fontsize=12, color='black',
                            verticalalignment='top',
                            xytext=(0, -10), textcoords='offset points', rotation=90)
        thisn += 1

    # Add gene boundaries
    for start_idx, end_idx, gene in gene_limits:
        #Don't draw boundaries for downstream, upstream or rs-SNPs
        ignoregenes=["ds", "us", "rs", "-D"]
        #Extract the annotated genes in the element_means list 
        mean_genes = [g[0] for g in [a.split("_") for a in [x[2] for x in element_means]]]
        if not(any(w in gene for w in ignoregenes)):
            meanlabel_names = gene
            #Index in element_means for the first occurrence of gene
            meanidx = mean_genes.index(gene)
            for meanax in meanaxs:
                meanax.axvline(meanidx, linestyle='solid', color='black')
                meanax.axvline(meanidx, linestyle='solid', color='black')
            for meanax in meanaxs[:3]:
                meanax.annotate(meanlabel_names, (meanidx, meanax.get_ylim()[1]), fontsize=14, color='black')
    
    #Marker annotation
    if ("CYP2D6" in [x[2] for x in gene_limits]) and ("CYP2D8P" in [x[2] for x in gene_limits]):
        markerlimit = np.round(cyp2d8cn, 2)
        meanaxs[0].annotate(f"red markers indicate deviation greater than 25% from CYP2D8p baseline (dotted line cn {markerlimit})",
                        xy=(0, 1.05),  # 5% above the top
                        xycoords='axes fraction',
                        fontsize=10, ha='left', va='bottom')

    # Final plot
    meanfig.suptitle(f"{sample_name} {plotlabel} exon means", fontsize=16)

    meanfig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], h_pad=2)

    # Save plot
    meanoutput_path = os.path.join(output_dir, f"{sample_name}.{plotlabel}.exon_means.png")
    meanfig.savefig(meanoutput_path, dpi=300)

    plt.close(meanfig)

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
    parser.add_argument("-cnv","--input_cnv", type=str, help="Path to the input CNV file")
    parser.add_argument("-vcf","--input_vcf", type=str, help="Path to the input VCF.gz file")
    parser.add_argument("-bed","--input_bed", type=str, help="Path to the input BED file")
    parser.add_argument("-z","--input_z_threshold", type=float, help="Threshold for the z-score")
    parser.add_argument("-cn","--input_cn_threshold", type=float, help="Threshold for the copy number deviation from 2")
    parser.add_argument("-p","--input_prop_threshold", type=float, help="Proportion threshold")
    parser.add_argument("-rc","--report_cyp2d6", type=str, help="true/false")
    parser.add_argument("-g","--selected_genes", type=str, help="comma-separated list of genes")
    parser.add_argument("-o","--output_folder", type=str, help="Folder to store the report output")
    
    args = parser.parse_args()

    print("****************************************************")
    print(f"Running Python CYP2D6 CNV Report Script with arguments: {args}")

    # Create report directory
    report_dir = args.output_folder
    if not os.path.exists(report_dir):
        print(f"\nCreating report directory: {report_dir}")
        os.makedirs(report_dir)
    else:
        print(f"\nUsing existing report directory: {report_dir}")

    log_file = os.path.join(report_dir, f"{args.input_sample_name}_CNV_report_log.txt")
    print(f"\nCreating log file: {log_file}")

    # Setup logging
    setup_logging(log_file)

    logging.info("****************************************************")
    logging.info(f"Running Python generate_CNV_reports.py {args.input_sample_name} {args.input_cnv} {args.input_vcf} {args.input_bed} {args.input_z_threshold} {args.input_cn_threshold} {args.input_prop_threshold} {args.report_cyp2d6} {args.selected_genes} {args.output_folder}")
    logging.info(f"Start time: {datetime.now()}")

    # Read input CNV file
    logging.info(f"\nReading input CNV file: {args.input_cnv}")
    df_cnv = read_cnv_file(args.input_cnv)

    # Read input VCF file
    logging.info(f"\nReading input VCF file and calculaing vafs: {args.input_vcf}")
    df_vaf = get_vafs(args.input_vcf, args.input_sample_name, report_dir)

    # Extract gene limits
    df_gene_limits = gene_starts_and_ends(args.input_bed)

    if args.report_cyp2d6 == "true":
        # Create CYP2D6 region plot
        logging.info("\nCreating CYP2D6 region plot...")
        plot_region_with_layers(df_cnv, df_vaf, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, "chr22", 42123143, 42155250, args.input_bed, "CYP2D6-D8", report_dir)

    # Create plots for selected genes
    cnvgenes = list(args.selected_genes.split(","))
    for index, row in df_gene_limits.iterrows():
        thischr, thisstart, thisend, thisgene = row["chr"], row["start"], row["end"], row["gene"]
        
        #if any(w in thisgene for w in cnvgenes):
        if thisgene in cnvgenes:
            print(f"Plotting gene: {thisgene} ({thischr}:{thisstart}-{thisend})")

            # Call the plot_region function for each gene
            plot_region_with_layers(df_cnv, df_vaf, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, thischr, thisstart, thisend, args.input_bed, thisgene, report_dir)

    # Final log statements
    logging.info("Analysis complete.")
    logging.info(f"End time: {datetime.now()}")

if __name__ == "__main__":
    main()
