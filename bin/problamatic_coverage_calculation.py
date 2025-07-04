#!/usr/bin/env python3
"""
Usage:
    python problematic_regions_overlap_reports_and_plots_annotated_bed_pgx.py <sample_name> <depths_file> <panel_name> <target_bed> <depth_threshold> <overlap_bed> <outfolder>

This script reads a depths file (with columns: chr, position, depth),
a panel bed file (columns: chr, start, end, gene), and a problematic regions (overlap) bed file
(columns: chr, start, end, name). It generates an overlap report and plots the read depth
with overlapping problematic regions shaded. It also calculates per‐gene overlap metrics
and creates per‐gene plots.
"""

import argparse, os, sys, datetime, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps

#---------------------------
def overlap_report(depths_df, sample_name, bed_name, bed_df, depth_threshold, shade_df, output_folder, report_file, report_html, log_file):
    # Get panel chromosome first occurrence indices.
    panel_chrs = depths_df['chr'].unique()
    panel_chrs_row_no = []
    for ch in panel_chrs:
        idx = depths_df.index[depths_df['chr'] == ch][0]
        panel_chrs_row_no.append((ch, idx))
    log_file.write("Panel chromosomes and first row indices:\n")
    for ch, idx in panel_chrs_row_no:
        log_file.write(f"  {ch} -> {idx}\n")
    
    # Find overlapping regions between the depths file and the shade (problematic regions) bed.
    shade_regions = []  # will hold tuples: (start_index, end_index, shade_name)
    overlap_array = np.zeros(len(depths_df), dtype=int)
    this_start = None
    this_name = ""
    
    for j in range(len(depths_df)):
        cur_chr = depths_df.iloc[j]['chr']
        cur_pos = depths_df.iloc[j]['position']
        # check if new exon (non-contiguous positions or chromosome change)
        if j > 0:
            prev_chr = depths_df.iloc[j-1]['chr']
            prev_pos = depths_df.iloc[j-1]['position']
            if (cur_chr != prev_chr) or (cur_pos != prev_pos + 1):
                if this_start is not None:
                    this_end = j - 1
                    shade_regions.append((this_start, this_end, this_name))
                    this_start = None
                    this_name = ""
        # Determine if the current position is overlapped by any problematic region
        match = shade_df[(shade_df['chr'] == cur_chr) &
                            (shade_df['start'] < cur_pos) &
                            (shade_df['end'] >= cur_pos)]
        if not match.empty:
            overlap_array[j] = 1
            if this_start is None:
                this_start = j
                this_name = match.iloc[0]['name']
        else:
            if this_start is not None:
                this_end = j - 1
                shade_regions.append((this_start, this_end, this_name))
                this_start = None
                this_name = ""
    # If an overlap region was started at the end of the file, close it.
    if this_start is not None:
        shade_regions.append((this_start, len(depths_df) - 1, this_name))
    
    # Create complete panel read depth plot.
    y_series_depth = depths_df['depth'].astype(int).values
    # Create a boolean series (True when below threshold)
    threshold_series = depths_df['depth'] < depth_threshold

    plt.figure(figsize=(12, 8))
    # For color mapping: we use 1 for below threshold, 0 otherwise.
    marker_colors = threshold_series.astype(int)
    # Use marker size = (boolean+2)*5 (similar to Julia: threshold_series .+ 2)
    marker_sizes = (threshold_series.astype(int) + 2) * 5
    plt.scatter(range(len(y_series_depth)), y_series_depth, c=marker_colors, cmap=colormaps.get_cmap('RdBu'), s=marker_sizes)
    plt.xticks([])
    plt.xlabel("positions covered by panel")
    plt.ylabel("read depth")
    plt.axhline(y=depth_threshold, linestyle='dashed', color='grey')
    y_max = np.max(y_series_depth)
    # Shade overlapping regions.
    for (s, e, shade_name) in shade_regions:
        plt.axvspan(s, e, color='darkgrey', alpha=0.5)
    plt.annotate(f"* dark shaded intervals indicate overlap with problematic regions, red dots indicate positions below threshold {depth_threshold}",(1, 1.075*y_max), fontsize=8, color='black')
    # Draw vertical dotted lines at chromosome boundaries.
    for ch, idx in panel_chrs_row_no:
        plt.axvline(x=idx, linestyle='dotted', color='grey')
        # For annotation: shift y slightly (alternating a little)
        y_pos = y_max + (panel_chrs_row_no.index((ch, idx)) % 2)*0.05*y_max
        plt.annotate(ch, (idx, y_pos), fontsize=8, color='black')
    # Save the plot.
    outplot = os.path.join(output_folder, f"{sample_name}.{bed_name}.read_depth_and_overlap_complete_panel.png")
    plt.savefig(outplot, bbox_inches='tight')
    plt.close()
    log_file.write(f"Plot saved as {outplot}\n")
    
    # Calculate total overlapping positions.
    tot_pos = len(depths_df)
    tot_overlap = int(np.sum(overlap_array))
    percent_overlap = round(100 * tot_overlap / tot_pos, 2)
    
    # Write summary to text report.
    report_file.write("Overlap of complete panel:\n----------------------------------------------------------------\n")
    report_file.write(f"No. of bases covered in depth file: {tot_pos}\n")
    report_file.write(f"No. of bases that overlap with overlap bed file: {tot_overlap}\n")
    report_file.write(f"% of input depth positions that overlap: {percent_overlap}%\n")
    report_file.write("----------------------------------------------------------------\n\n")
    report_file.write(f"Total panel read depth plot: {outplot}\n")
    report_file.write("----------------------------------------------------------------\n")
    report_file.write("Genes with overlapping regions\n")
    
    # Write header to html report.
    report_html.write("<h2>Overlap of complete panel</h2>\n")
    report_html.write(f"No. of bases covered in depth file: {tot_pos}<br>\n")
    report_html.write(f"No. of bases that overlap with overlap bed file: {tot_overlap}<br>\n")
    report_html.write(f"% of input depth positions that overlap: {percent_overlap}%<br>\n")
    report_html.write(f"<img src=\"{sample_name}.{bed_name}.read_depth_and_overlap_complete_panel.png\">\n")
    report_html.write("<h2 id=\"gene_table\">Genes with overlapping regions</h2>\n")
    
    # Initialize an overlap region matrix (as list of lists)
    # Each row: [chr, pos_start-1, pos_end, gene, overlap_length, region_length, percent_overlap, problem_region]
    overlap_reg_mx = []
    for (s, e, shade_name) in shade_regions:
        pos_start = depths_df.iloc[s]['position']
        pos_end = depths_df.iloc[e]['position']
        overlap_length = pos_end - pos_start + 1
        # Find matching gene from panel bed file.
        match = bed_df[(bed_df['chr'] == depths_df.iloc[s]['chr']) & (bed_df['start'] < pos_start) & (bed_df['end'] >= pos_start)]
        if match.empty:
            this_gene = "NA"
            reg_length = 0
            percent_reg_overlap = 0.0
        else:
            this_gene = match.iloc[0]['gene']
            reg_length = match.iloc[0]['end'] - match.iloc[0]['start']
            percent_reg_overlap = round(100 * overlap_length / reg_length, 2) if reg_length > 0 else 0.0
        # Modify gene name: if splitting by "_" gives more than 2 parts, join first two parts.
        parts = this_gene.split("_")
        gene_mod = f"{parts[0]}_{parts[1]}" if len(parts) > 2 else parts[0]
        overlap_reg_mx.append([depths_df.iloc[s]['chr'], pos_start - 1, pos_end, gene_mod, overlap_length, reg_length, percent_reg_overlap, shade_name[0] if shade_name else ""])
    
    # Get list of unique genes with overlapping regions.
    overlap_genes = sorted(list({row[3] for row in overlap_reg_mx}))
    overlap_genes_df = pd.DataFrame({'gene': overlap_genes, 'totpos': 0, 'overlappos': 0, 'overlappc': 0.0})
    
    # For each gene, calculate total positions (from panel bed) and total overlapping length.
    for gene in overlap_genes:
        this_gene_bed_df = bed_df[bed_df['gene'].str.contains(gene)]
        this_gene_mx = [row for row in overlap_reg_mx if row[3] == gene]
        this_gene_overlap = sum([row[4] for row in this_gene_mx])
        this_gene_sum = 0
        for idx, row in this_gene_bed_df.iterrows():
            this_gene_sum += row['end'] - row['start']
        overlap_genes_df.loc[overlap_genes_df['gene'] == gene, 'totpos'] = this_gene_sum
        overlap_genes_df.loc[overlap_genes_df['gene'] == gene, 'overlappos'] = this_gene_overlap
        overlap_genes_df.loc[overlap_genes_df['gene'] == gene, 'overlappc'] = round(100 * this_gene_overlap / this_gene_sum, 2) if this_gene_sum > 0 else 0.0

    # Write gene table to reports.
    report_html.write("<table>\n<tr>\n<th>gene</th>\n<th>no. of pos. in panel</th>\n<th>overlapping positions</th>\n<th>% overlapping positions</th>\n</tr>\n")
    report_file.write("gene\tno. of pos. in panel\toverlapping positions\t% overlapping positions\t\n")
    for i, row in overlap_genes_df.iterrows():
        report_html.write(f"<tr>\n<td><a href=\"#{row['gene']}\">{row['gene']}</a></td>\n<td>{row['totpos']}</td>\n<td>{row['overlappos']}</td>\n<td>{row['overlappc']}</td>\n</tr>\n")
        report_file.write(f"{row['gene']}\t{row['totpos']}\t{row['overlappos']}\t{row['overlappc']}\n")
    report_html.write("</table><br>\n")
    
    # For each gene, print detailed report and plot.
    for gene in overlap_genes_df['gene']:
        this_gene_bed_df = bed_df[bed_df['gene'].str.contains(gene)]
        report_html.write(f"<h3 id=\"{gene}\">{gene}</h3>\n")
        if not this_gene_bed_df.empty:
            chr_val = this_gene_bed_df.iloc[0]['chr']
            start_val = this_gene_bed_df.iloc[0]['start']
            end_val = this_gene_bed_df.iloc[-1]['end']
            report_html.write(f"{chr_val}:{start_val}-{end_val}<br><br>\n")
        totpos_gene = overlap_genes_df.loc[overlap_genes_df['gene'] == gene, 'totpos'].values[0]
        overlappos_gene = overlap_genes_df.loc[overlap_genes_df['gene'] == gene, 'overlappos'].values[0]
        overlappc_gene = overlap_genes_df.loc[overlap_genes_df['gene'] == gene, 'overlappc'].values[0]
        report_html.write(f"Tot. no. of positions in panel: {totpos_gene}<br>\n")
        report_html.write(f"No. of overlapping positions: {overlappos_gene}<br>\n")
        report_html.write(f"% of gene positions in panel overlapping with problematic regions: {overlappc_gene}%<br><br>\n")
        report_html.write("<table>\n<tr>\n<th>chr:start-end</th>\n<th>gene</th>\n<th>overlap length</th>\n<th>exon length</th>\n<th>% of exon</th>\n<th>problem region</th>\n</tr>\n")
        report_file.write(f"\n{gene}\n")
        if not this_gene_bed_df.empty:
            report_file.write(f"{this_gene_bed_df.iloc[0]['chr']}:{this_gene_bed_df.iloc[0]['start']}-{this_gene_bed_df.iloc[-1]['end']}\n")
        report_file.write(f"Tot. no. of positions in panel: {totpos_gene}\n")
        report_file.write(f"No. of overlapping positions: {overlappos_gene}\n")
        report_file.write(f"% of gene positions in panel overlapping with problematic regions: {overlappc_gene}%\n")
        report_file.write("chr:start-end\tgene\toverlap length\texon length\t% of exon\tproblem region\n")
        this_gene_mx = [row for row in overlap_reg_mx if row[3] == gene]
        for row in this_gene_mx:
            report_file.write(f"{row[0]}:{row[1]}-{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\t{row[6]}\t{row[7]}\n")
            report_html.write(f"<tr>\n<td>{row[0]}:{row[1]}-{row[2]}</td>\n<td>{row[3]}</td>\n<td>{row[4]}</td>\n<td>{row[5]}</td>\n<td>{row[6]}%</td>\n<td>{row[7]}</td>\n</tr>\n")
        report_html.write("</table><br>\n")
        # Create overlap plot for gene.
        plot_gene_overlap(depths_df, sample_name, gene, depth_threshold, this_gene_bed_df, this_gene_mx, output_folder, log_file)
        gene_short = gene.split("/")[0]
        report_file.write(f"Gene read depth plot: {output_folder}/{sample_name}.{gene_short}.overlap_problem_regions.png\n")
        report_html.write(f"<img src=\"{sample_name}.{gene_short}.overlap_problem_regions.png\">\n")
        report_html.write("<br>\n<a href=\"#gene_table\">back to genes table</a>\n<br>\n")
        
    # End of overlap_report

#---------------------------
def plot_gene_overlap(depths_df, sample_name, gene_name, depth_threshold, gene_bed_df, gene_overlap_mx, output_folder, log_file):
    # Filter depths file for the gene region.
    chr_val = gene_bed_df.iloc[0]['chr']
    start_val = gene_bed_df.iloc[0]['start']
    end_val = gene_bed_df.iloc[-1]['end']
    this_gene_depth_df = depths_df[(depths_df['chr'] == chr_val) & (depths_df['position'] > start_val) & (depths_df['position'] <= end_val)]
    y_series_depth = this_gene_depth_df['depth'].astype(int).values
    threshold_series = this_gene_depth_df['depth'] < depth_threshold

    plt.figure(figsize=(12, 8))
    marker_colors = threshold_series.astype(int)
    marker_sizes = (threshold_series.astype(int) + 2) * 5
    plt.scatter(range(len(y_series_depth)), y_series_depth, c=marker_colors, cmap=colormaps.get_cmap('RdBu'), s=marker_sizes)
    plt.xticks([])
    plt.xlabel("positions covered by panel")
    plt.ylabel("read depth")
    plt.axhline(y=depth_threshold, linestyle='dashed', color='grey')
    y_max = np.max(y_series_depth)
    # Map overlap regions for this gene into indices relative to this_gene_depth_df.
    shade_regions = []
    for ov in gene_overlap_mx:
        # ov is a list: [chr, pos_start, pos_end, gene, ...]
        start_idx = None
        end_idx = None
        for j in range(len(this_gene_depth_df)):
            pos = this_gene_depth_df.iloc[j]['position']
            if start_idx is None and pos == (ov[1] + 1):
                start_idx = j
            if pos == ov[2]:
                end_idx = j
                break
        if start_idx is not None and end_idx is not None:
            shade_regions.append((start_idx, end_idx))
    for (s, e) in shade_regions:
        plt.axvspan(s, e, color='darkgrey', alpha=0.5)
    plt.annotate(f"* dark shaded intervals indicate overlap with problematic regions, red dots indicate positions below threshold {depth_threshold}",(1, 1.075*y_max), fontsize=8, color='black')
    # Add vertical dotted lines for exon boundaries.
    bed_regions = []
    for i, row in gene_bed_df.iterrows():
        for j in range(len(this_gene_depth_df)):
            if this_gene_depth_df.iloc[j]['position'] == row['start'] + 1:
                bed_regions.append(j)
                break
    for br in bed_regions:
        plt.axvline(x=br, linestyle='dotted', color='grey')
    gene_short = gene_name.split("/")[0]
    outplot = os.path.join(output_folder, f"{sample_name}.{gene_short}.overlap_problem_regions.png")
    log_file.write(f"Plot saved as {outplot}\n")
    plt.savefig(outplot, bbox_inches='tight')
    plt.close()


def read_depth_file(depth_file: str) -> pd.DataFrame:
    """Read depth file and return as a DataFrame.
    
    If the first column does not begin with 'chr' (e.g. it's numeric),
    the "Chr" column is converted to int along with "Position" and "Depth".
    Otherwise, "Chr" remains as a string.
    """
    df = pd.read_csv(
        depth_file,
        sep="\t",
        header=None,
        names=["chr", "position", "depth"],
        dtype={"chr": str},
    )
    return df


def is_chr_format(df: pd.DataFrame) -> bool:
    """Check if the 'Chr' column in the DataFrame contains 'chr' prefixes.
    
    Returns:
        True if at least one entry in 'Chr' starts with 'chr'.
        False if all entries are purely numeric.
    """
    if df.empty:
        raise ValueError("DataFrame is empty.")

    # Convert first value to string and check if any starts with 'chr'
    return str(df["chr"].iloc[0]).lower().startswith("chr")


def read_bed_file(bed_file: str, is_chr: bool) -> pd.DataFrame:
    """Read BED file and return as a DataFrame.
    
    If `is_chr` is True and "Chr" values do not start with "chr", it adds "chr".
    If `is_chr` is False and "Chr" values start with "chr", it removes "chr".
    
    Args:
        bed_file (str): Path to the BED file.
        is_chr (bool): Indicates whether the chromosome column should have "chr" prefixes.
    
    Returns:
        pd.DataFrame: BED file as a DataFrame with adjusted "Chr" column.
    """
    df = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "gene"],
        dtype={"chr": str, "start": int, "end": int, "gene": str},
    )

    # Check if the first entry in "Chr" starts with "chr"
    first_val = str(df["chr"].iloc[0]).lower()
    has_chr_prefix = first_val.startswith("chr")

    if is_chr and not has_chr_prefix:
        # If `is_chr` is True but the file lacks "chr", add it.
        df["chr"] = "chr" + df["chr"].astype(str)
    elif not is_chr and has_chr_prefix:
        # If `is_chr` is False but the file has "chr", remove it.
        df["chr"] = df["chr"].str.replace("^chr", "", regex=True)

    return df

def read_overlap_bed(bed_file: str, is_chr: bool) -> pd.DataFrame:
    """Read BED file and return as a DataFrame.
    
    If `is_chr` is True and "Chr" values do not start with "chr", it adds "chr".
    If `is_chr` is False and "Chr" values start with "chr", it removes "chr".
    
    Args:
        bed_file (str): Path to the BED file.
        is_chr (bool): Indicates whether the chromosome column should have "chr" prefixes.
    
    Returns:
        pd.DataFrame: BED file as a DataFrame with adjusted "Chr" column.
    """
    df = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "name"],
        dtype={"chr": str, "start": int, "end": int, "name": str},
    )

    # Check if the first entry in "Chr" starts with "chr"
    first_val = str(df["chr"].iloc[0]).lower()
    has_chr_prefix = first_val.startswith("chr")

    if is_chr and not has_chr_prefix:
        # If `is_chr` is True but the file lacks "chr", add it.
        df["chr"] = "chr" + df["chr"].astype(str)
    elif not is_chr and has_chr_prefix:
        # If `is_chr` is False but the file has "chr", remove it.
        df["chr"] = df["chr"].str.replace("^chr", "", regex=True)

    return df



#---------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate overlap report and plots (Python version) for problematic regions in a panel.")
    parser.add_argument("--sample_name", help="Sample name")
    parser.add_argument("--depths_file", help="Depths file (tab-delimited, columns: chr, position, depth)")
    parser.add_argument("--panel_name", help="Panel name")
    parser.add_argument("--target_bed", help="Target bed file (columns: chr, start, end, gene)")
    parser.add_argument("--depth_threshold", type=float, help="Depth threshold")
    parser.add_argument("--overlap_bed", help="Overlap bed file (columns: chr, start, end, name)")
    parser.add_argument("--outfolder", help="Output folder")
    args = parser.parse_args()

    if not args.outfolder:
        args.outfolder = os.getcwd()

    os.makedirs(args.outfolder, exist_ok=True)

    # Open report and log files.
    report_filename = os.path.join(args.outfolder, f"{args.sample_name}.{args.panel_name}.overlap_report.txt")
    report_html_filename = os.path.join(args.outfolder, f"{args.sample_name}.{args.panel_name}.overlap_report.html")
    log_filename = os.path.join(args.outfolder, f"{args.sample_name}.{args.panel_name}.overlap_log.txt")
    report_file = open(report_filename, "w")
    report_html = open(report_html_filename, "w")
    log_file = open(log_filename, "w")
    now_str = datetime.datetime.now().isoformat()
    report_file.write(now_str + "\n")
    report_file.write(f"Overlap report for sample {args.sample_name} and panel {args.panel_name}\n")
    report_file.write(f"Input depths file: {args.depths_file}\n")
    report_file.write(f"Input panel bed file: {args.target_bed}\n")
    report_file.write(f"Input problematic regions panel bed file: {args.overlap_bed}\n")
    report_file.write(f"Input depth threshold: {args.depth_threshold}\n\n")
    report_html.write(f"<!doctype html>\n<html>\n<head>\n<title>Overlap report {args.sample_name} panel {args.panel_name}</title>\n</head>\n<body>\n")
    report_html.write(now_str + "<br>\n")
    report_html.write(f"<h1>Overlap report: {args.sample_name}</h1>\n")
    report_html.write(f"sample: {args.sample_name}<br>panel: {args.panel_name}<br>\n")
    report_html.write("<h2>Input files</h2>\n")
    report_html.write(f"Input depths file: {args.depths_file}<br>\n")
    report_html.write(f"Input panel bed file: {args.target_bed}<br>\n")
    report_html.write(f"Input problematic regions panel bed file: {args.overlap_bed}<br>\n")
    report_html.write(f"Input depth threshold: {args.depth_threshold}<br><br>\n")
    log_file.write(now_str + "\n")
    log_file.write(f"Overlap analysis {args.sample_name}\n")
    log_file.write(f"python {sys.argv[0]} {args.sample_name} {args.depths_file} {args.panel_name} {args.target_bed} {args.depth_threshold} {args.overlap_bed} {args.outfolder}\n")

    # Read input files.
    depths_df = read_depth_file(args.depths_file)
    is_chr = is_chr_format(depths_df)
    bed_df = read_bed_file(args.target_bed, is_chr)
    shade_df = read_overlap_bed(args.overlap_bed, is_chr)

    overlap_report(depths_df, args.sample_name, args.panel_name, bed_df, args.depth_threshold, shade_df, args.outfolder, report_file, report_html, log_file)

    report_html.write("\n</body>\n</html>")
    report_file.close()
    report_html.close()
    log_file.close()

if __name__ == "__main__":
    main()
