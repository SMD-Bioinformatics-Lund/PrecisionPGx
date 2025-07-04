#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict
import logging

# Supported variant callers
SUPPORTED_CALLERS = {
    'freebayes', 'mutect2', 'tnscope', 'vardict', 'pindel',
    'bcftools', 'gatk-haplotyper', 'sentieon-haplotyper'
}

# ---------------------------
# Parse Command-Line Arguments
# ---------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate multiple VCFs into one.")
    parser.add_argument('--vcfs', required=True, help='Comma-separated list of VCF files to aggregate')
    parser.add_argument('--tumor-id', help='Tumor sample ID')
    parser.add_argument('--normal-id', help='Normal sample ID')
    parser.add_argument('--fluffify-pindel', action='store_true', help='Modify Pindel REF/ALT fields')
    parser.add_argument('--sample-order', help='Comma-separated list of sample order')
    return parser.parse_args()

# ---------------------------
# Read VCF File
# ---------------------------
def read_vcf(file):
    metadata = []
    variants = []
    sample_names = []

    opener = gzip.open if file.endswith('.gz') else open
    with opener(file, 'rt') as f:
        for line in f:
            if line.startswith("##"):
                metadata.append(line.strip())
            elif line.startswith("#CHROM"):
                header = line.strip().split("\t")
                sample_names = header[9:]  # Extract sample names from header
            else:
                variants.append(line.strip().split("\t"))

    return metadata, sample_names, variants

# ---------------------------
# Determine Variant Caller
# ---------------------------
def which_variantcaller(metadata):
    source = None
    for line in metadata:
        if line.startswith("##GATKCommandLine="):
            source = "gatk-haplotyper"
            break
        elif line.startswith("##SentieonCommandLine.Haplotyper"):
            source = "sentieon-haplotyper"
            break
        elif line.startswith("##bcftoolsCommand=mpileup"):
            source = "bcftools"
            break
        elif line.startswith("##source="):
            source = line.split("##source=")[1].split(" ")[0].split("_")[0].lower()
            break
    for caller in SUPPORTED_CALLERS:
        if source and caller in source:
            return caller
    return "unknown"

# ---------------------------
# Summarize Filters
# ---------------------------
def summarize_filters(filters):
    non_pass = [f for f in filters if f not in {"PASS", "."}]
    return "PASS" if not non_pass else ";".join(non_pass)

# ---------------------------
# Add INFO Field
# ---------------------------
def add_info(var, key, val):
    """ Mimics Perl's add_info: adds a key-value pair to the INFO field """
    if 'INFO_order' not in var:
        var['INFO_order'] = []
    var['INFO_order'].append(key)
    if 'info' not in var:
        var['info'] = {}
    var['info'][key] = val

# ---------------------------
# Add Genotype Data
# ---------------------------
def add_gt(var, sample, key, val):
    """ Mimics Perl's add_gt: adds FORMAT key and assigns genotype values per sample """
    if 'FORMAT' not in var:
        var['FORMAT'] = []
    if key not in var['FORMAT']:
        var['FORMAT'].append(key)
    for format_value in var.get('FORMAT_VAL', []):
        if format_value.get('_sample_id') == sample:
            format_value[key] = val

# ---------------------------
# Fix Genotype (FORMAT) Fields
# ---------------------------
def fix_gt(var, info_dict, caller):
    var['FORMAT'] = []  # Reset FORMAT field

    for format_val in var['FORMAT_VAL']:
        gt = format_val.get("GT", "./.")
        dp = format_val.get("DP", ".")
        vd = format_val.get("AD", ".,.").split(",")[1]
        vaf = format_val.get("AF", ".")

        if caller in {"mutect2", "tnscope", "vardict", "pindel"}:
            ref_dp, alt_dp = 0, 0
            if 'AD' in format_val and format_val['AD']:
                parts = format_val['AD'].split(",")
                ref_dp = int(parts[0]) if parts[0] else 0
                alt_dp = int(parts[1]) if len(parts) > 1 and parts[1] else 0
            vaf = round(float(alt_dp / (alt_dp + ref_dp) if (alt_dp + ref_dp) > 0 else 0),4)
            dp = ref_dp + alt_dp
            vd = alt_dp
        elif caller in {"freebayes", "gatk-haplotyper", "sentieon-haplotyper"}:
            try:
                vaf = round(float(int(vd) / int(dp)),4)
            except:
                vaf = "."
        elif caller == "bcftools":
            dp4 = info_dict.get("DP4", "0,0,0,0").split(",")
            ref_dp = int(dp4[0]) + int(dp4[1])
            alt_dp = int(dp4[2]) + int(dp4[3])
            dp = info_dict.get("DP", "0")
            fAD = format_val.get("AD", "0,0,0,0").split(",")

            if alt_dp in fAD:
                vd = alt_dp
            else:
                vd = 0

            try:
                vaf = round(float(vd / (dp)),4)
            except:
                vaf = "."

        add_gt(var, format_val['_sample_id'], "DP", dp)
        add_gt(var, format_val['_sample_id'], "VD", str(vd))
        add_gt(var, format_val['_sample_id'], "VAF", str(vaf))
        add_gt(var, format_val['_sample_id'], "GT", str(gt))

# ---------------------------
# Aggregate VCF Files
# ---------------------------
def aggregate_vcfs(vcf_files):
    aggregated = {}
    filters = defaultdict(set)
    all_filters = set()
    headers = []
    global_samples = set()

    for vcf_file in vcf_files:
        metadata, samples, variants = read_vcf(vcf_file)
        headers.append(metadata)
        global_samples.update(samples)
        caller = which_variantcaller(metadata)
        logging.info(f"Processing {vcf_file} (caller: {caller})")

        for var in variants:
            chrom, pos, var_id, ref, alt, qual, filt, info, fmt, *sample_data = var
            info_dict = {item.split("=")[0]: item.split("=")[1] for item in info.split(";") if "=" in item}
            dp_value = info_dict.get("DP", ".")
            key = f"{chrom}_{pos}_{ref}_{alt}"

            var_dict = {
                'CHROM': chrom,
                'POS': pos,
                'ID': var_id,
                'REF': ref,
                'ALT': alt,
                'QUAL': qual,
                'filter': filt,
                'info': {},
                'INFO_order': []
            }

            fmt_keys = fmt.split(":")
            gt_array = [{fmt_key: sample_data[i].split(":")[j] for j, fmt_key in enumerate(fmt_keys)} | {'_sample_id': sample} for i, sample in enumerate(samples)]
            var_dict['FORMAT_VAL'] = gt_array

            if key in aggregated:
                aggregated[key]['info']['variant_callers'] += f",{caller}"
                aggregated[key]['info']['VC'] = int(aggregated[key]['info']['VC']) + 1
            else:
                add_info(var_dict, "variant_callers", caller)
                add_info(var_dict, "DP", dp_value)
                add_info(var_dict, "VC", 1)
                fix_gt(var_dict, info_dict, caller)
                aggregated[key] = var_dict

    for key in aggregated:
        aggregated[key]['filter'] = summarize_filters(filters[key])

    return aggregated, headers, sorted(all_filters), list(global_samples)

# ---------------------------
# Print VCF Header
# ---------------------------
def print_header(filters, vcf_files, global_samples, command_line=None):
    print("##fileformat=VCFv4.2")
    print(f"##origin={vcf_files[0]}")
    print(f"##sources={','.join(vcf_files)}")
    print('##INFO=<ID=variant_callers,Number=.,Type=String,Description="List of variant callers">')
    print('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
    print('##INFO=<ID=VC,Number=1,Type=Integer,Description="Number of variant callers that called this variant">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
    print('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="ALT allele observation fraction">')
    print('##FORMAT=<ID=VD,Number=1,Type=Integer,Description="ALT allele observation count">')
    print(f'##CommandLine="{command_line}">')
    for f in filters:
        print(f'##FILTER=<ID={f},Description="{f}">')
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sorted(global_samples)))

# ---------------------------
# Print VCF Variant Record
# ---------------------------
def vcfstr(v):
    fields = [
        v.get("CHROM", "."),
        v.get("POS", "."),
        v.get("ID", "."),
        v.get("REF", "."),
        v.get("ALT", "."),
        v.get("QUAL", "."),
        v.get("filter", "."),
        ";".join([f"{key}={v['info'][key]}" for key in v.get("INFO_order", [])]),
        ":".join(v.get("FORMAT", []))
    ]

    sample_fields = [":".join([gt.get(key, ".") for key in v.get("FORMAT", [])]) for gt in v.get("FORMAT_VAL", [])]
    fields.append("\t".join(sample_fields))

    print("\t".join(fields))

# ---------------------------
# Main Routine
# ---------------------------
def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO)
    command_line = " ".join(["aggregate_vcf.py"] + ["--" + k + " " + str(v) for k, v in vars(args).items() if v])
    logging.info(command_line)
    vcf_files = args.vcfs.split(',')
    aggregated_vcfs, headers, all_filters, global_samples = aggregate_vcfs(vcf_files)
    print_header(all_filters, vcf_files, global_samples, command_line)

    for key in aggregated_vcfs:
        vcfstr(aggregated_vcfs[key])

if __name__ == "__main__":
    main()
