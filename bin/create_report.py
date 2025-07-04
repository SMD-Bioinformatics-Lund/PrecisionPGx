#!/usr/bin/env python3

import os
import sys
import argparse
from report_class import Report
import argparse


# DEFAULTS
TEMPLATE_FOLDER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "resources", "templates"
)

# Create a Pharmacogenomics (PGx) Report
parser = argparse.ArgumentParser(description="Generate a Pharmacogenomics (PGx) report")

# Required Arguments
parser.add_argument("--group", type=str, required=True, help="Specify the sample group")
parser.add_argument(
    "--read_depth",
    type=int,
    default=100,
    help="Read depth threshold for target regions (default: 100)",
)
parser.add_argument(
    "--report_template",
    type=str,
    default=os.path.join(
        TEMPLATE_FOLDER,
        "report.html",
    ),
    help=f"Path to the report template file (default: {os.path.join(TEMPLATE_FOLDER,'report.html')})",
)
parser.add_argument(
    "--detected_variants",
    type=str,
    required=True,
    help="Path to detected variants file",
)
parser.add_argument(
    "--missing_annotated_depth",
    type=str,
    required=True,
    help="Path to missing annotated depth file",
)
parser.add_argument(
    "--haplotype_definitions",
    type=str,
    required=True,
    help="Path to haplotype definitions file",
)
parser.add_argument(
    "--possible_diplotypes",
    type=str,
    required=True,
    help="Path to possible diplotypes file",
)
parser.add_argument(
    "--possible_interactions",
    type=str,
    required=True,
    help="Path to possible interactions file",
)
parser.add_argument(
    "--target_bed", type=str, required=True, help="Path to the target BED file"
)
parser.add_argument(
    "--padded_baits_depth",
    type=str,
    required=True,
    help="Path to padded baits depth file",
)
parser.add_argument(
    "--target_rsids", type=str, required=True, help="Path to target rsIDs file"
)
parser.add_argument(
    "--annotated_vcf", type=str, required=True, help="Path to annotated VCF file"
)
parser.add_argument(
    "--dbSNP_version", type=str, required=True, help="Specify the dbSNP version"
)
parser.add_argument(
    "--genome_version", type=str, required=True, help="Specify the Genome version"
)
parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="Specify the output file path for the PGx report",
)
parser.add_argument(
    "--logo",
    type=str,
    default=os.path.join(
        TEMPLATE_FOLDER,
        "rs_logo_rgb.png",
    ),
    help=f"Specify the logo file path for the PGx report (default: {os.path.join(TEMPLATE_FOLDER,'rs_logo_rgb.png')})",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # Example usage:
    report_instance = Report(
        group=args.group,
        read_depth=args.read_depth,
        detected_variants=args.detected_variants,
        missing_annotated_depth=args.missing_annotated_depth,
        haplotype_definitions=args.haplotype_definitions,
        possible_diplotypes=args.possible_diplotypes,
        possible_interactions=args.possible_interactions,
        target_bed=args.target_bed,
        padded_baits_depth=args.padded_baits_depth,
        target_rsids=args.target_rsids,
        annotated_vcf=args.annotated_vcf,
        dbSNP_version=args.dbSNP_version,
        genome_version=args.genome_version,
        output=args.output,
        report_template=args.report_template,
        logo=args.logo,
    )

    report_instance.create_report()
