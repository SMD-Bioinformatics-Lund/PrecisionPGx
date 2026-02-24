#!/usr/bin/env bash


input="$1"  # Input VCF file path
output="$2" # Output VCF file path


# Fix Chromosome names in the input VCF using Perl, adding "chr" prefix to chromosome names that do not already have it and skipping header lines
perl -pe '/^((?!^chr).)*$/ && s/^([^#])/chr$1/gsi' $input > merged_output.chrfixed.vcf

