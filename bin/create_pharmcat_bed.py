#!/usr/bin/env python3
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="PharmCat Coverage Analysis")
    parser.add_argument("--pharmcat_vcf", type=str, help="Path to PharmCat positions VCF file")
    parser.add_argument("--output", type=str, help="Path to output PharmCat positions BED file")
    return parser.parse_args()


def create_bed_from_vcf(pharmcat_vcf, bed_path):
    with open(pharmcat_vcf, "r") as vcf_file, open(bed_path, "w") as bed_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos, variant_id, ref, alt, qual, filter, info = fields[:8]
            alt_alleles = alt.split(",")
            is_snv = len(ref) == 1 and all(len(a) == 1 for a in alt_alleles)
            
            annot = f"SNV_{variant_id}_" if is_snv else f"INDEL_{variant_id}_"
            try:
                info_gene = info.split('=')[1]
            except:
                info_gene = info.split('=')[0]

            bed_file.write(f"{chrom}\t{int(pos)-1}\t{pos}\t{annot}{info_gene}\t{chrom}:{pos}\n")
            
            if not is_snv:
                for allele in alt_alleles:
                    if len(allele) > 1:
                        bed_file.write(f"{chrom}\t{pos}\t{int(pos) + 1}\tafterINS_{variant_id}_{info_gene}\t{chrom}:{pos}\n")
                    elif len(ref) > 1:
                        bed_file.write(f"{chrom}\t{int(pos) + len(ref) - 1}\t{int(pos) + len(ref)}\tafterDEL_{variant_id}_{info_gene}\t{chrom}:{pos}\n")


def main():
    args = parse_arguments()
    create_bed_from_vcf(args.pharmcat_vcf, args.output)
    print(f"PharmCat BED file created at {args.output}")


if __name__ == "__main__":
    main()
