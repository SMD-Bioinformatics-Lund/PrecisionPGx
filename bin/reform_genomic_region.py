#!/usr/bin/env python3

from ast import arg
import pandas as pd
import argparse
import sys
from distutils.util import strtobool


def reform(target_bed, output_f, detected_variants, padding, file_format, addchr):
    targets = pd.read_csv(
        target_bed,
        sep="\t",
        names=["CHROM", "START", "END", "ID", "GENE"],
        dtype={"START": int, "END": int},
    )
    if detected_variants != "":
        detected_rsid = pd.read_csv(detected_variants, sep="\t").ID
        targets = targets[~targets.ID.isin(detected_rsid)]

    bed = file_format == "bed"
    with open(output_f, "w+") as f:
        for i, row in targets.iterrows():
            chrom, start, end, id = row[0:4]
            if start > end:
                end, start = start, end
            start -= padding
            end += padding
            if bed:
                if addchr:
                    f.write(f"chr{chrom}\t{start}\t{end}\t{id}\n")
                else:
                    f.write(f"{chrom}\t{start}\t{end}\t{id}\n")
            else:
                if addchr:
                    f.write(f"chr{chrom}:{start}-{end}\n")
                else:
                    f.write(f"{chrom}:{start}-{end}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Rewrite bed to chr:start-end list. Removing wt targets or adding padding"
    )
    parser.add_argument("--target_bed", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--detected_variants", type=str, default="")
    parser.add_argument("--padding", type=int, default=0)
    parser.add_argument("--format", type=str, default="list")
    parser.add_argument(
        "--addchr",
        type=lambda x: bool(strtobool(x)),
        nargs="?",
        const=True,
        default=False,
        help="add chr to the chromosomes",
    )

    args = parser.parse_args(sys.argv[1:])
    target_bed = args.target_bed
    output_file = args.output_file
    detected_variants = args.detected_variants
    padding = args.padding
    file_format = args.format
    addchr = args.addchr
    reform(target_bed, output_file, detected_variants, padding, file_format, addchr)


if __name__ == "__main__":
    main()
