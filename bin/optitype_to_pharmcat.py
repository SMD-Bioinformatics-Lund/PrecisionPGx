#!/usr/bin/env python

import pandas as pd
import sys
import os

# Prompt for the Sample Name, which will be used as the file name and input file path to the TSV file
sample_name = sys.argv[1]
input_path = sys.argv[2]

# Load TSV file as a dataframe
df = pd.read_csv(input_path, sep="\t")

# Reformat Data
a_values = df.iloc[0, 1:3].str.replace("[AB](?=\\*)", "", regex=True).str.cat(sep="/")
b_values = df.iloc[0, 3:5].str.replace("[AB](?=\\*)", "", regex=True).str.cat(sep="/")

# Create new dataframe with reformatted data
new_df = pd.DataFrame([["HLA-A", a_values], ["HLA-B", b_values]])

# Export the data frame into a Tab-delimited file
new_df.to_csv(f"{sample_name}_HLA.txt", sep="\t", index=False, header=False)
