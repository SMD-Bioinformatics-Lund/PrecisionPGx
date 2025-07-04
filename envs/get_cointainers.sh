#!/bin/bash

# Run me in envs folder to get sinuglarities needed
# Please adjust paths of simg location in appropriate config_file


sudo singularity build target_variants_python.sif recepies/get_target_variants
sudo singularity build jinja_report.sif recepies/jinja_report
sudo singularity build samtools.sif recepies/samtools
sudo singularity build gatk3.sif docker://broadinstitute/gatk3:3.8-1
sudo singularity build gatk4.sif docker://broadinstitute/gatk
sudo singularity build bcftools_1.20.sif recepies/bcftools
sudo singularity build pharmcat.sif docker://pgkb/pharmcat



