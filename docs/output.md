# PrecisionPGx: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Output directories

- [PrecisionPGx-outputdir](#precisionpgx-output)
  - [alignment](#alignment)
  - [CYP2D6](#cyp2d6)
  - [HLA](#hla)
  - [multiqc](#multiqc)
    - [multiqc data](#multiqc_data)
    - [multiqc plots](#multiqc_plots) 
  - [pharmcat](#pharmcat)
    - [complete report](#complete_report)
    - [selected genes](#selected_genes)
  - [pipeline info](#pipeline_info)
  - [qc bam](#qc_bam)
    - [qualimap](#qualimap)
  - [snv calls](#snv_calls)
  - [target depth](#target_depth)
  - [trimming](#trimming)

### `alignment`

Alignment output.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/alignment/`
  - `*_sorted_md.bam`: Alignment file in bam format.
  - `*.bai`: Index of the corresponding bam file.
  - `*_sorted_md.MarkDuplicates.metrics.txt`: Text file containing the deduplication metrics.
  </details>

### `CYP2D6`

CYP2D6 genotyping output.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/CYP2D6/`
  - `*.pharmcat.CYP2D6.match.html`: PharmCAT matcher report in html format.
  - `*.pharmcat.CYP2D6.match.json`: PharmCAT matcher report in json format.
  - `*.pharmcat.CYP2D6.phenotype.json`: PharmCAT phenotyper report in json format.
  </details>

### `HLA`

HLA genotyping output.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/HLA/`
  - `*.optitype.coverage_plot.pdf`: Optitype coverage report in pdf format.
  - `*.optitype.result.tsv`: Optitype HLA genotyping results in tsv format.
  </details>

### `multiqc`

MultiQC report. Underlying data and plots in subfolders.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/multiqc/`
  - `multiqc_report.html`: MultiQC report in html format.
</details>

#### `multiqc_data`

Files containing the data presented in the MultiQC report.

#### `multiqc_plots`

Plots included in the MultiQC report in subfolders for formats pdf, png and svg. 

### `pharmcat`

PharmCAT VCF processing output. Genotyping and reporting output in subfolders.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/pharmcat/`
  - `*.pharmcat.preprocessed.vcf.bgz`: VCF preprocessed by PharmCAT. All positions relevant for genotyping that were missing from input VCF have been set to reference. 
  - `*.pharmcat.preprocessed.pass.vcf.gz`: Filtered VCF preprocessed by PharmCAT. Positions where read depth was lower than config parameter `min_dp` have been removed.
  - `*.pharmcat.missing_pgx_var.vcf`: PharmCAT output VCF including all positions without a variant.
  - `*.tbi`: Index of corresponding VCF.
</details>

#### `complete_report`

PharmCAT genotyping and reporting output.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/pharmcat/complete_report`
  - `*.pharmcat.match.html`: PharmCAT matcher report in html format.
  - `*.pharmcat.match.json`: PharmCAT matcher report in json format.
  - `*.pharmcat.phenotype.json`: PharmCAT phenotyper report in json format.
  - `*.pharmcat.report.html`: PharmCAT final report in html format.
  - `*.pharmcat.report.json`: PharmCAT final report in json format. 
  - `*.pharmcat.report.tsv`: PharmCAT final report in tsv format.
</details>

#### `selected_genes`

PharmCAT genotyping and reporting output including only selected genes.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/pharmcat/selected_genes`
  - `*.pharmcat.selected_genes.match.html`: PharmCAT matcher report in html format.
  - `*.pharmcat.selected_genes.match.json`: PharmCAT matcher report in json format.
  - `*.pharmcat.selected_genes.phenotype.json`: PharmCAT phenotyper report in json format.
  - `*.pharmcat.selected_genes.report.html`: PharmCAT final report in html format.
  - `*.pharmcat.selected_genes.report.json`: PharmCAT final report in json format.
  - `*.pharmcat.selected_genes.report.tsv`: PharmCAT final report in tsv format.
</details>

### `pipeline_info`

Execution and pipeline reports in html and txt format.
Input parameters in json format and software versions in yml format.

### `qc_bam`

Output from qc analysis of bam files.
Qualimap output in sample subfolders.

### `snv_calls`

SNV calling output.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/snv_calls/`
  - `*.{variant_caller}.vcf.gz`: VCF from given variant caller.
  - `*.{variant_caller}.g.vcf.gz`: gVCF from given variant caller.
  - `*.{variant_caller}.norm.vcf.gz`: Normalised VCF from given variant caller.
  - `*.{variant_caller}.norm.target.vcf.gz`: Normalised and target region filtered VCF from given variant caller. Target regions are defined by config parameters `target_bed` and `bait_padding`.
  - `*.{variant_caller}.norm.target.filtered.vcf.gz`: Normalised, target region and quality filtered VCF from given variant caller. Quality filter parameters for each caller are given in `conf/modules/variant_calling.config`. This VCF is used as input to PharmCAT.
  - `*.tbi`: Index of corresponding VCF.
</details>

### `target_depth`
Target depth analysis output.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/target_depth/`
  - `*.target_depth.tsv`: Target depth output from samtools depth. Used for CNV calling and filtering of VCFs processed for PharmCAT genotyping.
  - `*.target.pass.bed`: Bed file describing the target regions where the read depth was at least the depth given with config parameter `min_dp`.
</details>

### `trimming`

Output from trimming of fastq files with fastp.

<details markdown="1">
<summary>Output files:</summary>

- `{outputdir}/trimming/`
  - `*.fastp.html`: Fastp report in html format.
  - `*.fastp.json`: Fastp report in json format.
  - `*.fastp.log`: Fastp log file.
  - `*.fastp.fastq.gz`: Trimmed fastq files. These files are used for alignment.
</details>
