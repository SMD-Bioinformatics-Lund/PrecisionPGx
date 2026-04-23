# PrecisionPGx: Usage

## Table of contents
- [Introduction](#introduction)
- [Requirements](#requirements)
- [Reference files](#reference-files)
  - [Reference genome files](#reference-genome-files)
    - [Reference genome](#reference-genome)
    - [Reference genome indices](#reference-genome-indices)
    - [Target definition file](#target-definition-file)
- [Configuration settings](#configuration-settings)
    - [Main config parameters](#main-config-parameters)
    - [Subworkflow config settings](#subworkflow-config-settings)
    - [Resource use config settings](#resource-use-config-settings)
- [Run the pipeline with test data](#run-the-pipeline-with-test-data)
- [Run the pipeline with input data](#run-the-pipeline-with-input-data)
    - [Input sample sheet](#input-sample-sheet)
    - [Run the pipeline with input parameters in parameters.yaml](#run-the-pipeline-with-input-parameters-in-parametersyaml)
    - [Run the pipeline with user profile and user_profile.config](#run-the-pipeline-with-user-profile-and-user_profileconfig)
    - [Nextflow logs](#nextflow-logs)
- [Pipeline overview](#pipeline-overview)
    - [Prepare references](#prepare-references)
    - [Decompress spring compressed fastq](#decompress-spring-compressed-fastq)
    - [Alignment](#alignment)
        - [Aligner](#aligner)
        - [Trim FASTQ](#trim-fastq)
        - [Duplicates](#duplicates)
        - [Output folder](#output-folder)
    - [QC BAM](#qc-bam)
        - [QC statistics](#qc-statistics)
        - [Contamination](#contamination)
        - [Output folder](#output-folder-1)
    - [Variant calling](#variant-calling)
        - [Variant caller](#variant-caller)
        - [DeepVariant](#deepvariant)
        - [Sentieon](#sentieon)
        - [GATK haplotype caller](#gatk-haplotype-caller)
        - [Output folder](#output-folder-2)
    - [Variant filtration](#variant-filtration)
        - [Normalise VCF](#normalise-vcf)
        - [Merge sample VCF with PharmCAT definition VCF](#merge-sample-vcf-with-pharmcat-definition-vcf)
        - [Create target bed filtered VCF](#create-target-bed-filtered-vcf)
        - [Create additionally filtered VCF](#create-additionally-filtered-vcf)
        - [Output folder](#output-folder-3)
    - [Target depth](#target-depth)
        - [Output folder](#output-folder-4)
    - [CNV calling](#cnv_calling)
        - [Output folder](#output-folder-5)
    - [HLA calling](#hla_calling)
        - [Output folder](#output-folder-6)
    - [PharmCAT VCF processing](#pharmcat-vcf-processing)
        - [Create PharmCAT formatted VCF](#create-pharmcat-formatted-vcf)
            - [Output folder](#output-folder-7)
        - [Filter PharmCAT VCF by target depth](#filter-pharmcat-vcf-by-target-depth)
            - [Output folder](#output-folder-8)
    - [PharmCAT genotyping and reporting](#pharmcat-genotyping-and-reporting)
        - [PharmCAT matcher](#pharmcat-matcher)
            - [Output folder](#output-folder-9)
        - [PharmCAT phenotyper](#pharmcat-phenotyper)
            - [Output folder](#output-folder-10)
        - [PharmCAT reporter](#pharmcat-reporter)
            - [Output folder](#output-folder-11)
    - [CYP2D6 calling](#cyp2d6_calling)
        - [Output folder](#output-folder-12)
    - [Collect software versions and multiqc](#collect-software-versions-and-multiqc)
        - [MultiQC](#multiqc)
        - [Output folder](#output-folder-13)

## Introduction

PrecisionPGx calls variants from targeted (panel) or whole genome sequencing data and does pharmacogenetic genotyping and annotation, primarily using the tool [PharmCAT](https://pharmcat.clinpgx.org/).

The pipeline is built using [Nextflow](https://www.nextflow.io/).

This document includes an overview of the pipeline and instructions for how to install and run the pipeline.

## Requirements

Nextflow >=25.04.0

Container technology: Docker, Singularity, Apptainer, Podman, Shifter or Charliecloud. *Only singularity tested*


## Reference files

Reference files are either provided or downloaded from iGenomes. *iGenomes download has not been tested yet*

Reference files that must be provided are described below with their corresponding config parameters.

### Reference genome files

*Note: PharmCAT only supports GRCh38*

The reference genome must be GRCh38. 

Inclusion of alternate contigs has not been tested but may cause problems for some genes, e.g. *CYP2D6*, where alternate contigs include some of the haplotype-defining variations.   

#### Reference genome
Set the reference genome with the config parameter `fasta`.

`fasta`: hg38_no_alt.fa

Paths to indices for the reference genome must also be provided with config parameters `fai`and `bwa`/`bwamem2`. *Automatic indexing if these are missing should be added*

Reference indexed with samtools:

`fai`: hg38_no_alt.fa.fai

Example command to index with samtools:

```
samtools faidx hg38_no_alt.fa
```

Reference indexed with bwa or bwa-mem2:

`bwa`: Path to directory for pre-built bwa index

`bwamem2`: Path to directory for pre-built bwa index

Example command to index with bwa-mem2: 

```
bwa-mem2 index hg38_no_alt.fa
```

Bwa-mem2 index includes the following files:

`hg38_no_alt.fa`

`hg38_no_alt.fa.0123`

`hg38_no_alt.fa.amb`

`hg38_no_alt.fa.ann`

`hg38_no_alt.fa.bwt.2bit.64`

`hg38_no_alt.fa.pac`

#### PharmCAT reference genome
A separate reference genome is used for PharmCAT. This is the reference genome that is automatically downloaded when running PharmCAT without a provided reference genome. 

PharmCAT provides these files on [zenodo](https://zenodo.org/records/7288118) with a description of how they were created.

Set PharmCAT reference genome and indices with config parameters `pharmcat_reference_fasta`, `pharmcat_reference_fasta_index` and `pharmcat_reference_fasta_fai`.

`pharmcat_reference_fasta`: /path/to/pharmcat_reference_genome/HG38_p13_pharmcat/reference.fna.bgz

`pharmcat_reference_fasta_index`: /path/to/pharmcat_reference_genome/HG38_p13_pharmcat/reference.fna.bgz.gzi

`pharmcat_reference_fasta_fai`: /path/to/pharmcat_reference_genome/HG38_p13_pharmcat/reference.fna.bgz.fai

#### Target definition file

The target definition file is used to define the regions where variants are called and where depth and copy numbers are calculated.

The target definition file is provided in the repo. *The target bed needs to be updated with the new design including new ID-SNP and chrX sites*

Set the target bed with the config parameter `target_bed`.

`target_bed`: assets/panpgx.TE-97619858.grch38.refseq_mane.whole_genes_snps.annot.v2.bed

The input target bed can be padded with a given number of bases to include broader regions.

Set the pad length with the config parameter `bait_padding`.

`bait_padding`: 50

## Configuration settings

### Main config parameters

Config parameters are set either in a `*.config` file or input `parameters.yaml` file, or a combination of both.

The main config file, `nextflow.config`, holds default or null values for input parameters. 

*To-do: Remove uny unused parameters from config. Add explanations where needed.*
```
// Global default params, used in configs
params {

    // TODO nf-core: Specify your pipeline's command line flags
    // Input options
    input                           = null

    // Reference genome information; iGenomes is effectively disabled but retained for linting
    genome                          = 'GRCh38'
    igenomes_base                   = 's3://ngi-igenomes/igenomes/'
    igenomes_ignore                 = false
    local_genomes                   = null
    save_reference                  = false

    // Main params
    analysis_type                   = "panel" // panel or wgs
    bait_padding                    = 100
    extract_alignments              = false
    restrict_to_contigs             = null
    save_mapped_as_cram             = false
    skip_tools                      = null
    skip_subworkflows               = null
    platform                        = 'illumina'

    // References files
    fasta                           = null
    fai                             = null
    sequence_dictionary             = null
    bwa                             = null
    bwamem2                         = null
    bwameme                         = null
    target_bed                      = null
    intervals_wgs                   = null

    // Alignment
    aligner                         = "bwamem2" // bwa, bwamem2, bwameme, sentieon
    bwa_as_fallback                 = false
    mbuffer_mem                     = 3072
    samtools_sort_threads           = 4
    min_trimmed_length              = 40
    rmdup                           = false

    // Haplotype/variant calling
    variant_caller                  = 'sentieon' // gatk4, sentieon, deepvariant
    min_dp                          = 50 // For Panel, change for WGS
    call_interval                   = null
    ml_model                        = null
    emit_mode                       = "confident" // variant, confident, all, gvcf: applicable only for sentieon
    filter_vcf                      = false

    // Verifybamid
    verifybamid_svd_bed             = null
    verifybamid_svd_mu              = null
    verifybamid_svd_ud              = null

    // CNV Calling
    // TODO:

    // HLA Calling
    // TODO:

    // Pharmcat
    //pharmcat                        = null
    //pharmcat_version                = null
    pharmcat_resource_dir           = null
    pharmcat_positions_vcf          = null
    pharmcat_positions              = null
    pharmcat_positions_index        = null
    pharmcat_uniallelic_pos         = null
    pharmcat_uniallelic_pos_index   = null
    pharmcat_reference_fasta        = null
    pharmcat_reference_fasta_index  = null
    pharmcat_reference_fasta_fai    = null
    pharmcat_reporter_sources       = 'CPIC,FDA,DPWG'
    // TODO:

    // MultiQC params
    multiqc_config                  = null
    multiqc_title                   = null
    multiqc_logo                    = null
    max_multiqc_email_size          = '25.MB'
    multiqc_methods_description     = null


    // Boilerplate options
    outdir                          = null
    publish_dir_mode                = 'copy'
    email                           = null
    email_on_fail                   = null
    plaintext_email                 = false
    monochrome_logs                 = false
    hook_url                        = null // TODO: add hook url to nextflow config file, System.getenv('HOOK_URL')
    version                         = false
    help                            = false
    help_full                       = false
    show_hidden                     = false
    trace_report_suffix             = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
    pipelines_testdata_base_path    = 'https://raw.githubusercontent.com/SMD-Bioinformatics-Lund/PrecisionPGx/master/test-datasets/'

    // Config options
    config_profile_name        = null
    config_profile_description = null

    custom_config_version      = 'master'
    // custom_config_base         = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
    custom_config_base         = null
    config_profile_contact     = null
    config_profile_url         = null

    // Enable parameter valiation against defined schema with nf-schema
    validate_params = true
}
```

### Subworkflow config settings

Additional configuration parameters for subworkflows are defined in `*.config` files in the folder `conf/modules/`.

For example, parameters for additional variant filtration based on call quality are defined in config file `conf/modules/variant_filtration.config`.

### Resource use config settings

Resource use settings for processes are defined in the config file `conf/base.config`.

Processes are labelled as:

`process_single`
`process_low`
`process_medium`
`process_high`
`process_long`
`process_medium`
`process_high_memory`

Each label has a defined resource usage in `conf/base.config` 

## Run the pipeline with test data

*Not available yet.*

## Run the pipeline with input data 

### Input samplesheet

The columns that are accepted in the samplesheet are defined in the file `assets/schema_input.json`. 

Required columns are: 

`sample`: sample name 

`type`: T (tumour) or N (normal)

`seq_type`: dna or rna. (required by nf-core modules for hla genotyping with optitype)

For analysis based on FASTQ files `fastq_1` and `fastq_2` must be provided. 

For analysis based on BAM files `bam`: path to BAM file, and `bai`: path to index BAI file, must be provided. 

To generate PharmCAT reports including only a pre-defined set of genes the column `genes` must be provided. `genes` should be a comma-separated list of gene names within citation marks, e.g. `"CYP2C9,CYP2C19,TPMT"`. The samples in the samplesheet can all have different values in the column `genes`.

Example of the contents of a minimal sample sheet:

```
sample,case_id,type,lane,fastq_1,fastq_2,seq_type,genes
NA12878,NA12878,N,L001,/path/to/NA12878-single1_S1_R1_001.fastq.gz,/path/to/NA12878-single1_S1_R2_001.fastq.gz,dna,""
NA12877,NA12878,N,L001,/path/to/NA12877_S3_R1_001.fastq.gz,/path/to/NA12877_S3_R2_001.fastq.gz,dna,"CYP2C9,DPYD,CYP3A4"
```

*To-do: Add table explaining all accepted additional columns*

### Run the pipeline with input parameters in `parameters.yaml`
The input parameters can be set to other values using a `parameters.yaml` file.

Run the pipeline with a `parameters.yaml` file:

```
nextflow run <path to project directory> -profile <container technology, e.g. singularity> --input <input_samplesheet.csv> -params-file <input_parameters.yaml> --outdir <output_dir>`
```

Note that nextflow parameters have a single dash "-", e.g. -profile, -params-file, while pipeline parameters have double dashes "--", e.g. --input, --outdir.


Example of a file `parameters.yaml` with required parameters:

```
#*************************
#Analysis options
analysis_type: panel
platform: illumina

#*************************
#Workflow options
#null/cyp2d6_calling,hla_calling,cnv_calling
skip_subworkflows: hla_calling,cnv_calling

#*************************
#Reference genome parameters
genome: "GRCh38"
#Do not download igenomes
igenomes_ignore: true
#Reference genome file
fasta: /path/to/reference_genome/hg38_no_alt.fa
fai: /path/to/reference_genome/hg38_no_alt.fa.fai

#*************************
#Alignment parameters
aligner: bwamem2
#Path to directory for pre-built bwa index
bwa: /path/to/reference_genome/bwa-idx/
bwamem2: /path/to/reference_genome/bwa-mem2-idx/
#Path to the interval list of the genome. This is used to calculate genome-wide coverage statistics.
intervals_wgs: /path/to/reference_genome/hg38_no_alt.interval_list
#Remove duplicates prior to variant calling
rmdup: true

#*************************
#Target parameters
#The target bed defines the regions where variants are called and where depth is calculated for CNV analysis.
#Applied also when input data is WGS or a broader panel
target_bed: /path/to/project/PrecisionPGx/assets/panpgx.TE-97619858.grch38.refseq_mane.whole_genes_snps.annot.v2.bed
#The amount to pad each end of the target intervals to create bait intervals.
bait_padding: 50

#*************************
#Variant calling parameters
variant_caller: deepvariant
#Filter VCFs by quality and read depth
filter_vcf: true
#Read depth filter threshold, used in additional VCF filtering and PharmCAT genotyping.
min_dp: 20

#*************************
#PharmCAT parameters
pharmcat: true

#Resources
pharmcat_resource_dir: /path/to/project/PrecisionPGx/assets/pharmcat_resources

#Reference genome
pharmcat_reference_fasta: /path/to/pharmcat_reference_genome/HG38_p13_pharmcat/reference.fna.bgz
pharmcat_reference_fasta_index: /path/to/pharmcat_reference_genome//HG38_p13_pharmcat/reference.fna.bgz.gzi
pharmcat_reference_fasta_fai: /path/to/pharmcat_reference_genome/HG38_p13_pharmcat/reference.fna.bgz.fai

#PharmCAT position files
pharmcat_positions: /path/to/project//PrecisionPGx/assets/pharmcat_resources/pharmcat_positions_3.1.1.vcf.gz
pharmcat_positions_index: /path/to/project//PrecisionPGx/assets/pharmcat_resources/pharmcat_positions_3.1.1.vcf.gz.tbi
pharmcat_uniallelic_pos: /path/to/project//PrecisionPGx/assets/pharmcat_resources/pharmcat_positions_3.1.1.uniallelic.vcf.bgz
pharmcat_uniallelic_pos_index: /path/to/project//PrecisionPGx/assets/pharmcat_resources/pharmcat_positions_3.1.1.uniallelic.vcf.bgz.csi

#PharmCAT report options
#Set the sources for prescription recommendations. Defaults to all three: CPIC,FDA,DPWG
pharmcat_reporter_sources: CPIC,FDA,DPWG

#Email
email: user.name@domain.se
email_on_fail: user.name@domain.se
```

### Run the pipeline with user profile and `user_profile.config`
The input parameters can also be set with a custom profile.

The main config file, `nextflow.config`, has a section `profiles` which includes the settings for container technologies (apptainer,docker etc.) and test and user profiles. To add a custom profile, add a line to this section.

Example, with a profile called `smd_lund` defined by a config file called `smd_lund.config`:
```
profiles {
    smd_lund            { includeConfig 'conf/smd_lund.config'              }
}
```

The custom config file should have the same formatting as the main config file, `nextflow.config`. 

Run the pipeline with a user profile:
```
nextflow run <path to project directory> -profile singularity,<user_profile_name> --input <input_samplesheet.csv> --outdir <output_dir>
```

### Nextflow logs
Several different logs are created when running nextflow and can be generated after a run.

See nextflow [report documentation](https://docs.seqera.io/nextflow/reports) for details.

#### Process logs
A general log including input parameters and executed processes is written to standard out by nextflow when running a workflow.

A more detailed process log is written to the file `.nextflow.log`, in the folder from where nextflow is launched. 

#### Command log
Use the command
```
nextflow log <run_name> -f <list of fields>
```
with the field `scripts` included in the list of fields, to generate a log file with the actual commands run in each step.

Example:
```
nextflow log chaotic_easley -f name,script,exit,status
```

To view the log of the latest run, parameter `last` can be used instead of the run name.

```
nextflow log last -f name,script,exit,status
``` 

The run name is given in the standard output or more detailed process log and can also be retrived from the run history file `.nextflow/history`, in the folder from where nextflow is launched, or by running the command
```
nextflow log
``` 

Use a template to specify order and format of the output fields when generating a command log.

Example of a template, *command_log_template.md*:
```
id: $task_id hash: [$hash]
name: $name
module: $module
container: $container

script:
$script

requests:
cpu: $cpus
disk: $disk
memory: $memory
time: $time

exit status: $exit
task status: $status
duration: $duration
task folder: $workdir
------------------------------------------------
```

Run nextflow log using a template:
```
nextflow log chaotic_easley -t config/command_log_template.md
```

With the example template presented here the resulting command log file will then include sections with the following information for each executed task:
```
id: 111 hash: [4a/59d4b1]
name: PRECISIONPGX_MAIN:PRECISIONPGX:CYP2D6_CALLING:PHARMCAT_CYP2D6_MATCHER (NA19109)
module: -
container: /path/to/containers/nxf_nfcore_apptainer_images/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-2b-2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a-data.img

script:

    pharmcat \
        -vcf NA19109.deepvariant.pharmcat.preprocessed.pass.vcf.gz \
        --base-filename NA19109.deepvariant.pharmcat.CYP2D6 \
        --output-dir . \
        -matcher \
        --samples NA19109 \
         --research-mode combinations,cyp2d6 --matcher-all-results --matcher-save-html  \
        --genes CYP2D6


requests:
cpu: 2
disk: -
memory: 12 GB
time: 4h

exit status: 0
task status: COMPLETED
duration: 3.1s
task folder: /labdata/KG/PrecisionPGx/PrecisionPGx_nfcore/work/4a/59d4b1ea36e339ae847d4016dc4b55
------------------------------------------------

id: 110 hash: [71/6c967a]
name: PRECISIONPGX_MAIN:PRECISIONPGX:PHARMCAT_GENOTYPING_REPORTING:PHARMCAT_MATCHER (NA19109)
module: -
container: /path/to/containers/nxf_nfcore_apptainer_images/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-2b-2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a-data.img

...
```

## Pipeline overview

### Prepare references
This subworkflow processes input references and initializes refererence file channels.

### Decompress spring compressed fastq
FASTQ-files compressed with spring, either interleaved or separate, are decompressed.

Decompressed FASTQ-files are not saved as output.

*Not tested*

### Alignment
This subworkflow aligns reads to reference genome with either bwa, bwa-meme,bwa-mem2 or sentieon (bwamem).

#### Aligner
Set aligner with config parameter àligner`.

`aligner`: bwa/bwameme/bwamem2/sentieon

*Bwamem2 tested.*

#### Trim FASTQ
Prior to alignment the subworkflow trims reads with fastp. This can be disabled by including fastp in optional config parameter `skip_tools`.

`skip_tools`: fastp

Fastp runs with default settings for trimming. 

From the [fastp documentation](https://github.com/OpenGene/fastp):
```
Parameter	Description
-q, --qualified_quality_phred	The quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
-u, --unqualified_percent_limit	How many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%
-y or --low_complexity_filter	The complexity is defined as the percentage of base that is different from its next base… its default value is 30, which means 30% complexity is required.
```

Set additional input parameters to fastp in `conf/modules/align.config`. 

Included default settings are:

`-—length_required` The default is set to 40 but can be changed with config parameter min_trimmed_length. 

`--correction` Base mismatches between overlapping paired end reads are corrected using the value and quality of the highest quality called base. Default settings for minimum overlap and quality differences are used. ”This function is based on overlapping detection, which has adjustable parameters overlap_len_require (default 30), overlap_diff_limit (default 5) and overlap_diff_percent_limit (default 20%). Please note that the reads should meet these three conditions simultaneously.”

`--overrepresentation_analysis` Generates counts and analysis of overrepresented sequence 

#### Duplicates
The subworkflow marks duplicate reads with picard MarkDuplicates or sentieon deduplication.

#### Output folder
`alignment`

### QC-BAM
Quality control by analysis of BAM-files.
#### QC statistics  

This subworkflow calculates statistics with `picard`, `tiddit`, `chromograph`, `qualimap` and `mosdepth`.

When input data is from whole genome sequencing the subworkflow calculates additional metrics.

Set input data type with config parameter `analysis_type`.

`analysis_type`: panel/wgs

#### Contamination
The subworkflow estimates contamination with `verifybamid2`.

*Not tested* 

#### Output folder
`qc_bam`

### Variant calling
This subworkflow calls variants with chosen caller.

Set variant caller with config parameter `variant_caller`.

`variant_caller`: deepvariant/sentieon/gatk4.

#### DeepVariant
Calls variants with DeepVariant only in regions specified by input parameter `target_bed`.

DeepVariant uses the WGS model if input data is WGS and WES model if input data is panel data, as set with config parameter `analysis_type`.

`analysis_type`: wgs/panel

(Optional re-genotyping of HET variants in chrX for chrX,chrY is possible with DeepVariant– karyotypes. Not implemented now, as only one gene, G6PD, included in the panel in chrX and none in chrY. G6PD is not covered in all PharmCAT positions)

#### Sentieon
Calls variants with sentieon haplotyper.

#### GATK haplotype caller
Calls variants with GATK haplotype caller.

#### Output folder
`snv_calls`

### Variant filtration
This subworkflow normalises and filters the VCF before input to PharmCAT. 

#### Normalise VCF
Normalise VCF with bcftools.
```
bcftools norm -m+ -c 
```

Default settings:
``` 
-m+: join biallelic sites into multiallelic records 
-c ws: check-references, warn (w) and set/fix (s) bad sites
```

output suffix: `*.norm.vcf.gz`

#### Merge sample VCF with PharmCAT definition VCF
*Removed in this branch. By merging with PharmCAT VCF PharmCAT sites are added and multiallelic sites are annotated as such. Instead of merging with PharmCAT VCF in this step the PharmCAT preprocessor can be run with the setting to add non-variant sites as reference.*

Create two-sample VCF (sample + PharmCAT) with bcftools.

``` 
bcftools merge -m both
```
 
Default settings:
```
-m both: both SNP and indel records can be multiallelic
```

This adds the PharmCAT positions to the VCF, with sample GT ”./.” in positions with no variant.

#### Create target bed filtered VCF
Filter VCF to retain only variants in the regions specified by the target bed with bcftools.
```
bcftools view -s
```
Settings:
```
-s ^PharmCAT: Exclude the PharmCAT sample (This setting is only necessary if the merge step is included)
```

This step is redundant when using deepvariant and variant calling already is limited to the target regions. However, even when the number of variants doesn’t change, this step adds INFO-field annotations for allele count and total number of alleles for each variant. 

output suffix: `*.norm.target.vcf.gz`

#### Create quality filtered VCF
Filter the VCF by call quality and read depth with bcftools.
```
bcftools filter
```

Set additional parameters for filtering on quality and read depth in the config file for the subworkflow: `conf/modules/variant_filtration.config`.

**Sentieon** settings:

QUAL > 10

INFO/DP > `min_dp`

**DeepVariant** settings:

FORMAT/GQ > 19

FORMAT/DP > `min_dp`

For quality filtering of DeepVariant calls FORMAT/GQ is used. 

” The best value to filter on is the GQ value of the variant call itself. We know from empirical investigation that the GQ value is very well-calibrated with the empirical error rate (See Figure 2 of - https://www.nature.com/articles/nbt.4235). This error probability is in the PHRED scale, so a GQ of 10 means DeepVariant indicates a 90% probability the call is correct, while a GQ of 20 indicates a 99% probability the call is correct.” [DeepVariant call filtering issue](https://github.com/google/deepvariant/issues/503)

**GATK haplotyper** settings:
*Not implemented*

output suffix: `*.norm.target.filtered.vcf.gz`

#### Output folder
`snv_calls`

### Target depth
This subworkflow calculates read depth in the regions specified by the input target bed file `target_bed`.

The calculated depths are used to determine the regions where coverage was sufficient to allow for reliable variant detection and for copy number analysis. These regions are included in an output BED-file which is used in subsequent analysis steps.

Set the target depth threshold with the config parameter `min_dp`.

#### Output folder
`target_depth`

### CNV calling
This subworkflow calculates copy numbers using the CNV-Z [algorithm](https://doi.org/10.1016/j.softx.2023.101530).

Expected read depths can be provided with the config parameters or calculated from the analysed cohort.

Skip this subworkflow by adding "cnv_calling" to the config parameter `skip_subworkflows`.

*In progress. Add/test other additional CNV detection tools?* 

#### Output folder
`cnv_calls`

### HLA calling
This subworkflow calls HLA genotypes with [optitype](https://github.com/FRED-2/Optitype).

Skip this subworkflow by adding "hla_calling" to the config parameter `skip_subworkflows`.

*In progress.*

#### Output folder
`hla_calls`

### PharmCAT VCF processing
#### Create PharmCAT formatted VCF

This subworkflow calls the PharmCAT preprocessor to format and filter the input VCF to only include the positions used by PharmCAT.

The PharmCAT preprocessor runs with the following options:
```
--missing-to-ref
```

With this option all positions that does not have a variant in the input VCF or have a variant with genotype ./. are added as homozygous reference 0/0. The option makes it possible to use standard VCF as input, where only variants are present, instead of gVCF, where all positions are included. 

However, with this option an additional filtering step is necessary to distinguish between sites where no variant was called and sites where variants could not be called due to insufficient coverage.

##### Output folder
`pharmcat`

output suffix: `*.pharmcat.preprocessed.vcf`, `*.pharmcat.missing_pgx_var.vcf`

#### Filter PharmCAT VCF by target depth
The subworkflow filters the PharmCAT formatted VCF so that reference PharmCAT positions where coverage was below input threshold are removed. These sites will then be treated by PharmCAT as missing.

Filtering is done using the BED-file generated in the `target depth` subworkflow.

Set the target depth threshold with the config parameter `min_dp`.

##### Output folder
`pharmcat`

Output suffix: `*.pharmcat.preprocessed.pass.vcf`

### PharmCAT genotyping and reporting
This subworkflow does pharmacogenetic genotyping and reporting using PharmCAT.

The subworkflow is called to generate reports including all genes genotyped by PharmCAT and to generate reports including only those genes stated in column `genes` in the samplesheet for a given sample. 

When called to generate a complete report the subworkflow is called as `PHARMCAT_GENOTYPING_REPORTING` and when called to generate a report including only selected genes the subworkflow is called as `PHARMCAT_GENOTYPING_REPORTING_SELECTED`.

Complete report is created by default and activate and inactivate the generation of a selected report is controlled by the column `genes` in the samplesheet for a given sample. 

The genes that are included in the selected report are set in column `genes` in the input samplesheet.

#### PharmCAT matcher
The PharmCAT matcher calls diplotypes from the processed and filtered VCF.

The PharmCAT matcher runs with the following options:
```
--report-all-matches
```
With this option all matching allele combinations for a given gene are included in the summary table in the final report.

##### Output folder
`pharmcat/complete_report` and/or `pharmcat/selected_genes`

Output suffix: `*.pharmcat.match.json` `*.pharmcat.match.html` or `*.pharmcat.selected_genes.match.json` `*.pharmcat.selected_genes.match.html`

The output html report shows the matching defined haplotype combinations and also reports any missing positions.

#### PharmCAT phenotyper
The PharmCAT phenotyper matches the diplotypes with annotated phenotypes.

##### Output folder
`pharmcat/complete_report` and/or `pharmcat/selected_genes`

Output suffix: `*.pharmcat.phenotype.json` or `*.pharmcat.selected_genes.phenotype.json`

#### PharmCAT reporter
The PharmCAT reporter creates reports including prescription recommendations for the obtained phenotypes.

Specify the sources for prescription recommendations with the config parameter `pharmcat_reporter_sources`. Default is all three CPIC,FDA,DPWG.

`pharmcat_reporter_sources`: CPIC,FDA,DPWG

Note: The gene CYP3A4 only have recommendations from DPWG.

##### Output folder
`pharmcat/complete_report` and/or `pharmcat/selected_genes`

Output suffix: `*.pharmcat.report.json` `*.pharmcat.report.tsv` `*.pharmcat.report.html`

### CYP2D6 calling
Currently only calls *CYP2D6* with PharmCAT.

This subworkflow calls the PharmCAT matcher and phenotyper to genotype *CYP2D6* based only on detected SNVs. This will fail to detect copy number variants and hybrid genes and cannot be used alone to call *CYP2D6* genotypes.

Skip this subworkflow by adding "cnv_calling" to the config parameter `skip_subworkflows`.

`skip_subworkflows`: cnv_calling

#### Output folder
`CYP2D6`

Output suffix: `*.pharmcat.CYP2D6.phenotype.json` `*.pharmcat.CYP2D6.match.json` `*.pharmcat.CYP2D6.match.html`

*In progress. Must also include CNV calls to appropriately detect copy number variants and hybrid genes.*

### Collect software versions & multiqc
Collects software versions for the used programs and creates the MultiQC report.

#### MultiQC
This subworkflow collects data and creates a MultiQC report.

#### Output folder
`multiqc`
