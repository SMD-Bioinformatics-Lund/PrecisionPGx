//
// Subworkflow with functionality specific to the nf-core/precisionpgx pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Validate parameters against schema at runtime
    monochrome_logs   // boolean: Disable ANSI colours
    nextflow_cli_args //   array: Positional nextflow CLI args
    outdir            //  string: Output directory
    input             //  string: Path to input samplesheet

    main:

    //
    // Print version and exit if required + dump params json + conda/mamba handling
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    
    // Validate parameters against schema at runtime
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        null,
        null,
        null,
        null,
        null,
        null,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(nextflow_cli_args)

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()
    checkRequiredParameters(params)


    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .tap { ch_original_input }
        .map { meta, fastq1, fastq2, spring1, spring2, bam, bai -> meta.id }
        .reduce([:]) { counts, sample -> //get counts of each sample in the samplesheet - for groupTuple
            counts[sample] = (counts[sample] ?: 0) + 1
            counts
        }
        .combine( ch_original_input )
        .map { counts, meta, fastq1, fastq2, spring1, spring2, bam, bai ->
            def new_meta = meta + [num_lanes:counts[meta.id]]
            if (fastq1 && fastq2) {
                new_meta += [read_group: generateReadGroupLine(fastq1, meta, params)]
                return [new_meta + [single_end: false, data_type: "fastq_gz"], [fastq1, fastq2]]
            } else if (fastq1 && !fastq2) {
                new_meta += [read_group: generateReadGroupLine(fastq1, meta, params)]
                return [new_meta + [single_end: true, data_type: "fastq_gz"], [fastq1]]
            } else if (spring1 && spring2) {
                new_meta += [read_group: generateReadGroupLine(spring1, meta, params)]
                return [new_meta + [single_end: false, data_type: "separate_spring"], [spring1, spring2]]
            } else if (spring1 && !spring2) {
                new_meta += [read_group: generateReadGroupLine(spring1, meta, params)]
                return [new_meta + [single_end: false, data_type: "interleaved_spring"], [spring1]]
            } else if (bam && bai) {
                new_meta += [read_group: generateReadGroupLine(bam, meta, params)]
                return [new_meta, [bam, bai]]
            }
        }
        .tap{ ch_input_counts }
        .map { meta, files -> files }
        .reduce([:]) { counts, files -> //get line number for each row to construct unique sample ids
            counts[files] = counts.size() + 1
            return counts
        }
        .combine( ch_input_counts )
        .map { lineno, meta, files -> //append line number to sampleid
            def new_meta = meta + [id:meta.id+"_LNUMBER"+lineno[files]]
            return [ new_meta, files ]
        }
        .tap { ch_samplesheet }
        .branch { meta, files  ->
            fastq: !files[0].toString().endsWith("bam")
                return [meta, files]
            align: files[0].toString().endsWith("bam")
                return [meta, files]
        }
        .set {ch_samplesheet_by_type}

    ch_samples  = ch_samplesheet.map { meta, files ->
                    def new_id = meta.sample
                    def new_meta = meta - meta.subMap('lane', 'read_group') + [id:new_id]
                    return new_meta
                    }.unique()


    emit:
    reads     = ch_samplesheet_by_type.fastq
    align     = ch_samplesheet_by_type.align
    samples   = ch_samples
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           // string: completion email
    email_on_fail   // string: email on fail
    plaintext_email // boolean: plain-text emails
    outdir          // path: output directory
    monochrome_logs // boolean: disable ANSI colour
    hook_url        // string: hook url for notifications
    multiqc_report  // path: multiqc report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    // Send completion email and summary
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.value,
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) imNotification(summary_params, hook_url)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def generateReadGroupLine(file, meta, params) {
    return "\'@RG\\tID:" + file.simpleName + "_" + meta.type + "\\tPL:" + params.platform.toUpperCase() + "\\tSM:" + meta.id + "\'"
}

def validateInputParameters() {
    genomeExistsError()

    // Basic required checks (schema likely enforces, but keep hard checks)
    if (!params.input)  error "Missing required parameter: --input"
    if (!params.outdir) error "Missing required parameter: --outdir"

    // If not using --genome, user must provide the core references explicitly
    // (you can adjust required set depending on Sentieon needs)
    if (!(params.genomes && params.genome)) {
        if (!params.fasta) error "Missing --fasta (or supply --genome with --genomes config)"
    }
}


def checkRequiredParameters(params) {
    def mandatoryParams = [
        "analysis_type",
        "fasta",
        "input",
        "outdir",
        "pharmcat_reporter_sources",
        "pharmcat_complete_report",
        "pharmcat_selected_report",
    ]

    def pharmcatParams = [
        "pharmcat_resource_dir",
        "pharmcat_positions",
        "pharmcat_positions_index",
        "pharmcat_uniallelic_pos",
        "pharmcat_uniallelic_pos_index",
        "pharmcat_reference_fasta",
        "pharmcat_reference_fasta_index",
        "pharmcat_reference_fasta_fai",
    ]

    mandatoryParams += pharmcatParams

    // Static requirements that are not influenced by user-defined skips
    def staticRequirements   = [
        analysis_type_panel      : ["target_bed"],
        variant_caller_sentieon  : ["ml_model"],
    ]

    // Requirements that can be modified by the user using either skip_tools or skip_subworkflows here
    def dynamicRequirements = [
        variant_calling              : ["genome"],
        variant_annotation           : ["genome"],
    ]

    def missingParamsCount = 0

    staticRequirements.each { condition, paramsList ->
        if ((condition == "analysis_type_panel" && params.analysis_type == "panel") ||
            (condition == "variant_caller_sentieon" && params.variant_caller.equals('sentieon'))) {
                mandatoryParams += paramsList
        }
    }

    all_skips = params.skip_subworkflows+","+params.skip_tools
    dynamicRequirements.each { condition, paramsList ->
        if (!all_skips.split(',').contains(condition)) {
                mandatoryParams += paramsList
        }
    }


    mandatoryParams.unique().each { param ->
        if (params[param] == null) {
            println("params." + param + " not set.")
            missingParamsCount += 1
        }
    }

    pharmcatParams.unique().each { param ->
        def p = file(params[param])
        if (!p.exists()) error "PharmCAT resource file does not exist: ${params[param]}"
    }

    if (missingParamsCount > 0) {
        error("\nSet missing parameters and restart the run. For more information please check usage documentation on github.")
    }
}


//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}


//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}


def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Available genome keys are:\n" +
            "  ${params.genomes.keySet().join(', ')}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def align_text                  = []
    def variant_annotation_text   = []
    def haplotype_calls_text        = []
    def qc_bam_text                 = []
    def preprocessing_text          = []
    def other_citation_text         = []

    align_text = [
        params.aligner.equals("bwa")      ? "BWA (Li, 2013),"                        :"",
        params.aligner.equals("bwamem2")  ? "BWA-MEM2 (Vasimuddin et al., 2019),"    : "",
        params.aligner.equals("bwameme")  ? "BWA-MEME (Jung et al., 2022),"          : "",
        params.aligner.equals("sentieon") ? "Sentieon DNASeq (Kendig et al., 2019)," : "",
        params.aligner.equals("sentieon") ? "Sentieon Tools (Freed et al., 2017),"   : ""
    ]

    // TODO:
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('variant_annotation'))) {
        variant_annotation_text = [
            "CADD (Rentzsch et al., 2019, 2021),",
            "Vcfanno (Pedersen et al., 2016),",
            "VEP (McLaren et al., 2016),",
            "Genmod (Magnusson et al., 2018),"
        ]
    }

    // TODO:
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('variant_calling'))) {
        haplotype_calls_text = [
            params.variant_caller.equals('gatk4') ? "GATK (McKenna et al., 2010),"      : "",
            params.variant_caller.equals('sentieon')    ? "Sentieon DNAscope (Freed et al., 2022)," : "",
        ]
    }

    // TODO:
    qc_bam_text = [
        "Picard (Broad Institute, 2023)",
        "Qualimap (Okonechnikov et al., 2016),",
        "TIDDIT (Eisfeldt et al., 2017),",
        "UCSC Bigwig and Bigbed (Kent et al., 2010),",
        (params.verifybamid_svd_bed && params.verifybamid_svd_mu && params.verifybamid_svd_ud) ? "VerifyBamID2 (Zhang et al., 2020)," : "",
        "Mosdepth (Pedersen & Quinlan, 2018),"
    ]

    // TODO: Seqtk
    preprocessing_text = [
        "FastQC (Andrews 2010),",
        (params.skip_tools && params.skip_tools.split(',').contains('seqtk')) ? "" : "Fastp (Chen, 2023),"
    ]

    // TODO:
    other_citation_text = [
        "BCFtools (Danecek et al., 2021),",
        "BEDTools (Quinlan & Hall, 2010),",
        "GATK (McKenna et al., 2010),",
        "MultiQC (Ewels et al. 2016),",
        "SAMtools (Li et al., 2009),",
        "Tabix (Li, 2011)",
        "."
    ]

    def concat_text = align_text +
                        variant_annotation_text   +
                        haplotype_calls_text        +
                        qc_bam_text                 +
                        preprocessing_text          +
                        other_citation_text

    def citation_text = [ "Tools used in the workflow included:" ] + concat_text.unique(false) { a, b -> a <=> b } - ""
    return citation_text.join(' ').trim()
}

def toolBibliographyText() {

    def align_text                  = []
    def variant_annotation_text   = []
    def haplotype_calls_text        = []
    def qc_bam_text                 = []
    def preprocessing_text          = []
    def other_citation_text         = []

    align_text = [
        params.aligner.equals("bwa") ? "<li>Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM (arXiv:1303.3997). arXiv. http://arxiv.org/abs/1303.3997</li>" :"",
        params.aligner.equals("bwamem2") ? "<li>Vasimuddin, Md., Misra, S., Li, H., & Aluru, S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), 314–324. https://doi.org/10.1109/IPDPS.2019.00041</li>" : "",
        params.aligner.equals("bwameme") ? "<li>Jung Y, Han D. BWA-MEME: BWA-MEM emulated with a machine learning approach. Bioinformatics. 2022;38(9):2404-2413. doi:10.1093/bioinformatics/btac137</li>" : "",
        params.aligner.equals("sentieon") ? "<li>Kendig, K. I., Baheti, S., Bockol, M. A., Drucker, T. M., Hart, S. N., Heldenbrand, J. R., Hernaez, M., Hudson, M. E., Kalmbach, M. T., Klee, E. W., Mattson, N. R., Ross, C. A., Taschuk, M., Wieben, E. D., Wiepert, M., Wildman, D. E., & Mainzer, L. S. (2019). Sentieon DNASeq Variant Calling Workflow Demonstrates Strong Computational Performance and Accuracy. Frontiers in Genetics, 10, 736. https://doi.org/10.3389/fgene.2019.00736</li>" : "",
        params.aligner.equals("sentieon") ? "<li>Freed, D., Aldana, R., Weber, J. A., & Edwards, J. S. (2017). The Sentieon Genomics Tools—A fast and accurate solution to variant calling from next-generation sequence data (p. 115717). bioRxiv. https://doi.org/10.1101/115717</li>" : ""
    ]

    // TODO:
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('variant_annotation'))) {
        variant_annotation_text = [
            "<li>Rentzsch, P., Schubach, M., Shendure, J., & Kircher, M. (2021). CADD-Splice—Improving genome-wide variant effect prediction using deep learning-derived splice scores. Genome Medicine, 13(1), 31. https://doi.org/10.1186/s13073-021-00835-9</li>",
            "<li>Rentzsch, P., Witten, D., Cooper, G. M., Shendure, J., & Kircher, M. (2019). CADD: Predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Research, 47(D1), D886–D894. https://doi.org/10.1093/nar/gky1016</li>",
            "<li>Pedersen, B. S., Layer, R. M., & Quinlan, A. R. (2016). Vcfanno: Fast, flexible annotation of genetic variants. Genome Biology, 17(1), 118. https://doi.org/10.1186/s13059-016-0973-5</li>",
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Magnusson, M., Hughes, T., Glabilloy, & Bitdeli Chef. (2018). genmod: Version 3.7.3 (3.7.3) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.3841142</li>"
        ]
    }
    // TODO:
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('variant_calling'))) {
        haplotype_calls_text = [
            params.variant_caller.equals('gatk4') ? "<li>Poplin, R., Chang, P.-C., Alexander, D., Schwartz, S., Colthurst, T., Ku, A., Newburger, D., Dijamco, J., Nguyen, N., Afshar, P. T., Gross, S. S., Dorfman, L., McLean, C. Y., & DePristo, M. A. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983–987. https://doi.org/10.1038/nbt.4235</li>" : "",
            params.variant_caller.equals('sentieon') ? "<li>Freed, D., Pan, R., Chen, H., Li, Z., Hu, J., & Aldana, R. (2022). DNAscope: High accuracy small variant calling using machine learning [Preprint]. Bioinformatics. https://doi.org/10.1101/2022.05.20.492556</li>" : ""
        ]
    }

    qc_bam_text = [
        "<li>Broad Institute. (2023). Picard Tools. In Broad Institute, GitHub repository. http://broadinstitute.github.io/picard/</li>",
        "<li>Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: Advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2), 292–294. https://doi.org/10.1093/bioinformatics/btv566</li>",
        "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
        "<li>Kent, W. J., Zweig, A. S., Barber, G., Hinrichs, A. S., & Karolchik, D. (2010). BigWig and BigBed: Enabling browsing of large distributed datasets. Bioinformatics, 26(17), 2204–2207. https://doi.org/10.1093/bioinformatics/btq351</li>",
        "<li>Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: Quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867–868. https://doi.org/10.1093/bioinformatics/btx699</li>"
    ]
    // TODO:
    preprocessing_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        (params.skip_tools && params.skip_tools.split(',').contains('seqtk')) ? "" : "<li>Chen, S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta, 2(2), e107. https://doi.org/10.1002/imt2.107</li>"
    ]

    // TODO:
    other_citation_text = [
        "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008</li>",
        "<li>McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354</li>",
        "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352</li>",
        "<li>Li, H. (2011). Tabix: Fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics, 27(5), 718–719. https://doi.org/10.1093/bioinformatics/btq671</li>",
        "<li>Quinlan, AR., Hall IM. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinfomatics, 26(6), 841-842. https://doi.org/10.1093/bioinformatics/btq033</li>"
    ]

    def concat_text = align_text +
                        variant_annotation_text   +
                        haplotype_calls_text        +
                        qc_bam_text                 +
                        preprocessing_text          +
                        other_citation_text

    def reference_text = concat_text.unique(false) { a, b -> a <=> b } - ""
    return reference_text.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}