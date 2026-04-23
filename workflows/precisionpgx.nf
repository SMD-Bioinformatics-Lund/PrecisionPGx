/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap;samplesheetToList } from 'plugin/nf-schema'
include { paramsSummaryMultiqc               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText             } from '../subworkflows/local/utils_nfcore_precisionpgx_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R1_FQ   } from '../modules/nf-core/spring/decompress/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R2_FQ   } from '../modules/nf-core/spring/decompress/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_FQ_PAIR } from '../modules/nf-core/spring/decompress/main'

//
// SUBWORKFLOWS
//

include { ALIGN                                                                     } from '../subworkflows/local/align'
include { PREPARE_REFERENCES                                                        } from '../subworkflows/local/prepare_references'
include { QC_BAM                                                                    } from '../subworkflows/local/qc_bam'
include { VARIANT_CALLING                                                           } from '../subworkflows/local/variant_calling'
include { VARIANT_FILTRATION                                                        } from '../subworkflows/local/variant_filtration'
include { VARIANT_FILTRATION as GVCF_FILTRATION                                     } from '../subworkflows/local/variant_filtration'
include { PHARMCAT_VCF_PROCESSING                                                   } from '../subworkflows/local/pharmcat_vcf_processing'
include { PHARMCAT_GENOTYPING_REPORTING                                             } from '../subworkflows/local/pharmcat_genotyping_reporting'
include { PHARMCAT_GENOTYPING_REPORTING as PHARMCAT_GENOTYPING_REPORTING_SELECTED   } from '../subworkflows/local/pharmcat_genotyping_reporting'
include { TARGET_DEPTH                                                              } from '../subworkflows/local/target_depth'
include { CYP2D6_CALLING                                                            } from '../subworkflows/local/cyp2d6_calling'
//include { CNV_CALLING                                                             } from '../subworkflows/local/cnv_calling'
//include { HLA_CALLING                                                             } from '../subworkflows/local/hla_calling'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PRECISIONPGX {

    take:
    ch_reads
    ch_alignments
    ch_samples

    main:

    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }



    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Initialize file channels for PREPARE_REFERENCES subworkflow
    //
    ch_genome_fasta              = Channel.fromPath(params.fasta).map { it -> [[id:it.simpleName], it] }.collect()
    ch_genome_fai                = params.fai                 ? Channel.fromPath(params.fai).map {it -> [[id:it.simpleName], it]}.collect()
                                                                : Channel.empty()
    ch_genome_dictionary         = params.sequence_dictionary ? Channel.fromPath(params.sequence_dictionary).map {it -> [[id:it.simpleName], it]}.collect()
                                                                : Channel.empty()
    ch_target_bed_unprocessed    = params.target_bed          ? Channel.fromPath(params.target_bed).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                : Channel.value([[],[]])

    //
    // Prepare references and indices.
    //
    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary,
        ch_target_bed_unprocessed,
    )
    .set { ch_references }

    //
    // Gather built indices or get them from the params
    //
    ch_bait_intervals               = ch_references.bait_intervals
    ch_target_bed                   = ch_references.target_bed
    ch_target_intervals             = ch_references.target_intervals
    ch_target_bed_uncompressed      = ch_references.target_bed_uncompressed
    ch_call_interval                = params.call_interval                      ? Channel.fromPath(params.call_interval).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])
    ch_genome_bwaindex              = params.bwa                                ? Channel.fromPath(params.bwa).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : ch_references.genome_bwa_index
    ch_genome_bwamem2index          = params.bwamem2                            ? Channel.fromPath(params.bwamem2).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : ch_references.genome_bwamem2_index
    ch_genome_bwamemeindex          = params.bwameme                            ? Channel.fromPath(params.bwameme).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : ch_references.genome_bwameme_index
    ch_genome_chrsizes              = ch_references.genome_chrom_sizes
    ch_genome_fai                   = ch_references.genome_fai
    ch_genome_dictionary            = ch_references.genome_dict
    ch_ml_model                     = params.variant_caller.equals('sentieon')  ? Channel.fromPath(params.ml_model).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])

    ch_intervals_wgs                = params.intervals_wgs                      ? Channel.fromPath(params.intervals_wgs).collect()
                                                                                : Channel.empty()
    ch_svd_bed                      = params.verifybamid_svd_bed                ? Channel.fromPath(params.verifybamid_svd_bed)
                                                                                : Channel.empty()
    ch_svd_mu                       = params.verifybamid_svd_mu                 ? Channel.fromPath(params.verifybamid_svd_mu)
                                                                                : Channel.empty()
    ch_svd_ud                       = params.verifybamid_svd_ud                 ? Channel.fromPath(params.verifybamid_svd_ud)
                                                                                : Channel.empty()

    ch_sentieon_emit_vcf            = params.emit_mode.equals('gvcf')           ? Channel.value(false) : Channel.value(params.emit_mode)
    ch_sentieon_emit_gvcf           = params.emit_mode.equals('gvcf')           ? Channel.value(params.emit_mode) : Channel.value(false)

    // Pharmcat 
    ch_pc_positions_vcf             = params.pharmcat_positions                 ? Channel.fromPath(params.pharmcat_positions).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])
    ch_pc_positions_vcf_index       = params.pharmcat_positions_index           ? Channel.fromPath(params.pharmcat_positions_index).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])

    ch_pc_uniallelic_pos_vcf        = params.pharmcat_uniallelic_pos            ? Channel.fromPath(params.pharmcat_uniallelic_pos).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])
    ch_pc_uniallelic_pos_vcf_index  = params.pharmcat_uniallelic_pos_index      ? Channel.fromPath(params.pharmcat_uniallelic_pos_index).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])

    ch_pc_reference_fasta           = params.pharmcat_reference_fasta           ? Channel.fromPath(params.pharmcat_reference_fasta).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])
    ch_pc_reference_fasta_index     = params.pharmcat_reference_fasta_index     ? Channel.fromPath(params.pharmcat_reference_fasta_index).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])
    ch_pc_reference_fasta_fai       = params.pharmcat_reference_fasta_fai       ? Channel.fromPath(params.pharmcat_reference_fasta_fai).map {it -> [[id:it.simpleName], it]}.collect()
                                                                                : Channel.value([[:],[]])


    //
    // Input QC (ch_reads will be empty if fastq input isn't provided so FASTQC won't run if input is not fastq)
    //

    ch_input_by_sample_type = ch_reads.branch{
        fastq_gz:           it[0].data_type == "fastq_gz"
        interleaved_spring: it[0].data_type == "interleaved_spring"
        separate_spring:    it[0].data_type == "separate_spring"
    }

    // Just one fastq.gz.spring-file with both R1 and R2
    ch_one_fastq_gz_pair_from_spring = SPRING_DECOMPRESS_TO_FQ_PAIR(ch_input_by_sample_type.interleaved_spring, false).fastq

    // Two fastq.gz.spring-files - one for R1 and one for R2
    ch_r1_fastq_gz_from_spring  = SPRING_DECOMPRESS_TO_R1_FQ(ch_input_by_sample_type.separate_spring.map{ meta, files -> [meta, files[0] ]}, true).fastq
    ch_r2_fastq_gz_from_spring  = SPRING_DECOMPRESS_TO_R2_FQ(ch_input_by_sample_type.separate_spring.map{ meta, files -> [meta, files[1] ]}, true).fastq
    ch_two_fastq_gz_from_spring = ch_r1_fastq_gz_from_spring.join(ch_r2_fastq_gz_from_spring).map{ meta, fastq_1, fastq_2 -> [meta, [fastq_1, fastq_2]]}

    ch_input_fastqs = ch_input_by_sample_type.fastq_gz.mix(ch_one_fastq_gz_pair_from_spring).mix(ch_two_fastq_gz_from_spring)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ALIGN & FETCH STATS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ALIGN (
        ch_input_fastqs,
        ch_alignments,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        ch_genome_bwamemeindex,
        ch_genome_dictionary,
        params.mbuffer_mem,
        params.platform,
        params.samtools_sort_threads
    )
    .set { ch_mapped }

    //
    // BAM QUALITY CHECK
    //
    QC_BAM (
        ch_mapped.genome_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary,
        ch_target_bed.map {
            meta, bed_path, bed_tbi -> [ meta, bed_path ]
        },
        ch_target_bed_uncompressed,
        ch_bait_intervals,
        ch_target_intervals,
        ch_intervals_wgs,
        ch_svd_bed,
        ch_svd_mu,
        ch_svd_ud,
    )


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        HAPLOTYPE CALLING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    //
    // VARIANT CALLING
    //

    VARIANT_CALLING (
        ch_mapped.genome_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary,
        ch_target_bed.map {
            meta, bed_path, bed_tbi -> [ meta, bed_path ]
        },
        ch_target_intervals,
        ch_sentieon_emit_vcf,
        ch_sentieon_emit_gvcf
    )
    .set { ch_haplotypes }


    ch_filter_vcf_input = Channel.empty().mix(
        ch_haplotypes.sentieon_vcf,
        ch_haplotypes.gatk_vcf,
        ch_haplotypes.deepvariant_vcf
        )
    
    ch_filter_vcf_input_tbi = Channel.empty().mix(
        ch_haplotypes.sentieon_vcf_tbi,
        ch_haplotypes.gatk_vcf_tbi,
        ch_haplotypes.deepvariant_vcf_tbi
        )


    ch_filter_vcf_input.join(
        ch_filter_vcf_input_tbi, 
        failOnMismatch:true, 
        failOnDuplicate:true
    ).set { ch_filter_vcf_input_joined }

    //Create gvcf channel
    //ch_haplotypes.deepvariant_gvcf.join(
    //    ch_haplotypes.deepvariant_gvcf_tbi,
    //    failOnMismatch:true,
    //    failOnDuplicate:true
    //).set { ch_filter_gvcf_input_joined }

    VARIANT_FILTRATION (
        ch_filter_vcf_input_joined,
        ch_genome_fasta,
        ch_genome_fai,
        ch_target_bed.map {
            meta, bed_path, bed_tbi -> [ meta, bed_path ]
        },
        ch_pc_positions_vcf.map {
            meta, pc_vcf -> [ pc_vcf ]
        },
        ch_pc_positions_vcf_index.map {
            meta, pc_vcf_tbi -> [ pc_vcf_tbi ]
        }
    )
    .set { ch_filtered_haplotypes }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CNV CALLING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('cnv_calling'))){
        // Input QC_BAM.out.target_depth_tsv
        /*
        CNV_CALLING (
            QC_BAM.out.target_depth_tsv
        )
        .set { ch_cnvcalls }
        */
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        HLA CALLING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('hla_calling'))){
        // Input ch_mapped.genome_bam_bai
        /*
        HLA_CALLING (
            ch_mapped.genome_bam_bai
        )
        .set { ch_hlacalls }
        */
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PHARMCAT VCF PROCESSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // PHARMCAT
    //
    ch_filtered_haplotypes.filtered_vcf.join(
        ch_filtered_haplotypes.filtered_vcf_tbi,
        failOnMismatch:true,
        failOnDuplicate:true
    ).set { ch_pharmcat_input_joined }

    PHARMCAT_VCF_PROCESSING(
        ch_pharmcat_input_joined, 
        ch_pc_reference_fasta, 
        ch_pc_reference_fasta_fai, 
        ch_pc_positions_vcf.join(
            ch_pc_positions_vcf_index, 
            failOnMismatch:true, 
            failOnDuplicate:true
        ),
        ch_pc_uniallelic_pos_vcf.join(
            ch_pc_uniallelic_pos_vcf_index, 
            failOnMismatch:true, 
            failOnDuplicate:true
        ),
        QC_BAM.out.target_pass_bed
    )
    .set { ch_pharmcat }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PHARMCAT GENOTYPING REPORTING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_pharmcat.preprocessed_vcf_pass.join(
        ch_pharmcat.preprocessed_vcf_pass_tbi,
        failOnMismatch:true,
        failOnDuplicate:true
    ).set { ch_pc_input }

    //Generate complete report
    PHARMCAT_GENOTYPING_REPORTING(
        ch_pc_input,
        [] // This is by defualt empty because we want the complete report
    )
    .set { ch_pharmcat_complete }

    //Generate report with selected genes
    PHARMCAT_GENOTYPING_REPORTING_SELECTED(
        ch_pc_input,
        ch_pc_input.map {
            meta, vcf, tbi -> meta.genes
        } // Here we send the meta.genes for selected genes report
    )
    .set { ch_pharmcat_selected }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CYP2D6 CALLING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('cyp2d6_calling'))){

        CYP2D6_CALLING(
            ch_pharmcat.preprocessed_vcf_pass.join(
                ch_pharmcat.preprocessed_vcf_pass_tbi,
                failOnMismatch:true,
                failOnDuplicate:true
            )
        )
        .set { ch_cyp2d6 }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT SOFTWARE VERSIONS & MultiQC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // TASK: Aggregate software versions

    //
    /*
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
    */

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'precisionpgx_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.fromPath("$projectDir/assets/PrecisionPGx-light.png", checkIfExists: true)


    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.fastp_json.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.markdup_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.hs_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.self_sm.map{it[1]}.collect().ifEmpty([]))


    ch_multiqc_configs_list = ch_multiqc_config
        .mix(ch_multiqc_custom_config)
        .collect()
        .ifEmpty([])
    ch_multiqc_logo_list = ch_multiqc_logo
        .collect()
        .ifEmpty([])

    ch_multiqc_files_keyed = ch_multiqc_files
        .collect()
        .map { files -> ['multiqc', files] }
    ch_multiqc_configs_keyed = ch_multiqc_configs_list
        .map { configs -> ['multiqc', configs] }
    ch_multiqc_logo_keyed = ch_multiqc_logo_list
        .map { logo -> ['multiqc', logo] }

    ch_multiqc_input = ch_multiqc_files_keyed
        .join(ch_multiqc_configs_keyed)
        .join(ch_multiqc_logo_keyed)
        .map { _, files, configs, logo ->
            [[id: 'multiqc'], files, configs, logo, [], []]
        }


    MULTIQC (
        ch_multiqc_input
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()        // channel: /path/to/multiqc_report.html
    versions       = ch_versions                        // channel: [ path(versions.yml) ]

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
