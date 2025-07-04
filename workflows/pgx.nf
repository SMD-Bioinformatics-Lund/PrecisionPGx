#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { SAMPLE                        } from '../subworkflows/local/sample'
include { ALIGN_SENTIEON                } from '../subworkflows/local/align_sentieon'
include { QC                            } from '../subworkflows/local/qc'
include { PADDED_INTERVALS              } from '../subworkflows/local/padded_intervals'
include { HAPLOTYPING                   } from '../subworkflows/local/haplotyping'
include { ONTARGET                      } from '../subworkflows/local/ontarget'
include { ANNOTATION                    } from '../subworkflows/local/annotations'
include { PHARMCAT                      } from '../subworkflows/local/pharmcat'
include { COVERAGE                      } from '../subworkflows/local/coverage'
include { CLINICAL_INFORMATION          } from '../subworkflows/local/clinical_information'
include { PGX_REPORT                    } from '../subworkflows/local/reports'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'
include { RUN_MULTIQC                   } from '../subworkflows/local/multiqc'

csv = file(params.csv)

workflow PGX_FULL {

    take:
        fastq_input
        samples

    main:
        ch_versions = Channel.empty()

        // ALIGNMENT
        ALIGN_SENTIEON ( fastq_input ).set { ch_mapped }
        ch_versions = ch_versions.mix(ch_mapped.versions)

        QC ( 
            ch_mapped.qc_out, 
            ch_mapped.bam_lowcov
        ).set { ch_qc }
        ch_versions = ch_versions.mix(ch_qc.versions)

        aligned_bam = ch_mapped.bam_umi.map { it -> it[0..3] }

        // HAPLOTYPING
        HAPLOTYPING ( aligned_bam ).set { ch_haplotypes }
        ch_versions = ch_versions.mix(ch_haplotypes.versions)

        // PHARMCAT
        PHARMCAT ( HAPLOTYPING.out.filtered_haplotypes ).set { ch_pharmcat }
        ch_versions = ch_versions.mix(ch_pharmcat.versions)

        // ONTARGET
        ONTARGET ( 
            aligned_bam, 
            ch_haplotypes.filtered_haplotypes
        ).set { ch_ontarget }


        if (params.ontarget) {
            ontarget_bams           = ch_ontarget.ontarget_bam
            ontarget_haplotypes     = ch_ontarget.ontarget_vcf
        } else {
            ontarget_bams           = aligned_bam
            ontarget_haplotypes     = ch_haplotypes.filtered_haplotypes
        }

        // ANNOTATION
        ANNOTATION ( ontarget_haplotypes ).set { ch_annotation }

        // COVERAGE
        COVERAGE ( 
            ontarget_bams, // CAN be ontarget or all, modify this to supprt both
        ).set { ch_coverage }

        // CLINICAL_INFORMATION
        CLINICAL_INFORMATION ( ch_annotation.detected_variants )

        // PGX REPORT
        PGX_REPORT ( 
            ch_annotation.annotated_haplotypes, 
            ch_annotation.detected_variants,
            ch_coverage.missing_targets_coverage_annotated,
            CLINICAL_INFORMATION.out.diplotypes,
            ch_coverage.baits_coverage,
            CLINICAL_INFORMATION.out.interactions
        ).set { ch_report }
        ch_versions = ch_versions.mix(ch_report.versions)

        // MODULE: Software Versions

        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml'),
            samples
        )


        // SUBWORKFLOW: MULTIQC
        fastqc_files = Channel.empty()
        // fastqc_files = fastqc_files.mix(ch_mapped.fastqc_out)

        RUN_MULTIQC (
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_report.targets_depth,
            fastqc_files,
        ) .set { ch_multiqc }


    emit:
        umi_bam                     = ch_mapped.bam_umi                                 // channel: [ tuple val(group), val(meta), file("aligned.bam"), file("aligned.bam.bai") ]
        report                      = ch_report.pgx_report                              // channel: [ val(group), val(meta), file(pgx-report) ]
        pharmcat_preprocessed       = ch_pharmcat.pharmcat_preprocessed                 // channel: [ tuple val(group), val(meta) file("pharmcat_preprocessed.vcf") ]
        pharmcat_report             = ch_pharmcat.pharmcat_report                       // channel: [ tuple val(group), val(meta) file("pharmcat.html") ]
        pharmcat_pheno_json         = ch_pharmcat.pharmcat_pheno_json                   // channel: [ tuple val(group), val(meta) file("pharmcat.phenotype.json") ]
        pharmcat_match_json         = ch_pharmcat.pharmcat_match_json                   // channel: [ tuple val(group), val(meta) file("pharmcat.match.json") ]
        pharmcat_match_html         = ch_pharmcat.pharmcat_match_html                   // channel: [ tuple val(group), val(meta) file("pharmcat.match.html") ]
        pharmcat_report_json        = ch_pharmcat.pharmcat_report_json                  // channel: [ tuple val(group), val(meta) file("pharmcat.match.html") ]
        multiqc_report              = ch_multiqc.multiqc_report                         // channel: [ tuple val(group), val(meta) file("multiqc.html") ]
        multiqc_data                = ch_multiqc.multiqc_data                           // channel: [ tuple val(group), val(meta) file("multiqc_data") ]
        multiqc_plots               = ch_multiqc.multiqc_data                           // channel: [ tuple val(group), val(meta) file("multiqc_plots") ]
        versions                    = ch_versions                                       // channel: [ file(versions) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	logFile = file("${params.crondir}/logs/" + base + "pgx.complete")
	logFile.text = msg
	logFile.append(error)
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def concatenate_values(channel) {
    channel.collect().join('_')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/