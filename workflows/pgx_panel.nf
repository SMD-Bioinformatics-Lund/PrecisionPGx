#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPROCESS_FASTQ              } from '../subworkflows/local/preprocess_fastq'
include { ALIGN                         } from '../subworkflows/local/alignment'
include { VARIANT_CALLING               } from '../subworkflows/local/variant_calling'
include { VARIANT_FILTRATION            } from '../subworkflows/local/variant_filtration'
include { COVERAGE                      } from '../subworkflows/local/coverage'
include { PHARMCAT                      } from '../subworkflows/local/pharmcat'
include { CYP2D6_CNVCALL                } from '../subworkflows/local/cnv'



// include { QC                            } from '../subworkflows/local/qc'
// include { PADDED_INTERVALS              } from '../subworkflows/local/padded_intervals'
// include { HAPLOTYPING                   } from '../subworkflows/local/haplotyping'
// include { ONTARGET                      } from '../subworkflows/local/ontarget'
// include { ANNOTATION                    } from '../subworkflows/local/annotations'
// include { PHARMCAT                      } from '../subworkflows/local/pharmcat'
// include { COVERAGE                      } from '../subworkflows/local/coverage'
// include { CLINICAL_INFORMATION          } from '../subworkflows/local/clinical_information'
// include { PGX_REPORT                    } from '../subworkflows/local/reports'
// include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'
// include { TEST_FASTQ   } from '../modules/local/test/main'
// include { TEST_BAM   } from '../modules/local/test/main'
include { RUN_MULTIQC                   } from '../subworkflows/local/multiqc'

csv = file(params.csv)

workflow PGX_PANEL {

    take:
        fastq_input                 // FASTQ CHANNEL for PGX_FULL
        bam_input
        samples


    main:
        ch_versions = Channel.empty()

        // TEST_FASTQ ( fastq_input )
        // TEST_BAM ( bam_input )

        // Preprocessing Fastq (subsample and trim)
        PREPROCESS_FASTQ ( fastq_input ).set { fastq_processed }
        ch_versions = ch_versions.mix(fastq_processed.versions)

        // Alignment with Sentieon
        ALIGN ( fastq_processed.proccessed_fastq ).set { bam_aligned }
        ch_versions = ch_versions.mix(bam_aligned.versions)

        // Basically it will either have input from csv or from align sub workflow
        bam_input.mix(bam_aligned.bam_dedup).set{bam}

        // Variant calling
        VARIANT_CALLING ( bam )
        ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

        // Varinat_Filtering
        VARIANT_FILTRATION ( VARIANT_CALLING.out.aggregate_vcf )

        // Coverage
        COVERAGE ( bam )
        ch_versions = ch_versions.mix(COVERAGE.out.versions)

        // Pharmcat
        PHARMCAT ( VARIANT_CALLING.out.aggregate_vcf_tbi, COVERAGE.out.pc_panel_depth )
        ch_versions = ch_versions.mix(PHARMCAT.out.versions)

        // CNV calling CYP2D6
        CYP2D6_CNVCALL ( COVERAGE.out.cnv_depth )
        ch_versions = ch_versions.mix(CYP2D6_CNVCALL.out.versions)

    
        // Multiqc

    // emit:
    //     fastq_trimmed   =   TEST_FASTQ.out
    //     bam_trimmed     =   TEST_BAM.out
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