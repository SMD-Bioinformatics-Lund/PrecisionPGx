#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_INPUT       } from './subworkflows/local/create_meta'
// include { PGX_FULL          } from './workflows/pgx'
// include { PGX_SHORT         } from './workflows/pgx_short'
include { PGX_PANEL         } from './workflows/pgx_panel'

println(params.genome_file)
csv = file(params.csv)
println(csv)



workflow {

    // SUBWORKFLOW: CHECK_INPUT
    CHECK_INPUT ( Channel.fromPath(csv) )
    samples = CHECK_INPUT.out.meta.collect { it[0] }.map { concatenate_values(it) }

    // FASTQ CHANNEL for PGX_FULL
    fastq_ch = CHECK_INPUT.out.bam_fq_ch
        .map { item -> item[2] ? item[0..3] : Channel.empty([[],[],[],[]]) }
    
    fastq_ch.view()
    // BAM CHANNEL for PGX_SHORT
    bam_ch = CHECK_INPUT.out.bam_fq_ch
        .map { item -> item[4] ? item[0,1,4,5] : Channel.empty([[],[],[],[]]) }

    bam_ch.view()

    // WORKFLOW: PGX_FULL
    // PGX_FULL(fastq_ch, samples)

    // WORKFLOW: PGX_SHORT
    // PGX_SHORT(bam_ch, samples)

    // WORKFLOW: PGX_PANEL
    PGX_PANEL(fastq_ch, bam_ch, samples)
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def concatenate_values(channel) {
    channel.collect().join('_')
}