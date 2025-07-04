#!/usr/bin/env nextflow

include { SENTIEON_BWA_UMI          } from '../../modules/local/sentieon/main'
include { SENTIEON_MARKDUP          } from '../../modules/local/sentieon/main'
include { SENTIEON_BQSR_UMI         } from '../../modules/local/sentieon/main'
include { SENTIEON_QC               } from '../../modules/local/sentieon/main'


workflow ALIGN {
    take: 
        fastq_input         // channel: [mandatory] [ val(group), val(meta), read1, read2  ]

    main:
        ch_versions = Channel.empty()

        SENTIEON_BWA_UMI ( fastq_input )
        ch_versions = ch_versions.mix(SENTIEON_BWA_UMI.out.versions)

        SENTIEON_MARKDUP ( SENTIEON_BWA_UMI.out.bam_umi_markdup )
        ch_versions = ch_versions.mix(SENTIEON_MARKDUP.out.versions)

        SENTIEON_BQSR_UMI ( SENTIEON_BWA_UMI.out.bam_umi )
        ch_versions = ch_versions.mix(SENTIEON_BQSR_UMI.out.versions)

        SENTIEON_QC ( SENTIEON_MARKDUP.out.bam_qc )
        ch_versions = ch_versions.mix(SENTIEON_QC.out.versions)

    emit:
        bam_lowcov              =   SENTIEON_MARKDUP.out.bam_qc             // channel: [ val(group), val(meta), file(bam), file(bai), file(dedup_metrics.txt) ]
        bam_umi                 =   SENTIEON_BQSR_UMI.out.bam_varcall       // channel: [ val(group), val(meta), file(bam), file(bai), file(bqsr.table) ]
        qc_out                  =   SENTIEON_QC.out.qc_cdm                  // channel: [ val(group), val(meta), file(QC) ]
        dedup_bam_is_metrics    =   SENTIEON_QC.out.dedup_bam_is_metrics    // channel: [ val(group), val(meta), file(is_metrics.txt) ]    
        bam_dedup               =   SENTIEON_MARKDUP.out.bam_bqsr           // channel: [ val(group), val(meta), file(bam), file(bai)] 
        versions                =   ch_versions                             // channel: [ file(versions) ]
}