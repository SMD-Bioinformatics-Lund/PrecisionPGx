#!/usr/bin/env nextflow

include { ONTARGET_BAM                          } from '../../modules/local/ontarget/main'
include { ONTARGET_VCF                          } from '../../modules/local/ontarget/main'

workflow ONTARGET {

    take:
        bam                     // channel: [ tuple val(group) val(meta) file(".bam") file(".bai") ]
        haplotypes_filtered     // channel: [ tuple val(group) val(meta) file(".vcf") ]

    main:
        ch_versions = Channel.empty()

        // Get the on-target bam file
        ONTARGET_BAM ( bam )
        ch_versions = ch_versions.mix(ONTARGET_BAM.out.versions)

        // ONtarget Variants
        ONTARGET_VCF ( haplotypes_filtered )
        ch_versions = ch_versions.mix(ONTARGET_VCF.out.versions)


    emit:
        ontarget_vcf                    = ONTARGET_VCF.out.vcf_ontarget                 // channel: [ tuple val(group), val(meta) file("ontarget.vcf") ]
        ontarget_bam                    = ONTARGET_BAM.out.bam_ontarget                 // channel: [ tuple val(group), val(meta) file("ontarget.bam") ]
        versions                        = ch_versions                                   // channel: [ path(versions.yml) ]
}