#!/usr/bin/env nextflow

include { BCFTOOLS_ANNOTATION                   } from '../../modules/local/annotation/main'
include { DETECTED_VARIANTS                     } from '../../modules/local/variant_detection/main'


workflow ANNOTATION {

    take:
        haplotypes             // channel: [ tuple val(meta) file(".bam") file(".bai") ]

    main:
        ch_versions = Channel.empty()

        // Annotate the haplotypes
        BCFTOOLS_ANNOTATION ( haplotypes )
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATION.out.versions)

        // Detect the variants from the target list
        DETECTED_VARIANTS ( BCFTOOLS_ANNOTATION.out.annotations )
        ch_versions = ch_versions.mix(DETECTED_VARIANTS.out.versions)

    emit:
        annotated_haplotypes            = BCFTOOLS_ANNOTATION.out.annotations           // channel: [ tuple val(group), val(meta) file("annotated.vcf") ]
        detected_variants               = DETECTED_VARIANTS.out.detected_tsv            // channel: [ tuple val(group), val(meta) file("detected_variants.tsv") ]
        versions                        = ch_versions                                   // channel: [ path(versions.yml) ]
}