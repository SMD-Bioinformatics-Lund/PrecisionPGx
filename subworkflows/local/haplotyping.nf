#!/usr/bin/env nextflow

include { GATK_HAPLOTYPING                      } from '../../modules/local/haplotyping/main'
include { SENTIEON_HAPLOTYPING                  } from '../../modules/local/haplotyping/main'
include { VARIANT_FILTRATION                    } from '../../modules/local/filtration/main'
include { ONTARGET_BAM                          } from '../../modules/local/ontarget/main'
include { ONTARGET_VCF                          } from '../../modules/local/ontarget/main'
include { BCFTOOLS_ANNOTATION                   } from '../../modules/local/annotation/main'
include { DETECTED_VARIANTS                     } from '../../modules/local/variant_detection/main'


workflow HAPLOTYPING {

    take:
        bam_input             // channel: [ tuple val(meta) file(".bam") file(".bai") ]

    main:
        ch_versions = Channel.empty()

        // Runs only when the haplotype caller is GATK
        GATK_HAPLOTYPING ( bam_input )

        // Runs only when the haplotype caller is SENTIEON
        SENTIEON_HAPLOTYPING ( bam_input )

        ch_haplotypes = Channel.empty().mix(GATK_HAPLOTYPING.out.haplotypes, SENTIEON_HAPLOTYPING.out.haplotypes)
        ch_versions   = Channel.empty().mix(GATK_HAPLOTYPING.out.versions, SENTIEON_HAPLOTYPING.out.versions)

        // Filter the haplotypes
        VARIANT_FILTRATION ( ch_haplotypes ) .set { ch_filtered_vcf }        
        ch_versions = ch_versions.mix(VARIANT_FILTRATION.out.versions)

    emit:
        haplotypes                      = ch_haplotypes                                 // channel: [ tuple val(group), val(meta) file("haplotypes.vcf.gz") file("haplotypes.vcf.gz.tbi") ]
        filtered_haplotypes             = ch_filtered_vcf.haplotypes_filtered           // channel: [ tuple val(group), val(meta) file("haplotypes_filtered.vcf.gz") file("haplotypes_filtered.vcf.gz.tbi") ]
        versions                        = ch_versions                                   // channel: [ path(versions.yml) ]
}