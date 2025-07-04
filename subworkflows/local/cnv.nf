#!/usr/bin/env nextflow

include { CALCULATE_CNV_CYP2D6                } from '../../modules/local/cnv_calling/main'
include { REPORT_CNV_CYP2D6                   } from '../../modules/local/cnv_calling/main'


workflow CYP2D6_CNVCALL {

    take:
        depth_file             // channel: [ tuple val(group), val(meta) file("haplotypes_filtered.vcf.gz") file("haplotypes_filtered.vcf.gz.tbi") ]

    main:
        ch_versions = Channel.empty()

        // CALCULATE THE CNV for cyped26 regions
        CALCULATE_CNV_CYP2D6 ( depth_file )
        ch_versions = ch_versions.mix(CALCULATE_CNV_CYP2D6.out.versions)

        // Annotate and filter CNV for CYPD26 regions
        REPORT_CNV_CYP2D6 ( CALCULATE_CNV_CYP2D6.out.cypd26_cnv )
        ch_versions = ch_versions.mix( REPORT_CNV_CYP2D6.out.versions)

    emit:
        reports   = REPORT_CNV_CYP2D6.out.cnvreport_cypd26      // channel: [ tuple val(group), val(meta), file("*.report_filtered.txt"), file ("*.CYP2D6_vs_CYP2D7_table.txt"), file("*.png") ]
        versions                = ch_versions                  // channel: [ path(versions.yml) ]
}