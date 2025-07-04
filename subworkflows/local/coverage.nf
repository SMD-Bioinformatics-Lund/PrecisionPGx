#!/usr/bin/env nextflow

include { DEPTH                                 } from '../../modules/local/samtools/main'
include { DEPTH as PHARMCAT_DEPTH               } from '../../modules/local/samtools/main'
include { DEPTH as CNV_DEPTH                    } from '../../modules/local/samtools/main'
include { COVERAGE_REPORTS                      } from '../../modules/local/coverage/main'
include { PROBLAMATIC_REGIONS_COVERAGE_REPORTS  } from '../../modules/local/coverage/main'

workflow COVERAGE {

    take:
        bam                     // channel: [ tuple val(meta) file(".bam") file(".bai") ]

    main:
        ch_versions = Channel.empty()

        // samtools depth
        DEPTH ( bam )
        ch_versions = ch_versions.mix(DEPTH.out.versions)

        // samtools depth for pharmcat
        PHARMCAT_DEPTH ( bam )
        ch_versions = ch_versions.mix(PHARMCAT_DEPTH.out.versions)

        // samtools depth for CNV
        CNV_DEPTH ( bam )
        ch_versions = ch_versions.mix(CNV_DEPTH.out.versions)

        // coverage reports
        COVERAGE_REPORTS ( DEPTH.out.depth )
        ch_versions = ch_versions.mix(COVERAGE_REPORTS.out.versions)

        // Coverage Problamatic Regions
        PROBLAMATIC_REGIONS_COVERAGE_REPORTS ( DEPTH.out.depth )
        ch_versions = ch_versions.mix(PROBLAMATIC_REGIONS_COVERAGE_REPORTS.out.versions)

    emit:
        panel_depth                 = DEPTH.out.depth
        cnv_depth                   = CNV_DEPTH.out.depth
        pc_panel_depth              = PHARMCAT_DEPTH.out.depth
        cov_stats                   = COVERAGE_REPORTS.out.cov_stats
        cov_annotated               = COVERAGE_REPORTS.out.cov_annotated
        cov_html                    = COVERAGE_REPORTS.out.cov_html
        cov_plots                   = COVERAGE_REPORTS.out.cov_plots
        overlap_cov_txt             = PROBLAMATIC_REGIONS_COVERAGE_REPORTS.out.overlap_cov_txt
        overlap_cov_html            = PROBLAMATIC_REGIONS_COVERAGE_REPORTS.out.overlap_cov_html
        overlap_cov_plots           = PROBLAMATIC_REGIONS_COVERAGE_REPORTS.out.overlap_cov_plots
        versions                    = ch_versions                                                       // channel: [ path(versions.yml) ]
}