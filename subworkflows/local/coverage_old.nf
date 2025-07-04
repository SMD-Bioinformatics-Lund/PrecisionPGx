#!/usr/bin/env nextflow

include { DEPTH_OF_TARGETS                      } from '../../modules/local/summary/main'
include { DEPTH_OF_BAITS                        } from '../../modules/local/summary/main'
include { APPEND_ID_TO_GDF                      } from '../../modules/local/summary/main'
include { DEPTH as SAMTOOLS_DEPTH               } from '../../modules/local/samtools/main'

workflow COVERAGE {

    take:
        bam                     // channel: [ tuple val(meta) file(".bam") file(".bai") ]

    main:
        ch_versions = Channel.empty()

        // Get the depth of targets
        DEPTH_OF_TARGETS ( bam ) 
        ch_versions = ch_versions.mix(DEPTH_OF_TARGETS.out.versions)

        // Get the depth of baits
        DEPTH_OF_BAITS ( bam )
        ch_versions = ch_versions.mix(DEPTH_OF_BAITS.out.versions)

        // Append the id to the gdf
        APPEND_ID_TO_GDF ( DEPTH_OF_TARGETS.out.pgx_depth_at_missing )
        ch_versions = ch_versions.mix(APPEND_ID_TO_GDF.out.versions)

        // samtools depth
        SAMTOOLS_DEPTH ( bam )
        ch_versions = ch_versions.mix(DEPTH.out.versions)


    emit:
        missing_targets_coverage               = DEPTH_OF_TARGETS.out.pgx_depth_at_missing             // channel: [ tuple val(group), val(meta) file("depth_at_missing.gdf") ]
        missing_targets_coverage_annotated     = APPEND_ID_TO_GDF.out.depth_at_missing_annotate_gdf    // channel: [ tuple val(group), val(meta) file("depth_at_missing_annotate.gdf") ]
        baits_coverage                         = DEPTH_OF_BAITS.out.padded_baits_list                  // channel: [ tuple val(group), val(meta) file("baits_depth.gdf") ]
        versions                               = ch_versions                                           // channel: [ path(versions.yml) ]
}