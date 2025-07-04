#!/usr/bin/env nextflow


include { GET_PADDED_BAITS                      } from '../../modules/local/summary/main'
include { PADDED_BED_INTERVALS                  } from '../../modules/local/summary/main'

workflow PADDED_INTERVALS {

    main:
        ch_versions = Channel.empty()

        // Get the padded bed intervals
        PADDED_BED_INTERVALS ()
        ch_versions = ch_versions.mix(PADDED_BED_INTERVALS.out.versions)

        // Get the padded baits
        GET_PADDED_BAITS ()
        ch_versions = ch_versions.mix(GET_PADDED_BAITS.out.versions)


    emit:
        padded_intervals_list   = PADDED_BED_INTERVALS.out.padded_baits_list             // channel: [ tuple val(group), val(meta) file("padded_bait_intervals.list") ]
        padded_intervals_bed    = GET_PADDED_BAITS.out.padded_bed_intervals              // channel: [ tuple val(group), val(meta) file("padded_bait_intervals.bed") ]
        versions                = ch_versions                                            // channel: [ path(versions.yml) ]
}