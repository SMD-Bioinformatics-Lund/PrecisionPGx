#!/usr/bin/env nextflow


include { SAMPLE_TARGET_LIST                    } from '../../modules/local/summary/main'
include { GET_CLIINICAL_GUIDELINES              } from '../../modules/local/clinical_information/main'
include { GET_INTERACTION_GUIDELINES            } from '../../modules/local/clinical_information/main'

workflow CLINICAL_INFORMATION {

    take:
        detected_variants             // channel: [ tuple val(meta) file(".bam") file(".bai") ]

    main:
        ch_versions = Channel.empty()


        // Get the sample target list
        // SAMPLE_TARGET_LIST ( detected_variants )
        // ch_versions = ch_versions.mix(SAMPLE_TARGET_LIST.out.versions)

        // Get the clinical guidelines
        GET_CLIINICAL_GUIDELINES ( detected_variants )
        ch_versions = ch_versions.mix(GET_CLIINICAL_GUIDELINES.out.versions)

        // Get the interaction guidelines
        GET_INTERACTION_GUIDELINES ( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes )
        ch_versions = ch_versions.mix(GET_INTERACTION_GUIDELINES.out.versions)

    emit:
        diplotypes              = GET_CLIINICAL_GUIDELINES.out.possible_diplotypes      // channel: [ tuple val(group), val(meta) file("possible_diplotypes.tsv") ]
        interactions            = GET_INTERACTION_GUIDELINES.out.possible_interactions  // channel: [ tuple val(group), val(meta) file("possible_interactions.tsv") ]
        versions                = ch_versions                                           // channel: [ path(versions.yml) ]
}