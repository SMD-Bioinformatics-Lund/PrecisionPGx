#!/usr/bin/env nextflow

include { GET_PGX_REPORT                        } from '../../modules/local/pgx_report/main'


workflow PGX_REPORT {

    take:
        annotated_vcf                       // channel: [ tuple val(meta) file(".bam") file(".bai") ]
        detected_variants                   // channel: [ tuple val(meta) file(".bam") file(".bai") ]
        depth_at_missing_annotate_gdf       // channel: [ tuple val(meta) file(".bam") file(".bai") ]
        possible_diplotypes                 // channel: [ tuple val(meta) file(".bam") file(".bai") ]
        baits_coverage                      // channel: [ tuple val(meta) file(".bam") file(".bai") ]
        possible_interactions               // channel: [ tuple val(meta) file(".bam") file(".bai") ]


    main:
        ch_versions = Channel.empty()


        // Get the PGx report
        GET_PGX_REPORT ( 
            annotated_vcf
            .join( detected_variants, by: [0,1])
            .join( depth_at_missing_annotate_gdf, by: [0,1] )
            .join( possible_diplotypes, by: [0,1] )
            .join( baits_coverage, by: [0,1] )
            .join( possible_interactions, by: [0,1] )
        )
        ch_versions = ch_versions.mix(GET_PGX_REPORT.out.versions)


    emit:
        pgx_report              = GET_PGX_REPORT.out.pgx_html               // channel: [ tuple val(group), val(meta) file("pgx.html") ]
        targets_depth           = GET_PGX_REPORT.out.targets_depth          // channel: [ tuple val(group), val(meta) file("targets_depth.tsv") ]
        versions                = ch_versions                               // channel: [ path(versions.yml) ]
}