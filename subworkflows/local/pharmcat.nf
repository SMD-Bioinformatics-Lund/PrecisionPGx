#!/usr/bin/env nextflow

// include { PHARMCAT_VCF_BED                      } from '../../modules/local/pharmcat/main'
include { PHARMCAT_PREPROCESSING                } from '../../modules/local/pharmcat/main'
include { PHARMCAT_RUN                          } from '../../modules/local/pharmcat/main'
include { PHARMCAT_RUN as PHARMCAT_RUN_CYP2D6   } from '../../modules/local/pharmcat/main'
include { EXTRACT_CYP2D6_HAPLOTYPES             } from '../../modules/local/summary/main'
include { UNCOVERED_PHARMCAT_POSTIONS           } from '../../modules/local/coverage/main'

workflow PHARMCAT {

    take:
        filtered_haplotypes             // channel: [ tuple val(group), val(meta) file("haplotypes_filtered.vcf.gz") file("haplotypes_filtered.vcf.gz.tbi") ]
        pc_pos_depth                    // channel: [ tuple val(group), val(meta) file("pharmcat_vcf.depth") ]

    main:
        ch_versions = Channel.empty()

        // Convert pharmcat vcf to bed
        // PHARMCAT_VCF_BED ()
        // ch_versions = ch_versions.mix(PHARMCAT_VCF_BED.out.versions)

        // Preprocess the pharmcat
        PHARMCAT_PREPROCESSING ( filtered_haplotypes )
        ch_versions = ch_versions.mix(PHARMCAT_PREPROCESSING.out.versions)

        // Run the pharmcat
        PHARMCAT_RUN ( PHARMCAT_PREPROCESSING.out.pharmcat_preprocessed_vcf )
        ch_versions = ch_versions.mix(PHARMCAT_RUN.out.versions)

        // uncovered pharmcat positions
        UNCOVERED_PHARMCAT_POSTIONS ( pc_pos_depth )
        ch_versions = ch_versions.mix(UNCOVERED_PHARMCAT_POSTIONS.out.versions)

        // Run the pharmcat
        PHARMCAT_RUN_CYP2D6 ( PHARMCAT_PREPROCESSING.out.pharmcat_preprocessed_vcf )
        ch_versions = ch_versions.mix(PHARMCAT_RUN.out.versions)

        // Extract CYP2D6 haplotypes
        EXTRACT_CYP2D6_HAPLOTYPES ( PHARMCAT_RUN_CYP2D6.out.pharmcat_match_json )


    emit:
        pc_pos_uncovered            = UNCOVERED_PHARMCAT_POSTIONS.out.lowcovered_pc_positions       // channel: [ tuple val(group), val(meta) file("*.pharmcat.pos.low.coverage.txt") ]
        // pharmcat_positions_bed      = PHARMCAT_VCF_BED.out.pharmcat_positions_bed                // channel: file("pharmcat_positions.bed") ]
        pharmcat_preprocessed       = PHARMCAT_PREPROCESSING.out.pharmcat_preprocessed_vcf          // channel: [ tuple val(group), val(meta) file("pharmcat_preprocessed.vcf") ]
        pharmcat_report             = PHARMCAT_RUN.out.pharmcat_report                              // channel: [ tuple val(group), val(meta) file("pharmcat.html") ]
        pharmcat_pheno_json         = PHARMCAT_RUN.out.pharmcat_pheno_json                          // channel: [ tuple val(group), val(meta) file("pharmcat.phenotype.json") ]
        pharmcat_match_json         = PHARMCAT_RUN.out.pharmcat_match_json                          // channel: [ tuple val(group), val(meta) file("pharmcat.match.json") ]
        // pharmcat_match_html         = PHARMCAT_RUN.out.pharmcat_match_html                       // channel: [ tuple val(group), val(meta) file("pharmcat.match.html") ]
        pharmcat_report_json        = PHARMCAT_RUN.out.pharmcat_report_json                         // channel: [ tuple val(group), val(meta) file("pharmcat.report.json") ]
        cy2d6_pharmcat_pheno_json   = PHARMCAT_RUN_CYP2D6.out.pharmcat_pheno_json                   // channel: [ tuple val(group), val(meta) file("pharmcat.phenotype.json") ]
        cy2d6_pharmcat_match_json   = PHARMCAT_RUN_CYP2D6.out.pharmcat_match_json                   // channel: [ tuple val(group), val(meta) file("pharmcat.match.json") ]
        cy2d6_pharmcat_report_json  = PHARMCAT_RUN_CYP2D6.out.pharmcat_report_json                  // channel: [ tuple val(group), val(meta) file("pharmcat.report.json") ]
        cyp2d6_haplotypes           = EXTRACT_CYP2D6_HAPLOTYPES.out.cy2d6_genotypes                 // channel: [ tuple val(group), val(meta) file("pharmcat.report.json") ]
        versions                    = ch_versions                                                   // channel: [ path(versions.yml) ]
}