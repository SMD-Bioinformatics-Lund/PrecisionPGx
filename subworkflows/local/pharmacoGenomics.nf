#!/usr/bin/env nextflow

include { ONTARGET_BAM                          } from '../../modules/local/ontarget/main'
include { ONTARGET_VCF                          } from '../../modules/local/ontarget/main'
include { GATK_HAPLOTYPING                      } from '../../modules/local/haplotyping/main'
include { SENTIEON_HAPLOTYPING                  } from '../../modules/local/haplotyping/main'
include { BCFTOOLS_ANNOTATION                   } from '../../modules/local/annotation/main'
include { VARIANT_FILTRATION                    } from '../../modules/local/filtration/main'
include { PHARMCAT_PREPROCESSING                } from '../../modules/local/pharmcat/main'
include { PHARMCAT                              } from '../../modules/local/pharmcat/main'
include { DETECTED_VARIANTS                     } from '../../modules/local/variant_detection/main'
include { SAMPLE_TARGET_LIST                    } from '../../modules/local/summary/main'
include { DEPTH_OF_TARGETS                      } from '../../modules/local/summary/main'
include { GET_PADDED_BAITS                      } from '../../modules/local/summary/main'
include { DEPTH_OF_BAITS                        } from '../../modules/local/summary/main'
include { APPEND_ID_TO_GDF                      } from '../../modules/local/summary/main'
include { PADDED_BED_INTERVALS                  } from '../../modules/local/summary/main'
include { GET_CLIINICAL_GUIDELINES              } from '../../modules/local/pgx_report/main'
include { GET_INTERACTION_GUIDELINES            } from '../../modules/local/pgx_report/main'
include { GET_PGX_REPORT                        } from '../../modules/local/pgx_report/main'


workflow PHARMACO_GENOMICS {

    take:
        bam_input             // channel: [ tuple val(meta) file(".bam") file(".bai") ]

    main:
        ch_versions = Channel.empty()

        // Get the padded bed intervals
        PADDED_BED_INTERVALS ()
        ch_versions = ch_versions.mix(PADDED_BED_INTERVALS.out.versions)

        // Get the padded baits
        GET_PADDED_BAITS ()
        ch_versions = ch_versions.mix(GET_PADDED_BAITS.out.versions)

        // Runs only when the haplotype caller is GATK
        GATK_HAPLOTYPING ( bam_input )

        // Runs only when the haplotype caller is SENTIEON
        SENTIEON_HAPLOTYPING ( bam_input )

        ch_haplotypes = Channel.empty().mix(GATK_HAPLOTYPING.out.haplotypes, SENTIEON_HAPLOTYPING.out.haplotypes)
        ch_versions   = Channel.empty().mix(GATK_HAPLOTYPING.out.versions, SENTIEON_HAPLOTYPING.out.versions)

        // Filter the haplotypes
        VARIANT_FILTRATION ( ch_haplotypes )        
        ch_versions = ch_versions.mix(VARIANT_FILTRATION.out.versions)

        // Preprocess the pharmcat
        PHARMCAT_PREPROCESSING ( VARIANT_FILTRATION.out.haplotypes_filtered )
        ch_versions = ch_versions.mix(PHARMCAT_PREPROCESSING.out.versions)

        // Run the pharmcat
        PHARMCAT ( PHARMCAT_PREPROCESSING.out.pharmcat_preprocessed_vcf )
        ch_versions = ch_versions.mix(PHARMCAT.out.versions)

        // ONtarget Variants
        ONTARGET_VCF ( VARIANT_FILTRATION.out.haplotypes_filtered, PADDED_BED_INTERVALS.out.padded_bed_intervals )
        ch_versions = ch_versions.mix(ONTARGET_VCF.out.versions)

        // Annotate the haplotypes
        BCFTOOLS_ANNOTATION ( ONTARGET_VCF.out.vcf_ontarget )
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATION.out.versions)

        // Detect the variants from the target list
        DETECTED_VARIANTS ( BCFTOOLS_ANNOTATION.out.annotations )
        ch_versions = ch_versions.mix(DETECTED_VARIANTS.out.versions)

        // Get the sample target list
        SAMPLE_TARGET_LIST ( DETECTED_VARIANTS.out.detected_tsv )
        ch_versions = ch_versions.mix(SAMPLE_TARGET_LIST.out.versions)

        // Get the on-target bam file
        ONTARGET_BAM ( 
            bam_input, 
            PADDED_BED_INTERVALS.out.padded_bed_intervals 
        )
        ch_versions = ch_versions.mix(ONTARGET_BAM.out.versions)

        // Get the depth of targets
        DEPTH_OF_TARGETS ( 
            ONTARGET_BAM.out.bam_ontarget
            .join( SAMPLE_TARGET_LIST.out.pgx_target_interval_list, by: [0,1] ) 
        )
        ch_versions = ch_versions.mix(DEPTH_OF_TARGETS.out.versions)

        // Get the depth of baits
        DEPTH_OF_BAITS ( 
            ONTARGET_BAM.out.bam_ontarget, 
            GET_PADDED_BAITS.out.padded_baits_list 
        )
        ch_versions = ch_versions.mix(DEPTH_OF_BAITS.out.versions)

        // Append the id to the gdf
        APPEND_ID_TO_GDF ( DEPTH_OF_TARGETS.out.pgx_depth_at_missing )
        ch_versions = ch_versions.mix(APPEND_ID_TO_GDF.out.versions)

        // Get the clinical guidelines
        GET_CLIINICAL_GUIDELINES ( DETECTED_VARIANTS.out.detected_tsv )
        ch_versions = ch_versions.mix(GET_CLIINICAL_GUIDELINES.out.versions)

        // Get the interaction guidelines
        GET_INTERACTION_GUIDELINES ( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes )
        ch_versions = ch_versions.mix(GET_INTERACTION_GUIDELINES.out.versions)

        // Get the PGx report
        GET_PGX_REPORT ( 
            BCFTOOLS_ANNOTATION.out.annotations
            .join( DETECTED_VARIANTS.out.detected_tsv, by: [0,1])
            .join( APPEND_ID_TO_GDF.out.depth_at_missing_annotate_gdf, by: [0,1] )
            .join( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes, by: [0,1] )
            .join( DEPTH_OF_BAITS.out.padded_baits_list, by: [0,1] )
            .join( GET_INTERACTION_GUIDELINES.out.possible_interactions, by: [0,1] )
        )
        ch_versions = ch_versions.mix(GET_PGX_REPORT.out.versions)


    emit:
        pgx_report              = GET_PGX_REPORT.out.pgx_html                           // channel: [ tuple val(group), val(meta) file("pgx.html") ]
        // pharmcat_pgx_regions    = PHARMCAT_PREPROCESSING.out.pharmcat_pgx_regions    // channel: [ tuple val(group), val(meta) file("pharmcat_pgx_regions.vcf") ]
        pharmcat_preprocessed   = PHARMCAT_PREPROCESSING.out.pharmcat_preprocessed_vcf  // channel: [ tuple val(group), val(meta) file("pharmcat_preprocessed.vcf") ]
        pharmcat_report         = PHARMCAT.out.pharmcat_report                          // channel: [ tuple val(group), val(meta) file("pharmcat.html") ]
        pharmcat_pheno_json     = PHARMCAT.out.pharmcat_pheno_json                      // channel: [ tuple val(group), val(meta) file("pharmcat.phenotype.json") ]
        pharmcat_match_json     = PHARMCAT.out.pharmcat_match_json                      // channel: [ tuple val(group), val(meta) file("pharmcat.match.json") ]
        pharmcat_match_html     = PHARMCAT.out.pharmcat_match_html                      // channel: [ tuple val(group), val(meta) file("pharmcat.match.html") ]
        pharmcat_report_json     = PHARMCAT.out.pharmcat_report_json                    // channel: [ tuple val(group), val(meta) file("pharmcat.report.json") ]
        targets_depth           = GET_PGX_REPORT.out.targets_depth                      // channel: [ tuple val(group), val(meta) file(".targets.depth.tsv") ]
        versions                = ch_versions                                           // channel: [ path(versions.yml) ]
}