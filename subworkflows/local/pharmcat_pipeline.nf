include { PHARMCAT_VCFPREPROCESSOR  } from '../../modules/nf-core/pharmcat/vcfpreprocessor/main'
include { PHARMCAT_MATCHER          } from '../../modules/nf-core/pharmcat/matcher/main'
include { PHARMCAT_PHENOTYPER       } from '../../modules/nf-core/pharmcat/phenotyper/main'
include { PHARMCAT_REPORTER         } from '../../modules/nf-core/pharmcat/reporter/main'
include { TABIX_TABIX               } from '../../modules/nf-core/tabix/tabix/main'


workflow PHARMCAT_PIPELINE {

    take:
    ch_vcf                      // channel: [ val(meta), [ vcf ], [ tbi ] ]
    ch_ref_fasta                // channel: [ val(meta), path(fasta) ]
    ch_ref_fasta_index          // channel: [ val(meta), path(fasta_index) ]
    ch_pc_pos                   // channel: [ val(meta), path(pharmcat_positions_vcf), path(pharmcat_positions_vcf_index) ]
    ch_pc_uniallelic_pos        // channel: [ val(meta), path(pharmcat_uniallelic_pos_vcf), path(pharmcat_uniallelic_pos_vcf_index) ]

    main:

    // VCF Preprocessing
    PHARMCAT_VCFPREPROCESSOR(ch_vcf, ch_ref_fasta, ch_ref_fasta_index, ch_pc_pos, ch_pc_uniallelic_pos).set { ch_preprocessed_vcf } 
    TABIX_TABIX(ch_preprocessed_vcf.preprocessed_vcf).set { ch_preprocessed_vcf_tbi }

    // Pharmcat Alelle Matching
    PHARMCAT_MATCHER(
        ch_preprocessed_vcf.preprocessed_vcf.join(
            ch_preprocessed_vcf_tbi.index, 
            failOnMismatch:true, 
            failOnDuplicate:true
            ),
            [],
        ).set { ch_pc_matches }

    // Pharmcat Phenotyping
    PHARMCAT_PHENOTYPER(
        ch_pc_matches.matcher_json.map {
            meta, matcher_json -> 
            [
                meta,
                matcher_json,
                []  // TODO: Outside Pheno calls, will be updated when we add CYP2D6/HLA calls
            ]
        }
    ).set { ch_phenotypes }

    // Pharmcat Reporting
    PHARMCAT_REPORTER(ch_phenotypes.phenotyper_json).set { ch_pc_reports }

    emit:
    preprocessed_vcf    = ch_preprocessed_vcf.preprocessed_vcf              // channel: [ val(meta), [ preprocessed.vcf.bgz ] ]
    missing_pgx_var_vcf = ch_preprocessed_vcf.missing_pgx_var               // channel: [ val(meta), [ missing_pgx_var.vcf ] ]
    matcher_json        = ch_pc_matches.matcher_json                        // channel: [ val(meta), [ match.json ] ]
    matcher_html        = ch_pc_matches.matcher_html                        // channel: [ val(meta), [ match.html ] ]
    phenotyper_json     = ch_phenotypes.phenotyper_json                     // channel: [ val(meta), [ phenotype.json ] ]
    report_json         = ch_pc_reports.report_json                         // channel: [ val(meta), [ report.json ] ]
    report_html         = ch_pc_reports.report_html                         // channel: [ val(meta), [ report.html ] ]
    report_tsv          = ch_pc_reports.report_tsv                          // channel: [ val(meta), [ report.tsv ] ]
}
