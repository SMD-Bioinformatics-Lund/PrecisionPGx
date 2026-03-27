include { PHARMCAT_MATCHER          } from '../../modules/nf-core/pharmcat/matcher/main'
include { PHARMCAT_PHENOTYPER       } from '../../modules/nf-core/pharmcat/phenotyper/main'
include { PHARMCAT_REPORTER         } from '../../modules/nf-core/pharmcat/reporter/main'

workflow PHARMCAT_GENOTYPING_REPORTING {

    take:
    ch_preprocessed_vcf_pass    // channel: [ val(meta), path(vcf), path(tbi) ]

    main:

    // Pharmcat Allele Matching
    PHARMCAT_MATCHER(
        ch_preprocessed_vcf_pass,
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
    matcher_json        = ch_pc_matches.matcher_json                        // channel: [ val(meta), [ match.json ] ]
    matcher_html        = ch_pc_matches.matcher_html                        // channel: [ val(meta), [ match.html ] ]
    phenotyper_json     = ch_phenotypes.phenotyper_json                     // channel: [ val(meta), [ phenotype.json ] ]
    report_json         = ch_pc_reports.report_json                         // channel: [ val(meta), [ report.json ] ]
    report_html         = ch_pc_reports.report_html                         // channel: [ val(meta), [ report.html ] ]
    report_tsv          = ch_pc_reports.report_tsv                          // channel: [ val(meta), [ report.tsv ] ]
}
