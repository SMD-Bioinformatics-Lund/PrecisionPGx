include { PHARMCAT_VCFPREPROCESSOR  } from '../../modules/nf-core/pharmcat/vcfpreprocessor/main'
include { PHARMCAT_MATCHER          } from '../../modules/nf-core/pharmcat/matcher/main'
include { PHARMCAT_PHENOTYPER       } from '../../modules/nf-core/pharmcat/phenotyper/main'
include { PHARMCAT_REPORTER         } from '../../modules/nf-core/pharmcat/reporter/main'


workflow PHARMCAT_PIPELINE {

    take:
    ch_vcf              // channel: [ val(meta), [ vcf ], [ tbi ] ]

    main:

    PHARMCAT_VCFPREPROCESSOR(ch_vcf)
    PHARMCAT_MATCHER(PHARMCAT_VCFPREPROCESSOR.out.vcf)
    PHARMCAT_PHENOTYPER(PHARMCAT_MATCHER.out.matches)
    PHARMCAT_REPORTER(PHARMCAT_PHENOTYPER.out.phenotypes)

    emit:
    preprocessed_vcf    = PHARMCAT_VCFPREPROCESSOR.out.preprocessed_vcf     // channel: [ val(meta), [ preprocessed.vcf.bgz ] ]
    missing_pgx_var_vcf = PHARMCAT_VCFPREPROCESSOR.out.missing_pgx_var      // channel: [ val(meta), [ missing_pgx_var.vcf ] ]
    matcher_json        = PHARMCAT_MATCHER.out.matcher_json                 // channel: [ val(meta), [ match.json ] ]
    matcher_html        = PHARMCAT_MATCHER.out.matcher_html                 // channel: [ val(meta), [ match.html ] ]
    phenotyper_json     = PHARMCAT_PHENOTYPER.out.phenotyper_json           // channel: [ val(meta), [ phenotype.json ] ]
    report_json         = PHARMCAT_REPORTER.out.report_json                 // channel: [ val(meta), [ report.json ] ]
    report_html         = PHARMCAT_REPORTER.out.report_html                 // channel: [ val(meta), [ report.html ] ]
    report_tsv          = PHARMCAT_REPORTER.out.report_tsv                  // channel: [ val(meta), [ report.tsv ] ]
}
