include { PHARMCAT_MATCHER as PHARMCAT_CYP2D6_MATCHER                } from '../../modules/nf-core/pharmcat/matcher/main'
include { PHARMCAT_PHENOTYPER as PHARMCAT_CYP2D6_PHENOTYPER          } from '../../modules/nf-core/pharmcat/phenotyper/main'

workflow CYP2D6_CALLING {

    take:
    ch_preprocessed_vcf_pass    // channel: [ val(meta), path(vcf), path(tbi) ]

    main:
    // Pharmcat Allele Matching
    PHARMCAT_CYP2D6_MATCHER(
        ch_preprocessed_vcf_pass,
        ["CYP2D6"]
        ).set { ch_pc_cyp2d6_matches }

    // Pharmcat Phenotyping
    PHARMCAT_CYP2D6_PHENOTYPER(
        ch_pc_cyp2d6_matches.matcher_json.map {
            meta, matcher_json ->
            [
                meta,
                matcher_json,
                []  // TODO: Outside Pheno calls, will be updated when we add CYP2D6/HLA calls
            ]
        }
    ).set { ch_cyp2d6_phenotypes }

    emit:
    cyp2d6_matcher_json        = ch_pc_cyp2d6_matches.matcher_json                // channel: [ val(meta), [ match.json ] ]
    cyp2d6_matcher_html        = ch_pc_cyp2d6_matches.matcher_html                // channel: [ val(meta), [ match.html ] ]
    cyp2d6_phenotyper_json     = ch_cyp2d6_phenotypes.phenotyper_json       // channel: [ val(meta), [ phenotype.json ] ]
}
