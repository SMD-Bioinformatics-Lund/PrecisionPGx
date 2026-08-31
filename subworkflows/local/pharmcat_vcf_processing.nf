include { PHARMCAT_VCFPREPROCESSOR  } from '../../modules/nf-core/pharmcat/vcfpreprocessor/main'
include { TABIX_TABIX               } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_VIEW             } from '../../modules/nf-core/bcftools/view/main'

workflow PHARMCAT_VCF_PROCESSING {

    take:
    ch_vcf                      // channel: [ val(meta), [ vcf ], [ tbi ] ]
    ch_ref_fasta                // channel: [ val(meta), path(fasta) ]
    ch_ref_fasta_indices        // channel: [ val(meta), [ path(fai), path(gzi) ] ]
    ch_pc_pos                   // channel: [ val(meta), path(pharmcat_positions_vcf), path(pharmcat_positions_vcf_index) ]
    ch_pc_uniallelic_pos        // channel: [ val(meta), path(pharmcat_uniallelic_pos_vcf), path(pharmcat_uniallelic_pos_vcf_index) ]
    ch_target_pass_bed          // channel: [ val(meta), path(target_pass_bed) ]

    main:

    // VCF Preprocessing
    PHARMCAT_VCFPREPROCESSOR(
        ch_vcf,
        ch_ref_fasta,
        ch_ref_fasta_indices,
        ch_pc_pos.first(),
        ch_pc_uniallelic_pos.first()
        ).set { ch_preprocessed_vcf }

    // VCF indexing
    TABIX_TABIX(ch_preprocessed_vcf.preprocessed_vcf).set { ch_preprocessed_vcf_tbi }

    // Filter VCF to only include sites with depth > input min_dp.
    // Join VCF, index and target-pass BED on meta so each sample's BED stays paired
    // with its own VCF; then split into the two channels BCFTOOLS_VIEW expects.
    ch_preprocessed_vcf.preprocessed_vcf
        .join(ch_preprocessed_vcf_tbi.index, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_target_pass_bed,            failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, vcf, tbi, pass_bed ->
            vcf_tbi:  [ meta, vcf, tbi ]
            pass_bed: pass_bed
        }
        .set { ch_view_input }

    BCFTOOLS_VIEW (
            ch_view_input.vcf_tbi,
            ch_view_input.pass_bed,
            [],
            []
        ).set { ch_preprocessed_vcf_pass }

    emit:
    preprocessed_vcf    = ch_preprocessed_vcf.preprocessed_vcf    // channel: [ val(meta), [ preprocessed.vcf.bgz ] ]
    missing_pgx_var_vcf = ch_preprocessed_vcf.missing_pgx_var               // channel: [ val(meta), [ missing_pgx_var.vcf ] ]
    preprocessed_vcf_pass    = ch_preprocessed_vcf_pass.vcf       // channel: [ val(meta), path(vcf) ]
    preprocessed_vcf_pass_tbi    = ch_preprocessed_vcf_pass.tbi   // channel: [ val(meta), path(vcf.tbi) ]
}
