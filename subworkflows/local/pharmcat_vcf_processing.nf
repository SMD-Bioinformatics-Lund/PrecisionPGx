include { PHARMCAT_VCFPREPROCESSOR  } from '../../modules/nf-core/pharmcat/vcfpreprocessor/main'
include { TABIX_TABIX               } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_VIEW             } from '../../modules/nf-core/bcftools/view/main'

workflow PHARMCAT_VCF_PROCESSING {

    take:
    ch_vcf                      // channel: [ val(meta), [ vcf ], [ tbi ] ]
    ch_ref_fasta                // channel: [ val(meta), path(fasta) ]
    ch_ref_fasta_index          // channel: [ val(meta), path(fasta_index) ]
    ch_pc_pos                   // channel: [ val(meta), path(pharmcat_positions_vcf), path(pharmcat_positions_vcf_index) ]
    ch_pc_uniallelic_pos        // channel: [ val(meta), path(pharmcat_uniallelic_pos_vcf), path(pharmcat_uniallelic_pos_vcf_index) ]
    ch_target_pass_bed          // channel: [ val(meta), path(target_pass_bed) ]

    main:

    // VCF Preprocessing
    PHARMCAT_VCFPREPROCESSOR(
        ch_vcf, 
        ch_ref_fasta, 
        ch_ref_fasta_index, 
        ch_pc_pos.first(), 
        ch_pc_uniallelic_pos.first()
        ).set { ch_preprocessed_vcf } 

    // VCF indexing
    TABIX_TABIX(ch_preprocessed_vcf.preprocessed_vcf).set { ch_preprocessed_vcf_tbi }

    // Filter VCF to only include sites with depth > input min_dp.
    //
    // The VCF+tbi tuple and the per-sample target-pass BED come from
    // different upstream subworkflows (PHARMCAT_VCFPREPROCESSOR vs QC_BAM)
    // and carry meta dicts that share `id` but differ in other fields. A
    // `.join` on the full meta therefore fails to match (issue #27). Join on
    // `meta.id` only, then `multiMap` back into the two channels BCFTOOLS_VIEW
    // expects (VCF tuple + plain BED path).
    ch_preprocessed_vcf.preprocessed_vcf
        .join(ch_preprocessed_vcf_tbi.index, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .join(
            ch_target_pass_bed.map { meta, bed -> [meta.id, bed] },
            failOnMismatch: true, failOnDuplicate: true
        )
        .multiMap { _id, meta, vcf, tbi, bed ->
            vcf_tbi:  [meta, vcf, tbi]
            pass_bed: bed
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
