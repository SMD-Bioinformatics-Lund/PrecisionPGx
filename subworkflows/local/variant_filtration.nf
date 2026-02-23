include { BCFTOOLS_NORM                 } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_MERGE                } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_VIEW                 } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_FILTER               } from '../../modules/nf-core/bcftools/filter/main'


workflow VARIANT_FILTRATION {

    take:
        vcf                         // channel: [ val(meta), path(vcf), path(vcf_tbi) ]
        ch_ref_fasta                // channel: [ val(meta), path(fasta) ]
        ch_ref_fasta_index          // channel: [ val(meta), path(fasta_index)
        ch_target_bed               // channel: [ path(target_bed) ]
        ch_pharmcat_positions_vcf   // channel: [ val(meta), path(pharmcat_positions_vcf) ]
    main:

        // Normalize the VCF using BCFTOOLS NORM
        BCFTOOLS_NORM ( vcf, ch_ref_fasta ).set({ ch_norm_vcf })

        // Merge the normalized VCF with the reference VCF using BCFTOOLS MERGE

        merge_input_vcfs = ch_norm_vcf.vcf.map{
            meta, vcf_path -> tuple(meta, [vcf_path, ch_pharmcat_positions_vcf] )
        }

        BCFTOOLS_MERGE ( 
            merge_input_vcfs.groupTuple(
                ch_norm_vcf.tbi
            ),
            ch_ref_fasta,
            ch_ref_fasta_index,
            ch_target_bed 
        ).set { ch_merged_vcf }

        // Filter the merged VCF using BCFTOOLS VIEW to retain only variants in the target regions
        BCFTOOLS_VIEW ( 
            ch_merged_vcf.vcf.groupTuple(
                ch_merged_vcf.index
            ),
            ch_target_bed,
            [],
            [] 
        ).set { ch_merged_view }

        // Apply additional filters to the merged VCF using BCFTOOLS FILTER
        BCFTOOLS_FILTER ( 
            ch_merged_view.vcf.groupTuple(
                ch_merged_view.tbi
            ) 
        ).set { ch_filtered_vcf }

    emit:
        merged_vcf      =   ch_merged_view.vcf.groupTuple(ch_merged_view.tbi)           // channel: [ val(meta), path(merged.vcf), path(merged.vcf.tbi) ]
        filtered_vcf    =   ch_filtered_vcf.vcf.groupTuple(ch_filtered_vcf.tbi)         // channel: [ val(meta), path(filtered.vcf), path(filtered.vcf.tbi) ]

}