include { BCFTOOLS_NORM                 } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_MERGE                } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_VIEW                 } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_FILTER               } from '../../modules/nf-core/bcftools/filter/main'


workflow VARIANT_FILTRATION {

    take:
        vcf                             // channel: [ val(meta), path(vcf), path(vcf_tbi) ]
        ch_ref_fasta                    // channel: [ val(meta), path(fasta) ]
        ch_ref_fasta_index              // channel: [ val(meta), path(fasta_index)
        ch_target_bed                   // channel: [ val(meta), path(target_bed) ]
        ch_pc_positions_vcf             // channel: [ path(pharmcat_positions_vcf) ]
        ch_pc_positions_vcf_tbi         // channel: [ path(pharmcat_positions_vcf.tbi) ]
    main:

        // Normalize the VCF using BCFTOOLS NORM
        BCFTOOLS_NORM ( vcf, ch_ref_fasta ).set{ ch_norm_vcf }


        // Merge the normalized VCF with the PharmCAT positions VCF using BCFTOOLS MERGE, ensuring that both the VCF and its index (TBI) are provided for merging
        merge_vcfs_input = ch_norm_vcf.vcf
            .combine(ch_pc_positions_vcf)
            .map { meta, vcf_path , pc_vcf ->  [ meta, [vcf_path, pc_vcf] ] }

        merge_tbis_input = ch_norm_vcf.tbi
            .combine(ch_pc_positions_vcf_tbi)
            .map { meta, tbi_path, pc_tbi ->  [ meta, [tbi_path, pc_tbi] ] }


        BCFTOOLS_MERGE ( 
            merge_vcfs_input.join(
                merge_tbis_input,
                failOnMismatch:true, 
                failOnDuplicate:true
            ),
            ch_ref_fasta,
            ch_ref_fasta_index,
            ch_target_bed
        ).set { ch_merged_vcf }


        // Filter the merged VCF using BCFTOOLS VIEW to retain only variants in the target regions
        BCFTOOLS_VIEW ( 
            ch_merged_vcf.vcf.join(
                ch_merged_vcf.index,
                failOnMismatch:true, 
                failOnDuplicate:true
            ),
            ch_target_bed.map { meta, bed_path -> [ bed_path ] },
            [],
            [] 
        ).set { ch_merged_view }

        // Apply additional filters to the merged VCF using BCFTOOLS FILTER
        BCFTOOLS_FILTER ( 
            ch_merged_view.vcf.join(
                ch_merged_view.tbi,
                failOnMismatch:true, 
                failOnDuplicate:true
            ) 
        ).set { ch_filtered_vcf }

    emit:
        merged_vcf          = ch_merged_view.vcf            // channel: [ val(meta), path(merged.vcf), path(merged.vcf.tbi) ]
        merged_vcf_tbi      = ch_merged_view.tbi        // channel: [ val(meta), path(merged.vcf), path(merged.vcf.tbi) ]
        filtered_vcf        = ch_filtered_vcf.vcf           // channel: [ val(meta), path(filtered.vcf), path(filtered.vcf.tbi) ]
        filtered_vcf_tbi    = ch_filtered_vcf.tbi       // channel: [ val(meta), path(filtered.vcf), path(filtered.vcf.tbi) ]

}