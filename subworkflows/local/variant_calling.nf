include { SENTIEON_HAPLOTYPER                       } from '../../modules/nf-core/sentieon/haplotyper'
include { GATK4_HAPLOTYPECALLER                     } from '../../modules/nf-core/gatk4/haplotypecaller'
include { DEEPVARIANT_RUNDEEPVARIANT                } from '../../modules/nf-core/deepvariant/rundeepvariant'



workflow VARIANT_CALLING {

    take:
        bam_bai_ch               // channel: [ val(meta), path(bam), path(bai) ]
        ch_ref_fasta             // channel: [ val(meta), path(fasta) ]
        ch_ref_fasta_index       // channel: [ val(meta), path(fasta_index) ]
        ch_ref_dict              // channel: [ val(meta), path(fasta_dict) ]
        ch_target_bed            // channel: [ val(meta), path(target_bed) ]
        ch_intervals             // channel: [ path(intervals) ]
        val_sentieon_emit_vcf    // boolean: [ optional ] applicable when the haplotyper is sentieon and emit mode is vcf
        val_sentieon_emit_gvcf   // boolean: [ optional ] applicable when the haplotyper is sentieon and emit mode is gvcf


    main:
        bam_bai_ch.combine(ch_intervals).set {ch_haplotyping_input }

        // Sentieon haplotyping
        ch_sentieon_hp_input = ch_haplotyping_input.map { meta, bam, bai, intervals ->
            return [ meta, bam, bai, intervals, [] ]
        }


        SENTIEON_HAPLOTYPER ( 
            ch_sentieon_hp_input, 
            ch_ref_fasta, 
            ch_ref_fasta_index,
            [[],[]],
            [[],[]],
            val_sentieon_emit_vcf,
            val_sentieon_emit_gvcf
        )

        if (val_sentieon_emit_vcf) {
            ch_sentieon_vcf = SENTIEON_HAPLOTYPER.out.vcf
            ch_sentieon_vcf_tbi = SENTIEON_HAPLOTYPER.out.vcf_tbi
        } else if (val_sentieon_emit_gvcf) {
            ch_sentieon_vcf = SENTIEON_HAPLOTYPER.out.gvcf
            ch_sentieon_vcf_tbi = SENTIEON_HAPLOTYPER.out.gvcf_tbi
        } else {
            ch_sentieon_vcf = channel.empty()
            ch_sentieon_vcf_tbi = channel.empty()
        }
        
        // GATK haplotyping
        ch_gatk_hp_input = ch_haplotyping_input.map { meta, bam, bai, intervals ->
                return [ meta, bam, bai, intervals, [] ]
            }
        GATK4_HAPLOTYPECALLER ( 
            ch_gatk_hp_input,
            ch_ref_fasta,
            ch_ref_fasta_index,
            ch_ref_dict,
            [[],[]],
            [[],[]],
        )

        // DEEP VARIANT
        // Provide the target bed to DEEPVARIANT
        bam_bai_ch.combine(ch_target_bed.map { meta, target_bed -> return target_bed}).set {ch_targeted_variant_calls_input}
        ch_deepvariant_input = ch_targeted_variant_calls_input.map { meta, bam, bai, target_bed ->
            return [ meta, bam, bai, target_bed ]
        }

        //ch_deepvariant_input = bam_bai_ch.map { meta, bam, bai -> return [ meta, bam, bai, [] ] }
        DEEPVARIANT_RUNDEEPVARIANT (
            ch_deepvariant_input,
            ch_ref_fasta,
            ch_ref_fasta_index,
            [[],[]],
            [[],[]],
        )

    emit:
        sentieon_vcf                            = ch_sentieon_vcf                           // channel: [ val(meta), path(vcf) ]
        sentieon_vcf_tbi                        = ch_sentieon_vcf_tbi                       // channel: [ val(meta), path(tbi) ]
        gatk_vcf                                = GATK4_HAPLOTYPECALLER.out.vcf             // channel: [ val(meta), path(vcf) ]
        gatk_vcf_tbi                            = GATK4_HAPLOTYPECALLER.out.tbi             // channel: [ val(meta), path(tbi) ]
        deepvariant_vcf                         = DEEPVARIANT_RUNDEEPVARIANT.out.vcf        // channel: [ val(meta), path(vcf) ]
        deepvariant_gvcf                        = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf       // channel: [ val(meta), path(vcf) ]
        deepvariant_vcf_tbi                     = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi    // channel: [ val(meta), path(tbi) ]
        deepvariant_gvcf_tbi                    = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_tbi   // channel: [ val(meta), path(tbi) ]

}
