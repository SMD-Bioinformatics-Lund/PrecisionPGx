include { SAMTOOLS_DEPTH                            } from '../../modules/nf-core/samtools/depth'
include { TARGET_PASS_REGIONS                       } from '../../modules/local/target_pass_regions'

workflow TARGET_DEPTH {

    take:
        ch_bam_bai               // channel: [ val(meta), path(bam), path(bai) ]
        ch_target_bed            // channel: [ val(meta), path(target_bed) ]

    main:
        // Samtools depth
        SAMTOOLS_DEPTH (
            ch_bam_bai,
            ch_target_bed
        )

        ch_tsv = SAMTOOLS_DEPTH.out.tsv

        // Create target pass regions bed
        TARGET_PASS_REGIONS (
            ch_tsv
        )

    emit:
        target_depth_tsv                        = SAMTOOLS_DEPTH.out.tsv                    // channel: [ val(meta), path(tsv) ]
        target_pass_bed                         = TARGET_PASS_REGIONS.out.pass_bed          // channel: [ val(meta), path(bed) ] 
        versions_samtools                       = SAMTOOLS_DEPTH.out.versions_samtools      // channel: [ val(process), val(samtools), val(version) ]
        versions                                = TARGET_PASS_REGIONS.out.versions          // channel: [ path(versions.yml) ]
}
