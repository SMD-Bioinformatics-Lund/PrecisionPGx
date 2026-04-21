//
// A quality check subworkflow for processed bams.
//

include { PICARD_COLLECTMULTIPLEMETRICS                            } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS                                  } from '../../modules/nf-core/picard/collecthsmetrics/main'
include { CHROMOGRAPH as CHROMOGRAPH_COV                           } from '../../modules/nf-core/chromograph/main'
include { TIDDIT_COV                                               } from '../../modules/nf-core/tiddit/cov/main'
include { MOSDEPTH                                                 } from '../../modules/nf-core/mosdepth/main'
include { VERIFYBAMID_VERIFYBAMID2                                 } from '../../modules/nf-core/verifybamid/verifybamid2/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_WG  } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { SENTIEON_WGSMETRICS as SENTIEON_WGSMETRICS_WG            } from '../../modules/nf-core/sentieon/wgsmetrics/main'

workflow QC_BAM {

    take:
        ch_bam_bai                  // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fasta             // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai               // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_dict              // channel: [mandatory] [ val(meta), path(dict) ]
        target_bed_uncompressed     // channel: [mandatory] [ path(target_bed) ]
        ch_bait_intervals           // channel: [mandatory] [ path(intervals_list) ]
        ch_target_intervals         // channel: [mandatory] [ path(intervals_list) ]
        ch_intervals_wgs            // channel: [mandatory] [ path(intervals) ]
        ch_svd_bed                  // channel: [optional] [ path(bed) ]
        ch_svd_mu                   // channel: [optional] [ path(meanpath) ]
        ch_svd_ud                   // channel: [optional] [ path(ud) ]

    main:
        ch_cov      = Channel.empty()
        ch_tiddit   = Channel.empty()
        ch_versions = Channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS (ch_bam_bai, ch_genome_fasta, ch_genome_fai)

        ch_bam_bai
            .combine(ch_bait_intervals)
            .combine(ch_target_intervals)
            .set { ch_hsmetrics_in}

        PICARD_COLLECTHSMETRICS (ch_hsmetrics_in, ch_genome_fasta, ch_genome_fai, ch_genome_dict, [[],[]])

        ch_bam_bai.combine(target_bed_uncompressed).set{ch_mosdepth_in}
        MOSDEPTH (ch_mosdepth_in, ch_genome_fasta)

        // COLLECT WGS METRICS
        if (!params.analysis_type.equals("panel")) {
            PICARD_COLLECTWGSMETRICS_WG ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs )
            SENTIEON_WGSMETRICS_WG ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs.map{ interval -> [[:], interval]} )
            ch_cov   = Channel.empty().mix(PICARD_COLLECTWGSMETRICS_WG.out.metrics, SENTIEON_WGSMETRICS_WG.out.wgs_metrics)

            TIDDIT_COV (ch_bam_bai, [[],[]]).set { tiddit_cov } // 2nd pos. arg is req. only for cram input
            ch_tiddit = Channel.empty().mix(tiddit_cov.wig)
            ch_versions = ch_versions.mix(tiddit_cov.versions.first())

            CHROMOGRAPH_COV([[:],[]], TIDDIT_COV.out.wig, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])
        }

        // Check contamination
        ch_svd_in = ch_svd_ud.combine(ch_svd_mu).combine(ch_svd_bed).collect()
        VERIFYBAMID_VERIFYBAMID2(ch_bam_bai, ch_svd_in, [], ch_genome_fasta.map {it-> it[1]})
        
    emit:
        multiple_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), path(metrics) ]
        tiddit_wig       = ch_tiddit                                 // channel: [ val(meta), path(wig) ]
        hs_metrics       = PICARD_COLLECTHSMETRICS.out.metrics       // channel: [ val(meta), path(metrics) ]
        d4               = MOSDEPTH.out.per_base_d4                  // channel: [ val(meta), path(d4) ]
        global_dist      = MOSDEPTH.out.global_txt                   // channel: [ val(meta), path(txt) ]
        self_sm          = VERIFYBAMID_VERIFYBAMID2.out.self_sm      // channel: [ val(meta), path(selfSM) ]
        cov              = ch_cov                                    // channel: [ val(meta), path(metrics) ]
        versions         = ch_versions                               // channel: [ path(versions.yml) ]
}