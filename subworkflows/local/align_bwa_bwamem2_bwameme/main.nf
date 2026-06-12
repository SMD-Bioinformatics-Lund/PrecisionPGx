//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

include { BWA_MEM as BWA                           } from '../../../modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWAMEM_FALLBACK               } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM                              } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEME_MEM                              } from '../../../modules/nf-core/bwameme/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGN   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_EXTRACT } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                           } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_MERGE                           } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_VIEW as EXTRACT_ALIGNMENTS      } from '../../../modules/nf-core/samtools/view/main'
include { PICARD_MARKDUPLICATES as MARKDUPLICATES  } from '../../../modules/nf-core/picard/markduplicates/main'


workflow ALIGN_BWA_BWAMEM2_BWAMEME {
    take:
        ch_reads_input   // channel: [mandatory] [ val(meta), path(reads_input) ]
        ch_bwa_index     // channel: [mandatory] [ val(meta), path(bwamem2_index) ]
        ch_bwamem2_index // channel: [mandatory] [ val(meta), path(bwamem2_index) ]
        ch_bwameme_index // channel: [mandatory] [ val(meta), path(bwamem2_index) ]
        ch_genome_fasta  // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai    // channel: [mandatory] [ val(meta), path(fai) ]
        val_mbuffer_mem  // integer: [mandatory] default: 3072
        val_platform     // string:  [mandatory] default: illumina
        val_sort_threads // integer: [mandatory] default: 4
    main:

        // Map, sort, and index
        if (params.aligner.equals("bwa")) {
            BWA ( ch_reads_input, ch_bwa_index, ch_genome_fasta, true )
            ch_align = BWA.out.bam
        } else if (params.aligner.equals("bwameme")) {
            BWAMEME_MEM ( ch_reads_input, ch_bwameme_index, ch_genome_fasta, true, val_mbuffer_mem, val_sort_threads )
            ch_align = BWAMEME_MEM.out.bam
        } else {
            BWAMEM2_MEM ( ch_reads_input, ch_bwamem2_index, ch_genome_fasta, true )
            ch_align    = BWAMEM2_MEM.out.bam

            if (params.bwa_as_fallback) {
                ch_reads_input
                    .join(BWAMEM2_MEM.out.bam, remainder: true)
                    .branch { it ->
                        ERROR: it[2].equals(null)
                            return [it[0], it[1]] // return reads
                        SUCCESS: !it[2].equals(null)
                            return [it[0], it[2]]  // return bam
                    }
                    .set { ch_fallback }

                BWAMEM_FALLBACK ( ch_fallback.ERROR, ch_bwa_index, ch_genome_fasta, true )
                ch_align = ch_fallback.SUCCESS.mix(BWAMEM_FALLBACK.out.bam)
            }
        }

        SAMTOOLS_INDEX_ALIGN ( ch_align )

        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = ch_align.join(SAMTOOLS_INDEX_ALIGN.out.index, failOnMismatch:true, failOnDuplicate:true)
        SAMTOOLS_STATS (
            bam_sorted_indexed,
            ch_genome_fasta.combine(ch_genome_fai.map { _meta, fai -> fai })
                           .map { meta, fasta, fai -> [meta, fasta, fai] }
        )

        // Merge multiple lane samples and index
        ch_align
            .map{ meta, bam ->
                    def new_id   = meta.sample
                    def new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"] - meta.subMap('lane','data_type')
                    [groupKey(new_meta, new_meta.num_lanes), bam]
                }
            .groupTuple()
            .branch{
                single: it[1].size() == 1
                multiple: it[1].size() > 1
                }
            .set{ bams }

        // If there are no samples to merge, skip the process.
        // The nf-core SAMTOOLS_MERGE module bundles (bams, indices) into one
        // tuple and (fasta, fai, gzi) into another, so build matching shapes.
        SAMTOOLS_MERGE (
            bams.multiple.map { meta, bam_list -> [meta, bam_list, []] },
            ch_genome_fasta.combine(ch_genome_fai.map { _meta, fai -> fai })
                           .map { meta, fasta, fai -> [meta, fasta, fai, []] }
        )
        prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)

        // GET ALIGNMENT FROM SELECTED CONTIGS
        if (params.extract_alignments) {
            SAMTOOLS_INDEX_EXTRACT ( prepared_bam )
            extract_bam_sorted_indexed = prepared_bam.join(SAMTOOLS_INDEX_EXTRACT.out.index, failOnMismatch:true, failOnDuplicate:true)
            EXTRACT_ALIGNMENTS( extract_bam_sorted_indexed, ch_genome_fasta, [])
            prepared_bam = EXTRACT_ALIGNMENTS.out.bam
        }

        // Marking duplicates
        MARKDUPLICATES ( prepared_bam , ch_genome_fasta, ch_genome_fai )
        SAMTOOLS_INDEX_MARKDUP ( MARKDUPLICATES.out.bam )
    emit:
        stats       = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), path(stats) ]
        metrics     = MARKDUPLICATES.out.metrics     // channel: [ val(meta), path(metrics) ]
        marked_bam  = MARKDUPLICATES.out.bam         // channel: [ val(meta), path(bam) ]
        marked_bai  = SAMTOOLS_INDEX_MARKDUP.out.index // channel: [ val(meta), path(bai) ]
}
