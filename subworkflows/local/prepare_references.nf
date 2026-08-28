//
// Prepare reference files
//

include { BEDTOOLS_SLOP as BEDTOOLS_PAD_TARGET_BED           } from '../../modules/nf-core/bedtools/slop/main'
include { BWA_INDEX as BWA_INDEX_GENOME                      } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME              } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEME_INDEX as BWAMEME_INDEX_GENOME              } from '../../modules/nf-core/bwameme/index/main'
include { CAT_CAT as CAT_CAT_BAIT                            } from '../../modules/nf-core/cat/cat/main'
include { GATK4_BEDTOINTERVALLIST as GATK_BILT               } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD          } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD_MT       } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS as GATK_ILT                } from '../../modules/nf-core/gatk4/intervallisttools/main'
include { GET_CHROM_SIZES                                    } from '../../modules/local/get_chrom_sizes'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME            } from '../../modules/nf-core/samtools/faidx/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_GENOME      } from '../../modules/nf-core/sentieon/bwaindex/main'
include { TABIX_BGZIPTABIX as TABIX_PBT                      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_PT                            } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPINDEX_PADDED_BED    } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow PREPARE_REFERENCES {
    take:
        ch_genome_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai                // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_dictionary         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_target_bed                // channel: [mandatory for WES] [ path(bed) ]

    main:
        ch_tbi           = Channel.empty()
        ch_bgzip_tbi     = Channel.empty()
        ch_bwa           = Channel.empty()
        ch_sentieonbwa   = Channel.empty()

        // Genome indices
        // TODO: Update this  to get the sizes from the fai file.
        SAMTOOLS_FAIDX_GENOME(
            ch_genome_fasta.join(ch_genome_fai), 
            false
        ) 
        GATK_SD(ch_genome_fasta)
        ch_fai  = Channel.empty().mix(ch_genome_fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        ch_dict = Channel.empty().mix(ch_genome_dictionary, GATK_SD.out.dict).collect()
        GET_CHROM_SIZES( ch_fai )

        // Genome alignment indices
        BWA_INDEX_GENOME(ch_genome_fasta).index.set{ch_bwa}
        BWAMEM2_INDEX_GENOME(ch_genome_fasta)
        BWAMEME_INDEX_GENOME(ch_genome_fasta)
        SENTIEON_BWAINDEX_GENOME(ch_genome_fasta).index.set{ch_sentieonbwa}

        // Index target bed file in case of gz input
        TABIX_PT(ch_target_bed)
        ch_target_bed
            .join(TABIX_PT.out.index)
            .set{ ch_trgt_bed_tbi }
        // Compress and index target bed file in case of uncompressed input
        TABIX_PBT(ch_target_bed).gz_index
            .set { ch_bgzip_tbi }
        ch_target_bed_gz_tbi = Channel.empty()
            .mix(ch_trgt_bed_tbi, ch_bgzip_tbi)

        // Pad bed file
        BEDTOOLS_PAD_TARGET_BED(
            ch_target_bed,
            ch_fai.map { _meta, fai -> return fai }
        )
        TABIX_BGZIPINDEX_PADDED_BED(BEDTOOLS_PAD_TARGET_BED.out.bed).gz_index
            .set { ch_target_bed_gz_tbi }

        // Generate bait and target intervals
        GATK_BILT(ch_target_bed, ch_dict).interval_list
        GATK_ILT(GATK_BILT.out.interval_list)
        GATK_ILT.out.interval_list
            .collect{ it[1] }
            .map { it ->
                def meta = it[0].toString().split("_split")[0].split("/")[-1] + "_bait.intervals_list"
                return [[id:meta], it]
            }
            .set { ch_bait_intervals_cat_in }
        CAT_CAT_BAIT ( ch_bait_intervals_cat_in )

    emit:
        genome_bwa_index        = Channel.empty().mix(ch_bwa, ch_sentieonbwa).collect()                                 // channel: [ val(meta), path(index) ]
        genome_bwamem2_index    = BWAMEM2_INDEX_GENOME.out.index.collect()                                              // channel: [ val(meta), path(index) ]
        genome_bwameme_index    = BWAMEME_INDEX_GENOME.out.index.collect()                                              // channel: [ val(meta), path(index) ]
        genome_chrom_sizes      = GET_CHROM_SIZES.out.sizes.collect()                                                   // channel: [ path(sizes) ]
        genome_fai              = ch_fai                                                                                // channel: [ val(meta), path(fai) ]
        genome_dict             = ch_dict                                                                               // channel: [ val(meta), path(dict) ]
        target_bed_uncompressed = BEDTOOLS_PAD_TARGET_BED.out.bed.map{ meta, inter -> inter}.collect().ifEmpty([[]])    // channel: [ val(meta), path(bed) ]
        target_bed              = ch_target_bed_gz_tbi.collect()                                                        // channel: [ val(meta), path(bed), path(tbi) ]
        bait_intervals          = CAT_CAT_BAIT.out.file_out.map{ meta, inter -> inter}.collect().ifEmpty([[]])          // channel: [ path(intervals) ]
        target_intervals        = GATK_BILT.out.interval_list.map{ meta, inter -> inter}.collect()                      // channel: [ path(interval_list) ]

}