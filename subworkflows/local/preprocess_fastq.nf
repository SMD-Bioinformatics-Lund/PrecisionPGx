#!/usr/bin/env nextflow

include { SUB_SAMPLE             } from '../../modules/local/seqtk/main'
include { TRIM                  } from '../../modules/local/fastp/main'


workflow PREPROCESS_FASTQ {
    take:
        fastq       // channel: [ val(meta), [ reads ] ]

    main:
        ch_versions = Channel.empty()

        // Only sub-sample if meta.sub == value
        SUB_SAMPLE ( fastq.filter{ item -> item[1].sub != false } )
        ch_versions = ch_versions.mix(SUB_SAMPLE.out.versions)


        //combine with any sample that does not get sub-sampled
        fastq_sample = SUB_SAMPLE.out.fastq_sub.mix( fastq.filter{ item -> item[1].sub == false } )

        TRIM ( fastq_sample )
        // fastq_done = FASTP.out.fastq_trimmed
        ch_versions = ch_versions.mix(TRIM.out.versions)


        if (params.trim ) {
            preprocessed_fastq = TRIM.out.fastq_trimmed
        } else {
            preprocessed_fastq = fastq_sample
        }

    emit:
        proccessed_fastq        =   preprocessed_fastq      // channel: [ val(meta), [ reads ] ] 
        versions                =   ch_versions             // channel: [ file(versions) ]

}