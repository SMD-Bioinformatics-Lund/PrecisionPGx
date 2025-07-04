#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK      } from '../../modules/local/check_input/main'


workflow CHECK_INPUT {
    take:
        csv     // file(csv)

    main:
        CSV_CHECK ( csv )
        checkedCsv = CSV_CHECK.out.csv.splitCsv( header:true, sep:',').set { csvmap }

        bam_fq_ch   = csvmap.map { create_input_channel(it) }
        meta        = csvmap.map { create_samples_channel(it) }

    emit:
        bam_fq_ch       // channel: [ val(group), val(meta), read1, read2. bam, bai ]
        meta            // channel: [ val(group), val(meta) ]

}

// Function to get list of [ group, meta, read1, read2, bam, bai  ]
def create_input_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.group              = row.group
    meta.id                 = row.id
    meta.type               = row.type
    meta.diagnosis          = (row.containsKey("diagnosis") ? row.diagnosis : '')
    meta.clarity_sample_id  = (row.containsKey("clarity_sample_id") ? row.clarity_sample_id : '')
    meta.sequencing_run     = (row.containsKey("sequencing_run") ? row.sequencing_run : '')
    meta.clarity_pool_id    = (row.containsKey("clarity_pool_id") ? row.clarity_pool_id : '')
    meta.ffpe               = (row.containsKey("ffpe") ? row.ffpe : false)
    meta.purity             = (row.containsKey("purity") ? row.purity : false)
	sub                     = false
	if (meta.reads && params.sample) {  
        if (meta.reads.toInteger() > params.sample_val) {
            sub = (params.sample_val / meta.reads.toInteger()).round(2)
            if (sub == 1.00){
                sub = 0.99
            }
        }
        else {
            sub = false
        }
	}
	meta.sub = sub
	// add path(s) of the fastq file(s) to the meta map
    // add path(s) of the bam file(s) to the meta map
    def input_ch = []
    if (row.containsKey("bam")) {
        input_ch = [row.group, meta, [], [], file(row.bam), file(row.bai)]
    } else if (row.containsKey("read1")) {
        input_ch = [row.group, meta, file(row.read1), file(row.read2), [], []]
    }

    return input_ch
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_samples_channel(LinkedHashMap row) {
    def meta                = [:]
    meta.id                 = row.id
    meta.group              = row.group
    meta.type               = row.type
    meta.diagnosis          = (row.containsKey("diagnosis") ? row.diagnosis : '')
    meta.clarity_sample_id  = (row.containsKey("clarity_sample_id") ? row.clarity_sample_id : '')
    meta.clarity_pool_id    = (row.containsKey("clarity_pool_id") ? row.clarity_pool_id : '')
    meta.sequencing_run     = (row.containsKey("sequencing_run") ? row.sequencing_run : '')
    meta.ffpe               = (row.containsKey("ffpe") ? row.ffpe : false)
    meta.purity             = (row.containsKey("purity") ? row.purity : false)
    meta.reads              = (row.containsKey("n_reads") ? row.n_reads : false)
    def sample_meta = [row.group, meta]
    return sample_meta
}
