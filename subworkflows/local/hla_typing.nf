include { YARA_INDEX                   } from '../../modules/nf-core/yara/index/main'
include { YARA_MAPPER                  } from '../../modules/nf-core/yara/mapper/main'
include { OPTITYPE                     } from '../../modules/nf-core/optitype/main'
include { REFORMAT_OPTITYPE_RESULTS    } from '../../modules/local/reformat_optitype_results/main'

workflow HLA_TYPING {
    take:
    ch_reads_input   // channel: [ val(meta), path(reads) ]
    hla_fasta     // path(fasta)

    main:
        ch_versions = channel.empty()
        ch_input_with_references = ch_reads_input.map { meta, _reads -> [ meta, file("$hla_fasta") ]
            }

        ch_input_with_references.view { "input to YARA_INDEX: $it "}        

        // MODULE: Run Yara indexing on HLA reference
        YARA_INDEX ( ch_input_with_references )
        ch_versions = ch_versions.mix(YARA_INDEX.out.versions_yara)
        // Map sample-specific reads and index
        ch_mapping_input = ch_reads_input.cross(YARA_INDEX.out.index).multiMap { reads, index -> 
            reads: reads
            index: index 
            }
        // MODULE: Run Yara mapping. Preparation Step - Pre-mapping against HLA
        YARA_MAPPER ( ch_mapping_input.reads, ch_mapping_input.index)
        ch_versions = ch_versions.mix(YARA_MAPPER.out.versions_yara)
        // MODULE: OptiType
        OPTITYPE ( YARA_MAPPER.out.bam.join(YARA_MAPPER.out.bai) )
        ch_versions = ch_versions.mix(OPTITYPE.out.versions_optitype)
        
        //Reformat OptiType output to match PharmCAT outside call format
        REFORMAT_OPTITYPE_RESULTS ( OPTITYPE.out.hla_type )
        ch_versions = ch_versions.mix(REFORMAT_OPTITYPE_RESULTS.out.versions)

    emit:
        hla_type         = OPTITYPE.out.hla_type                      // channel: [ val(meta), path(tsv) ]
        coverage_plot    = OPTITYPE.out.coverage_plot                 // channel: [ val(meta), path(pdf) ]
        hla_calls        = REFORMAT_OPTITYPE_RESULTS.out.hla_calls    // channel: [ val(meta), path(txt) ]
        versions         = ch_versions                                // channel: [ path(versions.yml) ]
}
