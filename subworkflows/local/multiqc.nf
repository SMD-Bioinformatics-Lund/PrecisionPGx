#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC                       } from '../../modules/local/multiqc/main'

// Config File Channels //

ch_multiqc_config                       = Channel.fromPath(params.multiqc_config, checkIfExists: true)
ch_multiqc_custom_config                = params.multiqc_extra_config           ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true )  : Channel.empty()
ch_multiqc_logo                         = params.multiqc_logo                   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true )          : Channel.empty()
ch_multiqc_custom_methods_description   = params.multiqc_methods_desc           ? file(params.multiqc_methods_desc, checkIfExists: true)                : Channel.empty()

workflow RUN_MULTIQC {

    take:
        software_versions
        coverage_files
        fastqc_files
    
    main:
        ch_multiqc_files    = Channel.empty()
        ch_multiqc_files    = ch_multiqc_files.mix(software_versions)
        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            coverage_files
        )
        multiqc_reports = MULTIQC.out.report.toList()
    
    emit:
        multiqc_report          = multiqc_reports                               // channel: [ tuple val(group), val(meta) file("multiqc.html") ]
        multiqc_data            = MULTIQC.out.data                              // channel: [ tuple val(group), val(meta) file("multiqc_data") ]
        multiqc_plots           = MULTIQC.out.plots                             // channel: [ tuple val(group), val(meta) file("multiqc_plots") ]
}