include { CALCULATE_CN                       } from '../../modules/local/cnvz/calculate_cn'
include { GENERATE_REPORT                    } from '../../modules/local/cnvz/generate_report'

workflow CNV_CALLING {

    take:
        target_depth_tsv   // channel: [ val(meta), path(tsv) ]
        vcf                // channel: [ val(meta), path(vcf) ]
        ch_target_bed      // channel: [ val(meta), path(target_bed) ]
        cn_means           // path
    main:
        CALCULATE_CN (
            target_depth_tsv,
            cn_means
        )

        ch_generate_report_input = CALCULATE_CN.out.target_cn.map { meta, csv -> [meta.id, csv]}
        .join( 
            vcf.map { meta, vcf -> [meta.id, meta, vcf] },
            failOnMismatch: true, failOnDuplicate: true
        )
        .map { _id, meta, vcf, csv -> return [meta, csv, vcf] }
       
        //Debug prints
        //ch_generate_report_input.view { "ch_generate_report_input: $it "}
        //ch_target_bed.view { "ch_target_bed: $it "}
        //ch_generate_report_input.csv.view { "ch_generate_report_input.csv: $it "}
        //ch_generate_report_input.vcf.view { "ch_generate_report_input.vcf: $it "}
 
        GENERATE_REPORT (
            ch_generate_report_input,
            ch_target_bed.map { meta, target_bed -> return target_bed}
        )
       
    emit:
        target_cn                 = CALCULATE_CN.out.target_cn              // channel: [ val(meta), path(csv) ]
        cn_report                 = GENERATE_REPORT.out.cn_report           // channel: [ val(meta), path(png), path(txt)]
        versions                  = CALCULATE_CN.out.versions               // path(yml)
}
