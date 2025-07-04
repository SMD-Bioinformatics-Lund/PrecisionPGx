process MULTIQC {
    label 'process_single'
    tag "$meta.group"

    input:
        path multiqc_files, stageAs: "?/*"
        path(multiqc_config)
        path(extra_multiqc_config)
        path(multiqc_logo)
        tuple val(group), val(meta), file(targets_depth)

    output:
        tuple val(group), val(meta), file("*multiqc_report.html"), emit: report
        tuple val(group), val(meta), file("*_data")              , emit: data
        tuple val(group), val(meta), file("*_plots")             , optional:true, emit: plots
        path "versions.yml"        , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def config = multiqc_config ? "--config $multiqc_config" : ''
        def prefix = task.ext.prefix ?: "${meta.group}"
        def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
        def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''
        def sampleId = "${meta.group}" ?: ''
        """
        sed -i "s/REPLACE_SAMPLE_ID/$sampleId/" $multiqc_config

        multiqc \\
            --force \\
            --filename ${prefix}_multiqc_report.html \\
            $args \\
            $extra_config \\
            $config \\
            $logo . 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        mkdir ${prefix}_multiqc_data
        mkdir ${prefix}_multiqc_plots
        touch ${prefix}_multiqc_report.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
}