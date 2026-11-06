process GENERATE_REPORT  {
    tag "${meta.id}"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'docker://community.wave.seqera.io/library/pip_argparse_matplotlib_numpy_pandas:c7b55a761ab349d1' :
    //    'seqeracontainers/pip_argparse_matplotlib_numpy_pandas' }"
    container "community.wave.seqera.io/library/pip_argparse_matplotlib_numpy_pandas:c7b55a761ab349d1"
    
    input:
        tuple path(cnv), path(vcf), val(meta)
        path(target_bed)

    output:
        tuple val(meta), path("*.png"), path("*.txt"), emit: cn_report
        path "versions.yml"                        , emit: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args     ?: ''

        """
        generate_cnv_reports.py -s ${prefix} -cnv ${cnv} -vcf ${vcf} -bed ${target_bed} ${args} -o . 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        python --version
        python3 --version
        which python3

        touch ${prefix}.png
        touch ${prefix}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
}
