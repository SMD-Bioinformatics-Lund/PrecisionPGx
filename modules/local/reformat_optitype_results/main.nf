process REFORMAT_OPTITYPE_RESULTS  {
    tag "${meta.id}"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'docker://community.wave.seqera.io/library/pip_argparse_matplotlib_numpy_pandas:c7b55a761ab349d1' :
    //    'seqeracontainers/pip_argparse_matplotlib_numpy_pandas' }"
    container "community.wave.seqera.io/library/pip_argparse_matplotlib_numpy_pandas:c7b55a761ab349d1"
    
    input:
        tuple val(meta), path(tsv)

    output:
        tuple val(meta), path("*.txt"), emit: hla_calls
        path "versions.yml"           , emit: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def args    = task.ext.args     ?: ''

        """
        optitype_to_pharmcat.py ${prefix} ${tsv} 

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

        touch ${prefix}_HLA.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """
}
