process CALCULATE_CN  {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/395d36614b9be92f68dc83cdf1658008dfcbe4e24eb51e91b772fc8e301cb66f/data' :
        'seqeracontainers/pandas-numpy-argparse' }"
    
    input:
        tuple val(meta), path(target_depths)
        path(mean_depths)

    output:
        tuple val(meta), path("*.csv"), emit: target_cn
        path "versions.yml", emit: versions

    script:
        
        def prefix = task.ext.prefix ?: "${meta.id}"

        """

        calculate_cn.py -s ${prefix} -d ${target_depths} -m ${mean_depths} -o . 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        python --version
        python3 --version
        which python3

        touch ${prefix}.cnvz.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

}
