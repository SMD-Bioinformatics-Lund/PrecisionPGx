process COHORT_MEAN_DEPTHS  {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/395d36614b9be92f68dc83cdf1658008dfcbe4e24eb51e91b772fc8e301cb66f/data' :
        'seqeracontainers/pandas-numpy-argparse' }"
    
    input:
        path(depthfiles)

    output:
        path("*.mean_proportional_depths.csv"),    emit: cohort_means
        path "versions.yml",        emit: versions

    script:
        def prefix = task.ext.prefix ?: "${params.trace_report_suffix}"

        """

        calculate_means.py -d ${depthfiles.join(',')} -p ${prefix} -o . 

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
        touch ${prefix}.mean_proportional_depths.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

}
