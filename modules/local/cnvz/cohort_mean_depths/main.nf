process COHORT_MEAN_DEPTHS  {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/395d36614b9be92f68dc83cdf1658008dfcbe4e24eb51e91b772fc8e301cb66f/data' :
        'seqeracontainers/pandas-numpy-argparse' }"
    
    input:
        val(depthfiles)

    output:
        path("*.mean_proportional_depths.csv"),    emit: cohort_means
        path "versions.yml",        emit: versions

    script:
        //Passing $depthfiles directly results in command:
        //calculate_means.py -d [/labdata/KG/PrecisionPGx/PrecisionPGx_260529/work/0e/2408e997971e2f9c26894961ceaadf/NA12877.target_depth.tsv, /labdata/KG/PrecisionPGx/PrecisionPGx_260529/work/13/6081d2cb672b97e50a834c2edb7737/NA12878.target_depth.tsv] -p <prefix> -o .
        //Which will be interpreted as -d having two arguments. Expected input is comma-separated list with no spaces: path1,path2,path3
        // Use sed to remove the square brackets and the spaces: sed 's/[][" "]//g'

        def prefix = task.ext.prefix ?: "${params.trace_report_suffix}"

        """
        depthslist=`echo ${depthfiles} | sed 's/[][" "]//g'`

        echo "\$depthslist"

        calculate_means.py -d "\$depthslist" -p ${prefix} -o . 

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
