process TARGET_PASS_REGIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.pass.bed"), emit: pass_bed
    tuple val("${task.process}"), val('coreutils'), eval('echo \$VERSION'), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = "8.31"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '{if (\$3 >= ${params.min_dp}) print \$1 "\t" \$2-1 "\t" \$2 "\t" \$3 }' ${tsv} > ${prefix}.target.pass.bed  
    """

    stub:
    def VERSION = "8.31"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.target.pass.bed
    """
}
