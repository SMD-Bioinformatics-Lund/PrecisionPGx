process DEPTH {
    tag "${meta.group}"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.depth"),  emit: depth
        path "versions.yml",                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        samtools depth $args ${bam} > ${prefix}.depth

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n1 | sed 's/.*\s//g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.depth

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n1 | sed 's/.*\s//g')
        END_VERSIONS
        """
}