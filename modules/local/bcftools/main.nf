process BCFTOOLS_CALL {
    tag "${meta.group}"
    label "process_medium"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.bcftools.raw.vcf.gz"), file("*.bcftools.raw.vcf.gz.tbi"),  emit: bcftools_vcf
        path "versions.yml",                                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def args2   = task.ext.args2 ?: ''
        def args3   = task.ext.args3 ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}.bcftools.raw"
        """
        bcftools mpileup --threads ${task.cpus} $args ${bam} | bcftools call --threads ${task.cpus} ${args2} | bcftools filter ${args3} -o ${prefix}.vcf.gz
        tabix -p vcf ${prefix}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | head -n1 | sed 's/.*\s//g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}.bcftools.raw"
        """
        touch ${prefix}.vcf.gz
        touch ${prefix}.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | head -n1 | sed 's/.*\s//g')
        END_VERSIONS
        """
}