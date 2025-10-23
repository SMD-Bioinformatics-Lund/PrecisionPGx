process BGZIP {
    label 'process_low'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.haplotypes.vcf.gz"), file("*.haplotypes.vcf.gz.tbi"),  emit: compressed_vcf
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def args2       = task.ext.args2  ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}"
        """
        bgzip -c $vcf > ${prefix}.haplotypes.vcf.gz
        tabix ${prefix}.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.haplotypes.vcf.gz ${prefix}.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bgzip: \$(bgzip --version | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
}