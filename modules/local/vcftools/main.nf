process VCFTOOLS_FILTER {
    tag "${meta.group}"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.filtered.recode.vcf"), emit: filtered_vcf
        tuple val(group), val(meta), file("*.filtered.tagged.vcf"), emit: tagged_vcf
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def args2   = task.ext.args2 ?: ''
        def args3   = task.ext.args3 ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}.filtered"
        """ 
        vcffilter $args ${vcf} | vcffilter $args2 | vcfglxgt > ${prefix}.tagged.vcf

        vcftools --vcf ${prefix}.tagged.vcf --out ${prefix} $args3

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.tagged.vcf
        touch ${prefix}.recode.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
        END_VERSIONS
        """
}