process FREEBAYES {
    label "process_medium"
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.freebayes.raw.vcf.gz"), file("*.freebayes.raw.vcf.gz.tbi"),    emit: freebayes_vcf
        path "versions.yml",                                                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args                ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}.freebayes.raw"
        """
        freebayes $args $bam > ${prefix}.vcf1

        filter_freebayes_unpaired.pl ${prefix}.vcf1 > ${prefix}.vcf

        bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
        tabix -p vcf ${prefix}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
    stub:
        def prefix      = task.ext.prefix ?: "${meta.group}.freebayes.raw"
        """
        touch ${prefix}.vcf.gz
        touch ${prefix}.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}